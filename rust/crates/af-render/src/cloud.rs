//! Out-of-core atom clouds (af_cloud .afcloud files) of up to billions of
//! atoms (the C++ graphics/CloudRenderer). Atoms are sphere impostors (lit,
//! with correct depth), 8 bytes each on the GPU. Every frame each chunk is
//! culled against the view and given a number of atoms from its size on
//! screen; as chunks are stored in random order, drawing a prefix draws a
//! uniform subsample. A total point budget keeps the frame rate, a memory
//! budget bounds the GPU memory, and a background thread streams the needed
//! prefixes from the file.
//!
//! The planning (culling, level of detail, budgets, streaming requests,
//! eviction) is [`Planner`], plain CPU code ported line by line; the drawing is
//! wgpu: instanced quads instead of GL point sprites, occlusion queries on the
//! chunk boxes (wgpu's are boolean, so a box counts as hidden when no sample
//! passed, where the C++ took fewer than 16), no GPU timer.

use std::collections::VecDeque;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;

use af_cloud::{CloudFile, CloudInfo, PackedAtom};
use bytemuck::{Pod, Zeroable};
use glam::{Mat4, Vec3, Vec4};

/// Atoms per load granule.
pub const BLOCK: u64 = 1 << 16;
/// Atoms per load request.
pub const MAX_LOAD: u64 = 4 << 20;
pub const MAX_IN_FLIGHT: usize = 12;
pub const UPLOAD_BYTES_PER_FRAME: f64 = 384.0 * 1024.0 * 1024.0;
const ATOM_BYTES: u64 = std::mem::size_of::<PackedAtom>() as u64;

const COLOUR_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;
const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CloudSettings {
    /// Atoms drawn per frame at most.
    pub point_budget: f64,
    /// GPU memory for atom data.
    pub memory_budget_mb: f64,
    /// Atoms per atom-sized screen area of a chunk (level of detail).
    pub density: f32,
    pub radius_scale: f32,
    /// false: every atom of visible chunks (within the budgets).
    pub level_of_detail: bool,
    pub show_box: bool,
    /// Skip chunks hidden behind drawn atoms.
    pub occlusion_culling: bool,
}

impl Default for CloudSettings {
    fn default() -> Self {
        Self { point_budget: 60e6, memory_budget_mb: 2560.0, density: 3.0, radius_scale: 1.0, level_of_detail: true, show_box: true, occlusion_culling: true }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CloudStats {
    pub total_atoms: u64,
    pub drawn_atoms: u64,
    pub resident_atoms: u64,
    pub chunks: usize,
    pub visible_chunks: usize,
    pub occluded_chunks: usize,
    pub loading_chunks: usize,
    pub resident_mb: f64,
    pub uploaded_mb_this_frame: f64,
    /// Every visible chunk has the atoms it wants.
    pub complete: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ChunkState {
    /// Atoms the GPU buffer holds.
    pub capacity: u64,
    /// Leading atoms uploaded.
    pub resident: u64,
    /// Atoms wanted this frame.
    pub wanted: u64,
    /// Prefix requested from the loader.
    pub requested: u64,
    /// Frame it was last drawn.
    pub last_used: u64,
}

/// A chunk to draw this frame: `cover` atoms fill its screen area; `want` is
/// what it gets (all atoms once they are about a pixel or larger on screen).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Visible {
    pub chunk: usize,
    pub want: f64,
    pub cover: f64,
    pub distance: f32,
    pub area: f32,
}

/// Atoms [first, first + count) of a chunk, to read from the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LoadRequest {
    pub chunk: usize,
    pub first: u64,
    pub count: u64,
}

/// The CPU side of the renderer: per-chunk residency, the frame's plan and
/// the streaming bookkeeping (CloudRenderer::plan, request, evict, upload).
#[derive(Clone, Debug, Default)]
pub struct Planner {
    pub chunks: Vec<ChunkState>,
    /// Hidden behind atoms drawn in an earlier frame (occlusion queries).
    pub occluded: Vec<bool>,
    /// The chunks to draw this frame, front to back.
    pub visible: Vec<Visible>,
    /// In the frustum and not behind the camera: occlusion-tested.
    pub tested: Vec<usize>,
    pub frame: u64,
    pub resident_bytes: u64,
    pub in_flight: usize,
    pub complete: bool,
    pub pixels_per_unit: f32,
    pub max_radius: f32,
}

impl Planner {
    pub fn new(chunks: usize) -> Self {
        Self { chunks: vec![ChunkState::default(); chunks], occluded: vec![false; chunks], pixels_per_unit: 1.0, max_radius: 0.3, ..Default::default() }
    }

    /// Requests a chunk's atoms up to `up_to`, in pieces, so that one large
    /// chunk does not hold up the others.
    pub fn request(&mut self, info: &CloudInfo, chunk: usize, up_to: u64, out: &mut Vec<LoadRequest>) {
        let count = info.chunks[chunk].count;
        let c = &mut self.chunks[chunk];
        let up_to = count.min(up_to);
        if up_to <= c.requested {
            return;
        }
        let mut first = c.requested;
        while first < up_to {
            out.push(LoadRequest { chunk, first, count: MAX_LOAD.min(up_to - first) });
            first += MAX_LOAD;
        }
        self.in_flight += (up_to - c.requested).div_ceil(MAX_LOAD) as usize;
        c.requested = up_to;
    }

    /// A finished load: the chunk's buffer capacity to upload it into, or
    /// None when it is out of order (after an eviction) or failed, and is
    /// forgotten (requested again later).
    pub fn accept_load(&mut self, info: &CloudInfo, load: &LoadRequest, ok: bool) -> Option<u64> {
        let c = &mut self.chunks[load.chunk];
        if !ok || load.first != c.resident {
            c.requested = c.resident;
            return None;
        }
        let end = load.first + load.count;
        let mut capacity = c.capacity;
        if end > c.capacity {
            let total = info.chunks[load.chunk].count;
            capacity = total.min(end.max(c.capacity * 2));
            self.resident_bytes += (capacity - c.capacity) * ATOM_BYTES;
        }
        Some(capacity)
    }

    /// The upload of an accepted load is done.
    pub fn commit_load(&mut self, load: &LoadRequest, capacity: u64) {
        let c = &mut self.chunks[load.chunk];
        c.capacity = capacity;
        c.resident = load.first + load.count;
    }

    /// An accepted load could not be uploaded: forgotten as a failed one.
    pub fn reject_load(&mut self, load: &LoadRequest, capacity: u64) {
        let c = &mut self.chunks[load.chunk];
        self.resident_bytes -= (capacity - c.capacity) * ATOM_BYTES;
        c.requested = c.resident;
    }

    /// Releases the least recently drawn chunks beyond the memory budget
    /// (never one with loads in flight); returns them.
    pub fn evict(&mut self, settings: &CloudSettings, frame: u64) -> Vec<usize> {
        let budget = settings.memory_budget_mb * 1024.0 * 1024.0;
        if self.resident_bytes as f64 <= budget {
            return Vec::new();
        }
        let mut candidates: Vec<usize> =
            (0..self.chunks.len()).filter(|&i| self.chunks[i].capacity != 0 && self.chunks[i].last_used < frame && self.chunks[i].requested == self.chunks[i].resident).collect();
        af_core::cxx::sort_by(&mut candidates, |&a, &b| self.chunks[a].last_used < self.chunks[b].last_used);
        let mut released = Vec::new();
        for i in candidates {
            if self.resident_bytes as f64 <= budget * 0.9 {
                break;
            }
            self.resident_bytes -= self.chunks[i].capacity * ATOM_BYTES;
            self.chunks[i] = ChunkState::default();
            released.push(i);
        }
        released
    }

    /// Culling, level of detail, budgets and streaming of a frame. Matrices
    /// are glm's (clip depth -1..1). Returns the chunks evicted.
    #[allow(clippy::too_many_arguments)]
    pub fn plan(&mut self, info: &CloudInfo, settings: &CloudSettings, projection: &Mat4, view: &Mat4, width: i32, height: i32, orthographic: bool, requests: &mut Vec<LoadRequest>) -> Vec<usize> {
        // Visibility and level of detail per chunk.
        let view_projection = *projection * *view;
        let pixels_per_unit = projection.y_axis.y * height as f32 * 0.5;
        let mut max_radius = 0.3f32;
        for s in &info.species {
            max_radius = max_radius.max(s.radius);
        }
        max_radius *= settings.radius_scale;
        self.pixels_per_unit = pixels_per_unit;
        self.max_radius = max_radius;
        let screen_area = width as f64 * height as f64;
        let mut visible: Vec<Visible> = Vec::with_capacity(self.chunks.len());
        let mut tested: Vec<usize> = Vec::with_capacity(self.chunks.len());
        let mut total_wanted = 0.0f64;
        for (i, chunk) in info.chunks.iter().enumerate() {
            let lo = Vec3::new(chunk.lower[0] as f32, chunk.lower[1] as f32, chunk.lower[2] as f32);
            let hi = Vec3::new(chunk.upper[0] as f32, chunk.upper[1] as f32, chunk.upper[2] as f32);
            let pad = Vec3::splat(max_radius);
            let (a, b) = (lo - pad, hi + pad);
            // Frustum test on the box corners (clip space), with screen extent.
            let mut outside = [0; 6];
            let mut behind = false;
            let (mut min_x, mut min_y, mut max_x, mut max_y) = (1e30f32, 1e30f32, -1e30f32, -1e30f32);
            for k in 0..8 {
                let corner = Vec4::new(if k & 1 != 0 { b.x } else { a.x }, if k & 2 != 0 { b.y } else { a.y }, if k & 4 != 0 { b.z } else { a.z }, 1.0);
                let clip = view_projection * corner;
                outside[0] += (clip.x < -clip.w) as i32;
                outside[1] += (clip.x > clip.w) as i32;
                outside[2] += (clip.y < -clip.w) as i32;
                outside[3] += (clip.y > clip.w) as i32;
                outside[4] += (clip.z < -clip.w) as i32;
                outside[5] += (clip.z > clip.w) as i32;
                if clip.w <= 1e-6 {
                    behind = true;
                    continue;
                }
                let (x, y) = (clip.x / clip.w, clip.y / clip.w);
                min_x = min_x.min(x);
                max_x = max_x.max(x);
                min_y = min_y.min(y);
                max_y = max_y.max(y);
            }
            if outside.iter().any(|&n| n == 8) {
                self.occluded[i] = false;
                continue;
            }
            if !behind {
                tested.push(i);
            } else {
                self.occluded[i] = false;
            }
            // Hidden behind atoms drawn last frame: neither drawn nor streamed.
            if settings.occlusion_culling && self.occluded[i] {
                continue;
            }
            let mut area = screen_area;
            if !behind {
                let w = (max_x.clamp(-1.0, 1.0) - min_x.clamp(-1.0, 1.0)) as f64 * 0.5 * width as f64;
                let h = (max_y.clamp(-1.0, 1.0) - min_y.clamp(-1.0, 1.0)) as f64 * 0.5 * height as f64;
                area = (w * h).max(1.0);
            }
            let centre = 0.5 * (lo + hi);
            let distance = (-(*view * Vec4::from((centre, 1.0))).z).max(1e-3);
            let count = chunk.count as f64;
            let (mut want, mut cover) = (count, count);
            if settings.level_of_detail && !behind {
                let half_diagonal = 0.5 * (hi - lo).length();
                let near = (distance - half_diagonal).max(1e-3);
                let far = distance + half_diagonal;
                let near_pixels = (max_radius * pixels_per_unit / if orthographic { 1.0 } else { near }) as f64;
                let far_pixels = (max_radius * pixels_per_unit / if orthographic { 1.0 } else { far }) as f64;
                // Enough atoms to cover the chunk's screen area `density` times over...
                let atom_area = 3.14159265 * near_pixels.max(0.5) * near_pixels.max(0.5);
                cover = count.min(256.0f64.max(settings.density as f64 * area / atom_area));
                // ...rising to every atom as atoms grow from 0.75 to 1.5 pixels in radius,
                // where the arrangement of atoms becomes visible.
                let resolved = ((far_pixels - 0.75) / 0.75).clamp(0.0, 1.0);
                want = cover + (count - cover) * resolved;
            }
            visible.push(Visible { chunk: i, want, cover, distance, area: area as f32 });
            total_wanted += want;
        }
        // Point and memory budgets: every chunk keeps its cover; the remaining
        // budget gives the nearest chunks full detail first.
        let memory_atoms = settings.memory_budget_mb * 1024.0 * 1024.0 / ATOM_BYTES as f64 * 0.85;
        let limit = settings.point_budget.min(memory_atoms);
        if total_wanted > limit {
            let covers: f64 = visible.iter().map(|v| v.want.min(v.cover)).sum();
            if covers > limit {
                let scale = limit / covers;
                for v in visible.iter_mut() {
                    v.want = 256.0f64.max(v.want.min(v.cover) * scale);
                }
            } else {
                let mut order: Vec<usize> = (0..visible.len()).collect();
                af_core::cxx::sort_by(&mut order, |&x, &y| visible[x].distance < visible[y].distance);
                let mut remaining = limit - covers;
                for k in order {
                    let v = &mut visible[k];
                    let base = v.want.min(v.cover);
                    let extra = (v.want - base).min(remaining);
                    v.want = base + extra;
                    remaining -= extra;
                }
            }
        }
        // Stream what is missing: the largest gaps on screen first.
        let mut missing: Vec<(f64, usize)> = Vec::new();
        let mut complete = true;
        for v in &visible {
            let c = &mut self.chunks[v.chunk];
            c.wanted = v.want as u64;
            c.last_used = self.frame;
            if c.resident < c.wanted {
                complete = false;
            }
            if c.requested < c.wanted {
                missing.push((v.area as f64 * (1.0 - c.resident as f64 / v.want.max(1.0)), v.chunk));
            }
        }
        af_core::cxx::sort_by(&mut missing, |x, y| x.0 > y.0);
        for &(_, chunk) in &missing {
            if self.in_flight >= MAX_IN_FLIGHT {
                break;
            }
            let target = self.chunks[chunk].wanted.div_ceil(BLOCK) * BLOCK;
            self.request(info, chunk, target, requests);
        }
        let evicted = self.evict(settings, self.frame);
        self.complete = complete;
        // Draw front to back (cheap depth rejection of hidden atoms).
        af_core::cxx::sort_by(&mut visible, |x, y| x.distance < y.distance);
        self.visible = visible;
        self.tested = tested;
        evicted
    }

    /// Atoms of each planned chunk to draw now: the resident prefix of its want.
    pub fn draw_count(&self, v: &Visible) -> u64 {
        self.chunks[v.chunk].resident.min(v.want.ceil() as u64)
    }
}

/// The cloud's outline: the cell when known, else the bounding box (24 line ends).
pub fn outline(info: &CloudInfo) -> Vec<[f32; 3]> {
    let v = |x: [f64; 3]| Vec3::new(x[0] as f32, x[1] as f32, x[2] as f32);
    let (o, a, b, c) = if info.has_cell {
        (v(info.cell_origin), v(info.cell[0]), v(info.cell[1]), v(info.cell[2]))
    } else {
        (
            v(info.lower),
            Vec3::new((info.upper[0] - info.lower[0]) as f32, 0.0, 0.0),
            Vec3::new(0.0, (info.upper[1] - info.lower[1]) as f32, 0.0),
            Vec3::new(0.0, 0.0, (info.upper[2] - info.lower[2]) as f32),
        )
    };
    let corner = [o, o + a, o + b, o + a + b, o + c, o + a + c, o + b + c, o + a + b + c];
    const EDGES: [[usize; 2]; 12] = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];
    EDGES.iter().flat_map(|e| [corner[e[0]].to_array(), corner[e[1]].to_array()]).collect()
}

/// The centre of a cloud's bounding box.
pub fn cloud_centre(info: &CloudInfo) -> Vec3 {
    Vec3::new((0.5 * (info.lower[0] + info.upper[0])) as f32, (0.5 * (info.lower[1] + info.upper[1])) as f32, (0.5 * (info.lower[2] + info.upper[2])) as f32)
}

/// Half the diagonal of a cloud's bounding box (at least 1).
pub fn cloud_radius(info: &CloudInfo) -> f32 {
    let d = Vec3::new((info.upper[0] - info.lower[0]) as f32, (info.upper[1] - info.lower[1]) as f32, (info.upper[2] - info.lower[2]) as f32);
    1.0f32.max(0.5 * d.length())
}

// ---------------------------------------------------------------- loader

#[derive(Default)]
struct Queues {
    queue: VecDeque<LoadRequest>,
    done: VecDeque<(LoadRequest, Vec<PackedAtom>)>,
    stop: bool,
}

/// The background thread reading requested chunk prefixes from the file.
struct Loader {
    shared: Arc<(Mutex<Queues>, Condvar)>,
    thread: Option<JoinHandle<()>>,
}

impl Loader {
    fn start(file: Arc<CloudFile>) -> Self {
        let shared: Arc<(Mutex<Queues>, Condvar)> = Arc::default();
        let worker = shared.clone();
        let thread = std::thread::spawn(move || loop {
            let load = {
                let (lock, wake) = &*worker;
                let mut q = lock.lock().unwrap();
                while !q.stop && q.queue.is_empty() {
                    q = wake.wait(q).unwrap();
                }
                if q.stop {
                    return;
                }
                q.queue.pop_front().unwrap()
            };
            // A read error drops the load: requested again later.
            let atoms = file.read(load.chunk, load.first, load.count).unwrap_or_default();
            worker.0.lock().unwrap().done.push_back((load, atoms));
        });
        Self { shared, thread: Some(thread) }
    }

    fn push(&self, loads: &[LoadRequest]) {
        if loads.is_empty() {
            return;
        }
        self.shared.0.lock().unwrap().queue.extend(loads.iter().copied());
        self.shared.1.notify_one();
    }

    fn pop_done(&self) -> Option<(LoadRequest, Vec<PackedAtom>)> {
        self.shared.0.lock().unwrap().done.pop_front()
    }
}

impl Drop for Loader {
    fn drop(&mut self) {
        self.shared.0.lock().unwrap().stop = true;
        self.shared.1.notify_all();
        if let Some(t) = self.thread.take() {
            let _ = t.join();
        }
    }
}

// ---------------------------------------------------------------- GPU

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct FrameUniform {
    view: [f32; 16],
    projection: [f32; 16],
    light: [f32; 4],
    params: [f32; 4],
    colour_range: [f32; 4],
    line_colour: [f32; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ChunkUniform {
    lower: [f32; 4],
    scale: [f32; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct BoxInstance {
    a: [f32; 3],
    b: [f32; 3],
}

/// Occlusion results being read back (one batch at a time).
struct Readback {
    chunks: Vec<usize>,
    ready: Arc<AtomicBool>,
    /// map_async was called.
    mapping: bool,
}

/// Renders an atom cloud offscreen (the desktop shows the texture).
pub struct CloudRenderer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    point_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,
    box_pipeline: wgpu::RenderPipeline,
    frame_layout: wgpu::BindGroupLayout,
    chunk_layout: wgpu::BindGroupLayout,
    frame_uniform: wgpu::Buffer,
    frame_group: Option<wgpu::BindGroup>,
    /// Per-chunk uniforms (dynamic offsets), their capacity in chunks.
    chunk_uniform: Option<(wgpu::Buffer, wgpu::BindGroup, usize)>,
    uniform_stride: u64,
    unit_box: wgpu::Buffer,
    line_buffer: Option<wgpu::Buffer>,
    queries: Option<(wgpu::QuerySet, wgpu::Buffer, wgpu::Buffer)>,
    readback: Option<Readback>,
    file: Option<Arc<CloudFile>>,
    path: PathBuf,
    buffers: Vec<Option<wgpu::Buffer>>,
    planner: Planner,
    loader: Option<Loader>,
    target: Option<(u32, u32, wgpu::Texture, wgpu::Texture)>,
    stats: CloudStats,
    pub settings: CloudSettings,
    /// Colours by species, or by height along an axis (0 none, 1 x, 2 y, 3 z).
    pub colour_axis: i32,
}

impl CloudRenderer {
    pub fn with_device(device: wgpu::Device, queue: wgpu::Queue) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("cloud"), source: wgpu::ShaderSource::Wgsl(include_str!("cloud.wgsl").into()) });
        let frame_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("cloud frame"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                    ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::VERTEX,
                    ty: wgpu::BindingType::Texture { sample_type: wgpu::TextureSampleType::Float { filterable: false }, view_dimension: wgpu::TextureViewDimension::D2, multisampled: false },
                    count: None,
                },
            ],
        });
        let chunk_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("cloud chunk"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: true,
                    min_binding_size: wgpu::BufferSize::new(std::mem::size_of::<ChunkUniform>() as u64),
                },
                count: None,
            }],
        });
        let points_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("cloud points"),
            bind_group_layouts: &[Some(&frame_layout), Some(&chunk_layout)],
            ..Default::default()
        });
        let frame_only = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("cloud lines"), bind_group_layouts: &[Some(&frame_layout)], ..Default::default() });
        let depth = |write: bool| {
            Some(wgpu::DepthStencilState {
                format: DEPTH_FORMAT,
                depth_write_enabled: Some(write),
                depth_compare: Some(wgpu::CompareFunction::Less),
                stencil: Default::default(),
                bias: Default::default(),
            })
        };
        let colour = [Some(wgpu::ColorTargetState { format: COLOUR_FORMAT, blend: None, write_mask: wgpu::ColorWrites::ALL })];
        let no_colour = [Some(wgpu::ColorTargetState { format: COLOUR_FORMAT, blend: None, write_mask: wgpu::ColorWrites::empty() })];
        let point_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("cloud points"),
            layout: Some(&points_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("point_vs"),
                buffers: &[Some(wgpu::VertexBufferLayout {
                    array_stride: ATOM_BYTES,
                    step_mode: wgpu::VertexStepMode::Instance,
                    attributes: &wgpu::vertex_attr_array![0 => Uint16x4],
                })],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState { module: &shader, entry_point: Some("point_fs"), targets: &colour, compilation_options: Default::default() }),
            primitive: wgpu::PrimitiveState { topology: wgpu::PrimitiveTopology::TriangleStrip, ..Default::default() },
            depth_stencil: depth(true),
            multisample: Default::default(),
            multiview_mask: None,
            cache: None,
        });
        let line_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("cloud box"),
            layout: Some(&frame_only),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("line_vs"),
                buffers: &[Some(wgpu::VertexBufferLayout { array_stride: 12, step_mode: wgpu::VertexStepMode::Vertex, attributes: &wgpu::vertex_attr_array![0 => Float32x3] })],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState { module: &shader, entry_point: Some("line_fs"), targets: &colour, compilation_options: Default::default() }),
            primitive: wgpu::PrimitiveState { topology: wgpu::PrimitiveTopology::LineList, ..Default::default() },
            depth_stencil: depth(true),
            multisample: Default::default(),
            multiview_mask: None,
            cache: None,
        });
        let box_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("cloud occlusion boxes"),
            layout: Some(&frame_only),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("box_vs"),
                buffers: &[
                    Some(wgpu::VertexBufferLayout { array_stride: 12, step_mode: wgpu::VertexStepMode::Vertex, attributes: &wgpu::vertex_attr_array![0 => Float32x3] }),
                    Some(wgpu::VertexBufferLayout {
                        array_stride: std::mem::size_of::<BoxInstance>() as u64,
                        step_mode: wgpu::VertexStepMode::Instance,
                        attributes: &wgpu::vertex_attr_array![1 => Float32x3, 2 => Float32x3],
                    }),
                ],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState { module: &shader, entry_point: Some("box_fs"), targets: &no_colour, compilation_options: Default::default() }),
            primitive: wgpu::PrimitiveState { topology: wgpu::PrimitiveTopology::TriangleList, ..Default::default() },
            depth_stencil: depth(false),
            multisample: Default::default(),
            multiview_mask: None,
            cache: None,
        });
        let frame_uniform = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("cloud frame"),
            size: std::mem::size_of::<FrameUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        // Unit cube (triangles) for the chunk occlusion queries.
        const CORNERS: [[f32; 3]; 8] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.], [1., 1., 1.]];
        const FACES: [[usize; 3]; 12] = [[0, 2, 1], [1, 2, 3], [4, 5, 6], [5, 7, 6], [0, 1, 4], [1, 5, 4], [2, 6, 3], [3, 6, 7], [0, 4, 2], [2, 4, 6], [1, 3, 5], [3, 7, 5]];
        let triangles: Vec<[f32; 3]> = FACES.iter().flat_map(|f| f.iter().map(|&v| CORNERS[v])).collect();
        let unit_box = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("unit box"),
            size: std::mem::size_of_val(triangles.as_slice()) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&unit_box, 0, bytemuck::cast_slice(&triangles));
        let alignment = device.limits().min_uniform_buffer_offset_alignment as u64;
        let uniform_stride = (std::mem::size_of::<ChunkUniform>() as u64).div_ceil(alignment) * alignment;
        Self {
            device,
            queue,
            point_pipeline,
            line_pipeline,
            box_pipeline,
            frame_layout,
            chunk_layout,
            frame_uniform,
            frame_group: None,
            chunk_uniform: None,
            uniform_stride,
            unit_box,
            line_buffer: None,
            queries: None,
            readback: None,
            file: None,
            path: PathBuf::new(),
            buffers: Vec::new(),
            planner: Planner::default(),
            loader: None,
            target: None,
            stats: CloudStats::default(),
            settings: CloudSettings::default(),
            colour_axis: 0,
        }
    }

    /// Opens a cloud for display.
    pub fn open(&mut self, path: &Path) -> Result<(), String> {
        self.close();
        let file = Arc::new(CloudFile::open(path)?);
        let info = file.info().clone();
        let n = info.chunks.len();
        self.planner = Planner::new(n);
        self.buffers = (0..n).map(|_| None).collect();
        self.stats = CloudStats { total_atoms: info.atoms, chunks: n, ..Default::default() };
        // Species colours and radii.
        let rows = info.species.len().div_ceil(256).max(1);
        let mut texels = vec![0.7f32; 256 * rows * 4];
        for (i, s) in info.species.iter().enumerate() {
            texels[i * 4..i * 4 + 4].copy_from_slice(&[s.color[0], s.color[1], s.color[2], s.radius]);
        }
        let species = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("cloud species"),
            size: wgpu::Extent3d { width: 256, height: rows as u32, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba32Float,
            usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
            view_formats: &[],
        });
        self.queue.write_texture(
            wgpu::TexelCopyTextureInfo { texture: &species, mip_level: 0, origin: wgpu::Origin3d::ZERO, aspect: wgpu::TextureAspect::All },
            bytemuck::cast_slice(&texels),
            wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(256 * 16), rows_per_image: Some(rows as u32) },
            wgpu::Extent3d { width: 256, height: rows as u32, depth_or_array_layers: 1 },
        );
        let species_view = species.create_view(&Default::default());
        self.frame_group = Some(self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("cloud frame"),
            layout: &self.frame_layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: self.frame_uniform.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: wgpu::BindingResource::TextureView(&species_view) },
            ],
        }));
        let lines = outline(&info);
        let line_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("cloud outline"),
            size: std::mem::size_of_val(lines.as_slice()) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.queue.write_buffer(&line_buffer, 0, bytemuck::cast_slice(&lines));
        self.line_buffer = Some(line_buffer);
        if n > 0 {
            let set = self.device.create_query_set(&wgpu::QuerySetDescriptor { label: Some("cloud occlusion"), ty: wgpu::QueryType::Occlusion, count: n as u32 });
            let size = (n as u64 * 8).max(8);
            let resolve = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("occlusion resolve"),
                size,
                usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            });
            let read = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("occlusion read"),
                size,
                usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
                mapped_at_creation: false,
            });
            self.queries = Some((set, resolve, read));
        }
        self.loader = Some(Loader::start(file.clone()));
        self.file = Some(file);
        self.path = path.to_path_buf();
        Ok(())
    }

    pub fn close(&mut self) {
        self.loader = None;
        self.readback = None;
        self.queries = None;
        self.buffers.clear();
        self.planner = Planner::default();
        self.file = None;
        self.frame_group = None;
        self.line_buffer = None;
        self.stats = CloudStats::default();
    }

    pub fn active(&self) -> bool {
        self.file.is_some()
    }

    pub fn info(&self) -> Option<&CloudInfo> {
        self.file.as_ref().map(|f| f.info())
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn stats(&self) -> &CloudStats {
        &self.stats
    }

    pub fn centre(&self) -> Vec3 {
        self.info().map(cloud_centre).unwrap_or(Vec3::ZERO)
    }

    pub fn radius(&self) -> f32 {
        self.info().map(cloud_radius).unwrap_or(1.0)
    }

    /// Finished loads, up to a per-frame upload volume, onto the GPU.
    fn consume_loads(&mut self, encoder: &mut wgpu::CommandEncoder) {
        let (Some(file), Some(loader)) = (self.file.clone(), self.loader.as_ref()) else { return };
        let info = file.info();
        let mut bytes = 0.0f64;
        let mut finished = Vec::new();
        while bytes <= UPLOAD_BYTES_PER_FRAME {
            let Some(done) = loader.pop_done() else { break };
            bytes += (done.0.count * ATOM_BYTES) as f64;
            finished.push(done);
        }
        for (load, atoms) in finished {
            self.planner.in_flight = self.planner.in_flight.saturating_sub(1);
            let ok = atoms.len() as u64 == load.count && load.count > 0;
            let Some(capacity) = self.planner.accept_load(info, &load, ok) else { continue };
            let old = self.planner.chunks[load.chunk].capacity;
            if capacity != old {
                let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("cloud chunk"),
                    size: capacity * ATOM_BYTES,
                    usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
                    mapped_at_creation: false,
                });
                let resident = self.planner.chunks[load.chunk].resident;
                if let Some(previous) = &self.buffers[load.chunk] {
                    if resident > 0 {
                        encoder.copy_buffer_to_buffer(previous, 0, &buffer, 0, resident * ATOM_BYTES);
                    }
                }
                self.buffers[load.chunk] = Some(buffer);
            }
            let Some(buffer) = &self.buffers[load.chunk] else {
                self.planner.reject_load(&load, capacity);
                continue;
            };
            // Plain repr(C) u16 records.
            let data = unsafe { std::slice::from_raw_parts(atoms.as_ptr() as *const u8, atoms.len() * ATOM_BYTES as usize) };
            self.queue.write_buffer(buffer, load.first * ATOM_BYTES, data);
            self.planner.commit_load(&load, capacity);
            self.stats.uploaded_mb_this_frame += (load.count * ATOM_BYTES) as f64 / (1024.0 * 1024.0);
        }
    }

    /// Occlusion results of an earlier frame, when already available (no stall).
    fn read_occlusion(&mut self) {
        let Some(r) = &self.readback else { return };
        let _ = self.device.poll(wgpu::PollType::Poll);
        if !r.ready.load(Ordering::Acquire) {
            return;
        }
        let r = self.readback.take().unwrap();
        let Some((_, _, read)) = &self.queries else { return };
        {
            let Ok(data) = read.get_mapped_range(..) else {
                read.unmap();
                return;
            };
            for (k, &chunk) in r.chunks.iter().enumerate() {
                let samples = u64::from_le_bytes(data[k * 8..k * 8 + 8].try_into().unwrap());
                if chunk < self.planner.occluded.len() {
                    self.planner.occluded[chunk] = samples == 0;
                }
            }
        }
        read.unmap();
    }

    fn ensure_target(&mut self, width: u32, height: u32) {
        if self.target.as_ref().is_some_and(|t| t.0 == width && t.1 == height) {
            return;
        }
        let extent = wgpu::Extent3d { width, height, depth_or_array_layers: 1 };
        let texture = |format, usage| {
            self.device.create_texture(&wgpu::TextureDescriptor {
                label: Some("cloud target"),
                size: extent,
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage,
                view_formats: &[],
            })
        };
        let colour = texture(COLOUR_FORMAT, wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC | wgpu::TextureUsages::TEXTURE_BINDING);
        let depth = texture(DEPTH_FORMAT, wgpu::TextureUsages::RENDER_ATTACHMENT);
        self.target = Some((width, height, colour, depth));
    }

    /// Plans, streams and draws one frame into the target texture; returns
    /// the encoder (not yet submitted). Matrices are glm's.
    #[allow(clippy::too_many_arguments)]
    fn encode(&mut self, projection: &Mat4, view: &Mat4, width: u32, height: u32, orthographic: bool, light_theme: bool, background: [f32; 4]) -> Option<wgpu::CommandEncoder> {
        let file = self.file.clone()?;
        let info = file.info();
        self.planner.frame += 1;
        self.stats.uploaded_mb_this_frame = 0.0;
        let mut encoder = self.device.create_command_encoder(&Default::default());
        self.consume_loads(&mut encoder);
        self.read_occlusion();
        let mut requests = Vec::new();
        let evicted = self.planner.plan(info, &self.settings, projection, view, width as i32, height as i32, orthographic, &mut requests);
        for i in evicted {
            self.buffers[i] = None;
        }
        if let Some(loader) = &self.loader {
            loader.push(&requests);
        }
        self.ensure_target(width, height);
        // Uniforms.
        let light = Vec3::new(0.45, 0.6, 0.65).normalize();
        let axis = self.colour_axis.clamp(0, 3);
        let range = if axis > 0 { [info.lower[axis as usize - 1] as f32, info.upper[axis as usize - 1] as f32] } else { [0.0, 1.0] };
        let line = if light_theme { 0.25 } else { 0.85 };
        let frame = FrameUniform {
            view: view.to_cols_array(),
            projection: projection.to_cols_array(),
            light: [light.x, light.y, light.z, 0.38],
            params: [self.planner.pixels_per_unit, if orthographic { 1.0 } else { 0.0 }, self.settings.radius_scale, axis as f32],
            colour_range: [range[0], range[1], 0.0, 0.0],
            line_colour: [line, line, line, 1.0],
        };
        self.queue.write_buffer(&self.frame_uniform, 0, bytemuck::bytes_of(&frame));
        // The chunks to draw and their uniforms.
        let mut draws: Vec<(usize, u64)> = Vec::new();
        let mut chunk_data: Vec<u8> = Vec::new();
        for v in &self.planner.visible {
            let n = self.planner.draw_count(v);
            if n == 0 || self.buffers[v.chunk].is_none() {
                continue;
            }
            let chunk = &info.chunks[v.chunk];
            // A subsample is drawn with larger points so that the chunk still looks solid.
            let fraction = n as f64 / chunk.count as f64;
            let fill = (1.0 / fraction).cbrt().clamp(1.0, 2.5) as f32;
            let u = ChunkUniform {
                lower: [chunk.lower[0] as f32, chunk.lower[1] as f32, chunk.lower[2] as f32, fill],
                scale: [
                    ((chunk.upper[0] - chunk.lower[0]) / 65535.0) as f32,
                    ((chunk.upper[1] - chunk.lower[1]) / 65535.0) as f32,
                    ((chunk.upper[2] - chunk.lower[2]) / 65535.0) as f32,
                    0.0,
                ],
            };
            chunk_data.extend_from_slice(bytemuck::bytes_of(&u));
            chunk_data.resize(chunk_data.len().next_multiple_of(self.uniform_stride as usize), 0);
            draws.push((v.chunk, n));
        }
        if !draws.is_empty() {
            let needed = draws.len();
            if self.chunk_uniform.as_ref().is_none_or(|c| c.2 < needed) {
                let capacity = needed.next_power_of_two().max(64);
                let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("cloud chunks"),
                    size: capacity as u64 * self.uniform_stride,
                    usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
                    mapped_at_creation: false,
                });
                let group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("cloud chunks"),
                    layout: &self.chunk_layout,
                    entries: &[wgpu::BindGroupEntry {
                        binding: 0,
                        resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding { buffer: &buffer, offset: 0, size: wgpu::BufferSize::new(std::mem::size_of::<ChunkUniform>() as u64) }),
                    }],
                });
                self.chunk_uniform = Some((buffer, group, capacity));
            }
            self.queue.write_buffer(&self.chunk_uniform.as_ref().unwrap().0, 0, &chunk_data);
        }
        // Occlusion queries: each tested chunk's box against the atoms drawn now.
        let query = self.settings.occlusion_culling && self.readback.is_none() && self.queries.is_some() && !self.planner.tested.is_empty();
        let mut boxes: Vec<BoxInstance> = Vec::new();
        let mut box_buffer = None;
        if query {
            let r = self.planner.max_radius;
            for &i in &self.planner.tested {
                let c = &info.chunks[i];
                boxes.push(BoxInstance {
                    a: [c.lower[0] as f32 - r, c.lower[1] as f32 - r, c.lower[2] as f32 - r],
                    b: [c.upper[0] as f32 + r, c.upper[1] as f32 + r, c.upper[2] as f32 + r],
                });
            }
            let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("occlusion boxes"),
                size: std::mem::size_of_val(boxes.as_slice()) as u64,
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.queue.write_buffer(&buffer, 0, bytemuck::cast_slice(&boxes));
            box_buffer = Some(buffer);
        }
        let (_, _, colour, depth) = self.target.as_ref().unwrap();
        let colour_view = colour.create_view(&Default::default());
        let depth_view = depth.create_view(&Default::default());
        let mut drawn = 0u64;
        {
            let bg = background;
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("cloud"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &colour_view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color { r: bg[0] as f64, g: bg[1] as f64, b: bg[2] as f64, a: bg[3] as f64 }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &depth_view,
                    depth_ops: Some(wgpu::Operations { load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Store }),
                    stencil_ops: None,
                }),
                occlusion_query_set: if query { self.queries.as_ref().map(|q| &q.0) } else { None },
                ..Default::default()
            });
            let frame_group = self.frame_group.as_ref().unwrap();
            pass.set_bind_group(0, frame_group, &[]);
            if let Some((_, group, _)) = &self.chunk_uniform {
                pass.set_pipeline(&self.point_pipeline);
                for (k, &(chunk, n)) in draws.iter().enumerate() {
                    let buffer = self.buffers[chunk].as_ref().unwrap();
                    pass.set_bind_group(1, group, &[(k as u64 * self.uniform_stride) as u32]);
                    pass.set_vertex_buffer(0, buffer.slice(0..n * ATOM_BYTES));
                    pass.draw(0..4, 0..n as u32);
                    drawn += n;
                }
            }
            if let Some(boxes_buffer) = &box_buffer {
                pass.set_pipeline(&self.box_pipeline);
                pass.set_vertex_buffer(0, self.unit_box.slice(..));
                pass.set_vertex_buffer(1, boxes_buffer.slice(..));
                for k in 0..boxes.len() as u32 {
                    pass.begin_occlusion_query(k);
                    pass.draw(0..36, k..k + 1);
                    pass.end_occlusion_query();
                }
            }
            if self.settings.show_box {
                if let Some(lines) = &self.line_buffer {
                    pass.set_pipeline(&self.line_pipeline);
                    pass.set_vertex_buffer(0, lines.slice(..));
                    pass.draw(0..24, 0..1);
                }
            }
        }
        if query {
            let (set, resolve, read) = self.queries.as_ref().unwrap();
            let n = boxes.len() as u32;
            encoder.resolve_query_set(set, 0..n, resolve, 0);
            encoder.copy_buffer_to_buffer(resolve, 0, read, 0, n as u64 * 8);
            self.readback = Some(Readback { chunks: self.planner.tested.clone(), ready: Arc::new(AtomicBool::new(false)), mapping: false });
        }
        let occluded = self.planner.tested.iter().filter(|&&i| self.planner.occluded[i]).count();
        let resident: u64 = self.planner.chunks.iter().map(|c| c.resident).sum();
        self.stats.drawn_atoms = drawn;
        self.stats.resident_atoms = resident;
        self.stats.visible_chunks = self.planner.visible.len();
        self.stats.occluded_chunks = if self.settings.occlusion_culling { occluded } else { 0 };
        self.stats.loading_chunks = self.planner.in_flight;
        self.stats.resident_mb = self.planner.resident_bytes as f64 / (1024.0 * 1024.0);
        self.stats.complete = self.planner.complete && self.planner.in_flight == 0;
        Some(encoder)
    }

    /// Starts reading the occlusion results of the frame just submitted.
    fn after_submit(&mut self) {
        if let (Some(r), Some((_, _, read))) = (&mut self.readback, &self.queries) {
            if !r.mapping {
                r.mapping = true;
                let ready = r.ready.clone();
                read.map_async(wgpu::MapMode::Read, .., move |result| {
                    if result.is_ok() {
                        ready.store(true, Ordering::Release);
                    }
                });
            }
        }
    }

    /// Draws a frame into the target texture and leaves it on the GPU (the
    /// desktop shows it). None when no cloud is open.
    #[allow(clippy::too_many_arguments)]
    pub fn render_texture(&mut self, projection: &Mat4, view: &Mat4, width: u32, height: u32, orthographic: bool, light_theme: bool, background: [f32; 4]) -> Option<wgpu::TextureView> {
        let (width, height) = (width.max(1), height.max(1));
        let encoder = self.encode(projection, view, width, height, orthographic, light_theme, background)?;
        self.queue.submit([encoder.finish()]);
        self.after_submit();
        Some(self.target.as_ref().unwrap().2.create_view(&Default::default()))
    }

    /// Image export: draws until the atoms have streamed in (at most
    /// `timeout`) and returns RGBA rows, top first (captureVulkan).
    #[allow(clippy::too_many_arguments)]
    pub fn capture(&mut self, projection: &Mat4, view: &Mat4, width: u32, height: u32, orthographic: bool, light_theme: bool, background: [f32; 4], timeout: std::time::Duration) -> Result<Vec<u8>, String> {
        let start = std::time::Instant::now();
        let (width, height) = (width.max(1), height.max(1));
        let mut pass = 0;
        loop {
            let encoder = self.encode(projection, view, width, height, orthographic, light_theme, background).ok_or("No atom cloud open")?;
            self.queue.submit([encoder.finish()]);
            self.after_submit();
            let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
            pass += 1;
            if (pass > 2 && self.stats.complete) || start.elapsed() > timeout {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        let (_, _, colour, _) = self.target.as_ref().unwrap();
        let padded_row = (width * 4).div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT) * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("cloud readback"),
            size: (padded_row * height) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let mut encoder = self.device.create_command_encoder(&Default::default());
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo { texture: colour, mip_level: 0, origin: wgpu::Origin3d::ZERO, aspect: wgpu::TextureAspect::All },
            wgpu::TexelCopyBufferInfo { buffer: &readback, layout: wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(padded_row), rows_per_image: Some(height) } },
            wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
        );
        self.queue.submit([encoder.finish()]);
        let slice = readback.slice(..);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
        let mut rgba = Vec::with_capacity((width * height * 4) as usize);
        {
            let data = slice.get_mapped_range().map_err(|e| format!("{e:?}"))?;
            for row in 0..height {
                let start = (row * padded_row) as usize;
                rgba.extend_from_slice(&data[start..start + (width * 4) as usize]);
            }
        }
        readback.unmap();
        Ok(rgba)
    }
}

impl Drop for CloudRenderer {
    fn drop(&mut self) {
        self.close();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_cloud::{Chunk, Species};

    /// A 4 x 1 x 1 row of 10 Angstrom chunks along x, 100 000 atoms each.
    fn row() -> CloudInfo {
        let chunks = (0..4)
            .map(|i| Chunk { offset: 64 + i as u64 * 800_000, count: 100_000, lower: [i as f64 * 10.0, 0.0, 0.0], upper: [i as f64 * 10.0 + 10.0, 10.0, 10.0] })
            .collect();
        CloudInfo {
            atoms: 400_000,
            lower: [0.0; 3],
            upper: [40.0, 10.0, 10.0],
            species: vec![Species { name: "Cu".into(), color: [0.8, 0.5, 0.2], radius: 1.28 }],
            chunks,
            ..Default::default()
        }
    }

    /// glm::perspective / lookAt (clip depth -1..1).
    fn camera(eye: Vec3, target: Vec3) -> (Mat4, Mat4) {
        (crate::scene::perspective(45f32.to_radians(), 1.0, 0.5, 20000.0), crate::scene::look_at(eye, target, Vec3::Y))
    }

    #[test]
    fn culls_chunks_outside_the_view() {
        let info = row();
        let mut p = Planner::new(4);
        let mut requests = Vec::new();
        // Looking at the first chunk from close by, along -z: the far chunks are outside.
        let (proj, view) = camera(Vec3::new(5.0, 5.0, 25.0), Vec3::new(5.0, 5.0, 5.0));
        p.plan(&info, &CloudSettings::default(), &proj, &view, 800, 800, false, &mut requests);
        let chunks: Vec<usize> = p.visible.iter().map(|v| v.chunk).collect();
        assert!(chunks.contains(&0) && !chunks.contains(&3), "{chunks:?}");
        assert_eq!(p.tested.len(), p.visible.len());
        // Looking away: nothing.
        let (proj, view) = camera(Vec3::new(5.0, 5.0, 25.0), Vec3::new(5.0, 5.0, 50.0));
        p.plan(&info, &CloudSettings::default(), &proj, &view, 800, 800, false, &mut Vec::new());
        assert!(p.visible.is_empty());
    }

    #[test]
    fn level_of_detail_and_streaming() {
        let info = row();
        let settings = CloudSettings::default();
        let mut p = Planner::new(4);
        let mut requests = Vec::new();
        // Far away: subsampled to the cover, every chunk at least 256 atoms.
        let (proj, view) = camera(Vec3::new(20.0, 5.0, 5000.0), Vec3::new(20.0, 5.0, 5.0));
        p.plan(&info, &settings, &proj, &view, 400, 400, false, &mut requests);
        assert_eq!(p.visible.len(), 4);
        for v in &p.visible {
            assert!(v.want < 1000.0 && v.want >= 256.0, "{v:?}");
            assert_eq!(v.want, v.cover);
        }
        assert!(!p.complete);
        // Requests are rounded up to whole blocks and counted in flight.
        assert_eq!(requests.len(), 4);
        assert!(requests.iter().all(|r| r.first == 0 && r.count == BLOCK));
        assert_eq!(p.in_flight, 4);
        assert!(p.chunks.iter().all(|c| c.requested == BLOCK));
        // The loads arrive: the plan is complete.
        for r in requests.clone() {
            let capacity = p.accept_load(&info, &r, true).unwrap();
            assert_eq!(capacity, BLOCK);
            p.commit_load(&r, capacity);
            p.in_flight -= 1;
        }
        assert_eq!(p.resident_bytes, 4 * BLOCK * 8);
        requests.clear();
        p.plan(&info, &settings, &proj, &view, 400, 400, false, &mut requests);
        assert!(p.complete && requests.is_empty());
        assert_eq!(p.draw_count(&p.visible[0]), p.visible[0].want.ceil() as u64);
        // Close up: every atom of the nearest chunk.
        let (proj, view) = camera(Vec3::new(5.0, 5.0, 14.0), Vec3::new(5.0, 5.0, 5.0));
        requests.clear();
        p.plan(&info, &settings, &proj, &view, 1000, 1000, false, &mut requests);
        let nearest = p.visible[0];
        assert_eq!((nearest.chunk, nearest.want), (0, 100_000.0));
        // Prefixes after the resident one; a large request comes in pieces of
        // at most MAX_LOAD.
        let r = requests.iter().find(|r| r.chunk == 0).unwrap();
        assert_eq!((r.first, r.count), (BLOCK, 100_000 - BLOCK));
    }

    #[test]
    fn budgets_scale_the_covers() {
        let info = row();
        let settings = CloudSettings { point_budget: 1000.0, level_of_detail: false, ..Default::default() };
        let mut p = Planner::new(4);
        let (proj, view) = camera(Vec3::new(20.0, 5.0, 200.0), Vec3::new(20.0, 5.0, 5.0));
        p.plan(&info, &settings, &proj, &view, 400, 400, false, &mut Vec::new());
        // Without level of detail every atom is wanted; covers exceed the budget,
        // so each is scaled down to the 256-atom minimum.
        assert_eq!(p.visible.len(), 4);
        assert!(p.visible.iter().all(|v| v.want == 256.0), "{:?}", p.visible);
        // A budget between the covers and the total: the nearest get more.
        let settings = CloudSettings { point_budget: 250_000.0, density: 0.01, ..Default::default() };
        let (proj, view) = camera(Vec3::new(-30.0, 5.0, 5.0), Vec3::new(20.0, 5.0, 5.0));
        p.plan(&info, &settings, &proj, &view, 400, 400, false, &mut Vec::new());
        let total: f64 = p.visible.iter().map(|v| v.want).sum();
        assert!(total <= 250_000.0 + 1.0, "{total}");
        let nearest = p.visible.iter().min_by(|a, b| a.distance.total_cmp(&b.distance)).unwrap();
        assert_eq!(nearest.chunk, 0);
        assert_eq!(nearest.want, 100_000.0);
    }

    #[test]
    fn eviction_and_out_of_order_loads() {
        let info = row();
        let mut p = Planner::new(4);
        for (i, used) in [(0usize, 5u64), (1, 2), (2, 9), (3, 1)] {
            p.chunks[i] = ChunkState { capacity: 100_000, resident: 100_000, requested: 100_000, wanted: 0, last_used: used };
        }
        p.resident_bytes = 4 * 100_000 * 8;
        // 2 MB budget: release the least recently drawn until below 90 % of it.
        let settings = CloudSettings { memory_budget_mb: 2.0, ..Default::default() };
        let released = p.evict(&settings, 10);
        assert_eq!(released, vec![3, 1]);
        assert_eq!(p.resident_bytes, 2 * 100_000 * 8);
        assert_eq!(p.chunks[3], ChunkState::default());
        // A load that no longer starts at the resident prefix is forgotten.
        p.chunks[0].requested = 120_000;
        assert_eq!(p.accept_load(&info, &LoadRequest { chunk: 0, first: 110_000, count: 10 }, true), None);
        assert_eq!(p.chunks[0].requested, 100_000);
        // Failed reads too.
        p.chunks[1].requested = BLOCK;
        assert_eq!(p.accept_load(&info, &LoadRequest { chunk: 1, first: 0, count: BLOCK }, false), None);
        assert_eq!(p.chunks[1].requested, 0);
        // Capacity doubles (up to the chunk's atoms).
        let r = LoadRequest { chunk: 1, first: 0, count: 1000 };
        assert_eq!(p.accept_load(&info, &r, true), Some(1000));
        p.commit_load(&r, 1000);
        let r = LoadRequest { chunk: 1, first: 1000, count: 10 };
        assert_eq!(p.accept_load(&info, &r, true), Some(2000));
    }

    #[test]
    fn outline_of_box_and_cell() {
        let mut info = row();
        let lines = outline(&info);
        assert_eq!(lines.len(), 24);
        assert_eq!(lines[0], [0.0, 0.0, 0.0]);
        assert_eq!(lines[1], [40.0, 0.0, 0.0]);
        info.has_cell = true;
        info.cell = [[2.0, 0.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        info.cell_origin = [1.0, 1.0, 1.0];
        let lines = outline(&info);
        assert_eq!(lines[23], [4.0, 3.0, 4.0]);
        assert_eq!(cloud_centre(&row()), Vec3::new(20.0, 5.0, 5.0));
        assert!((cloud_radius(&row()) - 0.5 * (1600f32 + 200.0).sqrt()).abs() < 1e-4);
    }

    /// Writes a small cloud and renders it on the default adapter (skipped without one).
    #[test]
    fn renders_a_generated_cloud() {
        let dir = std::env::temp_dir().join(format!("af-render-cloud-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cu.afcloud");
        af_cloud::generate_crystal(&path, "fcc", 3.615, 0.0, [8, 8, 8], "Cu", "", None, af_cloud::BuildOptions::default()).unwrap();
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let Ok(adapter) = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions::default())) else {
            eprintln!("no GPU adapter; skipped");
            return;
        };
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default())).unwrap();
        let mut r = CloudRenderer::with_device(device, queue);
        r.open(&path).unwrap();
        assert_eq!(r.stats().total_atoms, 2048);
        let centre = r.centre();
        let (proj, view) = camera(centre + Vec3::new(0.0, 0.0, 3.0 * r.radius()), centre);
        let rgba = r.capture(&proj, &view, 128, 128, false, false, [0.0, 0.0, 0.0, 1.0], std::time::Duration::from_secs(20)).unwrap();
        assert!(r.stats().complete, "{:?}", r.stats());
        assert_eq!(r.stats().drawn_atoms, 2048);
        // Copper atoms fill the middle; the corner stays background.
        let px = |x: usize, y: usize| &rgba[(y * 128 + x) * 4..(y * 128 + x) * 4 + 3];
        assert!(px(64, 64)[0] > 60, "{:?}", px(64, 64));
        assert_eq!(px(1, 1), &[0, 0, 0]);
        drop(r);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
