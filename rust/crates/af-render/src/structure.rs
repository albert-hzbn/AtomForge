//! The structure renderer (the C++ graphics/Renderer with app/ImageExport's
//! offscreen capture): atoms as instanced 40x40 sphere meshes, bonds as
//! 32-sided cylinders, a 2048x2048 PCF shadow map, the cell box as wide
//! lines, and the display settings' lighting, depth cueing and outlines.

use bytemuck::{Pod, Zeroable};
use glam::Mat4;

use af_core::display::DisplaySettings;

use crate::scene::{cylinder_mesh, sphere_mesh, FrameView, Scene};

const COLOUR_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;
const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
const SHADOW_SIZE: u32 = 2048;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct FrameUniform {
    projection: [f32; 16],
    view: [f32; 16],
    light_mvp: [f32; 16],
    light_pos: [f32; 4],
    view_pos: [f32; 4],
    material: [f32; 4],
    material2: [f32; 4],
    fog: [f32; 4],
    fog_color: [f32; 4],
    outline: [f32; 4],
    outline_color: [f32; 4],
    line: [f32; 4],
    viewport: [f32; 4],
    selection: [f32; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct HaloInstance {
    centre: [f32; 3],
    radius: f32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct AtomInstance {
    centre: [f32; 3],
    scale: f32,
    colour: [f32; 3],
    shininess: f32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct BondInstance {
    start: [f32; 3],
    radius: f32,
    end: [f32; 3],
    shininess_a: f32,
    colour_a: [f32; 3],
    shininess_b: f32,
    colour_b: [f32; 3],
    pad: f32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct LineInstance {
    a: [f32; 3],
    b: [f32; 3],
}

/// The lighting and finish of a frame (Renderer's members set by applyDisplaySettings).
#[derive(Clone, Copy, Debug)]
pub struct Look {
    pub ambient: f32,
    pub saturation: f32,
    pub contrast: f32,
    pub shadow_strength: f32,
    pub specular_intensity: f32,
    pub shininess_scale: f32,
    pub shininess_floor: f32,
    pub fog_strength: f32,
    pub fog_near: f32,
    pub fog_far: f32,
    pub fog_colour: [f32; 3],
    pub outline_width: f32,
    pub outline_colour: [f32; 3],
    pub box_line_width: f32,
}

impl Look {
    /// applyDisplaySettings: the display's lighting with the frame's depth range.
    pub fn new(d: &DisplaySettings, frame: &FrameView, background: [f32; 3]) -> Self {
        Self {
            ambient: d.ambient,
            saturation: d.saturation,
            contrast: d.contrast,
            shadow_strength: d.shadow_strength,
            specular_intensity: d.specular_intensity,
            shininess_scale: d.shininess_scale,
            shininess_floor: d.shininess_floor,
            fog_strength: d.depth_cueing,
            fog_near: frame.depth_near,
            fog_far: frame.depth_far,
            fog_colour: background,
            outline_width: d.outline_width,
            outline_colour: d.outline_colour,
            box_line_width: d.cell_line_width,
        }
    }
}

/// What to draw.
#[derive(Clone, Copy, Debug)]
pub struct DrawOptions {
    pub background: [f32; 4],
    pub show_bonds: bool,
    pub show_atoms: bool,
    pub show_box: bool,
    /// The cell lines' colour (None: dark on light backgrounds, light on dark).
    pub box_colour: Option<[f32; 3]>,
    /// The halo of selected atoms (DisplaySettings::selection_colour).
    pub selection_colour: [f32; 3],
}

struct Buffers {
    atoms: Option<(wgpu::Buffer, u32)>,
    bonds: Option<(wgpu::Buffer, u32)>,
    lines: Option<(wgpu::Buffer, u32)>,
    halos: Option<(wgpu::Buffer, u32)>,
}

/// Renders structure scenes offscreen and reads them back as RGBA.
pub struct StructureRenderer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    uniform: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    shadow_bind_group: wgpu::BindGroup,
    shadow_view: wgpu::TextureView,
    atom_pipeline: wgpu::RenderPipeline,
    bond_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,
    atom_shadow_pipeline: wgpu::RenderPipeline,
    bond_shadow_pipeline: wgpu::RenderPipeline,
    halo_pipeline: wgpu::RenderPipeline,
    halo_hidden_pipeline: wgpu::RenderPipeline,
    sphere_vertices: wgpu::Buffer,
    sphere_indices: wgpu::Buffer,
    sphere_index_count: u32,
    cylinder_vertices: wgpu::Buffer,
    cylinder_vertex_count: u32,
    buffers: Buffers,
    /// (width, height, colour, depth) of the last frame.
    target: Option<(u32, u32, wgpu::Texture, wgpu::Texture)>,
    pub adapter_name: String,
}

fn buffer_with(device: &wgpu::Device, queue: &wgpu::Queue, label: &str, data: &[u8], usage: wgpu::BufferUsages) -> wgpu::Buffer {
    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some(label),
        size: (data.len() as u64).max(16).next_multiple_of(4),
        usage: usage | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&buffer, 0, data);
    buffer
}

impl StructureRenderer {
    /// A renderer on the default (high-performance) adapter.
    pub fn new() -> Result<Self, String> {
        pollster::block_on(async {
            let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
            let adapter = instance
                .request_adapter(&wgpu::RequestAdapterOptions { power_preference: wgpu::PowerPreference::HighPerformance, ..Default::default() })
                .await
                .map_err(|e| format!("no graphics adapter: {e}"))?;
            let name = adapter.get_info().name;
            // Depth clamping in the colour pass, as the C++ GL_DEPTH_CLAMP (when the adapter has it).
            let features = adapter.features() & wgpu::Features::DEPTH_CLIP_CONTROL;
            let (device, queue) = adapter
                .request_device(&wgpu::DeviceDescriptor { required_features: features, ..Default::default() })
                .await
                .map_err(|e| format!("no graphics device: {e}"))?;
            Ok(Self::with_device(device, queue, name))
        })
    }

    pub fn with_device(device: wgpu::Device, queue: wgpu::Queue, adapter_name: String) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("structure"),
            source: wgpu::ShaderSource::Wgsl(include_str!("structure.wgsl").into()),
        });
        let uniform = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("frame"),
            size: std::mem::size_of::<FrameUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let shadow = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("shadow map"),
            size: wgpu::Extent3d { width: SHADOW_SIZE, height: SHADOW_SIZE, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let shadow_view = shadow.create_view(&Default::default());
        // A dummy depth texture bound while the shadow map itself is drawn.
        let dummy = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("no shadow"),
            size: wgpu::Extent3d { width: 1, height: 1, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH_FORMAT,
            usage: wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("shadow"),
            address_mode_u: wgpu::AddressMode::ClampToEdge,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            compare: Some(wgpu::CompareFunction::LessEqual),
            ..Default::default()
        });
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("frame"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                    ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Depth,
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Comparison),
                    count: None,
                },
            ],
        });
        let group = |view: &wgpu::TextureView| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("frame"),
                layout: &layout,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: uniform.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: wgpu::BindingResource::TextureView(view) },
                    wgpu::BindGroupEntry { binding: 2, resource: wgpu::BindingResource::Sampler(&sampler) },
                ],
            })
        };
        let bind_group = group(&shadow_view);
        let shadow_bind_group = group(&dummy.create_view(&Default::default()));
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("structure"),
            bind_group_layouts: &[Some(&layout)],
            ..Default::default()
        });
        let mesh_layout = wgpu::VertexBufferLayout {
            array_stride: 12,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &wgpu::vertex_attr_array![0 => Float32x3],
        };
        let atom_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<AtomInstance>() as u64,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &wgpu::vertex_attr_array![1 => Float32x3, 2 => Float32, 3 => Float32x3, 4 => Float32],
        };
        let bond_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<BondInstance>() as u64,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &wgpu::vertex_attr_array![1 => Float32x3, 2 => Float32, 3 => Float32x3, 4 => Float32, 5 => Float32x3, 6 => Float32, 7 => Float32x3],
        };
        let line_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<LineInstance>() as u64,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3],
        };
        let clamp = device.features().contains(wgpu::Features::DEPTH_CLIP_CONTROL);
        let colour_target = [Some(wgpu::ColorTargetState { format: COLOUR_FORMAT, blend: None, write_mask: wgpu::ColorWrites::ALL })];
        let make = |label: &str,
                    vs: &str,
                    fs: Option<&str>,
                    buffers: &[Option<wgpu::VertexBufferLayout>],
                    shadow_pass: bool| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some(label),
                layout: Some(&pipeline_layout),
                vertex: wgpu::VertexState { module: &shader, entry_point: Some(vs), buffers, compilation_options: Default::default() },
                fragment: fs.map(|fs| wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some(fs),
                    targets: &colour_target,
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleList,
                    front_face: wgpu::FrontFace::Ccw,
                    // Only back faces go into the shadow map (no self-shadow acne).
                    cull_mode: if shadow_pass { Some(wgpu::Face::Front) } else { None },
                    unclipped_depth: !shadow_pass && clamp,
                    ..Default::default()
                },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: DEPTH_FORMAT,
                    depth_write_enabled: Some(true),
                    depth_compare: Some(wgpu::CompareFunction::Less),
                    stencil: Default::default(),
                    // glPolygonOffset(3, 6) in the shadow pass.
                    bias: if shadow_pass { wgpu::DepthBiasState { constant: 6, slope_scale: 3.0, clamp: 0.0 } } else { Default::default() },
                }),
                multisample: Default::default(),
                multiview_mask: None,
                cache: None,
            })
        };
        let atom_pipeline = make("atoms", "atom_vs", Some("atom_fs"), &[Some(mesh_layout.clone()), Some(atom_layout.clone())], false);
        let bond_pipeline = make("bonds", "bond_vs", Some("bond_fs"), &[Some(mesh_layout.clone()), Some(bond_layout.clone())], false);
        let line_pipeline = make("lines", "line_vs", Some("line_fs"), &[Some(line_layout)], false);
        let atom_shadow_pipeline = make("atom shadows", "atom_shadow_vs", None, &[Some(mesh_layout.clone()), Some(atom_layout)], true);
        let bond_shadow_pipeline = make("bond shadows", "bond_shadow_vs", None, &[Some(mesh_layout), Some(bond_layout)], true);
        // The selection halo: blended, never writing depth; hidden parts first, faintly.
        let halo_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<HaloInstance>() as u64,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32],
        };
        let halo_target = [Some(wgpu::ColorTargetState { format: COLOUR_FORMAT, blend: Some(wgpu::BlendState::ALPHA_BLENDING), write_mask: wgpu::ColorWrites::ALL })];
        let halo = |label: &str, fs: &str, compare: wgpu::CompareFunction| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some(label),
                layout: Some(&pipeline_layout),
                vertex: wgpu::VertexState { module: &shader, entry_point: Some("halo_vs"), buffers: &[Some(halo_layout.clone())], compilation_options: Default::default() },
                fragment: Some(wgpu::FragmentState { module: &shader, entry_point: Some(fs), targets: &halo_target, compilation_options: Default::default() }),
                primitive: wgpu::PrimitiveState { topology: wgpu::PrimitiveTopology::TriangleList, cull_mode: None, unclipped_depth: clamp, ..Default::default() },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: DEPTH_FORMAT,
                    depth_write_enabled: Some(false),
                    depth_compare: Some(compare),
                    stencil: Default::default(),
                    bias: Default::default(),
                }),
                multisample: Default::default(),
                multiview_mask: None,
                cache: None,
            })
        };
        let halo_hidden_pipeline = halo("selection halo (hidden)", "halo_hidden_fs", wgpu::CompareFunction::Greater);
        let halo_pipeline = halo("selection halo", "halo_fs", wgpu::CompareFunction::LessEqual);
        let (sphere, indices) = sphere_mesh(40, 40);
        let cylinder = cylinder_mesh(32);
        let sphere_vertices = buffer_with(&device, &queue, "sphere", bytemuck::cast_slice(&sphere), wgpu::BufferUsages::VERTEX);
        let sphere_indices = buffer_with(&device, &queue, "sphere indices", bytemuck::cast_slice(&indices), wgpu::BufferUsages::INDEX);
        let cylinder_vertices = buffer_with(&device, &queue, "cylinder", bytemuck::cast_slice(&cylinder), wgpu::BufferUsages::VERTEX);
        Self {
            device,
            queue,
            uniform,
            bind_group,
            shadow_bind_group,
            shadow_view,
            atom_pipeline,
            bond_pipeline,
            line_pipeline,
            atom_shadow_pipeline,
            bond_shadow_pipeline,
            halo_pipeline,
            halo_hidden_pipeline,
            sphere_vertices,
            sphere_indices,
            sphere_index_count: indices.len() as u32,
            cylinder_vertices,
            cylinder_vertex_count: cylinder.len() as u32,
            buffers: Buffers { atoms: None, bonds: None, lines: None, halos: None },
            target: None,
            adapter_name,
        }
    }

    /// Uploads a scene's atoms, bonds and box (SceneBuffers::upload).
    pub fn upload(&mut self, scene: &Scene) {
        let d = &scene.instances;
        let atoms: Vec<AtomInstance> = (0..d.positions.len())
            .map(|i| AtomInstance { centre: d.positions[i].to_array(), scale: d.scales[i], colour: d.colors[i], shininess: d.shininess[i] })
            .collect();
        let b = &scene.bonds;
        let bonds: Vec<BondInstance> = (0..b.starts.len())
            .map(|i| BondInstance {
                start: b.starts[i].to_array(),
                radius: b.radii[i],
                end: b.ends[i].to_array(),
                shininess_a: b.shininess_a[i],
                colour_a: b.colors_a[i],
                shininess_b: b.shininess_b[i],
                colour_b: b.colors_b[i],
                pad: 0.0,
            })
            .collect();
        let lines: Vec<LineInstance> =
            d.box_lines.chunks_exact(2).map(|p| LineInstance { a: p[0].to_array(), b: p[1].to_array() }).collect();
        let halos: Vec<HaloInstance> =
            scene.selected.iter().filter(|&&i| i < d.positions.len()).map(|&i| HaloInstance { centre: d.positions[i].to_array(), radius: d.scales[i] }).collect();
        let make = |label: &str, bytes: &[u8], count: usize| {
            (count > 0).then(|| (buffer_with(&self.device, &self.queue, label, bytes, wgpu::BufferUsages::VERTEX), count as u32))
        };
        self.buffers = Buffers {
            atoms: make("atoms", bytemuck::cast_slice(&atoms), atoms.len()),
            bonds: make("bonds", bytemuck::cast_slice(&bonds), bonds.len()),
            lines: make("lines", bytemuck::cast_slice(&lines), lines.len()),
            halos: make("selection halos", bytemuck::cast_slice(&halos), halos.len()),
        };
    }

    fn write_uniform(&self, frame: &FrameView, look: &Look, line_colour: [f32; 3], selection: [f32; 3]) {
        let v4 = |v: glam::Vec3| [v.x, v.y, v.z, 0.0];
        let c4 = |c: [f32; 3]| [c[0], c[1], c[2], 0.0];
        let m = |m: &Mat4| m.to_cols_array();
        let u = FrameUniform {
            projection: m(&frame.projection),
            view: m(&frame.view),
            light_mvp: m(&frame.light_mvp),
            light_pos: v4(frame.light_position),
            view_pos: v4(frame.camera_position),
            material: [look.ambient, look.saturation, look.contrast, look.specular_intensity],
            material2: [look.shadow_strength, look.shininess_scale, look.shininess_floor, 0.0],
            fog: [look.fog_near, look.fog_far, look.fog_strength, 0.0],
            fog_color: c4(look.fog_colour),
            outline: [look.outline_width, 0.0, 0.0, 0.0],
            outline_color: c4(look.outline_colour),
            line: [line_colour[0], line_colour[1], line_colour[2], look.box_line_width],
            viewport: [frame.width as f32, frame.height as f32, 0.0, 0.0],
            selection: c4(selection),
        };
        self.queue.write_buffer(&self.uniform, 0, bytemuck::bytes_of(&u));
    }

    /// The colour and depth targets of this size (kept between frames).
    fn ensure_target(&mut self, width: u32, height: u32) {
        if self.target.as_ref().is_some_and(|t| t.0 == width && t.1 == height) {
            return;
        }
        let extent = wgpu::Extent3d { width, height, depth_or_array_layers: 1 };
        let texture = |format, usage, view_formats: &[wgpu::TextureFormat]| {
            self.device.create_texture(&wgpu::TextureDescriptor {
                label: None,
                size: extent,
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage,
                view_formats,
            })
        };
        let colour = texture(
            COLOUR_FORMAT,
            wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC | wgpu::TextureUsages::TEXTURE_BINDING,
            &[],
        );
        let depth = texture(DEPTH_FORMAT, wgpu::TextureUsages::RENDER_ATTACHMENT, &[]);
        self.target = Some((width, height, colour, depth));
    }

    /// Renders the uploaded scene into the target texture and leaves it on the
    /// GPU (the desktop shows it directly).
    pub fn render_texture(&mut self, frame: &FrameView, look: &Look, options: &DrawOptions) -> Result<wgpu::TextureView, String> {
        let encoder = self.encode(frame, look, options)?;
        self.queue.submit([encoder.finish()]);
        Ok(self.target.as_ref().expect("target").2.create_view(&Default::default()))
    }

    /// Renders the uploaded scene and reads it back as RGBA rows, top row first
    /// (captureSceneToRgba).
    pub fn render(&mut self, frame: &FrameView, look: &Look, options: &DrawOptions) -> Result<Vec<u8>, String> {
        let mut encoder = self.encode(frame, look, options)?;
        let (width, height) = (frame.width as u32, frame.height as u32);
        let extent = wgpu::Extent3d { width, height, depth_or_array_layers: 1 };
        let colour = &self.target.as_ref().expect("target").2;
        let padded_row = (width * 4).div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT) * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: (padded_row * height) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo { texture: colour, mip_level: 0, origin: wgpu::Origin3d::ZERO, aspect: wgpu::TextureAspect::All },
            wgpu::TexelCopyBufferInfo {
                buffer: &readback,
                layout: wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(padded_row), rows_per_image: Some(height) },
            },
            extent,
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

    fn encode(&mut self, frame: &FrameView, look: &Look, options: &DrawOptions) -> Result<wgpu::CommandEncoder, String> {
        let (width, height) = (frame.width, frame.height);
        if width <= 0 || height <= 0 {
            return Err("Invalid viewport size for image export.".into());
        }
        let (width, height) = (width as u32, height as u32);
        let bg = options.background;
        let line_colour = options.box_colour.unwrap_or(if bg[0] > 0.5 { [0.25; 3] } else { [0.85; 3] });
        self.write_uniform(frame, look, line_colour, options.selection_colour);
        self.ensure_target(width, height);
        let (_, _, colour, depth) = self.target.as_ref().expect("target");
        let mut encoder = self.device.create_command_encoder(&Default::default());
        let draw_bonds = options.show_bonds && options.show_atoms;
        // Shadow passes first.
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("shadow"),
                color_attachments: &[],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.shadow_view,
                    depth_ops: Some(wgpu::Operations { load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Store }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_bind_group(0, &self.shadow_bind_group, &[]);
            if options.show_atoms {
                if let Some((atoms, count)) = &self.buffers.atoms {
                    pass.set_pipeline(&self.atom_shadow_pipeline);
                    pass.set_vertex_buffer(0, self.sphere_vertices.slice(..));
                    pass.set_vertex_buffer(1, atoms.slice(..));
                    pass.set_index_buffer(self.sphere_indices.slice(..), wgpu::IndexFormat::Uint32);
                    pass.draw_indexed(0..self.sphere_index_count, 0, 0..*count);
                }
            }
            if draw_bonds {
                if let Some((bonds, count)) = &self.buffers.bonds {
                    pass.set_pipeline(&self.bond_shadow_pipeline);
                    pass.set_vertex_buffer(0, self.cylinder_vertices.slice(..));
                    pass.set_vertex_buffer(1, bonds.slice(..));
                    pass.draw(0..self.cylinder_vertex_count, 0..*count);
                }
            }
        }
        {
            let colour_view = colour.create_view(&Default::default());
            let depth_view = depth.create_view(&Default::default());
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("scene"),
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
                ..Default::default()
            });
            pass.set_bind_group(0, &self.bind_group, &[]);
            if draw_bonds {
                if let Some((bonds, count)) = &self.buffers.bonds {
                    pass.set_pipeline(&self.bond_pipeline);
                    pass.set_vertex_buffer(0, self.cylinder_vertices.slice(..));
                    pass.set_vertex_buffer(1, bonds.slice(..));
                    pass.draw(0..self.cylinder_vertex_count, 0..*count);
                }
            }
            if options.show_atoms {
                if let Some((atoms, count)) = &self.buffers.atoms {
                    pass.set_pipeline(&self.atom_pipeline);
                    pass.set_vertex_buffer(0, self.sphere_vertices.slice(..));
                    pass.set_vertex_buffer(1, atoms.slice(..));
                    pass.set_index_buffer(self.sphere_indices.slice(..), wgpu::IndexFormat::Uint32);
                    pass.draw_indexed(0..self.sphere_index_count, 0, 0..*count);
                }
            }
            if options.show_box {
                if let Some((lines, count)) = &self.buffers.lines {
                    pass.set_pipeline(&self.line_pipeline);
                    pass.set_vertex_buffer(0, lines.slice(..));
                    pass.draw(0..6, 0..*count);
                }
            }
            if options.show_atoms {
                if let Some((halos, count)) = &self.buffers.halos {
                    pass.set_vertex_buffer(0, halos.slice(..));
                    pass.set_pipeline(&self.halo_hidden_pipeline);
                    pass.draw(0..6, 0..*count);
                    pass.set_pipeline(&self.halo_pipeline);
                    pass.draw(0..6, 0..*count);
                }
            }
        }
        Ok(encoder)
    }
}
