//! Offscreen rendering of atomistic scenes with wgpu. The desktop app shows
//! the rendered frames in its egui window; the command line saves them.

pub mod cloud;
pub mod colouring;
pub mod grab;
pub mod measure;
pub mod polyhedra;
pub mod region;
pub mod scene;
pub mod structure;

use bytemuck::{Pod, Zeroable};
use glam::{Mat4, Vec3, Vec4};

/// One atom drawn as a sphere.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct SphereInstance {
    pub centre: [f32; 3],
    pub radius: f32,
    pub colour: [f32; 4],
}

/// Camera looking at `target` from `eye`.
#[derive(Clone, Copy, Debug)]
pub struct Camera {
    pub eye: Vec3,
    pub target: Vec3,
    pub up: Vec3,
    pub perspective: bool,
    pub fov_y_degrees: f32,
    /// Half the visible height (orthographic).
    pub ortho_half_height: f32,
    pub near: f32,
    pub far: f32,
}

impl Camera {
    pub fn view(&self) -> Mat4 {
        Mat4::look_at_rh(self.eye, self.target, self.up)
    }

    pub fn projection(&self, aspect: f32) -> Mat4 {
        if self.perspective {
            Mat4::perspective_rh(self.fov_y_degrees.to_radians(), aspect, self.near, self.far)
        } else {
            let h = self.ortho_half_height;
            Mat4::orthographic_rh(-h * aspect, h * aspect, -h, h, self.near, self.far)
        }
    }
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CameraUniform {
    view: [[f32; 4]; 4],
    projection: [[f32; 4]; 4],
    light: [f32; 4],
    params: [f32; 4],
}

/// A frame in BGRA8 (premultiplied, opaque), rows tightly packed.
pub struct Frame {
    pub width: u32,
    pub height: u32,
    pub bgra: Vec<u8>,
}

/// Headless renderer: a device, the sphere pipeline and a target texture.
pub struct OffscreenRenderer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::RenderPipeline,
    bind_group: wgpu::BindGroup,
    uniform: wgpu::Buffer,
    instances: Option<(wgpu::Buffer, u32)>,
    target: Option<Target>,
    pub background: [f64; 3],
    pub adapter_name: String,
}

struct Target {
    width: u32,
    height: u32,
    colour: wgpu::Texture,
    depth: wgpu::Texture,
    readback: wgpu::Buffer,
    padded_row: u32,
}

const COLOUR_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Bgra8Unorm;
const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

impl OffscreenRenderer {
    /// A renderer on the default (high-performance) adapter.
    pub fn new() -> Result<Self, String> {
        pollster::block_on(Self::create())
    }

    async fn create() -> Result<Self, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                ..Default::default()
            })
            .await
            .map_err(|e| format!("no graphics adapter: {e}"))?;
        let adapter_name = adapter.get_info().name;
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default())
            .await
            .map_err(|e| format!("no graphics device: {e}"))?;
        Ok(Self::with_device(device, queue, adapter_name))
    }

    /// A renderer on a device owned elsewhere (a UI toolkit's wgpu device),
    /// so frames can stay on the GPU.
    pub fn with_device(device: wgpu::Device, queue: wgpu::Queue, adapter_name: String) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("spheres"),
            source: wgpu::ShaderSource::Wgsl(include_str!("spheres.wgsl").into()),
        });
        let uniform = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("camera"),
            size: std::mem::size_of::<CameraUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("camera"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("camera"),
            layout: &layout,
            entries: &[wgpu::BindGroupEntry { binding: 0, resource: uniform.as_entire_binding() }],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("spheres"),
            bind_group_layouts: &[Some(&layout)],
            ..Default::default()
        });
        let attributes = wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32, 2 => Float32x4];
        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("spheres"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                buffers: &[Some(wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<SphereInstance>() as u64,
                    step_mode: wgpu::VertexStepMode::Instance,
                    attributes: &attributes,
                })],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: COLOUR_FORMAT,
                    blend: None,
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState { topology: wgpu::PrimitiveTopology::TriangleStrip, ..Default::default() },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: DEPTH_FORMAT,
                depth_write_enabled: Some(true),
                depth_compare: Some(wgpu::CompareFunction::Less),
                stencil: Default::default(),
                bias: Default::default(),
            }),
            multisample: Default::default(),
            multiview_mask: None,
            cache: None,
        });
        Self {
            device,
            queue,
            pipeline,
            bind_group,
            uniform,
            instances: None,
            target: None,
            background: [0.08, 0.09, 0.12],
            adapter_name,
        }
    }

    /// Uploads the spheres (kept until the next call).
    pub fn set_spheres(&mut self, spheres: &[SphereInstance]) {
        if spheres.is_empty() {
            self.instances = None;
            return;
        }
        let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("spheres"),
            size: std::mem::size_of_val(spheres) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.queue.write_buffer(&buffer, 0, bytemuck::cast_slice(spheres));
        self.instances = Some((buffer, spheres.len() as u32));
    }

    fn ensure_target(&mut self, width: u32, height: u32) {
        if let Some(t) = &self.target {
            if t.width == width && t.height == height {
                return;
            }
        }
        let extent = wgpu::Extent3d { width, height, depth_or_array_layers: 1 };
        let texture = |format, usage, label| {
            self.device.create_texture(&wgpu::TextureDescriptor {
                label: Some(label),
                size: extent,
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage,
                view_formats: &[],
            })
        };
        let colour = texture(
            COLOUR_FORMAT,
            wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC | wgpu::TextureUsages::TEXTURE_BINDING,
            "colour",
        );
        let depth = texture(DEPTH_FORMAT, wgpu::TextureUsages::RENDER_ATTACHMENT, "depth");
        let padded_row = (width * 4).div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT) * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: (padded_row * height) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        self.target = Some(Target { width, height, colour, depth, readback, padded_row });
    }

    /// Renders into the target texture and leaves it on the GPU.
    pub fn render_texture(&mut self, camera: &Camera, width: u32, height: u32) -> wgpu::TextureView {
        let width = width.max(1);
        let height = height.max(1);
        self.ensure_target(width, height);
        let encoder = self.encode(camera, width, height);
        self.queue.submit([encoder.finish()]);
        self.target.as_ref().unwrap().colour.create_view(&Default::default())
    }

    /// Renders the spheres and reads the frame back.
    pub fn render(&mut self, camera: &Camera, width: u32, height: u32) -> Frame {
        let width = width.max(1);
        let height = height.max(1);
        self.ensure_target(width, height);
        let mut encoder = self.encode(camera, width, height);
        let target = self.target.as_ref().unwrap();
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &target.colour,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &target.readback,
                layout: wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(target.padded_row), rows_per_image: Some(height) },
            },
            wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
        );
        self.queue.submit([encoder.finish()]);
        let slice = target.readback.slice(..);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
        let mut bgra = Vec::with_capacity((width * height * 4) as usize);
        {
            let data = slice.get_mapped_range().expect("mapped readback");
            for row in 0..height {
                let start = (row * target.padded_row) as usize;
                bgra.extend_from_slice(&data[start..start + (width * 4) as usize]);
            }
        }
        target.readback.unmap();
        Frame { width, height, bgra }
    }

    fn encode(&mut self, camera: &Camera, width: u32, height: u32) -> wgpu::CommandEncoder {
        let view = camera.view();
        let light = view.transform_vector3(Vec3::new(0.4, 0.6, 1.0).normalize());
        let uniform = CameraUniform {
            view: view.to_cols_array_2d(),
            projection: camera.projection(width as f32 / height as f32).to_cols_array_2d(),
            light: Vec4::from((light, 0.0)).to_array(),
            params: [if camera.perspective { 1.0 } else { 0.0 }, 0.0, 0.0, 0.0],
        };
        self.queue.write_buffer(&self.uniform, 0, bytemuck::bytes_of(&uniform));
        let target = self.target.as_ref().unwrap();
        let colour_view = target.colour.create_view(&Default::default());
        let depth_view = target.depth.create_view(&Default::default());
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let [r, g, b] = self.background;
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("scene"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &colour_view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations { load: wgpu::LoadOp::Clear(wgpu::Color { r, g, b, a: 1.0 }), store: wgpu::StoreOp::Store },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &depth_view,
                    depth_ops: Some(wgpu::Operations { load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Store }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            if let Some((buffer, count)) = &self.instances {
                pass.set_pipeline(&self.pipeline);
                pass.set_bind_group(0, &self.bind_group, &[]);
                pass.set_vertex_buffer(0, buffer.slice(..));
                pass.draw(0..4, 0..*count);
            }
        }
        encoder
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_a_sphere_in_the_middle() {
        let Ok(mut r) = OffscreenRenderer::new() else {
            eprintln!("no GPU adapter; skipped");
            return;
        };
        r.background = [0.0, 0.0, 0.0];
        r.set_spheres(&[SphereInstance { centre: [0.0, 0.0, 0.0], radius: 1.0, colour: [1.0, 0.0, 0.0, 1.0] }]);
        let camera = Camera {
            eye: Vec3::new(0.0, 0.0, 10.0),
            target: Vec3::ZERO,
            up: Vec3::Y,
            perspective: false,
            fov_y_degrees: 30.0,
            ortho_half_height: 2.0,
            near: 0.1,
            far: 100.0,
        };
        let f = r.render(&camera, 64, 64);
        let px = |x: u32, y: u32| {
            let i = ((y * f.width + x) * 4) as usize;
            [f.bgra[i], f.bgra[i + 1], f.bgra[i + 2]]
        };
        // Centre: red (B, G small; R large); corner: background.
        let centre = px(32, 32);
        assert!(centre[2] > 150 && centre[1] < 120 && centre[0] < 120, "centre {centre:?}");
        assert_eq!(px(1, 1), [0, 0, 0]);
        // The sphere's radius is a quarter of the height: 16 pixels.
        assert!(px(32 + 13, 32)[2] > 20);
        assert_eq!(px(32 + 19, 32), [0, 0, 0]);
    }
}
