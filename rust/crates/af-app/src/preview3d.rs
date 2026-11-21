//! The dialogs' 3D previews (the C++ preview FBOs: CSLGrainBoundaryDialog,
//! DislocationBuilderDialog, ...): the structure drawn by the real renderer,
//! with its periodic boundary images, bonds and cell, from the isometric
//! camera of the main view, orbited by dragging and zoomed by the wheel.
//!
//! A dialog calls `show` while it lays out; the app renders every preview
//! requested this frame (`Previews::render`) before egui paints, so the image
//! is the current frame's.

use std::cell::RefCell;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use af_core::display::DisplaySettings;
use af_core::structure::Structure;
use af_render::colouring::{apply_colour_mode, ColourMode};
use af_render::scene::{apply_default_view, build_frame_view, build_instance_data, BondStyle, OrbitCamera, Scene};
use af_render::structure::{DrawOptions, Look, StructureRenderer};
use eframe::egui;

use crate::document::IDENTITY;

/// The camera of one preview (kept in egui's memory under the preview's id).
#[derive(Clone, Copy, PartialEq)]
struct View {
    yaw: f32,
    pitch: f32,
    zoom: f32,
}

impl Default for View {
    fn default() -> Self {
        Self { yaw: 45.0, pitch: 35.0, zoom: 1.0 }
    }
}

struct Request {
    id: egui::Id,
    /// The structure when it changed since it was last drawn.
    structure: Option<Structure>,
    fingerprint: u64,
    size: [u32; 2],
    view: View,
    background: [f32; 4],
}

thread_local! {
    static REQUESTS: RefCell<Vec<Request>> = const { RefCell::new(Vec::new()) };
    /// The texture of each preview and the structure it shows.
    static TEXTURES: RefCell<HashMap<egui::Id, (egui::TextureId, u64)>> = RefCell::new(HashMap::new());
}

/// A cheap identity of a structure's contents (atoms, elements, cell).
fn fingerprint(s: &Structure) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    s.atoms.len().hash(&mut h);
    s.has_unit_cell.hash(&mut h);
    for v in s.cell_vectors.iter().flatten() {
        v.to_bits().hash(&mut h);
    }
    for a in &s.atoms {
        a.atomic_number.hash(&mut h);
        a.x.to_bits().hash(&mut h);
        a.y.to_bits().hash(&mut h);
        a.z.to_bits().hash(&mut h);
    }
    h.finish()
}

/// A preview of `s` filling `size`; `hint` is shown while it has no atoms.
pub fn show(ui: &mut egui::Ui, size: egui::Vec2, s: &Structure, hint: &[&str]) -> egui::Response {
    let id = ui.next_auto_id().with("preview3d");
    let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click_and_drag());
    let p = crate::theme::palette(ui.visuals().dark_mode);
    let background = p.surface;
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 4.0, background);
    if s.atoms.is_empty() {
        let line = ui.text_style_height(&egui::TextStyle::Body);
        for (i, text) in hint.iter().enumerate() {
            let y = rect.center().y + (i as f32 - (hint.len() as f32 - 1.0) * 0.5) * line * 1.2;
            painter.text(egui::pos2(rect.center().x, y), egui::Align2::CENTER_CENTER, *text, egui::FontId::proportional(13.0), p.text_muted);
        }
        painter.rect_stroke(rect, 4.0, egui::Stroke::new(1.0, p.border), egui::StrokeKind::Inside);
        return response;
    }
    let mut view: View = ui.ctx().data(|d| d.get_temp(id)).unwrap_or_default();
    if response.dragged() {
        let d = response.drag_delta();
        view.yaw -= d.x * 0.5;
        view.pitch = (view.pitch + d.y * 0.5).clamp(-89.0, 89.0);
    }
    if response.hovered() {
        let wheel = ui.input(|i| i.smooth_scroll_delta.y);
        if wheel != 0.0 {
            view.zoom = (view.zoom * (1.0 + wheel * 0.002)).clamp(0.1, 20.0);
        }
    }
    ui.ctx().data_mut(|d| d.insert_temp(id, view));
    let k = ui.ctx().pixels_per_point();
    let px = [(rect.width() * k).round().max(1.0) as u32, (rect.height() * k).round().max(1.0) as u32];
    let print = fingerprint(s);
    let known = TEXTURES.with(|t| t.borrow().get(&id).copied());
    let changed = known.is_none_or(|(_, f)| f != print);
    // The renderer takes the background as the main view's are given: sRGB components.
    let srgb = |c: u8| c as f32 / 255.0;
    REQUESTS.with(|r| {
        r.borrow_mut().push(Request { id, structure: changed.then(|| s.clone()), fingerprint: print, size: px, view, background: [srgb(background.r()), srgb(background.g()), srgb(background.b()), 1.0] })
    });
    if let Some((texture, _)) = known {
        painter.image(texture, rect, egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)), egui::Color32::WHITE);
    } else {
        ui.ctx().request_repaint();
    }
    painter.rect_stroke(rect, 4.0, egui::Stroke::new(1.0, p.border), egui::StrokeKind::Inside);
    response
}

struct Slot {
    renderer: StructureRenderer,
    scene: Scene,
    has_cell: bool,
    fingerprint: u64,
    texture: Option<egui::TextureId>,
    /// What the texture shows (view, size, structure): unchanged frames are not redrawn.
    drawn: Option<(View, [u32; 2], u64)>,
}

/// The renderers of the previews on screen.
#[derive(Default)]
pub struct Previews {
    slots: HashMap<egui::Id, Slot>,
}

impl Previews {
    /// Draws the previews requested this frame; drops the ones no longer shown.
    pub fn render(&mut self, state: &eframe::egui_wgpu::RenderState, display: &DisplaySettings, colours: &[[f32; 3]], radii: &[f32]) {
        let requests = REQUESTS.with(|r| std::mem::take(&mut *r.borrow_mut()));
        let shown: std::collections::HashSet<egui::Id> = requests.iter().map(|r| r.id).collect();
        let stale: Vec<egui::Id> = self.slots.keys().filter(|id| !shown.contains(id)).copied().collect();
        for id in stale {
            if let Some(slot) = self.slots.remove(&id) {
                if let Some(texture) = slot.texture {
                    state.renderer.write().free_texture(&texture);
                }
            }
            TEXTURES.with(|t| t.borrow_mut().remove(&id));
        }
        for request in requests {
            let slot = self.slots.entry(request.id).or_insert_with(|| Slot {
                renderer: StructureRenderer::with_device(state.device.clone(), state.queue.clone(), state.adapter.get_info().name),
                scene: Scene::default(),
                has_cell: false,
                fingerprint: 0,
                texture: None,
                drawn: None,
            });
            if let Some(mut s) = request.structure {
                apply_colour_mode(&mut s, colours, ColourMode::Element, None);
                let shininess = vec![32.0; 119];
                let mut data = build_instance_data(&s, false, &IDENTITY, radii, &shininess);
                // Smaller spheres, so the bonds between them show (the C++ previews).
                for scale in data.scales.iter_mut() {
                    *scale *= 0.85;
                }
                slot.scene = Scene::new(data, &BondStyle::from_display(display));
                slot.has_cell = s.has_unit_cell;
                slot.fingerprint = request.fingerprint;
                slot.renderer.upload(&slot.scene);
                slot.drawn = None;
            }
            let key = (request.view, request.size, slot.fingerprint);
            if slot.drawn == Some(key) && slot.texture.is_some() {
                continue;
            }
            let (w, h) = (request.size[0] as i32, request.size[1] as i32);
            let mut camera = OrbitCamera::default();
            apply_default_view(&mut camera, &slot.scene, w, h, true);
            camera.yaw = request.view.yaw;
            camera.pitch = request.view.pitch;
            camera.distance /= request.view.zoom;
            let frame = build_frame_view(&camera, &slot.scene, false, w, h, display);
            let bg = request.background;
            let look = Look::new(display, &frame, [bg[0], bg[1], bg[2]]);
            let options = DrawOptions { background: bg, show_bonds: true, show_atoms: true, show_box: slot.has_cell, box_colour: None, selection_colour: display.selection_colour };
            let Ok(view) = slot.renderer.render_texture(&frame, &look, &options) else { continue };
            let mut egui_renderer = state.renderer.write();
            let texture = match slot.texture {
                Some(id) => {
                    egui_renderer.update_egui_texture_from_wgpu_texture(&state.device, &view, eframe::wgpu::FilterMode::Linear, id);
                    id
                }
                None => egui_renderer.register_native_texture(&state.device, &view, eframe::wgpu::FilterMode::Linear),
            };
            slot.texture = Some(texture);
            slot.drawn = Some(key);
            TEXTURES.with(|t| t.borrow_mut().insert(request.id, (texture, slot.fingerprint)));
        }
    }
}
