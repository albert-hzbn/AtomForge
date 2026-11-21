//! What the View menu changes (the C++ FileBrowser's view state) and the
//! viewport overlays drawn over the 3D picture (ui/MeasurementOverlay): atom
//! labels, the element legend and the orientation axes.

use af_core::structure::Structure;
use af_render::colouring::ColourMode;
use af_render::scene::{FrameView, Scene};
use af_science::atom_properties::PropertyDisplay;
use eframe::egui;
use glam::{Mat3, Vec3, Vec4};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AtomStyle {
    #[default]
    Balls,
    BallAndStick,
    SpaceFilling,
    Polyhedral,
}

impl AtomStyle {
    pub const ALL: [AtomStyle; 4] = [AtomStyle::Balls, AtomStyle::BallAndStick, AtomStyle::SpaceFilling, AtomStyle::Polyhedral];

    pub fn label(self) -> &'static str {
        match self {
            AtomStyle::Balls => "Balls",
            AtomStyle::BallAndStick => "Ball and Stick",
            AtomStyle::SpaceFilling => "Space Filling",
            AtomStyle::Polyhedral => "Polyhedral",
        }
    }

    pub fn hint(self) -> &'static str {
        match self {
            AtomStyle::Balls => "Spheres at the element radii, no bonds.",
            AtomStyle::BallAndStick => "Smaller spheres joined by bonds.",
            AtomStyle::SpaceFilling => "Spheres 1.8 times the element radii, no bonds.",
            AtomStyle::Polyhedral => "Coordination polyhedra with bonds (needs atoms).",
        }
    }

    /// FileBrowser::getAtomRadiusScale.
    pub fn radius_scale(self) -> f32 {
        match self {
            AtomStyle::SpaceFilling => 1.8,
            AtomStyle::Polyhedral => 0.7,
            _ => 1.0,
        }
    }
}

/// What Label Atoms With shows (AtomLabelContent).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LabelContent {
    #[default]
    Element,
    Property,
    ElementAndProperty,
    Index,
}

impl LabelContent {
    pub const ALL: [LabelContent; 4] = [LabelContent::Element, LabelContent::Property, LabelContent::ElementAndProperty, LabelContent::Index];

    pub fn label(self) -> &'static str {
        match self {
            LabelContent::Element => "Element",
            LabelContent::Property => "Per-Atom Property Value",
            LabelContent::ElementAndProperty => "Element and Property Value",
            LabelContent::Index => "Atom Index",
        }
    }
}

/// The view state of the C++ FileBrowser (its defaults: balls, no bonds,
/// orthographic projection).
#[derive(Clone, Copy, Debug)]
pub struct ViewOptions {
    pub style: AtomStyle,
    pub colour_mode: ColourMode,
    pub property: PropertyDisplay,
    pub show_bonds: bool,
    pub show_atoms: bool,
    pub show_box: bool,
    pub show_labels: bool,
    pub label_content: LabelContent,
    pub show_legend: bool,
    pub orthographic: bool,
}

impl Default for ViewOptions {
    fn default() -> Self {
        Self {
            style: AtomStyle::Balls,
            colour_mode: ColourMode::Element,
            property: PropertyDisplay::default(),
            show_bonds: false,
            show_atoms: true,
            show_box: true,
            show_labels: false,
            label_content: LabelContent::Element,
            show_legend: true,
            orthographic: true,
        }
    }
}

impl ViewOptions {
    /// FileBrowser::isShowBondsEnabled: the bonded styles always draw bonds.
    pub fn bonds_drawn(&self) -> bool {
        matches!(self.style, AtomStyle::BallAndStick | AtomStyle::Polyhedral) || self.show_bonds
    }
}

pub fn has_property(s: &Structure) -> bool {
    !s.atoms.is_empty() && s.atom_property.len() == s.atoms.len()
}

/// The text of an atom's label (empty: no label).
pub fn label_text(s: &Structure, index: usize, content: LabelContent) -> String {
    let value = || -> String {
        if s.atom_property.len() != s.atoms.len() {
            return String::new();
        }
        let v = s.atom_property[index];
        if !v.is_finite() {
            return "-".into();
        }
        af_core::format::format_g(v, 4)
    };
    match content {
        LabelContent::Property => value(),
        LabelContent::ElementAndProperty => {
            let v = value();
            if v.is_empty() {
                s.atoms[index].symbol.clone()
            } else {
                format!("{} {v}", s.atoms[index].symbol)
            }
        }
        LabelContent::Index => index.to_string(),
        LabelContent::Element => s.atoms[index].symbol.clone(),
    }
}

/// A world point in framebuffer pixels (projectToScreen).
pub fn project(p: Vec3, frame: &FrameView) -> Option<(f32, f32)> {
    let clip = frame.projection * frame.view * Vec4::new(p.x, p.y, p.z, 1.0);
    if clip.w <= 1e-5 {
        return None;
    }
    let (nx, ny) = (clip.x / clip.w, clip.y / clip.w);
    Some(((nx * 0.5 + 0.5) * frame.width as f32, (1.0 - (ny * 0.5 + 0.5)) * frame.height as f32))
}

/// Where the labels go (drawElementLabelsOverlay): nearest atoms first, at
/// most 1200, none overlapping another. `size_of` measures a label in pixels.
pub fn place_labels(scene: &Scene, s: &Structure, frame: &FrameView, content: LabelContent, size_of: &mut dyn FnMut(&str) -> (f32, f32)) -> Vec<(f32, f32, String)> {
    const MAX_CANDIDATES: usize = 8000;
    const MAX_LABELS: usize = 1200;
    const PADDING: f32 = 2.0;
    let mut placed = Vec::new();
    if s.atoms.is_empty() {
        return placed;
    }
    let mut candidates: Vec<(usize, f32, f32, f32)> = Vec::new();
    for (i, &p) in scene.instances.positions.iter().enumerate() {
        if candidates.len() >= MAX_CANDIDATES {
            break;
        }
        let atom = scene.instances.atom_indices.get(i).copied().unwrap_or(-1);
        if atom < 0 || atom as usize >= s.atoms.len() || label_text(s, atom as usize, content).is_empty() {
            continue;
        }
        let view_pos = frame.view * Vec4::new(p.x, p.y, p.z, 1.0);
        if view_pos.z >= 0.0 {
            continue;
        }
        let Some((sx, sy)) = project(p, frame) else { continue };
        candidates.push((atom as usize, sx, sy, -view_pos.z));
    }
    af_core::cxx::sort_by(&mut candidates, |a, b| a.3 < b.3);
    let mut occupied: Vec<[f32; 4]> = Vec::new();
    for (atom, sx, sy, _) in candidates {
        let text = label_text(s, atom, content);
        let (w, h) = size_of(&text);
        let rect = [sx - w * 0.5 - PADDING, sy - h * 0.5 - PADDING, sx + w * 0.5 + PADDING, sy + h * 0.5 + PADDING];
        if rect[2] < 0.0 || rect[0] > frame.width as f32 || rect[3] < 0.0 || rect[1] > frame.height as f32 {
            continue;
        }
        if occupied.iter().any(|o| rect[0] < o[2] && rect[2] > o[0] && rect[1] < o[3] && rect[3] > o[1]) {
            continue;
        }
        placed.push((rect[0] + PADDING, rect[1] + PADDING, text));
        occupied.push(rect);
        if placed.len() >= MAX_LABELS {
            break;
        }
    }
    placed
}

/// The element legend's rows: symbol, count and colour, most common first.
pub fn legend_entries(s: &Structure) -> Vec<(String, usize, [f32; 3])> {
    let mut entries: Vec<(String, usize, [f32; 3])> = Vec::new();
    for a in &s.atoms {
        match entries.iter_mut().find(|e| e.0 == a.symbol) {
            Some(e) => e.1 += 1,
            None => entries.push((a.symbol.clone(), 1, [a.r, a.g, a.b])),
        }
    }
    af_core::cxx::sort_by(&mut entries, |a, b| a.1 > b.1);
    entries
}

fn colour(c: [f32; 3]) -> egui::Color32 {
    egui::Color32::from_rgb((c[0].clamp(0.0, 1.0) * 255.0).round() as u8, (c[1].clamp(0.0, 1.0) * 255.0).round() as u8, (c[2].clamp(0.0, 1.0) * 255.0).round() as u8)
}

/// Draws the labels over the viewport rectangle (pixels to points by `scale`).
pub fn draw_labels(painter: &egui::Painter, rect: egui::Rect, scale: f32, scene: &Scene, s: &Structure, frame: &FrameView, content: LabelContent, light: bool) {
    let font = egui::FontId::proportional(14.0);
    let mut size_of = |text: &str| {
        let galley = painter.layout_no_wrap(text.to_string(), font.clone(), egui::Color32::WHITE);
        (galley.size().x * scale, galley.size().y * scale)
    };
    let text_colour = if light { egui::Color32::from_rgb(15, 15, 15) } else { egui::Color32::WHITE };
    for (x, y, text) in place_labels(scene, s, frame, content, &mut size_of) {
        painter.text(rect.min + egui::vec2(x / scale, y / scale), egui::Align2::LEFT_TOP, text, font.clone(), text_colour);
    }
}

/// The element legend at the top left; a click on a row returns its symbol
/// (and whether Ctrl was held, to add to the selection).
pub fn draw_legend(ui: &mut egui::Ui, rect: egui::Rect, s: &Structure) -> Option<(String, bool)> {
    let entries = legend_entries(s);
    if entries.len() < 2 {
        return None;
    }
    let mut clicked = None;
    // The C++ legend window (drawElementLegend): rows of the text height, a
    // swatch, the symbol and the count at the right, sizes scaled by the font
    // (18 px / 16).
    let scale = 18.0 / 16.0;
    let p = crate::theme::palette(ui.visuals().dark_mode);
    let font = egui::TextStyle::Body.resolve(ui.style());
    let line = ui.fonts_mut(|f| f.row_height(&font));
    let shown = entries.len().min(12);
    let width = |t: &str| ui.painter().layout_no_wrap(t.to_string(), font.clone(), p.text).size().x;
    let count_width = entries[..shown].iter().map(|e| width(&e.1.to_string())).fold(0.0f32, f32::max);
    let symbol_width = entries[..shown].iter().map(|e| width(&e.0)).fold(0.0f32, f32::max);
    let swatch = 6.0 * scale;
    let row_width = swatch * 2.0 + 8.0 * scale + symbol_width + 16.0 * scale + count_width;
    let area = egui::Rect::from_min_max(rect.min + egui::vec2(14.0 * scale, 14.0 * scale), rect.max);
    ui.scope_builder(egui::UiBuilder::new().max_rect(area).layout(egui::Layout::top_down(egui::Align::Min)), |ui| {
        egui::Frame::new()
            .fill(p.popup.gamma_multiply(0.9))
            .stroke(egui::Stroke::new(1.0, p.border))
            .corner_radius(8.0 * scale)
            .inner_margin(8.0 * scale)
            .show(ui, |ui| {
                ui.spacing_mut().item_spacing = egui::vec2(6.0 * scale, 2.0 * scale);
                for (symbol, count, c) in &entries[..shown] {
                    let (row, response) = ui.allocate_exact_size(egui::vec2(row_width, line), egui::Sense::click());
                    if response.hovered() {
                        ui.painter().rect_filled(row.expand(2.0), 3.0, ui.visuals().selection.bg_fill.gamma_multiply(0.6));
                    }
                    let centre = row.center().y;
                    let painter = ui.painter();
                    painter.circle(egui::pos2(row.left() + swatch, centre), swatch, colour(*c), egui::Stroke::new(1.0, p.text_muted));
                    painter.text(egui::pos2(row.left() + swatch * 2.0 + 8.0 * scale, centre), egui::Align2::LEFT_CENTER, symbol, font.clone(), p.text);
                    painter.text(egui::pos2(row.right(), centre), egui::Align2::RIGHT_CENTER, count.to_string(), font.clone(), p.text_muted);
                    let response = response.on_hover_text(format!("Select the {count} {symbol} atoms (Ctrl: add them to the selection)"));
                    if response.clicked() {
                        clicked = Some((symbol.clone(), ui.input(|i| i.modifiers.command)));
                    }
                }
                if entries.len() > shown {
                    ui.label(egui::RichText::new(format!("+{} more", entries.len() - shown)).color(p.text_muted));
                }
            });
    });
    clicked
}

/// The orientation axes at the bottom left (drawOrientationAxesGizmo): the
/// cell vectors a, b, c (or x, y, z without a cell) as the camera sees them.
pub fn draw_axes(painter: &egui::Painter, rect: egui::Rect, s: &Structure, frame: &FrameView, light: bool) {
    if rect.width() < 120.0 || rect.height() < 120.0 {
        return;
    }
    let origin = egui::pos2(rect.min.x + 73.0, rect.max.y - 73.0);
    let radius = 53.0;
    let rgba = |r, g, b, a| egui::Color32::from_rgba_unmultiplied(r, g, b, a);
    painter.circle_filled(origin, radius, if light { rgba(235, 238, 242, 200) } else { rgba(18, 22, 28, 190) });
    painter.circle_stroke(origin, radius, egui::Stroke::new(1.0, if light { rgba(0, 0, 0, 55) } else { rgba(255, 255, 255, 55) }));
    painter.circle_stroke(origin, radius - 7.0, egui::Stroke::new(1.0, if light { rgba(0, 0, 0, 20) } else { rgba(255, 255, 255, 20) }));
    let rotation = Mat3::from_mat4(frame.view);
    let (world, labels) = if s.has_unit_cell {
        let v = |k: usize| Vec3::new(s.cell_vectors[k][0] as f32, s.cell_vectors[k][1] as f32, s.cell_vectors[k][2] as f32).normalize_or_zero();
        ([v(0), v(1), v(2)], ["a", "b", "c"])
    } else {
        ([Vec3::X, Vec3::Y, Vec3::Z], ["X", "Y", "Z"])
    };
    let colours = if light {
        [rgba(200, 50, 50, 255), rgba(40, 160, 50, 255), rgba(40, 110, 220, 255)]
    } else {
        [rgba(235, 92, 92, 255), rgba(110, 220, 120, 255), rgba(110, 175, 255, 255)]
    };
    let dirs: Vec<Vec3> = world.iter().map(|&w| (rotation * w).normalize_or_zero()).collect();
    let mut order = [0usize, 1, 2];
    af_core::cxx::sort_by(&mut order, |&l, &r| dirs[l].z > dirs[r].z);
    let (scale, length) = (53.0f32, 0.77f32);
    for k in order {
        let d = dirs[k];
        let depth = -d.z * length;
        let end = egui::pos2(origin.x + d.x * length * scale, origin.y - d.y * length * scale);
        let depth01 = ((depth + 1.0) * 0.5).clamp(0.0, 1.0);
        painter.line_segment([origin, end], egui::Stroke::new(1.8 + 1.2 * depth01, colours[k]));
        let r = 2.5 + 2.0 * depth01;
        painter.circle_filled(end, r, colours[k]);
        painter.circle_stroke(end, r + 1.0, egui::Stroke::new(1.0, if light { rgba(30, 30, 40, 200) } else { rgba(245, 245, 250, 200) }));
        painter.text(end + egui::vec2(5.5, 4.0), egui::Align2::LEFT_TOP, labels[k], egui::FontId::proportional(14.0), colours[k]);
    }
    painter.circle_filled(origin, 3.2, if light { rgba(20, 20, 25, 240) } else { rgba(240, 240, 245, 240) });
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::structure::AtomSite;

    fn atoms(symbols: &[&str]) -> Structure {
        let mut s = Structure::default();
        for (i, sym) in symbols.iter().enumerate() {
            s.atoms.push(AtomSite { symbol: sym.to_string(), x: i as f64 * 3.0, ..Default::default() });
        }
        s
    }

    #[test]
    fn labels_and_legend() {
        let mut s = atoms(&["Cu", "Ni", "Cu", "O"]);
        assert_eq!(label_text(&s, 1, LabelContent::Element), "Ni");
        assert_eq!(label_text(&s, 1, LabelContent::Property), "");
        assert_eq!(label_text(&s, 1, LabelContent::ElementAndProperty), "Ni");
        assert_eq!(label_text(&s, 3, LabelContent::Index), "3");
        s.atom_property = vec![1.0, 2.345678, f64::NAN, 1e-7];
        assert_eq!(label_text(&s, 1, LabelContent::Property), "2.346");
        assert_eq!(label_text(&s, 2, LabelContent::ElementAndProperty), "Cu -");
        assert_eq!(label_text(&s, 3, LabelContent::Property), "1e-07");
        let legend = legend_entries(&s);
        assert_eq!(legend[0].0, "Cu");
        assert_eq!(legend[0].1, 2);
        assert_eq!(legend.len(), 3);
    }

    #[test]
    fn styles() {
        assert_eq!(AtomStyle::SpaceFilling.radius_scale(), 1.8);
        let mut v = ViewOptions::default();
        assert!(!v.bonds_drawn());
        v.style = AtomStyle::BallAndStick;
        assert!(v.bonds_drawn());
    }
}
