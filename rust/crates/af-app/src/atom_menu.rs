//! The right-click menu of selected atoms (the C++ ui/AtomContextMenu):
//! substitute the selected atoms (S), insert an atom at their midpoint (I),
//! measure, Atom Info, delete and deselect; the element comes from a
//! periodic table.

use af_core::elements;
use af_core::structure::{AtomSite, Structure};
use af_render::scene::Scene;
use eframe::egui;


/// What the element picker is for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pick {
    Substitute,
    InsertMidpoint,
}

/// Gives the selected atoms element `z` (each atom once); returns how many changed.
pub fn substitute(structure: &mut Structure, scene: &Scene, selection: &[usize], z: i32, colours: &[[f32; 3]]) -> usize {
    let mut count = 0;
    for &i in selection {
        let Some(&base) = scene.instances.atom_indices.get(i) else { continue };
        let Some(atom) = (base >= 0).then(|| structure.atoms.get_mut(base as usize)).flatten() else { continue };
        atom.symbol = elements::symbol(z).to_string();
        atom.atomic_number = z;
        [atom.r, atom.g, atom.b] = colour_of(z, colours);
        count += 1;
    }
    count
}

fn colour_of(z: i32, colours: &[[f32; 3]]) -> [f32; 3] {
    if z >= 0 && (z as usize) < colours.len() {
        colours[z as usize]
    } else {
        elements::default_colour(z)
    }
}

/// Adds an atom of element `z` at the mean of the selected instances'
/// positions (as drawn, so periodic images count where they are); at least two.
pub fn insert_midpoint(structure: &mut Structure, scene: &Scene, selection: &[usize], z: i32, colours: &[[f32; 3]]) -> bool {
    let (mut sum, mut n) = ([0.0f64; 3], 0usize);
    for &i in selection {
        let Some(p) = scene.instances.positions.get(i) else { continue };
        sum[0] += p.x as f64;
        sum[1] += p.y as f64;
        sum[2] += p.z as f64;
        n += 1;
    }
    if n < 2 {
        return false;
    }
    let [r, g, b] = colour_of(z, colours);
    structure.atoms.push(AtomSite {
        symbol: elements::symbol(z).to_string(),
        atomic_number: z,
        x: sum[0] / n as f64,
        y: sum[1] / n as f64,
        z: sum[2] / n as f64,
        r,
        g,
        b,
        ..Default::default()
    });
    // Per-atom data stays aligned with the atoms.
    if !structure.grain_colors.is_empty() {
        structure.grain_colors.push([r, g, b]);
    }
    if !structure.grain_region_ids.is_empty() {
        structure.grain_region_ids.push(-1);
    }
    if !structure.atom_property.is_empty() {
        structure.atom_property.push(f64::NAN);
    }
    if !structure.frozen_property.is_empty() {
        structure.frozen_property.push(f64::NAN);
    }
    true
}

/// The element categories of the periodic table and their colours (kCatColors).
fn category(z: i32) -> usize {
    match z {
        1 => 0,
        3 | 11 | 19 | 37 | 55 | 87 => 1,
        4 | 12 | 20 | 38 | 56 | 88 => 2,
        57..=71 => 3,
        89..=103 => 4,
        13 | 31 | 49 | 50 | 81..=84 | 113..=116 => 6,
        5 | 14 | 32 | 33 | 51 | 52 => 7,
        6 | 7 | 8 | 15 | 16 | 34 => 8,
        9 | 17 | 35 | 53 | 85 | 117 => 9,
        2 | 10 | 18 | 36 | 54 | 86 | 118 => 10,
        _ => 5,
    }
}

const CATEGORY_COLOURS: [[f32; 3]; 11] = [
    [0.78, 0.78, 0.20],
    [0.92, 0.40, 0.25],
    [0.92, 0.76, 0.30],
    [1.00, 0.72, 0.77],
    [1.00, 0.78, 0.60],
    [0.72, 0.72, 0.95],
    [0.55, 0.85, 0.80],
    [0.72, 0.90, 0.55],
    [0.60, 0.90, 0.58],
    [0.55, 0.90, 0.70],
    [0.50, 0.78, 0.95],
];

fn rgb(c: [f32; 3], brighter: f32, alpha: f32) -> egui::Color32 {
    let f = |v: f32| ((v + brighter).min(1.0) * 255.0).round() as u8;
    egui::Color32::from_rgba_unmultiplied(f(c[0]), f(c[1]), f(c[2]), (alpha * 255.0).round() as u8)
}

/// The periodic table of the C++ picker (drawPeriodicTable): category colours,
/// the f-block below with markers where it belongs; a click selects, a
/// double-click applies. Returns (selection changed, double-clicked).
pub fn category_table(ui: &mut egui::Ui, selected: &mut i32) -> (bool, bool) {
    let font_size = crate::theme::BASE_FONT * 18.0 / 13.5;
    let (cell_w, cell_h) = ((font_size * 2.8).max(44.0), (font_size * 2.1).max(34.0));
    let (pad_x, gap) = (6.0, 10.0);
    let row_y = |row: usize| if row <= 6 { row as f32 * cell_h } else { 7.0 * cell_h + gap + (row - 7) as f32 * cell_h };
    let (rect, _) = ui.allocate_exact_size(egui::vec2(2.0 * pad_x + 18.0 * cell_w, row_y(8) + cell_h + 4.0), egui::Sense::hover());
    let origin = rect.min + egui::vec2(pad_x, 4.0);
    let symbol_font = egui::FontId::proportional(crate::theme::BASE_FONT);
    let (mut changed, mut applied) = (false, false);
    for z in 1..=118 {
        let (row, col) = crate::display_settings::table_position(z);
        let r = egui::Rect::from_min_size(origin + egui::vec2(col as f32 * cell_w, row_y(row)), egui::vec2(cell_w - 1.0, cell_h - 1.0));
        let response = ui.interact(r, ui.scope_id().with(("element", z)), egui::Sense::click());
        let is_selected = *selected == z;
        let c = CATEGORY_COLOURS[category(z)];
        let bg = rgb(c, if is_selected { 0.35 } else if response.hovered() { 0.15 } else { 0.0 }, 1.0);
        ui.painter().rect_filled(r, 3.0, bg);
        if is_selected {
            ui.painter().rect_stroke(r, 3.0, egui::Stroke::new(2.0, egui::Color32::from_rgba_unmultiplied(0, 0, 100, 200)), egui::StrokeKind::Inside);
        } else {
            ui.painter().rect_stroke(r, 3.0, egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(0, 0, 0, 55)), egui::StrokeKind::Inside);
        }
        ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, elements::symbol(z), symbol_font.clone(), egui::Color32::BLACK);
        let tip = if is_selected { format!("{z}  {}\n{}\n(Selected)", elements::symbol(z), elements::name(z)) } else { format!("{z}  {}\n{}", elements::symbol(z), elements::name(z)) };
        let response = response.on_hover_text(tip);
        if response.clicked() {
            *selected = z;
            changed = true;
        }
        if response.double_clicked() {
            *selected = z;
            applied = true;
        }
    }
    // The markers where the lanthanides and actinides belong.
    for (row, label, cat) in [(5usize, "57-71", 3usize), (6, "89-103", 4)] {
        let r = egui::Rect::from_min_size(origin + egui::vec2(2.0 * cell_w, row_y(row)), egui::vec2(cell_w - 1.0, cell_h - 1.0));
        ui.painter().rect_filled(r, 3.0, rgb(CATEGORY_COLOURS[cat], 0.0, 0.55));
        ui.painter().rect_stroke(r, 3.0, egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(0, 0, 0, 55)), egui::StrokeKind::Inside);
        ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, label, egui::FontId::proportional(crate::theme::BASE_FONT * 0.8), egui::Color32::from_rgb(40, 40, 40));
    }
    (changed, applied)
}

/// The element picker (Periodic Table): click or type a symbol to choose,
/// double-click, Enter or Apply to use it (drawPeriodicTable).
#[derive(Default)]
pub struct ElementPicker {
    pub pending: Option<Pick>,
    element: i32,
    typed: String,
    typed_at: f64,
}

impl ElementPicker {
    pub fn open(&mut self, pick: Pick) {
        self.pending = Some(pick);
        self.element = -1;
        self.typed.clear();
    }

    /// Typed letters pick the element with that symbol; letters more than a
    /// second apart start a new symbol (pickTypedSymbol).
    fn type_text(&mut self, text: &str, now: f64) {
        for c in text.chars().filter(|c| c.is_ascii_alphabetic()) {
            if now - self.typed_at > 1.0 || self.typed.len() >= 2 {
                self.typed.clear();
            }
            self.typed_at = now;
            self.typed.push(c);
            if let Some(z) = (1..=118).find(|&z| elements::symbol(z).eq_ignore_ascii_case(&self.typed)) {
                self.element = z;
            }
        }
    }

    /// UI scripts: `open substitute|insert`, `type LETTERS` (the app's `element-picker` command).
    #[allow(dead_code)]
    pub fn script(&mut self, command: &str, args: &str) -> Result<(), String> {
        match command {
            "open" => self.open(if args.trim() == "insert" { Pick::InsertMidpoint } else { Pick::Substitute }),
            "type" => {
                let now = self.typed_at;
                self.type_text(args.trim(), now);
            }
            _ => return Err(format!("unknown element-picker command {command}")),
        }
        Ok(())
    }

    pub fn show(&mut self, ctx: &egui::Context, _colours: &[[f32; 3]]) -> Option<(Pick, i32)> {
        let pick = self.pending?;
        let mut open = true;
        let mut chosen = None;
        let mut cancel = false;
        let now = ctx.input(|i| i.time);
        let typed: String = ctx.input(|i| i.events.iter().filter_map(|e| if let egui::Event::Text(t) = e { Some(t.clone()) } else { None }).collect());
        if !typed.is_empty() && !ctx.egui_wants_keyboard_input() {
            self.type_text(&typed, now);
        }
        crate::theme::modal_window("Periodic Table").id(egui::Id::unique("element picker")).open(&mut open).resizable(false).show(ctx, |ui| {
            let (_, applied) = category_table(ui, &mut self.element);
            ui.add_space(4.0);
            let picked = (1..=118).contains(&self.element);
            if picked {
                ui.label(format!("{} ({}, Z = {})", elements::name(self.element), elements::symbol(self.element), self.element));
            } else {
                let muted = crate::theme::palette(ui.visuals().dark_mode).text_muted;
                ui.label(egui::RichText::new("Click an element or type its symbol; double-click or Enter applies it.").color(muted));
            }
            let enter = picked && ui.input(|i| i.key_pressed(egui::Key::Enter));
            match crate::dialog_layout::footer(ui, "Apply", Some("Cancel"), picked) {
                1 => chosen = Some((pick, self.element)),
                2 => cancel = true,
                _ => {}
            }
            if (applied || enter) && picked {
                chosen = Some((pick, self.element));
            }
        });
        if !open || cancel || chosen.is_some() {
            self.pending = None;
        }
        chosen
    }
}

/// What the context menu asks for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MenuAction {
    Pick(Pick),
    Distance,
    Angle,
    Dihedral,
    AtomInfo,
    Delete,
    Deselect,
}

/// The menu's items for a selection of `n` instances.
pub fn menu(ui: &mut egui::Ui, n: usize) -> Option<MenuAction> {
    let items = [
        ("Substitute Atom", true, MenuAction::Pick(Pick::Substitute)),
        ("Insert Atom at Midpoint", n >= 2, MenuAction::Pick(Pick::InsertMidpoint)),
        ("Measure Distance", n == 2, MenuAction::Distance),
        ("Measure Angle", n == 3, MenuAction::Angle),
        ("Measure Dihedral", n == 4, MenuAction::Dihedral),
        ("Atom Info", n == 1, MenuAction::AtomInfo),
        ("Delete", true, MenuAction::Delete),
        ("Deselect", true, MenuAction::Deselect),
    ];
    let mut action = None;
    for (label, enabled, a) in items {
        if ui.add_enabled(enabled, egui::Button::new(label)).clicked() {
            action = Some(a);
            ui.close();
        }
    }
    action
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_render::scene::{build_instance_data, BondStyle};

    #[test]
    fn substitute_and_insert() {
        let mut s = Structure::default();
        for (i, x) in [0.0, 2.0, 4.0].iter().enumerate() {
            s.atoms.push(AtomSite { symbol: "Cu".into(), atomic_number: 29, x: *x, y: i as f64, ..Default::default() });
        }
        s.atom_property = vec![1.0, 2.0, 3.0];
        let look = af_core::display::default_look();
        let scene = Scene::new(build_instance_data(&s, false, &[[1, 0, 0], [0, 1, 0], [0, 0, 1]], &look.radii, &vec![32.0; 119]), &BondStyle::default());
        assert_eq!(substitute(&mut s, &scene, &[0, 2], 28, &look.colours), 2);
        assert_eq!((s.atoms[0].symbol.as_str(), s.atoms[1].symbol.as_str(), s.atoms[2].atomic_number), ("Ni", "Cu", 28));
        assert_eq!([s.atoms[2].r, s.atoms[2].g, s.atoms[2].b], look.colours[28]);
        assert!(!insert_midpoint(&mut s, &scene, &[1], 8, &look.colours));
        assert!(insert_midpoint(&mut s, &scene, &[0, 1, 2], 8, &look.colours));
        let o = &s.atoms[3];
        assert_eq!((o.symbol.as_str(), o.x, o.y), ("O", 2.0, 1.0));
        assert_eq!(s.atom_property.len(), 4);
        assert!(s.atom_property[3].is_nan());
    }
}
