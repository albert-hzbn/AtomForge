//! Edit > Merge Structures (the C++ ui/MergeStructuresDialog): structures
//! loaded from files (or the one in the view) are arranged in 3D, each with a
//! translation and an XYZ rotation about its centroid, by the numeric fields
//! or by dragging the gizmo's arrows and rings in the arrangement view; Merge
//! replaces the structure in the view with all of them (no cell), one undo step.

use af_core::elements;
use af_core::structure::Structure;
use eframe::egui;
use glam::{Mat3, Vec3};

use crate::file_browser::{FileBrowser, Purpose};
use af_core::json::Json;
use af_pipeline::pipeline::{find_modifier_type, parameter};

use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

const STRUCTURE_EXTENSIONS: &[&str] = &["cif", "vasp", "poscar", "contcar", "xyz", "extxyz", "pdb", "json", "lmp", "data", "dump", "xsf", "cfg", "gen", "mol", "sdf", "res", "in", "pwi", "cell", "gro", "mol2"];
const HIT_RADIUS: f32 = 8.0;
const ARC_SEGMENTS: usize = 32;

#[derive(Clone)]
pub struct MergeEntry {
    pub name: String,
    pub path: String,
    pub structure: Structure,
    pub pivot: Vec3,
    pub translation: Vec3,
    pub rotation_deg: Vec3,
    pub enabled: bool,
    /// The structure entering a pipeline step (left in place by the step).
    pub step_input: bool,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Translate,
    Rotate,
}

struct Drag {
    mode: Mode,
    axis: usize,
    start: egui::Pos2,
    translation: Vec3,
    rotation: Vec3,
}

pub struct MergeDialog {
    open: bool,
    pub entries: Vec<MergeEntry>,
    selected: Option<usize>,
    status: String,
    show_box: bool,
    mode: Mode,
    drag: Option<Drag>,
    yaw: f32,
    pitch: f32,
    distance: f32,
    centre: Vec3,
    fit: bool,
    browser: FileBrowser,
    /// A pipeline merge step being edited, and the menu arrangement kept aside meanwhile.
    step: Option<StepEdit>,
    saved: Option<(Vec<MergeEntry>, Option<usize>)>,
}

impl Default for MergeDialog {
    fn default() -> Self {
        Self {
            open: false,
            entries: Vec::new(),
            selected: None,
            status: String::new(),
            show_box: false,
            mode: Mode::Translate,
            drag: None,
            yaw: 0.0,
            pitch: 0.0,
            distance: 12.0,
            centre: Vec3::ZERO,
            fit: true,
            browser: FileBrowser::default(),
            step: None,
            saved: None,
        }
    }
}

/// Rz * Ry * Rx (degrees), as the C++ eulerRotationXYZ.
pub fn euler_rotation_xyz(deg: Vec3) -> Mat3 {
    Mat3::from_rotation_z(deg.z.to_radians()) * Mat3::from_rotation_y(deg.y.to_radians()) * Mat3::from_rotation_x(deg.x.to_radians())
}

pub fn centroid(s: &Structure) -> Vec3 {
    if s.atoms.is_empty() {
        return Vec3::ZERO;
    }
    let sum = s.atoms.iter().fold(Vec3::ZERO, |acc, a| acc + Vec3::new(a.x as f32, a.y as f32, a.z as f32));
    sum / s.atoms.len() as f32
}

fn base_name(path: &str) -> String {
    path.rsplit(['/', '\\']).next().unwrap_or(path).to_string()
}

/// Each enabled entry's atoms, rotated about its pivot and translated (no cell).
pub fn combined(entries: &[MergeEntry]) -> Structure {
    let mut out = Structure { has_unit_cell: false, ..Default::default() };
    for e in entries.iter().filter(|e| e.enabled) {
        let rot = euler_rotation_xyz(e.rotation_deg);
        for atom in &e.structure.atoms {
            let p = Vec3::new(atom.x as f32, atom.y as f32, atom.z as f32);
            let t = e.pivot + rot * (p - e.pivot) + e.translation;
            let mut a = atom.clone();
            a.set_position([t.x as f64, t.y as f64, t.z as f64]);
            out.atoms.push(a);
        }
    }
    out
}

/// The structure in the view with its atoms replaced by the merged ones and
/// no cell (Merge Structures).
pub fn merged_into(view: &Structure, entries: &[MergeEntry]) -> Structure {
    let merged = combined(entries);
    let mut s = view.clone();
    s.atoms = merged.atoms;
    s.has_unit_cell = false;
    s.cell_offset = [0.0; 3];
    s.cell_vectors = [[0.0; 3]; 3];
    s.pbc_boundary_tol = 0.0;
    s.grain_colors.clear();
    s.grain_region_ids.clear();
    s.ipf_load_status.clear();
    s.atom_property.clear();
    s.atom_property_name.clear();
    s.dislocation_loop_points.clear();
    s
}

fn axis_colour(axis: usize, highlight: bool) -> egui::Color32 {
    let normal = [egui::Color32::from_rgb(220, 50, 50), egui::Color32::from_rgb(50, 200, 50), egui::Color32::from_rgb(50, 80, 220)];
    let bright = [egui::Color32::from_rgb(255, 130, 130), egui::Color32::from_rgb(130, 255, 130), egui::Color32::from_rgb(130, 160, 255)];
    if highlight {
        bright[axis]
    } else {
        normal[axis]
    }
}

fn segment_distance(p: egui::Pos2, a: egui::Pos2, b: egui::Pos2) -> f32 {
    let ab = b - a;
    let ap = p - a;
    let len = ab.length_sq();
    if len < 1e-6 {
        return ap.length();
    }
    let t = (ap.dot(ab) / len).clamp(0.0, 1.0);
    (ap - ab * t).length()
}

const RING_U: [Vec3; 3] = [Vec3::Y, Vec3::X, Vec3::X];
const RING_V: [Vec3; 3] = [Vec3::Z, Vec3::Z, Vec3::Y];

impl MergeDialog {
    pub fn add_path(&mut self, path: &str) -> bool {
        match af_io::loader::load_structure_from_file(path) {
            Ok(s) => {
                let name = base_name(path);
                self.status = format!("Loaded {name} (atoms={})", s.atoms.len());
                self.entries.push(MergeEntry { pivot: centroid(&s), structure: s, name, path: path.to_string(), translation: Vec3::ZERO, rotation_deg: Vec3::ZERO, enabled: true, step_input: false });
                if self.selected.is_none() {
                    self.selected = Some(0);
                }
                self.fit = true;
                true
            }
            Err(e) => {
                self.status = format!("Failed to load: {}", if e.is_empty() { path.to_string() } else { e });
                false
            }
        }
    }

    pub fn add_scene(&mut self, s: &Structure) {
        if s.atoms.is_empty() {
            self.status = "Current scene is empty.".into();
            return;
        }
        self.entries.push(MergeEntry { pivot: centroid(s), structure: s.clone(), name: "scene".into(), path: String::new(), translation: Vec3::ZERO, rotation_deg: Vec3::ZERO, enabled: true, step_input: false });
        if self.selected.is_none() {
            self.selected = Some(0);
        }
        self.status = format!("Added current scene (atoms={})", s.atoms.len());
        self.fit = true;
    }

    /// editStep: the structure entering the step stays where it is; the
    /// step's file is added, shifted by its offset.
    pub fn begin_step(&mut self, edit: StepEdit) {
        self.finish_step();
        self.saved = Some((std::mem::take(&mut self.entries), self.selected));
        self.selected = None;
        self.drag = None;
        self.status.clear();
        if !edit.input.atoms.is_empty() {
            self.entries.push(MergeEntry {
                pivot: centroid(&edit.input),
                structure: edit.input.clone(),
                name: "step input".into(),
                path: String::new(),
                translation: Vec3::ZERO,
                rotation_deg: Vec3::ZERO,
                enabled: true,
                step_input: true,
            });
            self.selected = Some(0);
        }
        let read = || -> Result<(String, Json), String> {
            let t = find_modifier_type("merge").ok_or("merge step missing")?;
            let file = parameter(&edit.parameters, &t, "file")?.string().map_err(|e| e.to_string())?.to_string();
            Ok((file, parameter(&edit.parameters, &t, "offset")?))
        };
        match read() {
            Ok((file, offset)) => {
                if !file.is_empty() && self.add_path(&file) {
                    let last = self.entries.len() - 1;
                    if let Ok(items) = offset.items() {
                        for (k, v) in items.iter().take(3).enumerate() {
                            self.entries[last].translation[k] = v.number().unwrap_or(0.0) as f32;
                        }
                    }
                    self.selected = Some(last);
                }
            }
            Err(e) => self.status = format!("Step parameters not read: {e}"),
        }
        self.step = Some(edit);
        self.fit = true;
        self.open = true;
    }

    /// finishStepEdit: the menu arrangement comes back.
    pub fn finish_step(&mut self) {
        if self.step.take().is_none() {
            return;
        }
        if let Some((entries, selected)) = self.saved.take() {
            self.entries = entries;
            self.selected = selected;
        }
        self.drag = None;
        self.fit = true;
    }

    /// The first entry loaded from a file (the step's file).
    fn first_file_entry(&self) -> Option<usize> {
        self.entries.iter().position(|e| !e.path.is_empty() && !e.step_input)
    }

    /// stepParameters: the file and its offset relative to the step input.
    pub fn step_parameters(&self, base: &Json) -> Json {
        let mut p = base.clone();
        let Some(file) = self.first_file_entry() else { return p };
        let mut offset = self.entries[file].translation;
        for e in self.entries.iter().filter(|e| e.step_input) {
            offset -= e.translation;
        }
        p.set("file", self.entries[file].path.as_str());
        p.set("offset", Json::array(vec![Json::from(offset.x as f64), Json::from(offset.y as f64), Json::from(offset.z as f64)]));
        p
    }

    fn update_step(&mut self) -> Result<Outcome, String> {
        let edit = self.step.clone().ok_or("Merge Structures is not editing a pipeline step")?;
        if self.first_file_entry().is_none() {
            return Err("Add a structure file to merge.".into());
        }
        let parameters = self.step_parameters(&edit.parameters);
        self.finish_step();
        self.open = false;
        Ok(Outcome::UpdateStep(StepEdit { parameters, ..edit }))
    }

    fn remove(&mut self, i: usize) {
        if i >= self.entries.len() {
            return;
        }
        self.entries.remove(i);
        self.selected = if self.entries.is_empty() { None } else { self.selected.map(|s| s.min(self.entries.len() - 1)) };
        self.drag = None;
        self.fit = true;
    }

    fn merge(&mut self, view: &Structure) -> Option<Outcome> {
        let s = merged_into(view, &self.entries);
        if s.atoms.is_empty() {
            self.status = "No enabled structures to merge.".into();
            return None;
        }
        self.status = format!("Merged {} atoms.", s.atoms.len());
        Some(Outcome::Edited(s, self.status.clone()))
    }

    fn rotation(&self) -> Mat3 {
        Mat3::from_rotation_x(self.pitch.to_radians()) * Mat3::from_rotation_y(self.yaw.to_radians())
    }

    /// Screen position and depth (larger is nearer) of a world point.
    fn project(&self, p: Vec3, rect: egui::Rect, scale: f32) -> (egui::Pos2, f32) {
        let v = self.rotation() * (p - self.centre);
        (rect.center() + egui::vec2(v.x * scale, -v.y * scale), v.z)
    }

    fn auto_fit(&mut self, atoms: &Structure) {
        self.centre = centroid(atoms);
        let mut max_r = 0.0f32;
        for a in &atoms.atoms {
            let r = elements::covalent_radius(a.atomic_number).max(0.3);
            max_r = max_r.max((Vec3::new(a.x as f32, a.y as f32, a.z as f32) - self.centre).length() + r);
        }
        self.distance = max_r.max(1.0) * 1.2;
    }

    fn preview(&mut self, ui: &mut egui::Ui, colours: &[[f32; 3]]) {
        let combined = combined(&self.entries);
        // A fixed height: the panel around it grows with what it holds.
        let size = egui::vec2(ui.available_width(), 510.0);
        let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click_and_drag());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 4.0, ui.visuals().extreme_bg_color);
        if combined.atoms.is_empty() {
            painter.text(rect.center(), egui::Align2::CENTER_CENTER, "The 3D arrangement appears here once structures are added.", egui::FontId::proportional(crate::theme::BASE_FONT), ui.visuals().weak_text_color());
            return;
        }
        if std::mem::take(&mut self.fit) {
            self.auto_fit(&combined);
        }
        let scale = 0.5 * rect.height().min(rect.width()) / self.distance;
        let gizmo = self.distance * 0.3;
        // Atoms, back to front, coloured by element; the selected entry outlined.
        let mut atoms: Vec<(f32, egui::Pos2, f32, egui::Color32, bool)> = Vec::new();
        for (ei, e) in self.entries.iter().enumerate().filter(|(_, e)| e.enabled) {
            let rot = euler_rotation_xyz(e.rotation_deg);
            for a in &e.structure.atoms {
                let p = e.pivot + rot * (Vec3::new(a.x as f32, a.y as f32, a.z as f32) - e.pivot) + e.translation;
                let (at, depth) = self.project(p, rect, scale);
                let z = a.atomic_number;
                let c = if z >= 0 && (z as usize) < colours.len() { colours[z as usize] } else { [a.r, a.g, a.b] };
                let radius = (elements::covalent_radius(z).max(0.3) * 0.75 * scale).max(1.5);
                atoms.push((depth, at, radius, crate::theme::srgb(c[0], c[1], c[2]).into(), Some(ei) == self.selected));
            }
        }
        atoms.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (_, at, r, c, selected) in &atoms {
            painter.circle_filled(*at, *r, *c);
            let edge = if *selected { egui::Stroke::new(1.5, egui::Color32::WHITE) } else { egui::Stroke::new(0.8, egui::Color32::from_black_alpha(140)) };
            painter.circle_stroke(*at, *r, edge);
        }
        if self.show_box {
            let (lo, hi) = combined.atoms.iter().fold((Vec3::splat(f32::MAX), Vec3::splat(f32::MIN)), |(lo, hi), a| {
                let p = Vec3::new(a.x as f32, a.y as f32, a.z as f32);
                (lo.min(p), hi.max(p))
            });
            let corner = |i: usize| Vec3::new(if i & 1 != 0 { hi.x } else { lo.x }, if i & 2 != 0 { hi.y } else { lo.y }, if i & 4 != 0 { hi.z } else { lo.z });
            for i in 0..8usize {
                for bit in [1usize, 2, 4] {
                    if i & bit == 0 {
                        let (a, _) = self.project(corner(i), rect, scale);
                        let (b, _) = self.project(corner(i | bit), rect, scale);
                        painter.line_segment([a, b], egui::Stroke::new(1.0, egui::Color32::from_gray(170)));
                    }
                }
            }
        }
        // The gizmo of the selected entry: arrows translate, rings rotate.
        let selected = self.selected.filter(|&i| i < self.entries.len() && self.entries[i].enabled);
        if let Some(si) = selected {
            let e = &self.entries[si];
            let centre = e.pivot + e.translation;
            let (c2, _) = self.project(centre, rect, scale);
            let active = self.drag.as_ref().map(|d| d.axis);
            let ring = |axis: usize, s: usize| centre + (RING_U[axis] * (s as f32 / ARC_SEGMENTS as f32 * std::f32::consts::TAU).cos() + RING_V[axis] * (s as f32 / ARC_SEGMENTS as f32 * std::f32::consts::TAU).sin()) * gizmo;
            painter.circle_stroke(c2, (gizmo * 1.1 * scale).max(14.0), egui::Stroke::new(1.8, egui::Color32::from_white_alpha(96)));
            for axis in 0..3 {
                let col = axis_colour(axis, active == Some(axis));
                let (tip, _) = self.project(centre + Vec3::AXES[axis] * gizmo, rect, scale);
                painter.line_segment([c2, tip], egui::Stroke::new(if self.mode == Mode::Translate { 2.8 } else { 1.8 }, col));
                let dir = tip - c2;
                if dir.length() > 5.0 {
                    let head = dir.length() * 0.15;
                    let nd = dir.normalized();
                    let perp = egui::vec2(-nd.y, nd.x);
                    let base = tip - nd * head;
                    painter.add(egui::Shape::convex_polygon(vec![tip, base + perp * head * 0.5, base - perp * head * 0.5], col, egui::Stroke::NONE));
                }
                painter.text(tip + egui::vec2(4.0, -8.0), egui::Align2::LEFT_TOP, ["X", "Y", "Z"][axis], egui::FontId::proportional(12.0), col);
                let points: Vec<egui::Pos2> = (0..=ARC_SEGMENTS).map(|s| self.project(ring(axis, s), rect, scale).0).collect();
                painter.add(egui::Shape::line(points, egui::Stroke::new(if self.mode == Mode::Rotate { 2.8 } else { 1.8 }, col)));
                let (label, _) = self.project(centre + (RING_U[axis] + RING_V[axis]) * 0.707 * gizmo, rect, scale);
                painter.text(label + egui::vec2(4.0, -8.0), egui::Align2::LEFT_TOP, ["Rx", "Ry", "Rz"][axis], egui::FontId::proportional(12.0), col);
            }
            painter.circle_filled(c2, 4.0, egui::Color32::from_white_alpha(200));
            painter.circle_stroke(c2, 4.0, egui::Stroke::new(1.5, egui::Color32::from_black_alpha(200)));
        }
        // Interaction.
        let pointer = response.interact_pointer_pos().or(response.hover_pos());
        if response.drag_started_by(egui::PointerButton::Primary) {
            if let (Some(si), Some(p)) = (selected, response.interact_pointer_pos()) {
                let e = &self.entries[si];
                let centre = e.pivot + e.translation;
                let (c2, _) = self.project(centre, rect, scale);
                let mut best = (HIT_RADIUS, None);
                for axis in 0..3 {
                    let (tip, _) = self.project(centre + Vec3::AXES[axis] * gizmo, rect, scale);
                    let d = segment_distance(p, c2, tip);
                    if d < best.0 {
                        best = (d, Some((Mode::Translate, axis)));
                    }
                }
                for axis in 0..3 {
                    for s in 0..ARC_SEGMENTS {
                        let point = |s: usize| {
                            let t = s as f32 / ARC_SEGMENTS as f32 * std::f32::consts::TAU;
                            centre + (RING_U[axis] * t.cos() + RING_V[axis] * t.sin()) * gizmo
                        };
                        let d = segment_distance(p, self.project(point(s), rect, scale).0, self.project(point(s + 1), rect, scale).0);
                        if d < best.0 {
                            best = (d, Some((Mode::Rotate, axis)));
                        }
                    }
                }
                if let Some((mode, axis)) = best.1 {
                    self.drag = Some(Drag { mode, axis, start: p, translation: e.translation, rotation: e.rotation_deg });
                }
            }
        }
        if let (Some(drag), Some(p), Some(si)) = (&self.drag, pointer, selected) {
            if response.dragged_by(egui::PointerButton::Primary) {
                let d = p - drag.start;
                let (mode, axis, t0, r0) = (drag.mode, drag.axis, drag.translation, drag.rotation);
                let centre = self.entries[si].pivot + self.entries[si].translation;
                let (c2, _) = self.project(centre, rect, scale);
                let (a2, _) = self.project(centre + Vec3::AXES[axis], rect, scale);
                let e = &mut self.entries[si];
                if mode == Mode::Translate {
                    let screen_axis = a2 - c2;
                    if screen_axis.length() > 1e-3 {
                        let along = screen_axis.normalized().dot(d);
                        e.translation = t0;
                        e.translation[axis] += along / scale;
                    }
                } else {
                    e.rotation_deg = r0;
                    e.rotation_deg[axis] += (d.x - d.y) * 0.5;
                }
            }
        } else if response.dragged() {
            let delta = response.drag_delta();
            self.yaw -= delta.x * 0.45;
            self.pitch += delta.y * 0.45;
        }
        if response.drag_stopped() {
            self.drag = None;
        }
        if response.clicked() {
            if let Some(p) = response.interact_pointer_pos() {
                // The front-most atom under the cursor selects its entry.
                let mut best: Option<(f32, usize)> = None;
                for (ei, e) in self.entries.iter().enumerate().filter(|(_, e)| e.enabled) {
                    let rot = euler_rotation_xyz(e.rotation_deg);
                    for a in &e.structure.atoms {
                        let q = e.pivot + rot * (Vec3::new(a.x as f32, a.y as f32, a.z as f32) - e.pivot) + e.translation;
                        let (at, depth) = self.project(q, rect, scale);
                        let r = (elements::covalent_radius(a.atomic_number).max(0.3) * 0.75 * scale).max(1.5);
                        if (at - p).length() <= r && best.is_none_or(|b| depth > b.0) {
                            best = Some((depth, ei));
                        }
                    }
                }
                if let Some((_, ei)) = best {
                    self.selected = Some(ei);
                }
            }
        }
        if response.hovered() {
            let wheel = ui.input(|i| i.smooth_scroll_delta.y);
            if wheel != 0.0 {
                self.distance = (self.distance * (1.0 - wheel * 0.002)).clamp(0.5, 5000.0);
            }
        }
    }
}

impl StructureDialog for MergeDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        self.add_path(path);
        true
    }
    fn takes_all_drops(&self) -> bool {
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["merge"]
    }

    fn title(&self) -> &'static str {
        "Merge Structures"
    }

    fn open(&mut self, _step: &str) {
        self.finish_step();
        self.open = true;
        self.status.clear();
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        if edit.step != "merge" {
            return false;
        }
        self.begin_step(edit);
        true
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if let Some(path) = self.browser.show(ctx) {
            self.add_path(&path.to_string_lossy());
        }
        if !self.open {
            return None;
        }
        let mut result = None;
        let mut open = true;
        let mut close = false;
        let editing = self.step.as_ref().map(|s| s.step.clone());
        crate::theme::modal_window("Merge Structures").open(&mut open).default_width(1560.0).resizable(false).show(ctx, |ui| {
            ui.set_width(1530.0);
            if let Some(step) = &editing {
                crate::structure_dialog::step_banner(ui, step);
                ui.label("The step appends one structure file to the structure entering it, shifted by an offset: the first entry loaded from a file is that file and its Pos is the offset. Rotations and further entries are not part of the step.");
                match self.first_file_entry() {
                    None => {
                        ui.colored_label(egui::Color32::from_rgb(220, 80, 60), "Add a structure file to merge.");
                    }
                    Some(f) if self.entries[f].rotation_deg != Vec3::ZERO => {
                        ui.colored_label(egui::Color32::from_rgb(230, 170, 60), format!("The rotation of {} is not written to the step.", self.entries[f].name));
                    }
                    Some(_) => {}
                }
            }
            ui.horizontal(|ui| {
                if ui.button("Add Structure").clicked() {
                    self.browser.show_for(Purpose::Open, "Add structure", "", STRUCTURE_EXTENSIONS);
                }
                if ui.add_enabled(!input.source.atoms.is_empty(), egui::Button::new("Add Current Scene")).on_hover_text("Adds the structure in the view as one more entry").clicked() {
                    self.add_scene(input.source);
                }
            });
            ui.label("Arrange multiple structures in 3D before merging. Left-click atoms to select, drag gizmo handles to move/rotate, and left-drag empty space to orbit the preview.");
            crate::dialog_layout::muted(ui, "Shortcuts: T = Translate, R = Rotate, Del = Delete selected");
            ui.horizontal(|ui| {
                if ui.button("Clear All").clicked() {
                    self.entries.clear();
                    self.selected = None;
                    self.drag = None;
                    self.status = "Cleared.".into();
                    self.fit = true;
                }
                crate::dialog_layout::labelled_check(ui, &mut self.show_box, "Show Bounding Box").on_hover_text("Displays the assembly enclosure for inspection. This overlay is not a new periodic lattice definition.");
            });
            if !ctx.egui_wants_keyboard_input() && ui.ui_contains_pointer() {
                let (t, r, del) = ctx.input(|i| (i.key_pressed(egui::Key::T), i.key_pressed(egui::Key::R), i.key_pressed(egui::Key::Delete) || i.key_pressed(egui::Key::Backspace)));
                if t {
                    self.mode = Mode::Translate;
                }
                if r {
                    self.mode = Mode::Rotate;
                }
                if del {
                    if let Some(i) = self.selected {
                        self.remove(i);
                        self.status = "Deleted selected structure.".into();
                    }
                }
            }
            ui.separator();
            let height = 640.0;
            ui.horizontal_top(|ui| {
                let left = (ui.available_width() * 0.29).floor();
                crate::dialog_layout::panel(ui, egui::vec2(left, height), |ui| {
                    crate::dialog_layout::section(ui, &format!("Loaded Structures ({})", self.entries.len()));
                    let mut remove = None;
                    egui::ScrollArea::vertical().id_salt("merge entries").max_height(height - 30.0).show(ui, |ui| {
                        for (i, e) in self.entries.iter_mut().enumerate() {
                            ui.push_id(i, |ui| {
                                ui.horizontal(|ui| {
                                    if ui.selectable_label(self.selected == Some(i), &e.name).on_hover_text(if e.path.is_empty() { "The structure in the view" } else { e.path.as_str() }).clicked() {
                                        self.selected = Some(i);
                                    }
                                    ui.weak(format!("atoms={}", e.structure.atoms.len()));
                                });
                                ui.horizontal(|ui| {
                                    for v in [&mut e.translation.x, &mut e.translation.y, &mut e.translation.z] {
                                        ui.add(egui::DragValue::new(v).speed(0.05).range(-10000.0..=10000.0).fixed_decimals(3));
                                    }
                                    ui.label("Pos (Å)");
                                })
                                .response
                                .on_hover_text("The three values translate the selected component in Cartesian Angstrom. Use these fields for reproducible offsets after coarse gizmo placement.");
                                ui.horizontal(|ui| {
                                    for v in [&mut e.rotation_deg.x, &mut e.rotation_deg.y, &mut e.rotation_deg.z] {
                                        ui.add(egui::DragValue::new(v).speed(0.25).range(-3600.0..=3600.0).fixed_decimals(2));
                                    }
                                    ui.label("Rot (°)");
                                })
                                .response
                                .on_hover_text("Three rotation values orient the selected component. The preview reflects its component transform; orbiting the camera changes the view of all components together.");
                                ui.horizontal(|ui| {
                                    if ui.button("Reset").clicked() {
                                        e.translation = Vec3::ZERO;
                                        e.rotation_deg = Vec3::ZERO;
                                    }
                                    if ui.button("Remove").clicked() {
                                        remove = Some(i);
                                    }
                                });
                                ui.separator();
                            });
                        }
                        if self.entries.is_empty() {
                            crate::dialog_layout::muted(ui, "No structures loaded yet.");
                            crate::dialog_layout::muted(ui, "Add Structure loads CIF/XYZ/PDB/... files.");
                        }
                    });
                    if let Some(i) = remove {
                        self.remove(i);
                    }
                });
                let rest = ui.available_width() - 2.0;
                crate::dialog_layout::panel(ui, egui::vec2(rest, height), |ui| {
                    crate::dialog_layout::section(ui, "3D Arrangement View");
                    self.preview(ui, input.colours);
                    crate::dialog_layout::separator_text(ui, "Selected Structure Controls");
                    match self.selected.filter(|&i| i < self.entries.len()) {
                        Some(i) => {
                            ui.label(format!("{}  ({} atoms)", self.entries[i].name, self.entries[i].structure.atoms.len()));
                            crate::dialog_layout::muted(ui, "Use gizmo handles in the 3D view to move/rotate this structure.");
                            if ui.button("Delete Selected").clicked() {
                                self.remove(i);
                            }
                        }
                        None => {
                            crate::dialog_layout::muted(ui, "Click a structure in the 3D preview or list to select it.");
                        }
                    }
                });
            });
            ui.add_space(2.0);
            let (primary, enabled) = if editing.is_some() { ("Update Step", self.first_file_entry().is_some()) } else { ("Merge Structures", !self.entries.is_empty()) };
            match crate::dialog_layout::action_bar(ui, &self.status.clone(), crate::dialog_layout::Note::Info, primary, Some("Close"), enabled) {
                1 if editing.is_some() => result = self.update_step().ok(),
                1 => result = self.merge(input.source),
                2 => close = true,
                _ => {}
            }
        });
        if !self.open {
            // Update Step closed it.
            return result;
        }
        self.open = open && !close;
        if !self.open {
            self.finish_step();
        }
        result
    }

    /// add PATH | add-scene | move I X Y Z | rotate I RX RY RZ | remove I |
    /// clear | merge | expect-entries N
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let entry = |words: &[&str], n: usize| -> Result<usize, String> {
            let i: usize = words.first().and_then(|w| w.parse().ok()).ok_or(format!("{command} wants an entry number (from 0)"))?;
            if i >= n {
                return Err(format!("no entry {i}"));
            }
            Ok(i)
        };
        let vector = |words: &[&str]| -> Result<Vec3, String> {
            let v: Vec<f32> = words.iter().skip(1).map(|w| w.parse::<f32>()).collect::<Result<_, _>>().map_err(|_| format!("{command} wants three numbers"))?;
            (v.len() == 3).then(|| Vec3::new(v[0], v[1], v[2])).ok_or(format!("{command} wants three numbers"))
        };
        match command {
            "add" => {
                if !self.add_path(args.trim()) {
                    return Err(self.status.clone());
                }
            }
            "add-scene" => self.add_scene(input.source),
            "move" => {
                let i = entry(&words, self.entries.len())?;
                self.entries[i].translation = vector(&words)?;
                self.fit = true;
            }
            "rotate" => {
                let i = entry(&words, self.entries.len())?;
                self.entries[i].rotation_deg = vector(&words)?;
            }
            "remove" => {
                let i = entry(&words, self.entries.len())?;
                self.remove(i);
            }
            "clear" => {
                self.entries.clear();
                self.selected = None;
            }
            "merge" => return self.merge(input.source).map(Some).ok_or(self.status.clone()),
            "update-step" => return self.update_step().map(Some),
            "expect-entries" => {
                let n: usize = args.trim().parse().map_err(|_| "expect-entries N")?;
                if self.entries.len() != n {
                    return Err(format!("{} entries", self.entries.len()));
                }
            }
            other => return Err(format!("Merge Structures has no command {other}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::structure::AtomSite;

    #[test]
    fn merge_step_round_trip() {
        // Run from the repository root (as the app is), the example file loads.
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
        let file = root.join("docs/manual/examples/fe_bcc.cif").to_string_lossy().to_string();
        let mut input = Structure::default();
        input.atoms.push(AtomSite::new("Cu", [0.0, 0.0, 0.0]));
        let mut p = Json::object();
        p.set("file", file.as_str());
        p.set("offset", Json::array(vec![Json::from(1.0), Json::from(2.0), Json::from(3.0)]));
        p.set("select", false);
        let mut d = MergeDialog::default();
        d.begin_step(StepEdit { index: 2, step: "merge".into(), parameters: p.clone(), input, ..Default::default() });
        assert_eq!(d.entries.len(), 2);
        assert!(d.entries[0].step_input);
        assert_eq!(d.entries[1].translation, Vec3::new(1.0, 2.0, 3.0));
        // The offset is relative to the step input, which the merge leaves in place.
        d.entries[0].translation = Vec3::new(0.5, 0.0, 0.0);
        let back = d.step_parameters(&p);
        assert_eq!(back.find("offset").unwrap().dump(-1), "[0.5,2,3]");
        assert_eq!(back.find("file").unwrap().string().unwrap(), file);
        assert!(!back.find("select").unwrap().boolean().unwrap());
        // Ending the edit brings back the menu arrangement.
        d.finish_step();
        assert!(d.entries.is_empty());
    }

    #[test]
    fn rotation_about_pivot_and_translation() {
        let mut s = Structure::default();
        s.atoms.push(AtomSite::new("O", [1.0, 0.0, 0.0]));
        s.atoms.push(AtomSite::new("O", [-1.0, 0.0, 0.0]));
        let e = MergeEntry { pivot: centroid(&s), structure: s, name: "a".into(), path: String::new(), translation: Vec3::new(0.0, 0.0, 5.0), rotation_deg: Vec3::new(0.0, 0.0, 90.0), enabled: true, step_input: false };
        let mut off = e.clone();
        off.enabled = false;
        let c = combined(&[e, off]);
        assert_eq!(c.atoms.len(), 2);
        let p = c.atoms[0].position();
        assert!((p[0]).abs() < 1e-6 && (p[1] - 1.0).abs() < 1e-6 && (p[2] - 5.0).abs() < 1e-6, "{p:?}");
        assert!(!c.has_unit_cell);
    }
}
