//! Edit > Add Interstitial Atoms (the C++ ui/InterstitialAtomsDialog): void
//! detection on a source structure, voids grouped by volume, and placement of
//! one element into the voids at random (by percentage or count), inside an
//! OBJ/STL mesh region, or into hand-picked voids. The preview is a small
//! orbiting projection (click a void to select it, Ctrl+click for several);
//! the mesh is positioned with numeric controls.

use af_algorithms::interstitial_voids::{detect_interstitial_voids, VoidKind, VoidParams, VoidRegion};
use af_algorithms::mesh::{parse_obj_mesh, parse_stl_mesh, Mesh};
use af_core::cxx_random::{shuffle, Mt19937};
use af_core::structure::{AtomSite, Structure};
use eframe::egui;
use glam::{Mat3, Mat4, Vec3, Vec4};

use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

/// Void kinds of the insert-interstitials step, in the order of its "kind" choice.
const STEP_KINDS: [&str; 4] = ["any", "tetrahedral", "octahedral", "irregular"];
const STEP_KIND_LABELS: [&str; 4] = ["Any", "Tetrahedral", "Octahedral", "Irregular"];

const SEED_HINT: &str = "Seed controls stochastic selection; zero is a deterministic seed in this tool, unlike builders where zero explicitly means time.";
const GRID_HINT: &str = "Grid controls the search sampling (8-30 in the GUI); Max Voids caps the retained candidates. Increasing sampling can reveal sites missed by a coarse search and increases work.";
const CLEARANCE_HINT: &str = "These are distinct detection controls in Å: clearance filters available space around a candidate; separation prevents retaining nearby duplicate candidates. They are not the requested number of inserted atoms.";

/// A small orbiting orthographic preview drawn with the painter (the C++
/// dialogs render these into framebuffers).
pub(crate) struct OrbitPreview {
    pub yaw: f32,
    pub pitch: f32,
    /// Pixels per Å.
    pub zoom: f32,
    pub centre: Vec3,
    pub fit_pending: bool,
}

impl Default for OrbitPreview {
    fn default() -> Self {
        Self { yaw: 45.0, pitch: 35.0, zoom: 20.0, centre: Vec3::ZERO, fit_pending: true }
    }
}

impl OrbitPreview {
    fn rotation(&self) -> Mat3 {
        Mat3::from_rotation_x(self.pitch.to_radians()) * Mat3::from_rotation_y(-self.yaw.to_radians())
    }

    /// Screen position and depth (larger is nearer the viewer).
    pub fn project(&self, p: Vec3, rect: egui::Rect) -> (egui::Pos2, f32) {
        let v = self.rotation() * (p - self.centre);
        (rect.center() + egui::vec2(v.x * self.zoom, -v.y * self.zoom), v.z)
    }

    /// Frames the points (with their radii) in the rect.
    pub fn fit(&mut self, points: &[(Vec3, f32)], rect: egui::Rect) {
        self.yaw = 45.0;
        self.pitch = 35.0;
        if points.is_empty() {
            self.centre = Vec3::ZERO;
            self.zoom = 20.0;
            return;
        }
        let (lo, hi) = points.iter().fold((Vec3::splat(f32::MAX), Vec3::splat(f32::MIN)), |(lo, hi), (p, _)| (lo.min(*p), hi.max(*p)));
        self.centre = (lo + hi) * 0.5;
        let extent = points.iter().map(|(p, r)| (*p - self.centre).length() + r).fold(1.0f32, f32::max);
        self.zoom = 0.5 * rect.width().min(rect.height()) / (extent * 1.1);
    }

    /// Allocates the canvas: drag orbits, the wheel zooms.
    pub fn canvas(&mut self, ui: &mut egui::Ui, size: egui::Vec2) -> (egui::Response, egui::Painter) {
        let (response, painter) = ui.allocate_painter(size, egui::Sense::click_and_drag());
        if response.dragged() {
            let d = response.drag_delta();
            self.yaw -= d.x * 0.5;
            self.pitch = (self.pitch + d.y * 0.5).clamp(-85.0, 85.0);
        }
        if response.hovered() {
            let wheel = ui.input(|i| i.smooth_scroll_delta.y);
            if wheel != 0.0 {
                self.zoom = (self.zoom * (1.0 + wheel * 0.002)).clamp(0.5, 500.0);
            }
        }
        let rect = response.rect;
        painter.rect_filled(rect, 4.0, ui.visuals().extreme_bg_color);
        painter.rect_stroke(rect, 4.0, ui.visuals().widgets.noninteractive.bg_stroke, egui::StrokeKind::Inside);
        (response, painter)
    }

    /// The cell's twelve edges, dashed (the C++ preview's box lines); drawn
    /// before the atoms, which hide the edges behind them.
    pub fn draw_cell(&self, painter: &egui::Painter, rect: egui::Rect, s: &Structure, ui: &egui::Ui) {
        if !s.has_unit_cell {
            return;
        }
        let v = |i: usize| Vec3::new(s.cell_vectors[i][0] as f32, s.cell_vectors[i][1] as f32, s.cell_vectors[i][2] as f32);
        let o = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
        let corner = |k: usize| o + v(0) * (k & 1) as f32 + v(1) * ((k >> 1) & 1) as f32 + v(2) * ((k >> 2) & 1) as f32;
        let colour = if ui.visuals().dark_mode { egui::Color32::from_gray(200) } else { egui::Color32::from_gray(60) };
        for (i, j) in [(0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7)] {
            let (a, b) = (self.project(corner(i), rect).0, self.project(corner(j), rect).0);
            painter.extend(egui::Shape::dashed_line(&[a, b], egui::Stroke::new(1.4, colour), 7.0, 4.0));
        }
    }

    /// Atoms as shaded discs, back to front: (position, colour, radius Å, alpha).
    pub fn draw_atoms(&self, painter: &egui::Painter, rect: egui::Rect, atoms: &[(Vec3, [f32; 3], f32, f32)]) {
        let mut items: Vec<(f32, egui::Pos2, [f32; 3], f32, f32)> = atoms
            .iter()
            .map(|(p, c, r, a)| {
                let (s, z) = self.project(*p, rect);
                (z, s, *c, r * self.zoom, *a)
            })
            .filter(|i| rect.expand(i.3).contains(i.1))
            .collect();
        items.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (_, at, c, radius, alpha) in items {
            let radius = radius.max(1.0);
            let rgba = |c: [f32; 3], a: f32| -> egui::Color32 { crate::theme::srgba(c[0].clamp(0.0, 1.0), c[1].clamp(0.0, 1.0), c[2].clamp(0.0, 1.0), a).into() };
            painter.circle_filled(at, radius, rgba([c[0] * 0.5, c[1] * 0.5, c[2] * 0.5], alpha));
            painter.circle_filled(at, radius * 0.86, rgba(c, alpha));
            painter.circle_filled(at - egui::vec2(radius * 0.3, radius * 0.3), radius * 0.32, rgba([c[0] * 0.5 + 0.5, c[1] * 0.5 + 0.5, c[2] * 0.5 + 0.5], alpha * 0.7));
        }
    }
}

/// The preview atoms of a structure: (position, colour, radius, alpha).
pub(crate) fn preview_atoms(s: &Structure, alpha: f32) -> Vec<(Vec3, [f32; 3], f32, f32)> {
    s.atoms.iter().map(|a| (Vec3::new(a.x as f32, a.y as f32, a.z as f32), [a.r * 0.8, a.g * 0.8, a.b * 0.8], af_core::elements::covalent_radius(a.atomic_number) * 0.6, alpha)).collect()
}

/// A hint shown centred in an empty canvas.
pub(crate) fn canvas_hint(painter: &egui::Painter, rect: egui::Rect, text: &str, ui: &egui::Ui) {
    painter.text(rect.center(), egui::Align2::CENTER_CENTER, text, egui::FontId::proportional(14.0), ui.visuals().weak_text_color());
}

fn ray_intersects_triangle(orig: Vec3, dir: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> bool {
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let pvec = dir.cross(e2);
    let det = e1.dot(pvec);
    if det.abs() < 1e-8 {
        return false;
    }
    let inv_det = 1.0 / det;
    let tvec = orig - v0;
    let u = tvec.dot(pvec) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let qvec = tvec.cross(e1);
    let v = dir.dot(qvec) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    e2.dot(qvec) * inv_det > 1e-6
}

fn mesh_transform(scale: f32, translation: Vec3, rotation_deg: Vec3) -> Mat4 {
    Mat4::from_translation(translation)
        * Mat4::from_rotation_x(rotation_deg.x.to_radians())
        * Mat4::from_rotation_y(rotation_deg.y.to_radians())
        * Mat4::from_rotation_z(rotation_deg.z.to_radians())
        * Mat4::from_scale(Vec3::splat(scale))
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Placement {
    Random,
    MeshRegion,
    Manual,
}

const PLACEMENTS: [(Placement, &str); 3] = [(Placement::Random, "Random distribution"), (Placement::MeshRegion, "Inside OBJ/STL mesh region"), (Placement::Manual, "Manual (selected voids only)")];

/// One volume group of the detected voids.
struct VoidGroup {
    bin: i32,
    count: usize,
    label: String,
    enabled: bool,
}

pub struct InterstitialDialog {
    open: bool,
    source: Option<Structure>,
    source_label: String,
    voids: Vec<VoidRegion>,
    groups: Vec<VoidGroup>,
    selected: usize,
    picked: Vec<usize>,
    grid: i32,
    max_voids: i32,
    min_clearance: f32,
    min_separation: f32,
    show_polyhedra: bool,
    element: i32,
    picking_element: bool,
    placement: Placement,
    by_percent: bool,
    percent: f32,
    count: i32,
    seed: i32,
    min_insert: f32,
    mesh: Option<Mesh>,
    mesh_label: String,
    mesh_path: String,
    structure_path: String,
    mesh_scale: f32,
    mesh_translation: Vec3,
    mesh_rotation: Vec3,
    status: (String, bool),
    preview: OrbitPreview,
    /// The pipeline step being edited (Edit in Dialog), with the step-only
    /// fields: the void kind (index into STEP_KINDS) and random site order.
    step: Option<StepEdit>,
    step_kind: usize,
    step_random: bool,
}

impl Default for InterstitialDialog {
    fn default() -> Self {
        Self {
            open: false,
            source: None,
            source_label: String::new(),
            voids: Vec::new(),
            groups: Vec::new(),
            selected: 0,
            picked: Vec::new(),
            grid: 14,
            max_voids: 0,
            min_clearance: 0.6,
            min_separation: 0.5,
            show_polyhedra: true,
            element: 6,
            picking_element: false,
            placement: Placement::Random,
            by_percent: true,
            percent: 25.0,
            count: 20,
            seed: 42,
            min_insert: 0.6,
            mesh: None,
            mesh_label: String::new(),
            mesh_path: String::new(),
            structure_path: String::new(),
            mesh_scale: 1.0,
            mesh_translation: Vec3::ZERO,
            mesh_rotation: Vec3::ZERO,
            status: (String::new(), false),
            preview: OrbitPreview::default(),
            step: None,
            step_kind: 0,
            step_random: false,
        }
    }
}

impl InterstitialDialog {
    /// editStep: the controls from the step's parameters (defaults where absent)
    /// and the voids from the structure entering the step.
    fn begin_step(&mut self, edit: StepEdit) {
        let read = |edit: &StepEdit, d: &mut Self| -> Result<(), String> {
            let t = af_pipeline::pipeline::find_modifier_type("insert-interstitials").ok_or("no insert-interstitials step type")?;
            let value = |name: &str| af_pipeline::pipeline::parameter(&edit.parameters, &t, name);
            let z = af_core::elements::atomic_number(value("element")?.string()?);
            if (1..=118).contains(&z) {
                d.element = z;
            }
            let kind = value("kind")?;
            let kind = kind.string()?;
            d.step_kind = STEP_KINDS.iter().position(|k| *k == kind).unwrap_or(0);
            d.by_percent = false;
            d.count = (value("count")?.number()? as i32).max(0);
            d.step_random = value("order")?.string()? == "random";
            d.seed = value("seed")?.number()? as i32;
            d.grid = value("resolution")?.number()? as i32;
            d.min_clearance = value("clearance")?.number()? as f32;
            d.min_separation = value("separation")?.number()? as f32;
            Ok(())
        };
        self.step_kind = 0;
        self.status = match read(&edit, self) {
            Ok(()) => (String::new(), false),
            Err(e) => (format!("Step parameters not read: {e}"), true),
        };
        // The voids are those of the structure entering the step.
        let input = edit.input.clone();
        self.step = Some(edit);
        if input.atoms.is_empty() {
            self.source = None;
            self.source_label = "step input".into();
            self.voids.clear();
            self.picked.clear();
            self.selected = 0;
            self.rebuild_groups();
        } else {
            self.set_source(input, "step input");
        }
        self.open = true;
    }

    /// stepParameters: the step's parameters with the dialog's values.
    fn step_parameters(&self) -> Option<StepEdit> {
        let edit = self.step.as_ref()?;
        let mut out = edit.clone();
        let p = &mut out.parameters;
        p.set("element", af_core::elements::symbol(self.element.clamp(1, 118)));
        p.set("kind", STEP_KINDS[self.step_kind.min(3)]);
        p.set("count", self.count.max(0));
        p.set("order", if self.step_random { "random" } else { "largest" });
        p.set("seed", self.seed);
        p.set("resolution", self.grid);
        p.set("clearance", self.min_clearance as f64);
        p.set("separation", self.min_separation as f64);
        Some(out)
    }

    /// Update Step: the new parameters, and the edit ends.
    fn update_step(&mut self) -> Option<Outcome> {
        let out = self.step_parameters()?;
        self.step = None;
        self.open = false;
        Some(Outcome::UpdateStep(out))
    }

    fn set_source(&mut self, s: Structure, label: &str) {
        self.source = Some(s);
        self.source_label = label.into();
        self.voids.clear();
        self.picked.clear();
        self.selected = 0;
        self.rebuild_groups();
        self.preview.fit_pending = true;
    }

    fn use_scene(&mut self, scene: &Structure) {
        if scene.atoms.is_empty() {
            self.status = ("Current scene has no atoms.".into(), true);
        } else {
            self.set_source(scene.clone(), "scene");
            self.status = ("Source set from current scene.".into(), false);
        }
    }

    fn load_structure(&mut self, path: &str) {
        match af_io::loader::load_structure_from_file(path) {
            Err(e) => self.status = (format!("Failed to load structure: {}", if e.is_empty() { path.to_string() } else { e }), true),
            Ok(s) if s.atoms.is_empty() => self.status = ("Loaded structure has no atoms.".into(), true),
            Ok(s) => {
                let name = std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.into());
                self.set_source(s, &name);
                self.status = (format!("Structure loaded: {name}"), false);
            }
        }
    }

    fn load_mesh(&mut self, path: &str) {
        let lower = path.to_ascii_lowercase();
        let result = if lower.ends_with(".obj") {
            parse_obj_mesh(std::path::Path::new(path))
        } else if lower.ends_with(".stl") {
            parse_stl_mesh(std::path::Path::new(path))
        } else {
            Err(String::new())
        };
        match result {
            Ok(m) => {
                self.mesh = Some(m);
                self.mesh_label = std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.into());
                self.status = (format!("Mesh loaded: {}", self.mesh_label), false);
            }
            Err(e) => self.status = (format!("Failed to load mesh: {}", if e.is_empty() { path.to_string() } else { e }), true),
        }
    }

    /// detectVoids: true when any void was found.
    fn detect(&mut self) -> bool {
        self.voids.clear();
        self.picked.clear();
        self.selected = 0;
        let Some(source) = self.source.as_ref().filter(|s| !s.atoms.is_empty()) else {
            self.rebuild_groups();
            return false;
        };
        let params = VoidParams { grid_resolution: self.grid, max_voids: if self.step.is_some() { 0 } else { self.max_voids }, min_clearance: self.min_clearance, min_separation: self.min_separation, ..Default::default() };
        let detection = detect_interstitial_voids(source, &params);
        if detection.success {
            self.voids = detection.regions;
        }
        if !self.voids.is_empty() {
            self.picked.push(0);
        }
        self.rebuild_groups();
        let found = !self.voids.is_empty();
        self.status = if found {
            (format!("Detected {} void candidates{}.", self.voids.len(), if self.max_voids > 0 { " (capped)" } else { "" }), false)
        } else {
            ("No voids detected. Increase grid or lower clearance threshold.".into(), true)
        };
        found
    }

    /// rebuildVoidTypeOptions: groups by volume bin, keeping earlier choices.
    fn rebuild_groups(&mut self) {
        let previous: Vec<(i32, bool)> = self.groups.iter().map(|g| (g.bin, g.enabled)).collect();
        let mut bins: std::collections::BTreeMap<i32, (usize, f64)> = std::collections::BTreeMap::new();
        for v in &self.voids {
            let g = bins.entry(v.volume_bin).or_default();
            g.0 += 1;
            g.1 += v.volume as f64;
        }
        self.groups = bins
            .into_iter()
            .enumerate()
            .map(|(i, (bin, (count, total)))| VoidGroup {
                bin,
                count,
                label: format!("Void {} ({count}) [V {:.3} Å³]", i + 1, if count > 0 { total / count as f64 } else { 0.0 }),
                enabled: previous.iter().find(|p| p.0 == bin).is_none_or(|p| p.1),
            })
            .collect();
    }

    fn matches_groups(&self, v: &VoidRegion) -> bool {
        self.groups.is_empty() || self.groups.iter().any(|g| g.enabled && g.bin == v.volume_bin)
    }

    fn inside_mesh(&self, point: Vec3) -> bool {
        let Some(mesh) = self.mesh.as_ref() else { return false };
        if mesh.vertices.is_empty() || mesh.indices.is_empty() || self.mesh_scale <= 0.0 {
            return false;
        }
        let inv = mesh_transform(self.mesh_scale, self.mesh_translation, self.mesh_rotation).inverse();
        let local = (inv * Vec4::new(point.x, point.y, point.z, 1.0)).truncate();
        let dirs = [Vec3::new(1.0, 0.137, 0.319).normalize(), Vec3::new(0.173, 1.0, 0.271).normalize(), Vec3::new(0.259, 0.347, 1.0).normalize()];
        let n = mesh.vertices.len();
        let votes = dirs
            .iter()
            .filter(|&&dir| {
                let hits = mesh
                    .indices
                    .chunks_exact(3)
                    .filter(|t| t.iter().all(|&i| (i as usize) < n) && ray_intersects_triangle(local, dir, mesh.vertices[t[0] as usize], mesh.vertices[t[1] as usize], mesh.vertices[t[2] as usize]))
                    .count();
                hits % 2 == 1
            })
            .count();
        votes >= 2
    }

    /// The candidate voids for placement (before the target count).
    fn pool(&self) -> Vec<usize> {
        match self.placement {
            Placement::Manual => {
                let mut pool = Vec::new();
                for &i in &self.picked {
                    if i < self.voids.len() && !pool.contains(&i) {
                        pool.push(i);
                    }
                }
                pool
            }
            _ => (0..self.voids.len()).filter(|&i| self.matches_groups(&self.voids[i]) && (self.placement != Placement::MeshRegion || self.inside_mesh(self.voids[i].position))).collect(),
        }
    }

    /// applyPlacement: the source with the inserted atoms, and how many.
    fn place(&self, colours: &[[f32; 3]]) -> Option<(Structure, usize)> {
        let source = self.source.as_ref()?;
        let mut pool = self.pool();
        if pool.is_empty() {
            return None;
        }
        let desired = match self.placement {
            Placement::Manual => pool.len(),
            _ if self.by_percent => (pool.len() as f64 * self.percent.clamp(0.0, 100.0) as f64 / 100.0).round() as usize,
            _ => self.count.max(0) as usize,
        }
        .min(pool.len());
        if desired == 0 {
            return None;
        }
        if self.placement != Placement::Manual {
            let mut rng = Mt19937::new(self.seed.max(0) as u32);
            shuffle(&mut pool, &mut rng);
        }
        let mut out = source.clone();
        let cell = af_core::glm::cell_matrices(&out);
        let mut existing: Vec<Vec3> = out.atoms.iter().map(|a| Vec3::new(a.x as f32, a.y as f32, a.z as f32)).collect();
        let z = self.element.clamp(1, 118);
        let c = colours.get(z as usize).copied().unwrap_or([0.8, 0.8, 0.2]);
        let min_insert = self.min_insert.max(0.0);
        let mut added = 0;
        for &i in pool.iter().take(desired) {
            let p = self.voids[i].position;
            let clash = existing.iter().any(|&e| {
                let d = match &cell {
                    Some((m, inv)) => af_core::glm::minimum_image_delta(p - e, true, m, inv),
                    None => p - e,
                };
                af_core::glm::dot(d, d) < min_insert * min_insert
            });
            if clash {
                continue;
            }
            out.atoms.push(AtomSite { atomic_number: z, symbol: af_core::elements::symbol(z).into(), x: p.x as f64, y: p.y as f64, z: p.z as f64, r: c[0], g: c[1], b: c[2], ..Default::default() });
            existing.push(p);
            added += 1;
        }
        Some((out, added))
    }

    fn can_apply(&self) -> bool {
        self.source.is_some() && !self.voids.is_empty() && (self.placement != Placement::MeshRegion || self.mesh.is_some()) && (self.placement != Placement::Manual || !self.picked.is_empty())
    }

    /// Apply Interstitials.
    fn apply(&mut self, colours: &[[f32; 3]]) -> Option<Outcome> {
        match self.place(colours) {
            Some((out, added)) if added > 0 => {
                // The result becomes the source, as in the C++ dialog.
                self.source = Some(out.clone());
                self.preview.fit_pending = true;
                let text = format!("Added {added} interstitial atom(s).");
                self.status = (text.clone(), false);
                Some(Outcome::Edited(out, text))
            }
            _ => {
                self.status = ("No interstitials placed. Check target, type filter, or mesh transform.".into(), true);
                None
            }
        }
    }

    /// Clicking a void selects it (Ctrl+click adds or removes).
    fn pick(&mut self, index: usize, add: bool) {
        self.selected = index;
        if add {
            match self.picked.iter().position(|&i| i == index) {
                None => self.picked.push(index),
                Some(at) if self.picked.len() > 1 => {
                    self.picked.remove(at);
                }
                _ => {}
            }
        } else {
            self.picked = vec![index];
        }
    }

    fn draw_preview(&mut self, ui: &mut egui::Ui, size: egui::Vec2) {
        let (response, painter) = self.preview.canvas(ui, size);
        let rect = response.rect;
        let Some(source) = self.source.as_ref().filter(|s| !s.atoms.is_empty()) else {
            canvas_hint(&painter, rect, if self.source.is_some() { "No atoms to preview." } else { "Load a structure or press 'Use Current Scene'." }, ui);
            return;
        };
        let atoms = preview_atoms(source, 1.0);
        if self.preview.fit_pending {
            let points: Vec<(Vec3, f32)> = atoms.iter().map(|a| (a.0, a.2)).collect();
            self.preview.fit(&points, rect);
            self.preview.fit_pending = false;
        }
        self.preview.draw_cell(&painter, rect, source, ui);
        self.preview.draw_atoms(&painter, rect, &atoms);
        // The mesh region as a translucent wireframe.
        if let Some(mesh) = self.mesh.as_ref().filter(|_| self.placement == Placement::MeshRegion) {
            let t = mesh_transform(self.mesh_scale, self.mesh_translation, self.mesh_rotation);
            let stroke = egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(115, 166, 217, 140));
            for tri in mesh.indices.chunks_exact(3).take(20_000) {
                if tri.iter().any(|&i| i as usize >= mesh.vertices.len()) {
                    continue;
                }
                let p: Vec<egui::Pos2> = tri.iter().map(|&i| self.preview.project((t * mesh.vertices[i as usize].extend(1.0)).truncate(), rect).0).collect();
                painter.add(egui::Shape::closed_line(p, stroke));
            }
        }
        if self.voids.is_empty() {
            return;
        }
        let edge = |k: VoidKind| match k {
            VoidKind::Tetrahedral => crate::theme::srgba(0.2, 0.8, 1.0, 0.95),
            VoidKind::Octahedral => crate::theme::srgba(1.0, 0.68, 0.2, 0.95),
            VoidKind::Irregular => crate::theme::srgba(0.78, 0.4, 1.0, 0.92),
        };
        let fill = |k: VoidKind| match k {
            VoidKind::Tetrahedral => crate::theme::srgba(0.2, 0.8, 1.0, 0.18),
            VoidKind::Octahedral => crate::theme::srgba(1.0, 0.68, 0.2, 0.2),
            VoidKind::Irregular => crate::theme::srgba(0.78, 0.4, 1.0, 0.16),
        };
        let click = response.clicked().then(|| response.interact_pointer_pos()).flatten();
        let ctrl = ui.input(|i| i.modifiers.command);
        let mut best: Option<(f32, usize)> = None;
        if self.selected >= self.voids.len() {
            self.selected = 0;
        }
        for (i, v) in self.voids.iter().enumerate() {
            let (at, _) = self.preview.project(v.position, rect);
            let chosen = self.picked.contains(&i);
            let primary = i == self.selected;
            let (radius, width) = if primary { (5.4, 2.3) } else if chosen { (4.6, 1.9) } else { (3.6, 1.2) };
            painter.circle_stroke(at, radius, egui::Stroke::new(width, edge(v.kind)));
            painter.circle_filled(at, if chosen { 2.0 } else { 1.4 }, edge(v.kind));
            if let Some(m) = click {
                let d2 = (m - at).length_sq();
                if d2 < 64.0 && best.is_none_or(|b| d2 < b.0) {
                    best = Some((d2, i));
                }
            }
        }
        if let Some((_, i)) = best {
            self.pick(i, ctrl);
        }
        if self.show_polyhedra {
            let v = &self.voids[self.selected];
            for tri in &v.face_vertices {
                let p: Vec<egui::Pos2> = tri.iter().map(|&q| self.preview.project(q, rect).0).collect();
                painter.add(egui::Shape::convex_polygon(p.clone(), fill(v.kind), egui::Stroke::NONE));
                painter.add(egui::Shape::closed_line(p, egui::Stroke::new(1.0, edge(v.kind))));
            }
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui, input: &Input) {
        let editing = self.step.is_some();
        ui.strong("Source");
        ui.separator();
        if editing {
            ui.weak("The structure entering the step.");
        } else {
            ui.horizontal(|ui| {
                if ui.button("Use Current Scene").clicked() {
                    self.use_scene(input.source);
                }
                ui.weak("or load a file:");
            });
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut self.structure_path).hint_text("structure file path").desired_width(220.0));
                if ui.add_enabled(!self.structure_path.trim().is_empty(), egui::Button::new("Load")).clicked() {
                    let p = self.structure_path.trim().to_string();
                    self.load_structure(&p);
                }
            });
        }
        match &self.source {
            Some(s) => {
                ui.label(format!("Loaded: {}", self.source_label));
                ui.weak(format!("Atoms: {}", s.atoms.len()));
            }
            None => {
                ui.weak("No source loaded.");
            }
        }
        ui.add_space(6.0);
        ui.strong("Void Detection");
        ui.separator();
        ui.add(egui::Slider::new(&mut self.grid, 8..=30).text("Grid")).on_hover_text(GRID_HINT);
        if !editing {
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.max_voids).range(0..=i32::MAX));
                ui.label("Max Voids");
            })
            .response
            .on_hover_text(GRID_HINT);
            ui.weak("0 = detect all voids");
        }
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.min_clearance).speed(0.02).range(0.05..=5.0).fixed_decimals(2));
            ui.label("Min Clearance (Å)");
        })
        .response
        .on_hover_text(CLEARANCE_HINT);
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.min_separation).speed(0.02).range(0.05..=5.0).fixed_decimals(2));
            ui.label("Min Separation (Å)");
        })
        .response
        .on_hover_text(CLEARANCE_HINT);
        if ui.add_enabled(self.source.is_some(), egui::Button::new("Detect Voids").min_size(egui::vec2(140.0, 0.0))).clicked() {
            self.detect();
        }
        for g in &self.groups {
            ui.label(format!("• {}", g.label));
        }
        ui.add_space(6.0);
        ui.strong("Placement");
        ui.separator();
        if let Some(v) = self.voids.get(self.selected) {
            let kind = match v.kind {
                VoidKind::Tetrahedral => "Tetrahedral",
                VoidKind::Octahedral => "Octahedral",
                VoidKind::Irregular => "Irregular",
            };
            ui.label("Selected Void Info");
            ui.weak("Click to select one, Ctrl+click to add/remove multiple.");
            ui.label(format!("• Index: {} / {}", self.selected + 1, self.voids.len()));
            ui.label(format!("• Manual selection count: {}", self.picked.len()));
            ui.label(format!("• Type: {kind}"));
            ui.label(format!("• Volume: {:.3} Å³", v.volume));
            ui.label(format!("• Clearance: {:.3} Å", v.clearance));
            ui.label(format!("• Coordination: {}", v.coordination));
            ui.label(format!("• Center: ({:.3}, {:.3}, {:.3})", v.position.x, v.position.y, v.position.z));
            ui.weak("Only selected void polygon is shown.");
        }
        if editing {
            self.step_controls(ui, input);
            return;
        }
        ui.label("Void Types For Distribution");
        if self.groups.is_empty() {
            ui.weak("Detect voids to populate void types.");
        } else {
            ui.horizontal(|ui| {
                if ui.button("Select All Types").clicked() {
                    self.groups.iter_mut().for_each(|g| g.enabled = true);
                }
                if ui.button("Clear All Types").clicked() {
                    self.groups.iter_mut().for_each(|g| g.enabled = false);
                }
            });
            ui.weak("Use checkboxes to choose one or multiple void types.");
            for g in self.groups.iter_mut() {
                crate::dialog_widgets::checkbox(ui, &mut g.enabled, g.label.as_str());
            }
            if !self.groups.iter().any(|g| g.enabled) {
                ui.colored_label(egui::Color32::from_rgb(220, 80, 60), "No void type selected. Placement pool will be empty.");
            }
        }
        crate::dialog_widgets::checkbox(ui, &mut self.show_polyhedra, "Show Void Polyhedra in Scene");
        ui.horizontal(|ui| {
            if ui.button("Select Element").clicked() {
                self.picking_element = !self.picking_element;
            }
            ui.label(af_core::elements::symbol(self.element.clamp(1, 118)));
        });
        if self.picking_element && crate::display_settings::periodic_table(ui, &mut self.element, Some(input.colours)) {
            self.picking_element = false;
        }
        egui::ComboBox::from_label("Mode").width(260.0).selected_text(PLACEMENTS.iter().find(|p| p.0 == self.placement).map(|p| p.1).unwrap_or("")).show_ui(ui, |ui| {
            for (mode, label) in PLACEMENTS {
                ui.selectable_value(&mut self.placement, mode, label);
            }
        });
        if self.placement == Placement::Manual {
            ui.weak("Manual mode inserts atoms into the currently selected void centers.");
            ui.weak(format!("Selected centers: {}", self.picked.len()));
        } else {
            egui::ComboBox::from_label("Target").selected_text(if self.by_percent { "By percentage" } else { "By atom count" }).show_ui(ui, |ui| {
                ui.selectable_value(&mut self.by_percent, true, "By percentage");
                ui.selectable_value(&mut self.by_percent, false, "By atom count");
            });
            if self.by_percent {
                ui.add(egui::Slider::new(&mut self.percent, 0.0..=100.0).fixed_decimals(1).suffix("%").text("Percent"));
            } else {
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut self.count).range(0..=i32::MAX));
                    ui.label("Atom Count");
                });
            }
        }
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.seed));
            ui.label("Seed");
        })
        .response
        .on_hover_text(SEED_HINT);
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.min_insert).speed(0.02).range(0.0..=5.0).fixed_decimals(2));
            ui.label("Min Insert Distance (Å)");
        })
        .response
        .on_hover_text(SEED_HINT);
        if self.placement == Placement::MeshRegion {
            ui.add_space(6.0);
            ui.strong("3D Region (OBJ/STL)");
            ui.separator();
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut self.mesh_path).hint_text(".obj or .stl path").desired_width(220.0));
                if ui.add_enabled(!self.mesh_path.trim().is_empty(), egui::Button::new("Load Mesh")).clicked() {
                    let p = self.mesh_path.trim().to_string();
                    self.load_mesh(&p);
                }
            });
            if self.mesh.is_some() {
                ui.label(format!("Loaded: {}", self.mesh_label));
            } else {
                ui.weak("Load an OBJ/STL mesh to define the region.");
            }
            ui.weak("Position it with the numeric controls:");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.mesh_scale).speed(0.01).range(1e-4..=1000.0).fixed_decimals(3));
                ui.label("Scale");
            });
            ui.horizontal(|ui| {
                for c in [&mut self.mesh_translation.x, &mut self.mesh_translation.y, &mut self.mesh_translation.z] {
                    ui.add(egui::DragValue::new(c).speed(0.05).range(-10000.0..=10000.0).fixed_decimals(2));
                }
                ui.label("Translate");
            });
            ui.horizontal(|ui| {
                for c in [&mut self.mesh_rotation.x, &mut self.mesh_rotation.y, &mut self.mesh_rotation.z] {
                    ui.add(egui::DragValue::new(c).speed(0.5).range(-360.0..=360.0).fixed_decimals(1));
                }
                ui.label("Rotate (°)");
            });
            if ui.button("Reset Transform").clicked() {
                self.mesh_translation = Vec3::ZERO;
                self.mesh_rotation = Vec3::ZERO;
                self.mesh_scale = 1.0;
            }
        }
    }

    /// The step-mode placement controls (element, void kind, site order, count, seed).
    fn step_controls(&mut self, ui: &mut egui::Ui, input: &Input) {
        crate::dialog_widgets::checkbox(ui, &mut self.show_polyhedra, "Show Void Polyhedra in Scene");
        ui.horizontal(|ui| {
            if ui.button("Select Element").clicked() {
                self.picking_element = !self.picking_element;
            }
            ui.label(af_core::elements::symbol(self.element.clamp(1, 118)));
        });
        if self.picking_element && crate::display_settings::periodic_table(ui, &mut self.element, Some(input.colours)) {
            self.picking_element = false;
        }
        egui::ComboBox::from_label("Void Kind").width(180.0).selected_text(STEP_KIND_LABELS[self.step_kind.min(3)]).show_ui(ui, |ui| {
            for (k, label) in STEP_KIND_LABELS.iter().enumerate() {
                ui.selectable_value(&mut self.step_kind, k, *label);
            }
        });
        egui::ComboBox::from_label("Choose Sites").width(180.0).selected_text(if self.step_random { "Random" } else { "Largest voids first" }).show_ui(ui, |ui| {
            ui.selectable_value(&mut self.step_random, false, "Largest voids first");
            ui.selectable_value(&mut self.step_random, true, "Random");
        });
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.count).range(0..=i32::MAX));
            ui.label("Atom Count");
        });
        ui.weak("0 = fill every site of the kind");
        ui.add_enabled_ui(self.step_random, |ui| {
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.seed));
                ui.label("Seed");
            })
            .response
            .on_hover_text(SEED_HINT);
        });
    }

    fn set(&mut self, name: &str, value: &str) -> Result<(), String> {
        let v = value.trim();
        let f = || v.parse::<f32>().map_err(|_| format!("{name} wants a number, got {v:?}"));
        let i = || v.parse::<i32>().map_err(|_| format!("{name} wants an integer, got {v:?}"));
        let triple = || -> Result<Vec3, String> {
            let w: Vec<f32> = v.split_whitespace().map(str::parse).collect::<Result<_, _>>().map_err(|_| format!("{name} wants three numbers"))?;
            (w.len() == 3).then(|| Vec3::new(w[0], w[1], w[2])).ok_or(format!("{name} wants three numbers"))
        };
        match name {
            "grid" => self.grid = i()?.clamp(8, 30),
            "max-voids" => self.max_voids = i()?.max(0),
            "clearance" => self.min_clearance = f()?,
            "separation" => self.min_separation = f()?,
            "element" => {
                let z = af_core::elements::atomic_number_any_case(v);
                if z < 1 {
                    return Err(format!("no element {v}"));
                }
                self.element = z;
            }
            "mode" => {
                self.placement = match v {
                    "random" => Placement::Random,
                    "mesh" => Placement::MeshRegion,
                    "manual" => Placement::Manual,
                    _ => return Err(format!("mode is random, mesh or manual, not {v}")),
                }
            }
            "target" => self.by_percent = v == "percent",
            "percent" => {
                self.by_percent = true;
                self.percent = f()?;
            }
            "count" => {
                self.by_percent = false;
                self.count = i()?.max(0);
            }
            "seed" => self.seed = i()?,
            "kind" => self.step_kind = STEP_KINDS.iter().position(|k| *k == v).ok_or(format!("kind is any, tetrahedral, octahedral or irregular, not {v}"))?,
            "order" => {
                self.step_random = match v {
                    "random" => true,
                    "largest" => false,
                    _ => return Err(format!("order is largest or random, not {v}")),
                }
            }
            "min-distance" => self.min_insert = f()?,
            "mesh-scale" => self.mesh_scale = f()?,
            "mesh-translate" => self.mesh_translation = triple()?,
            "mesh-rotate" => self.mesh_rotation = triple()?,
            _ => return Err(format!("the interstitial dialog has no {name}")),
        }
        Ok(())
    }
}

impl StructureDialog for InterstitialDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".obj") || lower.ends_with(".stl") {
            self.load_mesh(path);
        } else {
            self.load_structure(path);
        }
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["insert-interstitials"]
    }

    fn title(&self) -> &'static str {
        "Add Interstitial Atoms"
    }

    fn open(&mut self, _step: &str) {
        // From the menu: the dialog inserts normally again.
        self.step = None;
        self.open = true;
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        self.begin_step(edit);
        true
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        let mut outcome = None;
        let mut open = self.open;
        let mut close = false;
        crate::theme::modal_window("Add Interstitial Atoms").open(&mut open).default_size([1000.0, 680.0]).show(ctx, |ui| {
            if let Some(edit) = &self.step {
                ui.colored_label(egui::Color32::from_rgb(110, 170, 240), format!("Editing pipeline step {}: {}", edit.index + 1, edit.step));
                ui.weak("The step fills voids of the structure entering it; mesh regions, manual picks and volume groups are not part of the step.");
            }
            let fit = crate::dialog_layout::fit_height(ui, "interstitials", 48.0, 300.0);
            let height = fit.height;
            ui.horizontal_top(|ui| {
                let width = (ui.available_width() * 0.5).clamp(260.0, 520.0);
                self.draw_preview(ui, egui::vec2(width, height));
                egui::ScrollArea::vertical().id_salt("interstitial controls").max_height(height).show(ui, |ui| {
                    ui.vertical(|ui| {
                        ui.set_min_width(340.0);
                        self.controls(ui, input);
                    });
                });
            });
            ui.separator();
            {
                let kind = if self.status.1 { crate::dialog_layout::Note::Error } else { crate::dialog_layout::Note::Success };
                let (label, can) = if self.step.is_some() { ("Update Step", true) } else { ("Apply Interstitials", self.can_apply()) };
                match crate::dialog_layout::action_bar(ui, &self.status.0.clone(), kind, label, Some("Close"), can) {
                    1 if self.step.is_some() => outcome = self.update_step(),
                    1 => outcome = self.apply(input.colours),
                    2 => close = true,
                    _ => {}
                }
            }
            fit.finish(ui);
        });
        if !(open && !close) && self.open {
            // Closed (Close, X): the menu path inserts normally again.
            self.open = false;
            self.step = None;
        }
        outcome
    }

    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let number = |i: usize| args.split_whitespace().nth(i).and_then(|w| w.parse::<f64>().ok());
        match command {
            "use-scene" => {
                self.use_scene(input.source);
                if self.source.is_none() {
                    return Err(self.status.0.clone());
                }
            }
            "load" => {
                self.load_structure(args.trim());
                if self.status.1 {
                    return Err(self.status.0.clone());
                }
            }
            "load-mesh" => {
                self.load_mesh(args.trim());
                if self.status.1 {
                    return Err(self.status.0.clone());
                }
            }
            "set" => {
                let (name, value) = args.trim().split_once(char::is_whitespace).unwrap_or((args.trim(), ""));
                self.set(name, value)?;
            }
            "detect" => {
                if !self.detect() {
                    return Err(self.status.0.clone());
                }
            }
            "types" => {
                // types all | none | 1,3 (1-based volume groups)
                match args.trim() {
                    "all" => self.groups.iter_mut().for_each(|g| g.enabled = true),
                    "none" => self.groups.iter_mut().for_each(|g| g.enabled = false),
                    list => {
                        let wanted: Vec<usize> = list.split([',', ' ']).filter(|w| !w.is_empty()).map(|w| w.parse::<usize>().map_err(|_| format!("bad group {w}"))).collect::<Result<_, _>>()?;
                        for (i, g) in self.groups.iter_mut().enumerate() {
                            g.enabled = wanted.contains(&(i + 1));
                        }
                    }
                }
            }
            "select" => {
                // select I [J ...]: picks voids (1-based), as clicks then Ctrl+clicks.
                for (n, w) in args.split_whitespace().enumerate() {
                    let i = w.parse::<usize>().ok().filter(|&i| i >= 1 && i <= self.voids.len()).ok_or(format!("no void {w}"))?;
                    self.pick(i - 1, n > 0);
                }
            }
            "update-step" => return self.update_step().map(Some).ok_or("no pipeline step is being edited".into()),
            "apply" => {
                if self.step.is_some() {
                    return Err("a pipeline step is being edited: use update-step".into());
                }
                if !self.can_apply() {
                    return Err("nothing to apply: load a source, detect voids (and load a mesh or pick voids)".into());
                }
                return self.apply(input.colours).map(Some).ok_or(self.status.0.clone());
            }
            "expect-voids" => {
                let n = number(0).ok_or("expect-voids N")? as usize;
                if self.voids.len() != n {
                    return Err(format!("{} voids detected", self.voids.len()));
                }
            }
            "expect-groups" => {
                let n = number(0).ok_or("expect-groups N")? as usize;
                if self.groups.len() != n {
                    return Err(format!("{} volume groups: {:?}", self.groups.len(), self.groups.iter().map(|g| (&g.label, g.count)).collect::<Vec<_>>()));
                }
            }
            "expect-pool" => {
                let n = number(0).ok_or("expect-pool N")? as usize;
                let pool = self.pool().len();
                if pool != n {
                    return Err(format!("{pool} voids in the placement pool"));
                }
            }
            _ => return Err(format!("unknown interstitial command {command}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_parameters_round_trip() {
        let t = af_pipeline::pipeline::find_modifier_type("insert-interstitials").unwrap();
        let input = Input { source: &Structure::default(), colours: &[] };
        // Missing parameters take the step's defaults.
        let mut d = InterstitialDialog::default();
        assert!(d.edit_step(StepEdit { index: 2, step: "insert-interstitials".into(), ..Default::default() }));
        assert_eq!((d.element, d.step_kind, d.count, d.step_random, d.seed, d.grid), (1, 0, 1, false, 1, 14));
        assert!(!d.by_percent);
        let Some(Outcome::UpdateStep(out)) = d.script("update-step", "", &input).unwrap() else { panic!("no update") };
        assert_eq!(out.index, 2);
        let value = |name: &str| af_pipeline::pipeline::parameter(&out.parameters, &t, name).unwrap();
        assert_eq!(value("element").string().unwrap(), "H");
        assert_eq!(value("order").string().unwrap(), "largest");
        assert_eq!(value("clearance").number().unwrap(), 0.6f32 as f64);
        assert!(d.step.is_none() && !d.open);
        // Changed values come back, and read back the same.
        let mut edit = StepEdit { step: "insert-interstitials".into(), ..Default::default() };
        edit.parameters.set("element", "Li");
        edit.parameters.set("kind", "octahedral");
        edit.parameters.set("count", 3);
        edit.parameters.set("order", "random");
        edit.parameters.set("seed", 7);
        d.edit_step(edit);
        assert_eq!((d.element, d.step_kind, d.count, d.step_random, d.seed), (3, 2, 3, true, 7));
        d.script("set", "kind tetrahedral", &input).unwrap();
        d.script("set", "separation 0.75", &input).unwrap();
        let Some(Outcome::UpdateStep(out)) = d.script("update-step", "", &input).unwrap() else { panic!("no update") };
        let mut again = InterstitialDialog::default();
        again.edit_step(out.clone());
        assert_eq!((again.element, again.step_kind, again.count, again.step_random, again.seed, again.min_separation), (3, 1, 3, true, 7, 0.75));
        // Opening from the menu ends the step edit.
        again.open("insert-interstitials");
        assert!(again.script("update-step", "", &input).is_err());
    }

    #[test]
    fn ray_counts_cube_crossings() {
        // A unit cube as 12 triangles: a ray from inside crosses once.
        let v: Vec<Vec3> = (0..8).map(|i| Vec3::new((i & 1) as f32, ((i >> 1) & 1) as f32, ((i >> 2) & 1) as f32)).collect();
        let faces = [[0, 1, 3], [0, 3, 2], [4, 6, 7], [4, 7, 5], [0, 4, 5], [0, 5, 1], [2, 3, 7], [2, 7, 6], [0, 2, 6], [0, 6, 4], [1, 5, 7], [1, 7, 3]];
        let mut d = InterstitialDialog { mesh: Some(Mesh { vertices: v, indices: faces.iter().flatten().map(|&i| i as u32).collect() }), ..Default::default() };
        assert!(d.inside_mesh(Vec3::splat(0.5)));
        assert!(!d.inside_mesh(Vec3::splat(1.5)));
        d.mesh_translation = Vec3::splat(1.0);
        assert!(d.inside_mesh(Vec3::splat(1.5)));
        d.mesh_scale = 0.2;
        assert!(!d.inside_mesh(Vec3::splat(1.5)));
    }
}
