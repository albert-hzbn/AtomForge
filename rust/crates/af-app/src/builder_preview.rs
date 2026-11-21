//! The small orbitable previews of the builder dialogs (the C++ preview FBOs
//! of the nanocrystal and polycrystal builders): a reference structure drawn
//! as shaded discs with its cell, and a convex polyhedron (the Wulff shape)
//! drawn as lit faces, both with the painter's algorithm.

use af_core::structure::Structure;
use eframe::egui;
use glam::{Mat3, Vec3};

/// Left-drag orbits, the wheel zooms (the C++ yaw/pitch/distance camera).
#[derive(Clone, Copy)]
pub struct Orbit {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
}

impl Default for Orbit {
    fn default() -> Self {
        Self { yaw: 45.0, pitch: 35.0, zoom: 1.0 }
    }
}

impl Orbit {
    fn rotation(&self) -> Mat3 {
        Mat3::from_rotation_x(self.pitch.to_radians()) * Mat3::from_rotation_y(-self.yaw.to_radians())
    }

    /// Drag and wheel over the preview's response.
    fn interact(&mut self, ui: &egui::Ui, response: &egui::Response) {
        if response.dragged() {
            let d = response.drag_delta();
            self.yaw -= d.x * 0.5;
            self.pitch = (self.pitch + d.y * 0.5).clamp(-89.0, 89.0);
        }
        if response.hovered() {
            let wheel = ui.input(|i| i.smooth_scroll_delta.y);
            if wheel != 0.0 {
                self.zoom = (self.zoom * (1.0 + wheel * 0.002)).clamp(0.1, 20.0);
            }
        }
    }
}

fn frame(ui: &mut egui::Ui, size: egui::Vec2) -> (egui::Response, egui::Painter) {
    let (response, painter) = ui.allocate_painter(size, egui::Sense::click_and_drag());
    painter.rect_stroke(response.rect, 2.0, ui.visuals().widgets.noninteractive.bg_stroke, egui::StrokeKind::Inside);
    (response, painter)
}

/// Centre and radius of a point cloud.
fn bounds(points: impl Iterator<Item = Vec3> + Clone) -> (Vec3, f32) {
    let (mut lo, mut hi) = (Vec3::splat(f32::MAX), Vec3::splat(f32::MIN));
    for p in points.clone() {
        lo = lo.min(p);
        hi = hi.max(p);
    }
    if lo.x > hi.x {
        return (Vec3::ZERO, 1.0);
    }
    let centre = (lo + hi) * 0.5;
    let radius = points.map(|p| (p - centre).length()).fold(0.0f32, f32::max).max(1.0);
    (centre, radius)
}

/// A reference structure, drawn by the structure renderer with its periodic
/// images, bonds and cell (crate::preview3d; the preview keeps its own camera).
pub fn structure(ui: &mut egui::Ui, size: egui::Vec2, s: &Structure, _colours: &[[f32; 3]], _orbit: &mut Orbit, empty_hint: &[&str]) -> egui::Response {
    crate::preview3d::show(ui, size, s, empty_hint)
}

/// One face of a convex polyhedron: its vertices (in order), outward normal and colour.
pub struct Face {
    pub vertices: Vec<Vec3>,
    pub normal: Vec3,
    pub colour: [f32; 3],
}

/// A convex polyhedron: lit front faces with dark edges.
pub fn polyhedron(ui: &mut egui::Ui, size: egui::Vec2, faces: &[Face], orbit: &mut Orbit) -> egui::Response {
    let (response, painter) = frame(ui, size);
    orbit.interact(ui, &response);
    let rect = response.rect;
    let painter = painter.with_clip_rect(rect.shrink(2.0));
    let (centre, radius) = bounds(faces.iter().flat_map(|f| f.vertices.iter().copied()).collect::<Vec<_>>().into_iter());
    let r = orbit.rotation();
    let scale = rect.width().min(rect.height()) * 0.45 / radius * orbit.zoom;
    let light = Vec3::new(-0.4, 0.6, 1.0).normalize();
    let mut drawn: Vec<(f32, Vec<egui::Pos2>, egui::Color32)> = Vec::new();
    for f in faces.iter().filter(|f| f.vertices.len() >= 3) {
        let n = r * f.normal;
        if n.z <= 0.0 {
            continue;
        }
        let rotated: Vec<Vec3> = f.vertices.iter().map(|p| r * (*p - centre)).collect();
        let depth = rotated.iter().map(|q| q.z).sum::<f32>() / rotated.len() as f32;
        let points = rotated.iter().map(|q| egui::pos2(rect.center().x + q.x * scale, rect.center().y - q.y * scale)).collect();
        let diffuse = n.normalize_or_zero().dot(light).max(0.34);
        let k = 0.30 + 0.70 * diffuse;
        let c = f.colour;
        drawn.push((depth, points, crate::theme::srgb((c[0] * k).min(1.0), (c[1] * k).min(1.0), (c[2] * k).min(1.0)).into()));
    }
    af_core::cxx::sort_by(&mut drawn, |a, b| a.0 < b.0);
    for (_, points, fill) in drawn {
        painter.add(egui::Shape::convex_polygon(points.clone(), fill, egui::Stroke::NONE));
        painter.add(egui::Shape::closed_line(points, egui::Stroke::new(1.2, egui::Color32::from_gray(25))));
    }
    response
}

/// The structure files a reference can be loaded from (as File > Open).
const STRUCTURE_EXTENSIONS: &[&str] = &["cif", "xyz", "extxyz", "vasp", "poscar", "contcar", "pdb", "mol", "sdf", "mol2", "gjf", "xsf", "pwi", "in", "lmp", "data", "cell", "gen"];

/// A builder's reference crystal: loaded from a file (Load Reference) or
/// taken from the structure in the view.
pub struct Reference {
    pub structure: Structure,
    /// The file's name, or "the structure in the view".
    pub name: String,
    /// "Loaded: N atoms" or why a load failed.
    pub status: String,
    /// The loaded file's path (written as --input by an edited step).
    pub path: String,
    /// The reference is the pipeline step's input (m_referenceFromStep).
    pub from_step: bool,
    browser: crate::file_browser::FileBrowser,
    browser_title: &'static str,
}

impl Reference {
    pub fn new(browser_title: &'static str) -> Self {
        Self { structure: Structure::default(), name: String::new(), status: String::new(), path: String::new(), from_step: false, browser: Default::default(), browser_title }
    }

    pub fn clear(&mut self) {
        self.structure = Structure::default();
        self.name.clear();
        self.status.clear();
        self.path.clear();
        self.from_step = false;
    }

    /// An edited step's reference (editStep): the file its --input names
    /// (kept even when it cannot be loaded now), else the step's input.
    pub fn from_step_edit(&mut self, edit: &crate::structure_dialog::StepEdit, needs_cell: bool, input_name: &str) {
        self.clear();
        let r = edit.reader();
        if r.has("--input") {
            let path = r.text("--input", "");
            let _ = self.load(&path, needs_cell);
            self.path = path;
        } else if edit.use_input() && !edit.input.atoms.is_empty() {
            self.structure = edit.input.clone();
            self.name = input_name.to_string();
            self.from_step = true;
            self.status = format!("Step input: {} atoms", self.structure.atoms.len());
        }
    }

    /// Loads a file; `needs_cell` refuses structures without a unit cell.
    pub fn load(&mut self, path: &str, needs_cell: bool) -> Result<(), String> {
        let loaded = match af_io::loader::load_structure_from_file(path) {
            Ok(s) => s,
            Err(e) => {
                self.status = format!("Load failed: {e}");
                return Err(self.status.clone());
            }
        };
        if needs_cell && !loaded.has_unit_cell {
            self.status = "Error: structure has no unit cell (required for polycrystal).".into();
            return Err(self.status.clone());
        }
        self.name = std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string());
        self.structure = loaded;
        self.path = path.to_string();
        self.from_step = false;
        self.status = format!("Loaded: {} atoms", self.structure.atoms.len());
        Ok(())
    }

    /// Takes the structure in the view.
    pub fn take_view(&mut self, source: &Structure, needs_cell: bool) -> Result<(), String> {
        if source.atoms.is_empty() {
            return Err("there is no structure in the view".into());
        }
        if needs_cell && !source.has_unit_cell {
            self.status = "Error: structure has no unit cell (required for polycrystal).".into();
            return Err(self.status.clone());
        }
        self.structure = source.clone();
        self.name = "the structure in the view".into();
        self.path.clear();
        self.from_step = false;
        self.status = format!("Loaded: {} atoms", self.structure.atoms.len());
        Ok(())
    }

    /// Load Reference and Use Structure in View.
    pub fn buttons(&mut self, ui: &mut egui::Ui, source: &Structure, needs_cell: bool) {
        ui.horizontal(|ui| {
            if ui.button("Load Reference").clicked() {
                self.browser.show_for(crate::file_browser::Purpose::Open, self.browser_title, "", STRUCTURE_EXTENSIONS);
            }
            if ui.add_enabled(!source.atoms.is_empty(), egui::Button::new("Use Structure in View")).on_hover_text("Takes the structure in the main view as the reference.").clicked() {
                let _ = self.take_view(source, needs_cell);
            }
        });
    }

    /// The file browser window (call every frame while the dialog is open).
    pub fn browser(&mut self, ctx: &egui::Context, needs_cell: bool) {
        if let Some(path) = self.browser.show(ctx) {
            let _ = self.load(&path.to_string_lossy(), needs_cell);
        }
    }
}

/// The lengths of a cell's vectors.
pub fn cell_lengths(s: &Structure) -> [f32; 3] {
    [0, 1, 2].map(|i| Vec3::new(s.cell_vectors[i][0] as f32, s.cell_vectors[i][1] as f32, s.cell_vectors[i][2] as f32).length())
}

#[cfg(test)]
mod golden {
    use crate::structure_dialog::{StepEdit, StructureDialog};

    /// editStep then stepOptions of the C++ Nanocrystal, Polycrystal and
    /// Amorphous dialogs (rust/tools/golden/golden_step_builders.cpp).
    #[test]
    fn builder_step_options_match_cpp() {
        let table = include_str!("../../../tests/golden/step_builders.tsv");
        let mut n = 0;
        for line in table.lines().filter(|l| !l.is_empty()) {
            let f: Vec<&str> = line.split('\t').collect();
            let (step, options, expected) = (f[0], f[1], f[2]);
            let mut parameters = af_core::json::Json::object();
            parameters.set("options", options);
            let edit = StepEdit { index: 0, step: step.into(), parameters, input: Default::default(), ..Default::default() };
            let got = match step {
                "build-nano" => {
                    let mut d = crate::nano_dialog::NanoDialog::default();
                    d.edit_step(edit);
                    d.step_options()
                }
                "build-poly" => {
                    let mut d = crate::poly_dialog::PolyDialog::default();
                    d.edit_step(edit);
                    d.step_options()
                }
                _ => {
                    let mut d = crate::amorphous_dialog::AmorphousDialog::default();
                    d.edit_step(edit);
                    d.step_options()
                }
            };
            assert_eq!(got, expected, "{step} {options}");
            n += 1;
        }
        assert_eq!(n, 26);
    }
}
