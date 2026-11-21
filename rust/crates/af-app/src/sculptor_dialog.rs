//! Edit > Cell Sculptor (the C++ ui/CellSculptorDialog): a supercell of the
//! source cut by (hkl) slabs, each either between two distances along its
//! normal or (with a unit cell) a whole number of structural repeats from a
//! start plane; three non-parallel slabs close the region and define the new
//! cell. Two previews: the supercell with the slab faces (atoms outside the
//! region faded), and the cut result.

use af_algorithms::cell_sculptor::{
    csc_apply_slabs, csc_atom_layer_spacing, csc_build_supercell, csc_centroid, csc_clip_halfspace, csc_d1, csc_d2, csc_is_bounded, csc_normal, csc_slab_ref, csc_structural_dhkl, CellSlabPlane, SLAB_TOL,
};
use af_core::glm;
use af_core::structure::Structure;
use eframe::egui;
use glam::{Mat3, Vec3};

use crate::interstitial_dialog::{canvas_hint, preview_atoms, OrbitPreview};
use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

const SLAB_COLOURS: [[f32; 3]; 5] = [[0.40, 0.70, 1.00], [1.00, 0.50, 0.30], [0.40, 1.00, 0.50], [1.00, 0.90, 0.30], [0.90, 0.40, 1.00]];
const PERIODIC_HINT: &str = "Periodic mode uses a starting plane and an integer number of repeat intervals. Start can be negative; N is positive. Inspect the displayed spacing and retained sites because one crystallographic interval need not equal one chemically identical atomic layer.";

pub struct SculptorDialog {
    open: bool,
    source: Option<Structure>,
    source_name: String,
    /// The file the source was loaded from (empty: scene or step input).
    source_path: String,
    /// The path typed into the Load field.
    path_text: String,
    /// The pipeline step being edited (Edit in Dialog).
    step: Option<StepEdit>,
    nx: i32,
    ny: i32,
    nz: i32,
    slabs: Vec<CellSlabPlane>,
    new_h: i32,
    new_k: i32,
    new_l: i32,
    new_d1: f32,
    new_d2: f32,
    ghost_alpha: f32,
    status: String,
    /// The tiled source, the slab faces, the inside mask and the cut, rebuilt when dirty.
    supercell: Structure,
    faces: Vec<(Vec<Vec3>, [f32; 3])>,
    inside: Vec<bool>,
    result: Structure,
    dirty: bool,
    source_view: OrbitPreview,
    result_view: OrbitPreview,
}

impl Default for SculptorDialog {
    fn default() -> Self {
        Self {
            open: false,
            source: None,
            source_name: String::new(),
            source_path: String::new(),
            path_text: String::new(),
            step: None,
            nx: 1,
            ny: 1,
            nz: 1,
            slabs: Vec::new(),
            new_h: 1,
            new_k: 0,
            new_l: 0,
            new_d1: -5.0,
            new_d2: 5.0,
            ghost_alpha: 0.28,
            status: String::new(),
            supercell: Structure::default(),
            faces: Vec::new(),
            inside: Vec::new(),
            result: Structure::default(),
            dirty: true,
            source_view: OrbitPreview::default(),
            result_view: OrbitPreview::default(),
        }
    }
}

/// The slab faces (rebuildFaces): a parallelepiped for three slabs, else each
/// face clipped by the other slabs.
fn slab_faces(source: &Structure, supercell: &Structure, slabs: &[CellSlabPlane]) -> Vec<(Vec<Vec3>, [f32; 3])> {
    let mut faces = Vec::new();
    if slabs.is_empty() {
        return faces;
    }
    let centroid = csc_centroid(supercell);
    let data: Vec<(Vec3, Vec3, f32, f32)> = slabs.iter().map(|s| (csc_normal(s, source), csc_slab_ref(s, source, centroid), csc_d1(s, source), csc_d2(s, source))).collect();
    if slabs.len() == 3 {
        let a_mat = Mat3::from_cols(data[0].0, data[1].0, data[2].0).transpose();
        if glm::determinant(&a_mat).abs() > 1e-6 {
            let rhs = Vec3::new(data[0].2 + glm::dot(data[0].1, data[0].0), data[1].2 + glm::dot(data[1].1, data[1].0), data[2].2 + glm::dot(data[2].1, data[2].0));
            let p0 = glm::mul(&glm::inverse(&a_mat), rhs);
            let a = data[0].0 * (data[0].3 - data[0].2);
            let b = data[1].0 * (data[1].3 - data[1].2);
            let c = data[2].0 * (data[2].3 - data[2].2);
            let colour = |i: usize| slabs[i].color;
            faces.push((vec![p0, p0 + b, p0 + b + c, p0 + c], colour(0)));
            faces.push((vec![p0 + a, p0 + a + b, p0 + a + b + c, p0 + a + c], colour(0)));
            faces.push((vec![p0, p0 + c, p0 + a + c, p0 + a], colour(1)));
            faces.push((vec![p0 + b, p0 + a + b, p0 + a + b + c, p0 + b + c], colour(1)));
            faces.push((vec![p0, p0 + a, p0 + a + b, p0 + b], colour(2)));
            faces.push((vec![p0 + c, p0 + a + c, p0 + a + b + c, p0 + b + c], colour(2)));
            return faces;
        }
    }
    let max_extent = supercell.atoms.iter().map(|a| glm::length(Vec3::new(a.x as f32, a.y as f32, a.z as f32) - centroid)).fold(5.0f32, f32::max) * 1.4 + 3.0;
    for (si, &(n, reference, d1, d2)) in data.iter().enumerate() {
        for (face_d, lower) in [(d1, true), (d2, false)] {
            let up = if glm::dot(n, Vec3::Y).abs() < 0.9 { Vec3::Y } else { Vec3::X };
            let u = glm::normalize(up - glm::dot(up, n) * n);
            let v = glm::cross(n, u);
            let origin = reference + n * face_d;
            let mut poly = vec![origin + (u + v) * max_extent, origin + (u - v) * max_extent, origin + (-u - v) * max_extent, origin + (-u + v) * max_extent];
            for (sj, &(nj, rj, d1j, d2j)) in data.iter().enumerate() {
                if poly.is_empty() {
                    break;
                }
                if sj == si {
                    poly = if lower { csc_clip_halfspace(&poly, -n, -d2 - glm::dot(reference, n)) } else { csc_clip_halfspace(&poly, n, d1 + glm::dot(reference, n)) };
                } else {
                    poly = csc_clip_halfspace(&poly, nj, d1j + glm::dot(rj, nj));
                    if !poly.is_empty() {
                        poly = csc_clip_halfspace(&poly, -nj, -d2j - glm::dot(rj, nj));
                    }
                }
            }
            if poly.len() >= 3 {
                faces.push((poly, slabs[si].color));
            }
        }
    }
    faces
}

impl SculptorDialog {
    fn has_periodic(&self) -> bool {
        self.source.as_ref().is_some_and(|s| s.has_unit_cell)
    }

    fn set_source(&mut self, s: Structure, name: &str, path: &str) {
        self.source = Some(s);
        self.source_name = name.into();
        self.source_path = path.into();
        self.dirty = true;
        self.source_view.fit_pending = true;
        self.result_view.fit_pending = true;
    }

    fn use_scene(&mut self, scene: &Structure) {
        if scene.atoms.is_empty() {
            self.status = "Scene is empty.".into();
            return;
        }
        self.set_source(scene.clone(), "scene", "");
        self.status = format!("Using scene ({} atoms)", scene.atoms.len());
    }

    fn load(&mut self, path: &str) -> bool {
        match af_io::loader::load_structure_from_file(path) {
            Ok(s) => {
                let name = std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.into());
                self.status = format!("Loaded {name} ({} atoms)", s.atoms.len());
                self.set_source(s, &name, path);
                true
            }
            Err(e) => {
                self.status = format!("Failed to load: {}", if e.is_empty() { path.to_string() } else { e });
                false
            }
        }
    }

    fn clear(&mut self) {
        *self = Self { open: self.open, step: self.step.take(), slabs: std::mem::take(&mut self.slabs), nx: self.nx, ny: self.ny, nz: self.nz, ..Default::default() };
    }

    /// editStep: the source from --input (else the structure entering the
    /// step), the supercell and the slabs from the step's options.
    fn begin_step(&mut self, edit: StepEdit) {
        self.status.clear();
        let options = match af_pipeline::options::Reader::new(&edit.options()) {
            Ok(r) => r,
            Err(e) => {
                self.status = format!("Options not read: {e}");
                af_pipeline::options::Reader::default()
            }
        };
        let input_file = options.text("--input", "");
        if input_file.is_empty() || !self.load(&input_file) {
            if !edit.input.atoms.is_empty() {
                self.set_source(edit.input.clone(), "step input", "");
            } else {
                self.source = None;
                self.source_path.clear();
                self.dirty = true;
            }
        }
        self.nx = options.integer("--nx", 1).max(1);
        self.ny = options.integer("--ny", 1).max(1);
        self.nz = options.integer("--nz", 1).max(1);
        // --slabs "h k l lower upper periodic;..." (periodic bounds are start plane and count).
        self.slabs.clear();
        for row in options.text("--slabs", "").split(';') {
            let w: Vec<&str> = row.split_whitespace().collect();
            let ints = |i: usize| w.get(i).and_then(|v| v.parse::<i32>().ok());
            let floats = |i: usize| w.get(i).and_then(|v| v.parse::<f64>().ok());
            let (Some(h), Some(k), Some(l), Some(lower), Some(upper), Some(periodic)) = (ints(0), ints(1), ints(2), floats(3), floats(4), ints(5)) else { continue };
            let mut slab = CellSlabPlane { h, k, l, use_periodic: periodic != 0, ..Default::default() };
            if slab.use_periodic {
                slab.start_plane = lower.round() as i32;
                slab.n_periods = (upper.round() as i32).max(1);
            } else {
                slab.d1 = lower as f32;
                slab.d2 = upper as f32;
            }
            slab.color = SLAB_COLOURS[self.slabs.len() % 5];
            self.slabs.push(slab);
        }
        self.step = Some(edit);
        self.dirty = true;
        self.source_view.fit_pending = true;
        self.result_view.fit_pending = true;
        self.open = true;
    }

    /// stepOptions: the slabs (periodic ones as start and count), the supercell
    /// and a source file that replaces the step input.
    fn step_options(&self) -> String {
        let number = |v: f32| af_core::format::format_g(v as f64, 7);
        let has_periodic = self.has_periodic();
        let slabs: Vec<String> = self
            .slabs
            .iter()
            .map(|slab| {
                let bounds = if slab.use_periodic && has_periodic { format!("{} {} 1", slab.start_plane, slab.n_periods.max(1)) } else { format!("{} {} 0", number(slab.d1), number(slab.d2)) };
                format!("{} {} {} {bounds}", slab.h, slab.k, slab.l)
            })
            .collect();
        let mut out = af_pipeline::options::Writer::new();
        out.add("--slabs", &slabs.join(";"));
        if self.nx != 1 {
            out.add_integer("--nx", self.nx as i64);
        }
        if self.ny != 1 {
            out.add_integer("--ny", self.ny as i64);
        }
        if self.nz != 1 {
            out.add_integer("--nz", self.nz as i64);
        }
        if !self.source_path.is_empty() {
            out.add("--input", &self.source_path);
        }
        out.str().to_string()
    }

    /// Update Step (needs a slab; the region need not be closed).
    fn update_step(&mut self) -> Option<Outcome> {
        if self.slabs.is_empty() {
            return None;
        }
        let edit = self.step.take()?;
        self.open = false;
        Some(Outcome::UpdateStep(edit.with_options(&self.step_options())))
    }

    fn add_slab(&mut self, h: i32, k: i32, l: i32, d1: f32, d2: f32) {
        let mut slab = CellSlabPlane { h, k, l, d1, d2: if d2 < d1 + 0.05 { d1 + 0.05 } else { d2 }, use_periodic: self.has_periodic(), start_plane: 0, n_periods: 1, ..Default::default() };
        slab.color = SLAB_COLOURS[self.slabs.len() % 5];
        self.slabs.push(slab);
        self.dirty = true;
    }

    /// The supercell, slab faces, inside mask and cut (rebuildFaces, rebuildResult).
    fn rebuild(&mut self) {
        if !self.dirty {
            return;
        }
        self.dirty = false;
        let Some(source) = self.source.as_ref() else {
            self.supercell = Structure::default();
            self.faces.clear();
            self.inside.clear();
            self.result = Structure::default();
            return;
        };
        self.supercell = csc_build_supercell(source, self.nx, self.ny, self.nz);
        self.faces = slab_faces(source, &self.supercell, &self.slabs);
        let centroid = csc_centroid(&self.supercell);
        let data: Vec<(Vec3, Vec3, f32, f32, bool)> = self.slabs.iter().map(|s| (csc_normal(s, source), csc_slab_ref(s, source, centroid), csc_d1(s, source), csc_d2(s, source), s.use_periodic && source.has_unit_cell)).collect();
        self.inside = self
            .supercell
            .atoms
            .iter()
            .map(|a| {
                let p = Vec3::new(a.x as f32, a.y as f32, a.z as f32);
                data.iter().all(|&(n, r, d1, d2, periodic)| {
                    let proj = glm::dot(p - r, n);
                    proj >= d1 - SLAB_TOL && if periodic { proj < d2 - SLAB_TOL } else { proj <= d2 + SLAB_TOL }
                })
            })
            .collect();
        self.result = if self.slabs.is_empty() { Structure::default() } else { csc_apply_slabs(&self.supercell, &self.slabs, source) };
        self.result_view.fit_pending = true;
    }

    fn bounded(&self) -> bool {
        self.source.as_ref().is_some_and(|s| csc_is_bounded(&self.slabs, s))
    }

    fn generate(&mut self) -> Option<Outcome> {
        self.rebuild();
        if !self.bounded() || self.result.atoms.is_empty() {
            return None;
        }
        let n = self.result.atoms.len();
        Some(Outcome::Edited(self.result.clone(), format!("Cell sculptor: {n} atoms")))
    }

    fn source_preview(&mut self, ui: &mut egui::Ui, size: egui::Vec2) {
        let (response, painter) = self.source_view.canvas(ui, size);
        let rect = response.rect;
        if self.source.is_none() {
            canvas_hint(&painter, rect, "Load a structure: 'Use Scene' above", ui);
            return;
        }
        let mut atoms = preview_atoms(&self.supercell, 1.0);
        if !self.slabs.is_empty() {
            for (a, inside) in atoms.iter_mut().zip(&self.inside) {
                if !inside {
                    a.3 = self.ghost_alpha;
                }
            }
        }
        if self.source_view.fit_pending {
            let points: Vec<(Vec3, f32)> = atoms.iter().map(|a| (a.0, a.2)).collect();
            self.source_view.fit(&points, rect);
            self.source_view.fit_pending = false;
        }
        self.source_view.draw_atoms(&painter, rect, &atoms);
        for (poly, c) in &self.faces {
            let points: Vec<egui::Pos2> = poly.iter().map(|&p| self.source_view.project(p, rect).0).collect();
            let fill: egui::Color32 = crate::theme::srgba(c[0], c[1], c[2], 0.18).into();
            let edge: egui::Color32 = crate::theme::srgba(c[0], c[1], c[2], 0.9).into();
            painter.add(egui::Shape::convex_polygon(points.clone(), fill, egui::Stroke::NONE));
            painter.add(egui::Shape::closed_line(points, egui::Stroke::new(1.5, edge)));
        }
        painter.text(rect.left_bottom() + egui::vec2(6.0, -6.0), egui::Align2::LEFT_BOTTOM, "Drag=orbit  Scroll=zoom", egui::FontId::proportional(12.0), ui.visuals().weak_text_color());
    }

    fn result_preview(&mut self, ui: &mut egui::Ui, size: egui::Vec2) {
        let (response, painter) = self.result_view.canvas(ui, size);
        let rect = response.rect;
        if self.result.atoms.is_empty() {
            canvas_hint(&painter, rect, if self.source.is_some() { "Add cutting slabs to see the result" } else { "Load a structure first" }, ui);
            return;
        }
        let atoms = preview_atoms(&self.result, 1.0);
        if self.result_view.fit_pending {
            let points: Vec<(Vec3, f32)> = atoms.iter().map(|a| (a.0, a.2)).collect();
            self.result_view.fit(&points, rect);
            self.result_view.fit_pending = false;
        }
        self.result_view.draw_atoms(&painter, rect, &atoms);
        if self.result.has_unit_cell {
            let o = Vec3::new(self.result.cell_offset[0] as f32, self.result.cell_offset[1] as f32, self.result.cell_offset[2] as f32);
            let v = |r: usize| Vec3::new(self.result.cell_vectors[r][0] as f32, self.result.cell_vectors[r][1] as f32, self.result.cell_vectors[r][2] as f32);
            let (a, b, c) = (v(0), v(1), v(2));
            let stroke = egui::Stroke::new(1.2, egui::Color32::from_rgb(230, 200, 60));
            for (p, q) in [(o, o + a), (o, o + b), (o, o + c), (o + a, o + a + b), (o + a, o + a + c), (o + b, o + a + b), (o + b, o + b + c), (o + c, o + a + c), (o + c, o + b + c), (o + a + b, o + a + b + c), (o + a + c, o + a + b + c), (o + b + c, o + a + b + c)] {
                painter.line_segment([self.result_view.project(p, rect).0, self.result_view.project(q, rect).0], stroke);
            }
        }
        painter.text(rect.left_bottom() + egui::vec2(6.0, -6.0), egui::Align2::LEFT_BOTTOM, "Drag=orbit  Scroll=zoom", egui::FontId::proportional(12.0), ui.visuals().weak_text_color());
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        if !self.status.is_empty() {
            ui.label(&self.status);
        }
        ui.add_space(4.0);
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.set_min_width(ui.available_width());
            ui.horizontal(|ui| {
                ui.label("Supercell");
                ui.weak("(tiles +/- around origin)");
            });
            ui.separator();
            ui.horizontal(|ui| {
                let mut changed = false;
                for (label, v) in [("x", &mut self.nx), ("y", &mut self.ny), ("z", &mut self.nz)] {
                    ui.weak(label);
                    changed |= ui.add(egui::DragValue::new(v).speed(0.1).range(1..=30)).changed();
                }
                if changed {
                    self.dirty = true;
                    self.source_view.fit_pending = true;
                }
            });
        });
        ui.add_space(4.0);
        let has_periodic = self.has_periodic();
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.set_min_width(ui.available_width());
            ui.horizontal(|ui| {
                ui.label("Cutting Slabs");
                ui.weak("(atoms kept within [d1,d2] along hkl)");
            });
            ui.separator();
            let mut remove = None;
            for si in 0..self.slabs.len() {
                ui.push_id(si, |ui| {
                    let source = self.source.clone().unwrap_or_default();
                    let slab = &mut self.slabs[si];
                    let mut changed = false;
                    ui.horizontal(|ui| {
                        if ui.color_edit_button_rgb(&mut slab.color).changed() {
                            changed = true;
                        }
                        ui.label(format!("Slab {}", si + 1));
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui.button("Remove").clicked() {
                                remove = Some(si);
                            }
                        });
                    });
                    ui.horizontal(|ui| {
                        for (label, v) in [("h", &mut slab.h), ("k", &mut slab.k), ("l", &mut slab.l)] {
                            ui.weak(label);
                            changed |= ui.add(egui::DragValue::new(v).speed(0.15).range(-12..=12)).changed();
                        }
                    });
                    if !has_periodic {
                        slab.use_periodic = false;
                    }
                    ui.horizontal(|ui| {
                        let r = ui.add_enabled(has_periodic, egui::Checkbox::new(&mut slab.use_periodic, "Periodic")).on_hover_text(PERIODIC_HINT).on_disabled_hover_text("Load a crystal with a unit cell to enable periodic mode.");
                        changed |= r.changed();
                        if slab.use_periodic && has_periodic {
                            let d_rep = csc_structural_dhkl(slab, &source);
                            let d_atom = csc_atom_layer_spacing(slab, &source);
                            if (d_rep - d_atom).abs() > d_rep * 0.01 {
                                ui.weak(format!("d={} Å (layer {} Å)", af_core::format::format_g(d_rep as f64, 6), af_core::format::format_g(d_atom as f64, 6)));
                            } else {
                                ui.weak(format!("d={} Å", af_core::format::format_g(d_rep as f64, 6)));
                            }
                        }
                    });
                    if slab.use_periodic && has_periodic {
                        ui.horizontal(|ui| {
                            changed |= ui.add(egui::DragValue::new(&mut slab.start_plane).speed(0.1).range(-200..=200)).on_hover_text(PERIODIC_HINT).changed();
                            ui.label("Start");
                            changed |= ui.add(egui::DragValue::new(&mut slab.n_periods).speed(0.1).range(1..=200)).on_hover_text(PERIODIC_HINT).changed();
                            ui.label("N");
                        });
                        let (e1, e2) = (csc_d1(slab, &source), csc_d2(slab, &source));
                        ui.weak(format!("d1={e1:.3} Å   d2={e2:.3} Å   t={:.3} Å", e2 - e1));
                    } else {
                        ui.horizontal(|ui| {
                            changed |= ui.add(egui::DragValue::new(&mut slab.d1).speed(0.05).range(-500.0..=500.0).fixed_decimals(2)).changed();
                            ui.label("d1");
                            changed |= ui.add(egui::DragValue::new(&mut slab.d2).speed(0.05).range(-500.0..=500.0).fixed_decimals(2)).changed();
                            ui.label("d2");
                        });
                    }
                    if changed {
                        if !slab.use_periodic && slab.d2 < slab.d1 + 0.05 {
                            slab.d2 = slab.d1 + 0.05;
                        }
                        self.dirty = true;
                    }
                    ui.separator();
                });
            }
            if let Some(i) = remove {
                self.slabs.remove(i);
                self.dirty = true;
            }
            ui.horizontal(|ui| {
                for (label, v) in [("h", &mut self.new_h), ("k", &mut self.new_k), ("l", &mut self.new_l)] {
                    ui.weak(label);
                    ui.add(egui::DragValue::new(v).speed(0.15).range(-12..=12));
                }
            });
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.new_d1).speed(0.1).range(-500.0..=500.0).fixed_decimals(1));
                ui.label("d1");
                ui.add(egui::DragValue::new(&mut self.new_d2).speed(0.1).range(-500.0..=500.0).fixed_decimals(1));
                ui.label("d2");
            });
            if ui.add(egui::Button::new("+ Add Slab").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                if self.new_d2 < self.new_d1 + 0.05 {
                    self.new_d2 = self.new_d1 + 0.05;
                }
                let (h, k, l, d1, d2) = (self.new_h, self.new_k, self.new_l, self.new_d1, self.new_d2);
                self.add_slab(h, k, l, d1, d2);
            }
        });
        ui.add_space(4.0);
        ui.separator();
        ui.label("Outside atom opacity");
        ui.add(egui::Slider::new(&mut self.ghost_alpha, 0.0..=1.0).fixed_decimals(2));
        ui.separator();
        if self.slabs.is_empty() {
            ui.weak("Add >=3 non-parallel slabs to close the region.");
        } else if !self.bounded() {
            ui.colored_label(egui::Color32::from_rgb(230, 170, 60), "Region not closed — add more slabs.");
        } else {
            ui.colored_label(egui::Color32::from_rgb(70, 170, 90), "Closed region.");
            if !self.result.atoms.is_empty() {
                ui.weak(format!("{} atoms in result", self.result.atoms.len()));
            }
        }
    }

    fn slab_mut(&mut self, index: &str) -> Result<&mut CellSlabPlane, String> {
        let n = self.slabs.len();
        index.parse::<usize>().ok().filter(|&i| i >= 1 && i <= n).and_then(|i| self.slabs.get_mut(i - 1)).ok_or(format!("no slab {index} ({n} slabs)"))
    }
}

impl StructureDialog for SculptorDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        self.load(path);
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-sculpt"]
    }

    fn title(&self) -> &'static str {
        "Cell Sculptor"
    }

    fn open(&mut self, _step: &str) {
        // From the menu: the dialog builds normally again.
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
        self.rebuild();
        let mut outcome = None;
        let mut open = self.open;
        let mut close = false;
        crate::theme::modal_window("Cell Sculptor").open(&mut open).default_size([1100.0, 780.0]).show(ctx, |ui| {
            if let Some(edit) = &self.step {
                ui.colored_label(egui::Color32::from_rgb(110, 170, 240), format!("Editing pipeline step {}: {}", edit.index + 1, edit.step));
            }
            let fit = crate::dialog_layout::fit_height(ui, "sculptor", 8.0, 420.0);
            let height = fit.height;
            ui.horizontal_top(|ui| {
                let width = (ui.available_width() * 0.5).clamp(280.0, 580.0);
                ui.vertical(|ui| {
                    ui.set_width(width);
                    ui.horizontal(|ui| {
                        ui.label("Supercell + Cutting Planes");
                        if let Some(edit) = &self.step {
                            // The step works on the structure entering it, or on a file (--input).
                            if ui.button("Use Step Input").clicked() {
                                let s = edit.input.clone();
                                if s.atoms.is_empty() {
                                    self.status = "Scene is empty.".into();
                                } else {
                                    self.status = format!("Using scene ({} atoms)", s.atoms.len());
                                    self.set_source(s, "step input", "");
                                }
                            }
                        } else if ui.button("Use Scene").clicked() {
                            self.use_scene(input.source);
                        }
                        if self.source.is_some() {
                            if ui.button("Clear").clicked() {
                                self.clear();
                            }
                            ui.weak(format!("| {}  {} × {} × {}", self.source_name, self.nx, self.ny, self.nz));
                        }
                    });
                    self.rebuild();
                    let top = height * 0.52 - 30.0;
                    self.source_preview(ui, egui::vec2(width, top.max(120.0)));
                    ui.horizontal(|ui| {
                        ui.label("Cut Result");
                        if !self.result.atoms.is_empty() {
                            ui.weak(format!("| {} atoms", self.result.atoms.len()));
                        }
                    });
                    self.result_preview(ui, egui::vec2(width, (height - top - 60.0).max(120.0)));
                });
                ui.vertical(|ui| {
                    ui.label("Cell Sculptor — Controls");
                    ui.separator();
                    egui::ScrollArea::vertical().id_salt("sculptor controls").max_height(height - 70.0).show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.add(egui::TextEdit::singleline(&mut self.path_text).hint_text("or a structure file path").desired_width(220.0));
                            if ui.add_enabled(!self.path_text.trim().is_empty(), egui::Button::new("Load")).clicked() {
                                let p = self.path_text.trim().to_string();
                                self.load(&p);
                            }
                        });
                        self.controls(ui);
                    });
                    ui.separator();
                    {
                        // The step keeps the atoms between the slabs; the region need not be closed.
                        let (label, can) = if self.step.is_some() { ("Update Step", !self.slabs.is_empty()) } else { ("Generate", self.bounded()) };
                        match crate::dialog_layout::action_bar(ui, "", crate::dialog_layout::Note::Info, label, Some("Close"), can) {
                            1 if self.step.is_some() => outcome = self.update_step(),
                            1 => {
                                outcome = self.generate();
                                if outcome.is_some() {
                                    close = true;
                                }
                            }
                            2 => close = true,
                            _ => {}
                        }
                    }
                });
            });
            fit.finish(ui);
        });
        if !(open && !close) && self.open {
            // Closed (Close, X): the menu path builds normally again.
            self.open = false;
            self.step = None;
        }
        outcome
    }

    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let int = |i: usize| words.get(i).and_then(|w| w.parse::<i32>().ok());
        let float = |i: usize| words.get(i).and_then(|w| w.parse::<f32>().ok());
        match command {
            "use-scene" => {
                self.use_scene(input.source);
                if self.source.is_none() {
                    return Err(self.status.clone());
                }
            }
            "load" => {
                if !self.load(args.trim()) {
                    return Err(self.status.clone());
                }
            }
            "clear" => self.clear(),
            "supercell" => {
                // supercell NX NY NZ
                let (x, y, z) = (int(0), int(1), int(2));
                let (Some(x), Some(y), Some(z)) = (x, y, z) else { return Err("supercell NX NY NZ".into()) };
                (self.nx, self.ny, self.nz) = (x.clamp(1, 30), y.clamp(1, 30), z.clamp(1, 30));
                self.source_view.fit_pending = true;
                self.dirty = true;
            }
            "add" => {
                // add H K L [D1 D2]: periodic (start 0, N 1) when the source has a cell.
                let (Some(h), Some(k), Some(l)) = (int(0), int(1), int(2)) else { return Err("add H K L [D1 D2]".into()) };
                self.add_slab(h, k, l, float(3).unwrap_or(self.new_d1), float(4).unwrap_or(self.new_d2));
            }
            "slab" => {
                // slab I periodic on|off | start S | n N | d1 V | d2 V | hkl H K L
                let (index, what) = (words.first().copied().unwrap_or(""), words.get(1).copied().unwrap_or(""));
                let value = |i: usize| words.get(i).copied().ok_or(format!("slab {index} {what} wants a value"));
                let periodic = self.has_periodic();
                let slab = self.slab_mut(index)?;
                let parse_i = |w: &str| w.parse::<i32>().map_err(|_| format!("{w} is not an integer"));
                let parse_f = |w: &str| w.parse::<f32>().map_err(|_| format!("{w} is not a number"));
                match what {
                    "periodic" => slab.use_periodic = periodic && value(2)? == "on",
                    "start" => slab.start_plane = parse_i(value(2)?)?,
                    "n" => slab.n_periods = parse_i(value(2)?)?.max(1),
                    "d1" => slab.d1 = parse_f(value(2)?)?,
                    "d2" => slab.d2 = parse_f(value(2)?)?,
                    "hkl" => {
                        slab.h = parse_i(value(2)?)?;
                        slab.k = parse_i(value(3)?)?;
                        slab.l = parse_i(value(4)?)?;
                    }
                    _ => return Err(format!("unknown slab setting {what}")),
                }
                if !slab.use_periodic && slab.d2 < slab.d1 + 0.05 {
                    slab.d2 = slab.d1 + 0.05;
                }
                self.dirty = true;
            }
            "remove" => {
                let i = int(0).filter(|&i| i >= 1 && i as usize <= self.slabs.len()).ok_or(format!("no slab {args}"))?;
                self.slabs.remove(i as usize - 1);
                self.dirty = true;
            }
            "update-step" => {
                if self.step.is_none() {
                    return Err("no pipeline step is being edited".into());
                }
                return self.update_step().map(Some).ok_or("Add at least one slab.".into());
            }
            "generate" => {
                if self.step.is_some() {
                    return Err("a pipeline step is being edited: use update-step".into());
                }
                self.rebuild();
                if !self.bounded() {
                    return Err("Region not closed — add more slabs.".into());
                }
                return self.generate().map(Some).ok_or("the cut keeps no atoms".into());
            }
            "expect-result" => {
                self.rebuild();
                let n = int(0).ok_or("expect-result N")? as usize;
                if self.result.atoms.len() != n {
                    return Err(format!("{} atoms in the result", self.result.atoms.len()));
                }
            }
            "expect-bounded" => {
                let want = args.trim() != "no";
                if self.bounded() != want {
                    return Err(format!("region closed: {}", self.bounded()));
                }
            }
            _ => return Err(format!("unknown sculptor command {command}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::structure::AtomSite;

    fn cubic(a: f64) -> Structure {
        let mut s = Structure { has_unit_cell: true, ..Default::default() };
        s.cell_vectors = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
        for p in [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]] {
            s.atoms.push(AtomSite { symbol: "Cu".into(), atomic_number: 29, x: p[0] * a, y: p[1] * a, z: p[2] * a, ..Default::default() });
        }
        s
    }

    #[test]
    fn step_options_round_trip() {
        let input = Input { source: &Structure::default(), colours: &[] };
        let mut edit = StepEdit { index: 1, step: "build-sculpt".into(), input: cubic(3.61), ..Default::default() };
        // The step's default options (af-cli build_modifiers).
        edit.parameters.set("options", "--slabs \"1 0 0 0 1 1;0 1 0 0 1 1;0 0 1 0 1 1\"");
        let mut d = SculptorDialog::default();
        assert!(d.edit_step(edit.clone()));
        assert_eq!(d.slabs.len(), 3);
        assert!(d.slabs.iter().all(|s| s.use_periodic && s.start_plane == 0 && s.n_periods == 1));
        assert_eq!((d.nx, d.ny, d.nz), (1, 1, 1));
        assert_eq!(d.source_name, "step input");
        assert_eq!(d.step_options(), "--slabs \"1 0 0 0 1 1;0 1 0 0 1 1;0 0 1 0 1 1\"");
        // A manual slab and a supercell, written as the C++ stepOptions.
        d.script("supercell", "2 3 1", &input).unwrap();
        d.script("slab", "2 periodic off", &input).unwrap();
        d.script("slab", "2 d1 -1.25", &input).unwrap();
        d.script("slab", "2 d2 0.1", &input).unwrap();
        d.script("slab", "3 n 2", &input).unwrap();
        let Some(Outcome::UpdateStep(out)) = d.script("update-step", "", &input).unwrap() else { panic!("no update") };
        let text = out.options();
        assert_eq!(text, "--slabs \"1 0 0 0 1 1;0 1 0 -1.25 0.1 0;0 0 1 0 2 1\" --nx 2 --ny 3");
        assert_eq!(out.index, 1);
        assert!(d.step.is_none() && !d.open);
        // Read back: the same dialog state.
        let mut again = SculptorDialog::default();
        again.edit_step(out);
        assert_eq!((again.nx, again.ny, again.nz), (2, 3, 1));
        assert!(!again.slabs[1].use_periodic && again.slabs[1].d1 == -1.25 && again.slabs[1].d2 == 0.1);
        assert_eq!(again.slabs[2].n_periods, 2);
        assert_eq!(again.step_options(), text);
        // Opening from the menu ends the step edit.
        again.open("build-sculpt");
        assert!(again.script("update-step", "", &input).is_err());
    }

    #[test]
    fn three_periodic_slabs_close_the_cell() {
        let mut d = SculptorDialog::default();
        d.set_source(cubic(3.61), "fcc", "");
        d.script("supercell", "3 3 3", &Input { source: &Structure::default(), colours: &[] }).unwrap();
        for hkl in ["1 0 0", "0 1 0", "0 0 1"] {
            d.script("add", hkl, &Input { source: &Structure::default(), colours: &[] }).unwrap();
        }
        assert!(d.bounded());
        d.rebuild();
        // One conventional cell of fcc: 4 atoms, and 6 faces drawn.
        assert_eq!(d.result.atoms.len(), 4);
        assert_eq!(d.faces.len(), 6);
        assert_eq!(d.inside.iter().filter(|&&i| i).count(), 4);
    }
}
