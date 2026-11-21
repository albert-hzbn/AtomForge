//! Build > Nanocrystal (the C++ ui/NanoCrystalBuilderDialog): a reference
//! crystal carved by a shape (sphere, ellipsoid, box, cylinder, octahedra,
//! cuboctahedron) or by a Wulff construction from facet families and their
//! surface energies; the carving centre, the supercell tiling and the output
//! cell. The reference and the Wulff shape are previewed side by side.

use af_algorithms::nanocrystal::{
    build_nanocrystal, compute_bounding_radius, compute_shape_half_extents, compute_wulff_preview, shape_label, NanoBuildResult, NanoGenerationMode, NanoParams, NanoShape, WulffPlaneInput, WulffPreview,
};
use af_core::structure::Structure;
use eframe::egui;
use glam::Vec3;

use crate::builder_preview::{cell_lengths, Face, Orbit, Reference};
use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

const SHAPES: [NanoShape; 7] = [NanoShape::Sphere, NanoShape::Ellipsoid, NanoShape::Box, NanoShape::Cylinder, NanoShape::Octahedron, NanoShape::TruncatedOctahedron, NanoShape::Cuboctahedron];
const CENTRE_HINT: &str = "Automatic centring follows the reference atoms. Disable it to specify a Cartesian centre in Å. A shift smaller than a lattice spacing can still change the surface termination.";
const REPLICATE_HINT: &str = "Automatic mode sizes the source coverage. Manual mode uses the three integer repetitions; increase all necessary directions if the shape is clipped by insufficient source coverage.";
const CELL_HINT: &str = "Choose a rectangular enclosing output cell where needed and add the requested surrounding vacuum in Å. Inspect the saved cell in a periodic format; a plain XYZ consumer may ignore the cell.";
const WULFF_HINT: &str = "Max radius sets the distance of the highest-energy family, with others proportional to energy. Family colour changes the plane preview, not surface energy. Build the atoms again after changing family energies or geometry.";
const SUCCESS: egui::Color32 = egui::Color32::from_rgb(70, 170, 90);
const WARN: egui::Color32 = egui::Color32::from_rgb(230, 170, 60);
const ERROR: egui::Color32 = egui::Color32::from_rgb(220, 80, 60);

/// The C++ wulffFamilyColor palette.
fn family_colour(i: usize) -> [f32; 3] {
    const PALETTE: [[f32; 3]; 10] = [
        [0.92, 0.33, 0.29],
        [0.18, 0.63, 0.58],
        [0.24, 0.45, 0.86],
        [0.96, 0.73, 0.20],
        [0.72, 0.38, 0.82],
        [0.13, 0.64, 0.28],
        [0.89, 0.52, 0.16],
        [0.35, 0.73, 0.87],
        [0.82, 0.23, 0.47],
        [0.62, 0.65, 0.19],
    ];
    PALETTE[i % PALETTE.len()]
}

fn mode_label(m: NanoGenerationMode) -> &'static str {
    match m {
        NanoGenerationMode::Shape => "Shape Cut",
        NanoGenerationMode::WulffConstruction => "Wulff Construction",
    }
}

fn length(ui: &mut egui::Ui, label: &str, v: &mut f32, max: f32) {
    ui.horizontal(|ui| {
        ui.add(egui::DragValue::new(v).speed(0.5).range(0.1..=max).fixed_decimals(2));
        ui.label(label);
    })
    .response
    .on_hover_text("Drag to adjust, double-click to type");
}

pub struct NanoDialog {
    open: bool,
    params: NanoParams,
    last: NanoBuildResult,
    reference: Reference,
    /// Bumped on every reference change (the Wulff preview's cache key).
    reference_version: u64,
    family_colours: Vec<[f32; 3]>,
    wulff: Option<(NanoParams, u64, WulffPreview)>,
    orbit: Orbit,
    wulff_orbit: Orbit,
    /// The pipeline step being edited (Edit in Dialog).
    step: Option<StepEdit>,
}

impl Default for NanoDialog {
    fn default() -> Self {
        // The C++ dialog opens on the Wulff construction.
        let params = NanoParams { generation_mode: NanoGenerationMode::WulffConstruction, wulff_planes: vec![WulffPlaneInput::default()], ..NanoParams::default() };
        Self {
            open: false,
            params,
            last: NanoBuildResult::default(),
            reference: Reference::new("Load reference (Nanocrystal)"),
            reference_version: 0,
            family_colours: Vec::new(),
            wulff: None,
            orbit: Orbit::default(),
            wulff_orbit: Orbit::default(),
            step: None,
        }
    }
}

impl NanoDialog {
    /// Fills the dialog from a build-nano step (editStep); missing flags take
    /// runNano's defaults, the reference is --input or the step's input.
    pub fn load_step(&mut self, edit: &StepEdit) {
        let r = edit.reader();
        let mut p = NanoParams { generation_mode: NanoGenerationMode::Shape, ..NanoParams::default() };
        match r.text("--shape", "sphere").to_ascii_lowercase().as_str() {
            "ellipsoid" => p.shape = NanoShape::Ellipsoid,
            "box" => p.shape = NanoShape::Box,
            "cylinder" => p.shape = NanoShape::Cylinder,
            "octahedron" => p.shape = NanoShape::Octahedron,
            "truncated-octahedron" => p.shape = NanoShape::TruncatedOctahedron,
            "cuboctahedron" => p.shape = NanoShape::Cuboctahedron,
            "wulff" => p.generation_mode = NanoGenerationMode::WulffConstruction,
            _ => p.shape = NanoShape::Sphere,
        }
        for text in r.values("--facet") {
            // istream >> h >> k >> l >> energy
            let w: Vec<&str> = text.split_whitespace().collect();
            let int = |i: usize| w.get(i).and_then(|v| v.parse::<i32>().ok());
            if let (Some(h), Some(k), Some(l), Some(e)) = (int(0), int(1), int(2), w.get(3).and_then(|v| v.parse::<f32>().ok())) {
                p.wulff_planes.push(WulffPlaneInput { h, k, l, surface_energy: e });
            }
        }
        p.wulff_max_radius = r.number("--radius", 20.0) as f32;
        let radius = r.number("--radius", 15.0) as f32;
        p.sphere_radius = radius;
        p.oct_radius = radius;
        p.cubo_radius = radius;
        p.trunc_oct_radius = r.number("--radius", 18.0) as f32;
        p.ellip_rx = r.number("--rx", 15.0) as f32;
        p.ellip_ry = r.number("--ry", 12.0) as f32;
        p.ellip_rz = r.number("--rz", 10.0) as f32;
        p.box_hx = r.number("--hx", 15.0) as f32;
        p.box_hy = r.number("--hy", 15.0) as f32;
        p.box_hz = r.number("--hz", 15.0) as f32;
        p.trunc_oct_trunc = r.number("--trunc", 12.0) as f32;
        p.cyl_radius = r.number("--cylradius", 12.0) as f32;
        p.cyl_height = r.number("--cylheight", 30.0) as f32;
        p.cyl_axis = r.integer("--cylaxis", 2).clamp(0, 2);
        p.vacuum_padding = r.number("--vacuum", 5.0) as f32;
        p.set_output_cell = !r.has("--no-cell");
        let centre = r.numbers("--center");
        p.auto_center_from_atoms = centre.len() < 3;
        if !p.auto_center_from_atoms {
            (p.cx, p.cy, p.cz) = (centre[0] as f32, centre[1] as f32, centre[2] as f32);
        }
        let (a, b, c) = (r.integer("--repa", 0), r.integer("--repb", 0), r.integer("--repc", 0));
        p.auto_replicate = !(a > 0 || b > 0 || c > 0);
        if !p.auto_replicate {
            let or5 = |n: i32| if n > 0 { n } else { 5 };
            (p.rep_a, p.rep_b, p.rep_c) = (or5(a), or5(b), or5(c));
        }
        self.params = p;
        self.last = NanoBuildResult::default();
        self.family_colours.clear();
        self.reference.from_step_edit(edit, false, "(pipeline step input)");
        self.reference_changed();
    }

    /// The options of the matching pipeline step (stepOptions): only the
    /// selected shape's sizes.
    pub fn step_options(&self) -> String {
        let mut w = af_pipeline::options::Writer::new();
        if !self.reference.from_step && !self.reference.path.is_empty() {
            w.add("--input", &self.reference.path);
        }
        let p = &self.params;
        let f = |v: f32| v as f64;
        if p.generation_mode == NanoGenerationMode::WulffConstruction {
            w.add("--shape", "wulff");
            w.add_number("--radius", f(p.wulff_max_radius));
            for plane in &p.wulff_planes {
                w.add_numbers("--facet", &[plane.h as f64, plane.k as f64, plane.l as f64, f(plane.surface_energy)]);
            }
        } else {
            match p.shape {
                NanoShape::Ellipsoid => w.add("--shape", "ellipsoid").add_number("--rx", f(p.ellip_rx)).add_number("--ry", f(p.ellip_ry)).add_number("--rz", f(p.ellip_rz)),
                NanoShape::Box => w.add("--shape", "box").add_number("--hx", f(p.box_hx)).add_number("--hy", f(p.box_hy)).add_number("--hz", f(p.box_hz)),
                NanoShape::Cylinder => w.add("--shape", "cylinder").add_number("--cylradius", f(p.cyl_radius)).add_number("--cylheight", f(p.cyl_height)).add_integer("--cylaxis", p.cyl_axis as i64),
                NanoShape::Octahedron => w.add("--shape", "octahedron").add_number("--radius", f(p.oct_radius)),
                NanoShape::TruncatedOctahedron => w.add("--shape", "truncated-octahedron").add_number("--radius", f(p.trunc_oct_radius)).add_number("--trunc", f(p.trunc_oct_trunc)),
                NanoShape::Cuboctahedron => w.add("--shape", "cuboctahedron").add_number("--radius", f(p.cubo_radius)),
                _ => w.add("--shape", "sphere").add_number("--radius", f(p.sphere_radius)),
            };
        }
        if !p.auto_center_from_atoms {
            w.add_numbers("--center", &[f(p.cx), f(p.cy), f(p.cz)]);
        }
        if !p.auto_replicate {
            w.add_integer("--repa", p.rep_a.max(1) as i64).add_integer("--repb", p.rep_b.max(1) as i64).add_integer("--repc", p.rep_c.max(1) as i64);
        }
        if p.set_output_cell {
            w.add_number("--vacuum", f(p.vacuum_padding));
        } else {
            w.flag("--no-cell");
        }
        w.str().to_string()
    }

    /// Update Step: the step's new options (and use_input when the reference
    /// is the step's input); ends the step edit.
    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        self.open = false;
        self.last = NanoBuildResult::default();
        let mut updated = edit.with_options(&self.step_options());
        if self.reference.from_step {
            updated = updated.with_parameter("use_input", true);
        }
        Some(Outcome::UpdateStep(updated))
    }

    /// The Wulff shape of the current reference and parameters (recomputed on change).
    fn wulff(&mut self) -> Option<&WulffPreview> {
        if self.reference.structure.atoms.is_empty() || self.params.generation_mode != NanoGenerationMode::WulffConstruction {
            return None;
        }
        let stale = self.wulff.as_ref().is_none_or(|(p, v, _)| *p != self.params || *v != self.reference_version);
        if stale {
            let preview = compute_wulff_preview(&self.reference.structure, &self.params);
            self.wulff = Some((self.params.clone(), self.reference_version, preview));
        }
        self.wulff.as_ref().map(|w| &w.2)
    }

    fn can_build(&mut self) -> bool {
        !self.reference.structure.atoms.is_empty() && (self.params.generation_mode != NanoGenerationMode::WulffConstruction || self.wulff().is_some_and(|w| w.success))
    }

    fn build(&mut self, colours: &[[f32; 3]]) -> Option<Outcome> {
        let mut s = Structure::default();
        self.last = build_nanocrystal(&mut s, &self.reference.structure, &self.params, colours, &[], &[]);
        self.last.success.then(|| Outcome::Built(s, "Built a nanocrystal".into()))
    }

    fn reference_changed(&mut self) {
        self.reference_version += 1;
        self.orbit = Orbit::default();
        self.wulff_orbit = Orbit::default();
    }

    fn reference_panel(&mut self, ui: &mut egui::Ui, input: &Input) {
        ui.strong("Reference Structure");
        ui.separator();
        if !self.reference.structure.atoms.is_empty() {
            let mut clear = false;
            ui.horizontal(|ui| {
                ui.label(&self.reference.status);
                clear = ui.button("Clear").clicked();
            });
            if clear {
                self.reference.clear();
                self.reference_changed();
            }
        }
        let height = (ui.available_height() - 2.0 * ui.text_style_height(&egui::TextStyle::Body) - 16.0).max(60.0);
        crate::builder_preview::structure(ui, egui::vec2(ui.available_width(), height), &self.reference.structure, input.colours, &mut self.orbit, &["Load a structure file", "(.cif .xyz .pdb .sdf .mol .vasp ...)"]);
        let r = &self.reference.structure;
        if !r.atoms.is_empty() {
            ui.weak("Drag to orbit, scroll to zoom");
            if r.has_unit_cell {
                let [a, b, c] = cell_lengths(r);
                ui.weak(format!("{} atoms | a={a:.3}  b={b:.3}  c={c:.3} Å | periodic", r.atoms.len()));
            } else {
                ui.weak(format!("{} atoms | no unit cell | shape carving only", r.atoms.len()));
            }
        }
    }

    fn wulff_panel(&mut self, ui: &mut egui::Ui) {
        ui.strong("Wulff Plane View");
        ui.separator();
        if self.params.generation_mode != NanoGenerationMode::WulffConstruction {
            ui.label("Switch generation mode to Wulff Construction to preview symmetry-expanded planes.");
            return;
        }
        if self.reference.structure.atoms.is_empty() {
            ui.label("Load a periodic reference structure to generate the Wulff plane preview.");
            return;
        }
        let colours = self.family_colours.clone();
        let Some(w) = self.wulff().cloned() else { return };
        if !w.success {
            ui.label(&w.message);
            return;
        }
        ui.weak(format!("{} faces from {} symmetry-expanded planes", w.faces.len(), w.planes.len()));
        let faces: Vec<Face> = w.faces.iter().map(|f| Face { vertices: f.vertices.clone(), normal: f.normal, colour: colours.get(f.family_index.max(0) as usize).copied().unwrap_or_else(|| family_colour(f.family_index.max(0) as usize)) }).collect();
        let height = (ui.available_height() - ui.text_style_height(&egui::TextStyle::Body) - 10.0).max(80.0);
        crate::builder_preview::polyhedron(ui, egui::vec2(ui.available_width(), height), &faces, &mut self.wulff_orbit);
        ui.weak("Left-drag = orbit   Scroll = zoom");
    }

    fn shape_parameters(p: &mut NanoParams, ui: &mut egui::Ui) {
        match p.shape {
            NanoShape::Sphere => length(ui, "Radius (Å)", &mut p.sphere_radius, 500.0),
            NanoShape::Ellipsoid => {
                length(ui, "Semi-axis a (Å)", &mut p.ellip_rx, 500.0);
                length(ui, "Semi-axis b (Å)", &mut p.ellip_ry, 500.0);
                length(ui, "Semi-axis c (Å)", &mut p.ellip_rz, 500.0);
            }
            NanoShape::Box => {
                length(ui, "Half-size X (Å)", &mut p.box_hx, 500.0);
                length(ui, "Half-size Y (Å)", &mut p.box_hy, 500.0);
                length(ui, "Half-size Z (Å)", &mut p.box_hz, 500.0);
                ui.weak(format!("Full box: {:.2} × {:.2} × {:.2} Å", 2.0 * p.box_hx, 2.0 * p.box_hy, 2.0 * p.box_hz));
            }
            NanoShape::Cylinder => {
                length(ui, "Radius (Å)", &mut p.cyl_radius, 500.0);
                length(ui, "Height (Å)", &mut p.cyl_height, 1000.0);
                ui.horizontal(|ui| {
                    for (i, axis) in ["X", "Y", "Z"].into_iter().enumerate() {
                        crate::dialog_widgets::radio_value(ui, &mut p.cyl_axis, i as i32, axis);
                    }
                    ui.label("Axis");
                });
            }
            NanoShape::Octahedron => {
                length(ui, "Radius (Å)", &mut p.oct_radius, 500.0);
                ui.weak("Condition: |x|+|y|+|z| <= R");
            }
            NanoShape::TruncatedOctahedron => {
                length(ui, "Octahedron R (Å)", &mut p.trunc_oct_radius, 500.0);
                length(ui, "Truncation R (Å)", &mut p.trunc_oct_trunc, 500.0);
                ui.weak("Cond: |x|+|y|+|z|<=R_oct  AND  max(|x|,|y|,|z|)<=R_trunc");
            }
            NanoShape::Cuboctahedron => {
                length(ui, "Radius (Å)", &mut p.cubo_radius, 500.0);
                ui.weak("Cond: |x|+|y|, |y|+|z|, |x|+|z| <= R");
            }
            NanoShape::MeshModel => {}
        }
    }

    fn wulff_parameters(&mut self, ui: &mut egui::Ui) {
        let p = &mut self.params;
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut p.wulff_max_radius).speed(0.25).range(0.1..=1000.0).fixed_decimals(2));
            ui.label("Max radius (Å)");
        })
        .response
        .on_hover_text(WULFF_HINT);
        if p.wulff_planes.is_empty() {
            p.wulff_planes.push(WulffPlaneInput::default());
        }
        while self.family_colours.len() < p.wulff_planes.len() {
            self.family_colours.push(family_colour(self.family_colours.len()));
        }
        self.family_colours.truncate(p.wulff_planes.len());
        ui.add_space(4.0);
        ui.label("Facet families");
        ui.weak("Distances follow energy ratios; the highest-energy family sits at the max radius.");
        let mut remove = None;
        let single = p.wulff_planes.len() == 1;
        for (i, plane) in p.wulff_planes.iter_mut().enumerate() {
            ui.separator();
            ui.horizontal(|ui| {
                ui.color_edit_button_rgb(&mut self.family_colours[i]);
                ui.label(format!("Family {}", i + 1));
                if ui.button("Remove").clicked() && !single {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                for (v, l) in [(&mut plane.h, "h"), (&mut plane.k, "k"), (&mut plane.l, "l")] {
                    ui.add(egui::DragValue::new(v).speed(0.2).range(-12..=12));
                    ui.label(l);
                }
                ui.add(egui::DragValue::new(&mut plane.surface_energy).speed(0.01).range(0.001..=1000.0).fixed_decimals(4));
                ui.label("Energy");
            });
        }
        if let Some(i) = remove {
            p.wulff_planes.remove(i);
            self.family_colours.remove(i);
        }
        if ui.add(egui::Button::new("Add Facet Family").min_size(egui::vec2(150.0, 0.0))).clicked() {
            p.wulff_planes.push(WulffPlaneInput::default());
            self.family_colours.push(family_colour(self.family_colours.len()));
        }
    }

    /// The tiling the build will use, and whether it is refused (drawReplicationPreview).
    fn replication_preview(&self, ui: &mut egui::Ui) {
        let r = &self.reference.structure;
        if !r.has_unit_cell {
            return;
        }
        let [la, lb, lc] = cell_lengths(r);
        if la < 1e-8 || lb < 1e-8 || lc < 1e-8 {
            return;
        }
        let p = &self.params;
        let max_r = compute_bounding_radius(p);
        let mut clamped = false;
        let mut rep = |l: f32| {
            let n = (max_r / l).ceil() as i32 + 2;
            if n > 40 {
                clamped = true;
                40
            } else {
                n
            }
        };
        let (a, b, c) = if p.auto_replicate { (rep(la), rep(lb), rep(lc)) } else { (p.rep_a.max(1), p.rep_b.max(1), p.rep_c.max(1)) };
        let tested = (2 * a as i64 + 1) * (2 * b as i64 + 1) * (2 * c as i64 + 1) * r.atoms.len() as i64;
        ui.weak(format!("Tiling: {} × {} × {}  ({tested} atoms tested)", 2 * a + 1, 2 * b + 1, 2 * c + 1));
        if clamped {
            ui.colored_label(WARN, "Warning: replication clamped to 40.");
        }
        if tested > 8_000_000 {
            ui.colored_label(ERROR, "Exceeds 8M limit -- build will be refused.");
        }
    }

    fn result(&self, ui: &mut egui::Ui) {
        let r = &self.last;
        if r.message.is_empty() {
            return;
        }
        ui.colored_label(if r.success { SUCCESS } else { ERROR }, format!("Status: {}", r.message));
        if !r.success {
            return;
        }
        ui.separator();
        egui::Grid::new("nano result").num_columns(2).spacing([16.0, 2.0]).show(ui, |ui| {
            let mut row = |label: &str, value: String| {
                ui.weak(label);
                ui.label(value);
                ui.end_row();
            };
            row("Mode", mode_label(r.mode).into());
            if r.mode == NanoGenerationMode::Shape {
                row("Shape", shape_label(r.shape).into());
            }
            row("Input atoms", r.input_atoms.to_string());
            row("Output atoms", r.output_atoms.to_string());
            row("Est. diameter", format!("{:.2} A", r.estimated_diameter));
            if r.mode == NanoGenerationMode::WulffConstruction {
                row("Wulff faces", r.wulff_face_count.to_string());
            }
            if r.tiling_used {
                row("Tiling", format!("{}x{}x{} ({})", 2 * r.rep_a + 1, 2 * r.rep_b + 1, 2 * r.rep_c + 1, if r.rep_clamped { "clamped" } else { "auto" }));
            } else {
                row("Tiling", "none (no unit cell in reference)".into());
            }
        });
    }

    fn options(&mut self, ui: &mut egui::Ui) {
        ui.strong("Builder Options");
        ui.separator();
        let before = self.params.generation_mode;
        egui::ComboBox::from_id_salt("nano generation mode").width(ui.available_width()).selected_text(mode_label(self.params.generation_mode)).show_ui(ui, |ui| {
            for m in [NanoGenerationMode::Shape, NanoGenerationMode::WulffConstruction] {
                ui.selectable_value(&mut self.params.generation_mode, m, mode_label(m));
            }
        });
        if self.params.generation_mode != before {
            self.last = NanoBuildResult::default();
            self.wulff_orbit = Orbit::default();
        }
        ui.separator();
        if self.params.generation_mode == NanoGenerationMode::Shape {
            ui.label("Shape");
            if !SHAPES.contains(&self.params.shape) {
                self.params.shape = NanoShape::Sphere;
            }
            egui::ComboBox::from_id_salt("nano shape").width(ui.available_width()).selected_text(shape_label(self.params.shape)).show_ui(ui, |ui| {
                for s in SHAPES {
                    ui.selectable_value(&mut self.params.shape, s, shape_label(s));
                }
            });
            ui.separator();
        }
        let scroll_height = (ui.available_height() - 140.0).max(160.0);
        egui::ScrollArea::vertical().id_salt("nano parameters").max_height(scroll_height).show(ui, |ui| {
            if self.params.generation_mode == NanoGenerationMode::Shape {
                Self::shape_parameters(&mut self.params, ui);
            } else {
                self.wulff_parameters(ui);
            }
            ui.separator();
            ui.label("Carving center");
            let p = &mut self.params;
            crate::dialog_widgets::checkbox(ui, &mut p.auto_center_from_atoms, "Auto-center from atoms").on_hover_text(CENTRE_HINT);
            if !p.auto_center_from_atoms {
                for (v, l) in [(&mut p.cx, "X (Å)"), (&mut p.cy, "Y (Å)"), (&mut p.cz, "Z (Å)")] {
                    ui.horizontal(|ui| {
                        ui.add(egui::DragValue::new(v).speed(0.1).range(-1000.0..=1000.0).fixed_decimals(3));
                        ui.label(l);
                    })
                    .response
                    .on_hover_text(CENTRE_HINT);
                }
            }
            if self.reference.structure.has_unit_cell {
                ui.separator();
                ui.label("Supercell replication");
                let p = &mut self.params;
                crate::dialog_widgets::checkbox(ui, &mut p.auto_replicate, "Auto-replicate").on_hover_text(format!("{REPLICATE_HINT}\nAutomatically tile the unit cell to fully cover the generated shape"));
                if !p.auto_replicate {
                    for (v, l) in [(&mut p.rep_a, "Along a"), (&mut p.rep_b, "Along b"), (&mut p.rep_c, "Along c")] {
                        ui.horizontal(|ui| {
                            ui.add(egui::DragValue::new(v).speed(0.2).range(1..=40));
                            ui.label(l);
                        });
                    }
                }
                self.replication_preview(ui);
            }
            ui.separator();
            ui.label("Output cell");
            let p = &mut self.params;
            crate::dialog_widgets::checkbox(ui, &mut p.set_output_cell, "Rectangular cell").on_hover_text(CELL_HINT);
            if p.set_output_cell {
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut p.vacuum_padding).speed(0.1).range(0.0..=f32::MAX).fixed_decimals(2));
                    ui.label("Vacuum (Å)");
                })
                .response
                .on_hover_text(CELL_HINT);
                let vacuum = p.vacuum_padding;
                if p.generation_mode == NanoGenerationMode::Shape {
                    let he = compute_shape_half_extents(p);
                    ui.weak(format!("{:.2} × {:.2} × {:.2} Å", 2.0 * (he.hx + vacuum), 2.0 * (he.hy + vacuum), 2.0 * (he.hz + vacuum)));
                } else if let Some(w) = self.wulff().filter(|w| w.success) {
                    let e: Vec3 = w.max_point - w.min_point;
                    ui.weak(format!("{:.2} × {:.2} × {:.2} Å", e.x + 2.0 * vacuum, e.y + 2.0 * vacuum, e.z + 2.0 * vacuum));
                }
            }
            if let Some(w) = self.wulff().filter(|w| !w.message.is_empty()) {
                let message = w.message.clone();
                ui.add_space(4.0);
                ui.label(format!("Preview: {message}"));
            }
        });
        ui.separator();
        ui.label("Result");
        self.result(ui);
    }
}

impl StructureDialog for NanoDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        let _ = self.reference.load(path, false);
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-nano"]
    }

    fn title(&self) -> &'static str {
        "Build Nanocrystal"
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
        self.step = None;
        self.last = NanoBuildResult::default();
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        let mut built = None;
        let mut open = self.open;
        let version = (self.reference.structure.atoms.len(), self.reference.name.clone());
        crate::theme::modal_window("Build Nanocrystal").open(&mut open).default_size([1100.0, 760.0]).show(ctx, |ui| {
            if let Some(edit) = &self.step {
                crate::structure_dialog::step_banner(ui, &edit.step);
            }
            self.reference.buttons(ui, input.source, false);
            ui.label("Load a reference crystal (or use the structure in the view), then build the nanocrystal.");
            ui.separator();
            let column = 700.0f32.min(ctx.content_rect().height() - 220.0).max(360.0);
            let stacked = ui.available_width() < 800.0;
            let left_width = if stacked { ui.available_width() } else { 460.0 };
            let left = |me: &mut Self, ui: &mut egui::Ui| {
                ui.vertical(|ui| {
                    let top = column * 0.52;
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        ui.set_width(left_width - 14.0);
                        ui.set_height(top);
                        me.reference_panel(ui, input);
                    });
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        ui.set_width(left_width - 14.0);
                        ui.set_height(column - top - 20.0);
                        me.wulff_panel(ui);
                    });
                });
            };
            let right = |me: &mut Self, ui: &mut egui::Ui| {
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.set_height(column);
                    me.options(ui);
                });
            };
            if stacked {
                left(self, ui);
                right(self, ui);
            } else {
                ui.horizontal_top(|ui| {
                    ui.vertical(|ui| left(self, ui));
                    ui.vertical(|ui| right(self, ui));
                });
            }
            // A step may be updated before its input exists; building needs a reference.
            let editing = self.step.is_some();
            let can = editing || self.can_build();
            let label = if editing { "Update Step" } else { "Build Nanocrystal" };
            match crate::dialog_layout::footer(ui, label, Some("Close"), can) {
                1 => built = if editing { self.update_step() } else { self.build(input.colours) },
                2 => self.open = false,
                _ => {}
            }
        });
        self.open = open && self.open;
        if !self.open {
            self.step = None;
        }
        self.reference.browser(ctx, false);
        if version != (self.reference.structure.atoms.len(), self.reference.name.clone()) {
            self.reference_changed();
        }
        built
    }

    /// reference PATH|view|clear | mode shape|wulff | shape NAME |
    /// size A [B C] (the shape's radius / semi-axes / half-sizes / radius height) |
    /// axis x|y|z | truncation R | wulff-radius R | facets "H K L E; ..." |
    /// centre auto | centre X Y Z | replicate auto | replicate A B C |
    /// cell on|off | vacuum V | build | expect-faces N | expect-output MIN [MAX]
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let num = |i: usize| words.get(i).and_then(|w| w.parse::<f64>().ok()).ok_or(format!("{command}: expected a number in {args:?}"));
        let p = &mut self.params;
        match command {
            "reference" => {
                match args.trim() {
                    "view" => self.reference.take_view(input.source, false)?,
                    "clear" => self.reference.clear(),
                    path => self.reference.load(path, false)?,
                }
                self.reference_changed();
            }
            "mode" => {
                p.generation_mode = match args.trim() {
                    "shape" => NanoGenerationMode::Shape,
                    "wulff" => NanoGenerationMode::WulffConstruction,
                    other => return Err(format!("no generation mode {other}")),
                }
            }
            "shape" => p.shape = *SHAPES.iter().find(|s| shape_label(**s).eq_ignore_ascii_case(args.trim())).ok_or(format!("no shape {args}"))?,
            "size" => {
                let a = num(0)? as f32;
                let (b, c) = (num(1).map(|v| v as f32).unwrap_or(a), num(2).map(|v| v as f32).unwrap_or(a));
                match p.shape {
                    NanoShape::Sphere => p.sphere_radius = a,
                    NanoShape::Ellipsoid => (p.ellip_rx, p.ellip_ry, p.ellip_rz) = (a, b, c),
                    NanoShape::Box => (p.box_hx, p.box_hy, p.box_hz) = (a, b, c),
                    NanoShape::Cylinder => (p.cyl_radius, p.cyl_height) = (a, num(1)? as f32),
                    NanoShape::Octahedron => p.oct_radius = a,
                    NanoShape::TruncatedOctahedron => p.trunc_oct_radius = a,
                    NanoShape::Cuboctahedron => p.cubo_radius = a,
                    NanoShape::MeshModel => return Err("no mesh models here".into()),
                }
            }
            "axis" => p.cyl_axis = ["x", "y", "z"].iter().position(|a| a.eq_ignore_ascii_case(args.trim())).ok_or("axis x|y|z")? as i32,
            "truncation" => p.trunc_oct_trunc = num(0)? as f32,
            "wulff-radius" => p.wulff_max_radius = num(0)? as f32,
            "facets" => {
                let mut planes = Vec::new();
                for part in args.split(';').filter(|s| !s.trim().is_empty()) {
                    let w: Vec<f64> = part.split_whitespace().map(|x| x.parse::<f64>()).collect::<Result<_, _>>().map_err(|_| format!("a facet is H K L ENERGY, not {part:?}"))?;
                    if w.len() != 4 {
                        return Err(format!("a facet is H K L ENERGY, not {part:?}"));
                    }
                    planes.push(WulffPlaneInput { h: w[0] as i32, k: w[1] as i32, l: w[2] as i32, surface_energy: w[3] as f32 });
                }
                p.wulff_planes = planes;
            }
            "centre" if args.trim() == "auto" => p.auto_center_from_atoms = true,
            "centre" => {
                p.auto_center_from_atoms = false;
                (p.cx, p.cy, p.cz) = (num(0)? as f32, num(1)? as f32, num(2)? as f32);
            }
            "replicate" if args.trim() == "auto" => p.auto_replicate = true,
            "replicate" => {
                p.auto_replicate = false;
                (p.rep_a, p.rep_b, p.rep_c) = (num(0)? as i32, num(1)? as i32, num(2)? as i32);
            }
            "cell" => p.set_output_cell = args.trim() == "on",
            "vacuum" => p.vacuum_padding = num(0)? as f32,
            "close" => {
                self.open = false;
                self.step = None;
            }
            "update-step" => return self.update_step().map(Some).ok_or_else(|| "no pipeline step is being edited".to_string()),
            "expect-options" => {
                if self.step_options() != args.trim() {
                    return Err(format!("step options are {:?}", self.step_options()));
                }
            }
            "build" => {
                if !self.can_build() {
                    let why = self.wulff().map(|w| w.message.clone()).filter(|m| !m.is_empty()).unwrap_or_else(|| "load a reference first".into());
                    return Err(why);
                }
                let out = self.build(input.colours);
                return if out.is_some() { Ok(out) } else { Err(self.last.message.clone()) };
            }
            "expect-faces" => {
                let n = num(0)? as usize;
                let w = self.wulff().ok_or("no Wulff preview")?;
                if w.faces.len() != n {
                    return Err(format!("{} Wulff faces: {}", w.faces.len(), w.message));
                }
            }
            "expect-output" => {
                let n = self.last.output_atoms as f64;
                if n < num(0)? || num(1).is_ok_and(|max| n > max) {
                    return Err(format!("{n} atoms: {}", self.last.message));
                }
            }
            other => return Err(format!("the nanocrystal dialog has no {other}")),
        }
        Ok(None)
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        self.load_step(&edit);
        self.step = Some(edit);
        self.open = true;
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edit(options: &str) -> StepEdit {
        let mut parameters = af_core::json::Json::object();
        parameters.set("options", options);
        StepEdit { index: 0, step: "build-nano".into(), parameters, input: Default::default(), ..Default::default() }
    }

    #[test]
    fn step_options_round_trip() {
        for options in [
            "--shape sphere --radius 10 --vacuum 5",
            "--shape wulff --radius 20 --facet \"1 0 0 1\" --facet \"1 1 1 1.200000048\" --vacuum 5",
            "--shape ellipsoid --rx 15 --ry 12 --rz 10 --center \"1 2 3\" --repa 4 --repb 5 --repc 6 --no-cell",
            "--shape box --hx 15 --hy 15 --hz 15 --vacuum 5",
            "--shape cylinder --cylradius 12 --cylheight 30 --cylaxis 0 --vacuum 5",
            "--shape octahedron --radius 15 --vacuum 5",
            "--shape truncated-octahedron --radius 18 --trunc 12 --vacuum 5",
            "--shape cuboctahedron --radius 15 --vacuum 2.5",
        ] {
            let mut d = NanoDialog::default();
            assert!(d.edit_step(edit(options)));
            assert_eq!(d.step_options(), options);
        }
        // Missing flags take runNano's defaults; one --radius sizes the truncated octahedron.
        let mut d = NanoDialog::default();
        d.edit_step(edit("--shape truncated-octahedron"));
        assert_eq!(d.step_options(), "--shape truncated-octahedron --radius 18 --trunc 12 --vacuum 5");
        d.edit_step(edit("--input cu.cif --repb 3"));
        assert_eq!(d.step_options(), "--input cu.cif --shape sphere --radius 15 --repa 5 --repb 3 --repc 5 --vacuum 5");
    }
}
