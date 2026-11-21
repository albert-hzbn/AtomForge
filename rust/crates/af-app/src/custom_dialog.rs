//! Build > Custom Structure (the C++ ui/CustomStructureDialog): a crystal
//! reference (the scene or a file) fills a 3D model (OBJ or STL) scaled to
//! Ångström, optionally oriented by angles or a Miller direction, with the
//! reference replicated to cover it and a rectangular output cell.

use af_algorithms::mesh::{parse_obj_mesh, parse_stl_mesh};
use af_algorithms::nanocrystal::{build_nanocrystal, NanoBuildResult, NanoParams, NanoShape};
use af_core::structure::Structure;
use eframe::egui;
use glam::Vec3;

use crate::file_browser::{FileBrowser, Purpose};
use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

const WARN: egui::Color32 = egui::Color32::from_rgb(230, 170, 60);
const GOOD: egui::Color32 = egui::Color32::from_rgb(70, 170, 90);
const INFO: egui::Color32 = egui::Color32::from_rgb(110, 170, 230);

pub struct CustomDialog {
    open: bool,
    params: NanoParams,
    reference: Structure,
    reference_name: String,
    vertices: Vec<Vec3>,
    indices: Vec<u32>,
    model_name: String,
    half_extents: Vec3,
    /// 1: rotation angles, 2: Miller direction along +z.
    orientation_type: i32,
    rotation: [f32; 3],
    miller: [i32; 3],
    result: NanoBuildResult,
    status: String,
    browser: FileBrowser,
    /// Where the reference and the model came from (for the step's options).
    reference_path: String,
    reference_from_step_input: bool,
    model_path: String,
    /// The pipeline step being edited ("Update Step" instead of Build Fill).
    step: Option<StepEdit>,
}

impl Default for CustomDialog {
    fn default() -> Self {
        Self {
            open: false,
            params: NanoParams { shape: NanoShape::MeshModel, ..Default::default() },
            reference: Structure::default(),
            reference_name: String::new(),
            vertices: Vec::new(),
            indices: Vec::new(),
            model_name: String::new(),
            half_extents: Vec3::ZERO,
            orientation_type: 1,
            rotation: [0.0; 3],
            miller: [0, 0, 1],
            result: NanoBuildResult::default(),
            status: String::new(),
            browser: FileBrowser::default(),
            reference_path: String::new(),
            reference_from_step_input: false,
            model_path: String::new(),
            step: None,
        }
    }
}

fn file_name(path: &str) -> String {
    std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string())
}

impl CustomDialog {
    /// A dropped or chosen file: OBJ/STL to the model, anything else to the reference.
    fn load(&mut self, path: &str) -> bool {
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".obj") || lower.ends_with(".stl") {
            self.load_model(path)
        } else {
            self.load_reference(path)
        }
    }

    fn load_reference(&mut self, path: &str) -> bool {
        match af_io::loader::load_structure_from_file(path) {
            Ok(s) => {
                self.reference = s;
                self.reference_name = file_name(path);
                self.reference_path = path.to_string();
                self.reference_from_step_input = false;
                self.status = format!("Crystal reference loaded: {} (atoms={}, unitCell={})", self.reference_name, self.reference.atoms.len(), if self.reference.has_unit_cell { "yes" } else { "no" });
                true
            }
            Err(e) => {
                self.status = format!("Failed to load crystal structure: {e}");
                false
            }
        }
    }

    fn load_model(&mut self, path: &str) -> bool {
        let p = std::path::Path::new(path);
        let lower = path.to_ascii_lowercase();
        let mesh = if lower.ends_with(".obj") {
            parse_obj_mesh(p)
        } else if lower.ends_with(".stl") {
            parse_stl_mesh(p)
        } else {
            Err(String::new())
        };
        let mesh = match mesh {
            Ok(m) => m,
            Err(e) => {
                self.status = format!("Failed to load model: {}", if e.is_empty() { "Unsupported model format. Use OBJ or STL.".to_string() } else { e });
                return false;
            }
        };
        let (mut lo, mut hi) = (Vec3::splat(1e30), Vec3::splat(-1e30));
        for v in &mesh.vertices {
            lo = lo.min(*v);
            hi = hi.max(*v);
        }
        let centre = 0.5 * (lo + hi);
        self.half_extents = 0.5 * (hi - lo);
        self.vertices = mesh.vertices.iter().map(|v| *v - centre).collect();
        self.indices = mesh.indices;
        self.model_name = file_name(path);
        self.model_path = path.to_string();
        let h = self.half_extents;
        self.status = format!("Model loaded: {} (triangles={}, bbox={}x{}x{})", self.model_name, self.indices.len() / 3, (2.0 * h.x) as i32, (2.0 * h.y) as i32, (2.0 * h.z) as i32);
        true
    }

    /// The dialog as the step's options (stepOptions).
    pub fn step_options(&self) -> String {
        let p = &self.params;
        let mut out = af_pipeline::options::Writer::new();
        if !self.reference_from_step_input && !self.reference_path.is_empty() {
            out.add("--input", &self.reference_path);
        }
        if !self.model_path.is_empty() {
            out.add("--mesh", &self.model_path);
        }
        out.add_number("--scale", p.model_scale as f64);
        if p.set_output_cell {
            out.add_number("--vacuum", p.vacuum_padding as f64);
        } else {
            out.flag("--no-cell");
        }
        if !p.auto_center_from_atoms {
            out.add_numbers("--center", &[p.cx as f64, p.cy as f64, p.cz as f64]);
        }
        if !p.auto_replicate {
            out.add_integer("--repa", p.rep_a.max(1) as i64).add_integer("--repb", p.rep_b.max(1) as i64).add_integer("--repc", p.rep_c.max(1) as i64);
        }
        if p.apply_crystal_orientation {
            if self.orientation_type == 2 {
                out.add("--miller", &format!("{} {} {}", self.miller[0], self.miller[1], self.miller[2]));
            } else {
                out.add_number("--rotx", self.rotation[0] as f64).add_number("--roty", self.rotation[1] as f64).add_number("--rotz", self.rotation[2] as f64);
            }
        }
        out.str().to_string()
    }

    /// "Update Step": the step's new options (use_input set when the reference
    /// is the step's input; without --input the pipeline passes it).
    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        let mut updated = edit.with_options(&self.step_options());
        if self.reference_from_step_input {
            updated = updated.with_parameter("use_input", true);
        }
        self.open = false;
        Some(Outcome::UpdateStep(updated))
    }

    fn use_step_input(&mut self) {
        let Some(input) = self.step.as_ref().map(|e| e.input.clone()) else { return };
        if input.atoms.is_empty() {
            self.status = "The step has no input structure. Load a structure first.".into();
            return;
        }
        self.reference = input;
        self.reference_name = "pipeline step input".into();
        self.reference_path.clear();
        self.reference_from_step_input = true;
        self.status = "Pipeline step input used as crystal reference.".into();
    }

    fn use_scene(&mut self, scene: &Structure) {
        self.reference_path.clear();
        self.reference_from_step_input = false;
        if scene.atoms.is_empty() {
            self.status = "Scene is empty. Load a structure first.".into();
            return;
        }
        self.reference = scene.clone();
        if scene.has_unit_cell {
            self.reference_name = "scene (unit cell)".into();
            self.status = "Current scene copied as crystal reference.".into();
        } else {
            self.reference_name = "scene (no unit cell)".into();
            self.status = "Warning: scene has no unit cell - tiling disabled.".into();
        }
    }

    fn can_build(&self) -> bool {
        !self.reference.atoms.is_empty() && !self.vertices.is_empty() && !self.indices.is_empty()
    }

    fn build(&mut self, colours: &[[f32; 3]]) -> Option<Outcome> {
        let p = &mut self.params;
        p.shape = NanoShape::MeshModel;
        p.model_hx = self.half_extents.x;
        p.model_hy = self.half_extents.y;
        p.model_hz = self.half_extents.z;
        if p.apply_crystal_orientation {
            p.use_miller_orientation = self.orientation_type == 2;
            if p.use_miller_orientation {
                (p.miller_h, p.miller_k, p.miller_l) = (self.miller[0], self.miller[1], self.miller[2]);
            } else {
                (p.orient_x_deg, p.orient_y_deg, p.orient_z_deg) = (self.rotation[0], self.rotation[1], self.rotation[2]);
            }
        } else {
            p.use_miller_orientation = false;
        }
        let mut s = Structure::default();
        self.result = build_nanocrystal(&mut s, &self.reference, &self.params, colours, &self.vertices, &self.indices);
        self.status = self.result.message.clone();
        self.result.success.then(|| Outcome::Built(s, "Built a custom structure".into()))
    }

    fn controls(&mut self, ui: &mut egui::Ui, input: &Input) -> Option<Outcome> {
        if ui.button("Load Structure or Model").clicked() {
            self.browser.show_for(Purpose::Open, "Load structure or model", "", &[]);
        }
        ui.label("Crystal structure files (CIF, XYZ, VASP, PDB, ...) go to the left panel. 3D model files (OBJ, STL) go to the right panel automatically.");
        ui.separator();
        ui.columns(2, |cols| {
            let ui = &mut cols[0];
            ui.strong("Crystal Structure");
            if self.reference.atoms.is_empty() {
                ui.weak("Load a crystal structure file");
                ui.weak("(CIF, XYZ, PDB, VASP, ...)");
            } else {
                ui.colored_label(GOOD, &self.reference_name);
                ui.weak(format!("{} atoms  |  cell: {}", self.reference.atoms.len(), if self.reference.has_unit_cell { "yes" } else { "no" }));
                if self.reference.has_unit_cell {
                    let len = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                    let c = &self.reference.cell_vectors;
                    ui.weak(format!("Cell: a={:.2}  b={:.2}  c={:.2} Å", len(c[0]), len(c[1]), len(c[2])));
                }
            }
            ui.horizontal(|ui| {
                if !self.reference.atoms.is_empty() && ui.button("Clear").clicked() {
                    self.reference = Structure::default();
                    self.reference_name.clear();
                    self.reference_path.clear();
                    self.reference_from_step_input = false;
                    self.status = "Crystal reference cleared.".into();
                }
                if self.step.is_some() {
                    // The step's reference is the structure entering it, not the scene.
                    if ui.button("Use Step Input").clicked() {
                        self.use_step_input();
                    }
                } else if ui.button("Use Current Scene").clicked() {
                    self.use_scene(input.source);
                }
            });
            let ui = &mut cols[1];
            ui.strong("3D Object");
            if self.vertices.is_empty() {
                ui.weak("Load a 3D model file");
                ui.weak("(OBJ or STL)");
            } else {
                let h = self.half_extents;
                let s = self.params.model_scale;
                ui.colored_label(GOOD, &self.model_name);
                ui.weak(format!("{} triangles", self.indices.len() / 3));
                ui.weak(format!("Model size: {:.1} × {:.1} × {:.1} (model units)", 2.0 * h.x, 2.0 * h.y, 2.0 * h.z));
                ui.colored_label(INFO, format!("= {:.1} × {:.1} × {:.1} Å", 2.0 * h.x * s, 2.0 * h.y * s, 2.0 * h.z * s));
                if ui.button("Clear").clicked() {
                    self.vertices.clear();
                    self.indices.clear();
                    self.model_name.clear();
                    self.model_path.clear();
                    self.status = "Model cleared.".into();
                }
            }
        });
        ui.separator();
        ui.strong("Fill Options");
        let p = &mut self.params;
        let hint = "Load the enclosing mesh and enter Å per model unit. A unit-radius mesh at scale 12 represents a 12 Å radius along the corresponding model axes. Verify the preview's physical extent before building.";
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut p.model_scale).speed(0.1).range(1e-4..=1e6).fixed_decimals(3)).on_hover_text(hint);
            ui.label("Scale (Å per model unit)");
            if !self.vertices.is_empty() {
                let h = self.half_extents * p.model_scale;
                ui.weak(format!("→ {:.1} × {:.1} × {:.1} Å", 2.0 * h.x, 2.0 * h.y, 2.0 * h.z));
            }
        });
        let centre_hint = "Automatic centring derives placement from the reference; manual values specify the centre. Inspect reference and mesh together to avoid filling an unintended region.";
        ui.checkbox(&mut p.auto_center_from_atoms, "Auto-center from reference atoms").on_hover_text(centre_hint);
        if !p.auto_center_from_atoms {
            ui.horizontal(|ui| {
                for (v, l) in [(&mut p.cx, "Center X"), (&mut p.cy, "Center Y"), (&mut p.cz, "Center Z")] {
                    ui.add(egui::DragValue::new(v).speed(0.1).fixed_decimals(3)).on_hover_text(centre_hint);
                    ui.label(l);
                }
            });
        }
        ui.separator();
        ui.strong("Crystal Orientation");
        ui.checkbox(&mut p.apply_crystal_orientation, "Apply crystal orientation").on_hover_text("Enables the orientation controls. The current radio labels 1 and 2 select Euler-angle and Miller-direction modes respectively. In mode 1 enter Rot X/Y/Z in degrees. In mode 2 enter h,k,l for direction alignment.");
        if p.apply_crystal_orientation {
            ui.horizontal(|ui| {
                ui.label("Orientation Type");
                ui.radio_value(&mut self.orientation_type, 1, "1");
                ui.radio_value(&mut self.orientation_type, 2, "2");
            });
            if self.orientation_type == 1 {
                ui.horizontal(|ui| {
                    for (v, l) in self.rotation.iter_mut().zip(["Rot X (°)", "Rot Y (°)", "Rot Z (°)"]) {
                        ui.add(egui::DragValue::new(v).speed(1.0).fixed_decimals(2));
                        ui.label(l);
                    }
                });
                ui.weak("Type 1: manual rotation about X/Y/Z axes.");
            } else {
                ui.horizontal(|ui| {
                    for (v, l) in self.miller.iter_mut().zip(["h", "k", "l"]) {
                        ui.add(egui::DragValue::new(v).speed(0.1));
                        ui.label(l);
                    }
                });
                ui.weak("Type 2: align crystal direction [h k l] with model +Z axis.");
            }
        }
        if self.reference.has_unit_cell {
            ui.checkbox(&mut p.auto_replicate, "Auto-replicate reference cell").on_hover_text("Ensure enough periodic sites to cover the whole mesh. Manual undercoverage can look like a broken mesh cut.");
            if !p.auto_replicate {
                ui.horizontal(|ui| {
                    for (v, l) in [(&mut p.rep_a, "Reps a"), (&mut p.rep_b, "Reps b"), (&mut p.rep_c, "Reps c")] {
                        ui.add(egui::DragValue::new(v).speed(0.1).range(1..=1000));
                        ui.label(l);
                    }
                });
            }
        }
        ui.horizontal(|ui| {
            ui.checkbox(&mut p.set_output_cell, "Set rectangular output cell");
            if p.set_output_cell {
                ui.add(egui::DragValue::new(&mut p.vacuum_padding).speed(0.1).range(0.0..=1e4).fixed_decimals(2));
                ui.label("Vacuum (Å)");
            }
        });
        ui.add_space(4.0);
        let can = self.can_build();
        if !can {
            let text = if self.reference.atoms.is_empty() && self.vertices.is_empty() {
                "Load a crystal structure and a 3D model to begin."
            } else if self.reference.atoms.is_empty() {
                "Load a crystal structure file (CIF, XYZ, ...) to continue."
            } else {
                "Load a 3D model (OBJ/STL) to continue."
            };
            ui.colored_label(WARN, text);
        }
        let mut built = None;
        let status = if self.status.is_empty() { String::new() } else { format!("Status: {}", self.status) };
        let label = if self.step.is_some() { "Update Step" } else { "Build Fill" };
        ui.separator();
        match crate::dialog_layout::action_bar(ui, &status, crate::dialog_layout::Note::Info, label, Some("Close"), can) {
            1 => built = if self.step.is_some() { self.update_step() } else { self.build(input.colours) },
            2 => self.open = false,
            _ => {}
        }
        built
    }
}

impl StructureDialog for CustomDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        self.load(path);
        true
    }
    fn takes_all_drops(&self) -> bool {
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-custom"]
    }

    fn title(&self) -> &'static str {
        "Custom Structure"
    }

    fn open(&mut self, _step: &str) {
        self.step = None;
        self.open = true;
        self.result = NanoBuildResult::default();
        self.status.clear();
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        let reader = edit.reader();
        self.result = NanoBuildResult::default();
        self.status.clear();
        // Reference crystal: --input FILE, else the structure entering the step.
        self.reference = Structure::default();
        self.reference_name.clear();
        self.reference_path.clear();
        self.reference_from_step_input = false;
        if reader.has("--input") {
            self.load_reference(&reader.text("--input", ""));
        } else if edit.use_input() && !edit.input.atoms.is_empty() {
            self.reference = edit.input.clone();
            self.reference_name = "pipeline step input".into();
            self.reference_from_step_input = true;
        }
        self.vertices.clear();
        self.indices.clear();
        self.model_name.clear();
        self.model_path.clear();
        self.half_extents = Vec3::splat(15.0);
        if reader.has("--mesh") {
            self.load_model(&reader.text("--mesh", ""));
        }
        // Missing flags take the defaults of AtomForge --build custom.
        let p = &mut self.params;
        p.model_scale = reader.number("--scale", 1.0) as f32;
        p.vacuum_padding = reader.number("--vacuum", 5.0) as f32;
        p.set_output_cell = !reader.has("--no-cell");
        let centre = reader.numbers("--center");
        p.auto_center_from_atoms = centre.len() < 3;
        if centre.len() >= 3 {
            (p.cx, p.cy, p.cz) = (centre[0] as f32, centre[1] as f32, centre[2] as f32);
        }
        let (a, b, c) = (reader.integer("--repa", 0), reader.integer("--repb", 0), reader.integer("--repc", 0));
        p.auto_replicate = !(a > 0 || b > 0 || c > 0);
        p.rep_a = if a > 0 { a } else { 5 };
        p.rep_b = if b > 0 { b } else { 5 };
        p.rep_c = if c > 0 { c } else { 5 };
        self.rotation = [reader.number("--rotx", 0.0) as f32, reader.number("--roty", 0.0) as f32, reader.number("--rotz", 0.0) as f32];
        if reader.has("--miller") {
            // --miller overrides the rotations, as in the command line.
            let m = reader.numbers("--miller");
            p.apply_crystal_orientation = true;
            self.orientation_type = 2;
            self.miller = [m.first().map_or(1, |v| *v as i32), m.get(1).map_or(0, |v| *v as i32), m.get(2).map_or(0, |v| *v as i32)];
        } else {
            self.orientation_type = 1;
            p.apply_crystal_orientation = self.rotation.iter().any(|&r| r != 0.0);
        }
        self.step = Some(edit);
        self.open = true;
        true
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if let Some(path) = self.browser.show(ctx) {
            self.load(&path.to_string_lossy());
        }
        if !self.open {
            return None;
        }
        let mut out = None;
        let mut open = self.open;
        let editing = self.step.as_ref().map(|e| e.step.clone());
        crate::theme::modal_window("Custom Structure").open(&mut open).default_size([720.0, 640.0]).min_width(640.0).show(ctx, |ui| {
            if let Some(step) = &editing {
                ui.colored_label(INFO, format!("Editing pipeline step: {step}"));
            }
            egui::ScrollArea::vertical().min_scrolled_height(420.0).show(ui, |ui| out = self.controls(ui, input));
        });
        self.open = open && self.open;
        if !self.open {
            // Closed: the next menu build is normal.
            self.step = None;
        }
        out
    }

    /// load PATH (structure or OBJ/STL) | use-scene | scale S | centre auto|X Y Z |
    /// orientation off|euler X Y Z|miller H K L | replicate auto|A B C |
    /// cell off|VACUUM | build | use-step-input | update-step | expect-status TEXT |
    /// expect-options TEXT.
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let nums = |from: usize| -> Result<Vec<f32>, String> { words[from.min(words.len())..].iter().map(|w| w.parse::<f32>().map_err(|e| format!("{w}: {e}"))).collect() };
        let p = &mut self.params;
        match command {
            "load" => {
                if !self.load(args.trim()) {
                    return Err(self.status.clone());
                }
            }
            "use-scene" => self.use_scene(input.source),
            "use-step-input" => self.use_step_input(),
            "update-step" => {
                if !self.can_build() {
                    return Err("load a crystal structure and a 3D model first".into());
                }
                return self.update_step().map(Some).ok_or_else(|| "no pipeline step is being edited".to_string());
            }
            "expect-options" => {
                if self.step_options() != args.trim() {
                    return Err(format!("options are {:?}", self.step_options()));
                }
            }
            "scale" => p.model_scale = args.trim().parse::<f32>().map_err(|e| e.to_string())?.max(1e-4),
            "centre" | "center" => {
                if args.trim() == "auto" {
                    p.auto_center_from_atoms = true;
                } else {
                    let v = nums(0)?;
                    if v.len() != 3 {
                        return Err("centre X Y Z".into());
                    }
                    p.auto_center_from_atoms = false;
                    (p.cx, p.cy, p.cz) = (v[0], v[1], v[2]);
                }
            }
            "orientation" => match words.first().copied() {
                Some("off") => p.apply_crystal_orientation = false,
                Some("euler") => {
                    let v = nums(1)?;
                    if v.len() != 3 {
                        return Err("orientation euler X Y Z".into());
                    }
                    p.apply_crystal_orientation = true;
                    self.orientation_type = 1;
                    self.rotation = [v[0], v[1], v[2]];
                }
                Some("miller") => {
                    let v = nums(1)?;
                    if v.len() != 3 {
                        return Err("orientation miller H K L".into());
                    }
                    p.apply_crystal_orientation = true;
                    self.orientation_type = 2;
                    self.miller = [v[0] as i32, v[1] as i32, v[2] as i32];
                }
                _ => return Err("orientation off|euler X Y Z|miller H K L".into()),
            },
            "replicate" => {
                if args.trim() == "auto" {
                    p.auto_replicate = true;
                } else {
                    let v = nums(0)?;
                    if v.len() != 3 {
                        return Err("replicate A B C".into());
                    }
                    p.auto_replicate = false;
                    (p.rep_a, p.rep_b, p.rep_c) = ((v[0] as i32).max(1), (v[1] as i32).max(1), (v[2] as i32).max(1));
                }
            }
            "cell" => {
                if args.trim() == "off" {
                    p.set_output_cell = false;
                } else {
                    p.set_output_cell = true;
                    p.vacuum_padding = args.trim().parse::<f32>().map_err(|e| e.to_string())?.max(0.0);
                }
            }
            "build" => {
                if !self.can_build() {
                    return Err("load a crystal structure and a 3D model first".into());
                }
                let o = self.build(input.colours);
                return if o.is_some() { Ok(o) } else { Err(self.status.clone()) };
            }
            "expect-status" => {
                if !self.status.contains(args.trim()) {
                    return Err(format!("status is {:?}", self.status));
                }
            }
            other => return Err(format!("the custom structure dialog has no command {other}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::json::Json;
    use af_core::structure::AtomSite;

    fn edit(options: &str) -> StepEdit {
        let mut parameters = Json::object();
        parameters.set("options", options);
        parameters.set("use_input", true);
        let mut input = Structure::default();
        input.atoms = vec![AtomSite::new("Cu", [0.0; 3])];
        StepEdit { index: 2, step: "build-custom".into(), parameters, input, ..Default::default() }
    }

    /// editStep then stepOptions, against the C++ (golden_step_options_composition);
    /// the tests run in the crate's folder, so the golden mesh path resolves.
    #[test]
    fn step_options_match_the_cpp() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/step_options_composition.json");
        let golden = Json::parse(&std::fs::read_to_string(path).unwrap()).unwrap();
        let cases = golden.find("custom").unwrap().items().unwrap();
        assert!(cases.len() >= 10);
        for case in cases {
            let options = case.find("options").unwrap().string().unwrap();
            let mut d = CustomDialog::default();
            assert!(d.edit_step(edit(options)));
            assert!(d.reference_from_step_input);
            assert_eq!(d.step_options(), case.find("out").unwrap().string().unwrap(), "options {options:?}");
        }
    }

    #[test]
    fn options_round_trip() {
        let mut d = CustomDialog::default();
        d.edit_step(edit("--scale 2.5 --no-cell --center \"1 2 3\" --repa 2 --repb 3 --repc 4 --rotx 30 --roty -15 --rotz 5"));
        let text = d.step_options();
        assert_eq!(text, "--scale 2.5 --no-cell --center \"1 2 3\" --repa 2 --repb 3 --repc 4 --rotx 30 --roty -15 --rotz 5");
        let mut again = CustomDialog::default();
        again.edit_step(edit(&text));
        assert_eq!(again.step_options(), text);
        match again.update_step() {
            Some(Outcome::UpdateStep(e)) => {
                assert_eq!(e.options(), text);
                assert_eq!(e.index, 2);
                assert!(e.use_input());
            }
            _ => panic!("no step update"),
        }
        // Opening from the menu ends a step edit.
        again.edit_step(edit(""));
        again.open("build-custom");
        assert!(again.update_step().is_none());
    }
}
