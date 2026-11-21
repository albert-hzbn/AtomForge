//! Edit > Insert Dislocation (the C++ ui/DislocationBuilderDialog): lattice
//! detection, character and region shape, crystallographic directions, the
//! line point, isotropic or anisotropic (Stroh) elasticity, dipoles and the
//! shape's parameters; the dislocation is built in the background and its
//! loop outline drawn over the view (View > Show > Dislocation Lines).

use std::sync::mpsc::{channel, Receiver};

use af_algorithms::dislocation::{
    build_dislocation, detect_dislocation_lattice, dislocation_lattice_family_name, DislocationCharacter, DislocationDetectionResult, DislocationParams, DislocationResult, DislocationShape, ElasticSymmetry,
};
use af_core::structure::Structure;
use eframe::egui;
use glam::{IVec3, Vec2, Vec3};

use af_pipeline::options::Writer;

use crate::structure_dialog::StepEdit;

const CHARACTERS: [(&str, DislocationCharacter); 3] = [("Edge", DislocationCharacter::Edge), ("Screw", DislocationCharacter::Screw), ("Mixed", DislocationCharacter::Mixed)];
const SHAPES: [(&str, DislocationShape); 5] = [
    ("Half-plane", DislocationShape::HalfPlane),
    ("Cylinder", DislocationShape::Cylinder),
    ("Sphere", DislocationShape::Sphere),
    ("Ellipsoid", DislocationShape::Ellipsoid),
    ("Freeform 2D", DislocationShape::Freeform2D),
];

pub struct DislocationDialog {
    pub open: bool,
    pub params: DislocationParams,
    detection: Option<DislocationDetectionResult>,
    running: Option<Receiver<DislocationResult>>,
    pub result: Option<DislocationResult>,
    pub status: (String, bool),
    /// Editing a pipeline step (Edit in Dialog).
    step_edit: Option<StepEdit>,
    /// The step's --input file (kept in its options).
    source_path: String,
    /// In a step edit: the structure the dislocation goes into (the step's
    /// input or its --input file) instead of the view's.
    step_source: Option<Structure>,
    /// Update Step was pressed: the step with its new options.
    pub updated: Option<StepEdit>,
}

impl Default for DislocationDialog {
    fn default() -> Self {
        Self {
            open: false,
            params: DislocationParams::default(),
            detection: None,
            running: None,
            result: None,
            status: (String::new(), false),
            step_edit: None,
            source_path: String::new(),
            step_source: None,
            updated: None,
        }
    }
}

fn ivec(ui: &mut egui::Ui, label: &str, v: &mut IVec3) {
    ui.horizontal(|ui| {
        for c in [&mut v.x, &mut v.y, &mut v.z] {
            ui.add(egui::DragValue::new(c).speed(0.2).range(-6..=6));
        }
        ui.label(label);
    });
}

fn vec3(ui: &mut egui::Ui, label: &str, v: &mut Vec3, speed: f64, range: std::ops::RangeInclusive<f32>) {
    ui.horizontal(|ui| {
        for c in [&mut v.x, &mut v.y, &mut v.z] {
            ui.add(egui::DragValue::new(c).speed(speed).range(range.clone()).fixed_decimals(3));
        }
        ui.label(label);
    });
}

fn number(ui: &mut egui::Ui, label: &str, v: &mut f32, speed: f64, range: std::ops::RangeInclusive<f32>, decimals: usize, hint: &str) {
    ui.horizontal(|ui| {
        ui.add(egui::DragValue::new(v).speed(speed).range(range).fixed_decimals(decimals));
        ui.label(label);
    })
    .response
    .on_hover_text(hint);
}

const CHARACTER_NAMES: [(&str, DislocationCharacter); 3] = [("edge", DislocationCharacter::Edge), ("screw", DislocationCharacter::Screw), ("mixed", DislocationCharacter::Mixed)];
const SHAPE_NAMES: [(&str, DislocationShape); 5] = [
    ("halfplane", DislocationShape::HalfPlane),
    ("cylinder", DislocationShape::Cylinder),
    ("sphere", DislocationShape::Sphere),
    ("ellipsoid", DislocationShape::Ellipsoid),
    ("freeform", DislocationShape::Freeform2D),
];

impl DislocationDialog {
    /// Insert Dislocation from the menu: inserts into the view again (ends a step edit).
    pub fn open_from_menu(&mut self) {
        self.step_edit = None;
        self.step_source = None;
        self.source_path.clear();
        self.open = true;
    }

    /// Opens the dialog on a build-dislocation step (editStep); missing flags
    /// take the defaults of AtomForge --build dislocation.
    pub fn edit_step(&mut self, edit: StepEdit) -> bool {
        let options = edit.reader();
        let ivec3 = |flag: &str, fallback: IVec3| {
            let v = options.numbers(flag);
            if v.len() == 3 {
                IVec3::new(v[0].round() as i32, v[1].round() as i32, v[2].round() as i32)
            } else {
                fallback
            }
        };
        let vec3 = |flag: &str, fallback: Vec3| {
            let v = options.numbers(flag);
            if v.len() == 3 {
                Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32)
            } else {
                fallback
            }
        };
        let mut p = DislocationParams::default();
        let character = options.text("--character", "edge").to_ascii_lowercase();
        p.character = CHARACTER_NAMES.iter().find(|(n, _)| *n == character).map(|(_, c)| *c).unwrap_or(DislocationCharacter::Edge);
        let shape = options.text("--shape", "halfplane").to_ascii_lowercase();
        p.shape = SHAPE_NAMES.iter().find(|(n, _)| *n == shape).map(|(_, s)| *s).unwrap_or(DislocationShape::HalfPlane);
        p.auto_directions = !options.has("--manual-vectors");
        p.plane_hkl = ivec3("--plane", p.plane_hkl);
        p.burgers_uvw = ivec3("--burgers", p.burgers_uvw);
        p.line_uvw = ivec3("--line", p.line_uvw);
        p.line_point_fractional = vec3("--line-frac", p.line_point_fractional);
        p.use_fractional_line_point = !options.has("--line-center");
        p.line_point_cartesian_offset = vec3("--line-offset", p.line_point_cartesian_offset);
        p.burgers_scale = options.number("--bscale", p.burgers_scale as f64) as f32;
        p.burgers_override_magnitude = options.number("--bmag", p.burgers_override_magnitude as f64) as f32;
        p.mixed_character_angle_deg = options.number("--mixed-angle", p.mixed_character_angle_deg as f64) as f32;
        p.poisson_ratio = options.number("--nu", p.poisson_ratio as f64) as f32;
        p.core_radius = options.number("--core", p.core_radius as f64) as f32;
        p.cutoff_radius = options.number("--cutoff", p.cutoff_radius as f64) as f32;
        p.line_half_length = options.number("--line-half", p.line_half_length as f64) as f32;
        p.cylinder_radius = options.number("--cyl-radius", p.cylinder_radius as f64) as f32;
        p.sphere_radius = options.number("--sphere-radius", p.sphere_radius as f64) as f32;
        p.ellipsoid_radii = vec3("--ellipsoid", p.ellipsoid_radii);
        if options.has("--poly2d") {
            // "x1 y1;x2 y2;..."
            let points: Vec<Vec2> = options
                .text("--poly2d", "")
                .split(';')
                .filter_map(|point| {
                    let mut w = point.split_whitespace().map(|v| v.parse::<f32>());
                    match (w.next(), w.next()) {
                        (Some(Ok(x)), Some(Ok(y))) => Some(Vec2::new(x, y)),
                        _ => None,
                    }
                })
                .collect();
            if points.len() >= 3 {
                p.freeform_points = points;
            }
        }
        p.dipole = options.has("--dipole");
        let offset = options.numbers("--dipole-offset");
        if offset.len() == 2 {
            p.dipole_offset = Vec2::new(offset[0] as f32, offset[1] as f32);
        }
        p.anisotropic_elasticity = options.has("--anisotropic");
        p.elastic_symmetry = if options.text("--elastic-symmetry", "cubic").to_ascii_lowercase() == "hexagonal" { ElasticSymmetry::Hexagonal } else { ElasticSymmetry::Cubic };
        p.elastic_c11 = options.number("--elastic-c11", 0.0);
        p.elastic_c12 = options.number("--elastic-c12", 0.0);
        p.elastic_c44 = options.number("--elastic-c44", 0.0);
        p.elastic_c13 = options.number("--elastic-c13", 0.0);
        p.elastic_c33 = options.number("--elastic-c33", 0.0);
        p.elastic_noise_amplitude = options.number("--elastic-noise", p.elastic_noise_amplitude);
        self.params = p;
        // Source: a file named in the options, else the structure entering the step.
        self.detection = None;
        self.result = None;
        self.status = (String::new(), false);
        self.source_path = options.text("--input", "");
        self.step_source = None;
        if !self.source_path.is_empty() {
            match af_io::loader::load_structure_from_file(&self.source_path) {
                Ok(s) => self.step_source = Some(s),
                Err(e) => self.status = (format!("Cannot load {}: {e}", self.source_path), true),
            }
        } else if edit.use_input() && !edit.input.atoms.is_empty() {
            self.step_source = Some(edit.input.clone());
            self.status = (format!("Step input: {} atoms", edit.input.atoms.len()), false);
        }
        self.step_edit = Some(edit);
        self.open = true;
        true
    }

    /// The options of the build-dislocation step (stepOptions).
    pub fn step_options(&self) -> String {
        let p = &self.params;
        let ivec3 = |v: IVec3| [v.x as f64, v.y as f64, v.z as f64];
        let vec3 = |v: Vec3| [v.x as f64, v.y as f64, v.z as f64];
        fn name<T: PartialEq + Copy>(list: &[(&'static str, T)], v: T) -> &'static str {
            list.iter().find(|(_, x)| *x == v).map(|(n, _)| *n).unwrap_or("")
        }
        let mut out = Writer::new();
        if !self.source_path.is_empty() {
            out.add("--input", &self.source_path);
        }
        out.add("--character", name(&CHARACTER_NAMES, p.character));
        out.add("--shape", name(&SHAPE_NAMES, p.shape));
        if !p.auto_directions {
            out.flag("--manual-vectors");
        }
        out.add_numbers("--plane", &ivec3(p.plane_hkl)).add_numbers("--burgers", &ivec3(p.burgers_uvw)).add_numbers("--line", &ivec3(p.line_uvw));
        if p.use_fractional_line_point {
            out.add_numbers("--line-frac", &vec3(p.line_point_fractional));
        } else {
            out.flag("--line-center");
        }
        out.add_numbers("--line-offset", &vec3(p.line_point_cartesian_offset));
        out.add_number("--bscale", p.burgers_scale as f64);
        if p.burgers_override_magnitude > 0.0 {
            out.add_number("--bmag", p.burgers_override_magnitude as f64);
        }
        if p.character == DislocationCharacter::Mixed {
            out.add_number("--mixed-angle", p.mixed_character_angle_deg as f64);
        }
        out.add_number("--nu", p.poisson_ratio as f64).add_number("--core", p.core_radius as f64).add_number("--cutoff", p.cutoff_radius as f64);
        out.add_number("--line-half", p.line_half_length as f64);
        match p.shape {
            DislocationShape::Cylinder => {
                out.add_number("--cyl-radius", p.cylinder_radius as f64);
            }
            DislocationShape::Sphere => {
                out.add_number("--sphere-radius", p.sphere_radius as f64);
            }
            DislocationShape::Ellipsoid => {
                out.add_numbers("--ellipsoid", &vec3(p.ellipsoid_radii));
            }
            DislocationShape::Freeform2D if p.freeform_points.len() >= 3 => {
                let polygon = p.freeform_points.iter().map(|q| format!("{} {}", Writer::number(q.x as f64), Writer::number(q.y as f64))).collect::<Vec<_>>().join(";");
                out.add("--poly2d", &polygon);
            }
            _ => {}
        }
        if p.dipole {
            out.flag("--dipole").add_numbers("--dipole-offset", &[p.dipole_offset.x as f64, p.dipole_offset.y as f64]);
        }
        if p.anisotropic_elasticity {
            let hexagonal = p.elastic_symmetry == ElasticSymmetry::Hexagonal;
            out.flag("--anisotropic").add("--elastic-symmetry", if hexagonal { "hexagonal" } else { "cubic" });
            out.add_number("--elastic-c11", p.elastic_c11).add_number("--elastic-c12", p.elastic_c12).add_number("--elastic-c44", p.elastic_c44);
            if hexagonal {
                out.add_number("--elastic-c13", p.elastic_c13).add_number("--elastic-c33", p.elastic_c33);
            }
            out.add_number("--elastic-noise", p.elastic_noise_amplitude);
        }
        out.str().to_string()
    }

    /// Update Step: hands back the step with these options and ends the edit
    /// (the step input stays the source unless the options name a file).
    pub fn update_step(&mut self) -> Option<StepEdit> {
        let edit = self.step_edit.take()?;
        self.step_source = None;
        self.open = false;
        let mut out = edit.with_options(&self.step_options());
        if self.source_path.is_empty() && !edit.input.atoms.is_empty() {
            out.parameters.set("use_input", true);
        }
        Some(out)
    }

    /// The structure the dialog works on: the step's in a step edit, else the view's.
    pub fn source<'a>(&'a self, view: &'a Structure) -> &'a Structure {
        self.step_source.as_ref().unwrap_or(view)
    }

    /// Sets one option by name (UI scripts): character, shape, burgers-scale,
    /// poisson, core, cutoff, half-length, mixed-angle, plane/burgers/line H K L
    /// (turns auto directions off), dipole on|off, radius R.
    pub fn set(&mut self, name: &str, value: &str) -> Result<(), String> {
        let p = &mut self.params;
        let f = || value.trim().parse::<f32>().map_err(|_| format!("{name} wants a number, got {value:?}"));
        let triple = || -> Result<IVec3, String> {
            let v: Vec<i32> = value.split_whitespace().map(|w| w.parse::<i32>()).collect::<Result<_, _>>().map_err(|_| format!("{name} wants three integers"))?;
            (v.len() == 3).then(|| IVec3::new(v[0], v[1], v[2])).ok_or(format!("{name} wants three integers"))
        };
        match name {
            "character" => p.character = CHARACTERS.iter().find(|(l, _)| l.eq_ignore_ascii_case(value.trim())).ok_or(format!("unknown character {value}"))?.1,
            "shape" => p.shape = SHAPES.iter().find(|(l, _)| l.eq_ignore_ascii_case(value.trim())).ok_or(format!("unknown shape {value}"))?.1,
            "burgers-scale" => p.burgers_scale = f()?,
            "poisson" => p.poisson_ratio = f()?,
            "core" => p.core_radius = f()?,
            "cutoff" => p.cutoff_radius = f()?,
            "half-length" => p.line_half_length = f()?,
            "mixed-angle" => p.mixed_character_angle_deg = f()?,
            "radius" => {
                let r = f()?;
                p.cylinder_radius = r;
                p.sphere_radius = r;
            }
            "plane" | "burgers" | "line" => {
                let v = triple()?;
                p.auto_directions = false;
                *match name {
                    "plane" => &mut p.plane_hkl,
                    "burgers" => &mut p.burgers_uvw,
                    _ => &mut p.line_uvw,
                } = v;
            }
            "dipole" => p.dipole = value.trim() == "on",
            "anisotropic" => p.anisotropic_elasticity = value.trim() == "on",
            _ => return Err(format!("unknown dislocation option {name}")),
        }
        Ok(())
    }

    /// Detects the lattice on this thread (UI scripts).
    pub fn detect_now(&mut self, source: &Structure) -> &DislocationDetectionResult {
        let d = detect_dislocation_lattice(source, true);
        self.status = (d.message.clone(), !d.success);
        self.detection = Some(d);
        self.detection.as_ref().unwrap()
    }

    /// Builds on this thread (UI scripts).
    pub fn generate_now(&mut self, source: &Structure) -> &DislocationResult {
        let r = build_dislocation(source, &self.params);
        self.status = (r.message.clone(), !r.success);
        self.result = Some(r);
        self.result.as_ref().unwrap()
    }

    /// Draws the dialog; returns a structure (with its loop) for Replace Main Scene.
    pub fn show(&mut self, ctx: &egui::Context, view: &Structure) -> Option<Structure> {
        let step_source = self.step_source.clone();
        let source = step_source.as_ref().unwrap_or(view);
        if let Some(rx) = &self.running {
            if let Ok(r) = rx.try_recv() {
                self.status = (r.message.clone(), !r.success);
                self.result = Some(r);
                self.running = None;
            } else {
                ctx.request_repaint_after(std::time::Duration::from_millis(100));
            }
        }
        if !self.open {
            return None;
        }
        let mut replace = None;
        let mut open = self.open;
        let mut close = false;
        // The C++ layout: input and output previews at the left, the options at
        // the right, the status and Close below.
        crate::theme::modal_window("Insert Dislocation").open(&mut open).resizable(false).show(ctx, |ui| {
            ui.set_width(1020.0);
            if let Some(edit) = &self.step_edit {
                crate::structure_dialog::step_banner(ui, &edit.step);
            }
            let validated = self.result.as_ref().is_some_and(|r| r.success);
            let good = egui::Color32::from_rgb(70, 170, 90);
            ui.horizontal(|ui| {
                ui.weak("Isotropic-elastic dislocation insertion with automatic FCC/HCP/BCC detection and shape-controlled application region.");
                if validated {
                    ui.colored_label(good, "Validated");
                } else {
                    ui.weak("Not Generated");
                }
            });
            ui.separator();
            let height = 600.0;
            ui.horizontal_top(|ui| {
                let left = 560.0;
                ui.allocate_ui_with_layout(egui::vec2(left, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                    ui.set_width(left);
                    let half = (height - 2.0 * 34.0 - 12.0) * 0.5;
                    ui.horizontal(|ui| {
                        ui.label("Input Preview");
                        if source.atoms.is_empty() {
                            ui.weak("(open or build a crystal first)");
                        } else {
                            ui.weak(format!("the structure in the view ({} atoms)", source.atoms.len()));
                        }
                    });
                    crate::preview3d::show(ui, egui::vec2(left, half), source, &["Open or build a crystal first:", "the dislocation is inserted into the structure in the view"]);
                    ui.horizontal(|ui| {
                        ui.label("Output Preview");
                        if validated {
                            ui.colored_label(good, "Validated");
                        }
                    });
                    let empty = Structure::default();
                    let output = self.result.as_ref().filter(|r| r.success).map(|r| &r.output).unwrap_or(&empty);
                    crate::preview3d::show(ui, egui::vec2(left, half), output, &["Generate dislocation structure", "to preview and validate"]);
                });
                ui.add_sized([9.0, height], egui::Separator::default().vertical());
                ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), height), egui::Layout::top_down(egui::Align::Min), |ui| {
                    egui::ScrollArea::vertical().max_height(height).show(ui, |ui| self.controls(ui, source, &mut replace));
                });
            });
            ui.separator();
            ui.horizontal(|ui| {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if crate::theme::secondary_button(ui, "Close", Some(egui::vec2(crate::dialog_layout::footer_button_width(), crate::theme::FRAME_HEIGHT))).clicked() {
                        close = true;
                    }
                    if !self.status.0.is_empty() {
                        let colour = if self.status.1 { egui::Color32::from_rgb(220, 80, 60) } else { good };
                        ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                            ui.add(egui::Label::new(egui::RichText::new(&self.status.0).color(colour)).truncate());
                        });
                    }
                });
            });
        });
        if close {
            open = false;
        }
        if !open {
            self.step_edit = None;
            self.step_source = None;
        }
        self.open = open && self.open;
        replace
    }

    fn controls(&mut self, ui: &mut egui::Ui, source: &Structure, replace: &mut Option<Structure>) {
        let has = !source.atoms.is_empty();
        ui.strong("Dislocation Options");
        if self.running.is_some() {
            ui.weak("Generation running in background...");
        }
        if ui.add_enabled(has, egui::Button::new("Detect Lattice").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
            let d = detect_dislocation_lattice(source, true);
            self.status = (d.message.clone(), !d.success);
            self.detection = Some(d);
        }
        match self.detection.as_ref().filter(|d| d.success) {
            Some(d) => {
                ui.label(format!("Detected lattice: {}", dislocation_lattice_family_name(d.family)));
                ui.weak(format!("Recognized atoms: {}", d.recognized_count));
            }
            None => {
                ui.weak("Detected lattice: not identified yet");
            }
        }
        let p = &mut self.params;
        ui.separator();
        ui.weak("Character / Region");
        fn label<T: PartialEq + Copy>(list: &[(&'static str, T)], v: T) -> &'static str {
            list.iter().find(|(_, x)| *x == v).map(|(l, _)| *l).unwrap_or("")
        }
        egui::ComboBox::from_label("Character").selected_text(label(&CHARACTERS, p.character)).show_ui(ui, |ui| {
            for (name, c) in CHARACTERS {
                ui.selectable_value(&mut p.character, c, name);
            }
        });
        egui::ComboBox::from_label("Shape").selected_text(label(&SHAPES, p.shape)).show_ui(ui, |ui| {
            for (name, s) in SHAPES {
                ui.selectable_value(&mut p.shape, s, name);
            }
        });
        ui.separator();
        ui.weak("Crystallographic Directions");
        crate::dialog_widgets::checkbox(ui, &mut p.auto_directions, "Auto lattice vectors").on_hover_text("Derives the construction directions. Disable it to enter integer Plane (hkl), Burgers [uvw], and Line [uvw] triplets. A direction and a plane normal are different objects in a general lattice.");
        if !p.auto_directions {
            ivec(ui, "Plane hkl", &mut p.plane_hkl);
            ivec(ui, "Burgers uvw", &mut p.burgers_uvw);
            ivec(ui, "Line uvw", &mut p.line_uvw);
        }
        crate::dialog_widgets::checkbox(ui, &mut p.use_fractional_line_point, "Line point in fractional coordinates").on_hover_text("Enables Line point frac, a position relative to the cell. Line point offset adds a Cartesian offset in Å.");
        if p.use_fractional_line_point {
            vec3(ui, "Line point frac", &mut p.line_point_fractional, 0.01, -1.0..=2.0);
        }
        vec3(ui, "Line point offset (Å)", &mut p.line_point_cartesian_offset, 0.05, -100.0..=100.0);
        ui.separator();
        ui.weak("Elastic Parameters");
        number(ui, "Burgers scale", &mut p.burgers_scale, 0.02, 0.05..=3.0, 3, "Multiplies the selected Burgers magnitude. A continuously adjustable value can produce a geometric displacement that is not a crystallographically allowed perfect dislocation.");
        ui.add_enabled_ui(!p.anisotropic_elasticity, |ui| {
            number(ui, "Poisson ratio", &mut p.poisson_ratio, 0.005, 0.05..=0.49, 3, "Elastic input to the displacement model (0.05-0.49); not fitted from the structure.");
        });
        number(ui, "Core radius (Å)", &mut p.core_radius, 0.05, 0.05..=20.0, 3, "Regularisation length in Å near the line.");
        number(ui, "Cutoff radius (Å)", &mut p.cutoff_radius, 0.2, 0.0..=100.0, 2, "Radial attenuation in Å; zero turns attenuation off.");
        number(ui, "Line half length (Å)", &mut p.line_half_length, 0.5, 0.1..=1_000_000.0, 1, "Limits the extent along the dislocation line.");
        if p.character == DislocationCharacter::Mixed {
            ui.add(egui::Slider::new(&mut p.mixed_character_angle_deg, 0.0..=90.0).fixed_decimals(1).text("Mixed angle (°)"));
        }
        ui.separator();
        ui.weak("Anisotropic Elasticity (Stroh formalism)");
        crate::dialog_widgets::checkbox(ui, &mut p.anisotropic_elasticity, "Use anisotropic elasticity").on_hover_text("Replaces the isotropic Poisson-ratio displacement field with the Stroh sextic formalism, given the crystal's elastic constants.");
        if p.anisotropic_elasticity {
            ui.horizontal(|ui| {
                ui.label("Elastic symmetry");
                crate::dialog_widgets::radio_value(ui, &mut p.elastic_symmetry, ElasticSymmetry::Cubic, "Cubic");
                crate::dialog_widgets::radio_value(ui, &mut p.elastic_symmetry, ElasticSymmetry::Hexagonal, "Hexagonal");
            });
            let constant = |ui: &mut egui::Ui, label: &str, v: &mut f64| {
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(v).speed(1.0).range(0.0..=2000.0).fixed_decimals(1));
                    ui.label(label);
                });
            };
            constant(ui, "C11 (GPa)", &mut p.elastic_c11);
            constant(ui, "C12 (GPa)", &mut p.elastic_c12);
            constant(ui, "C44 (GPa)", &mut p.elastic_c44);
            if p.elastic_symmetry == ElasticSymmetry::Hexagonal {
                constant(ui, "C13 (GPa)", &mut p.elastic_c13);
                constant(ui, "C33 (GPa)", &mut p.elastic_c33);
                ui.label("The structure's THIRD cell vector must be the c-axis.");
            }
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut p.elastic_noise_amplitude).speed(0.0001).range(0.0..=0.1).fixed_decimals(4));
                ui.label("Noise amplitude");
            });
            ui.label("Breaks the sextic formalism's degeneracy at isotropic or high-symmetry orientations (e.g. line = plane normal).");
            if p.elastic_c11 <= 0.0 || p.elastic_c44 <= 0.0 {
                ui.colored_label(egui::Color32::from_rgb(230, 170, 60), "C11 and C44 must be positive to generate.");
            }
        }
        ui.separator();
        ui.weak("Dipole");
        crate::dialog_widgets::checkbox(ui, &mut p.dipole, "Insert as dipole (adds an opposite-Burgers-vector partner)");
        if p.dipole {
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut p.dipole_offset.x).speed(0.2).range(-500.0..=500.0).fixed_decimals(2));
                ui.add(egui::DragValue::new(&mut p.dipole_offset.y).speed(0.2).range(-500.0..=500.0).fixed_decimals(2));
                ui.label("Dipole offset (Å)");
            });
            ui.label("Offset of the second (opposite-sign) dislocation, in the local slip-plane axes; a dipole's net Burgers vector is zero, so it is periodicity-compatible.");
        }
        ui.separator();
        ui.weak("Shape Parameters");
        match p.shape {
            DislocationShape::Cylinder => number(ui, "Cylinder radius (Å)", &mut p.cylinder_radius, 0.2, 0.1..=200.0, 2, ""),
            DislocationShape::Sphere => number(ui, "Sphere radius (Å)", &mut p.sphere_radius, 0.2, 0.1..=200.0, 2, ""),
            DislocationShape::Ellipsoid => vec3(ui, "Ellipsoid radii (Å)", &mut p.ellipsoid_radii, 0.2, 0.1..=200.0),
            DislocationShape::Freeform2D => {
                ui.weak("Freeform polygon in local (x, y) around the dislocation line.");
                ui.horizontal(|ui| {
                    if ui.button("Add Point").clicked() {
                        let next = p.freeform_points.last().map(|q| *q + Vec2::new(1.0, 0.0)).unwrap_or(Vec2::ZERO);
                        p.freeform_points.push(next);
                    }
                    if ui.button("Remove Last").clicked() {
                        p.freeform_points.pop();
                    }
                });
                for q in p.freeform_points.iter_mut() {
                    ui.horizontal(|ui| {
                        ui.add(egui::DragValue::new(&mut q.x).speed(0.05).range(-500.0..=500.0).fixed_decimals(2));
                        ui.add(egui::DragValue::new(&mut q.y).speed(0.05).range(-500.0..=500.0).fixed_decimals(2));
                    });
                }
            }
            DislocationShape::HalfPlane => {}
        }
        ui.add_space(6.0);
        if self.step_edit.is_some() {
            ui.add_space(6.0);
            if ui.add(egui::Button::new(egui::RichText::new("Update Step").strong()).min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                self.updated = self.update_step();
            }
            return;
        }
        let can = has && self.running.is_none() && (!p.anisotropic_elasticity || (p.elastic_c11 > 0.0 && p.elastic_c44 > 0.0));
        if ui.add_enabled_ui(can, |ui| crate::theme::primary_button(ui, "Generate Dislocation", Some(egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT)))).inner.clicked() {
            let (snapshot, params) = (source.clone(), self.params.clone());
            let (tx, rx) = channel();
            std::thread::spawn(move || {
                let _ = tx.send(build_dislocation(&snapshot, &params));
            });
            self.running = Some(rx);
            self.status = ("Generating dislocation...".into(), false);
        }
        if let Some(r) = self.result.as_ref().filter(|r| r.success) {
            ui.separator();
            ui.label("Result summary");
            let v = &r.validation;
            ui.weak(format!("Shifted atoms: {}", r.shifted_atom_count));
            ui.weak(format!("|b|: {:.4} Å", r.burgers_magnitude));
            ui.weak(format!("RMS displacement: {:.4} Å", v.rms_displacement));
            ui.weak(format!("Max displacement: {:.4} Å", v.max_displacement));
            ui.weak(format!("Min distance: {:.4} Å", v.min_interatomic_distance));
            ui.weak(format!("Family before/after: {} → {}", dislocation_lattice_family_name(v.family_before), dislocation_lattice_family_name(v.family_after)));
            if ui.button(egui::RichText::new("Replace Main Scene").strong()).clicked() {
                *replace = Some(result_structure(r));
                self.status = ("Main scene replaced with generated dislocation structure.".into(), false);
            }
        }
    }
}

/// The generated structure with its loop outline (for the overlay).
pub fn result_structure(r: &DislocationResult) -> Structure {
    let mut s = r.output.clone();
    s.dislocation_loop_points = r.loop_points.iter().map(|p| [p.x, p.y, p.z]).collect();
    s
}
