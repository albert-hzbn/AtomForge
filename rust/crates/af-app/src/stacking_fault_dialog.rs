//! Build > Stacking Fault (the C++ ui/StackingFaultBuilderDialog): a crystal
//! identified as FCC-, HCP- or BCC-like, the fault plane of its family, the
//! cell mode, layers and the displacement steps; the generated sliding
//! sequence is browsed frame by frame, one frame put in the view, or all of
//! them exported to a folder.

use af_algorithms::stacking_fault::{
    build_stacking_fault_sequence, detect_stacking_fault_family, stacking_fault_family_name, StackingFaultCellMode, StackingFaultDetectionResult, StackingFaultFamily, StackingFaultParams,
    StackingFaultPlane, StackingFaultResult,
};
use af_core::structure::Structure;
use eframe::egui;

use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

/// (label, extension, format) of the export formats.
const EXPORT_FORMATS: [(&str, &str, &str); 5] =
    [("XYZ (.xyz)", ".xyz", "xyz"), ("CIF (.cif)", ".cif", "cif"), ("VASP (.vasp)", ".vasp", "vasp"), ("LAMMPS Data (.data)", ".data", "data"), ("Quantum ESPRESSO (.pwi)", ".pwi", "pwscf")];

fn plane_options(family: StackingFaultFamily) -> &'static [(&'static str, StackingFaultPlane)] {
    match family {
        StackingFaultFamily::Fcc => &[("(112) <111>", StackingFaultPlane::Fcc111)],
        StackingFaultFamily::Hcp => &[("Basal", StackingFaultPlane::HcpBasal), ("Prismatic", StackingFaultPlane::HcpPrismatic), ("Pyramidal", StackingFaultPlane::HcpPyramidal)],
        StackingFaultFamily::Bcc => &[("(110) <111>", StackingFaultPlane::Bcc110), ("(112) <111>", StackingFaultPlane::Bcc112)],
        StackingFaultFamily::Unknown => &[("Auto", StackingFaultPlane::Auto)],
    }
}

/// The --plane index (the C++ enum order).
const PLANES: [StackingFaultPlane; 7] = [
    StackingFaultPlane::Auto,
    StackingFaultPlane::Fcc111,
    StackingFaultPlane::HcpBasal,
    StackingFaultPlane::HcpPrismatic,
    StackingFaultPlane::HcpPyramidal,
    StackingFaultPlane::Bcc110,
    StackingFaultPlane::Bcc112,
];

fn cell_mode_label(mode: StackingFaultCellMode) -> &'static str {
    match mode {
        StackingFaultCellMode::OrthogonalCell => "Orthogonal Cell",
        StackingFaultCellMode::SmallestUnitCell => "Smallest Unit Cell",
    }
}

fn base_name(path: &str) -> String {
    path.rsplit(['/', '\\']).next().unwrap_or(path).to_string()
}

fn stem(path: &str) -> String {
    let base = base_name(path);
    match base.rfind('.') {
        Some(i) => base[..i].to_string(),
        None => base,
    }
}

fn parent(path: &str) -> String {
    std::path::Path::new(path).parent().map(|p| p.to_string_lossy().replace('\\', "/")).filter(|p| !p.is_empty()).unwrap_or_else(|| ".".into())
}

/// Letters, digits, - and _ kept, spaces made _ (the C++ sanitizeName).
fn sanitize(value: &str) -> String {
    let out: String = value
        .chars()
        .filter_map(|c| match c {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '-' | '_' => Some(c),
            ' ' => Some('_'),
            _ => None,
        })
        .collect();
    if out.is_empty() { "stacking_fault".into() } else { out }
}

/// A small orthographic picture of a structure (atoms back to front), in
/// place of the C++ offscreen previews.
fn mini_view(ui: &mut egui::Ui, s: Option<&Structure>, size: egui::Vec2, hint: (&str, &str)) {
    let (rect, _) = ui.allocate_exact_size(size, egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let visuals = ui.visuals();
    painter.rect_filled(rect, 3.0, visuals.extreme_bg_color);
    painter.rect_stroke(rect, 3.0, visuals.widgets.noninteractive.bg_stroke, egui::StrokeKind::Inside);
    let Some(s) = s.filter(|s| !s.atoms.is_empty()) else {
        let c = rect.center();
        painter.text(c - egui::vec2(0.0, 9.0), egui::Align2::CENTER_CENTER, hint.0, egui::FontId::proportional(13.0), visuals.weak_text_color());
        painter.text(c + egui::vec2(0.0, 9.0), egui::Align2::CENTER_CENTER, hint.1, egui::FontId::proportional(13.0), visuals.weak_text_color());
        return;
    };
    let (yaw, pitch) = (45f64.to_radians(), 35f64.to_radians());
    let project = |p: [f64; 3]| {
        let x = p[0] * yaw.cos() - p[1] * yaw.sin();
        let y0 = p[0] * yaw.sin() + p[1] * yaw.cos();
        let y = p[2] * pitch.cos() - y0 * pitch.sin();
        let depth = p[2] * pitch.sin() + y0 * pitch.cos();
        (x, y, depth)
    };
    let points: Vec<(f64, f64, f64)> = s.atoms.iter().map(|a| project([a.x, a.y, a.z])).collect();
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for p in &points {
        lo = [lo[0].min(p.0), lo[1].min(p.1)];
        hi = [hi[0].max(p.0), hi[1].max(p.1)];
    }
    let span = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1.0);
    // Atoms 0.5 A in radius, the picture fitted with a margin of one atom.
    let fit = rect.width().min(rect.height()) as f64 - 16.0;
    let scale = (fit / (span + 1.0)).max(0.1);
    let centre = [(lo[0] + hi[0]) * 0.5, (lo[1] + hi[1]) * 0.5];
    let radius = (scale * 0.5).clamp(1.5, 30.0) as f32;
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by(|&a, &b| points[a].2.partial_cmp(&points[b].2).unwrap_or(std::cmp::Ordering::Equal));
    for i in order {
        let (x, y, _) = points[i];
        let at = rect.center() + egui::vec2(((x - centre[0]) * scale) as f32, (-(y - centre[1]) * scale) as f32);
        let a = &s.atoms[i];
        let colour = crate::theme::srgb(a.r, a.g, a.b);
        painter.circle(at, radius, colour, egui::Stroke::new(0.8, egui::Color32::from_black_alpha(140)));
    }
}

pub struct StackingFaultDialog {
    open: bool,
    take_view: bool,
    /// Editing a pipeline step: Update Step instead of Use Selected Structure.
    step: Option<StepEdit>,
    /// The step's --frame while no sequence is generated, and its --sequence folder (kept as is).
    step_frame: usize,
    step_sequence: String,
    source: Structure,
    source_loaded: bool,
    source_path: String,
    source_label: String,
    path_text: String,
    detection: StackingFaultDetectionResult,
    params: StackingFaultParams,
    result: StackingFaultResult,
    selected: usize,
    status: String,
    status_error: bool,
    export_dir: String,
    export_folder: String,
    export_format: usize,
}

impl Default for StackingFaultDialog {
    fn default() -> Self {
        Self {
            open: false,
            take_view: false,
            step: None,
            step_frame: 0,
            step_sequence: String::new(),
            source: Structure::default(),
            source_loaded: false,
            source_path: String::new(),
            source_label: String::new(),
            path_text: String::new(),
            detection: StackingFaultDetectionResult::default(),
            params: StackingFaultParams::default(),
            result: StackingFaultResult::default(),
            selected: 0,
            status: String::new(),
            status_error: false,
            export_dir: ".".into(),
            export_folder: "scene_stacking_fault".into(),
            export_format: 0,
        }
    }
}

impl StackingFaultDialog {
    fn set_status(&mut self, text: impl Into<String>, error: bool) {
        self.status = text.into();
        self.status_error = error;
    }

    fn load_file(&mut self, path: &str) -> bool {
        let loaded = match af_io::loader::load_structure_from_file(path) {
            Ok(s) => s,
            Err(e) => {
                self.set_status(format!("Load failed: {}", if e.is_empty() { path.to_string() } else { e }), true);
                return false;
            }
        };
        if !loaded.has_unit_cell {
            self.set_status("Stacking faults builder requires a structure with a unit cell.", true);
            return false;
        }
        self.source = loaded;
        self.source_loaded = true;
        self.source_path = path.into();
        self.source_label = base_name(path);
        self.export_dir = parent(path);
        self.export_folder = format!("{}_stacking_fault", sanitize(&stem(path)));
        self.analyze();
        self.result = StackingFaultResult::default();
        self.selected = 0;
        true
    }

    fn load_scene(&mut self, scene: &Structure) {
        if scene.atoms.is_empty() {
            self.set_status("Current scene has no atoms.", true);
            return;
        }
        if !scene.has_unit_cell {
            self.set_status("Current scene has no unit cell.", true);
            return;
        }
        self.source = scene.clone();
        self.source_loaded = true;
        self.source_path.clear();
        self.source_label = "scene".into();
        self.export_dir = ".".into();
        self.export_folder = "scene_stacking_fault".into();
        self.analyze();
        self.result = StackingFaultResult::default();
        self.selected = 0;
    }

    fn analyze(&mut self) {
        self.detection = detect_stacking_fault_family(&self.source, true);
        let (m, ok) = (self.detection.message.clone(), self.detection.success);
        self.set_status(m, !ok);
        if ok {
            self.params.plane = plane_options(self.detection.family)[0].1;
        }
    }

    fn generate(&mut self) {
        if !self.source_loaded {
            self.set_status("Load a structure first.", true);
            return;
        }
        self.result = build_stacking_fault_sequence(&self.source, &self.params);
        let (m, ok) = (self.result.message.clone(), self.result.success);
        self.set_status(m, !ok);
        self.selected = 0;
    }

    fn has_sequence(&self) -> bool {
        self.result.success && !self.result.sequence.is_empty()
    }

    fn use_selected(&mut self) -> Option<Outcome> {
        if !self.has_sequence() {
            return None;
        }
        let item = &self.result.sequence[self.selected.min(self.result.sequence.len() - 1)];
        let (label, s) = (item.label.clone(), item.structure.clone());
        self.set_status(format!("Applied {label} to the main scene."), false);
        Some(Outcome::Built(s, format!("Applied stacking fault {label}")))
    }

    fn export(&mut self) -> Result<String, String> {
        if !self.has_sequence() {
            self.set_status("Generate a stacking-fault sequence before exporting.", true);
            return Err(self.status.clone());
        }
        let (_, extension, format) = EXPORT_FORMATS[self.export_format];
        let target = std::path::Path::new(&self.export_dir.replace('\\', "/")).join(sanitize(&self.export_folder));
        if let Err(e) = std::fs::create_dir_all(&target) {
            self.set_status(format!("Could not create export folder: {e}"), true);
            return Err(self.status.clone());
        }
        for (i, item) in self.result.sequence.iter().enumerate() {
            let filename = format!("sf_{i:03}_{}{extension}", sanitize(&item.label));
            let path = target.join(&filename).to_string_lossy().replace('\\', "/");
            if af_io::loader::save_structure(&item.structure, &path, format).is_err() {
                self.set_status(format!("Failed to export: {filename}"), true);
                return Err(self.status.clone());
            }
        }
        let text = format!("Exported {} structures to {}", self.result.sequence.len(), target.to_string_lossy().replace('\\', "/"));
        self.set_status(text.clone(), false);
        Ok(text)
    }

    /// The dialog's settings as the step's options (the flags of AtomForge
    /// --build stacking-fault; the C++ stepOptions).
    pub fn step_options(&self) -> String {
        use af_pipeline::options::Writer;
        let mut out = Writer::new();
        if !self.source_path.is_empty() {
            out.add("--input", &self.source_path);
        }
        let plane = PLANES.iter().position(|p| *p == self.params.plane).unwrap_or(0);
        out.add_integer("--plane", plane as i64);
        out.add_integer("--layers", self.params.layer_count as i64);
        out.add_number("--interval", self.params.interval as f64);
        out.add_number("--maximum", self.params.max_displacement_factor as f64);
        out.add_integer("--frame", if self.has_sequence() { self.selected } else { self.step_frame } as i64);
        if self.params.cell_mode == StackingFaultCellMode::OrthogonalCell {
            out.flag("--orthogonal");
        }
        if !self.step_sequence.is_empty() {
            out.add("--sequence", &self.step_sequence);
        }
        out.str().to_string()
    }

    /// Fills the dialog from a pipeline step (the C++ editStep and the step
    /// setup its drawDialog does once the source is loaded): missing flags take
    /// the defaults of AtomForge --build stacking-fault.
    fn fill_from_step(&mut self, edit: &StepEdit) {
        let o = edit.reader();
        self.params = StackingFaultParams::default();
        let step_plane = PLANES[o.integer("--plane", 0).clamp(0, 6) as usize];
        self.params.plane = step_plane;
        self.params.layer_count = o.integer("--layers", 9);
        self.params.interval = o.number("--interval", 0.1) as f32;
        self.params.max_displacement_factor = o.number("--maximum", 2.0) as f32;
        self.params.cell_mode = if o.has("--orthogonal") { StackingFaultCellMode::OrthogonalCell } else { StackingFaultCellMode::SmallestUnitCell };
        self.step_frame = o.integer("--frame", 0).max(0) as usize;
        self.step_sequence = o.text("--sequence", "");
        // Source: a file named in the options, else the structure entering the step.
        self.source = Structure::default();
        self.source_loaded = false;
        self.source_label.clear();
        self.detection = StackingFaultDetectionResult::default();
        self.result = StackingFaultResult::default();
        self.selected = 0;
        self.set_status("", false);
        self.source_path = o.text("--input", "");
        if !self.source_path.is_empty() {
            let path = self.source_path.clone();
            if !self.load_file(&path) {
                self.source_path = path;
            }
        } else if !edit.input.atoms.is_empty() {
            self.load_scene(&edit.input);
            if self.source_loaded {
                self.source_label = "step input".into();
            }
        }
        // Its plane (detection picks one) and its sequence, at its frame.
        self.params.plane = step_plane;
        if self.source_loaded {
            self.generate();
            if self.has_sequence() {
                self.selected = self.step_frame.min(self.result.sequence.len() - 1);
            }
        }
    }

    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        self.open = false;
        Some(Outcome::UpdateStep(edit.with_options(&self.step_options())))
    }

    fn source_panel(&mut self, ui: &mut egui::Ui, input: &Input, height: f32) {
        ui.horizontal(|ui| {
            ui.label("Input Structure");
            if let Some(edit) = &self.step {
                if ui.small_button("Use Step Input").clicked() {
                    let structure = edit.input.clone();
                    self.load_scene(&structure);
                    if self.source_loaded && self.source_path.is_empty() {
                        self.source_label = "step input".into();
                    }
                }
            } else if ui.small_button("Use Current Scene").clicked() {
                self.load_scene(input.source);
            }
        });
        if self.source_loaded {
            ui.horizontal(|ui| {
                ui.colored_label(egui::Color32::from_rgb(70, 170, 90), &self.source_label);
                ui.weak(format!("({} atoms)", self.source.atoms.len()));
            });
        } else {
            ui.weak("Load a structure file or use the current scene.");
        }
        ui.horizontal(|ui| {
            ui.add(egui::TextEdit::singleline(&mut self.path_text).hint_text("Structure file path").desired_width(ui.available_width() - 60.0));
            if ui.add_enabled(!self.path_text.trim().is_empty(), egui::Button::new("Load")).clicked() {
                let path = self.path_text.trim().to_string();
                self.load_file(&path);
            }
        });
        let w = ui.available_width();
        let h = (height - 90.0).max(120.0);
        mini_view(ui, self.source_loaded.then_some(&self.source), egui::vec2(w, h), ("Load a structure file", "or use the current scene"));
    }

    fn controls(&mut self, ui: &mut egui::Ui) -> Option<Outcome> {
        let mut out = None;
        crate::dialog_layout::section(ui, "Builder");
        if self.detection.success {
            ui.label(format!("Detected phase: {}", stacking_fault_family_name(self.detection.family)));
            crate::dialog_layout::muted(ui, format!("Recognized atoms: {}", self.detection.recognized_count));
        } else {
            crate::dialog_layout::muted(ui, "Detected phase: not identified yet");
        }
        let full = ui.available_width();
        if ui.add_sized([full, crate::theme::FRAME_HEIGHT], egui::Button::new("Detect Structure")).clicked() {
            if self.source_loaded {
                self.analyze();
            } else {
                self.set_status("Load a structure first.", true);
            }
        }
        ui.add_space(4.0);
        // ImGui: the widget fills the column, its label after it.
        let label = 125.0f32.min(full * 0.4);
        let w = full - label - ui.spacing().item_spacing.x;
        let options = plane_options(self.detection.family);
        let current = options.iter().find(|o| o.1 == self.params.plane).map(|o| o.0).unwrap_or(options[0].0);
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("sf plane").width(w - 8.0).selected_text(current).show_ui(ui, |ui| {
                for (label, plane) in options {
                    ui.selectable_value(&mut self.params.plane, *plane, *label);
                }
            });
            ui.label("Plane");
        });
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("sf cell mode").width(w - 8.0).selected_text(cell_mode_label(self.params.cell_mode)).show_ui(ui, |ui| {
                for mode in [StackingFaultCellMode::SmallestUnitCell, StackingFaultCellMode::OrthogonalCell] {
                    ui.selectable_value(&mut self.params.cell_mode, mode, cell_mode_label(mode));
                }
            });
            ui.label("Cell Mode");
        });
        ui.horizontal(|ui| {
            let mut layers = self.params.layer_count as i32;
            crate::dialog_layout::slider_i32(ui, &mut layers, 2..=96, w, |v| v.to_string());
            self.params.layer_count = layers as _;
            ui.label("Layers");
        });
        ui.horizontal(|ui| {
            ui.add_sized([w, crate::theme::FRAME_HEIGHT], egui::DragValue::new(&mut self.params.interval).speed(0.01).range(0.01..=1.0).fixed_decimals(2));
            ui.label("Interval");
        });
        ui.horizontal(|ui| {
            ui.add_sized([w, crate::theme::FRAME_HEIGHT], egui::DragValue::new(&mut self.params.max_displacement_factor).speed(0.05).range(0.10..=4.0).fixed_decimals(2));
            ui.label("Max Displacement");
        });
        ui.add_space(4.0);
        if crate::dialog_layout::wide_primary(ui, "Generate Sequence", crate::theme::FRAME_HEIGHT, true).clicked() {
            self.generate();
        }
        // Without a sequence the step keeps its frame.
        if self.step.is_some() && !self.has_sequence() && ui.add(egui::Button::new("Update Step").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
            out = self.update_step();
        }
        if self.has_sequence() {
            crate::dialog_layout::section(ui, "Sequence");
            let last = self.result.sequence.len() - 1;
            self.selected = self.selected.min(last);
            ui.add(egui::Slider::new(&mut self.selected, 0..=last).text("Structure"));
            ui.label(format!("Label: {}", self.result.sequence[self.selected].label));
            ui.weak(format!("Partial displacement: {:.4} Å", self.result.partial_displacement));
            ui.weak(format!("Shifted atoms: {}", self.result.shifted_atom_count));
            if self.step.is_some() {
                if ui.add(egui::Button::new("Update Step").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                    out = self.update_step();
                }
            } else if ui.add(egui::Button::new("Use Selected Structure").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                out = self.use_selected();
                self.open = false;
            }
            crate::dialog_layout::section(ui, "Export");
            let field = (ui.available_width() - 80.0).max(120.0);
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut self.export_dir).desired_width(field));
                ui.label("Folder");
            });
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut self.export_folder).desired_width(field));
                ui.label("Subfolder");
            });
            if self.source_loaded && !self.source_path.is_empty() && ui.small_button("Use Source Folder").clicked() {
                self.export_dir = parent(&self.source_path);
            }
            egui::ComboBox::from_label("Format").selected_text(EXPORT_FORMATS[self.export_format].0).show_ui(ui, |ui| {
                for (i, f) in EXPORT_FORMATS.iter().enumerate() {
                    ui.selectable_value(&mut self.export_format, i, f.0);
                }
            });
            if ui.add(egui::Button::new("Export All Structures").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                let _ = self.export();
            }
        }
        out
    }

    fn output_panel(&mut self, ui: &mut egui::Ui, height: f32) {
        ui.horizontal(|ui| {
            ui.label("Generated Structures");
            if self.has_sequence() {
                ui.weak(&self.result.sequence[self.selected.min(self.result.sequence.len() - 1)].label);
            }
        });
        let shown = self.has_sequence().then(|| &self.result.sequence[self.selected.min(self.result.sequence.len() - 1)].structure);
        let w = ui.available_width();
        mini_view(ui, shown, egui::vec2(w, (height - 30.0).max(120.0)), ("Generate a stacking-fault sequence", "to preview intermediate structures"));
    }
}

impl StructureDialog for StackingFaultDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        self.load_file(path);
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-stacking-fault"]
    }

    fn title(&self) -> &'static str {
        "Stacking Faults"
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        self.fill_from_step(&edit);
        self.step = Some(edit);
        self.take_view = false;
        self.open = true;
        true
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
        self.step = None; // from the menu: build normally
        if !self.source_loaded {
            self.take_view = true;
        }
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        if std::mem::take(&mut self.take_view) && input.source.has_unit_cell && !input.source.atoms.is_empty() {
            self.load_scene(input.source);
        }
        let mut out = None;
        let mut open = self.open;
        crate::theme::modal_window(self.title()).open(&mut open).collapsible(false).default_size([1220.0, 760.0]).min_size([900.0, 560.0]).show(ctx, |ui| {
            if let Some(edit) = &self.step {
                ui.colored_label(
                    egui::Color32::from_rgb(90, 150, 220),
                    format!("Editing the pipeline step {}. Update step writes these settings into the pipeline; the structure changes through the pipeline.", edit.step),
                );
            }
            ui.weak("Load a crystal, identify whether it is FCC-like, HCP-like, or BCC-like, and build a sliding sequence across the detected stacking-fault family.");
            ui.separator();
            let height = (ui.ctx().content_rect().height() - 300.0).clamp(380.0, 640.0);
            let total = ui.available_width();
            ui.horizontal_top(|ui| {
                ui.allocate_ui_with_layout(egui::vec2(total * 0.34, height), egui::Layout::top_down(egui::Align::Min), |ui| self.source_panel(ui, input, height));
                ui.allocate_ui_with_layout(egui::vec2(total * 0.31, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                    egui::ScrollArea::vertical().id_salt("sf controls").max_height(height).show(ui, |ui| out = self.controls(ui));
                });
                ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), height), egui::Layout::top_down(egui::Align::Min), |ui| self.output_panel(ui, height));
            });
            ui.separator();
            ui.horizontal(|ui| {
                if !self.status.is_empty() {
                    let colour = if self.status_error { egui::Color32::from_rgb(220, 80, 60) } else { egui::Color32::from_rgb(70, 170, 90) };
                    ui.colored_label(colour, &self.status);
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.add(egui::Button::new("Close").min_size(egui::vec2(110.0, 0.0))).clicked() {
                        self.open = false;
                    }
                });
            });
        });
        self.open = self.open && open;
        if !self.open {
            self.step = None;
        }
        out
    }

    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let a = args.trim();
        let num = || a.parse::<f64>().map_err(|_| format!("{command}: {a:?} is not a number"));
        let failed = |d: &Self| if d.status_error { Err(d.status.clone()) } else { Ok(None) };
        match command {
            "close" => {
                self.open = false;
                self.step = None;
            }
            "update-step" => return self.update_step().map(Some).ok_or("not editing a pipeline step".into()),
            "expect-options" => {
                if self.step_options() != a {
                    return Err(format!("options: {}", self.step_options()));
                }
            }
            "use-view" => {
                self.load_scene(input.source);
                return failed(self);
            }
            "load" => {
                self.load_file(a);
                return failed(self);
            }
            "detect" => {
                if !self.source_loaded {
                    return Err("Load a structure first.".into());
                }
                self.analyze();
                return failed(self);
            }
            "plane" => {
                let options = plane_options(self.detection.family);
                self.params.plane = options.iter().find(|o| o.0.eq_ignore_ascii_case(a)).ok_or(format!("plane {a:?} is not one of {:?}", options.iter().map(|o| o.0).collect::<Vec<_>>()))?.1;
            }
            "cell" => self.params.cell_mode = if a == "orthogonal" { StackingFaultCellMode::OrthogonalCell } else { StackingFaultCellMode::SmallestUnitCell },
            "layers" => self.params.layer_count = (num()? as i32).clamp(2, 96),
            "interval" => self.params.interval = num()? as f32,
            "maximum" => self.params.max_displacement_factor = num()? as f32,
            "generate" => {
                self.generate();
                return failed(self);
            }
            "frame" => {
                let i = num()? as usize;
                if i >= self.result.sequence.len() {
                    return Err(format!("no frame {i} of {}", self.result.sequence.len()));
                }
                self.selected = i;
            }
            "use" if self.step.is_some() => return Err("editing a pipeline step: use update-step".into()),
            "use" => return self.use_selected().map(Some).ok_or("Generate a stacking-fault sequence first.".into()),
            "export" => {
                // export DIR [FORMAT INDEX]
                let w: Vec<&str> = a.split_whitespace().collect();
                if let Some(dir) = w.first() {
                    self.export_dir = dir.to_string();
                }
                if let Some(f) = w.get(1) {
                    self.export_format = f.parse::<usize>().map_err(|_| "export DIR [FORMAT]")?.min(EXPORT_FORMATS.len() - 1);
                }
                return self.export().map(|t| Some(Outcome::Status(t)));
            }
            "expect-family" => {
                let got = stacking_fault_family_name(self.detection.family);
                if !got.eq_ignore_ascii_case(a) {
                    return Err(format!("detected {got}"));
                }
            }
            "expect-sequence" => {
                // expect-sequence FRAMES [SHIFTED]
                let w: Vec<usize> = a.split_whitespace().filter_map(|w| w.parse().ok()).collect();
                let r = &self.result;
                if !r.success || w.first() != Some(&r.sequence.len()) || w.get(1).is_some_and(|&n| n != r.shifted_atom_count as usize) {
                    return Err(format!("{} frames, {} shifted ({})", r.sequence.len(), r.shifted_atom_count, r.message));
                }
            }
            _ => return Err(format!("unknown stacking fault dialog command {command}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edit(options: &str) -> StepEdit {
        let mut parameters = af_core::json::Json::object();
        parameters.set("options", options);
        StepEdit { index: 0, step: "build-stacking-fault".into(), parameters, input: Structure::default(), ..Default::default() }
    }

    #[test]
    fn step_options_round_trip() {
        // The C++ stepOptions for the same state (no source: no sequence, the frame kept).
        for text in [
            "--plane 0 --layers 9 --interval 0.1000000015 --maximum 2 --frame 0",
            "--input crystal.cif --plane 6 --layers 12 --interval 0.25 --maximum 1.5 --frame 3 --orthogonal --sequence out/sf",
        ] {
            let mut d = StackingFaultDialog::default();
            assert!(d.edit_step(edit(text)));
            assert_eq!(d.step_options(), text);
        }
        let mut d = StackingFaultDialog::default();
        d.edit_step(edit(""));
        assert_eq!(d.step_options(), "--plane 0 --layers 9 --interval 0.1000000015 --maximum 2 --frame 0");
    }
}
