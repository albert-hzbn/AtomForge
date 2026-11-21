//! Build > Shape Modeller (CAD) (the C++ ui/ShapeModellerDialog): the open
//! crystal carved into a shape scripted as in CAD (primitives, booleans,
//! transforms, patterns, regions), with templates, an insert menu, live
//! checking with a volume and atom-count estimate, a section view through the
//! shape and OBJ export of its surface. The result is the build-shape step.

use af_algorithms::shape::{commands, save_obj, templates, Model};
use af_core::elements;
use af_core::json::Json;
use af_core::structure::Structure;
use af_pipeline::pipeline::{find_modifier_type, make_modifier, PipelineData};
use eframe::egui;
use glam::{DMat3, DVec3};

use crate::params_form::{parameter_form, FormTexts};
use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

const STEP: &str = "build-shape";
const GOOD: egui::Color32 = egui::Color32::from_rgb(70, 170, 90);
const BAD: egui::Color32 = egui::Color32::from_rgb(220, 80, 60);
const WARN: egui::Color32 = egui::Color32::from_rgb(230, 170, 60);

pub struct ShapeDialog {
    open: bool,
    script: String,
    parameters: Json,
    texts: FormTexts,
    /// The parsed script and the cell it used.
    parsed_key: Option<(String, [f64; 9])>,
    model: Option<Model>,
    parse_error: String,
    volume: f64,
    /// Section view: plane (0 xy, 1 xz, 2 yz) and depth (0..1 of the extent).
    plane: usize,
    section_at: f32,
    raster: Vec<egui::Color32>,
    raster_key: String,
    mesh_path: String,
    mesh_browser: crate::file_browser::FileBrowser,
    message: (String, bool),
    /// The build-shape step being edited ("Update Step" instead of Carve).
    step: Option<StepEdit>,
}

impl Default for ShapeDialog {
    fn default() -> Self {
        Self {
            open: false,
            script: templates()[1].1.to_string(),
            parameters: make_modifier(STEP).map(|m| m.parameters).unwrap_or_else(|_| Json::object()),
            texts: FormTexts::default(),
            parsed_key: None,
            model: None,
            parse_error: String::new(),
            volume: 0.0,
            plane: 1,
            section_at: 0.5,
            raster: Vec::new(),
            raster_key: String::new(),
            mesh_path: String::new(),
            mesh_browser: crate::file_browser::FileBrowser::default(),
            message: (String::new(), false),
            step: None,
        }
    }
}

fn element_colour(symbol: &str, alpha: f32) -> egui::Color32 {
    let z = elements::atomic_number(symbol);
    let c = if z > 0 { elements::default_colour(z) } else { [0.62, 0.66, 0.72] };
    crate::theme::srgba(c[0], c[1], c[2], alpha).into()
}

fn cell_volume(s: &Structure) -> f64 {
    if !s.has_unit_cell {
        return 0.0;
    }
    let c = &s.cell_vectors;
    (c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0])).abs()
}

impl ShapeDialog {
    fn mode(&self) -> String {
        self.parameters.find("mode").and_then(|m| m.string().ok()).unwrap_or("carve").to_string()
    }

    /// Operations take the shape after them (or a { } block).
    fn insert(&mut self, snippet: &str) {
        if !self.script.is_empty() && !self.script.ends_with('\n') {
            self.script.push('\n');
        }
        self.script += snippet;
        self.script += if snippet.contains('{') || snippet.contains(';') { "\n" } else { " " };
    }

    fn reparse(&mut self, s: &Structure) {
        let mut cell = [0.0; 9];
        if s.has_unit_cell {
            for r in 0..3 {
                for c in 0..3 {
                    cell[3 * r + c] = s.cell_vectors[r][c];
                }
            }
        }
        if self.parsed_key.as_ref().is_some_and(|(t, c)| *t == self.script && *c == cell) {
            return;
        }
        self.parsed_key = Some((self.script.clone(), cell));
        self.model = None;
        self.parse_error.clear();
        self.volume = 0.0;
        self.raster_key.clear();
        let v = |r: usize| DVec3::new(s.cell_vectors[r][0], s.cell_vectors[r][1], s.cell_vectors[r][2]);
        let lattice = DMat3::from_cols(v(0), v(1), v(2));
        match Model::parse(&self.script, if s.has_unit_cell { Some(&lattice) } else { None }) {
            Ok(m) => {
                let b = m.bounds();
                if b.finite {
                    let e = b.high - b.low;
                    let spacing = 0.25f64.max((e.x * e.y * e.z).cbrt() / 60.0);
                    self.volume = m.volume(spacing);
                }
                self.model = Some(m);
            }
            Err(e) => self.parse_error = e.to_string(),
        }
    }

    /// The status line under the script.
    fn summary(&self, crystal: &Structure) -> Option<(String, egui::Color32)> {
        if !self.parse_error.is_empty() {
            return Some((self.parse_error.clone(), BAD));
        }
        let m = self.model.as_ref()?;
        let b = m.bounds();
        if !b.finite {
            return Some(("The shape is unbounded (a half-space or slab not intersected with a finite shape): it can delete or select atoms, but carving needs a finite shape.".into(), WARN));
        }
        let e = b.high - b.low;
        let n = m.primitives();
        let plural = if n == 1 { "" } else { "s" };
        let volume = cell_volume(crystal);
        let text = if volume > 0.0 && !crystal.atoms.is_empty() {
            format!("{n} primitive{plural}; {:.1} x {:.1} x {:.1} A; volume {:.0} A^3; about {:.0} atoms of this crystal", e.x, e.y, e.z, self.volume, self.volume / volume * crystal.atoms.len() as f64)
        } else {
            format!("{n} primitive{plural}; {:.1} x {:.1} x {:.1} A; volume {:.0} A^3", e.x, e.y, e.z, self.volume)
        };
        Some((text, GOOD))
    }

    fn section(&mut self, ui: &mut egui::Ui, crystal: &Structure, size: f32) {
        let (rect, _) = ui.allocate_exact_size(egui::vec2(size, size), egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 3.0, crate::theme::palette(ui.visuals().dark_mode).surface);
        let Some(model) = &self.model else {
            painter.text(rect.min + egui::vec2(8.0, 8.0), egui::Align2::LEFT_TOP, "No valid shape", egui::FontId::proportional(crate::theme::BASE_FONT), ui.visuals().weak_text_color());
            return;
        };
        let b = model.bounds();
        let (centre, extent) = if b.finite { (0.5 * (b.low + b.high), b.high - b.low) } else { (DVec3::ZERO, DVec3::splat(60.0)) };
        const AXES: [[usize; 3]; 3] = [[0, 1, 2], [0, 2, 1], [1, 2, 0]];
        let [u, v, w] = AXES[self.plane];
        let span = 1.1 * extent[u].max(extent[v]).max(1.0);
        let depth = centre[w] + (self.section_at as f64 - 0.5) * extent[w];
        let n = 96usize;
        let key = format!("{} {:.4} {}{}", self.plane, depth, n, self.script);
        if key != self.raster_key {
            self.raster_key = key;
            self.raster = vec![egui::Color32::TRANSPARENT; n * n];
            let body = crystal.atoms.first().map(|a| a.symbol.clone()).unwrap_or_default();
            for j in 0..n {
                for i in 0..n {
                    let mut q = DVec3::ZERO;
                    q[u] = centre[u] + ((i as f64 + 0.5) / n as f64 - 0.5) * span;
                    q[v] = centre[v] + (0.5 - (j as f64 + 0.5) / n as f64) * span;
                    q[w] = depth;
                    if !model.contains(q) {
                        continue;
                    }
                    self.raster[j * n + i] = match model.material(q) {
                        Some("none") => crate::theme::srgba(0.9, 0.2, 0.2, 0.35).into(),
                        Some(region) => element_colour(region, 1.0),
                        None => element_colour(&body, 1.0),
                    };
                }
            }
        }
        let cell = size / n as f32;
        for j in 0..n {
            let mut i = 0;
            while i < n {
                let colour = self.raster[j * n + i];
                let mut run = i + 1;
                while run < n && self.raster[j * n + run] == colour {
                    run += 1;
                }
                if colour != egui::Color32::TRANSPARENT {
                    painter.rect_filled(egui::Rect::from_min_max(rect.min + egui::vec2(i as f32 * cell, j as f32 * cell), rect.min + egui::vec2(run as f32 * cell, (j + 1) as f32 * cell)), 0.0, colour);
                }
                i = run;
            }
        }
        let to_screen = |a: f64, c: f64| rect.min + egui::vec2(((a - centre[u]) / span + 0.5) as f32 * size, (0.5 - (c - centre[v]) / span) as f32 * size);
        let o = to_screen(0.0, 0.0);
        let axis = egui::Stroke::new(1.0, ui.visuals().weak_text_color());
        painter.line_segment([o - egui::vec2(6.0, 0.0), o + egui::vec2(6.0, 0.0)], axis);
        painter.line_segment([o - egui::vec2(0.0, 6.0), o + egui::vec2(0.0, 6.0)], axis);
        let names = ['x', 'y', 'z'];
        let text = ui.visuals().text_color();
        painter.text(rect.min + egui::vec2(6.0, 4.0), egui::Align2::LEFT_TOP, format!("{}-{} section at {} = {:.1} Å", names[u], names[v], names[w], depth), egui::FontId::proportional(crate::theme::BASE_FONT), text);
        let mut bar = 10f64.powf((span / 3.0).log10().floor());
        if span / bar > 8.0 {
            bar *= 2.0;
        }
        let length = (bar / span) as f32 * size;
        let b0 = egui::pos2(rect.max.x - length - 8.0, rect.max.y - 10.0);
        painter.line_segment([b0, b0 + egui::vec2(length, 0.0)], egui::Stroke::new(2.0, text));
        painter.text(b0 - egui::vec2(0.0, 2.0), egui::Align2::LEFT_BOTTOM, format!("{bar} Å"), egui::FontId::proportional(crate::theme::BASE_FONT), text);
    }

    fn export_surface(&mut self) {
        let Some(model) = &self.model else { return };
        let spacing = self.parameters.find("mesh_spacing").and_then(|v| v.number().ok()).unwrap_or(0.5);
        let path = if self.mesh_path.is_empty() { "shape.obj".to_string() } else { self.mesh_path.clone() };
        self.message = match model.mesh(spacing).and_then(|mesh| save_obj(&mesh, &path).map(|_| mesh.vertices.len() / 3)) {
            Ok(n) => (format!("Wrote {n} triangles to {path}"), false),
            Err(e) => (e, true),
        };
    }

    /// The step's parameters: the settings and the script (stepParameters).
    pub fn step_parameters(&self) -> Json {
        let mut p = match self.parameters.members() {
            Ok(m) if !m.is_empty() => self.parameters.clone(),
            _ => make_modifier(STEP).map(|m| m.parameters).unwrap_or_else(|_| Json::object()),
        };
        p.set("script", self.script.as_str());
        p
    }

    /// "Update Step": the step's new parameters (the dialog stays on the step, as in the C++).
    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.as_ref()?;
        let mut updated = edit.clone();
        updated.parameters = self.step_parameters();
        self.step = Some(updated.clone());
        self.message = ("Step updated".into(), false);
        Some(Outcome::UpdateStep(updated))
    }

    /// Add to Pipeline: a build-shape step with these settings and script.
    fn add_to_pipeline(&mut self) -> Option<Outcome> {
        let mut step = make_modifier(STEP).ok()?;
        let mut parameters = self.parameters.clone();
        parameters.set("script", self.script.as_str());
        step.parameters = parameters;
        self.message = ("Added to the structure pipeline".into(), false);
        Some(Outcome::AddSteps(vec![step], "Added to the structure pipeline".into()))
    }

    /// Carve (a new structure), keep/delete (an edit) or select.
    fn apply(&mut self, structure: &Structure) -> Option<Outcome> {
        let t = find_modifier_type(STEP)?;
        let mut parameters = self.parameters.clone();
        parameters.set("script", self.script.as_str());
        let mut data = PipelineData::new(structure.clone());
        data.clear_selection();
        let mode = self.mode();
        match (t.apply)(&parameters, &mut data) {
            Ok(()) => {
                let note = data.notes.first().cloned().unwrap_or_else(|| "Done".into());
                self.message = (note.clone(), false);
                Some(match mode.as_str() {
                    "carve" => Outcome::Built(data.structure, "Carved a shape".into()),
                    "select" => Outcome::Select(data.selected.iter().enumerate().filter(|(_, &v)| v != 0).map(|(i, _)| i).collect(), note),
                    _ => Outcome::Edited(data.structure, note),
                })
            }
            Err(e) => {
                self.message = (e, true);
                None
            }
        }
    }

    fn script_panel(&mut self, ui: &mut egui::Ui, crystal: &Structure) {
        crate::dialog_layout::section(ui, "Shape");
        ui.horizontal(|ui| {
            ui.menu_button("Templates", |ui| {
                for (name, script) in templates() {
                    if ui.button(*name).clicked() {
                        self.script = script.to_string();
                        ui.close();
                    }
                }
            })
            .response
            .on_hover_text("Replace the script with a ready-made shape");
            ui.menu_button("Insert", |ui| {
                let mut groups: Vec<&str> = Vec::new();
                for c in commands() {
                    if !groups.contains(&c.group) {
                        groups.push(c.group);
                    }
                }
                for group in groups {
                    ui.menu_button(group, |ui| {
                        for c in commands().iter().filter(|c| c.group == group) {
                            if ui.button(c.name).on_hover_text(format!("{}\n{}", c.help, c.snippet)).clicked() {
                                self.insert(c.snippet);
                                ui.close();
                            }
                        }
                    });
                }
            })
            .response
            .on_hover_text("Add a primitive, boolean, transform, modifier, pattern or region to the end of the script");
            crate::dialog_layout::muted(ui, "Operations act on the shape after them or on a { } block.");
        });
        ui.add(egui::TextEdit::multiline(&mut self.script).font(egui::TextStyle::Body).desired_rows(11).desired_width(f32::INFINITY).margin(egui::vec2(8.0, 6.0)));
        self.reparse(crystal);
        if let Some((text, colour)) = self.summary(crystal) {
            let kind = if colour == BAD {
                crate::dialog_layout::Note::Error
            } else if colour == WARN {
                crate::dialog_layout::Note::Warning
            } else {
                crate::dialog_layout::Note::Success
            };
            crate::dialog_layout::note(ui, kind, &text);
        }
        if let Some(m) = &self.model {
            let regions = m.region_elements();
            if !regions.is_empty() {
                crate::dialog_layout::muted(ui, format!("Regions: {}", regions.join(", ")));
            }
        }
        crate::dialog_layout::section(ui, "Settings");
        if let Some(t) = find_modifier_type(STEP) {
            let mut settings = (*t).clone();
            settings.parameters.retain(|p| p.name != "script" && p.name != "mesh_file");
            for (name, value) in parameter_form(ui, "shape settings", &settings, &self.parameters, &mut self.texts) {
                self.parameters.set(&name, value);
            }
        }
        egui::CollapsingHeader::new("Command reference").show(ui, |ui| {
            egui::Grid::new("shape commands").num_columns(2).striped(true).show(ui, |ui| {
                for c in commands() {
                    if ui.selectable_label(false, c.name).on_hover_text(c.snippet).clicked() {
                        self.insert(c.snippet);
                    }
                    ui.label(c.help);
                    ui.end_row();
                }
            });
            ui.label("Variables: R = 20; then sphere(r = R / 2). Arithmetic + - * / ^, functions sqrt, sin, cos, tan (degrees), abs, min, max, pow. Operations act on the shape after them or on a { } block. Comments: // or #.");
        });
    }

    fn side_panel(&mut self, ui: &mut egui::Ui, crystal: &Structure) {
        crate::dialog_layout::section(ui, "Section view");
        ui.horizontal(|ui| {
            ui.radio_value(&mut self.plane, 0, "x-y");
            ui.radio_value(&mut self.plane, 1, "x-z");
            ui.radio_value(&mut self.plane, 2, "y-z");
            ui.spacing_mut().slider_width = (ui.available_width() - 120.0).max(60.0);
            ui.add(egui::Slider::new(&mut self.section_at, 0.0..=1.0).text("depth").fixed_decimals(2)).on_hover_text("Where the section cuts through the shape (0: one side of its bounds, 1: the other)");
        });
        let side = ui.available_width().min(400.0);
        self.section(ui, crystal, side);
        crate::dialog_layout::section(ui, "Surface mesh");
        ui.horizontal(|ui| {
            let w = ui.available_width() - 90.0;
            crate::dialog_layout::text_field(ui, &mut self.mesh_path, w, "shape.obj");
            if ui.button("Browse...").on_hover_text("Choose where the surface mesh is written").clicked() {
                self.mesh_browser.show_for(crate::file_browser::Purpose::SaveAs, "Surface mesh", "shape.obj", &["obj"]);
            }
        });
        if ui.add_enabled(self.model.is_some(), egui::Button::new("Export Surface (OBJ)")).on_hover_text("The shape's surface as a triangle mesh, in the shape's own coordinates (Angstrom)").clicked() {
            self.export_surface();
        }
        if !self.message.0.is_empty() {
            let kind = if self.message.1 { crate::dialog_layout::Note::Error } else { crate::dialog_layout::Note::Success };
            crate::dialog_layout::note(ui, kind, &self.message.0.clone());
        }
    }
}

impl StructureDialog for ShapeDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn steps(&self) -> &'static [&'static str] {
        &[STEP]
    }

    fn title(&self) -> &'static str {
        "Shape Modeller"
    }

    fn open(&mut self, _step: &str) {
        self.step = None;
        self.open = true;
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        let mut parameters = make_modifier(STEP).map(|m| m.parameters).unwrap_or_else(|_| Json::object());
        if let Ok(members) = edit.parameters.members() {
            for (name, value) in members {
                parameters.set(name, value.clone());
            }
        }
        if let Some(script) = parameters.find("script").and_then(|v| v.string().ok()) {
            self.script = script.to_string();
        }
        self.parameters = parameters;
        self.texts.clear();
        self.message = (String::new(), false);
        self.step = Some(edit);
        self.open = true;
        true
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        let mut out = None;
        let mut open = self.open;
        // While editing a step, the shape is checked against the structure entering it.
        let step = self.step.take();
        let crystal = step.as_ref().map_or(input.source, |e| &e.input);
        let title = if step.is_some() { "Shape Modeller (pipeline step)" } else { "Shape Modeller" };
        let mut update = false;
        crate::theme::window(title).id(egui::Id::unique("shape modeller")).open(&mut open).default_size([980.0, 700.0]).collapsible(false).show(ctx, |ui| {
            if let Some(e) = &step {
                crate::structure_dialog::step_banner(ui, &e.step);
            }
            let fit = crate::dialog_layout::fit_height(ui, "shape", 52.0, 200.0);
            let height = fit.height;
            ui.horizontal_top(|ui| {
                let left = (ui.available_width() * 0.56).floor();
                let down = egui::Layout::top_down(egui::Align::LEFT);
                ui.allocate_ui_with_layout(egui::vec2(left, height), down, |ui| {
                    egui::ScrollArea::vertical().id_salt("shape script").max_height(height).show(ui, |ui| {
                        ui.set_width(left - 14.0);
                        self.script_panel(ui, crystal);
                    });
                });
                ui.add_sized([9.0, height], egui::Separator::default().vertical());
                let rest = ui.available_width();
                ui.allocate_ui_with_layout(egui::vec2(rest, height), down, |ui| {
                    egui::ScrollArea::vertical().id_salt("shape side").max_height(height).show(ui, |ui| {
                        ui.set_width(rest - 14.0);
                        self.side_panel(ui, crystal)
                    });
                });
            });
            ui.add_space(2.0);
            ui.separator();
            if step.is_some() {
                match crate::dialog_layout::action_bar(ui, "", crate::dialog_layout::Note::Info, "Update Step", Some("Close"), self.model.is_some()) {
                    1 => update = true,
                    2 => self.open = false,
                    _ => {}
                }
                fit.finish(ui);
                return;
            }
            let mode = self.mode();
            let action = match mode.as_str() {
                "carve" => "Carve",
                "select" => "Select",
                _ => "Apply",
            };
            let can = self.model.is_some() && !input.source.atoms.is_empty();
            match crate::dialog_layout::action_bar(ui, "", crate::dialog_layout::Note::Info, action, Some("Add to Pipeline"), can) {
                1 => out = self.apply(input.source),
                2 => out = self.add_to_pipeline(),
                _ => {}
            }
            fit.finish(ui);
        });
        if let Some(path) = self.mesh_browser.show(ctx) {
            self.mesh_path = path.to_string_lossy().into_owned();
        }
        self.open = open && self.open;
        if self.open {
            self.step = step;
        }
        if update {
            out = self.update_step();
        }
        out
    }

    /// script TEXT (\n for new lines) | template NAME | insert COMMAND | set NAME VALUE |
    /// plane xy|xz|yz | depth F | export PATH | apply | update-step | expect-status TEXT |
    /// expect-volume MIN MAX | expect-script TEXT.
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        match command {
            "script" => self.script = args.replace("\\n", "\n"),
            "template" => self.script = templates().iter().find(|(n, _)| n.eq_ignore_ascii_case(args.trim())).ok_or(format!("no template {args}"))?.1.to_string(),
            "insert" => {
                let c = commands().iter().find(|c| c.name == args.trim()).ok_or(format!("no command {args}"))?;
                self.insert(c.snippet);
            }
            "set" => {
                let (name, value) = args.trim().split_once(' ').unwrap_or((args.trim(), ""));
                let t = find_modifier_type(STEP).ok_or("no build-shape step")?;
                let v = self.texts.set(&t, name, value.trim())?;
                self.parameters.set(name, v);
            }
            "plane" => self.plane = ["xy", "xz", "yz"].iter().position(|p| *p == args.trim()).ok_or("plane xy|xz|yz")?,
            "depth" => self.section_at = args.trim().parse::<f32>().map_err(|e| e.to_string())?.clamp(0.0, 1.0),
            "export" => {
                self.reparse(input.source);
                self.mesh_path = args.trim().to_string();
                self.export_surface();
                if self.message.1 || self.model.is_none() {
                    return Err(if self.model.is_none() { self.parse_error.clone() } else { self.message.0.clone() });
                }
            }
            "add-to-pipeline" => return Ok(self.add_to_pipeline()),
            "update-step" => {
                let crystal = self.step.as_ref().map_or_else(|| input.source.clone(), |e| e.input.clone());
                self.reparse(&crystal);
                if self.model.is_none() {
                    return Err(self.parse_error.clone());
                }
                return self.update_step().map(Some).ok_or_else(|| "no pipeline step is being edited".to_string());
            }
            "expect-script" => {
                if self.script.trim() != args.replace("\\n", "\n").trim() {
                    return Err(format!("script is {:?}", self.script));
                }
            }
            "apply" => {
                self.reparse(input.source);
                if self.model.is_none() {
                    return Err(self.parse_error.clone());
                }
                let o = self.apply(input.source);
                return if o.is_some() { Ok(o) } else { Err(self.message.0.clone()) };
            }
            "expect-status" => {
                self.reparse(input.source);
                let text = self.summary(input.source).map(|s| s.0).unwrap_or_default();
                if !text.contains(args.trim()) && !self.message.0.contains(args.trim()) {
                    return Err(format!("status is {text:?}, message {:?}", self.message.0));
                }
            }
            "expect-volume" => {
                self.reparse(input.source);
                let w: Vec<f64> = args.split_whitespace().filter_map(|w| w.parse().ok()).collect();
                if w.len() != 2 || self.volume < w[0] || self.volume > w[1] {
                    return Err(format!("shape volume {:.1} A^3", self.volume));
                }
            }
            other => return Err(format!("the shape modeller has no command {other}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// editStep fills the script and settings from the step's parameters; Update
    /// Step hands them back (the step's defaults filled in, as stepParameters).
    #[test]
    fn step_parameters_round_trip() {
        let mut parameters = Json::object();
        parameters.set("script", "box(size = [10, 10, 10]);");
        parameters.set("mode", "delete");
        let edit = StepEdit { index: 1, step: STEP.into(), parameters, input: Structure::default(), ..Default::default() };
        let mut d = ShapeDialog::default();
        assert!(d.edit_step(edit));
        assert_eq!(d.script, "box(size = [10, 10, 10]);");
        assert_eq!(d.mode(), "delete");
        let p = d.step_parameters();
        assert_eq!(p.find("script").and_then(|v| v.string().ok()), Some("box(size = [10, 10, 10]);"));
        assert_eq!(p.find("anchor").and_then(|v| v.string().ok()), Some("auto"));
        d.script = "sphere(r = 3);".into();
        match d.update_step() {
            Some(Outcome::UpdateStep(e)) => {
                assert_eq!(e.index, 1);
                assert_eq!(e.parameters.find("script").and_then(|v| v.string().ok()), Some("sphere(r = 3);"));
                assert_eq!(e.parameters.find("mode").and_then(|v| v.string().ok()), Some("delete"));
            }
            _ => panic!("no step update"),
        }
        // The dialog stays on the step after an update (as the C++ does); the menu ends it.
        assert!(d.step.is_some());
        d.open(STEP);
        assert!(d.update_step().is_none());
    }
}
