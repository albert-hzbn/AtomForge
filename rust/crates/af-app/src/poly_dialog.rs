//! Build > Polycrystal (the C++ ui/PolyCrystalBuilderDialog): a periodic
//! reference single crystal, the box, the number of grains, the seed and the
//! grains' orientations (all random, all given or some given as Bunge Euler
//! angles); a Voronoi tessellation filled with rotated copies of the crystal.

use af_algorithms::polycrystal::{build_polycrystal, GrainOrientation, GrainOrientationMode, PolyBuildResult, PolyParams};
use af_core::structure::Structure;
use eframe::egui;

use crate::builder_preview::{cell_lengths, Orbit, Reference};
use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

const MODES: [(&str, GrainOrientationMode); 3] =
    [("All Random", GrainOrientationMode::AllRandom), ("All Specified", GrainOrientationMode::AllSpecified), ("Partial (specify some, rest random)", GrainOrientationMode::PartialSpecified)];
const BOX_HINT: &str = "Supply a periodic reference and the three target box lengths in Å. The dimensions define the total specimen, not the size of every grain.";

pub struct PolyDialog {
    open: bool,
    params: PolyParams,
    last: PolyBuildResult,
    new_grain: i32,
    reference: Reference,
    orbit: Orbit,
    /// The pipeline step being edited (Edit in Dialog).
    step: Option<StepEdit>,
}

impl Default for PolyDialog {
    fn default() -> Self {
        Self { open: false, params: PolyParams::default(), last: PolyBuildResult::default(), new_grain: 1, reference: Reference::new("Load reference (Polycrystal)"), orbit: Orbit::default(), step: None }
    }
}

fn euler(ui: &mut egui::Ui, o: &mut GrainOrientation) {
    for (v, l) in [(&mut o.phi1, "phi1"), (&mut o.phi, "Phi"), (&mut o.phi2, "phi2")] {
        ui.add(egui::DragValue::new(v).speed(0.5).fixed_decimals(1));
        ui.label(l);
    }
}

/// Three or four whitespace-separated numbers (istream >> float ...).
fn floats(text: &str, n: usize) -> Option<Vec<f32>> {
    let v: Vec<f32> = text.split_whitespace().take(n).map(|w| w.parse::<f32>().ok()).collect::<Option<_>>()?;
    (v.len() == n).then_some(v)
}

impl PolyDialog {
    /// Fills the dialog from a build-poly step (editStep); missing flags take
    /// runPoly's defaults, the reference is --input or the step's input.
    pub fn load_step(&mut self, edit: &StepEdit) {
        let r = edit.reader();
        let mut p = PolyParams {
            size_x: r.number("--sizex", 50.0) as f32,
            size_y: r.number("--sizey", 50.0) as f32,
            size_z: r.number("--sizez", 50.0) as f32,
            num_grains: r.integer("--grains", 8),
            seed: r.integer("--seed", 42),
            ..PolyParams::default()
        };
        let eulers = r.values("--euler");
        for (i, text) in eulers.iter().enumerate() {
            if let Some(v) = floats(text, 3) {
                p.specified_orientations.push(GrainOrientation { grain_index: i as i32, phi1: v[0], phi: v[1], phi2: v[2] });
            }
        }
        let grain_eulers = r.values("--grain-euler");
        for text in &grain_eulers {
            let mut w = text.split_whitespace();
            let grain = w.next().and_then(|g| g.parse::<i32>().ok());
            let rest: Vec<&str> = w.collect();
            if let (Some(grain), Some(v)) = (grain, floats(&rest.join(" "), 3)) {
                if grain >= 1 {
                    p.specified_orientations.push(GrainOrientation { grain_index: grain - 1, phi1: v[0], phi: v[1], phi2: v[2] });
                }
            }
        }
        p.orientation_mode = if grain_eulers.is_empty() && !eulers.is_empty() && eulers.len() as i32 >= p.num_grains {
            GrainOrientationMode::AllSpecified
        } else if !p.specified_orientations.is_empty() {
            GrainOrientationMode::PartialSpecified
        } else {
            GrainOrientationMode::AllRandom
        };
        self.params = p;
        self.last = PolyBuildResult::default();
        self.reference.from_step_edit(edit, true, "(pipeline step input)");
        self.orbit = Orbit::default();
    }

    /// The options of the matching pipeline step (stepOptions).
    pub fn step_options(&self) -> String {
        let mut w = af_pipeline::options::Writer::new();
        if !self.reference.from_step && !self.reference.path.is_empty() {
            w.add("--input", &self.reference.path);
        }
        let p = &self.params;
        w.add_number("--sizex", p.size_x as f64).add_number("--sizey", p.size_y as f64).add_number("--sizez", p.size_z as f64);
        w.add_integer("--grains", p.num_grains as i64).add_integer("--seed", p.seed as i64);
        match p.orientation_mode {
            GrainOrientationMode::AllSpecified => {
                // One --euler per grain, in grain order.
                let n = (p.num_grains.max(0) as usize).min(p.specified_orientations.len());
                for o in &p.specified_orientations[..n] {
                    w.add_numbers("--euler", &[o.phi1 as f64, o.phi as f64, o.phi2 as f64]);
                }
            }
            GrainOrientationMode::PartialSpecified => {
                // Grain numbers are explicit: the specified grains need not be the first ones.
                for o in &p.specified_orientations {
                    w.add_numbers("--grain-euler", &[(o.grain_index + 1) as f64, o.phi1 as f64, o.phi as f64, o.phi2 as f64]);
                }
            }
            GrainOrientationMode::AllRandom => {}
        }
        w.str().to_string()
    }

    /// Update Step: the step's new options (and use_input when the reference
    /// is the step's input); ends the step edit.
    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        self.open = false;
        self.last = PolyBuildResult::default();
        let mut updated = edit.with_options(&self.step_options());
        if self.reference.from_step {
            updated = updated.with_parameter("use_input", true);
        }
        Some(Outcome::UpdateStep(updated))
    }

    fn clamp(&mut self) {
        let p = &mut self.params;
        p.size_x = p.size_x.max(1.0);
        p.size_y = p.size_y.max(1.0);
        p.size_z = p.size_z.max(1.0);
        p.num_grains = p.num_grains.clamp(1, 500);
        if p.orientation_mode == GrainOrientationMode::AllSpecified && p.specified_orientations.len() != p.num_grains as usize {
            p.specified_orientations.resize(p.num_grains as usize, GrainOrientation::default());
            for (i, o) in p.specified_orientations.iter_mut().enumerate() {
                o.grain_index = i as i32;
            }
        }
    }

    fn can_build(&self) -> bool {
        !self.reference.structure.atoms.is_empty() && self.reference.structure.has_unit_cell
    }

    fn build(&mut self, colours: &[[f32; 3]]) -> Option<Outcome> {
        self.clamp();
        let mut s = Structure::default();
        self.last = build_polycrystal(&mut s, &self.reference.structure, &self.params, colours);
        self.last.success.then(|| Outcome::Built(s, "Built a polycrystal".into()))
    }

    fn reference_panel(&mut self, ui: &mut egui::Ui, input: &Input) {
        ui.strong("Reference Single Crystal");
        ui.separator();
        let r = &self.reference.structure;
        if r.atoms.is_empty() {
            ui.weak("(none loaded)");
        } else {
            ui.label(format!("File: {}", self.reference.name));
            ui.label(format!("{} atoms", r.atoms.len()));
            if r.has_unit_cell {
                let [a, b, c] = cell_lengths(r);
                ui.weak(format!("a={a:.3}  b={b:.3}  c={c:.3} Å"));
            }
        }
        let height = (ui.available_height() - 60.0).max(120.0);
        crate::builder_preview::structure(ui, egui::vec2(ui.available_width(), height), &self.reference.structure, input.colours, &mut self.orbit, &["Load a structure file", "or use the structure in the view"]);
        if !self.reference.structure.atoms.is_empty() {
            ui.weak("Left-drag = orbit   Scroll = zoom");
        }
        if !self.reference.status.is_empty() {
            ui.label(&self.reference.status);
        }
    }

    fn options(&mut self, ui: &mut egui::Ui) {
        ui.strong("Builder Options");
        ui.separator();
        let p = &mut self.params;
        ui.label("Simulation Box Size (Å)");
        for (v, l) in [(&mut p.size_x, "X"), (&mut p.size_y, "Y"), (&mut p.size_z, "Z")] {
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(v).speed(0.5).range(1.0..=f32::MAX).fixed_decimals(2));
                ui.label(l);
            })
            .response
            .on_hover_text(BOX_HINT);
        }
        ui.separator();
        ui.label("Number of Grains");
        ui.add(egui::DragValue::new(&mut p.num_grains).range(1..=500));
        ui.separator();
        ui.label("Random Seed");
        ui.add(egui::DragValue::new(&mut p.seed));
        ui.separator();
        ui.label("Grain Orientations");
        let label = MODES.iter().find(|m| m.1 == p.orientation_mode).map(|m| m.0).unwrap_or("");
        egui::ComboBox::from_id_salt("poly orientation mode").width(ui.available_width()).selected_text(label).show_ui(ui, |ui| {
            for (l, m) in MODES {
                ui.selectable_value(&mut p.orientation_mode, m, l);
            }
        });
        match p.orientation_mode {
            GrainOrientationMode::AllSpecified => {
                if p.specified_orientations.len() != p.num_grains as usize {
                    p.specified_orientations.resize(p.num_grains as usize, GrainOrientation::default());
                    for (i, o) in p.specified_orientations.iter_mut().enumerate() {
                        o.grain_index = i as i32;
                    }
                }
                ui.weak("Bunge Euler angles (phi1, Phi, phi2) in degrees");
                egui::ScrollArea::vertical().id_salt("poly all orientations").max_height(180.0).show(ui, |ui| {
                    for (i, o) in p.specified_orientations.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(format!("Grain {}:", i + 1));
                            euler(ui, o);
                        });
                    }
                });
            }
            GrainOrientationMode::PartialSpecified => {
                ui.weak("Specify orientations for selected grains. Unspecified grains get random orientations.");
                ui.weak("Bunge Euler angles (phi1, Phi, phi2) in degrees");
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut self.new_grain).range(1..=p.num_grains.max(1)));
                    ui.label("Grain #");
                    self.new_grain = self.new_grain.clamp(1, p.num_grains.max(1));
                    if ui.button("Add Grain Orientation").clicked() && !p.specified_orientations.iter().any(|o| o.grain_index == self.new_grain - 1) {
                        p.specified_orientations.push(GrainOrientation { grain_index: self.new_grain - 1, ..Default::default() });
                    }
                });
                let mut remove = None;
                egui::ScrollArea::vertical().id_salt("poly partial orientations").max_height(180.0).show(ui, |ui| {
                    for (i, o) in p.specified_orientations.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(format!("Grain {}:", o.grain_index + 1));
                            euler(ui, o);
                            if ui.button("X").clicked() {
                                remove = Some(i);
                            }
                        });
                    }
                });
                if let Some(i) = remove {
                    p.specified_orientations.remove(i);
                }
            }
            GrainOrientationMode::AllRandom => {
                ui.weak("All grains will receive uniformly random orientations.");
            }
        }
        ui.separator();
        ui.label("Result");
        if self.last.success {
            ui.colored_label(egui::Color32::from_rgb(70, 170, 90), &self.last.message);
        } else if !self.last.message.is_empty() {
            ui.colored_label(egui::Color32::from_rgb(220, 80, 60), &self.last.message);
        } else {
            ui.weak("(not built yet)");
        }
    }
}

impl StructureDialog for PolyDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        let _ = self.reference.load(path, true);
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-poly"]
    }

    fn title(&self) -> &'static str {
        "Build Polycrystal"
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
        self.step = None;
        self.last = PolyBuildResult::default();
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
        crate::theme::modal_window("Build Polycrystal").open(&mut open).default_size([820.0, 640.0]).show(ctx, |ui| {
            if let Some(edit) = &self.step {
                crate::structure_dialog::step_banner(ui, &edit.step);
            }
            self.reference.buttons(ui, input.source, true);
            ui.add_space(4.0);
            let height = 560.0f32.min(ctx.content_rect().height() - 200.0).max(300.0);
            let stacked = ui.available_width() < 640.0;
            let left = |me: &mut Self, ui: &mut egui::Ui| {
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    ui.set_width(if stacked { ui.available_width() } else { 340.0 });
                    ui.set_height(height);
                    me.reference_panel(ui, input);
                });
            };
            let right = |me: &mut Self, ui: &mut egui::Ui| {
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.set_height(height);
                    egui::ScrollArea::vertical().id_salt("poly options").show(ui, |ui| me.options(ui));
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
            let label = if editing { "Update Step" } else { "Build" };
            match crate::dialog_layout::footer(ui, label, Some("Close"), editing || self.can_build()) {
                1 => built = if editing { self.update_step() } else { self.build(input.colours) },
                2 => self.open = false,
                _ => {}
            }
        });
        self.open = open && self.open;
        if !self.open {
            self.step = None;
        }
        self.reference.browser(ctx, true);
        built
    }

    /// reference PATH | reference view | size X Y Z | grains N | seed N |
    /// mode random|specified|partial | orientation GRAIN PHI1 PHI PHI2 |
    /// build | expect-grains N
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let num = |i: usize| words.get(i).and_then(|w| w.parse::<f64>().ok()).ok_or(format!("{command}: expected a number in {args:?}"));
        match command {
            "reference" if args.trim() == "view" => self.reference.take_view(input.source, true)?,
            "reference" => self.reference.load(args.trim(), true)?,
            "size" => (self.params.size_x, self.params.size_y, self.params.size_z) = (num(0)? as f32, num(1)? as f32, num(2)? as f32),
            "grains" => self.params.num_grains = num(0)? as i32,
            "seed" => self.params.seed = num(0)? as i32,
            "mode" => {
                self.params.orientation_mode = match args.trim() {
                    "random" => GrainOrientationMode::AllRandom,
                    "specified" => GrainOrientationMode::AllSpecified,
                    "partial" => GrainOrientationMode::PartialSpecified,
                    other => return Err(format!("no orientation mode {other}")),
                };
                self.clamp();
            }
            "orientation" => {
                // 1-based grain, Bunge angles in degrees.
                let grain = num(0)? as i32 - 1;
                let o = GrainOrientation { grain_index: grain, phi1: num(1)? as f32, phi: num(2)? as f32, phi2: num(3)? as f32 };
                let list = &mut self.params.specified_orientations;
                match list.iter_mut().find(|x| x.grain_index == grain) {
                    Some(x) => *x = o,
                    None => list.push(o),
                }
            }
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
                    return Err("load a periodic reference first".into());
                }
                let out = self.build(input.colours);
                return if out.is_some() { Ok(out) } else { Err(self.last.message.clone()) };
            }
            "expect-grains" => {
                if self.last.num_grains != num(0)? as i32 {
                    return Err(format!("{} grains: {}", self.last.num_grains, self.last.message));
                }
            }
            other => return Err(format!("the polycrystal dialog has no {other}")),
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
        StepEdit { index: 0, step: "build-poly".into(), parameters, input: Default::default(), ..Default::default() }
    }

    #[test]
    fn step_options_round_trip() {
        for options in [
            "--sizex 30 --sizey 30 --sizez 30 --grains 4 --seed 7",
            "--sizex 40 --sizey 35.5 --sizez 30 --grains 2 --seed 3 --euler \"0 0 0\" --euler \"45 30 10\"",
            "--sizex 50 --sizey 50 --sizez 50 --grains 8 --seed 42 --grain-euler \"3 10 20 30\" --grain-euler \"1 0 45 0\"",
        ] {
            let mut d = PolyDialog::default();
            assert!(d.edit_step(edit(options)));
            assert_eq!(d.step_options(), options);
        }
        let mut d = PolyDialog::default();
        d.edit_step(edit("--euler \"1 2 3\""));
        assert_eq!(d.params.orientation_mode, GrainOrientationMode::PartialSpecified);
        assert_eq!(d.step_options(), "--sizex 50 --sizey 50 --sizez 50 --grains 8 --seed 42 --grain-euler \"1 1 2 3\"");
    }
}
