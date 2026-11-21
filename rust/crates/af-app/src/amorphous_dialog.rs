//! Build > Amorphous Structure (the C++ ui/AmorphousBuilderDialog): the
//! composition, a box from a target density or given lengths, the cell scale
//! factor, and advanced placement settings (covalent tolerance, seed,
//! attempts, element-pair minimum distances); random sequential addition.

use af_algorithms::amorphous::{build_amorphous_structure, AmorphousBoxMode, AmorphousElementSpec, AmorphousParams, AmorphousPairDist, AmorphousResult};
use af_core::elements;
use eframe::egui;

use crate::display_settings::periodic_table;
use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

#[derive(Clone, Copy)]
struct PairRow {
    z1: i32,
    z2: i32,
    min_dist: f32,
    enabled: bool,
}

pub struct AmorphousDialog {
    open: bool,
    elements: Vec<AmorphousElementSpec>,
    /// 0 manual, 1 from density (the C++ m_boxMode).
    box_mode: i32,
    box_a: f32,
    box_b: f32,
    box_c: f32,
    target_density: f32,
    cell_scale_factor: f32,
    cov_tolerance: f32,
    seed: i32,
    max_attempts: i32,
    periodic: bool,
    pair_rows: Vec<PairRow>,
    picking: Option<usize>,
    pick_element: i32,
    status: (String, bool),
    last: Option<AmorphousResult>,
    radii: Vec<f32>,
    /// The pipeline step being edited (Edit in Dialog).
    step: Option<StepEdit>,
    /// The pair rows are rebuilt when the window next draws (as the C++ on open).
    rows_pending: bool,
}

impl Default for AmorphousDialog {
    fn default() -> Self {
        Self {
            open: false,
            elements: vec![AmorphousElementSpec { atomic_number: 14, count: 100 }],
            box_mode: 1,
            box_a: 20.0,
            box_b: 20.0,
            box_c: 20.0,
            target_density: 2.0,
            cell_scale_factor: 1.0,
            cov_tolerance: 0.75,
            seed: 42,
            max_attempts: 1000,
            periodic: true,
            pair_rows: Vec::new(),
            picking: None,
            pick_element: 14,
            status: (String::new(), false),
            last: None,
            radii: elements::literature_covalent_radii(),
            step: None,
            rows_pending: false,
        }
    }
}

const SUCCESS: egui::Color32 = egui::Color32::from_rgb(70, 170, 90);
const ERROR: egui::Color32 = egui::Color32::from_rgb(220, 80, 60);

impl AmorphousDialog {
    fn total(&self) -> i32 {
        self.elements.iter().map(|e| e.count.max(0)).sum()
    }

    fn default_min_dist(&self, z1: i32, z2: i32) -> f32 {
        let r = |z: i32| if z >= 1 && (z as usize) < self.radii.len() { self.radii[z as usize] } else { 1.5 };
        (r(z1) + r(z2)) * self.cov_tolerance
    }

    /// Every element pair of the composition, keeping existing overrides.
    fn rebuild_pair_rows(&mut self) {
        let mut zs: Vec<i32> = Vec::new();
        for e in &self.elements {
            if e.count > 0 && !zs.contains(&e.atomic_number) {
                zs.push(e.atomic_number);
            }
        }
        zs.sort();
        let mut rows = Vec::new();
        for i in 0..zs.len() {
            for j in i..zs.len() {
                let (z1, z2) = (zs[i], zs[j]);
                match self.pair_rows.iter().find(|p| p.z1.min(p.z2) == z1 && p.z1.max(p.z2) == z2) {
                    Some(p) => rows.push(*p),
                    None => rows.push(PairRow { z1, z2, min_dist: self.default_min_dist(z1, z2), enabled: false }),
                }
            }
        }
        self.pair_rows = rows;
    }

    /// The box the builder will use (summary line and density hint).
    fn box_lengths(&self) -> [f32; 3] {
        let total = self.total();
        if self.box_mode == 1 && total > 0 {
            let mass: f64 = self.elements.iter().filter(|e| (1..=118).contains(&e.atomic_number)).map(|e| e.count as f64 * elements::atomic_mass(e.atomic_number)).sum();
            let volume = mass * 1.66054 / (self.target_density as f64).max(1e-3);
            [volume.cbrt() as f32 * self.cell_scale_factor; 3]
        } else {
            [self.box_a * self.cell_scale_factor, self.box_b * self.cell_scale_factor, self.box_c * self.cell_scale_factor]
        }
    }

    fn params(&self) -> AmorphousParams {
        AmorphousParams {
            elements: self.elements.iter().filter(|e| e.count > 0).copied().collect(),
            box_mode: if self.box_mode == 0 { AmorphousBoxMode::Manual } else { AmorphousBoxMode::AutoFromDensity },
            box_a: self.box_a,
            box_b: self.box_b,
            box_c: self.box_c,
            target_density: self.target_density,
            cell_scale_factor: self.cell_scale_factor,
            pair_distances: self.pair_rows.iter().filter(|p| p.enabled).map(|p| AmorphousPairDist { z1: p.z1, z2: p.z2, min_dist: p.min_dist }).collect(),
            covalent_tolerance: self.cov_tolerance,
            seed: self.seed.max(0) as u32,
            max_attempts: self.max_attempts.max(1),
            periodic: self.periodic,
        }
    }

    /// Build: the structure when atoms were placed, with the status line.
    fn build(&mut self, colours: &[[f32; 3]]) -> Option<Outcome> {
        let r = build_amorphous_structure(&self.params(), &self.radii, colours);
        let mut message = r.message.clone();
        let mut out = None;
        if r.placed_atoms > 0 {
            message += &format!("  |  density: {:.3} g/cm³  |  box: {:.3} × {:.3} × {:.3} Å", r.actual_density, r.box_a, r.box_b, r.box_c);
            out = Some(Outcome::Built(r.output.clone(), "Built an amorphous structure".into()));
        }
        self.status = (message, !r.success);
        self.last = Some(r);
        out
    }

    fn element_list(text: &str) -> Result<Vec<AmorphousElementSpec>, String> {
        // "Si 64; O 128"
        let mut out = Vec::new();
        for part in text.split(';').filter(|p| !p.trim().is_empty()) {
            let w: Vec<&str> = part.split_whitespace().collect();
            let z = w.first().map(|s| elements::atomic_number(s)).unwrap_or(0);
            let count = w.get(1).and_then(|c| c.parse::<i32>().ok());
            match (z > 0, count, w.len()) {
                (true, Some(count), 2) => out.push(AmorphousElementSpec { atomic_number: z, count }),
                _ => return Err(format!("an element is SYMBOL COUNT, not {part:?}")),
            }
        }
        Ok(out)
    }

    /// Fills the dialog from a build-amorphous step's options (editStep);
    /// missing flags take runAmorphous's defaults.
    pub fn load_step(&mut self, edit: &StepEdit) {
        let r = edit.reader();
        self.elements.clear();
        for text in r.values("--element") {
            let w: Vec<&str> = text.split_whitespace().collect();
            let (Some(symbol), Some(count)) = (w.first(), w.get(1).and_then(|c| c.parse::<i32>().ok())) else { continue };
            let z = elements::atomic_number(symbol);
            if z > 0 {
                self.elements.push(AmorphousElementSpec { atomic_number: z, count });
            }
        }
        if self.elements.is_empty() {
            self.elements.push(AmorphousElementSpec { atomic_number: 14, count: 100 });
        }
        let manual = r.has("--boxa") || r.has("--boxb") || r.has("--boxc");
        self.box_mode = if manual { 0 } else { 1 };
        self.box_a = r.number("--boxa", 20.0) as f32;
        self.box_b = r.number("--boxb", 20.0) as f32;
        self.box_c = r.number("--boxc", 20.0) as f32;
        self.target_density = r.number("--density", 2.0) as f32;
        self.cell_scale_factor = r.number("--scale", 1.0) as f32;
        self.cov_tolerance = r.number("--covtol", 0.75) as f32;
        self.seed = r.integer("--seed", 42);
        self.max_attempts = r.integer("--attempts", 1000);
        self.periodic = !r.has("--no-periodic");
        self.pair_rows.clear();
        for text in r.values("--mindist") {
            let w: Vec<&str> = text.split_whitespace().collect();
            let (Some(s1), Some(s2), Some(d)) = (w.first(), w.get(1), w.get(2).and_then(|d| d.parse::<f32>().ok())) else { continue };
            let (z1, z2) = (elements::atomic_number(s1), elements::atomic_number(s2));
            if z1 <= 0 || z2 <= 0 {
                continue;
            }
            self.pair_rows.push(PairRow { z1: z1.min(z2), z2: z1.max(z2), min_dist: d, enabled: true });
        }
    }

    /// The options of the matching pipeline step (stepOptions).
    pub fn step_options(&self) -> String {
        let mut w = af_pipeline::options::Writer::new();
        for e in self.elements.iter().filter(|e| e.count > 0) {
            w.add("--element", &format!("{} {}", elements::symbol(e.atomic_number), e.count));
        }
        if self.box_mode == 0 {
            w.add_number("--boxa", self.box_a as f64).add_number("--boxb", self.box_b as f64).add_number("--boxc", self.box_c as f64);
        } else {
            w.add_number("--density", self.target_density as f64);
        }
        w.add_number("--scale", self.cell_scale_factor as f64);
        w.add_number("--covtol", self.cov_tolerance as f64);
        w.add_integer("--seed", self.seed.max(0) as i64);
        w.add_integer("--attempts", self.max_attempts.max(1) as i64);
        for p in self.pair_rows.iter().filter(|p| p.enabled) {
            w.add("--mindist", &format!("{} {} {}", elements::symbol(p.z1), elements::symbol(p.z2), af_pipeline::options::Writer::number(p.min_dist as f64)));
        }
        if !self.periodic {
            w.flag("--no-periodic");
        }
        w.str().to_string()
    }

    /// Update Step: the step with this dialog's options; ends the step edit.
    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        if std::mem::take(&mut self.rows_pending) {
            self.rebuild_pair_rows();
        }
        self.open = false;
        Some(Outcome::UpdateStep(edit.with_options(&self.step_options())))
    }

    fn controls(&mut self, ui: &mut egui::Ui, colours: &[[f32; 3]], built: &mut Option<Outcome>) {
        if let Some(edit) = &self.step {
            crate::structure_dialog::step_banner(ui, &edit.step);
        }
        ui.strong("Composition");
        let total = self.total();
        let mut remove = None;
        egui::Grid::new("amorphous composition").num_columns(5).spacing([12.0, 4.0]).show(ui, |ui| {
            for h in ["Element", "Name", "Count", "at%", ""] {
                ui.weak(h);
            }
            ui.end_row();
            for (i, e) in self.elements.iter_mut().enumerate() {
                if ui.add(egui::Button::new(elements::symbol(e.atomic_number)).min_size(egui::vec2(48.0, 0.0))).clicked() {
                    self.picking = Some(i);
                    self.pick_element = e.atomic_number;
                }
                if ui.add(egui::Button::new(elements::name(e.atomic_number)).frame(false)).on_hover_text("Click to change element").clicked() {
                    self.picking = Some(i);
                    self.pick_element = e.atomic_number;
                }
                ui.add(egui::DragValue::new(&mut e.count).range(0..=100_000).suffix(" atoms"));
                let pct = if total > 0 { 100.0 * e.count as f32 / total as f32 } else { 0.0 };
                ui.label(format!("{pct:.1}"));
                if ui.button("X").clicked() {
                    remove = Some(i);
                }
                ui.end_row();
            }
        });
        if let Some(i) = remove.filter(|_| self.elements.len() > 1) {
            self.elements.remove(i);
            self.rebuild_pair_rows();
        }
        ui.horizontal(|ui| {
            if ui.button("+ Add element").clicked() {
                self.elements.push(AmorphousElementSpec { atomic_number: 14, count: 10 });
                self.rebuild_pair_rows();
            }
            ui.weak(format!("Total: {total} atoms"));
        });

        ui.add_space(6.0);
        ui.strong("Simulation Box");
        ui.horizontal(|ui| {
            ui.radio_value(&mut self.box_mode, 1, "Auto from target density")
                .on_hover_text("Enter mass density in g cm^-3. The builder estimates a cubic box using the total atomic mass. The displayed side length includes the Cell scale factor.");
            ui.radio_value(&mut self.box_mode, 0, "Manual dimensions").on_hover_text("Enter a,b,c in Å to set box lengths explicitly. This mode replaces density-derived sizing.");
        });
        if self.box_mode == 1 {
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.target_density).speed(0.01).range(0.01..=30.0).fixed_decimals(3));
                ui.label("Target density (g/cm³)");
            });
            if total > 0 {
                ui.weak(format!("  → cubic box {:.2} Å", self.box_lengths()[0]));
            }
        } else {
            ui.horizontal(|ui| {
                for (v, l) in [(&mut self.box_a, "a (Å)"), (&mut self.box_b, "b (Å)"), (&mut self.box_c, "c (Å)")] {
                    ui.add(egui::DragValue::new(v).speed(0.5).range(1.0..=5000.0).fixed_decimals(2));
                    ui.label(l);
                }
            });
        }
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.cell_scale_factor).speed(0.01).range(0.5..=3.0).fixed_decimals(3));
            ui.label("Cell scale factor");
        })
        .response
        .on_hover_text("Multiplies box lengths before packing. Scaling lengths by s multiplies volume by s^3, so a density-derived box has effective density /s^3. For example, s=1.1 yields approximately 0.751.\nValues > 1 expand the box before packing (lower effective density). This mirrors the 1.1x cell-expansion trick used in MD melt-quench workflows to help the crystal disorder before quenching.");
        ui.checkbox(&mut self.periodic, "Periodic boundary conditions (attach unit cell)")
            .on_hover_text("Controls whether the output carries the unit cell. The current packing implementation checks minimum-image separations across the box even when this output-cell option is off.");

        ui.add_space(6.0);
        egui::CollapsingHeader::new("Advanced Settings").id_salt("amorphous advanced").show(ui, |ui| {
            ui.weak("Placement");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.cov_tolerance).speed(0.01).range(0.3..=2.0).fixed_decimals(2));
                ui.label("Covalent radius tolerance");
            })
            .response
            .on_hover_text("The fallback minimum separation is (r_i^cov+r_j^cov)t. Smaller t permits closer contacts; it does not calculate bonding.\nDefault min separation = (r_cov1 + r_cov2) * tolerance. 0.75 prevents overlaps without enforcing full bonding. Decrease for denser random packings.");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.seed).range(0..=i32::MAX));
                ui.label("RNG seed (0 = time)");
            })
            .response
            .on_hover_text("Use a nonzero integer to repeat the random placement. Zero requests a time-dependent seed.");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.max_attempts).speed(10.0).range(1..=100_000));
                ui.label("Max attempts per atom");
            })
            .response
            .on_hover_text("Bounds trials for each requested atom. Raising it increases work but cannot guarantee a feasible packing.\nIncrease if atoms are skipped; decrease to limit build time.");
            ui.weak("Element-pair minimum distances");
            ui.weak("Check a row to override the covalent-radii default for that pair.");
            if self.pair_rows.is_empty() {
                ui.weak("(no pairs – add elements above)");
            } else {
                let defaults: Vec<f32> = self.pair_rows.iter().map(|p| self.default_min_dist(p.z1, p.z2)).collect();
                egui::Grid::new("amorphous pairs").num_columns(3).spacing([16.0, 4.0]).show(ui, |ui| {
                    ui.weak("  Pair");
                    ui.weak("Default (Å)");
                    ui.weak("Override (Å)");
                    ui.end_row();
                    for (p, d) in self.pair_rows.iter_mut().zip(defaults) {
                        ui.checkbox(&mut p.enabled, format!("{}–{}", elements::symbol(p.z1), elements::symbol(p.z2)));
                        ui.weak(format!("{d:.3}"));
                        if p.enabled {
                            ui.add(egui::DragValue::new(&mut p.min_dist).speed(0.01).range(0.1..=10.0).fixed_decimals(3).suffix(" Å"));
                        } else {
                            ui.weak("(using default)");
                        }
                        ui.end_row();
                    }
                });
            }
        });

        ui.separator();
        let [a, b, c] = self.box_lengths();
        ui.weak(format!("Box: {a:.2} × {b:.2} × {c:.2} Å   |   Volume: {:.1} Å³   |   Atoms: {total}", a * b * c));
        if !self.status.0.is_empty() {
            ui.colored_label(if self.status.1 { ERROR } else { SUCCESS }, &self.status.0);
        }
        let label = if self.step.is_some() { "Update Step" } else { "Build" };
        match crate::dialog_layout::footer(ui, label, Some("Close"), total > 0) {
            1 => *built = if self.step.is_some() { self.update_step() } else { self.build(colours) },
            2 => {
                self.open = false;
                self.step = None;
            }
            _ => {}
        }
    }
}

impl StructureDialog for AmorphousDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-amorphous"]
    }

    fn title(&self) -> &'static str {
        "Amorphous Structure Builder"
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
        self.step = None;
        self.rows_pending = true;
        self.status = (String::new(), false);
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        if std::mem::take(&mut self.rows_pending) {
            self.rebuild_pair_rows();
        }
        let mut built = None;
        let mut open = self.open;
        crate::theme::modal_window("Amorphous Structure Builder").open(&mut open).default_size([680.0, 600.0]).show(ctx, |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).max_height(ctx.content_rect().height() - 120.0).show(ui, |ui| self.controls(ui, input.colours, &mut built));
        });
        self.open = open && self.open;
        if !self.open {
            self.step = None;
        }
        if let Some(i) = self.picking {
            let mut keep = true;
            let mut chosen = false;
            crate::theme::modal_window("Element").id(egui::Id::unique("amorphous element")).collapsible(false).open(&mut keep).show(ctx, |ui| {
                ui.weak("Click an element.");
                chosen = periodic_table(ui, &mut self.pick_element, Some(input.colours));
            });
            if chosen {
                if let Some(e) = self.elements.get_mut(i) {
                    e.atomic_number = self.pick_element;
                }
                self.rebuild_pair_rows();
            }
            if chosen || !keep {
                self.picking = None;
            }
        }
        built
    }

    /// elements "Si 64; O 128" | density D | box A B C | scale S | covtol T |
    /// seed N | attempts N | periodic on|off | mindist S1 S2 D | build |
    /// expect-placed N | expect-density MIN MAX
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let num = |i: usize| words.get(i).and_then(|w| w.parse::<f64>().ok()).ok_or(format!("{command}: expected a number in {args:?}"));
        match command {
            "elements" => {
                self.elements = Self::element_list(args)?;
                self.rebuild_pair_rows();
            }
            "density" => {
                self.box_mode = 1;
                self.target_density = num(0)? as f32;
            }
            "box" => {
                self.box_mode = 0;
                (self.box_a, self.box_b, self.box_c) = (num(0)? as f32, num(1)? as f32, num(2)? as f32);
            }
            "scale" => self.cell_scale_factor = num(0)? as f32,
            "covtol" => self.cov_tolerance = num(0)? as f32,
            "seed" => self.seed = num(0)? as i32,
            "attempts" => self.max_attempts = num(0)? as i32,
            "periodic" => self.periodic = args.trim() == "on",
            "mindist" => {
                self.rebuild_pair_rows();
                let (z1, z2) = (elements::atomic_number(words.first().unwrap_or(&"")), elements::atomic_number(words.get(1).unwrap_or(&"")));
                let d = num(2)? as f32;
                let row = self.pair_rows.iter_mut().find(|p| p.z1 == z1.min(z2) && p.z2 == z1.max(z2)).ok_or(format!("no {args} pair in the composition"))?;
                row.enabled = true;
                row.min_dist = d;
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
                let out = self.build(input.colours);
                if self.status.1 && out.is_none() {
                    return Err(self.status.0.clone());
                }
                return Ok(out);
            }
            "expect-placed" => {
                let r = self.last.as_ref().ok_or("nothing built")?;
                if r.placed_atoms != num(0)? as i32 {
                    return Err(format!("placed {} of {} atoms: {}", r.placed_atoms, r.requested_atoms, r.message));
                }
            }
            "expect-density" => {
                let r = self.last.as_ref().ok_or("nothing built")?;
                let d = r.actual_density as f64;
                if d < num(0)? || d > num(1)? {
                    return Err(format!("density {d}"));
                }
            }
            other => return Err(format!("the amorphous dialog has no {other}")),
        }
        Ok(None)
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        self.load_step(&edit);
        self.step = Some(edit);
        self.open = true;
        self.rows_pending = true;
        self.status = (String::new(), false);
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edit(options: &str) -> StepEdit {
        let mut parameters = af_core::json::Json::object();
        parameters.set("options", options);
        StepEdit { index: 0, step: "build-amorphous".into(), parameters, input: Default::default(), ..Default::default() }
    }

    #[test]
    fn step_options_round_trip() {
        for options in [
            "--element \"Si 64\" --element \"O 128\" --density 2.2 --scale 1 --covtol 0.75 --seed 42 --attempts 1000",
            "--element \"Cu 50\" --boxa 12 --boxb 13.5 --boxc 14 --scale 1.1 --covtol 0.8 --seed 7 --attempts 500 --mindist \"Cu Cu 2.1\" --no-periodic",
        ] {
            let mut d = AmorphousDialog::default();
            assert!(d.edit_step(edit(options)));
            let a = d.step_options();
            let mut again = AmorphousDialog::default();
            again.edit_step(edit(&a));
            assert_eq!(again.step_options(), a, "{options}");
            // Floats print as the C++ float -> double (%.10g).
            assert!(a.contains("--element \"Si 64\"") || a.contains("--element \"Cu 50\""), "{a}");
        }
        let mut d = AmorphousDialog::default();
        d.edit_step(edit("--element \"Si 64\" --density 2.2"));
        assert_eq!(d.step_options(), "--element \"Si 64\" --density 2.200000048 --scale 1 --covtol 0.75 --seed 42 --attempts 1000");
    }
}
