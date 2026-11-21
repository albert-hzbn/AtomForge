//! Build > Substitutional Solid Solution (the C++
//! ui/SubstitutionalSolidSolutionDialog): a host structure (the scene or a
//! file) keeps its positions and cell while its sites get random elements in
//! a target composition. Percentages are edited per element, the last one
//! taking up the rest; counts come by largest remainder.

use af_algorithms::solid_solution::{build_substitutional_solid_solution, ElementFraction, SssParams};
use af_core::elements;
use af_core::structure::Structure;
use eframe::egui;

use crate::display_settings::periodic_table;
use crate::file_browser::{FileBrowser, Purpose};
use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

#[derive(Clone, Copy, Debug, PartialEq)]
struct Entry {
    atomic_number: i32,
    /// Atomic percent (any scale; normalised).
    percent: f32,
    /// From the percents, by largest remainder.
    count: i32,
}

pub struct SssDialog {
    open: bool,
    source: Structure,
    loaded: bool,
    label: String,
    entries: Vec<Entry>,
    picker: Option<usize>,
    seed: i32,
    status: (String, bool),
    browser: FileBrowser,
    /// The file the host came from (empty for the scene or the step's input).
    source_path: String,
    /// The pipeline step being edited ("Update Step" instead of Build).
    step: Option<StepEdit>,
}

/// std::stof: leading spaces, then the longest number prefix.
fn stof(text: &str) -> Option<f32> {
    let t = text.trim_start();
    (1..=t.len()).rev().filter(|&n| t.is_char_boundary(n)).find_map(|n| {
        let head = &t[..n];
        if head.ends_with(|c: char| c.is_ascii_alphabetic() && !matches!(c, 'e' | 'E')) {
            return None;
        }
        head.parse::<f32>().ok()
    })
}

impl Default for SssDialog {
    fn default() -> Self {
        Self { open: false, source: Structure::default(), loaded: false, label: String::new(), entries: Vec::new(), picker: None, seed: 12345, status: (String::new(), false), browser: FileBrowser::default(), source_path: String::new(), step: None }
    }
}

const GOOD: egui::Color32 = egui::Color32::from_rgb(70, 170, 90);
const BAD: egui::Color32 = egui::Color32::from_rgb(220, 80, 60);

impl SssDialog {
    fn load_scene(&mut self, scene: &Structure) {
        if scene.atoms.is_empty() {
            self.status = ("Current scene has no atoms.".into(), true);
            return;
        }
        self.source = scene.clone();
        self.loaded = true;
        self.label = "scene".into();
        self.source_path.clear();
        self.loaded_message();
    }

    fn load_path(&mut self, path: &str) -> bool {
        match af_io::loader::load_structure_from_file(path) {
            Ok(s) => {
                self.source = s;
                self.loaded = true;
                self.source_path = path.to_string();
                self.label = std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string());
                self.loaded_message();
                true
            }
            Err(e) => {
                self.status = (format!("Load failed: {}", if e.is_empty() { path.to_string() } else { e }), true);
                false
            }
        }
    }

    /// The composition of the host (initCompositionFromSource).
    fn loaded_message(&mut self) {
        let mut counts: std::collections::BTreeMap<i32, i32> = Default::default();
        for a in &self.source.atoms {
            *counts.entry(a.atomic_number).or_default() += 1;
        }
        let total = self.source.atoms.len() as f32;
        self.entries = counts.into_iter().map(|(z, c)| Entry { atomic_number: z, count: c, percent: if total > 0.0 { c as f32 / total * 100.0 } else { 0.0 } }).collect();
        self.status = (format!("Loaded '{}'  ({} atoms)", self.label, self.source.atoms.len()), false);
    }

    /// The dialog as the step's options (stepOptions): --input when the host
    /// came from a file, --frac with each element's share of the atoms, --seed.
    pub fn step_options(&self) -> String {
        use af_pipeline::options::Writer;
        let total = self.source.atoms.len();
        let percent_sum: f32 = self.entries.iter().map(|e| e.percent.max(0.0)).sum();
        let mut frac = String::new();
        for e in &self.entries {
            let fraction = if total > 0 {
                e.count as f64 / total as f64
            } else if percent_sum > 0.0 {
                e.percent.max(0.0) as f64 / percent_sum as f64
            } else {
                0.0
            };
            if !frac.is_empty() {
                frac.push(',');
            }
            frac += &format!("{}={}", elements::symbol(e.atomic_number), Writer::number(fraction));
        }
        let mut out = Writer::new();
        if !self.source_path.is_empty() {
            out.add("--input", &self.source_path);
        }
        if !frac.is_empty() {
            out.add("--frac", &frac);
        }
        out.add_integer("--seed", self.seed as i64);
        out.str().to_string()
    }

    /// "Update Step": the step with the dialog's options; the edit ends.
    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        self.open = false;
        Some(Outcome::UpdateStep(edit.with_options(&self.step_options())))
    }

    fn recompute_counts(&mut self) {
        let total = self.source.atoms.len() as i32;
        let n = self.entries.len();
        if n == 0 || total == 0 {
            return;
        }
        let sum: f32 = self.entries.iter().map(|e| e.percent.max(0.0)).sum();
        if sum <= 0.0 {
            let (base, rem) = (total / n as i32, total % n as i32);
            for (i, e) in self.entries.iter_mut().enumerate() {
                e.count = base + if (i as i32) < rem { 1 } else { 0 };
            }
            return;
        }
        let real: Vec<f32> = self.entries.iter().map(|e| e.percent.max(0.0) / sum * total as f32).collect();
        let floored: Vec<i32> = real.iter().map(|r| r.floor() as i32).collect();
        let remainders: Vec<f32> = real.iter().zip(&floored).map(|(r, f)| r - *f as f32).collect();
        let assigned: i32 = floored.iter().sum();
        let mut order: Vec<usize> = (0..n).collect();
        af_core::cxx::sort_by(&mut order, |&a, &b| remainders[a] > remainders[b]);
        let mut counts = floored;
        for &k in order.iter().take((total - assigned).max(0) as usize) {
            counts[k] += 1;
        }
        for (e, c) in self.entries.iter_mut().zip(counts) {
            e.count = c;
        }
    }

    fn snap_percents(&mut self) {
        let total = self.source.atoms.len() as f32;
        if total > 0.0 {
            for e in &mut self.entries {
                e.percent = e.count as f32 / total * 100.0;
            }
        }
    }

    /// An entry's percent changed: the last entry takes up the rest
    /// (redistributePercentLastFree); percents snap to whole atoms.
    fn set_percent_last_free(&mut self, changed: usize, new_pct: f32) {
        let n = self.entries.len();
        if n == 0 {
            return;
        }
        if n == 1 {
            self.entries[0].percent = 100.0;
            self.recompute_counts();
            self.snap_percents();
            return;
        }
        let mut new_pct = new_pct.max(0.0);
        let mut sum_except_last = 0.0f32;
        for i in 0..n - 1 {
            sum_except_last += if i == changed { new_pct } else { self.entries[i].percent.max(0.0) };
        }
        if sum_except_last > 100.0 {
            new_pct = (new_pct - (sum_except_last - 100.0)).max(0.0);
            sum_except_last = 100.0;
        }
        self.entries[changed].percent = new_pct;
        self.entries[n - 1].percent = 100.0 - sum_except_last;
        self.recompute_counts();
        self.snap_percents();
    }

    /// An entry's percent set, the others scaled to keep their ratios (redistributePercent).
    fn set_percent_scaled(&mut self, changed: usize, new_pct: f32) {
        let n = self.entries.len();
        if n == 0 {
            return;
        }
        let new_pct = new_pct.clamp(0.0, 100.0);
        if n == 1 {
            self.entries[0].percent = 100.0;
            self.recompute_counts();
            return;
        }
        let other_sum: f32 = self.entries.iter().enumerate().filter(|(i, _)| *i != changed).map(|(_, e)| e.percent.max(0.0)).sum();
        let remaining = 100.0 - new_pct;
        self.entries[changed].percent = new_pct;
        for (i, e) in self.entries.iter_mut().enumerate() {
            if i == changed {
                continue;
            }
            e.percent = if remaining <= 0.0 {
                0.0
            } else if other_sum <= 0.0 {
                remaining / (n - 1) as f32
            } else {
                e.percent.max(0.0) * (remaining / other_sum)
            };
        }
        self.recompute_counts();
    }

    fn add_entry(&mut self) {
        let z = (1..=118).find(|z| !self.entries.iter().any(|e| e.atomic_number == *z)).unwrap_or(1);
        self.entries.push(Entry { atomic_number: z, percent: 0.0, count: 0 });
        self.set_percent_scaled(self.entries.len() - 1, 10.0);
    }

    fn remove_entry(&mut self, i: usize) {
        if i < self.entries.len() {
            self.entries.remove(i);
            self.recompute_counts();
        }
    }

    fn build(&mut self) -> Option<Outcome> {
        let total = self.source.atoms.len();
        let params = SssParams {
            seed: self.seed.max(0) as u32,
            composition: self.entries.iter().map(|e| ElementFraction { atomic_number: e.atomic_number, fraction: if total > 0 { e.count as f32 / total as f32 } else { 0.0 } }).collect(),
        };
        let r = build_substitutional_solid_solution(&self.source, &params);
        self.status = (r.message.clone(), !r.success);
        r.success.then(|| Outcome::Built(r.output, "Built a substitutional solid solution".into()))
    }

    fn host_panel(&mut self, ui: &mut egui::Ui, scene: &Structure) {
        ui.horizontal(|ui| {
            ui.strong("Host Structure");
            if let Some(input) = self.step.as_ref().map(|e| e.input.clone()) {
                // The step's host is the structure entering it, not the scene.
                if ui.button("Use Step Input").clicked() {
                    self.load_scene(&input);
                    if self.loaded && self.source_path.is_empty() {
                        self.label = "pipeline step input".into();
                    }
                }
            } else if ui.button("Use Current Scene").clicked() {
                self.load_scene(scene);
            }
            if ui.button("Load File...").clicked() {
                self.browser.show_for(Purpose::Open, "Load a host structure", "", &[]);
            }
        });
        if self.loaded {
            ui.horizontal(|ui| {
                ui.colored_label(GOOD, &self.label);
                ui.weak(format!("({} atoms)", self.source.atoms.len()));
            });
            if self.source.has_unit_cell {
                let len = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                let c = &self.source.cell_vectors;
                ui.weak(format!("Cell: a={:.2}  b={:.2}  c={:.2} Å", len(c[0]), len(c[1]), len(c[2])));
            }
        } else {
            ui.weak("No structure loaded yet.");
        }
    }

    fn composition_panel(&mut self, ui: &mut egui::Ui, colours: &[[f32; 3]]) -> Option<Outcome> {
        if !self.loaded {
            ui.add_space(4.0);
            ui.weak("Load a host structure above to");
            ui.weak("configure the target composition.");
            return None;
        }
        let mut built = None;
        ui.horizontal(|ui| {
            ui.strong("Target Composition");
            ui.weak(format!("({} total atoms)", self.source.atoms.len()));
        });
        let last = self.entries.len().saturating_sub(1);
        let mut remove = None;
        let mut changed = None;
        egui::Grid::new("sss composition").num_columns(4).striped(true).show(ui, |ui| {
            ui.weak("El.");
            ui.weak("at%");
            ui.weak("atoms");
            ui.end_row();
            let many = self.entries.len() > 1;
            for (i, e) in self.entries.iter_mut().enumerate() {
                if ui.add(egui::Button::new(elements::symbol(e.atomic_number)).min_size(egui::vec2(48.0, 0.0))).on_hover_text(elements::name(e.atomic_number)).clicked() {
                    self.picker = Some(i);
                }
                let mut pct = e.percent;
                let drag = egui::DragValue::new(&mut pct).speed(0.15).range(0.0..=100.0).fixed_decimals(2).suffix(" at%");
                if i == last {
                    ui.add_enabled(false, drag).on_disabled_hover_text("Auto-computed: 100 - sum of all other entries.");
                } else if ui.add(drag).on_hover_text("Drag or click to edit.\nThe last element adjusts automatically.").changed() {
                    changed = Some((i, pct));
                }
                ui.label(format!("{} atoms", e.count));
                if many && ui.button("X").clicked() {
                    remove = Some(i);
                }
                ui.end_row();
            }
        });
        if let Some((i, pct)) = changed {
            self.set_percent_last_free(i, pct);
        }
        if let Some(i) = remove {
            self.remove_entry(i);
        }
        if self.entries.len() < 118 && ui.button("+ Add Element").clicked() {
            self.add_entry();
        }
        if let Some(i) = self.picker {
            let mut z = self.entries.get(i).map_or(1, |e| e.atomic_number);
            let mut keep = true;
            crate::theme::modal_window("Select Element").id(egui::Id::unique("sss element picker")).open(&mut keep).collapsible(false).resizable(false).show(ui.ctx(), |ui| {
                if periodic_table(ui, &mut z, Some(colours)) {
                    if let Some(e) = self.entries.get_mut(i) {
                        e.atomic_number = z;
                    }
                    self.picker = None;
                }
            });
            if !keep {
                self.picker = None;
            }
        }
        ui.separator();
        ui.strong("Options");
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.seed).range(0..=i32::MAX)).on_hover_text("A positive fixed integer makes the site assignment reproducible for the same input ordering. Zero selects a random/time-dependent assignment; negative input is clamped to zero.");
            ui.label("RNG Seed");
            ui.weak("(0 = random)");
        });
        ui.separator();
        if self.step.is_some() {
            if ui.add(egui::Button::new(egui::RichText::new("Update Step").strong()).min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                built = self.update_step();
            }
        } else if ui.add(egui::Button::new(egui::RichText::new("Build Solid Solution").strong()).min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
            built = self.build();
        }
        built
    }
}

impl StructureDialog for SssDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        self.load_path(path);
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-sss"]
    }

    fn title(&self) -> &'static str {
        "Substitutional Solid Solution"
    }

    fn open(&mut self, _step: &str) {
        self.step = None;
        self.open = true;
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        let reader = edit.reader();
        self.status = (String::new(), false);
        // Host: --input FILE, else the structure entering the step.
        self.source = Structure::default();
        self.loaded = false;
        self.label.clear();
        self.source_path.clear();
        let mut failed = None;
        if reader.has("--input") {
            let path = reader.text("--input", "");
            if !self.load_path(&path) {
                self.source = Structure::default();
                failed = Some(self.status.clone());
            }
        } else if !edit.input.atoms.is_empty() {
            self.source = edit.input.clone();
            self.loaded = true;
            self.label = "pipeline step input".into();
        }
        self.loaded_message();
        self.status = failed.unwrap_or_default();
        // --frac "Cu=0.7,Zn=0.3" sets the composition (as atomic percents).
        let mut entries = Vec::new();
        for token in reader.text("--frac", "").split(',') {
            let Some((symbol, value)) = token.split_once('=') else { continue };
            let symbol: String = symbol.chars().filter(|c| !c.is_whitespace()).collect();
            let z = elements::atomic_number(&symbol);
            if z <= 0 {
                continue;
            }
            let Some(v) = stof(value) else { continue };
            entries.push(Entry { atomic_number: z, percent: v * 100.0, count: 0 });
        }
        if !entries.is_empty() {
            self.entries = entries;
            self.recompute_counts();
        }
        self.seed = reader.integer("--seed", 12345).max(0);
        self.step = Some(edit);
        self.open = true;
        true
    }


    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if let Some(path) = self.browser.show(ctx) {
            self.browser.open = false;
            self.load_path(&path.to_string_lossy());
        }
        if !self.open {
            return None;
        }
        let mut out = None;
        let mut open = self.open;
        let editing = self.step.as_ref().map(|e| e.step.clone());
        crate::theme::modal_window("Substitutional Solid Solution").open(&mut open).default_size([560.0, 560.0]).min_width(480.0).show(ctx, |ui| {
            if let Some(step) = &editing {
                ui.colored_label(egui::Color32::from_rgb(110, 170, 230), format!("Editing pipeline step: {step}"));
            }
            ui.weak("Substitute elements on host lattice sites to match a target composition.  Positions and unit cell are preserved.");
            ui.separator();
            self.host_panel(ui, input.source);
            ui.separator();
            egui::ScrollArea::vertical().max_height(420.0).min_scrolled_height(260.0).show(ui, |ui| {
                out = self.composition_panel(ui, input.colours);
            });
            ui.separator();
            if !self.status.0.is_empty() {
                ui.colored_label(if self.status.1 { BAD } else { GOOD }, &self.status.0);
            }
        });
        self.open = open && self.open;
        if out.is_some() {
            self.open = false;
        }
        if !self.open {
            // Closed: the next menu build is normal.
            self.step = None;
        }
        out
    }

    /// use-scene | load PATH | element ROW SYMBOL | percent ROW PCT | add |
    /// remove ROW | seed N | build | update-step | expect-counts SYMBOL=N ... |
    /// expect-options TEXT (rows from 1).
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let row = |i: usize| -> Result<usize, String> { words.get(i).and_then(|w| w.parse::<usize>().ok()).filter(|&r| r >= 1).map(|r| r - 1).ok_or("a row number (from 1) is needed".to_string()) };
        match command {
            "use-scene" => {
                self.load_scene(input.source);
                if !self.loaded {
                    return Err(self.status.0.clone());
                }
            }
            "load" => {
                if !self.load_path(args.trim()) {
                    return Err(self.status.0.clone());
                }
            }
            "element" => {
                let i = row(0)?;
                let z = elements::atomic_number(words.get(1).copied().unwrap_or(""));
                let e = self.entries.get_mut(i).ok_or("no such row")?;
                if z <= 0 {
                    return Err(format!("no element {args}"));
                }
                e.atomic_number = z;
            }
            "percent" => {
                let i = row(0)?;
                let pct = words.get(1).and_then(|w| w.parse::<f32>().ok()).ok_or("percent ROW PCT")?;
                if i >= self.entries.len() {
                    return Err("no such row".into());
                }
                self.set_percent_last_free(i, pct);
            }
            "add" => self.add_entry(),
            "remove" => self.remove_entry(row(0)?),
            "seed" => self.seed = args.trim().parse::<i32>().map_err(|e| e.to_string())?.max(0),
            "update-step" => return self.update_step().map(Some).ok_or_else(|| "no pipeline step is being edited".to_string()),
            "expect-options" => {
                if self.step_options() != args.trim() {
                    return Err(format!("options are {:?}", self.step_options()));
                }
            }
            "build" => {
                let o = self.build();
                return if o.is_some() { Ok(o) } else { Err(self.status.0.clone()) };
            }
            "expect-counts" => {
                for w in &words {
                    let (sym, n) = w.split_once('=').ok_or("expect-counts SYMBOL=N ...")?;
                    let z = elements::atomic_number(sym);
                    let got = self.entries.iter().find(|e| e.atomic_number == z).map_or(0, |e| e.count);
                    if got.to_string() != n {
                        return Err(format!("{sym} has {got} atoms, expected {n}"));
                    }
                }
            }
            "expect-status" => {
                if !self.status.0.contains(args.trim()) {
                    return Err(format!("status is {:?}", self.status.0));
                }
            }
            other => return Err(format!("the solid solution dialog has no command {other}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::json::Json;
    use af_core::structure::AtomSite;

    fn edit(options: &str, total: usize) -> StepEdit {
        let mut parameters = Json::object();
        parameters.set("options", options);
        let mut input = Structure::default();
        input.atoms = (0..total).map(|i| AtomSite::new("Cu", [i as f64, 0.0, 0.0])).collect();
        StepEdit { index: 0, step: "build-sss".into(), parameters, input, ..Default::default() }
    }

    /// editStep then stepOptions, against the C++ (golden_step_options_composition).
    #[test]
    fn step_options_match_the_cpp() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/step_options_composition.json");
        let golden = Json::parse(&std::fs::read_to_string(path).unwrap()).unwrap();
        let cases = golden.find("sss").unwrap().items().unwrap();
        assert!(cases.len() >= 8);
        for case in cases {
            let options = case.find("options").unwrap().string().unwrap();
            let total = case.find("total").unwrap().number().unwrap() as usize;
            let mut d = SssDialog::default();
            assert!(d.edit_step(edit(options, total)));
            assert_eq!(d.step_options(), case.find("out").unwrap().string().unwrap(), "options {options:?}, {total} atoms");
        }
    }

    /// The options written read back to the same dialog (and Update Step ends the edit).
    #[test]
    fn options_round_trip() {
        let mut d = SssDialog::default();
        d.edit_step(edit("--frac Cu=0.6,Ni=0.25,Fe=0.15 --seed 9", 40));
        let text = d.step_options();
        assert_eq!(text, "--frac Cu=0.6,Ni=0.25,Fe=0.15 --seed 9");
        let mut again = SssDialog::default();
        again.edit_step(edit(&text, 40));
        assert_eq!(again.entries.iter().map(|e| e.count).collect::<Vec<_>>(), vec![24, 10, 6]);
        assert_eq!(again.step_options(), text);
        match again.update_step() {
            Some(Outcome::UpdateStep(e)) => assert_eq!(e.options(), text),
            _ => panic!("no step update"),
        }
        assert!(again.step.is_none() && again.update_step().is_none());
    }
}
