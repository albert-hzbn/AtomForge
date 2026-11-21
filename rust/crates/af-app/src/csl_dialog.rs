//! Build > CSL Grain Boundary (the C++ ui/CSLGrainBoundaryDialog): a
//! reference crystal, the rotation axis and the coincidence (Σ) candidates up
//! to a maximum Σ, one of the three CSL boundary planes, the grain
//! thicknesses, vacuum, gap and overlap removal; the bicrystal replaces the
//! structure in the view.

use af_algorithms::csl::{build_bicrystal_from_candidate, classify_boundary_type, compute_gb_info, SigmaCandidate};
use af_core::structure::Structure;
use eframe::egui;

/// The structure files a reference can be loaded from (as File > Open).
const STRUCTURE_EXTENSIONS: &[&str] = &["cif", "xyz", "extxyz", "vasp", "poscar", "contcar", "pdb", "mol", "sdf", "mol2", "gjf", "xsf", "pwi", "in", "lmp", "data", "cell", "gen", "json"];

use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

/// The C++ CSLBuildResult: what the last Build made, or why not.
#[derive(Default)]
struct BuildResult {
    success: bool,
    message: String,
    boundary_type: String,
    sigma: i32,
    theta_deg: f32,
    input_atoms: usize,
    generated_atoms: usize,
    removed_overlap: i32,
    axis: [i32; 3],
    plane: [i32; 3],
    csl: [[i32; 3]; 3],
}

pub struct CslDialog {
    open: bool,
    /// The input structure's preview (drag to orbit, scroll to zoom).
    orbit: crate::builder_preview::Orbit,
    /// Load Reference: the file picker.
    browser: crate::file_browser::FileBrowser,
    /// Take the structure in the view as reference the next frame.
    take_view: bool,
    /// Editing a pipeline step: Update Step instead of Build.
    step: Option<StepEdit>,
    reference: Structure,
    /// The file the reference came from ("" = none, the view or the step input).
    reference_path: String,
    status: String,
    file_name: String,
    axis: [i32; 3],
    sigma_max: i32,
    candidates: Vec<SigmaCandidate>,
    sigma_selection: usize,
    plane_selection: usize,
    last_axis: [i32; 3],
    last_sigma_max: i32,
    uc_a: i32,
    uc_b: i32,
    vacuum: f32,
    gap: f32,
    overlap: f32,
    conventional: bool,
    result: BuildResult,
}

impl Default for CslDialog {
    fn default() -> Self {
        Self {
            open: false,
            take_view: false,
            orbit: crate::builder_preview::Orbit::default(),
            browser: crate::file_browser::FileBrowser::default(),
            step: None,
            reference: Structure::default(),
            reference_path: String::new(),
            status: "(no structure loaded)".into(),
            file_name: "(none)".into(),
            axis: [0, 0, 1],
            sigma_max: 200,
            candidates: Vec::new(),
            sigma_selection: 0,
            plane_selection: 0,
            last_axis: [0, 0, 0],
            last_sigma_max: 0,
            uc_a: 1,
            uc_b: 1,
            vacuum: 0.0,
            gap: 0.0,
            overlap: 0.0,
            conventional: false,
            result: BuildResult::default(),
        }
    }
}

const THICKNESS_HINT: &str = "Set repetitions used to thicken each grain along the stacking construction. Use both controls to study symmetric or unequal thicknesses. These are integer cell repetitions, not distances in Å.";
const SPACING_HINT: &str = "Vacuum adds the specified outside spacing in the construction; Gap separates the stacked grains. Both alter geometry. Neither is a relaxation or a search over rigid translations.";

fn candidate_label(c: &SigmaCandidate) -> String {
    format!("Σ{}  m={} n={}  {:.2} deg", c.sigma, c.m, c.n, c.theta_deg)
}

impl CslDialog {
    fn set_reference(&mut self, s: Structure, file_name: &str, status: String) {
        if !s.has_unit_cell {
            self.status = "Error: structure has no unit cell".into();
            self.file_name = "(none)".into();
            return;
        }
        self.reference = s;
        self.status = status;
        self.file_name = file_name.into();
    }

    fn use_view(&mut self, source: &Structure) {
        if source.atoms.is_empty() {
            self.status = "(no structure loaded)".into();
            return;
        }
        self.set_reference(source.clone(), "(structure in the view)", format!("Structure in view: {} atoms", source.atoms.len()));
        self.reference_path.clear();
    }

    fn load(&mut self, path: &str) -> Result<(), String> {
        match af_io::loader::load_structure_from_file(path) {
            Ok(s) => {
                let name = std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.into());
                let n = s.atoms.len();
                self.set_reference(s, &name, format!("Loaded: {n} atoms"));
                if self.reference.atoms.is_empty() || self.file_name == "(none)" {
                    return Err(self.status.clone());
                }
                self.reference_path = path.into();
                Ok(())
            }
            Err(e) => {
                self.status = format!("Error: {e}");
                self.file_name = "(none)".into();
                Err(self.status.clone())
            }
        }
    }

    /// The dialog's settings as the step's options (the flags of AtomForge
    /// --build gb; the C++ stepOptions).
    pub fn step_options(&self) -> String {
        use af_pipeline::options::Writer;
        let mut out = Writer::new();
        if !self.reference_path.is_empty() {
            out.add("--input", &self.reference_path);
        }
        out.add_numbers("--axis", &[self.axis[0] as f64, self.axis[1] as f64, self.axis[2] as f64]);
        if let Some(c) = self.candidates.get(self.sigma_selection) {
            out.add_integer("--sigma", c.sigma as i64);
        }
        out.add_integer("--sigmamax", self.sigma_max as i64);
        out.add_integer("--plane", self.plane_selection as i64);
        out.add_integer("--uca", self.uc_a as i64);
        out.add_integer("--ucb", self.uc_b as i64);
        out.add_number("--vacuum", self.vacuum as f64);
        out.add_number("--gap", self.gap as f64);
        out.add_number("--overlap", self.overlap as f64);
        if self.conventional {
            out.flag("--conventional");
        }
        out.str().to_string()
    }

    /// Fills the dialog from a pipeline step (the C++ editStep): missing flags
    /// take the defaults of AtomForge --build gb.
    fn fill_from_step(&mut self, edit: &StepEdit) {
        let o = edit.reader();
        let axis = o.numbers("--axis");
        for i in 0..3 {
            self.axis[i] = if axis.len() == 3 { axis[i].round() as i32 } else if i == 2 { 1 } else { 0 };
        }
        self.sigma_max = o.integer("--sigmamax", 100).max(3);
        self.plane_selection = o.integer("--plane", 0).clamp(0, 2) as usize;
        self.uc_a = o.integer("--uca", 1).max(1);
        self.uc_b = o.integer("--ucb", 1).max(1);
        self.vacuum = o.number("--vacuum", 0.0) as f32;
        self.gap = o.number("--gap", 0.0) as f32;
        self.overlap = o.number("--overlap", 0.0) as f32;
        self.conventional = o.has("--conventional");
        // The Σ list of the axis, with the step's Σ selected (none = the smallest).
        self.candidates = compute_gb_info(&self.axis, self.sigma_max);
        self.last_axis = self.axis;
        self.last_sigma_max = self.sigma_max;
        let sigma = o.integer("--sigma", 0);
        self.sigma_selection = self.candidates.iter().position(|c| c.sigma == sigma).unwrap_or(0);
        // Reference: a file named in the options, else the structure entering the step.
        self.reference = Structure::default();
        self.reference_path = o.text("--input", "");
        self.status = "(no structure loaded)".into();
        self.file_name = "(none)".into();
        if !self.reference_path.is_empty() {
            let path = self.reference_path.clone();
            if self.load(&path).is_err() {
                self.reference_path = path;
            }
        } else if edit.use_input() && !edit.input.atoms.is_empty() {
            self.reference = edit.input.clone();
            self.status = format!("Step input: {} atoms", self.reference.atoms.len());
            self.file_name = "(structure entering the step)".into();
        }
        self.result = BuildResult::default();
    }

    /// Update Step: the options, and whether the pipeline passes its input
    /// (no --input) as the reference.
    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        self.open = false;
        Some(Outcome::UpdateStep(edit.with_options(&self.step_options()).with_parameter("use_input", self.reference_path.is_empty())))
    }

    /// The Σ list is cleared when the axis or the maximum changes (as in C++).
    fn sync_candidates(&mut self) {
        self.sigma_max = self.sigma_max.max(3);
        if self.axis != self.last_axis || self.sigma_max != self.last_sigma_max {
            self.candidates.clear();
            self.sigma_selection = 0;
            self.plane_selection = 0;
            self.last_axis = self.axis;
            self.last_sigma_max = self.sigma_max;
        }
    }

    fn find(&mut self) {
        self.sync_candidates();
        self.candidates = compute_gb_info(&self.axis, self.sigma_max);
        self.sigma_selection = 0;
        self.plane_selection = 0;
    }

    /// Build: the bicrystal, or the reason in the result.
    fn build(&mut self) -> Option<Structure> {
        self.sync_candidates();
        self.result = BuildResult::default();
        if !self.reference.has_unit_cell || self.reference.atoms.is_empty() {
            self.result.message = "Please load a structure with a unit cell first.".into();
            return None;
        }
        if self.candidates.is_empty() {
            self.result.message = "Generate and select a Sigma candidate first.".into();
            return None;
        }
        let sel = self.candidates[self.sigma_selection.min(self.candidates.len() - 1)];
        let plane_index = self.plane_selection.min(2);
        let b = build_bicrystal_from_candidate(&self.reference, &sel, plane_index, self.uc_a.max(1), self.uc_b.max(1), self.vacuum, self.gap, self.overlap, self.conventional);
        let plane = sel.plane[plane_index];
        self.result = BuildResult {
            success: true,
            message: "CSL grain boundary generated.".into(),
            boundary_type: classify_boundary_type(&self.axis, &plane).into(),
            sigma: sel.sigma,
            theta_deg: sel.theta_deg,
            input_atoms: self.reference.atoms.len(),
            generated_atoms: b.structure.atoms.len(),
            removed_overlap: b.removed,
            axis: self.axis,
            plane,
            csl: sel.csl,
        };
        Some(b.structure)
    }

    fn outcome(s: Structure) -> Outcome {
        Outcome::Built(s, "Built a CSL grain boundary".into())
    }

    fn result_summary(&self, ui: &mut egui::Ui) {
        let r = &self.result;
        if r.message.is_empty() {
            return;
        }
        ui.separator();
        ui.strong("Result");
        if !r.success {
            ui.colored_label(egui::Color32::from_rgb(220, 80, 60), format!("Error: {}", r.message));
            return;
        }
        ui.colored_label(egui::Color32::from_rgb(70, 170, 90), format!("OK: {}", r.message));
        egui::Grid::new("csl result").num_columns(2).striped(false).show(ui, |ui| {
            let mut row = |label: &str, value: String| {
                ui.weak(label);
                ui.label(value);
                ui.end_row();
            };
            row("Input atoms", r.input_atoms.to_string());
            row("Output atoms", r.generated_atoms.to_string());
            row("Removed overlaps", r.removed_overlap.to_string());
            row("Sigma", r.sigma.to_string());
            row("Misorientation angle", format!("{:.4} deg", r.theta_deg));
            row("Boundary type", if r.boundary_type.is_empty() { "unknown".into() } else { r.boundary_type.clone() });
            row("Axis [u v w]", format!("[{} {} {}]", r.axis[0], r.axis[1], r.axis[2]));
            row("GB plane (h k l)", format!("({} {} {})", r.plane[0], r.plane[1], r.plane[2]));
            let c = &r.csl;
            row("CSL matrix", format!("[{} {} {}]  [{} {} {}]  [{} {} {}]", c[0][0], c[0][1], c[0][2], c[1][0], c[1][1], c[1][2], c[2][0], c[2][1], c[2][2]));
        });
    }

    fn controls(&mut self, ui: &mut egui::Ui, input: &Input) -> Option<Outcome> {
        use crate::dialog_layout as dl;
        let mut built = None;
        if let Some(edit) = &self.step {
            crate::structure_dialog::step_banner(ui, &edit.step);
        }
        ui.horizontal(|ui| {
            if ui.button("Load Reference").clicked() {
                self.browser.show_for(crate::file_browser::Purpose::Open, "Load reference", "", STRUCTURE_EXTENSIONS);
            }
            if ui.add_enabled(!input.source.atoms.is_empty(), egui::Button::new("Use Structure in View")).clicked() {
                self.use_view(input.source);
            }
        });
        ui.add_space(4.0);
        let height = 700.0f32.min(ui.available_height() - 60.0).max(360.0);
        ui.horizontal_top(|ui| {
            // Left: the input structure and its preview.
            dl::panel(ui, egui::vec2(380.0, height), |ui| {
                ui.label("Input Structure");
                ui.separator();
                ui.label(format!("Status: {}", self.status));
                dl::muted(ui, format!("File: {}", self.file_name));
                ui.add_space(4.0);
                let size = egui::vec2(ui.available_width(), (height - 190.0).max(80.0));
                let reference = if self.reference.has_unit_cell { self.reference.clone() } else { Structure::default() };
                crate::builder_preview::structure(ui, size, &reference, input.colours, &mut self.orbit, &["Drop a structure file here"]);
                if self.reference.has_unit_cell && !self.reference.atoms.is_empty() {
                    let len = |r: [f64; 3]| (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
                    let cv = &self.reference.cell_vectors;
                    dl::muted(ui, format!("{} atoms  a={:.2} b={:.2} c={:.2}", self.reference.atoms.len(), len(cv[0]), len(cv[1]), len(cv[2])));
                    dl::muted(ui, "Drag=orbit  Scroll=zoom");
                }
            });
            // Right: the settings.
            let rest = ui.available_width() - 2.0;
            dl::panel(ui, egui::vec2(rest, height), |ui| {
                egui::ScrollArea::vertical().id_salt("csl settings").show(ui, |ui| {
                    dl::separator_text(ui, "Misorientation");
                    ui.horizontal(|ui| {
                        dl::int_fields(ui, "csl axis", &mut self.axis, 220.0, -50..=50);
                        ui.label("Axis [u v w]").on_hover_text("Crystallographic rotation axis in Miller indices.");
                        dl::int_input(ui, "csl sigma max", &mut self.sigma_max, 110.0, 3..=10_000);
                        ui.label("Max sigma");
                        if ui.button("Find").clicked() {
                            self.find();
                        }
                    });
                    self.sync_candidates();
                    if self.candidates.is_empty() {
                        dl::muted(ui, "Set axis and press Find.");
                    } else {
                        if self.sigma_selection >= self.candidates.len() {
                            self.sigma_selection = 0;
                        }
                        egui::ComboBox::from_id_salt("csl sigma").width(ui.available_width() - 8.0).selected_text(candidate_label(&self.candidates[self.sigma_selection])).show_ui(ui, |ui| {
                            for (i, c) in self.candidates.iter().enumerate() {
                                ui.selectable_value(&mut self.sigma_selection, i, candidate_label(c));
                            }
                        });
                    }
                    dl::separator_text(ui, "Boundary Geometry");
                    if let Some(sel) = self.candidates.get(self.sigma_selection).copied() {
                        let labels: Vec<String> = (0..3).map(|i| format!("({} {} {})  {}", sel.plane[i][0], sel.plane[i][1], sel.plane[i][2], classify_boundary_type(&self.axis, &sel.plane[i]))).collect();
                        if self.plane_selection >= 3 {
                            self.plane_selection = 0;
                        }
                        ui.horizontal(|ui| {
                            egui::ComboBox::from_id_salt("csl gb plane").width(260.0).selected_text(&labels[self.plane_selection]).show_ui(ui, |ui| {
                                for (i, l) in labels.iter().enumerate() {
                                    ui.selectable_value(&mut self.plane_selection, i, l);
                                }
                            })
                            .response
                            .on_hover_text("Select a coincidence candidate, then a boundary plane from the generated choices. The plane selector contains candidate planes; the CLI index 0, 1, or 2 selects one of these choices. It is not the value of h, k, or l.\n\nGB plane normals derived from CSL matrix columns.\nTwist = plane parallel to axis, Tilt = plane perpendicular to axis.");
                            ui.label("GB plane");
                        });
                        dl::muted(ui, format!("Type: {}", classify_boundary_type(&self.axis, &sel.plane[self.plane_selection])));
                        let c = &sel.csl;
                        dl::muted(ui, format!("CSL: [{} {} {}] [{} {} {}] [{} {} {}]", c[0][0], c[0][1], c[0][2], c[1][0], c[1][1], c[1][2], c[2][0], c[2][1], c[2][2]));
                    } else {
                        dl::muted(ui, "Select a Sigma candidate first.");
                    }
                    ui.horizontal(|ui| {
                        dl::int_input(ui, "csl uc a", &mut self.uc_a, 110.0, 1..=100);
                        ui.label("Grain A (uc)").on_hover_text(THICKNESS_HINT);
                        dl::int_input(ui, "csl uc b", &mut self.uc_b, 110.0, 1..=100);
                        ui.label("Grain B (uc)").on_hover_text(format!("{THICKNESS_HINT}\n\nNumber of unit cells for each grain along the GB normal."));
                    });
                    dl::separator_text(ui, "Options");
                    ui.horizontal(|ui| {
                        let mut v = self.vacuum as f64;
                        if dl::float_input(ui, "csl vacuum", &mut v, 0.5, 2, 140.0, 0.0..=50.0) {
                            self.vacuum = v as _;
                        }
                        ui.label("Vacuum (Å)").on_hover_text(format!("{SPACING_HINT}\n\nVacuum space at the end of the bicrystal cell."));
                    });
                    ui.horizontal(|ui| {
                        let mut v = self.gap as f64;
                        if dl::float_input(ui, "csl gap", &mut v, 0.1, 2, 140.0, 0.0..=10.0) {
                            self.gap = v as _;
                        }
                        ui.label("Gap (Å)").on_hover_text(format!("{SPACING_HINT}\n\nExtra distance between the two grains at the interface."));
                    });
                    ui.horizontal(|ui| {
                        let mut v = self.overlap as f64;
                        if dl::float_input(ui, "csl overlap", &mut v, 0.1, 3, 140.0, 0.0..=5.0) {
                            self.overlap = v as _;
                        }
                        ui.label("Overlap dist (Å)").on_hover_text("Remove interfacial close contacts within the chosen cutoff. Inspect the removed-atom count and resulting stoichiometry. Set a conservative cutoff and compare before raising it.\n\nAtoms closer than this distance will be removed (0 = off).");
                    });
                    dl::labelled_check(ui, &mut self.conventional, "Conventional cell").on_hover_text("Retain the conventional construction when enabled; otherwise the builder attempts its primitive grain-boundary cell reduction. This can change the final cell representation and atom count.\n\nIf unchecked (default), reduce to the primitive cell.\nIf checked, keep the conventional (unreduced) cell.");
                    self.result_summary(ui);
                });
            });
        });
        let primary = if self.step.is_some() { "Update Step" } else { "Build" };
        match dl::footer(ui, primary, Some("Close"), true) {
            1 if self.step.is_some() => built = self.update_step(),
            1 => built = self.build().map(Self::outcome),
            2 => {
                self.open = false;
                self.step = None;
            }
            _ => {}
        }
        built
    }
}

impl StructureDialog for CslDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        let _ = self.load(path);
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-gb"]
    }

    fn title(&self) -> &'static str {
        "CSL Grain Boundary Builder"
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
        if self.reference.atoms.is_empty() {
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
        if std::mem::take(&mut self.take_view) && input.source.has_unit_cell {
            self.use_view(input.source);
        }
        let mut built = None;
        let mut open = self.open;
        crate::theme::modal_window(self.title()).open(&mut open).default_size([1100.0, 820.0]).show(ctx, |ui| {
            built = self.controls(ui, input);
        });
        if let Some(path) = self.browser.show(ctx) {
            let _ = self.load(&path.to_string_lossy());
        }
        self.open = self.open && open;
        if !self.open {
            self.step = None;
        }
        built
    }

    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let a = args.trim();
        let int = |w: &str| w.parse::<i32>().map_err(|_| format!("{command}: {w:?} is not an integer"));
        let float = || a.parse::<f32>().map_err(|_| format!("{command}: {a:?} is not a number"));
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
                self.use_view(input.source);
                if self.reference.atoms.is_empty() {
                    return Err(self.status.clone());
                }
            }
            "load" => self.load(a)?,
            "axis" => {
                let v: Vec<i32> = a.split_whitespace().map(int).collect::<Result<_, _>>()?;
                if v.len() != 3 {
                    return Err("axis U V W".into());
                }
                self.axis = [v[0], v[1], v[2]];
            }
            "sigma-max" => self.sigma_max = int(a)?,
            "find" => self.find(),
            "sigma" => {
                // sigma N [K]: the K-th (from 1) candidate with Σ = N.
                let w: Vec<&str> = a.split_whitespace().collect();
                let sigma = int(w.first().copied().unwrap_or(""))?;
                let k = w.get(1).map(|w| int(w)).transpose()?.unwrap_or(1).max(1) as usize;
                self.sigma_selection = self.candidates.iter().enumerate().filter(|(_, c)| c.sigma == sigma).nth(k - 1).map(|(i, _)| i).ok_or(format!("no candidate Σ{sigma} (#{k}); press find first"))?;
            }
            "plane" => self.plane_selection = (int(a)?).clamp(0, 2) as usize,
            "uca" => self.uc_a = int(a)?.max(1),
            "ucb" => self.uc_b = int(a)?.max(1),
            "vacuum" => self.vacuum = float()?,
            "gap" => self.gap = float()?,
            "overlap" => self.overlap = float()?,
            "conventional" => self.conventional = a == "on",
            "build" => {
                if self.step.is_some() {
                    return Err("editing a pipeline step: use update-step".into());
                }
                return match self.build() {
                    Some(s) => Ok(Some(Self::outcome(s))),
                    None => Err(self.result.message.clone()),
                };
            }
            "expect-candidates" => {
                let n = int(a)? as usize;
                if self.candidates.len() != n {
                    return Err(format!("{} Σ candidates: {:?}", self.candidates.len(), self.candidates.iter().map(|c| c.sigma).collect::<Vec<_>>()));
                }
            }
            "expect-result" => {
                // expect-result SIGMA TYPE REMOVED: the last build's Σ, boundary type and removed overlaps.
                let w: Vec<&str> = a.split_whitespace().collect();
                let r = &self.result;
                let got = format!("{} {} {}", r.sigma, r.boundary_type, r.removed_overlap);
                if !r.success || w.len() != 3 || got != w.join(" ") {
                    return Err(format!("last build: {got} ({})", r.message));
                }
            }
            _ => return Err(format!("unknown CSL dialog command {command}")),
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
        StepEdit { index: 0, step: "build-gb".into(), parameters, input: Structure::default(), ..Default::default() }
    }

    #[test]
    fn step_options_round_trip() {
        // The C++ stepOptions for the same state (see the golden program).
        for text in [
            "--axis \"0 0 1\" --sigma 5 --sigmamax 100 --plane 0 --uca 1 --ucb 1 --vacuum 0 --gap 0 --overlap 0",
            "--axis \"1 1 1\" --sigma 7 --sigmamax 50 --plane 2 --uca 2 --ucb 3 --vacuum 1.5 --gap 0.25 --overlap 0.5 --conventional",
        ] {
            let mut d = CslDialog::default();
            assert!(d.edit_step(edit(text)));
            assert_eq!(d.step_options(), text);
        }
        // Missing flags take the CLI defaults; a float keeps its float rounding.
        let mut d = CslDialog::default();
        d.edit_step(edit("--vacuum 0.1"));
        assert_eq!(d.step_options(), "--axis \"0 0 1\" --sigma 5 --sigmamax 100 --plane 0 --uca 1 --ucb 1 --vacuum 0.1000000015 --gap 0 --overlap 0");
    }
}
