//! Analysis > Structure and Defects > Dislocation Post-Processing (the C++
//! ui/DislocationAnalysisDialog): analyses comparing a deformed structure
//! with its undeformed reference, each a native reimplementation of a BABEL
//! program: the Nye tensor, the Vitek differential-displacement map,
//! crystallographic pattern matching, drag (migration-barrier) preparation
//! and the dislocation position/Burgers-vector fit. Results are exported to
//! files rather than written back into the scene.

use std::fmt::Write as _;
use std::sync::mpsc::{channel, Receiver};

use af_algorithms::nano_pattern::{build_pattern, detect_pattern, CrystalPattern, PatternDetectResult};
use af_algorithms::nye::{compute_nye_tensor, compute_vitek_map, fit_dislocation_from_nye, DislocationFitResult, NyeTensorResult, VitekResult};
use af_algorithms::sqs_drag::{prepare_drag, DragPrepResult};
use af_core::structure::Structure;
use eframe::egui;
use glam::DVec3;

use crate::analysis_panels::{empty_state, file_name, g, section, two_panels};
use crate::file_browser::{FileBrowser, Purpose};
use crate::structure_dialog::{Input, Outcome, StructureDialog};

/// What the file picker was opened for.
#[derive(Clone, Copy, PartialEq)]
enum Pick {
    Reference,
    Deformed,
    NyeCsv,
    VitekCsv,
    PatternCsv,
    DragStructure,
}

/// A background calculation and its error.
struct Task<T> {
    rx: Option<Receiver<T>>,
}

impl<T: Send + 'static> Task<T> {
    fn running(&self) -> bool {
        self.rx.is_some()
    }
    fn start(&mut self, f: impl FnOnce() -> T + Send + 'static) {
        let (tx, rx) = channel();
        std::thread::spawn(move || {
            let _ = tx.send(f());
        });
        self.rx = Some(rx);
    }
    fn poll(&mut self) -> Option<T> {
        let r = self.rx.as_ref()?.try_recv().ok()?;
        self.rx = None;
        Some(r)
    }
}

impl<T> Default for Task<T> {
    fn default() -> Self {
        Self { rx: None }
    }
}

pub struct DislocationAnalysisDialog {
    open: bool,
    picker: FileBrowser,
    pick: Pick,
    reference_path: String,
    reference_error: String,
    reference: Option<Structure>,
    deformed_path: String,
    deformed_error: String,
    deformed: Option<Structure>,
    status: String,

    nye_cutoff: f32,
    nye_no_pbc: bool,
    nye_task: Task<NyeTensorResult>,
    nye: NyeTensorResult,

    vitek_line: [f32; 3],
    vitek_burgers: f32,
    vitek_cutoff: f32,
    vitek_no_pbc: bool,
    vitek_task: Task<VitekResult>,
    vitek: VitekResult,

    pattern_cutoff: f32,
    pattern_no_pbc: bool,
    pattern_angle: f32,
    pattern_task: Task<(CrystalPattern, PatternDetectResult)>,
    pattern: CrystalPattern,
    pattern_result: PatternDetectResult,

    drag_zeta: f32,
    drag_clip: bool,
    drag_no_pbc: bool,
    drag_task: Task<DragPrepResult>,
    drag: DragPrepResult,

    fit_line: [f32; 3],
    fit_cutoff: f32,
    fit_area: f32,
    fit_no_pbc: bool,
    fit_task: Task<DislocationFitResult>,
    fit: DislocationFitResult,
}

impl Default for DislocationAnalysisDialog {
    fn default() -> Self {
        Self {
            open: false,
            picker: FileBrowser::default(),
            pick: Pick::Reference,
            reference_path: String::new(),
            reference_error: String::new(),
            reference: None,
            deformed_path: String::new(),
            deformed_error: String::new(),
            deformed: None,
            status: String::new(),
            nye_cutoff: 3.0,
            nye_no_pbc: false,
            nye_task: Task::default(),
            nye: NyeTensorResult::default(),
            vitek_line: [0.0, 0.0, 1.0],
            vitek_burgers: 0.0,
            vitek_cutoff: 0.0,
            vitek_no_pbc: false,
            vitek_task: Task::default(),
            vitek: VitekResult::default(),
            pattern_cutoff: 3.2,
            pattern_no_pbc: false,
            pattern_angle: 10.0,
            pattern_task: Task::default(),
            pattern: CrystalPattern::default(),
            pattern_result: PatternDetectResult::default(),
            drag_zeta: 0.5,
            drag_clip: false,
            drag_no_pbc: false,
            drag_task: Task::default(),
            drag: DragPrepResult::default(),
            fit_line: [0.0, 0.0, 1.0],
            fit_cutoff: 3.0,
            fit_area: 1.0,
            fit_no_pbc: false,
            fit_task: Task::default(),
            fit: DislocationFitResult::default(),
        }
    }
}

fn dvec(v: [f32; 3]) -> DVec3 {
    DVec3::new(v[0] as f64, v[1] as f64, v[2] as f64)
}

/// dialogLayout::primaryButton across the panel.
fn primary(ui: &mut egui::Ui, enabled: bool, text: &str) -> bool {
    ui.add_enabled_ui(enabled, |ui| crate::dialog_widgets::primary_button(ui, text, egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT))).inner.clicked()
}

/// ImGui::DragFloat at 55% of the panel, the label at its right.
fn drag(ui: &mut egui::Ui, label: &str, v: &mut f32, speed: f64, range: std::ops::RangeInclusive<f32>, decimals: usize, hint: &str) {
    let width = crate::analysis_panels::item_width(ui, 0.55);
    crate::dialog_widgets::drag_f32(ui, v, width, speed, range, decimals, label).on_hover_text(hint);
}

/// ImGui::DragFloat3 at 55% of the panel.
fn vector(ui: &mut egui::Ui, label: &str, v: &mut [f32; 3], hint: &str) {
    let width = crate::analysis_panels::item_width(ui, 0.55);
    ui.horizontal(|ui| {
        crate::dialog_widgets::field_style(ui);
        ui.spacing_mut().item_spacing.x = 4.0;
        let each = (width - 8.0) / 3.0;
        for c in v.iter_mut() {
            ui.add_sized([each, crate::theme::FRAME_HEIGHT], egui::DragValue::new(c).speed(0.01).fixed_decimals(3));
        }
        ui.spacing_mut().item_spacing.x = 8.0;
        ui.label(label);
    })
    .response
    .on_hover_text(hint);
}

const NO_PBC_HINT: &str = "Treat both structures as finite: neighbours are not found through periodic images.";

/// The Fit Dislocation calculation: the Nye tensor, then the weighted fit.
fn fit(reference: &Structure, deformed: &Structure, line: DVec3, cutoff: f32, area: f64, use_pbc: bool) -> DislocationFitResult {
    let nye = compute_nye_tensor(reference, deformed, cutoff, use_pbc);
    if !nye.success {
        return DislocationFitResult { message: format!("Nye tensor computation failed: {}", nye.message), ..Default::default() };
    }
    let positions: Vec<DVec3> = reference.atoms.iter().map(|a| DVec3::new(a.x, a.y, a.z)).collect();
    fit_dislocation_from_nye(&positions, &nye, line, area)
}

impl DislocationAnalysisDialog {
    pub fn load_reference(&mut self, path: &str) {
        self.reference_error.clear();
        match af_io::loader::load_structure_from_file(path) {
            Ok(s) => {
                self.reference = Some(s);
                self.reference_path = path.to_string();
            }
            Err(e) => {
                self.reference_error = e;
                self.reference = None;
                self.reference_path.clear();
            }
        }
    }

    pub fn load_deformed(&mut self, path: &str) {
        self.deformed_error.clear();
        match af_io::loader::load_structure_from_file(path) {
            Ok(s) => {
                self.deformed = Some(s);
                self.deformed_path = path.to_string();
            }
            Err(e) => {
                self.deformed_error = e;
                self.deformed = None;
                self.deformed_path.clear();
            }
        }
    }

    fn pair(&self) -> Option<(Structure, Structure)> {
        Some((self.reference.clone()?, self.deformed.clone()?))
    }

    pub fn nye_csv(&self) -> String {
        let mut out = String::from("index,symbol,alpha_xx,alpha_xy,alpha_xz,alpha_yx,alpha_yy,alpha_yz,alpha_zx,alpha_zy,alpha_zz,norm\n");
        let atoms = self.deformed.as_ref().map(|d| d.atoms.as_slice()).unwrap_or(&[]);
        for (i, alpha) in self.nye.alpha.iter().enumerate() {
            let _ = write!(out, "{i},{}", atoms.get(i).map(|a| a.symbol.as_str()).unwrap_or(""));
            for row in alpha {
                for v in row {
                    let _ = write!(out, ",{}", g(*v));
                }
            }
            let _ = writeln!(out, ",{}", g(self.nye.alpha_norm[i]));
        }
        out
    }

    pub fn vitek_csv(&self) -> String {
        let mut out = String::from("index_i,index_j,xi_A,yi_A,zi_A,xj_A,yj_A,zj_A,screw_A,edge_1_A,edge_2_A\n");
        for p in &self.vitek.pairs {
            let _ = writeln!(
                out,
                "{},{},{},{},{},{},{},{},{},{},{}",
                p.index_i,
                p.index_j,
                g(p.position_i.x),
                g(p.position_i.y),
                g(p.position_i.z),
                g(p.position_j.x),
                g(p.position_j.y),
                g(p.position_j.z),
                g(p.screw_component),
                g(p.edge_component.x),
                g(p.edge_component.y)
            );
        }
        out
    }

    pub fn pattern_csv(&self) -> String {
        let mut out = String::from("index,matched,neighbor_count,max_angle_deviation_deg\n");
        for r in &self.pattern_result.rows {
            let _ = writeln!(out, "{},{},{},{}", r.index, r.matched as i32, r.neighbor_count, g(r.max_angle_deviation_deg));
        }
        out
    }

    fn save(&mut self, pick: Pick, path: &str) {
        let result = match pick {
            Pick::NyeCsv => std::fs::write(path, self.nye_csv()).map_err(|e| e.to_string()),
            Pick::VitekCsv => std::fs::write(path, self.vitek_csv()).map_err(|e| e.to_string()),
            Pick::PatternCsv => std::fs::write(path, self.pattern_csv()).map_err(|e| e.to_string()),
            Pick::DragStructure => {
                let ext = std::path::Path::new(path).extension().map(|e| e.to_string_lossy().into_owned()).unwrap_or_default();
                let format = if ext.is_empty() { "vasp".to_string() } else { ext };
                af_io::loader::save_structure(&self.drag.interpolated, path, &format)
            }
            _ => Ok(()),
        };
        self.status = match result {
            Ok(()) => format!("Saved {}", file_name(path)),
            Err(e) => e,
        };
    }

    fn browse(&mut self, pick: Pick, title: &str, name: &str) {
        self.pick = pick;
        let purpose = if matches!(pick, Pick::Reference | Pick::Deformed) { Purpose::Open } else { Purpose::SaveAs };
        self.picker.show_for(purpose, title, name, &[]);
    }

    fn poll(&mut self) {
        if let Some(r) = self.nye_task.poll() {
            self.nye = r;
        }
        if let Some(r) = self.vitek_task.poll() {
            self.vitek = r;
        }
        if let Some((p, r)) = self.pattern_task.poll() {
            self.pattern = p;
            self.pattern_result = r;
        }
        if let Some(r) = self.drag_task.poll() {
            self.drag = r;
        }
        if let Some(r) = self.fit_task.poll() {
            self.fit = r;
        }
    }

    fn busy(&self) -> bool {
        self.nye_task.running() || self.vitek_task.running() || self.pattern_task.running() || self.drag_task.running() || self.fit_task.running()
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        section(ui, "Structures");
        ui.label("Both structures need the same atom count and ordering (e.g. a perfect crystal and the same crystal after inserting a dislocation).");
        if crate::dialog_widgets::wide_button(ui, "Open Reference (Undeformed)").clicked() {
            let name = file_name(&self.reference_path);
            self.browse(Pick::Reference, "Open reference structure", &name);
        }
        if let Some(r) = &self.reference {
            ui.label(format!("Reference: {} ({} atoms)", file_name(&self.reference_path), r.atoms.len()));
        }
        if !self.reference_error.is_empty() {
            ui.label(&self.reference_error);
        }
        if crate::dialog_widgets::wide_button(ui, "Open Deformed").clicked() {
            let name = file_name(&self.deformed_path);
            self.browse(Pick::Deformed, "Open deformed structure", &name);
        }
        if let Some(d) = &self.deformed {
            ui.label(format!("Deformed: {} ({} atoms)", file_name(&self.deformed_path), d.atoms.len()));
        }
        if !self.deformed_error.is_empty() {
            ui.label(&self.deformed_error);
        }
        let pair = self.pair();
        ui.add_enabled_ui(pair.is_some(), |ui| {
            section(ui, "Nye tensor");
            ui.label("Dislocation density tensor per atom (Hartley & Mishin lattice-correspondence method, matching BABEL's nyeTensor.f90).");
            drag(ui, "Cutoff (Å)", &mut self.nye_cutoff, 0.05, 0.1..=20.0, 2, "Neighbour cutoff (Å) for matching each atom's neighbours between the reference and the deformed structure; the first shell is usually enough.");
            crate::dialog_widgets::checkbox(ui, &mut self.nye_no_pbc, "Ignore periodicity").on_hover_text(NO_PBC_HINT);
            if primary(ui, !self.nye_task.running(), "Compute Nye tensor") {
                if let Some((r, d)) = pair.clone() {
                    let (cutoff, pbc) = (self.nye_cutoff, !self.nye_no_pbc);
                    self.nye_task.start(move || compute_nye_tensor(&r, &d, cutoff, pbc));
                }
            }
            if self.nye_task.running() {
                ui.label("Computing...");
            }

            section(ui, "Vitek differential-displacement map");
            ui.label("Displacement difference between neighbor pairs, split into screw (along line) and edge components (matching BABEL's vitek.f90).");
            vector(ui, "Line direction", &mut self.vitek_line, "Dislocation line direction (Cartesian). Displacement differences along it form the screw component, those across it the edge component.");
            drag(ui, "Burgers magnitude (Å)", &mut self.vitek_burgers, 0.01, 0.0..=20.0, 3, "Length of the Burgers vector (Å); the map's arrows are scaled relative to it.");
            drag(ui, "Cutoff (Å, 0=auto)", &mut self.vitek_cutoff, 0.05, 0.0..=20.0, 2, "Neighbour pairs closer than this (Å) are compared; 0 uses the first-neighbour distance.");
            crate::dialog_widgets::checkbox(ui, &mut self.vitek_no_pbc, "Ignore periodicity ").on_hover_text(NO_PBC_HINT);
            if primary(ui, !self.vitek_task.running(), "Compute Vitek map") {
                if let Some((r, d)) = pair.clone() {
                    let (line, b, cutoff, pbc) = (dvec(self.vitek_line), self.vitek_burgers as f64, self.vitek_cutoff, !self.vitek_no_pbc);
                    self.vitek_task.start(move || compute_vitek_map(&r, &d, line, b, cutoff, pbc));
                }
            }
            if self.vitek_task.running() {
                ui.label("Computing...");
            }

            section(ui, "Pattern matching");
            ui.label("Flags atoms whose local neighbor environment no longer matches the reference lattice (matching BABEL's patternInit/patternDetect).");
            drag(ui, "Cutoff (Å) ", &mut self.pattern_cutoff, 0.05, 0.1..=20.0, 2, "Neighbour cutoff (Å) for the local environments compared with the reference lattice.");
            drag(ui, "Angle threshold (°)", &mut self.pattern_angle, 0.2, 0.5..=90.0, 1, "Largest bond-direction deviation (degrees) that still matches the reference pattern; atoms beyond it are flagged.");
            crate::dialog_widgets::checkbox(ui, &mut self.pattern_no_pbc, "Ignore periodicity  ").on_hover_text(NO_PBC_HINT);
            if primary(ui, !self.pattern_task.running(), "Build Pattern + Detect") {
                if let Some((r, d)) = pair.clone() {
                    let (cutoff, angle, pbc) = (self.pattern_cutoff, self.pattern_angle as f64, !self.pattern_no_pbc);
                    self.pattern_result = PatternDetectResult::default();
                    self.pattern_task.start(move || {
                        let pattern = build_pattern(&r, cutoff, pbc);
                        let detected = if pattern.success { detect_pattern(&d, &pattern, cutoff, angle, pbc) } else { PatternDetectResult::default() };
                        (pattern, detected)
                    });
                }
            }
            if self.pattern_task.running() {
                ui.label("Building pattern and detecting...");
            }
            if !self.pattern.success && !self.pattern.message.is_empty() {
                ui.label(&self.pattern.message);
            }

            section(ui, "Prepare drag (migration barrier)");
            ui.label("Interpolates reference → deformed at reaction coordinate zeta for a constrained-minimization calculation (matching BABEL's prepareDrag). AtomForge does not itself run the constrained minimization.");
            ui.add(crate::dialog_widgets::Slider::new(&mut self.drag_zeta, 0.0..=1.0).fixed_decimals(3).text("Zeta")).on_hover_text("Reaction coordinate: 0 is the reference, 1 the deformed structure; positions in between are interpolated linearly.");
            crate::dialog_widgets::checkbox(ui, &mut self.drag_clip, "Clip displacement to minimum image").on_hover_text("Take each atom's displacement through the nearest periodic image, so atoms that crossed the cell boundary do not travel across the cell.");
            crate::dialog_widgets::checkbox(ui, &mut self.drag_no_pbc, "Ignore periodicity   ").on_hover_text(NO_PBC_HINT);
            if primary(ui, !self.drag_task.running(), "Interpolate") {
                if let Some((r, d)) = pair.clone() {
                    let (zeta, clip, pbc) = (self.drag_zeta as f64, self.drag_clip, !self.drag_no_pbc);
                    self.drag_task.start(move || prepare_drag(&r, &d, zeta, clip, pbc));
                }
            }
            if self.drag_task.running() {
                ui.label("Interpolating...");
            }
            if self.drag.success && crate::dialog_widgets::wide_button(ui, "Save Interpolated Structure").clicked() {
                self.browse(Pick::DragStructure, "Save interpolated structure", "interpolated.vasp");
            }

            section(ui, "Fit dislocation (position/Burgers vector)");
            ui.label("Recovers a dislocation's position and Burgers vector from the Nye tensor field (weighted centroid/integral), the same purpose as BABEL's displacementFit/vitekFit.");
            vector(ui, "Line direction ", &mut self.fit_line, "");
            drag(ui, "Cutoff (Å)  ", &mut self.fit_cutoff, 0.05, 0.1..=20.0, 2, "Neighbour cutoff (Å) for the Nye tensor the fit is computed from.");
            drag(ui, "Area per atom (Å²)", &mut self.fit_area, 0.05, 0.001..=1000.0, 3, "Cross-sectional area each atom represents perpendicular to the line (Å squared); it scales the integrated Burgers vector.");
            crate::dialog_widgets::checkbox(ui, &mut self.fit_no_pbc, "Ignore periodicity    ").on_hover_text(NO_PBC_HINT);
            if primary(ui, !self.fit_task.running(), "Fit Dislocation") {
                if let Some((r, d)) = pair.clone() {
                    let (line, cutoff, area, pbc) = (dvec(self.fit_line), self.fit_cutoff, self.fit_area as f64, !self.fit_no_pbc);
                    self.fit_task.start(move || fit(&r, &d, line, cutoff, area, pbc));
                }
            }
            if self.fit_task.running() {
                ui.label("Fitting (computes the Nye tensor internally)...");
            }
        });
        if pair.is_none() {
            ui.label("Load both a reference and a deformed structure to enable these analyses.");
        }
        if !self.status.is_empty() {
            ui.weak(&self.status);
        }
    }

    fn results(&mut self, ui: &mut egui::Ui) {
        if !self.nye.success && !self.vitek.success && !self.pattern_result.success && !self.drag.success && !self.fit.success {
            empty_state(
                ui,
                if self.reference.is_some() && self.deformed.is_some() {
                    "Compute any of the analyses on the left; the results appear here."
                } else {
                    "Open a reference (undeformed) and a deformed structure with the same atoms, then compute; the results appear here."
                },
            );
        }
        if !self.nye.success && !self.nye.message.is_empty() {
            ui.label(&self.nye.message);
        }
        if !self.vitek.success && !self.vitek.message.is_empty() {
            ui.label(&self.vitek.message);
        }
        if !self.pattern_result.success && !self.pattern_result.message.is_empty() {
            ui.label(&self.pattern_result.message);
        }
        if !self.drag.success && !self.drag.message.is_empty() {
            ui.label(&self.drag.message);
        }
        let full = egui::vec2(ui.available_width(), 0.0);
        if self.nye.success {
            section(ui, "Nye tensor result");
            ui.label(&self.nye.message);
            ui.label("Top atoms by |alpha| (dislocation content):");
            for i in top_atoms(&self.nye.alpha_norm, 10) {
                ui.label(format!("  atom {i}: |alpha| = {}", g4(self.nye.alpha_norm[i])));
            }
            if ui.add(egui::Button::new("Export Full Nye Tensor CSV").min_size(full)).clicked() {
                self.browse(Pick::NyeCsv, "Export Nye tensor", "nye_tensor.csv");
            }
        }
        if self.vitek.success {
            section(ui, "Vitek map result");
            ui.label(&self.vitek.message);
            if ui.add(egui::Button::new("Export Full Vitek Map CSV").min_size(full)).clicked() {
                self.browse(Pick::VitekCsv, "Export Vitek map", "vitek_map.csv");
            }
        }
        if self.pattern_result.success {
            section(ui, "Pattern matching result");
            ui.label(format!("Pattern: {}", self.pattern.message));
            ui.label(format!("Detection: {}", self.pattern_result.message));
            if ui.add(egui::Button::new("Export Full Pattern-match CSV").min_size(full)).clicked() {
                self.browse(Pick::PatternCsv, "Export pattern match", "pattern_match.csv");
            }
        }
        if self.drag.success {
            section(ui, "Drag preparation result");
            ui.label(&self.drag.message);
            ui.label(format!("Interpolated structure ready ({} atoms); use \"Save interpolated structure\" above.", self.drag.interpolated.atoms.len()));
        }
        if self.fit.success {
            section(ui, "Dislocation fit result");
            let (p, b) = (self.fit.line_position, self.fit.burgers_vector);
            ui.label(format!("Line position: ({:.4}, {:.4}, {:.4}) A", p.x, p.y, p.z));
            ui.label(format!("Burgers vector (uncalibrated by area per atom): ({:.4}, {:.4}, {:.4})", b.x, b.y, b.z));
            ui.label(format!("Magnitude: {:.4}", b.length()));
        } else if !self.fit.message.is_empty() {
            ui.label(&self.fit.message);
        }
    }
}

/// %.4g.
fn g4(v: f64) -> String {
    af_core::format::format_g(v, 4)
}

/// The indices of the `n` largest values, largest first (std::partial_sort with a stable tie order).
fn top_atoms(values: &[f64], n: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[b].total_cmp(&values[a]));
    order.truncate(n);
    order
}

fn flag(args: &str, name: &str) -> bool {
    args.split_whitespace().any(|w| w == name)
}

fn numbers(args: &str) -> Vec<f64> {
    args.split_whitespace().filter_map(|w| w.parse().ok()).collect()
}

impl StructureDialog for DislocationAnalysisDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn steps(&self) -> &'static [&'static str] {
        &["dislocation-analysis"]
    }

    fn title(&self) -> &'static str {
        "Dislocation Post-Processing"
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, _input: &Input) -> Option<Outcome> {
        self.poll();
        if self.busy() {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
        }
        if !self.open {
            return None;
        }
        let mut open = self.open;
        crate::theme::window("Dislocation Post-Processing").id(egui::Id::unique("dislocation post-processing")).open(&mut open).collapsible(false).default_size([1180.0, 860.0]).show(ctx, |ui| {
            two_panels(ui, "dislocation analysis", self, |s, ui| s.controls(ui), |s, ui| s.results(ui));
        });
        self.open = open;
        if let Some(path) = self.picker.show(ctx) {
            let path = path.to_string_lossy().into_owned();
            match self.pick {
                Pick::Reference => self.load_reference(&path),
                Pick::Deformed => self.load_deformed(&path),
                pick => self.save(pick, &path),
            }
        }
        None
    }

    /// reference PATH | deformed PATH | nye [CUTOFF] [nopbc] | vitek LX LY LZ B [CUTOFF] [nopbc] |
    /// pattern [CUTOFF] [ANGLE] [nopbc] | drag ZETA [clip] [nopbc] | fit LX LY LZ [CUTOFF] [AREA] [nopbc] |
    /// export nye|vitek|pattern|drag PATH | expect-nye-top ATOM | expect-nye-max VALUE |
    /// expect-vitek-pairs MIN [MAX] | expect-unmatched MIN [MAX] | expect-drag-atoms N |
    /// expect-fit-burgers MIN MAX | expect-fit-position X Y Z | expect-error TEXT
    fn script(&mut self, command: &str, args: &str, _input: &Input) -> Result<Option<Outcome>, String> {
        let args = args.trim();
        let n = numbers(args);
        let need = |s: &Self| s.pair().ok_or("load a reference and a deformed structure first".to_string());
        match command {
            "reference" => self.load_reference(args),
            "deformed" => self.load_deformed(args),
            "nye" => {
                let (r, d) = need(self)?;
                if let Some(&c) = n.first() {
                    self.nye_cutoff = c as f32;
                }
                self.nye_no_pbc = flag(args, "nopbc");
                self.nye = compute_nye_tensor(&r, &d, self.nye_cutoff, !self.nye_no_pbc);
                if !self.nye.success {
                    return Err(self.nye.message.clone());
                }
            }
            "vitek" => {
                let (r, d) = need(self)?;
                if n.len() < 4 {
                    return Err("vitek LX LY LZ BURGERS [CUTOFF] [nopbc]".into());
                }
                self.vitek_line = [n[0] as f32, n[1] as f32, n[2] as f32];
                self.vitek_burgers = n[3] as f32;
                self.vitek_cutoff = n.get(4).copied().unwrap_or(0.0) as f32;
                self.vitek_no_pbc = flag(args, "nopbc");
                self.vitek = compute_vitek_map(&r, &d, dvec(self.vitek_line), self.vitek_burgers as f64, self.vitek_cutoff, !self.vitek_no_pbc);
                if !self.vitek.success {
                    return Err(self.vitek.message.clone());
                }
            }
            "pattern" => {
                let (r, d) = need(self)?;
                if let Some(&c) = n.first() {
                    self.pattern_cutoff = c as f32;
                }
                if let Some(&a) = n.get(1) {
                    self.pattern_angle = a as f32;
                }
                self.pattern_no_pbc = flag(args, "nopbc");
                self.pattern = build_pattern(&r, self.pattern_cutoff, !self.pattern_no_pbc);
                if !self.pattern.success {
                    return Err(self.pattern.message.clone());
                }
                self.pattern_result = detect_pattern(&d, &self.pattern, self.pattern_cutoff, self.pattern_angle as f64, !self.pattern_no_pbc);
                if !self.pattern_result.success {
                    return Err(self.pattern_result.message.clone());
                }
            }
            "drag" => {
                let (r, d) = need(self)?;
                self.drag_zeta = n.first().copied().ok_or("drag ZETA [clip] [nopbc]")? as f32;
                self.drag_clip = flag(args, "clip");
                self.drag_no_pbc = flag(args, "nopbc");
                self.drag = prepare_drag(&r, &d, self.drag_zeta as f64, self.drag_clip, !self.drag_no_pbc);
                if !self.drag.success {
                    return Err(self.drag.message.clone());
                }
            }
            "fit" => {
                let (r, d) = need(self)?;
                if n.len() < 3 {
                    return Err("fit LX LY LZ [CUTOFF] [AREA] [nopbc]".into());
                }
                self.fit_line = [n[0] as f32, n[1] as f32, n[2] as f32];
                if let Some(&c) = n.get(3) {
                    self.fit_cutoff = c as f32;
                }
                if let Some(&a) = n.get(4) {
                    self.fit_area = a as f32;
                }
                self.fit_no_pbc = flag(args, "nopbc");
                self.fit = fit(&r, &d, dvec(self.fit_line), self.fit_cutoff, self.fit_area as f64, !self.fit_no_pbc);
                if !self.fit.success {
                    return Err(self.fit.message.clone());
                }
            }
            "export" => {
                let (what, path) = args.split_once(' ').ok_or("export nye|vitek|pattern|drag PATH")?;
                let pick = match what {
                    "nye" => Pick::NyeCsv,
                    "vitek" => Pick::VitekCsv,
                    "pattern" => Pick::PatternCsv,
                    "drag" => Pick::DragStructure,
                    _ => return Err(format!("nothing to export called {what}")),
                };
                self.save(pick, path.trim());
                if !self.status.starts_with("Saved") {
                    return Err(self.status.clone());
                }
            }
            "expect-nye-top" => {
                let want = n.first().copied().ok_or("expect-nye-top ATOM")? as usize;
                let top = top_atoms(&self.nye.alpha_norm, 1);
                if top.first() != Some(&want) {
                    return Err(format!("the largest |alpha| is at atom {top:?}"));
                }
            }
            "expect-nye-max" => {
                let want = n.first().copied().ok_or("expect-nye-max VALUE")?;
                let max = self.nye.alpha_norm.iter().copied().fold(0.0f64, f64::max);
                if (max - want).abs() > 1e-6 * want.abs().max(1.0) {
                    return Err(format!("max |alpha| is {max}"));
                }
            }
            "expect-vitek-pairs" => {
                let (min, max) = (n.first().copied().unwrap_or(1.0) as usize, n.get(1).map(|&v| v as usize).unwrap_or(usize::MAX));
                if self.vitek.pairs.len() < min || self.vitek.pairs.len() > max {
                    return Err(format!("{} Vitek pairs", self.vitek.pairs.len()));
                }
            }
            "expect-unmatched" => {
                // expect-unmatched MIN [MAX]
                let count = self.pattern_result.rows.iter().filter(|r| !r.matched).count();
                let (min, max) = (n.first().copied().unwrap_or(1.0) as usize, n.get(1).map(|&v| v as usize).unwrap_or(usize::MAX));
                if count < min || count > max {
                    return Err(format!("{count} atoms do not match the pattern"));
                }
            }
            "expect-drag-atoms" => {
                let want = n.first().copied().unwrap_or(0.0) as usize;
                if !self.drag.success || self.drag.interpolated.atoms.len() != want {
                    return Err(format!("interpolated structure has {} atoms ({})", self.drag.interpolated.atoms.len(), self.drag.message));
                }
            }
            "expect-fit-burgers" => {
                let mag = self.fit.burgers_vector.length();
                if n.len() != 2 || !self.fit.success || mag < n[0] || mag > n[1] {
                    return Err(format!("fitted |b| = {mag} ({})", self.fit.message));
                }
            }
            "expect-fit-position" => {
                let p = self.fit.line_position;
                if n.len() != 3 || !self.fit.success || (p - DVec3::new(n[0], n[1], n[2])).length() > 1e-6 {
                    return Err(format!("fitted line position {p} ({})", self.fit.message));
                }
            }
            "expect-error" => {
                let error = format!("{}{}", self.reference_error, self.deformed_error);
                if !error.contains(args) || error.is_empty() {
                    return Err(format!("error is {error:?}"));
                }
            }
            other => return Err(format!("unknown dislocation-analysis command {other}")),
        }
        Ok(None)
    }
}
