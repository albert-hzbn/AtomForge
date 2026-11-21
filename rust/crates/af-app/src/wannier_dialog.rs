//! Analysis > Electronic Structure > Wannier Analysis (the C++
//! ui/WannierAnalysisDialog): bands along a piecewise k-path, the Kubo Berry
//! curvature at a k-point and Fukui-Hatsugai-Suzuki Chern numbers from a
//! Wannier90 seedname_hr.dat tight-binding Hamiltonian. The calculations run
//! in the background (synchronously from UI scripts).

use std::sync::mpsc::{channel, Receiver};

use af_electronic::wannier::{berry_curvature, chern_number, read_wannier_hr, wannier_bands, WannierHamiltonian};
use eframe::egui;
use glam::DVec3;

use crate::analysis_panels::{empty_state, file_name, line_plot, section, two_panels};
use crate::file_browser::{FileBrowser, Purpose};
use crate::structure_dialog::{Input, Outcome, StructureDialog};

const PLANES: [&str; 3] = ["kx-ky", "kx-kz", "ky-kz"];
const PLANE_A: [i32; 3] = [0, 0, 1];
const PLANE_B: [i32; 3] = [1, 2, 2];

#[derive(Clone, Default)]
pub struct BandsResult {
    /// [k][band].
    pub energies: Vec<f64>,
    pub num_wann: i32,
    pub num_k: i32,
}

/// One k-point (fractional reciprocal) per line; at least two, interpolated
/// linearly with `samples` points per segment plus the final anchor.
pub fn interpolate_path(text: &str, samples: i32) -> Result<Vec<DVec3>, String> {
    let mut anchors = Vec::new();
    for line in text.lines() {
        let v: Vec<f64> = line.split_whitespace().map_while(|w| w.parse().ok()).take(3).collect();
        if v.len() == 3 {
            anchors.push(DVec3::new(v[0], v[1], v[2]));
        }
    }
    if anchors.len() < 2 {
        return Err("Supply at least two k-points (one per line, fractional reciprocal)".into());
    }
    let mut path = Vec::new();
    for w in anchors.windows(2) {
        for s in 0..samples {
            let t = s as f64 / samples as f64;
            path.push(w[0] * (1.0 - t) + w[1] * t);
        }
    }
    path.push(*anchors.last().unwrap());
    Ok(path)
}

pub fn compute_bands(model: &WannierHamiltonian, text: &str, samples: i32) -> Result<BandsResult, String> {
    let path = interpolate_path(text, samples)?;
    let energies = wannier_bands(model, &path)?;
    Ok(BandsResult { energies, num_wann: model.num_wann, num_k: path.len() as i32 })
}

/// A background calculation: its result arrives on the channel.
struct Task<T> {
    rx: Option<Receiver<Result<T, String>>>,
    error: String,
}

impl<T: Send + 'static> Task<T> {
    fn new() -> Self {
        Self { rx: None, error: String::new() }
    }
    fn running(&self) -> bool {
        self.rx.is_some()
    }
    fn start(&mut self, f: impl FnOnce() -> Result<T, String> + Send + 'static) {
        let (tx, rx) = channel();
        std::thread::spawn(move || {
            let _ = tx.send(f());
        });
        self.rx = Some(rx);
        self.error.clear();
    }
    /// The finished result, once.
    fn poll(&mut self) -> Option<T> {
        let r = self.rx.as_ref()?.try_recv().ok()?;
        self.rx = None;
        match r {
            Ok(v) => Some(v),
            Err(e) => {
                self.error = e;
                None
            }
        }
    }
}

pub struct WannierDialog {
    open: bool,
    picker: FileBrowser,
    loaded_path: String,
    load_error: String,
    model: WannierHamiltonian,
    kpath: String,
    samples_per_segment: i32,
    bands_task: Task<BandsResult>,
    bands: BandsResult,
    curvature_k: [f32; 3],
    curvature_plane: usize,
    curvature_task: Task<Vec<f64>>,
    curvature: Vec<f64>,
    chern_band: i32,
    chern_plane: usize,
    chern_grid: i32,
    chern_fixed: f32,
    chern_task: Task<f64>,
    chern: Option<f64>,
}

impl Default for WannierDialog {
    fn default() -> Self {
        Self {
            open: false,
            picker: FileBrowser::default(),
            loaded_path: String::new(),
            load_error: String::new(),
            model: WannierHamiltonian::default(),
            kpath: "0 0 0\n0.5 0 0\n0.5 0.5 0\n0 0 0".into(),
            samples_per_segment: 40,
            bands_task: Task::new(),
            bands: BandsResult::default(),
            curvature_k: [0.0; 3],
            curvature_plane: 0,
            curvature_task: Task::new(),
            curvature: Vec::new(),
            chern_band: 0,
            chern_plane: 0,
            chern_grid: 30,
            chern_fixed: 0.0,
            chern_task: Task::new(),
            chern: None,
        }
    }
}

fn plane_combo(ui: &mut egui::Ui, label: &str, plane: &mut usize) {
    let width = crate::analysis_panels::item_width(ui, 0.55);
    egui::ComboBox::new(label, label.trim()).width(width).selected_text(PLANES[*plane]).show_ui(ui, |ui| {
        for (i, name) in PLANES.iter().enumerate() {
            ui.selectable_value(plane, i, *name);
        }
    });
}

impl WannierDialog {
    pub fn load(&mut self, path: &str) {
        self.load_error.clear();
        match read_wannier_hr(path) {
            Ok(m) => {
                self.model = m;
                self.loaded_path = path.to_string();
            }
            Err(e) => {
                self.load_error = e;
                self.model = WannierHamiltonian::default();
                self.loaded_path.clear();
            }
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        if ui.button("Open seedname_hr.dat").clicked() {
            let name = if self.loaded_path.is_empty() { "wannier90_hr.dat".to_string() } else { file_name(&self.loaded_path) };
            self.picker.show_for(Purpose::Open, "Open Wannier90 Hamiltonian", &name, &[]);
        }
        if !self.loaded_path.is_empty() {
            ui.label(format!("Loaded: {} ({} Wannier functions)", file_name(&self.loaded_path), self.model.num_wann));
        }
        if !self.load_error.is_empty() {
            ui.label(&self.load_error);
        }
        ui.add(egui::Label::new("H(k) = sum_R exp(2*pi*i*k.R) H(R) / ndegen(R), k fractional reciprocal.").wrap());
        ui.add_enabled_ui(self.model.num_wann > 0, |ui| {
            use crate::dialog_widgets as w;
            section(ui, "Bands");
            ui.add(egui::Label::new("One fractional k-point (kx ky kz) per line; at least two, forming a piecewise path.").wrap());
            ui.add_sized([ui.available_width(), 90.0], egui::TextEdit::multiline(&mut self.kpath));
            let width = crate::analysis_panels::item_width(ui, 0.55);
            let mut samples = self.samples_per_segment as i64;
            if w::input_int(ui, &mut samples, width, i32::MIN as i64..=i32::MAX as i64, "Samples per segment").changed() {
                self.samples_per_segment = samples as i32;
            }
            if ui.add_enabled_ui(!self.bands_task.running(), |ui| w::primary_button(ui, "Compute Bands", egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT))).inner.clicked() {
                let (text, samples, model) = (self.kpath.clone(), self.samples_per_segment.clamp(1, 2000), self.model.clone());
                self.bands_task.start(move || compute_bands(&model, &text, samples));
            }
            if self.bands_task.running() {
                ui.label("Computing...");
            }
            if !self.bands_task.error.is_empty() {
                ui.add(egui::Label::new(&self.bands_task.error).wrap());
            }

            section(ui, "Berry curvature (Hamiltonian-gauge/Kubo term)");
            let width = crate::analysis_panels::item_width(ui, 0.55);
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 4.0;
                let each = (width - 8.0) / 3.0;
                w::field_style(ui);
                for c in &mut self.curvature_k {
                    ui.add_sized([each, crate::theme::FRAME_HEIGHT], egui::DragValue::new(c).speed(0.0).fixed_decimals(3));
                }
                ui.spacing_mut().item_spacing.x = 8.0;
                ui.label("k-point");
            });
            plane_combo(ui, "Plane", &mut self.curvature_plane);
            if ui.add_enabled_ui(!self.curvature_task.running(), |ui| w::primary_button(ui, "Compute Curvature", egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT))).inner.clicked() {
                let k = DVec3::new(self.curvature_k[0] as f64, self.curvature_k[1] as f64, self.curvature_k[2] as f64);
                let (a, b, model) = (PLANE_A[self.curvature_plane], PLANE_B[self.curvature_plane], self.model.clone());
                self.curvature_task.start(move || berry_curvature(&model, k, a, b));
            }
            if !self.curvature_task.error.is_empty() {
                ui.add(egui::Label::new(&self.curvature_task.error).wrap());
            }
            ui.add(egui::Label::new("Diverges at exact band degeneracies; use Chern number below for those k-points instead.").wrap());

            section(ui, "Chern number (Fukui-Hatsugai-Suzuki)");
            let width = crate::analysis_panels::item_width(ui, 0.55);
            let mut band = self.chern_band as i64;
            if w::input_int(ui, &mut band, width, i32::MIN as i64..=i32::MAX as i64, "Band (0-based)").changed() {
                self.chern_band = band as i32;
            }
            plane_combo(ui, "Plane ", &mut self.chern_plane);
            w::slider_i32(ui, &mut self.chern_grid, 4..=120, width, "Mesh resolution");
            w::drag_f32(ui, &mut self.chern_fixed, width, 0.0, f32::MIN..=f32::MAX, 3, "Fixed fractional coordinate");
            if ui.add_enabled_ui(!self.chern_task.running(), |ui| w::primary_button(ui, "Compute Chern number", egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT))).inner.clicked() {
                let (a, b) = (PLANE_A[self.chern_plane], PLANE_B[self.chern_plane]);
                let (band, grid, fixed, model) = (self.chern_band, self.chern_grid.clamp(2, 2000), self.chern_fixed as f64, self.model.clone());
                self.chern = None;
                self.chern_task.start(move || chern_number(&model, band, a, b, grid, fixed));
            }
            if !self.chern_task.error.is_empty() {
                ui.add(egui::Label::new(&self.chern_task.error).wrap());
            }
        });
    }

    fn results(&self, ui: &mut egui::Ui) {
        if self.bands.num_wann <= 0 && self.curvature.is_empty() && self.chern.is_none() {
            empty_state(
                ui,
                if self.loaded_path.is_empty() {
                    "Open a Wannier90 seedname_hr.dat, then compute bands, Berry curvature or a Chern number; the results appear here."
                } else {
                    "Compute bands, Berry curvature or a Chern number on the left; the results appear here."
                },
            );
        }
        if self.bands.num_wann > 0 {
            let (nw, nk) = (self.bands.num_wann as usize, self.bands.num_k as usize);
            ui.label(format!("{nw} bands, {nk} k-points along the path (index axis; no reciprocal-lattice distance is stored in seedname_hr.dat)."));
            let x: Vec<f64> = (0..nk).map(|k| k as f64).collect();
            let series = (0..nw).map(|b| (format!("Band {b}"), (0..nk).map(|k| self.bands.energies[k * nw + b]).collect())).collect();
            line_plot(ui, "Bands", "k-point index", "eV", &x, series, 260.0);
        }
        if !self.curvature.is_empty() {
            section(ui, "Berry curvature per band");
            for (i, v) in self.curvature.iter().enumerate() {
                ui.label(format!("Band {i}: {}", af_core::format::format_g(*v, 6)));
            }
        }
        if let Some(c) = self.chern {
            section(ui, "Chern number result");
            ui.label(format!("C = {c:.6}"));
            if (c - c.round()).abs() > 0.05 {
                ui.label("Far from an integer: refine the mesh resolution or check the band/plane for a degeneracy.");
            }
        }
    }

    fn poll(&mut self) {
        if let Some(b) = self.bands_task.poll() {
            self.bands = b;
        }
        if let Some(c) = self.curvature_task.poll() {
            self.curvature = c;
        }
        if let Some(c) = self.chern_task.poll() {
            self.chern = Some(c);
        }
    }
}

fn numbers(args: &str) -> Vec<f64> {
    args.split_whitespace().filter_map(|w| w.parse().ok()).collect()
}

fn plane_index(name: &str) -> Result<usize, String> {
    PLANES.iter().position(|p| *p == name.trim()).ok_or(format!("plane is one of {PLANES:?}"))
}

impl StructureDialog for WannierDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn steps(&self) -> &'static [&'static str] {
        &["wannier"]
    }

    fn title(&self) -> &'static str {
        "Wannier Analysis"
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, _input: &Input) -> Option<Outcome> {
        self.poll();
        if self.bands_task.running() || self.curvature_task.running() || self.chern_task.running() {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
        }
        if !self.open {
            return None;
        }
        let mut open = self.open;
        crate::theme::window("Wannier Analysis").id(egui::Id::unique("wannier analysis")).open(&mut open).collapsible(false).default_size([1180.0, 820.0]).show(ctx, |ui| {
            two_panels(ui, "wannier", self, |s, ui| s.controls(ui), |s, ui| s.results(ui));
        });
        self.open = open;
        if let Some(path) = self.picker.show(ctx) {
            self.load(&path.to_string_lossy());
        }
        None
    }

    /// load PATH | kpath K1; K2; ... | samples N | bands | curvature KX KY KZ [PLANE] |
    /// chern BAND [PLANE] [GRID] [FIXED] | expect-wann N | expect-bands NBANDS NK |
    /// expect-band-range BAND LOW HIGH | expect-chern C | expect-curvature BAND VALUE | expect-error TEXT
    fn script(&mut self, command: &str, args: &str, _input: &Input) -> Result<Option<Outcome>, String> {
        let args = args.trim();
        let loaded = |s: &Self| if s.model.num_wann > 0 { Ok(()) } else { Err("load a seedname_hr.dat first".to_string()) };
        match command {
            "load" => self.load(args),
            "kpath" => self.kpath = args.split(';').map(str::trim).collect::<Vec<_>>().join("\n"),
            "samples" => self.samples_per_segment = args.parse().map_err(|_| "samples N")?,
            "bands" => {
                loaded(self)?;
                self.bands = compute_bands(&self.model, &self.kpath, self.samples_per_segment.clamp(1, 2000))?;
            }
            "curvature" => {
                loaded(self)?;
                let w: Vec<&str> = args.split_whitespace().collect();
                let n = numbers(&w[..w.len().min(3)].join(" "));
                if n.len() != 3 {
                    return Err("curvature KX KY KZ [PLANE]".into());
                }
                self.curvature_k = [n[0] as f32, n[1] as f32, n[2] as f32];
                if let Some(p) = w.get(3) {
                    self.curvature_plane = plane_index(p)?;
                }
                let k = DVec3::new(self.curvature_k[0] as f64, self.curvature_k[1] as f64, self.curvature_k[2] as f64);
                self.curvature = berry_curvature(&self.model, k, PLANE_A[self.curvature_plane], PLANE_B[self.curvature_plane])?;
            }
            "chern" => {
                loaded(self)?;
                let w: Vec<&str> = args.split_whitespace().collect();
                self.chern_band = w.first().and_then(|v| v.parse().ok()).ok_or("chern BAND [PLANE] [GRID] [FIXED]")?;
                if let Some(p) = w.get(1) {
                    self.chern_plane = plane_index(p)?;
                }
                if let Some(g) = w.get(2) {
                    self.chern_grid = g.parse().map_err(|_| "grid is an integer")?;
                }
                if let Some(f) = w.get(3) {
                    self.chern_fixed = f.parse().map_err(|_| "fixed is a number")?;
                }
                let (a, b) = (PLANE_A[self.chern_plane], PLANE_B[self.chern_plane]);
                self.chern = Some(chern_number(&self.model, self.chern_band, a, b, self.chern_grid.clamp(2, 2000), self.chern_fixed as f64)?);
            }
            "expect-wann" => {
                let n: i32 = args.parse().map_err(|_| "expect-wann N")?;
                if self.model.num_wann != n {
                    return Err(format!("{} Wannier functions ({})", self.model.num_wann, self.load_error));
                }
            }
            "expect-bands" => {
                let n = numbers(args);
                if n.len() != 2 || self.bands.num_wann != n[0] as i32 || self.bands.num_k != n[1] as i32 {
                    return Err(format!("{} bands at {} k-points", self.bands.num_wann, self.bands.num_k));
                }
            }
            "expect-band-range" => {
                let n = numbers(args);
                let nw = self.bands.num_wann as usize;
                if n.len() != 3 || n[0] as usize >= nw {
                    return Err("expect-band-range BAND LOW HIGH".into());
                }
                let values = self.bands.energies.iter().skip(n[0] as usize).step_by(nw);
                let (lo, hi) = values.fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
                if (lo - n[1]).abs() > 1e-4 || (hi - n[2]).abs() > 1e-4 {
                    return Err(format!("band {} spans {lo} to {hi}", n[0]));
                }
            }
            "expect-chern" => {
                let want: f64 = args.parse().map_err(|_| "expect-chern C")?;
                match self.chern {
                    Some(c) if (c - want).abs() < 1e-3 => {}
                    other => return Err(format!("Chern number {other:?}")),
                }
            }
            "expect-curvature" => {
                let n = numbers(args);
                let got = n.first().and_then(|&b| self.curvature.get(b as usize)).copied();
                match (got, n.get(1)) {
                    (Some(v), Some(&want)) if (v - want).abs() <= 1e-6 * want.abs().max(1.0) => {}
                    _ => return Err(format!("curvature {:?}", self.curvature)),
                }
            }
            "expect-error" => {
                if self.load_error.is_empty() || !self.load_error.contains(args) {
                    return Err(format!("load error is {:?}", self.load_error));
                }
            }
            other => return Err(format!("unknown wannier command {other}")),
        }
        Ok(None)
    }
}
