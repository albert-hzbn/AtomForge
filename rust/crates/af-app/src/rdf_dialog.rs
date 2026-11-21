//! Analysis > Structure and Defects > Radial Distribution Function (the C++
//! ui/RadialDistributionAnalysis): species pair, periodicity, g(r) or
//! counts, radius range, bins and smoothing, raw-count and cumulative
//! overlays, and first-shell distortion metrics; computed in the background.

use std::collections::BTreeSet;
use std::sync::mpsc::{channel, Receiver};

use af_algorithms::rdf::{compute_rdf, RdfParams, RdfResult};
use af_core::elements;
use af_core::structure::Structure;
use eframe::egui;

pub struct RdfDialog {
    pub open: bool,
    pub params: RdfParams,
    reference: i32,
    target: i32,
    show_raw: bool,
    show_cumulative: bool,
    running: Option<Receiver<RdfResult>>,
    pub result: Option<RdfResult>,
}

impl Default for RdfDialog {
    fn default() -> Self {
        Self { open: false, params: RdfParams::default(), reference: 0, target: 0, show_raw: false, show_cumulative: false, running: None, result: None }
    }
}

/// "All elements" and each element of the structure (atomic number, label).
fn species(structure: &Structure) -> Vec<(i32, String)> {
    let mut options = vec![(0, "All elements".to_string())];
    let mut seen = BTreeSet::new();
    for a in &structure.atoms {
        if seen.insert(a.atomic_number) {
            options.push((a.atomic_number, format!("{} ({})", elements::symbol(a.atomic_number), a.atomic_number)));
        }
    }
    options
}

impl RdfDialog {
    /// Sets the species pair by element symbol ("all" for every element; UI scripts).
    pub fn set_pair(&mut self, reference: &str, target: &str) -> Result<(), String> {
        let z = |s: &str| if s.eq_ignore_ascii_case("all") { Ok(0) } else { Some(elements::atomic_number(s)).filter(|&z| z > 0).ok_or(format!("no element {s}")) };
        self.reference = z(reference)?;
        self.target = z(target)?;
        Ok(())
    }

    /// Runs the analysis on this thread (UI scripts).
    pub fn run_now(&mut self, structure: &Structure) -> &RdfResult {
        let mut p = self.params.clone();
        p.ref_atomic_number_filter = self.reference;
        p.target_atomic_number_filter = self.target;
        self.result = Some(compute_rdf(structure, &p));
        self.result.as_ref().unwrap()
    }

    pub fn show(&mut self, ctx: &egui::Context, structure: &Structure) {
        if let Some(rx) = &self.running {
            if let Ok(r) = rx.try_recv() {
                self.result = Some(r);
                self.running = None;
            } else {
                ctx.request_repaint_after(std::time::Duration::from_millis(100));
            }
        }
        if !self.open {
            return;
        }
        let options = species(structure);
        if !options.iter().any(|o| o.0 == self.reference) {
            self.reference = 0;
        }
        if !options.iter().any(|o| o.0 == self.target) {
            self.target = 0;
        }
        let mut open = self.open;
        crate::theme::modal_window("Radial Distribution Function").open(&mut open).fixed_size([1180.0, 820.0]).show(ctx, |ui| {
            // Controls (43 %) and results (57 %) in bordered panels above the footer.
            let footer = crate::theme::FRAME_HEIGHT + 3.0 * ui.spacing().item_spacing.y + 4.0;
            let height = (ui.available_height() - footer).max(300.0);
            let width = ui.available_width();
            let gap = ui.spacing().item_spacing.x;
            let left = (width - gap) * 0.43;
            ui.horizontal_top(|ui| {
                panel(ui, egui::vec2(left, height), "rdf controls", |ui| self.controls(ui, structure, &options));
                panel(ui, egui::vec2(width - left - gap, height), "rdf results", |ui| self.results(ui));
            });
            if crate::dialog_widgets::footer(ui, "Close", None, true) == 1 {
                self.open = false;
            }
        });
        self.open = open && self.open;
    }

    fn controls(&mut self, ui: &mut egui::Ui, structure: &Structure, options: &[(i32, String)]) {
        use crate::dialog_widgets as w;
        let p = &mut self.params;
        let iw = (ui.available_width() * 0.65).max(80.0);
        w::separator_text(ui, "Scope and Pair Selection");
        for (label, value, hint) in [
            ("Reference species", &mut self.reference, "Element of the centre atoms whose neighbours are counted (all elements: every atom)."),
            ("Target species", &mut self.target, "Element of the neighbours counted around the centre atoms; a pair of elements gives a partial g(r)."),
        ] {
            let text = options.iter().find(|o| o.0 == *value).map(|o| o.1.clone()).unwrap_or_default();
            egui::ComboBox::from_label(label).selected_text(text).width(iw).show_ui(ui, |ui| {
                for (z, name) in options {
                    ui.selectable_value(value, *z, name);
                }
            })
            .response
            .on_hover_text(hint);
        }
        w::checkbox(ui, &mut p.use_pbc, "Use PBC when unit cell is available").on_hover_text("Count neighbours through the periodic images of the cell (minimum image) when the structure has one.");
        w::checkbox(ui, &mut p.normalize, "Normalize to g(r)").on_hover_text("Divide the pair counts by those of an ideal gas of the same density, so g(r) tends to 1 at long range. Off: pair counts per shell.");
        w::separator_text(ui, "Sampling");
        range2(ui, iw, &mut p.r_min, &mut p.r_max, 0.01, 2, "Radius range");
        w::slider_i32(ui, &mut p.bin_count, 32..=2000, iw, "Bins").on_hover_text("Number of histogram bins between the minimum and maximum radius.");
        w::slider_i32(ui, &mut p.smoothing_passes, 0..=8, iw, "Smoothing passes").on_hover_text("Moving-average passes applied to the displayed curve. The raw histogram is not smoothed.");
        w::separator_text(ui, "Plot Overlays");
        w::checkbox(ui, &mut self.show_raw, "Overlay raw counts").on_hover_text("Also plot the unnormalised pair counts per bin.");
        w::checkbox(ui, &mut self.show_cumulative, "Overlay cumulative coordination").on_hover_text("Also plot the running coordination number: the average number of target atoms within r of a centre atom.");
        w::separator_text(ui, "Distortion Metrics (Alloy/Compound)");
        w::checkbox(ui, &mut p.enable_distortion_analysis, "Enable distortion analysis").on_hover_text("Report bond-length statistics of the first shell (mean, spread, coordination) for alloys and compounds.");
        if p.enable_distortion_analysis {
            w::checkbox(ui, &mut p.auto_distortion_window, "Auto shell window from first-shell minima").on_hover_text("Take the first-shell distance window from the g(r) minima around the first peak instead of fixed limits.");
            if !p.auto_distortion_window {
                range2(ui, iw, &mut p.manual_distortion_min, &mut p.manual_distortion_max, 0.005, 3, "Manual shell window");
            }
            ui.add(egui::Label::new("Use pair-resolved RDF (for example A-B) to track bond-length disorder and local distortion in multicomponent systems.").wrap());
        }
        ui.add_space(4.0);
        let running = self.running.is_some();
        let width = ui.available_width();
        if ui.add_enabled_ui(!running && !structure.atoms.is_empty(), |ui| w::primary_button(ui, "Run RDF", egui::vec2(width, crate::theme::FRAME_HEIGHT))).inner.clicked() {
            let mut params = self.params.clone();
            params.ref_atomic_number_filter = self.reference;
            params.target_atomic_number_filter = self.target;
            let snapshot = structure.clone();
            let (tx, rx) = channel();
            std::thread::spawn(move || {
                let _ = tx.send(compute_rdf(&snapshot, &params));
            });
            self.running = Some(rx);
        }
        if running {
            ui.label(egui::RichText::new("Computing...").color(ui.visuals().warn_fg_color));
        } else if self.result.is_some() {
            ui.label(egui::RichText::new("Done").color(w::note_colour(ui, w::Note::Success)));
        }
    }

    fn results(&mut self, ui: &mut egui::Ui) {
        // Empty until the first run, as the C++ panel.
        let Some(r) = &self.result else { return };
        ui.strong("RDF Summary");
        ui.label(format!("Status: {}", r.message));
        ui.label(format!("Atoms: {}    Ref: {}    Target: {}", r.atom_count, r.ref_count, r.target_count));
        ui.label(format!("PBC used: {}", if r.pbc_used { "Yes" } else { "No" }));
        ui.label(format!("Volume: {:.6} Å³    Target density: {:.6} Å⁻³", r.volume, r.density));
        ui.label(format!("Bin width: {:.5} Å", r.bin_width));
        ui.label(if r.has_first_peak { format!("First peak: r = {:.4} Å, value = {:.4}", r.first_peak_r, r.first_peak_value) } else { "First peak: not detected".into() });
        ui.label(if r.has_first_minimum { format!("First minimum after peak: r = {:.4} Å, value = {:.4}", r.first_minimum_r, r.first_minimum_value) } else { "First minimum after peak: not detected".into() });
        if r.has_first_peak_fwhm {
            ui.label(format!("First-peak FWHM: {:.4} Å", r.first_peak_fwhm));
        }
        draw_plot(ui, r, self.show_raw, self.show_cumulative);
        ui.strong("Alloy/Compound Distortion");
        if !r.valid {
            ui.weak("Run RDF to compute distortion metrics.");
        } else if !r.has_distortion_metrics {
            ui.label("No valid shell window was found for distortion analysis. Use manual shell bounds or adjust radius/bin settings.");
        } else {
            ui.label(format!("Shell window: [{:.4}, {:.4}] Å ({})", r.distortion_window_min, r.distortion_window_max, if r.distortion_window_auto { "auto" } else { "manual" }));
            ui.label(format!("Shell mean bond length: {:.5} Å", r.shell_mean_distance));
            ui.label(format!("Shell sigma (bond spread): {:.5} Å", r.shell_std_dev));
            ui.label(format!("Relative distortion sigma/r_mean: {:.3} %", r.shell_relative_distortion));
            ui.label(format!("Shell coordination (pair-resolved): {:.4}", r.shell_coordination));
        }
    }
}

/// g(r) with optional raw counts and cumulative coordination (drawPlot).
fn draw_plot(ui: &mut egui::Ui, r: &RdfResult, show_raw: bool, show_cumulative: bool) {
    let light = !ui.visuals().dark_mode;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().max(50.0), 260.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let rgba = egui::Color32::from_rgba_unmultiplied;
    painter.rect_filled(rect, 0.0, if light { rgba(245, 246, 248, 255) } else { rgba(20, 24, 31, 255) });
    painter.rect_stroke(rect, 0.0, egui::Stroke::new(1.0, if light { rgba(140, 145, 160, 255) } else { rgba(110, 125, 150, 255) }), egui::StrokeKind::Inside);
    let font = egui::FontId::proportional(12.0);
    let label = if light { rgba(50, 50, 65, 255) } else { rgba(200, 200, 200, 255) };
    if !r.valid || r.bins.is_empty() {
        painter.text(rect.min + egui::vec2(10.0, 10.0), egui::Align2::LEFT_TOP, "No RDF data.", font, label);
        return;
    }
    let g = egui::Rect::from_min_max(rect.min + egui::vec2(52.0, 16.0), rect.max - egui::vec2(16.0, 28.0));
    painter.rect_stroke(g, 0.0, egui::Stroke::new(1.0, if light { rgba(150, 155, 170, 255) } else { rgba(80, 90, 110, 255) }), egui::StrokeKind::Middle);
    let mut max_y = 0.0f32;
    for b in &r.bins {
        max_y = max_y.max(b.g);
        if show_raw {
            max_y = max_y.max(b.raw_count);
        }
        if show_cumulative {
            max_y = max_y.max(b.cumulative);
        }
    }
    let max_y = max_y.max(1.0);
    let grid = if light { rgba(170, 175, 190, 180) } else { rgba(50, 60, 75, 255) };
    for k in 0..=4 {
        let t = k as f32 / 4.0;
        let y = g.max.y + (g.min.y - g.max.y) * t;
        painter.line_segment([egui::pos2(g.min.x, y), egui::pos2(g.max.x, y)], egui::Stroke::new(1.0, grid));
        painter.text(egui::pos2(rect.min.x + 4.0, y), egui::Align2::LEFT_CENTER, format!("{:.2}", max_y * t), font.clone(), label);
        let x = g.min.x + g.width() * t;
        painter.line_segment([egui::pos2(x, g.min.y), egui::pos2(x, g.max.y)], egui::Stroke::new(1.0, grid));
        painter.text(egui::pos2(x - 10.0, g.max.y + 6.0), egui::Align2::LEFT_TOP, format!("{:.2}", r.r_min + (r.r_max - r.r_min) * t), font.clone(), label);
    }
    let span = (r.r_max - r.r_min).max(1e-6);
    let curve = |value: &dyn Fn(&af_algorithms::rdf::RdfBin) -> f32| -> Vec<egui::Pos2> {
        r.bins.iter().map(|b| egui::pos2(g.min.x + g.width() * ((b.r_center - r.r_min) / span), g.max.y - g.height() * (value(b) / max_y))).collect()
    };
    let (main, raw, cumulative) = if light { (rgba(20, 120, 210, 255), rgba(200, 110, 20, 220), rgba(30, 150, 40, 220)) } else { (rgba(80, 220, 255, 255), rgba(255, 180, 80, 220), rgba(140, 255, 140, 220)) };
    painter.add(egui::Shape::line(curve(&|b| b.g), egui::Stroke::new(2.0, main)));
    if show_raw {
        painter.add(egui::Shape::line(curve(&|b| b.raw_count), egui::Stroke::new(2.0, raw)));
    }
    if show_cumulative {
        painter.add(egui::Shape::line(curve(&|b| b.cumulative), egui::Stroke::new(2.0, cumulative)));
    }
    painter.text(g.min + egui::vec2(8.0, 6.0), egui::Align2::LEFT_TOP, if r.normalized { "g(r)" } else { "Histogram" }, font.clone(), main);
    if show_raw {
        painter.text(g.min + egui::vec2(64.0, 6.0), egui::Align2::LEFT_TOP, "Counts", font.clone(), raw);
    }
    if show_cumulative {
        painter.text(g.min + egui::vec2(130.0, 6.0), egui::Align2::LEFT_TOP, "Cumulative CN", font, cumulative);
    }
}

/// A bordered child panel of a fixed size, scrolling its content (the C++
/// responsive::beginChild with a border).
fn panel(ui: &mut egui::Ui, size: egui::Vec2, id: &str, add: impl FnOnce(&mut egui::Ui)) {
    let p = crate::theme::palette(ui.visuals().dark_mode);
    ui.allocate_ui_with_layout(size, egui::Layout::top_down(egui::Align::Min), |ui| {
        egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| {
            ui.set_width(size.x - 22.0);
            ui.set_height(size.y - 18.0);
            egui::ScrollArea::vertical().id_salt(id).auto_shrink([false, false]).show(ui, add);
        });
    });
}

/// ImGui::DragFloatRange2: two framed drags across `width` ("r_min = 1.00",
/// "r_max = 8.00"), the label after them.
fn range2(ui: &mut egui::Ui, width: f32, low: &mut f32, high: &mut f32, speed: f64, decimals: usize, label: &str) {
    ui.horizontal(|ui| {
        crate::dialog_widgets::field_style(ui);
        let each = ((width - ui.spacing().item_spacing.x) / 2.0).max(40.0);
        let (lo, hi) = (*low, *high);
        ui.add_sized([each, crate::theme::FRAME_HEIGHT], egui::DragValue::new(low).speed(speed).range(0.0..=hi).prefix("r_min = ").fixed_decimals(decimals));
        ui.add_sized([each, crate::theme::FRAME_HEIGHT], egui::DragValue::new(high).speed(speed).range(lo..=50.0).prefix("r_max = ").fixed_decimals(decimals));
        ui.label(label);
    });
}
