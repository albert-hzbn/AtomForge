//! Analysis > Structure and Defects > Angular Distribution Function (the C++
//! ui/AngularDistributionAnalysis): bond angles j-i-k within a cutoff, for
//! all centres, one element or a triplet of elements; smoothing and
//! normalisation; the distribution with reference geometries and peak
//! markers, the peak table, coordination statistics and a guide.

use std::sync::mpsc::{channel, Receiver};

use af_algorithms::adf::{compute_adf, AdfCentreMode, AdfParams, AdfResult};
use af_core::structure::Structure;
use eframe::egui;

pub struct AdfDialog {
    pub open: bool,
    pub params: AdfParams,
    show_peaks: bool,
    show_raw: bool,
    x_min: f32,
    x_max: f32,
    tab: usize,
    dirty: bool,
    running: Option<Receiver<AdfResult>>,
    pub result: Option<AdfResult>,
}

impl Default for AdfDialog {
    fn default() -> Self {
        Self { open: false, params: AdfParams::default(), show_peaks: true, show_raw: false, x_min: 0.0, x_max: 180.0, tab: 0, dirty: true, running: None, result: None }
    }
}

/// Reference bond angles and their abbreviations.
const REFERENCES: [(f32, &str, &str, &str); 7] = [
    (60.0, "60°", "tri", "Close-packed triangular layer"),
    (70.5, "70.5°", "fcc", "FCC / HCP second peak"),
    (90.0, "90°", "oct", "Octahedral (NaCl, MgO)"),
    (109.47, "109.5°", "tet", "Tetrahedral (diamond, ZnS)"),
    (120.0, "120°", "trpl", "Trigonal planar (graphene)"),
    (150.0, "150°", "bcc", "BCC second-shell"),
    (180.0, "180°", "lin", "Linear (CO₂, triple bonds)"),
];

fn elements_of(s: &Structure) -> Vec<String> {
    let mut list: Vec<String> = Vec::new();
    for a in &s.atoms {
        if !list.contains(&a.symbol) {
            list.push(a.symbol.clone());
        }
    }
    list.sort();
    list
}

impl AdfDialog {
    /// Runs the analysis on this thread (UI scripts).
    pub fn run_now(&mut self, structure: &Structure) -> &AdfResult {
        self.result = Some(compute_adf(structure, &self.params));
        self.dirty = false;
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
        let mut open = self.open;
        crate::theme::modal_window("Angular Distribution Function").open(&mut open).fixed_size([920.0, 660.0]).show(ctx, |ui| {
            // The settings in a bordered 300-wide child, the tabs beside it, a
            // Close footer (the C++ ADF window).
            let pal = crate::theme::palette(ui.visuals().dark_mode);
            let height = (ui.available_height() - 36.0 - ui.spacing().item_spacing.y).max(300.0);
            ui.horizontal_top(|ui| {
                ui.allocate_ui_with_layout(egui::vec2(273.0, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                    egui::Frame::new().fill(pal.surface).stroke(egui::Stroke::new(1.0, pal.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| {
                        ui.set_width(253.0);
                        ui.set_height(height - 18.0);
                        crate::dialog_widgets::separator_text(ui, "ADF Settings");
                        let compute_area = crate::theme::FRAME_HEIGHT * 2.4 + 16.0;
                        let settings_height = (ui.available_height() - compute_area).max(60.0);
                        egui::ScrollArea::vertical().id_salt("adf settings").max_height(settings_height).auto_shrink([false, false]).show(ui, |ui| self.settings(ui, structure));
                        ui.separator();
                        self.compute_area(ui, structure);
                    });
                });
                let rest = ui.available_width();
                ui.allocate_ui_with_layout(egui::vec2(rest, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                    ui.set_width(rest);
                    crate::dialog_layout::tab_bar(ui, &mut self.tab, &["Distribution", "Peaks", "Coordination", "Guide"]);
                    egui::ScrollArea::vertical().id_salt("adf tab").auto_shrink([false, false]).show(ui, |ui| match self.tab {
                        0 => {
                            self.plot(ui);
                            self.summary(ui);
                        }
                        1 => self.peak_table(ui),
                        2 => self.coordination(ui),
                        _ => guide(ui),
                    });
                });
            });
            if crate::dialog_widgets::footer(ui, "Close", None, true) == 1 {
                self.open = false;
            }
        });
        self.open = open && self.open;
    }

    fn settings(&mut self, ui: &mut egui::Ui, structure: &Structure) {
        let busy = self.running.is_some();
        let before = self.params.clone();
        let elements = elements_of(structure);
        ui.add_enabled_ui(!busy, |ui| {
            let p = &mut self.params;
            use crate::dialog_widgets as w;
            let label_w = 82.0;
            let field = |ui: &mut egui::Ui, label: &str, add: &mut dyn FnMut(&mut egui::Ui, f32) -> egui::Response| -> egui::Response {
                ui.horizontal(|ui| {
                    ui.add_sized([label_w, crate::theme::FRAME_HEIGHT], egui::Label::new(label).halign(egui::Align::Min));
                    let width = ui.available_width();
                    add(ui, width)
                })
                .inner
            };
            w::separator_text(ui, "Computation");
            field(ui, "Cutoff", &mut |ui, width| {
                w::field_style(ui);
                ui.add_sized([width, crate::theme::FRAME_HEIGHT], egui::DragValue::new(&mut p.r_cutoff).speed(0.05).range(0.5..=15.0).suffix(" Å").fixed_decimals(2))
            })
            .on_hover_text("Maximum neighbour distance (Å) for finding bonded triplets.\nDrag or Ctrl+click to type a value.");
            let bins_tip = format!("{} bins, {:.2}° per bin.", p.bin_count, 180.0 / p.bin_count.max(1) as f32);
            field(ui, "Bins", &mut |ui, width| {
                w::field_style(ui);
                ui.add_sized([width, crate::theme::FRAME_HEIGHT], egui::DragValue::new(&mut p.bin_count).range(30..=360))
            })
            .on_hover_text(bins_tip);
            let has_cell = structure.has_unit_cell;
            field(ui, "Use PBC", &mut |ui, _| ui.add_enabled_ui(has_cell, |ui| w::checkbox(ui, &mut p.use_pbc, "")).inner)
                .on_hover_text("Apply minimum-image PBC when searching for neighbours.")
                .on_disabled_hover_text("No unit cell: PBC unavailable.");
            w::separator_text(ui, "Atom Selection");
            w::radio_value(ui, &mut p.centre_mode, AdfCentreMode::All, "All").on_hover_text("Every atom is a valid triplet centre.");
            w::radio_value(ui, &mut p.centre_mode, AdfCentreMode::ByElement, "Element").on_hover_text("Only atoms of the chosen element type are centres.");
            w::radio_value(ui, &mut p.centre_mode, AdfCentreMode::ByPair, "Triplet").on_hover_text("Specify centre and both neighbour element types for j-i-k triplets.");
            let combo = |ui: &mut egui::Ui, id: &str, value: &mut String| {
                egui::ComboBox::from_id_salt(id).selected_text(if value.is_empty() { "(any)".to_string() } else { value.clone() }).show_ui(ui, |ui| {
                    ui.selectable_value(value, String::new(), "(any)");
                    for e in &elements {
                        ui.selectable_value(value, e.clone(), e);
                    }
                });
            };
            if p.centre_mode != AdfCentreMode::All {
                egui::Grid::new("adf filter").num_columns(2).show(ui, |ui| {
                    ui.label("Centre:");
                    combo(ui, "adf centre", &mut p.centre_symbol);
                    ui.end_row();
                    if p.centre_mode == AdfCentreMode::ByPair {
                        ui.label("Neigh j:");
                        combo(ui, "adf j", &mut p.neigh_symbol1);
                        ui.end_row();
                        ui.label("Neigh k:");
                        combo(ui, "adf k", &mut p.neigh_symbol2);
                        ui.end_row();
                    }
                });
            }
            w::separator_text(ui, "Post-processing");
            field(ui, "Smooth", &mut |ui, width| {
                let mut v = p.smooth_passes as f64;
                let text = format!("{} pass", p.smooth_passes);
                let r = w::slider(ui, &mut v, 0.0..=8.0, width, &text, "", true);
                p.smooth_passes = v.round() as _;
                r
            })
            .on_hover_text("3-point box smoothing passes (0 = off). 2-4 recommended.");
            field(ui, "Normalise", &mut |ui, _| w::checkbox(ui, &mut p.normalize, "")).on_hover_text("Scale histogram so tallest bin = 1. Disable to compare raw counts across runs.");
        });
        if self.params != before {
            self.dirty = true;
        }
        crate::dialog_widgets::separator_text(ui, "Display");
        crate::dialog_widgets::checkbox(ui, &mut self.show_peaks, "Peak markers");
        crate::dialog_widgets::checkbox(ui, &mut self.show_raw, "Raw counts");
        crate::dialog_widgets::disabled(ui, "Angle range:");
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut self.x_min).speed(0.5).range(0.0..=180.0).suffix("°").fixed_decimals(0));
            ui.add(egui::DragValue::new(&mut self.x_max).speed(0.5).range(0.0..=180.0).suffix("°").fixed_decimals(0)).on_hover_text("Zoom the plot to an angle sub-range. Does not affect the computation.");
            if self.x_max - self.x_min < 5.0 {
                self.x_max = (self.x_min + 5.0).min(180.0);
                self.x_min = self.x_min.min(self.x_max - 5.0);
            }
            if ui.button("All").on_hover_text("Reset to the full 0-180° view.").clicked() {
                self.x_min = 0.0;
                self.x_max = 180.0;
            }
        });
        if let Some(r) = self.result.as_ref().filter(|r| !r.element_counts.is_empty()) {
            crate::dialog_widgets::separator_text(ui, "Composition");
            let total = r.n_atoms.max(1) as f32;
            for (symbol, &count) in &r.element_counts {
                let fraction = count as f32 / total;
                ui.label(format!("{symbol:<2}  {:.1}%", fraction * 100.0));
                let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), 6.0), egui::Sense::hover());
                let fill = egui::Rect::from_min_size(rect.min, egui::vec2(rect.width() * fraction, 6.0));
                ui.painter().rect_filled(fill, 0.0, egui::Color32::from_rgba_unmultiplied(80, 200, 255, 200));
                ui.painter().rect_stroke(rect, 0.0, egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(90, 100, 120, 180)), egui::StrokeKind::Inside);
            }
        }
    }

    fn compute_area(&mut self, ui: &mut egui::Ui, structure: &Structure) {
        let busy = self.running.is_some();
        let width = ui.available_width();
        let response = if busy {
            let spin = ['|', '/', '-', '\\'][((ui.input(|i| i.time) * 4.0) as usize) % 4];
            ui.add_enabled(false, egui::Button::new(format!("Computing {spin}")).min_size(egui::vec2(width, crate::theme::FRAME_HEIGHT)))
        } else {
            ui.add_enabled_ui(!structure.atoms.is_empty(), |ui| crate::dialog_layout::wide_primary(ui, "Compute ADF", crate::theme::FRAME_HEIGHT, true)).inner
        };
        let response = response.on_hover_text(if self.dirty { "Parameters changed: click to recompute." } else { "Run the Angular Distribution Function analysis." });
        if response.clicked() {
            let (params, snapshot) = (self.params.clone(), structure.clone());
            let (tx, rx) = channel();
            std::thread::spawn(move || {
                let _ = tx.send(compute_adf(&snapshot, &params));
            });
            self.running = Some(rx);
            self.dirty = false;
        }
        match (&self.result, busy) {
            (_, true) => {
                ui.weak("Running on background thread...");
            }
            (Some(r), false) if r.valid => {
                ui.horizontal(|ui| {
                    ui.colored_label(egui::Color32::from_rgb(70, 170, 90), "Complete");
                    ui.weak(format!("{} triplets", r.n_triplets));
                });
            }
            (Some(r), false) => {
                ui.horizontal(|ui| {
                    ui.colored_label(egui::Color32::from_rgb(220, 80, 60), "Failed:");
                    ui.weak(&r.message);
                });
            }
            (None, false) => {
                ui.weak("Click Compute ADF to start.");
            }
        }
    }

    fn plot(&self, ui: &mut egui::Ui) {
        let light = !ui.visuals().dark_mode;
        let rgba = egui::Color32::from_rgba_unmultiplied;
        let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width().max(60.0), 270.0), egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, if light { rgba(245, 246, 248, 255) } else { rgba(20, 24, 31, 255) });
        painter.rect_stroke(rect, 0.0, egui::Stroke::new(1.0, if light { rgba(140, 145, 160, 255) } else { rgba(110, 125, 150, 255) }), egui::StrokeKind::Inside);
        let font = egui::FontId::proportional(12.0);
        let label = if light { rgba(50, 50, 65, 255) } else { rgba(200, 200, 200, 255) };
        let Some(r) = self.result.as_ref().filter(|r| !r.bins.is_empty()) else {
            let text = if self.running.is_some() { "Computing, please wait..." } else { "No data yet: set the parameters, then Compute ADF." };
            painter.text(rect.center(), egui::Align2::CENTER_CENTER, text, font, label);
            return;
        };
        let g = egui::Rect::from_min_max(rect.min + egui::vec2(54.0, 22.0), rect.max - egui::vec2(16.0, 34.0));
        painter.rect_stroke(g, 0.0, egui::Stroke::new(1.0, if light { rgba(150, 155, 170, 255) } else { rgba(80, 90, 110, 255) }), egui::StrokeKind::Middle);
        let range = (self.x_max - self.x_min).max(1.0);
        let visible = |a: f32| a >= self.x_min && a <= self.x_max;
        let value = |b: &af_algorithms::adf::AdfBin| if self.show_raw { b.count } else { b.value };
        let max_y = r.bins.iter().filter(|b| visible(b.angle_deg)).map(value).fold(0.0f32, f32::max).max(1e-6);
        let grid = if light { rgba(170, 175, 190, 180) } else { rgba(50, 60, 75, 255) };
        for k in 0..=4 {
            let t = k as f32 / 4.0;
            let y = g.max.y + (g.min.y - g.max.y) * t;
            painter.line_segment([egui::pos2(g.min.x, y), egui::pos2(g.max.x, y)], egui::Stroke::new(1.0, grid));
            painter.text(egui::pos2(rect.min.x + 3.0, y), egui::Align2::LEFT_CENTER, format!("{:.2}", max_y * t), font.clone(), label);
        }
        let px = |a: f32| g.min.x + g.width() * (a - self.x_min) / range;
        let (ref_line, ref_label) = if light { (rgba(130, 142, 160, 130), rgba(110, 122, 145, 210)) } else { (rgba(65, 75, 95, 200), rgba(120, 138, 165, 220)) };
        let mut last_right = f32::MIN;
        for (angle, _, abbr, _) in REFERENCES {
            if !visible(angle) {
                continue;
            }
            let x = px(angle);
            for s in [0, 2, 4] {
                let y0 = g.min.y + g.height() * s as f32 / 5.0;
                let y1 = g.min.y + g.height() * (s + 1) as f32 / 5.0;
                painter.line_segment([egui::pos2(x, y0), egui::pos2(x, y1)], egui::Stroke::new(1.0, ref_line));
            }
            let width = painter.layout_no_wrap(abbr.to_string(), font.clone(), ref_label).size().x;
            if x - width / 2.0 > last_right + 2.0 {
                painter.text(egui::pos2(x - width / 2.0, g.min.y + 2.0), egui::Align2::LEFT_TOP, abbr, font.clone(), ref_label);
                last_right = x + width / 2.0;
            }
        }
        let points: Vec<egui::Pos2> = r.bins.iter().filter(|b| visible(b.angle_deg)).map(|b| egui::pos2(px(b.angle_deg), g.max.y - g.height() * (value(b) / max_y).clamp(0.0, 1.0))).collect();
        let (line, fill) = if light { (rgba(20, 120, 210, 255), rgba(20, 120, 210, 40)) } else { (rgba(80, 220, 255, 255), rgba(80, 220, 255, 28)) };
        for pair in points.windows(2) {
            painter.add(egui::Shape::convex_polygon(vec![egui::pos2(pair[0].x, g.max.y), egui::pos2(pair[1].x, g.max.y), pair[1], pair[0]], fill, egui::Stroke::NONE));
        }
        painter.add(egui::Shape::line(points, egui::Stroke::new(2.0, line)));
        if self.show_peaks {
            let (peak, peak_line) = if light { (rgba(210, 35, 35, 220), rgba(210, 35, 35, 70)) } else { (rgba(255, 100, 100, 230), rgba(255, 100, 100, 55)) };
            for p in r.peaks.iter().filter(|p| visible(p.angle_deg)) {
                let x = px(p.angle_deg);
                let y = g.max.y - g.height() * (p.value / max_y).clamp(0.0, 1.0);
                painter.line_segment([egui::pos2(x, g.min.y + 2.0), egui::pos2(x, y - 6.0)], egui::Stroke::new(1.5, peak_line));
                painter.circle(egui::pos2(x, y), 5.0, peak, egui::Stroke::new(1.0, if light { rgba(170, 20, 20, 200) } else { rgba(255, 170, 170, 200) }));
                painter.text(egui::pos2(x, (y - 18.0).max(g.min.y + 14.0)), egui::Align2::CENTER_TOP, format!("{:.1}°", p.angle_deg), font.clone(), peak);
            }
        }
        let y_label = if self.show_raw { "Counts" } else if r.normalized { "g(θ) [norm]" } else { "g(θ)" };
        painter.text(egui::pos2(g.min.x + 5.0, rect.min.y + 4.0), egui::Align2::LEFT_TOP, y_label, font.clone(), line);
        painter.text(egui::pos2(g.center().x, g.max.y + 20.0), egui::Align2::CENTER_TOP, "Angle (°)", font.clone(), label);
        if let Some(m) = response.hover_pos().filter(|m| m.x >= g.min.x && m.x <= g.max.x) {
            let angle = self.x_min + (m.x - g.min.x) / g.width() * range;
            let index = ((angle / r.bin_width.max(1e-4)) as usize).min(r.bins.len() - 1);
            let bin = &r.bins[index];
            let v = value(bin);
            let cross = if light { rgba(80, 88, 108, 160) } else { rgba(185, 195, 215, 150) };
            painter.line_segment([egui::pos2(m.x, g.min.y), egui::pos2(m.x, g.max.y)], egui::Stroke::new(1.0, cross));
            painter.circle_filled(egui::pos2(m.x, g.max.y - g.height() * (v / max_y).clamp(0.0, 1.0)), 3.5, cross);
            let mut text = format!("θ = {:.2}°\ng(θ) = {:.5}", bin.angle_deg, v);
            if r.normalized && bin.count > 0.0 {
                text += &format!("\nRaw count: {:.0}", bin.count);
            }
            response.on_hover_text(text);
        }
    }

    fn summary(&self, ui: &mut egui::Ui) {
        let Some(r) = &self.result else { return };
        ui.separator();
        ui.weak(format!("Atoms: {}  |  Centres: {}  |  Triplets: {}  |  Bins: {}  ({:.2}° / bin)", r.n_atoms, r.n_centre_atoms, r.n_triplets, r.bin_count, r.bin_width));
        ui.weak(format!("{}  |  Cutoff: {:.2} Å", if r.pbc_used { "PBC applied" } else { "No PBC" }, r.r_cutoff));
    }

    fn peak_table(&self, ui: &mut egui::Ui) {
        let Some(r) = self.result.as_ref().filter(|r| !r.peaks.is_empty()) else {
            ui.weak("No peaks detected. Try reducing smooth passes or recompute with finer bins.");
            return;
        };
        ui.weak(format!("{} peak(s) detected:", r.peaks.len()));
        egui::Grid::new("adf peaks").num_columns(3).striped(true).show(ui, |ui| {
            ui.strong("Angle");
            ui.strong("Intensity");
            ui.strong("Geometry");
            ui.end_row();
            for p in &r.peaks {
                ui.label(format!("{:.2}°", p.angle_deg));
                ui.horizontal(|ui| {
                    ui.label(format!("{:.4}", p.value));
                    let (rect, _) = ui.allocate_exact_size(egui::vec2(120.0, 12.0), egui::Sense::hover());
                    ui.painter().rect_filled(egui::Rect::from_min_size(rect.min, egui::vec2(120.0 * p.value.min(1.0), 12.0)), 0.0, egui::Color32::from_rgba_unmultiplied(80, 220, 255, 160));
                });
                ui.label(if p.label.is_empty() { "—" } else { &p.label });
                ui.end_row();
            }
        });
    }

    fn coordination(&self, ui: &mut egui::Ui) {
        let Some(r) = self.result.as_ref().filter(|r| !r.coord_stats.is_empty()) else {
            ui.weak("No coordination data. Run a computation first.");
            return;
        };
        ui.weak(format!("Mean coordination within cutoff {:.2} Å:", r.r_cutoff));
        let max_cn = r.coord_stats.values().map(|c| c.mean).fold(1.0f32, f32::max);
        egui::Grid::new("adf coordination").num_columns(5).striped(true).show(ui, |ui| {
            for h in ["Elem", "Mean CN", "", "± Std", "Centres"] {
                ui.strong(h);
            }
            ui.end_row();
            for (symbol, c) in &r.coord_stats {
                ui.label(symbol);
                ui.label(format!("{:.2}", c.mean));
                let (rect, _) = ui.allocate_exact_size(egui::vec2(120.0, 12.0), egui::Sense::hover());
                ui.painter().rect_filled(egui::Rect::from_min_size(rect.min, egui::vec2(120.0 * c.mean / max_cn, 12.0)), 0.0, egui::Color32::from_rgba_unmultiplied(80, 230, 120, 140));
                ui.label(format!("{:.2}", c.stddev));
                ui.label(c.count.to_string());
                ui.end_row();
            }
        });
    }
}

fn guide(ui: &mut egui::Ui) {
    ui.strong("Method");
    ui.label("The ADF g(θ) accumulates the frequency of bond angles θ formed by atom triplets (j-i-k), where i is the centre atom and both j and k lie within the cutoff radius.");
    ui.strong("Reference Geometries");
    egui::Grid::new("adf guide").num_columns(3).striped(true).show(ui, |ui| {
        for h in ["Angle", "Abbrev.", "Geometry / Example"] {
            ui.strong(h);
        }
        ui.end_row();
        for (_, angle, abbr, description) in REFERENCES {
            ui.label(angle);
            ui.label(abbr);
            ui.label(description);
            ui.end_row();
        }
    });
    ui.strong("Tips");
    ui.label(
        "Cutoff: set to the first-shell peak of the RDF for best results.\n\nBins: 90-180 is a good default. More bins give finer resolution but noisier histograms.\n\nSmooth passes: 2-4 passes reduce noise without hiding real peaks. Set to 0 to see the raw histogram.\n\nAtom Selection: use Element or Triplet mode to isolate the ADF around a specific species in multi-element structures.\n\nAngle range: drag the range values to zoom into a region of interest without recomputing.",
    );
}
