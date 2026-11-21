//! Analysis > Structure and Defects > Short-Range Order (the C++
//! ui/ShortRangeOrderDialog): Warren-Cowley and Rao-Curtin parameters per
//! neighbour shell, the shells, the composition and a guide; alpha shown on
//! a blue (ordering) - grey (random) - red (clustering) scale.

use std::sync::mpsc::{channel, Receiver};

use af_algorithms::sro::{analyze_sro, SroAnalysisResult};
use af_core::structure::Structure;
use eframe::egui;

pub struct SroDialog {
    pub open: bool,
    pub shells: i32,
    pub tolerance: f64,
    tab: usize,
    running: Option<Receiver<SroAnalysisResult>>,
    pub result: Option<SroAnalysisResult>,
}

impl Default for SroDialog {
    fn default() -> Self {
        Self { open: false, shells: 2, tolerance: 0.05, tab: 0, running: None, result: None }
    }
}

/// The diverging colour of alpha in [-1, 1] (sroAlphaColor).
fn alpha_colour(alpha: f64, light: bool) -> egui::Color32 {
    let t = alpha.clamp(-1.0, 1.0) as f32;
    let (r, g, b) = if light {
        if t < 0.0 {
            let s = -t;
            (0.60 - s * 0.40, 0.60 - s * 0.20, 0.60 + s * 0.40)
        } else {
            (0.60 + t * 0.40, 0.60 - t * 0.40, 0.60 - t * 0.50)
        }
    } else if t < 0.0 {
        let s = -t;
        (1.0 - s * 0.55, 1.0 - s * 0.45, 1.0)
    } else {
        (1.0, 1.0 - t * 0.72, 1.0 - t * 0.88)
    };
    crate::theme::srgb(r, g, b).into()
}

/// A horizontal blue-grey-red gradient in `rect`.
fn gradient(painter: &egui::Painter, rect: egui::Rect, alpha: u8) {
    let (blue, grey, red) = (egui::Color32::from_rgba_unmultiplied(55, 95, 215, alpha), egui::Color32::from_rgba_unmultiplied(205, 205, 205, alpha), egui::Color32::from_rgba_unmultiplied(205, 55, 55, alpha));
    let mut mesh = egui::Mesh::default();
    let mid = rect.center().x;
    for (x0, x1, c0, c1) in [(rect.min.x, mid, blue, grey), (mid, rect.max.x, grey, red)] {
        let base = mesh.vertices.len() as u32;
        mesh.colored_vertex(egui::pos2(x0, rect.min.y), c0);
        mesh.colored_vertex(egui::pos2(x1, rect.min.y), c1);
        mesh.colored_vertex(egui::pos2(x1, rect.max.y), c1);
        mesh.colored_vertex(egui::pos2(x0, rect.max.y), c0);
        mesh.add_triangle(base, base + 1, base + 2);
        mesh.add_triangle(base, base + 2, base + 3);
    }
    painter.add(egui::Shape::mesh(mesh));
}

/// Where alpha sits in [-1, 1] (drawAlphaBar).
fn alpha_bar(ui: &mut egui::Ui, alpha: f64, light: bool) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(84.0, 8.0), egui::Sense::hover());
    let painter = ui.painter();
    gradient(painter, rect, 170);
    painter.rect_stroke(rect, 1.5, egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(110, 110, 110, 150)), egui::StrokeKind::Inside);
    let x = rect.min.x + ((alpha.clamp(-1.0, 1.0) as f32 + 1.0) * 0.5) * rect.width();
    let needle = if light { egui::Color32::from_rgb(30, 30, 30) } else { egui::Color32::WHITE };
    painter.line_segment([egui::pos2(x, rect.min.y - 1.0), egui::pos2(x, rect.max.y + 1.0)], egui::Stroke::new(2.0, needle));
    painter.add(egui::Shape::convex_polygon(vec![egui::pos2(x - 3.0, rect.min.y - 4.0), egui::pos2(x + 3.0, rect.min.y - 4.0), egui::pos2(x, rect.min.y - 1.0)], needle, egui::Stroke::NONE));
}

fn legend(ui: &mut egui::Ui, light: bool) {
    ui.horizontal(|ui| {
        ui.weak("αᵢⱼ colormap:");
        let (rect, _) = ui.allocate_exact_size(egui::vec2(260.0, 32.0), egui::Sense::hover());
        let bar = egui::Rect::from_min_size(rect.min, egui::vec2(260.0, 13.0));
        gradient(ui.painter(), bar, 235);
        let tick = if light { egui::Color32::from_rgb(30, 30, 30) } else { egui::Color32::from_rgb(215, 215, 215) };
        for (f, label) in [(0.0, "-1"), (0.25, "-0.5"), (0.5, "0"), (0.75, "+0.5"), (1.0, "+1")] {
            let x = bar.min.x + f * bar.width();
            ui.painter().line_segment([egui::pos2(x, bar.max.y), egui::pos2(x, bar.max.y + 4.0)], egui::Stroke::new(1.0, tick));
            ui.painter().text(egui::pos2(x, bar.max.y + 5.0), egui::Align2::CENTER_TOP, label, egui::FontId::proportional(11.0), tick);
        }
    });
    ui.horizontal(|ui| {
        ui.colored_label(egui::Color32::from_rgb(64, 115, 255), "Ordering");
        ui.weak("Random");
        ui.colored_label(egui::Color32::from_rgb(255, 71, 71), "Clustering");
    });
}

impl SroDialog {
    pub fn run_now(&mut self, structure: &Structure) -> &SroAnalysisResult {
        self.result = Some(analyze_sro(structure, self.shells, self.tolerance));
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
        let light = ctx.theme() == egui::Theme::Light;
        let periodic = !structure.atoms.is_empty() && structure.has_unit_cell;
        let mut open = self.open;
        crate::theme::modal_window("Short-Range Order Analysis").open(&mut open).fixed_size([920.0, 680.0]).show(ctx, |ui| {
            // Settings in a bordered 240-wide child, the result tabs beside it,
            // a Close footer (the C++ SRO window).
            let pal = crate::theme::palette(ui.visuals().dark_mode);
            let footer = crate::theme::FRAME_HEIGHT + 3.0 * ui.spacing().item_spacing.y + 4.0;
            let height = (ui.available_height() - footer).max(260.0);
            ui.horizontal_top(|ui| {
                ui.allocate_ui_with_layout(egui::vec2(218.0, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                    egui::Frame::new().fill(pal.surface).stroke(egui::Stroke::new(1.0, pal.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| {
                        ui.set_width(198.0);
                        ui.set_height(height - 18.0);
                        egui::ScrollArea::vertical().id_salt("sro left").auto_shrink([false, false]).show(ui, |ui| self.left(ui, structure, periodic));
                    });
                });
                let rest = ui.available_width();
                ui.allocate_ui_with_layout(egui::vec2(rest, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                    ui.set_width(rest);
                    ui.set_height(height);
                    self.right(ui, periodic, light);
                });
            });
            ui.add_space(2.0);
            ui.separator();
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if crate::dialog_widgets::primary_button(ui, "Close", egui::vec2(99.0, crate::theme::FRAME_HEIGHT)).clicked() {
                    self.open = false;
                }
            });
        });
        self.open = open && self.open;
    }

    fn left(&mut self, ui: &mut egui::Ui, structure: &Structure, periodic: bool) {
        use crate::dialog_widgets as w;
        let busy = self.running.is_some();
        let accent = crate::theme::palette(ui.visuals().dark_mode).accent;
        ui.label(egui::RichText::new("Settings").color(accent));
        ui.separator();
        let width = ui.available_width();
        ui.add_enabled_ui(!busy, |ui| {
            w::disabled(ui, "Neighbor shells");
            let mut shells = self.shells as i32;
            w::slider_i32(ui, &mut shells, 1..=8, width, "").on_hover_text("Number of nearest-neighbor shells to include. Higher values are slower on large structures.");
            self.shells = shells.max(1) as _;
            ui.add_space(4.0);
            w::disabled(ui, "Shell tolerance (Å)");
            let mut tolerance = self.tolerance as f32;
            w::slider_f32(ui, &mut tolerance, 0.001..=0.20, width, 3, "").on_hover_text("Distance window for merging atoms into the same shell. Increase if shells appear split.");
            self.tolerance = tolerance as _;
        });
        ui.add_space(4.0);
        ui.separator();
        let run = egui::Button::new(egui::RichText::new("Run Analysis").color(egui::Color32::WHITE)).fill(egui::Color32::from_rgb(33, 122, 209)).min_size(egui::vec2(width, 32.0));
        if ui.add_enabled(periodic && !busy, run).clicked() {
            let (snapshot, shells, tolerance) = (structure.clone(), self.shells, self.tolerance);
            let (tx, rx) = channel();
            std::thread::spawn(move || {
                let _ = tx.send(analyze_sro(&snapshot, shells, tolerance));
            });
            self.running = Some(rx);
            self.result = None;
        }
        let warn = egui::Color32::from_rgb(230, 170, 60);
        if busy {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.colored_label(warn, "Computing...");
            });
            ui.weak("Analysis running in background. UI remains responsive.");
        } else if !periodic {
            ui.colored_label(warn, "No unit-cell structure loaded.");
        } else if let Some(r) = self.result.as_ref().filter(|r| r.success) {
            ui.colored_label(egui::Color32::from_rgb(46, 224, 71), "Complete");
            egui::Grid::new("sro counts").show(ui, |ui| {
                ui.weak("Atoms");
                ui.label(r.warren_cowley.n_atoms.to_string());
                ui.end_row();
                ui.weak("Shells");
                ui.label(r.warren_cowley.n_shells.to_string());
                ui.end_row();
            });
            ui.label(egui::RichText::new("Composition").color(accent));
            ui.separator();
            for (element, &c) in &r.warren_cowley.composition {
                ui.horizontal(|ui| {
                    ui.label(format!("{element:<2}  {:.1}%", c * 100.0));
                    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width() - 10.0, 10.0), egui::Sense::hover());
                    ui.painter().rect_filled(egui::Rect::from_min_size(rect.min, egui::vec2(rect.width() * c as f32, 10.0)), 2.0, egui::Color32::from_rgba_unmultiplied(55, 175, 255, 140));
                    ui.painter().rect_stroke(rect, 2.0, egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(75, 85, 105, 180)), egui::StrokeKind::Inside);
                });
            }
        } else if let Some(r) = &self.result {
            ui.colored_label(warn, &r.message);
        }
    }

    fn right(&mut self, ui: &mut egui::Ui, periodic: bool, light: bool) {
        let Some(r) = self.result.as_ref().filter(|r| r.success) else {
            if self.running.is_some() {
                crate::dialog_layout::empty_state(ui, "SRO analysis running...");
            } else {
                crate::dialog_layout::empty_state(ui, if periodic { "Choose the neighbour shells and tolerance, then Run Analysis." } else { "Short-range order needs a periodic unit cell; this structure has none." });
            }
            return;
        };
        crate::dialog_layout::tab_bar(ui, &mut self.tab, &[" Warren-Cowley ", " Rao-Curtin (2022) ", " Shells ", " Guide "]);
        egui::ScrollArea::both().id_salt("sro results").show(ui, |ui| match self.tab {
            0 => {
                ui.weak("αᵢⱼ(s) = 1 − P(j|i,s) / cⱼ   [Cowley, Phys. Rev. 138 (1965) A1384]");
                let wc = &r.warren_cowley;
                for shell in 0..wc.n_shells.max(wc.shell_distances.len() as i32) {
                    let entries: Vec<_> = wc.entries.iter().filter(|e| e.shell == shell).collect();
                    let Some(first) = entries.first() else { continue };
                    crate::dialog_layout::collapsing(ui, ("wc", shell), &format!(" Shell {}  —  {:.4} Å", shell + 1, first.shell_distance), true, |ui| {
                        crate::dialog_layout::table(ui, &["Center", "Neighbor", "Nᵢⱼ", "Zᵢ", "P(j|i)", "αᵢⱼ"], &[62.0, 62.0, 58.0, 52.0, 70.0], entries.len(), |ui, row, col| {
                            let e = entries[row];
                            match col {
                                0 => {
                                    ui.label(&e.center_element);
                                }
                                1 => {
                                    ui.label(&e.neighbor_element);
                                }
                                2 => {
                                    ui.label(e.nij.to_string());
                                }
                                3 => {
                                    ui.label(e.zi.to_string());
                                }
                                4 => {
                                    ui.label(format!("{:.4}", e.pij));
                                }
                                _ => {
                                    ui.colored_label(alpha_colour(e.alpha, light), format!("{:>7.4}", e.alpha));
                                    alpha_bar(ui, e.alpha, light);
                                }
                            }
                        });
                    });
                }
                ui.separator();
                legend(ui, light);
            }
            1 => {
                ui.weak("αᵢⱼ(s) = 1 − Pᵢⱼ(s) / (2·cᵢ·cⱼ)   [Rao & Curtin, Acta Mater. 226 (2022) 117621]");
                let rc = &r.rao_curtin;
                for shell in 0..rc.n_shells.max(rc.shell_distances.len() as i32) {
                    let entries: Vec<_> = rc.entries.iter().filter(|e| e.shell == shell).collect();
                    let Some(first) = entries.first() else { continue };
                    let z = rc.shell_coordination.get(shell as usize).copied().unwrap_or(0.0);
                    crate::dialog_layout::collapsing(ui, ("rc", shell), &format!(" Shell {}  —  {:.4} Å   (avg. coord. Zₛ = {z:.2})", shell + 1, first.shell_distance), true, |ui| {
                        crate::dialog_layout::table(ui, &["Pair i–j", "cᵢ", "cⱼ", "Nᵢⱼ", "Nⱼᵢ", "Pᵢⱼ (obs.)", "2cᵢcⱼ (rand)", "αᵢⱼ"], &[62.0, 58.0, 58.0, 52.0, 52.0, 78.0, 88.0], entries.len(), |ui, row, col| {
                            let e = entries[row];
                            let text = match col {
                                0 => format!("{}–{}", e.elem_i, e.elem_j),
                                1 => format!("{:.4}", e.ci),
                                2 => format!("{:.4}", e.cj),
                                3 => e.n_ij.to_string(),
                                4 => e.n_ji.to_string(),
                                5 => format!("{:.5}", e.p_ij),
                                6 => format!("{:.5}", 2.0 * e.ci * e.cj),
                                _ => {
                                    ui.colored_label(alpha_colour(e.alpha, light), format!("{:>7.4}", e.alpha));
                                    alpha_bar(ui, e.alpha, light);
                                    return;
                                }
                            };
                            ui.label(text);
                        });
                    });
                }
                ui.separator();
                legend(ui, light);
            }
            2 => {
                let (wc, rc) = (&r.warren_cowley, &r.rao_curtin);
                crate::dialog_layout::table(ui, &["Shell", "Distance (Å)", "Avg. coord. Z", "WC pairs analyzed"], &[55.0, 115.0, 115.0], wc.shell_distances.len(), |ui, s, col| {
                    match col {
                        0 => {
                            ui.label((s + 1).to_string());
                        }
                        1 => {
                            ui.label(format!("{:.4}", wc.shell_distances[s]));
                        }
                        2 => {
                            match rc.shell_coordination.get(s) {
                                Some(z) if *z >= 0.0 => ui.label(format!("{z:.3}")),
                                _ => ui.weak("—"),
                            };
                        }
                        _ => {
                            ui.label(wc.entries.iter().filter(|e| e.shell == s as i32).count().to_string());
                        }
                    }
                });
            }
            _ => guide(ui),
        });
    }
}

fn guide(ui: &mut egui::Ui) {
    let accent = crate::theme::palette(ui.visuals().dark_mode).accent;
    ui.label(egui::RichText::new("Interpreting SRO Parameters").color(accent));
    ui.separator();
    for (colour, label, description) in [
        (egui::Color32::from_rgb(64, 115, 255), "α < 0  (blue)", "Ordering tendency. Atoms i and j prefer unlike neighbors (attractive effective pair interaction, heterogeneous bonding)."),
        (egui::Color32::from_rgb(166, 166, 166), "α ≈ 0  (grey)", "Near-random distribution. Ideal solid solution, no short-range preference."),
        (egui::Color32::from_rgb(255, 71, 71), "α > 0  (red)", "Clustering tendency. Like atoms prefer like neighbors (repulsive effective pair interaction, homogeneous bonding)."),
    ] {
        ui.horizontal(|ui| {
            let (rect, _) = ui.allocate_exact_size(egui::vec2(14.0, 14.0), egui::Sense::hover());
            ui.painter().rect_filled(rect, 3.0, colour);
            ui.colored_label(colour, label);
        });
        ui.label(description);
    }
    ui.separator();
    ui.label(egui::RichText::new("Method Comparison").color(accent));
    ui.separator();
    ui.label("Warren-Cowley: directional conditional probability P(j|i,s), counting j-type neighbors around a central i atom; alpha_ij and alpha_ji can differ for each pair type.");
    ui.label("Rao-Curtin (2022): symmetric joint probability P_ij(s), combining both directions so alpha_ij = alpha_ji by construction; better aligned with free-energy perturbation theory and more consistent for multicomponent alloys.");
    ui.weak("Both methods are accurate for weak-to-moderate interactions satisfying 0.5 < exp(-V_ij/kT) < 1.5.");
}
