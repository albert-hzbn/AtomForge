//! Analysis > Structure and Defects > Common Neighbour Analysis (the C++
//! ui/CommonNeighbourAnalysis): bonds by scaled covalent radii or one fixed
//! cutoff, then the pair signatures 1-j-k-l, the environments and every
//! atom's dominant signature; computed in the background.

use std::sync::mpsc::{channel, Receiver};

use af_algorithms::cna::{compute_cna, CnaParams, CnaResult, Signature};
use af_core::structure::Structure;
use eframe::egui;

pub struct CnaDialog {
    pub open: bool,
    pub params: CnaParams,
    running: Option<Receiver<CnaResult>>,
    pub result: Option<CnaResult>,
}

impl Default for CnaDialog {
    fn default() -> Self {
        Self { open: false, params: CnaParams::default(), running: None, result: None }
    }
}

fn signature(s: &Signature) -> String {
    format!("1-{}-{}-{}", s.common, s.bonds, s.chain)
}

impl CnaDialog {
    pub fn run_now(&mut self, structure: &Structure) -> &CnaResult {
        self.result = Some(compute_cna(structure, &self.params));
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
        crate::theme::modal_window("Common Neighbour Analysis").open(&mut open).fixed_size([900.0, 720.0]).show(ctx, |ui| {
            use crate::dialog_layout as dl;
            let before = self.params.clone();
            let p = &mut self.params;
            dl::labelled_check(ui, &mut p.use_pbc, "Use PBC when unit cell is available").on_hover_text("Find neighbours through the periodic images of the cell when the structure has one.");
            let width = dl::item_width(ui);
            ui.horizontal(|ui| {
                dl::slider_f32(ui, &mut p.cutoff_scale, 1.0..=1.6, width, |v| format!("{v:.2}"));
                ui.label("Bond cutoff scale");
            })
            .response
            .on_hover_text("Two atoms are neighbours when closer than this factor times the sum of their covalent radii (used when the fixed cutoff is 0).");
            ui.horizontal(|ui| {
                let mut v = p.cutoff_radius as f64;
                if dl::number_field(ui, "cna cutoff", &mut v, 3, width).changed() {
                    p.cutoff_radius = v.max(0.0) as f32;
                }
                ui.label("Fixed cutoff (Å; 0 = covalent radii)");
            })
            .response
            .on_hover_text("One neighbour cutoff in Å for every pair. Choose it between the first and second neighbour shells (fcc/hcp) or second and third (bcc). 0 uses the scaled covalent radii.");
            dl::note(ui, dl::Note::Info, "For crystalline CNA, choose a cutoff between neighbour shells: first and second for FCC/HCP; second and third for BCC. A single unit cell works: every periodic image of a neighbour is counted. Complex graphs with more than 16 common neighbours have chain index -1.");
            if self.params != before && self.running.is_none() {
                self.result = None;
            }
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                let busy = self.running.is_some();
                if dl::accent_button(ui, "Run CNA", !busy && !structure.atoms.is_empty()).clicked() {
                    let (params, snapshot) = (self.params.clone(), structure.clone());
                    let (tx, rx) = channel();
                    std::thread::spawn(move || {
                        let _ = tx.send(compute_cna(&snapshot, &params));
                    });
                    self.running = Some(rx);
                }
                if busy {
                    ui.label(egui::RichText::new("Computing...").color(ui.visuals().warn_fg_color));
                } else if self.result.is_some() {
                    ui.label(egui::RichText::new("Done").color(dl::note_colour(ui, dl::Note::Success)));
                }
            });
            match &self.result {
                Some(r) => egui::ScrollArea::vertical().id_salt("cna results").auto_shrink([false, false]).show(ui, |ui| details(ui, r)).inner,
                None if self.running.is_none() => dl::empty_state(ui, "Run CNA to classify every atom by its common-neighbour signature (fcc, hcp, bcc, icosahedral or other)."),
                None => {}
            }
        });
        self.open = open;
    }
}

fn details(ui: &mut egui::Ui, r: &CnaResult) {
    use crate::dialog_layout as dl;
    ui.separator();
    ui.label(format!("Status: {}", r.message));
    ui.label(format!("Atoms: {}", r.atom_count));
    ui.label(format!("Bonded pairs analyzed: {}", r.pair_count));
    ui.label(format!("PBC used in analysis: {}", if r.pbc_used { "Yes" } else { "No" }));
    if !r.valid {
        return;
    }
    let fraction = |count: i32, total: i32| if total > 0 { count as f32 / total as f32 } else { 0.0 };
    ui.separator();
    ui.label("Pair Signature Distribution (Honeycutt-Andersen form: 1-j-k-l)");
    let signatures: Vec<(String, i32)> = r.signature_counts.iter().map(|(s, &c)| (signature(s), c)).collect();
    child(ui, 164.0, |ui| {
        dl::table(ui, &["Signature", "Pair Count", "Fraction"], &[120.0, 120.0, 120.0], signatures.len(), |ui, row, col| {
            let (s, c) = &signatures[row];
            ui.label(match col {
                0 => s.clone(),
                1 => c.to_string(),
                _ => format!("{:.4}", fraction(*c, r.pair_count)),
            });
        });
    });
    ui.separator();
    ui.label("Per-Atom Environment Summary");
    let environments: Vec<(String, i32)> = r.environment_counts.iter().map(|(e, &c)| (e.clone(), c)).collect();
    child(ui, 118.0, |ui| {
        dl::table(ui, &["Environment", "Atom Count", "Fraction"], &[160.0, 120.0, 120.0], environments.len(), |ui, row, col| {
            let (e, c) = &environments[row];
            ui.label(match col {
                0 => e.clone(),
                1 => c.to_string(),
                _ => format!("{:.4}", fraction(*c, r.atom_count)),
            });
        });
    });
    ui.separator();
    ui.label("Per-Atom CNA Details");
    let widths = [52.0, 56.0, 56.0, 56.0, 160.0, 90.0, 0.0];
    let headers = ["#", "El", "Z", "CN", "Dominant Signature", "Sig Count", "Environment"];
    child(ui, 218.0, |ui| {
        dl::table(ui, &headers, &widths, 0, |_, _, _| {});
        let row_h = crate::theme::FRAME_HEIGHT - 6.0;
        egui::ScrollArea::vertical().id_salt("cna atoms").auto_shrink([false, false]).show_rows(ui, row_h, r.atom_rows.len(), |ui, range| {
            ui.spacing_mut().item_spacing.y = 0.0;
            let first = range.start;
            let rows = &r.atom_rows[range];
            dl::table(ui, &["", "", "", "", "", "", ""], &widths, rows.len(), |ui, row, col| {
                let a = &rows[row];
                let _ = first;
                ui.label(match col {
                    0 => a.index.to_string(),
                    1 => a.symbol.clone(),
                    2 => a.atomic_number.to_string(),
                    3 => a.coordination.to_string(),
                    4 => if a.dominant_signature_count > 0 { signature(&a.dominant_signature) } else { "N/A".into() },
                    5 => a.dominant_signature_count.to_string(),
                    _ => a.environment.clone(),
                });
            });
        });
    });
}

/// A bordered child of a fixed height (responsive::beginChild with a border).
fn child(ui: &mut egui::Ui, height: f32, add: impl FnOnce(&mut egui::Ui)) {
    let p = crate::theme::palette(ui.visuals().dark_mode);
    egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| {
        ui.set_width(ui.available_width());
        ui.set_height(height - 16.0);
        egui::ScrollArea::vertical().id_salt(("cna child", height as i32)).auto_shrink([false, false]).show(ui, |ui| {
            ui.spacing_mut().item_spacing.y = 0.0;
            add(ui);
        });
    });
}
