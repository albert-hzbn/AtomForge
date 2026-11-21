//! The welcome panel shown when nothing is open (FileBrowser::drawWelcome):
//! the name and version, Open Structure / Build Crystal / Open Project /
//! Read Manual, the example structures and the recent files.

use std::path::{Path, PathBuf};

use eframe::egui;

use crate::theme::{self, icon};

/// What the panel asks for.
pub enum WelcomeAction {
    OpenStructure,
    BuildCrystal,
    OpenProject,
    ReadManual,
    Open(String),
}

struct Example {
    file: &'static str,
    label: &'static str,
    about: &'static str,
}

const EXAMPLES: [Example; 8] = [
    Example { file: "cu_fcc.cif", label: "Copper (fcc)", about: "The four-atom conventional cell of copper" },
    Example { file: "fe_bcc.cif", label: "Iron (bcc)", about: "The two-atom conventional cell of iron" },
    Example { file: "cu_ni_alloy.cif", label: "Cu-Ni alloy", about: "A 500-atom random Cu70Ni30 solid solution" },
    Example { file: "cu_sigma5.cif", label: "Grain boundary", about: "A copper Sigma-5 [001] tilt grain boundary" },
    Example { file: "cu_polycrystal.cif", label: "Polycrystal", about: "A four-grain copper polycrystal" },
    Example { file: "cu_dislocation.cif", label: "Dislocation", about: "A copper crystal containing an edge dislocation" },
    Example { file: "cu_sphere.xyz", label: "Nanoparticle", about: "A spherical copper nanoparticle" },
    Example { file: "sio2_pack.xyz", label: "Amorphous SiO₂", about: "A 120-atom amorphous silica packing" },
];

/// The example structures: `examples` beside the program (a package), or
/// docs/manual/examples of the source tree above it (findExamplesDirectory).
pub fn examples_directory() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    let mut dir = exe.parent()?.to_path_buf();
    for candidate in [dir.join("examples"), dir.join("..").join("examples")] {
        if candidate.join("cu_fcc.cif").is_file() {
            return Some(candidate);
        }
    }
    loop {
        let candidate = dir.join("docs").join("manual").join("examples");
        if candidate.join("cu_fcc.cif").is_file() {
            return Some(candidate);
        }
        if !dir.pop() {
            return None;
        }
    }
}

fn file_name(path: &str) -> String {
    Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string())
}

/// Draws the panel centred in `area`.
pub fn show(ctx: &egui::Context, area: egui::Rect, recent: &[String]) -> Option<WelcomeAction> {
    let p = theme::palette(ctx.theme() == egui::Theme::Dark);
    let width = 560.0f32.min(area.width() - 40.0).max(200.0);
    let mut action = None;
    egui::Area::new(egui::Id::unique("welcome panel"))
        .order(egui::Order::Middle)
        .pivot(egui::Align2::CENTER_CENTER)
        .fixed_pos(egui::pos2(area.center().x, area.top() + area.height() * 0.52))
        .show(ctx, |ui| {
            egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(12).inner_margin(egui::Margin::symmetric(28, 24)).show(ui, |ui| {
                ui.set_width(width - 56.0);
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new("AtomForge").font(theme::heading(1.9)).color(p.text));
                    ui.label(egui::RichText::new(format!("Version {}", env!("CARGO_PKG_VERSION_MAJOR").to_string() + "." + env!("CARGO_PKG_VERSION_MINOR") + "." + env!("CARGO_PKG_VERSION_PATCH"))).color(p.text_muted));
                });
                ui.label(egui::RichText::new("Build, edit, analyse and render atomic structures.").color(p.text_muted));
                ui.add_space(12.0);
                let half = (ui.available_width() - ui.spacing().item_spacing.x) * 0.5;
                let big = egui::vec2(half, theme::FRAME_HEIGHT * 1.6);
                ui.horizontal(|ui| {
                    if theme::primary_button(ui, theme::with_icon(icon::OPEN, "Open Structure"), Some(big)).on_hover_text("CIF, POSCAR, XYZ, PDB, LAMMPS and more (Ctrl+O)").clicked() {
                        action = Some(WelcomeAction::OpenStructure);
                    }
                    if theme::secondary_button(ui, theme::with_icon(icon::CUBE, "Build Crystal"), Some(big)).on_hover_text("Build > Bulk Crystal: a crystal from its space group and lattice").clicked() {
                        action = Some(WelcomeAction::BuildCrystal);
                    }
                });
                ui.horizontal(|ui| {
                    if theme::secondary_button(ui, theme::with_icon(icon::FOLDER, "Open Project"), Some(big)).clicked() {
                        action = Some(WelcomeAction::OpenProject);
                    }
                    if theme::secondary_button(ui, theme::with_icon(icon::HELP, "Read Manual"), Some(big)).clicked() {
                        action = Some(WelcomeAction::ReadManual);
                    }
                });
                if let Some(examples) = examples_directory() {
                    theme::section(ui, "Examples");
                    ui.horizontal_wrapped(|ui| {
                        for e in &EXAMPLES {
                            let path = examples.join(e.file);
                            if !path.is_file() {
                                continue;
                            }
                            if ui.button(e.label).on_hover_text(format!("{}\n{}", e.about, e.file)).clicked() {
                                action = Some(WelcomeAction::Open(path.to_string_lossy().replace('\\', "/")));
                            }
                        }
                    });
                }
                if !recent.is_empty() {
                    theme::section(ui, "Recent files");
                    let font = egui::TextStyle::Body.resolve(ui.style());
                    for path in recent.iter().take(6) {
                        let name = file_name(path);
                        let missing = !Path::new(path).is_file();
                        ui.horizontal(|ui| {
                            let column = ui.available_width() * 0.42;
                            let (rect, response) = ui.allocate_exact_size(egui::vec2(column, theme::FRAME_HEIGHT), egui::Sense::click());
                            let response = if missing { response.on_hover_text(format!("{path}\n(the file is no longer there)")) } else { response.on_hover_text(path) };
                            if response.hovered() && !missing {
                                ui.painter().rect_filled(rect, 4.0, ui.visuals().selection.bg_fill.gamma_multiply(0.5));
                            }
                            ui.painter().text(rect.left_center(), egui::Align2::LEFT_CENTER, &name, font.clone(), if missing { p.text_muted } else { p.text });
                            if response.clicked() && !missing {
                                action = Some(WelcomeAction::Open(path.clone()));
                            }
                            // The folder, shortened from the left to fit.
                            let room = ui.available_width();
                            let mut folder = path[..path.len() - name.len().min(path.len())].to_string();
                            let fits = |ui: &egui::Ui, t: &str| ui.painter().layout_no_wrap(t.to_string(), font.clone(), p.text_muted).size().x <= room;
                            if !fits(ui, &folder) {
                                while folder.chars().count() > 4 && !fits(ui, &format!("…{folder}")) {
                                    folder.remove(0);
                                }
                                folder = format!("…{folder}");
                            }
                            ui.label(egui::RichText::new(folder).color(p.text_muted));
                        });
                    }
                }
                ui.add_space(4.0);
                ui.separator();
                ui.label(egui::RichText::new("Drop structure files here.").color(p.text_muted));
            });
        });
    action
}
