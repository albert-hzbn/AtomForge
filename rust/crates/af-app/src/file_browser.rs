//! A file browser window for opening and saving (the C++ Open Structure,
//! Save As and Export Image popups of ui/FileBrowser and ui/PathPicker):
//! back, forward and up with an editable path, a sidebar of locations and
//! drives beside the folder's listing (folder icons), the file name, the
//! format (Save As, Export Image), the export options and a footer.

use std::path::{Path, PathBuf};

use eframe::egui;

use crate::dialog_layout::{self, Note};
use crate::theme::FRAME_HEIGHT;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Purpose {
    Open,
    SaveAs,
    ExportImage,
    OpenProject,
    SaveProject,
    SavePipeline,
    LoadPipeline,
    OpenTrajectory,
    ExportGif,
}

/// Save As formats (label, extension), as the C++ kSaveFormats.
pub const SAVE_FORMATS: [(&str, &str); 13] = [
    ("XYZ (.xyz)", "xyz"),
    ("CIF (.cif)", "cif"),
    ("VASP POSCAR (.vasp)", "vasp"),
    ("LAMMPS Data (.data)", "data"),
    ("PDB (.pdb)", "pdb"),
    ("SDF (.sdf)", "sdf"),
    ("Mol2 (.mol2)", "mol2"),
    ("Quantum ESPRESSO (.pwi)", "pwi"),
    ("CASTEP Cell (.cell)", "cell"),
    ("DFTB+ (.gen)", "gen"),
    ("XCrySDen (.xsf)", "xsf"),
    ("CP2K SUBSYS (.cp2k)", "cp2k"),
    ("Gaussian Input (.gjf)", "gjf"),
];

/// Export Image formats (label, extension), as the C++ kImageExportFormats.
pub const IMAGE_FORMATS: [(&str, &str); 3] = [("PNG (.png)", "png"), ("JPEG (.jpg)", "jpg"), ("SVG (.svg)", "svg")];

pub struct FileBrowser {
    pub open: bool,
    pub purpose: Purpose,
    pub title: String,
    pub folder: PathBuf,
    pub name: String,
    /// Lower-case extensions shown (empty: all files).
    pub extensions: Vec<String>,
    pub error: String,
    /// Export Image: the background, the orientation gizmo and the resolution
    /// multiple (1-8) the image is drawn with.
    pub export_background: bool,
    pub export_gizmo: bool,
    pub export_scale: i32,
    /// The view's size in pixels (for the resolution text).
    pub view_size: (i32, i32),
    entries: Vec<(String, bool)>,
    listed: Option<PathBuf>,
    history: Vec<PathBuf>,
    history_index: usize,
    path_text: String,
    sidebar_width: f32,
}

impl Default for FileBrowser {
    fn default() -> Self {
        Self {
            open: false,
            purpose: Purpose::Open,
            title: String::new(),
            folder: std::env::current_dir().unwrap_or_default(),
            name: String::new(),
            extensions: Vec::new(),
            error: String::new(),
            export_background: true,
            export_gizmo: false,
            export_scale: 1,
            view_size: (1920, 1080),
            entries: Vec::new(),
            listed: None,
            history: Vec::new(),
            history_index: 0,
            path_text: String::new(),
            sidebar_width: 150.0,
        }
    }
}

/// A path with forward slashes, as the C++ shows them.
fn display(path: &Path) -> String {
    path.to_string_lossy().replace('\\', "/")
}

fn home() -> PathBuf {
    std::env::var_os(if cfg!(windows) { "USERPROFILE" } else { "HOME" }).map(PathBuf::from).unwrap_or_else(|| PathBuf::from("."))
}

/// The drives (Windows letters that exist; "/" elsewhere) and the home folder.
fn drive_roots() -> Vec<PathBuf> {
    let mut out = Vec::new();
    if cfg!(windows) {
        for letter in b'A'..=b'Z' {
            let root = PathBuf::from(format!("{}:/", letter as char));
            if root.is_dir() {
                out.push(root);
            }
        }
    } else {
        out.push(PathBuf::from("/"));
    }
    out.push(home());
    out
}

/// "CIF structure, 42.1 KB, modified 2026-09-09 14:02" (describeFile), "" if not a file.
fn describe_file(path: &Path) -> String {
    let Ok(meta) = std::fs::metadata(path) else { return String::new() };
    if !meta.is_file() {
        return String::new();
    }
    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let ext = path.extension().map(|e| e.to_string_lossy().to_ascii_lowercase()).unwrap_or_default();
    let kind = if ext == "afcloud" {
        "Large dataset (opens in the large-data viewer)".to_string()
    } else if name.starts_with("POSCAR") || name.starts_with("CONTCAR") || ext == "vasp" {
        "VASP structure".to_string()
    } else {
        format!("{} structure", if ext.is_empty() { name.to_ascii_uppercase() } else { ext.to_ascii_uppercase() })
    };
    let size = meta.len() as f64;
    let size_text = if size < 1024.0 {
        format!("{size:.0} bytes")
    } else if size < 1024.0 * 1024.0 {
        format!("{:.1} KB", size / 1024.0)
    } else if size < 1024.0 * 1024.0 * 1024.0 {
        format!("{:.1} MB", size / (1024.0 * 1024.0))
    } else {
        format!("{:.2} GB", size / (1024.0 * 1024.0 * 1024.0))
    };
    let mut text = format!("{kind}, {size_text}");
    if let Ok(modified) = meta.modified() {
        if let Ok(since) = modified.duration_since(std::time::UNIX_EPOCH) {
            // Local time is not available without a time-zone library: UTC.
            let secs = since.as_secs() as i64;
            let (days, rest) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
            let (y, m, d) = civil_from_days(days);
            text += &format!(", modified {y:04}-{m:02}-{d:02} {:02}:{:02}", rest / 3600, rest % 3600 / 60);
        }
    }
    text
}

/// The calendar date of a day count since 1970-01-01 (Howard Hinnant's algorithm).
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (yoe + era * 400 + if m <= 2 { 1 } else { 0 }, m, d)
}

/// `name` with its extension replaced (updateFilenameWithExtension).
fn with_extension(name: &str, ext: &str) -> String {
    let stem = match name.rfind('.') {
        Some(dot) if dot > 0 => &name[..dot],
        _ if name.is_empty() => "structure",
        _ => name,
    };
    format!("{stem}.{ext}")
}

impl FileBrowser {
    pub fn show_for(&mut self, purpose: Purpose, title: &str, name: &str, extensions: &[&str]) {
        self.open = true;
        self.purpose = purpose;
        // The C++ window titles of the three main file dialogs.
        self.title = match purpose {
            Purpose::Open => "Open Structure".to_string(),
            Purpose::SaveAs => "Save As".to_string(),
            Purpose::ExportImage => "Export Image".to_string(),
            _ => title.to_string(),
        };
        self.name = name.to_string();
        self.extensions = extensions.iter().map(|e| e.to_ascii_lowercase()).collect();
        self.error.clear();
        self.listed = None;
        if self.history.is_empty() {
            self.history.push(self.folder.clone());
            self.history_index = 0;
        }
    }

    /// Goes to a folder and records it in the history.
    fn navigate(&mut self, folder: PathBuf) {
        self.history.truncate(self.history_index + 1);
        self.history.push(folder.clone());
        self.history_index = self.history.len() - 1;
        self.folder = folder;
        self.error.clear();
    }

    fn opening(&self) -> bool {
        matches!(self.purpose, Purpose::Open | Purpose::OpenProject | Purpose::LoadPipeline | Purpose::OpenTrajectory)
    }

    fn refresh(&mut self) {
        if self.listed.as_ref() == Some(&self.folder) {
            return;
        }
        self.entries.clear();
        if let Ok(read) = std::fs::read_dir(&self.folder) {
            for entry in read.flatten() {
                let name = entry.file_name().to_string_lossy().into_owned();
                let is_dir = entry.path().is_dir();
                if !is_dir && !self.extensions.is_empty() {
                    let ext = Path::new(&name).extension().map(|e| e.to_string_lossy().to_ascii_lowercase()).unwrap_or_default();
                    let upper = name.to_ascii_uppercase();
                    let vasp = upper.contains("POSCAR") || upper.contains("CONTCAR");
                    if !self.extensions.contains(&ext) && !vasp {
                        continue;
                    }
                }
                self.entries.push((name, is_dir));
            }
        }
        self.entries.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.to_lowercase().cmp(&b.0.to_lowercase())));
        self.listed = Some(self.folder.clone());
    }

    /// Draws the window; returns the chosen path when the user confirms.
    pub fn show(&mut self, ctx: &egui::Context) -> Option<PathBuf> {
        if !self.open {
            return None;
        }
        self.refresh();
        let mut chosen = None;
        let mut open = self.open;
        let size = if matches!(self.purpose, Purpose::Open) { [720.0, 510.0] } else { [720.0, 560.0] };
        crate::theme::modal_window(self.title.clone()).id(egui::Id::unique("file browser")).open(&mut open).default_size(size).min_size([520.0, 400.0]).show(ctx, |ui| {
            chosen = self.contents(ui);
        });
        self.open = open && chosen.is_none();
        if let Some(path) = &chosen {
            if let Some(parent) = path.parent() {
                self.folder = parent.to_path_buf();
            }
        }
        chosen
    }

    fn contents(&mut self, ui: &mut egui::Ui) -> Option<PathBuf> {
        let mut chosen = None;
        let p = crate::theme::palette(ui.visuals().dark_mode);
        // --- The navigation bar with an editable path. ---
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 2.0;
            let nav = egui::vec2(32.0, FRAME_HEIGHT);
            if ui.add_enabled(self.history_index > 0, egui::Button::new("←").min_size(nav)).on_hover_text("Back").clicked() {
                self.history_index -= 1;
                self.folder = self.history[self.history_index].clone();
            }
            if ui.add_enabled(self.history_index + 1 < self.history.len(), egui::Button::new("→").min_size(nav)).on_hover_text("Forward").clicked() {
                self.history_index += 1;
                self.folder = self.history[self.history_index].clone();
            }
            if ui.add(egui::Button::new("↑").min_size(nav)).on_hover_text("Up one level").clicked() {
                if let Some(parent) = self.folder.parent().map(Path::to_path_buf) {
                    self.navigate(parent);
                }
            }
            ui.add_space(6.0);
            let shown = display(&self.folder);
            let id = ui.make_persistent_id("file browser path");
            let editing = ui.memory(|m| m.has_focus(id));
            if !editing {
                self.path_text = shown;
            }
            let response = ui.add_sized(egui::vec2(ui.available_width(), FRAME_HEIGHT), egui::TextEdit::singleline(&mut self.path_text).id(id).align(egui::Align2::LEFT_CENTER).margin(egui::vec2(8.0, 4.0)));
            if response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                let mut text = self.path_text.trim().replace('\\', "/");
                if text.len() == 2 && text.ends_with(':') {
                    text.push('/');
                }
                let folder = PathBuf::from(&text);
                if folder.is_dir() {
                    self.navigate(folder);
                } else {
                    self.error = format!("{text} is not a folder");
                }
            }
        });
        ui.separator();
        self.refresh();
        // --- The sidebar and the listing take the height the fields below leave. ---
        let below = match self.purpose {
            Purpose::ExportImage => 5.0 * (FRAME_HEIGHT + 8.0) + 60.0,
            Purpose::SaveAs => 2.0 * (FRAME_HEIGHT + 8.0) + 60.0,
            _ => FRAME_HEIGHT + 8.0 + 84.0,
        };
        let list_height = (ui.available_height() - below).max(300.0);
        let mut go_to: Option<PathBuf> = None;
        let mut double_clicked = false;
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 0.0;
            let sidebar = self.sidebar_width;
            ui.allocate_ui_with_layout(egui::vec2(sidebar, list_height), egui::Layout::top_down(egui::Align::Min), |ui| {
                list_frame(ui, list_height, |ui| {
                    if let Some(place) = locations(ui, &p) {
                        go_to = Some(place);
                    }
                });
            });
            // The splitter between them.
            let (rect, response) = ui.allocate_exact_size(egui::vec2(6.0, list_height), egui::Sense::drag());
            if response.hovered() || response.dragged() {
                ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeHorizontal);
            }
            if response.dragged() {
                self.sidebar_width = (self.sidebar_width + response.drag_delta().x).clamp(80.0, 300.0);
            }
            let hot = response.hovered() || response.dragged();
            ui.painter().vline(rect.center().x, (rect.top() + 4.0)..=(rect.bottom() - 4.0), egui::Stroke::new(1.5, if hot { egui::Color32::from_rgba_unmultiplied(150, 190, 230, 200) } else { egui::Color32::from_rgba_unmultiplied(120, 120, 120, 100) }));
            let width = ui.available_width();
            ui.allocate_ui_with_layout(egui::vec2(width, list_height), egui::Layout::top_down(egui::Align::Min), |ui| {
                list_frame(ui, list_height, |ui| {
                    if std::fs::read_dir(&self.folder).is_err() {
                        ui.label(egui::RichText::new("Unable to open folder").color(p.text_muted));
                        return;
                    }
                    let font = egui::TextStyle::Body.resolve(ui.style());
                    let row = crate::theme::BASE_FONT * 18.0 / 13.5 + 8.0;
                    ui.spacing_mut().item_spacing.y = 0.0;
                    for (name, is_dir) in &self.entries {
                        let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width(), row), egui::Sense::click());
                        let selected = !*is_dir && self.name == *name;
                        if selected || response.hovered() {
                            let a = if selected { if ui.visuals().dark_mode { 61 } else { 36 } } else if ui.visuals().dark_mode { 87 } else { 56 };
                            ui.painter().rect_filled(rect, 3.0, egui::Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), a));
                        }
                        let th = crate::theme::BASE_FONT * 18.0 / 13.5;
                        if *is_dir {
                            // A two-rectangle folder icon before the name.
                            let iw = th * 1.1;
                            let (x0, y0) = (rect.left() + 2.0, rect.center().y - th * 0.5);
                            let col = egui::Color32::from_rgba_unmultiplied(100, 175, 255, 225);
                            ui.painter().rect_filled(egui::Rect::from_min_max(egui::pos2(x0, y0 + th * 0.14), egui::pos2(x0 + iw * 0.48, y0 + th * 0.36)), 2.0, col);
                            ui.painter().rect_filled(egui::Rect::from_min_max(egui::pos2(x0, y0 + th * 0.28), egui::pos2(x0 + iw, y0 + th * 0.92)), 2.0, col);
                            let folder_text = if ui.visuals().dark_mode { egui::Color32::from_rgb(140, 199, 255) } else { egui::Color32::from_rgb(41, 56, 77) };
                            ui.painter().text(egui::pos2(x0 + iw + 11.0, rect.center().y), egui::Align2::LEFT_CENTER, name, font.clone(), folder_text);
                            if response.clicked() {
                                go_to = Some(self.folder.join(name));
                            }
                        } else {
                            ui.painter().text(egui::pos2(rect.left() + 8.0, rect.center().y), egui::Align2::LEFT_CENTER, name, font.clone(), p.text);
                            if response.clicked() {
                                self.name = name.clone();
                                self.error.clear();
                            }
                            if response.double_clicked() {
                                self.name = name.clone();
                                double_clicked = true;
                            }
                        }
                    }
                });
            });
        });
        if let Some(folder) = go_to {
            self.navigate(folder);
        }
        // --- The file name, the format and the export options. ---
        let label_width = if self.purpose == Purpose::Open { 66.0 } else { 119.0 };
        let hint = if self.purpose == Purpose::Open { "CIF, POSCAR/.vasp, XYZ, PDB, SDF, MOL2, PWI, GJF, LAMMPS dump or .afcloud" } else { "" };
        let mut entered = false;
        dialog_layout::row(ui, "File name", label_width, |ui| {
            let response = ui.add_sized(egui::vec2(ui.available_width(), FRAME_HEIGHT), egui::TextEdit::singleline(&mut self.name).hint_text(hint).align(egui::Align2::LEFT_CENTER).margin(egui::vec2(8.0, 4.0)));
            entered = response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        });
        match self.purpose {
            Purpose::SaveAs => self.format_row(ui, label_width, &SAVE_FORMATS, "save format"),
            Purpose::ExportImage => {
                self.format_row(ui, label_width, &IMAGE_FORMATS, "image format");
                dialog_layout::labelled_check(ui, &mut self.export_background, "Include background");
                if self.name.to_ascii_lowercase().ends_with(".jpg") && !self.export_background {
                    ui.label(egui::RichText::new("JPEG has no transparency; transparent areas will be blended on white.").color(p.text_muted));
                }
                dialog_layout::labelled_check(ui, &mut self.export_gizmo, "Include orientation gizmo");
                let (w, h) = (self.view_size.0 * self.export_scale, self.view_size.1 * self.export_scale);
                let (vw, vh) = self.view_size;
                dialog_layout::row(ui, "Resolution", label_width, |ui| {
                    let width = ui.available_width();
                    dialog_layout::slider_i32(ui, &mut self.export_scale, 1..=8, width, |s| format!("{s}x   ({} × {} px)", vw * s, vh * s));
                });
                let _ = (w, h);
            }
            _ => {}
        }
        // The selected file: format, size and date.
        if self.opening() && !self.name.is_empty() {
            let details = describe_file(&self.folder.join(&self.name));
            if !details.is_empty() {
                ui.label(egui::RichText::new(details).color(p.text_muted));
            }
        }
        if !self.error.is_empty() {
            dialog_layout::note(ui, Note::Error, &self.error.clone());
        }
        let verb = match self.purpose {
            Purpose::ExportImage | Purpose::ExportGif => "Export",
            _ if self.opening() => "Open",
            _ => "Save",
        };
        let mut confirm = double_clicked || (entered && !self.name.is_empty());
        match dialog_layout::footer(ui, verb, Some("Cancel"), !self.name.is_empty()) {
            1 => confirm = true,
            2 => self.open = false,
            _ => {}
        }
        if confirm && !self.name.is_empty() {
            let path = if Path::new(&self.name).is_absolute() { PathBuf::from(&self.name) } else { self.folder.join(&self.name) };
            if self.opening() && !path.is_file() {
                self.error = format!("{} does not exist", display(&path));
            } else {
                chosen = Some(path);
            }
        }
        chosen
    }

    /// The Format row: choosing a format gives the file name its extension.
    fn format_row(&mut self, ui: &mut egui::Ui, label_width: f32, formats: &[(&str, &str)], salt: &str) {
        let ext = Path::new(&self.name).extension().map(|e| e.to_string_lossy().to_ascii_lowercase()).unwrap_or_default();
        let current = formats.iter().position(|(_, e)| *e == ext).unwrap_or(0);
        let mut chosen = current;
        dialog_layout::row(ui, "Format", label_width, |ui| {
            let width = ui.available_width();
            egui::ComboBox::from_id_salt(salt).width(width - 8.0).selected_text(formats[current].0).show_ui(ui, |ui| {
                for (k, (label, _)) in formats.iter().enumerate() {
                    ui.selectable_value(&mut chosen, k, *label);
                }
            });
        });
        if chosen != current || ext.is_empty() {
            self.name = with_extension(&self.name, formats[chosen].1);
            self.error.clear();
        }
    }
}

/// A bordered list panel of a fixed height that scrolls (the C++ child windows).
fn list_frame(ui: &mut egui::Ui, height: f32, add: impl FnOnce(&mut egui::Ui)) {
    let p = crate::theme::palette(ui.visuals().dark_mode);
    egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(8, 8)).show(ui, |ui| {
        ui.set_width(ui.available_width());
        egui::ScrollArea::vertical().id_salt(ui.next_auto_id()).max_height(height - 18.0).min_scrolled_height(height - 18.0).auto_shrink([false, false]).show(ui, add);
    });
}

/// The sidebar: common folders, the bundled examples and the drives (drawLocations).
fn locations(ui: &mut egui::Ui, p: &crate::theme::Palette) -> Option<PathBuf> {
    let mut chosen = None;
    let font = egui::TextStyle::Body.resolve(ui.style());
    let row = crate::theme::BASE_FONT * 18.0 / 13.5 + 8.0;
    ui.spacing_mut().item_spacing.y = 0.0;
    let header = |ui: &mut egui::Ui, text: &str| {
        ui.add_space(4.0);
        ui.label(egui::RichText::new(text).color(p.text_muted));
        ui.add_space(4.0);
        let r = ui.available_rect_before_wrap();
        ui.painter().hline(r.x_range(), r.top(), egui::Stroke::new(1.0, p.separator));
        ui.add_space(6.0);
    };
    let mut item = |ui: &mut egui::Ui, label: &str, path: PathBuf, tip: Option<&str>| {
        let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width(), row), egui::Sense::click());
        if response.hovered() {
            ui.painter().rect_filled(rect, 3.0, egui::Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), if ui.visuals().dark_mode { 87 } else { 56 }));
        }
        ui.painter().text(egui::pos2(rect.left() + 4.0, rect.center().y), egui::Align2::LEFT_CENTER, label, font.clone(), p.text);
        let response = match tip {
            Some(t) => response.on_hover_text(t),
            None => response,
        };
        if response.clicked() {
            chosen = Some(path);
        }
    };
    header(ui, "Locations");
    let home = home();
    item(ui, "Home", home.clone(), None);
    for name in ["Desktop", "Documents", "Downloads"] {
        let folder = home.join(name);
        if folder.is_dir() {
            item(ui, name, folder, None);
        }
    }
    if let Some(examples) = crate::welcome::examples_directory() {
        item(ui, "Examples", examples, Some("The example structures shipped with AtomForge"));
    }
    ui.add_space(4.0);
    header(ui, "Drives");
    for root in drive_roots() {
        let label = display(&root);
        item(ui, &label, root.clone(), Some(&label));
    }
    chosen
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers() {
        assert_eq!(with_extension("cu.cif", "xyz"), "cu.xyz");
        assert_eq!(with_extension("", "png"), "structure.png");
        assert_eq!(with_extension("noext", "svg"), "noext.svg");
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(20_000), (2024, 10, 4));
    }
}
