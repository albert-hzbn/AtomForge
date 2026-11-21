//! Help > Find a Command (Ctrl+K, the C++ ui/CommandPalette) and Help >
//! Keyboard Shortcuts (F1): every menu command by name, ranked as the C++
//! ranks them, and the list of key and mouse bindings.

use eframe::egui;

use crate::measure_ui::Kind;
use crate::view::AtomStyle;
use af_render::colouring::ColourMode;

/// What a command does (the app runs it).
#[derive(Clone, Debug, PartialEq)]
pub enum Action {
    Open,
    Save,
    SaveAs,
    ExportImage,
    Close,
    Quit,
    Undo,
    Redo,
    Step(&'static str),
    Tool(&'static str),
    ColourMode(ColourMode),
    Style(AtomStyle),
    Toggle(&'static str),
    Projection(bool),
    ResetView,
    ViewAxis(usize),
    ViewLattice(usize),
    StructureInfo,
    AtomInfo,
    Measure(Kind),
    Theme(bool),
    SystemTheme,
    Manual,
    InterfaceScale(f32),
    AtomicSizes,
    DisplaySettings,
    SelectAll,
    ClearSelection,
    Shortcuts,
    About,
}

#[derive(Clone, Debug)]
pub struct Command {
    pub label: String,
    pub path: String,
    pub action: Action,
    pub enabled: bool,
}

impl Command {
    pub fn new(path: &str, label: impl Into<String>, action: Action, enabled: bool) -> Self {
        Self { label: label.into(), path: path.to_string(), action, enabled }
    }
}

/// Every word of the query must occur in the label or the menu path;
/// matches at the start of the label, then of a word in it, then anywhere in
/// it, rank first. None: no match.
pub fn score(command: &Command, words: &[String]) -> Option<i32> {
    let label = command.label.to_lowercase();
    let path = command.path.to_lowercase();
    let mut total = 0;
    for word in words {
        match label.find(word.as_str()) {
            Some(0) => total += 30,
            Some(i) => total += if label.as_bytes()[i - 1] == b' ' { 20 } else { 10 },
            None if path.contains(word.as_str()) => total += 2,
            None => return None,
        }
    }
    Some(total + i32::from(command.enabled))
}

/// The matches of a query, best first (a stable sort keeps the menu order).
pub fn ranked<'a>(commands: &'a [Command], query: &str) -> Vec<&'a Command> {
    let words: Vec<String> = query.to_lowercase().split_whitespace().map(str::to_string).collect();
    let mut matches: Vec<(i32, &Command)> = commands.iter().filter_map(|c| score(c, &words).map(|s| (s, c))).collect();
    matches.sort_by(|a, b| b.0.cmp(&a.0));
    matches.into_iter().map(|(_, c)| c).collect()
}

#[derive(Default)]
pub struct Palette {
    pub open: bool,
    query: String,
    selected: usize,
    focus: bool,
}

impl Palette {
    /// Types a query (UI scripts).
    pub fn set_query(&mut self, query: &str) {
        self.query = query.to_string();
        self.selected = 0;
    }

    pub fn show_palette(&mut self) {
        self.open = true;
        self.query.clear();
        self.selected = 0;
        self.focus = true;
    }

    /// Draws the palette; returns the chosen command's action.
    pub fn show(&mut self, ctx: &egui::Context, commands: &[Command]) -> Option<Action> {
        if !self.open {
            return None;
        }
        let mut chosen = None;
        let screen = ctx.content_rect();
        let width = 620.0f32.min(screen.width() - 32.0);
        let mut close = false;
        egui::Area::new(egui::Id::unique("command palette")).order(egui::Order::Foreground).fixed_pos(egui::pos2(screen.center().x - width / 2.0, screen.min.y + screen.height() * 0.12)).show(ctx, |ui| {
            let p = crate::theme::palette(ui.visuals().dark_mode);
            egui::Frame::popup(ui.style()).fill(p.popup).inner_margin(egui::Margin::symmetric(12, 10)).show(ui, |ui| {
                ui.set_width(width - 24.0);
                let edit = ui.add_sized(
                    egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT + 4.0),
                    egui::TextEdit::singleline(&mut self.query).hint_text("Type a command, e.g. surface slab, rdf, dark theme").align(egui::Align2::LEFT_CENTER).margin(egui::vec2(10.0, 4.0)),
                );
                if self.focus {
                    edit.request_focus();
                    self.focus = false;
                }
                if edit.changed() {
                    self.selected = 0;
                }
                let matches = ranked(commands, &self.query);
                let (down, up, enter, escape) = ui.input(|i| (i.key_pressed(egui::Key::ArrowDown), i.key_pressed(egui::Key::ArrowUp), i.key_pressed(egui::Key::Enter), i.key_pressed(egui::Key::Escape)));
                if down {
                    self.selected = (self.selected + 1).min(matches.len().saturating_sub(1));
                }
                if up {
                    self.selected = self.selected.saturating_sub(1);
                }
                self.selected = self.selected.min(matches.len().saturating_sub(1));
                if matches.is_empty() {
                    ui.add_space(4.0);
                    ui.label(egui::RichText::new(format!("No command matches \"{}\".", self.query)).color(p.text_muted));
                } else {
                    // Rows a text line high, the chosen one highlighted across the width.
                    let row = crate::theme::BASE_FONT * 18.0 / 13.5 + 8.0;
                    let font = egui::TextStyle::Body.resolve(ui.style());
                    ui.spacing_mut().item_spacing.y = 0.0;
                    egui::ScrollArea::vertical().max_height(matches.len().min(10) as f32 * row).auto_shrink([false, true]).show(ui, |ui| {
                        for (i, c) in matches.iter().enumerate() {
                            let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width(), row), if c.enabled { egui::Sense::click() } else { egui::Sense::hover() });
                            let selected = i == self.selected;
                            if selected || (response.hovered() && c.enabled) {
                                let a = if selected { if ui.visuals().dark_mode { 61 } else { 36 } } else if ui.visuals().dark_mode { 87 } else { 56 };
                                ui.painter().rect_filled(rect, 4.0, egui::Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), a));
                            }
                            let colour = if c.enabled { p.text } else { p.text_muted };
                            ui.painter().text(egui::pos2(rect.left() + 8.0, rect.center().y), egui::Align2::LEFT_CENTER, &c.label, font.clone(), colour);
                            ui.painter().text(egui::pos2(rect.right() - 4.0, rect.center().y), egui::Align2::RIGHT_CENTER, &c.path, font.clone(), p.text_muted);
                            if selected && (up || down) {
                                response.scroll_to_me(None);
                            }
                            if response.clicked() {
                                chosen = Some(c.action.clone());
                            }
                        }
                    });
                    if enter && matches[self.selected].enabled {
                        chosen = Some(matches[self.selected].action.clone());
                    }
                }
                ui.add_space(10.0);
                ui.label(egui::RichText::new("Up/Down to choose, Enter to run, Esc to close").color(p.text_muted));
                if escape {
                    close = true;
                }
            });
        });
        if chosen.is_some() || close {
            self.open = false;
        }
        chosen
    }
}

/// Help > Keyboard Shortcuts (F1): every key and mouse binding of the main view (FileBrowser::drawShortcuts).
pub fn shortcuts_window(ctx: &egui::Context, open: &mut bool) {
    let groups: [(&str, &[(&str, &str)]); 5] = [
        ("Files", &[("Ctrl+O", "Open a structure"), ("Ctrl+S", "Save"), ("Ctrl+Shift+S", "Save as"), ("Ctrl+Alt+S", "Export an image"), ("Ctrl+W", "Close the structure")]),
        (
            "Editing",
            &[
                ("Ctrl+Z", "Undo"),
                ("Ctrl+Y, Ctrl+Shift+Z", "Redo"),
                ("Delete", "Delete the selected atoms"),
                ("S", "Substitute the selected atoms"),
                ("I", "Insert an atom at the midpoint of the selection"),
                ("G", "Move the selected atoms; X, Y or Z locks an axis, click or Enter confirms, Esc cancels"),
            ],
        ),
        (
            "Selection",
            &[("Hover", "Name the atom under the mouse (element, index, position)"), ("Click", "Select an atom"), ("Ctrl+click", "Add or remove an atom"), ("Ctrl+A", "Select all atoms"), ("Esc, Ctrl+D", "Clear the selection")],
        ),
        (
            "View",
            &[
                ("Left drag", "Orbit"),
                ("Right drag", "Pan (when no box or lasso selection is active)"),
                ("Wheel", "Zoom"),
                ("R", "Reset the view"),
                ("1, 2, 3", "Look along x, y, z"),
                ("Shift+1, 2, 3", "Look along the lattice vectors a, b, c"),
                ("F1", "This list"),
                ("Ctrl+K", "Find any command by name"),
            ],
        ),
        ("Dialogs", &[("Esc", "Close the dialog (in a field: cancel the edit)"), ("Tab, Shift+Tab", "Move between fields")]),
    ];
    let mut close = false;
    let height = (ctx.content_rect().height() - 120.0).max(300.0);
    crate::theme::modal_window("Keyboard Shortcuts").open(open).resizable(false).fixed_size([560.0, height]).show(ctx, |ui| {
        let p = crate::theme::palette(ui.visuals().dark_mode);
        let font = egui::TextStyle::Body.resolve(ui.style());
        egui::ScrollArea::vertical().max_height(height - 70.0).auto_shrink([false, true]).show(ui, |ui| {
            ui.spacing_mut().item_spacing.y = 0.0;
            for (title, rows) in groups {
                ui.add_space(8.0);
                crate::theme::section(ui, title);
                for (k, (keys, action)) in rows.iter().enumerate() {
                    // Two stretched columns (0.38 / 0.62); every other row shaded.
                    let width = ui.available_width();
                    let keys_width = width * 0.38;
                    let action_galley = ui.painter().layout(action.to_string(), font.clone(), p.text_muted, width - keys_width - 16.0);
                    let height = action_galley.size().y.max(18.0) + 10.0;
                    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, height), egui::Sense::hover());
                    if k % 2 == 1 {
                        let shade = if ui.visuals().dark_mode { egui::Color32::from_rgba_unmultiplied(255, 255, 255, 8) } else { egui::Color32::from_rgba_unmultiplied(0, 0, 0, 6) };
                        ui.painter().rect_filled(rect, 0.0, shade);
                    }
                    ui.painter().text(egui::pos2(rect.left() + 8.0, rect.top() + 5.0), egui::Align2::LEFT_TOP, *keys, font.clone(), p.text);
                    ui.painter().galley(egui::pos2(rect.left() + keys_width + 8.0, rect.top() + 5.0), action_galley, p.text_muted);
                }
            }
        });
        ui.add_space(6.0);
        if crate::dialog_layout::footer(ui, "Close", None, true) == 1 {
            close = true;
        }
    });
    if close {
        *open = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(path: &str, label: &str, enabled: bool) -> Command {
        Command::new(path, label, Action::About, enabled)
    }

    #[test]
    fn ranking_follows_the_cpp_scores() {
        let commands = vec![
            c("Build", "Bulk Crystal", true),
            c("Build > Low-Dimensional", "Surface Slab", true),
            c("Analysis > Structure and Defects", "Radial Distribution Function", false),
            c("View > Theme", "Dark theme", true),
            c("Edit", "Slab Thickness", true),
        ];
        let words = |q: &str| q.split_whitespace().map(str::to_string).collect::<Vec<_>>();
        assert_eq!(score(&commands[1], &words("slab")), Some(21));
        assert_eq!(score(&commands[4], &words("slab")), Some(31));
        assert_eq!(score(&commands[1], &words("surface slab")), Some(51));
        assert_eq!(score(&commands[1], &words("low")), Some(3));
        assert_eq!(score(&commands[2], &words("distribution")), Some(20));
        assert_eq!(score(&commands[2], &words("xyz")), None);
        assert_eq!(score(&commands[3], &words("eme")), Some(11));
        let ranked = ranked(&commands, "Slab");
        assert_eq!(ranked.iter().map(|c| c.label.as_str()).collect::<Vec<_>>(), ["Slab Thickness", "Surface Slab"]);
        assert_eq!(super::ranked(&commands, "").len(), 5);
        assert_eq!(super::ranked(&commands, "")[0].label, "Bulk Crystal");
    }
}
