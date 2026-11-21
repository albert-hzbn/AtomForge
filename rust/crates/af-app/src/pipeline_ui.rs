//! Edit > Structure Pipeline (the C++ ui/PipelineDialog): non-destructive
//! modifiers on a copy of the structure. The steps in order with their
//! results, reordering (drag, arrows, Alt+Up/Down), bypassing, replacing,
//! inserting and removing; the selected step's parameters; the pipeline as a
//! command; Bake and Discard; pipeline files.

use af_core::json::Json;
use af_core::structure::Structure;
use af_pipeline::editor::PipelineEditor;
use af_pipeline::pipeline::{find_modifier_type, modifier_types, Pipeline};
use eframe::egui;

use crate::dialog_layout;
use crate::params_form::{parameter_form, FormTexts};

/// What the panel asks the app for.
#[derive(Default)]
pub struct Requests {
    /// The structure shown changed (pipeline output; not an undo step).
    pub output: Option<Structure>,
    /// Bake: the output as an undoable edit.
    pub bake: Option<Structure>,
    pub save_file: bool,
    pub load_file: bool,
    /// Colour by the pipeline's selection (Atom Property).
    pub show_property: bool,
    /// "Edit in Dialog": open this step (0-based) in its Build or Edit dialog.
    pub edit_in_dialog: Option<usize>,
}

pub struct PipelinePanel {
    pub open: bool,
    pub selected: Option<usize>,
    pub auto_update: bool,
    /// Docked at the right of the view instead of floating.
    pub docked: bool,
    pub colour_selection: bool,
    pub message: String,
    show_text: bool,
    text: String,
    texts: FormTexts,
    /// The step whose form texts are cached.
    texts_for: Option<(usize, String)>,
    /// A step to open in its dialog (handed to the app with the next Requests).
    pub edit_request: Option<usize>,
}

impl Default for PipelinePanel {
    fn default() -> Self {
        Self { open: false, selected: None, auto_update: true, docked: false, colour_selection: false, message: String::new(), show_text: false, text: String::new(), texts: FormTexts::default(), texts_for: None, edit_request: None }
    }
}

fn title_of(type_id: &str, label: &str) -> String {
    let title = find_modifier_type(type_id).map(|t| t.title.clone()).unwrap_or_else(|| type_id.to_string());
    if label.is_empty() {
        title
    } else {
        format!("{label}  ({title})")
    }
}

/// Every modifier by category; returns the type chosen.
fn modifier_menu(ui: &mut egui::Ui, exclude: &str) -> Option<String> {
    let mut chosen = None;
    let types = modifier_types();
    let mut categories: Vec<String> = Vec::new();
    for t in &types {
        if !categories.contains(&t.category) {
            categories.push(t.category.clone());
        }
    }
    for category in categories {
        ui.menu_button(&category, |ui| {
            for t in types.iter().filter(|t| t.category == category) {
                let help = t.help.split('\n').next().unwrap_or("").to_string();
                if ui.add_enabled(t.id != exclude, egui::Button::new(&t.title)).on_hover_text(help).clicked() {
                    chosen = Some(t.id.clone());
                    ui.close();
                }
            }
        });
    }
    chosen
}

impl PipelinePanel {
    /// The output after an edit (also while the panel is closed): the
    /// structure to show, with the selection as a property when asked.
    pub fn update(&mut self, editor: &mut PipelineEditor, force: bool) -> Requests {
        let mut r = Requests::default();
        if !editor.active() || !(force || (editor.dirty() && self.auto_update)) {
            return r;
        }
        match editor.evaluate() {
            Ok(out) => {
                let mut s = out.structure.clone();
                if self.colour_selection && !s.atoms.is_empty() {
                    s.atom_property = out.selected.iter().map(|&v| v as f64).collect();
                    s.atom_property_name = "Pipeline selection (1 = selected)".into();
                    r.show_property = true;
                }
                r.output = Some(s);
            }
            Err(e) => self.message = e,
        }
        r
    }

    fn select(&mut self, index: Option<usize>) {
        self.selected = index;
    }

    /// Selects a step after Update Step; its form shows the new parameters.
    pub fn select_step(&mut self, index: usize) {
        self.select(Some(index));
        self.texts_for = None;
        self.texts.clear();
    }

    /// Moves a step; the selection follows it (and the steps it passed).
    pub fn move_step(&mut self, editor: &mut PipelineEditor, from: usize, to: usize) {
        if from == to || editor.move_step(from, to).is_err() {
            return;
        }
        self.selected = self.selected.map(|s| {
            if s == from {
                to
            } else if from < to && s > from && s <= to {
                s - 1
            } else if from > to && s >= to && s < from {
                s + 1
            } else {
                s
            }
        });
    }

    pub fn remove_step(&mut self, editor: &mut PipelineEditor, index: usize) {
        if editor.remove(index).is_err() {
            return;
        }
        self.selected = match self.selected {
            Some(s) if s > index => Some(s - 1),
            Some(s) if s >= editor.size() => editor.size().checked_sub(1),
            other => other,
        };
    }

    /// Adds a step after the selected one (else at the end), starting the pipeline when needed.
    /// Steps appended after the last (PipelineDialog::addSteps); the panel opens on the last.
    pub fn add_steps(&mut self, editor: &mut PipelineEditor, structure: &Structure, steps: Vec<af_pipeline::pipeline::Modifier>) {
        if !editor.active() && !structure.atoms.is_empty() {
            editor.set_input(structure);
        }
        let mut pipeline = editor.pipeline().clone();
        pipeline.modifiers.extend(steps);
        editor.set_pipeline(pipeline);
        self.select(editor.size().checked_sub(1));
        self.open = true;
    }

    pub fn add_step(&mut self, editor: &mut PipelineEditor, structure: &Structure, type_id: &str, at: Option<usize>) {
        if !editor.active() && !structure.atoms.is_empty() {
            editor.set_input(structure);
        }
        let at = at.unwrap_or_else(|| self.selected.map(|s| s + 1).unwrap_or(editor.size()));
        match editor.add(type_id, at) {
            Ok(i) => self.select(Some(i)),
            Err(e) => self.message = e,
        }
    }

    /// The floating window (when not docked at the right of the view).
    pub fn show(&mut self, ctx: &egui::Context, editor: &mut PipelineEditor, structure: &Structure) -> Requests {
        let mut requests = Requests::default();
        if !self.open || self.docked {
            return requests;
        }
        let mut open = self.open;
        crate::theme::window("Structure Pipeline").id(egui::Id::unique("structure pipeline")).open(&mut open).default_size([560.0, 760.0]).min_size([420.0, 480.0]).show(ctx, |ui| {
            self.contents(ui, editor, structure, &mut requests);
        });
        self.open = open;
        requests
    }

    /// Docked at the right of the view (View > Dock Structure Pipeline at the Right):
    /// the panel's tab (as an ImGui dock tab) with its close cross, and the
    /// close button of the dock node at the right; right-click the tab to undock.
    pub fn show_docked(&mut self, ui: &mut egui::Ui, editor: &mut PipelineEditor, structure: &Structure) -> Requests {
        let mut requests = Requests::default();
        let p = crate::theme::palette(ui.visuals().dark_mode);
        let font = egui::TextStyle::Body.resolve(ui.style());
        ui.horizontal(|ui| {
            let galley = ui.painter().layout_no_wrap("Structure Pipeline".to_string(), font.clone(), p.text);
            let size = egui::vec2(galley.size().x + 44.0, crate::theme::FRAME_HEIGHT - 4.0);
            let (tab, response) = ui.allocate_exact_size(size, egui::Sense::click());
            ui.painter().rect_filled(tab, egui::CornerRadius { nw: 6, ne: 6, sw: 0, se: 0 }, p.surface);
            ui.painter().line_segment([tab.left_top() + egui::vec2(3.0, 1.0), tab.right_top() + egui::vec2(-3.0, 1.0)], egui::Stroke::new(2.0, p.accent));
            ui.painter().galley(egui::pos2(tab.left() + 9.0, tab.center().y - galley.size().y / 2.0), galley, p.text);
            let cross = egui::Rect::from_center_size(egui::pos2(tab.right() - 14.0, tab.center().y), egui::vec2(16.0, 16.0));
            let cross_response = ui.interact(cross, ui.scope_id().with("dock tab close"), egui::Sense::click());
            if cross_response.hovered() {
                ui.painter().rect_filled(cross, 3.0, p.button_hover);
            }
            let (c, d) = (cross.center(), 4.0);
            let stroke = egui::Stroke::new(1.2, p.text);
            ui.painter().line_segment([c - egui::vec2(d, d), c + egui::vec2(d, d)], stroke);
            ui.painter().line_segment([c + egui::vec2(-d, d), c + egui::vec2(d, -d)], stroke);
            if cross_response.clicked() {
                self.open = false;
            }
            response.context_menu(|ui| {
                if ui.button("Undock").clicked() {
                    self.docked = false;
                    ui.close();
                }
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let (rect, response) = ui.allocate_exact_size(egui::vec2(18.0, 18.0), egui::Sense::click());
                if response.hovered() {
                    ui.painter().rect_filled(rect, 3.0, p.button_hover);
                }
                let (c, d) = (rect.center(), 5.0);
                ui.painter().line_segment([c - egui::vec2(d, d), c + egui::vec2(d, d)], stroke);
                ui.painter().line_segment([c + egui::vec2(-d, d), c + egui::vec2(d, -d)], stroke);
                if response.on_hover_text("Close").clicked() {
                    self.open = false;
                }
            });
        });
        ui.add_space(4.0);
        self.contents(ui, editor, structure, &mut requests);
        requests
    }

    fn contents(&mut self, ui: &mut egui::Ui, editor: &mut PipelineEditor, structure: &Structure, requests: &mut Requests) {
        let muted = crate::theme::palette(ui.visuals().dark_mode).text_muted;
        ui.add(egui::Label::new(egui::RichText::new("Modifiers change a copy of the structure in order; the input is kept, so steps can be edited, bypassed, reordered or removed at any time.").color(muted)).wrap());
        if !editor.active() {
            if ui.add_enabled(!structure.atoms.is_empty(), egui::Button::new("Start from the Active Structure").min_size(egui::vec2(ui.available_width(), 34.0))).clicked() {
                editor.set_input(structure);
                self.message.clear();
            }
            ui.label(egui::RichText::new(if structure.atoms.is_empty() { "Load or build a structure first." } else { "Adding a modifier also starts the pipeline from the active structure." }).color(muted));
        } else {
            // Edits made outside the pipeline would be overwritten by its next update.
            if !editor.dirty() && !editor.is_current_output(structure) && !self.colour_selection {
                let red = dialog_layout::status_error(ui);
                ui.add(egui::Label::new(egui::RichText::new("The structure was changed outside the pipeline; its next update would replace those changes.").color(red)).wrap());
                ui.horizontal_wrapped(|ui| {
                    if ui.button("Use It as the Pipeline Input").clicked() {
                        editor.set_input(structure);
                    }
                    if ui.button("Restore the Pipeline Output").clicked() {
                        *requests = self.update(editor, true);
                    }
                });
            }
            ui.horizontal_wrapped(|ui| {
                if dialog_layout::accent_button(ui, "Bake", true).on_hover_text("Make the output the structure and end the pipeline.").clicked() {
                    requests.bake = Some(editor.bake());
                    self.select(None);
                }
                if dialog_layout::danger_button(ui, "Discard").on_hover_text("End the pipeline and show the input structure again.").clicked() {
                    let input = editor.input().cloned().unwrap_or_default();
                    editor.reset();
                    self.select(None);
                    requests.output = Some(input);
                }
                if ui.button("Use the Active Structure as Input").on_hover_text("Restart the pipeline from the structure now in the view.").clicked() {
                    editor.set_input(structure);
                }
            });
        }
        ui.separator();
        // The steps: exactly as tall as the input, the steps and the hint,
        // between four rows and 45% of the panel; it scrolls beyond that.
        let row = crate::theme::FRAME_HEIGHT + ui.spacing().item_spacing.y;
        let maximum = (4.0 * row).max(ui.available_height() * 0.45);
        dialog_layout::child_frame(ui, |ui| {
            egui::ScrollArea::vertical().id_salt("pipeline steps").max_height(maximum).auto_shrink([false, true]).show(ui, |ui| self.list(ui, editor, structure));
        });
        ui.horizontal_wrapped(|ui| {
            ui.menu_button("Add Modifier", |ui| {
                if let Some(t) = modifier_menu(ui, "") {
                    self.add_step(editor, structure, &t, None);
                }
            });
            dialog_layout::labelled_check(ui, &mut self.auto_update, "Update automatically").on_hover_text("Recompute the output after every edit");
            if !self.auto_update && ui.add_enabled(editor.dirty(), egui::Button::new("Update")).clicked() {
                *requests = self.update(editor, true);
            }
            if dialog_layout::labelled_check(ui, &mut self.colour_selection, "Colour atoms by selection").on_hover_text("Colour the output's selected atoms (1) apart from the others (0)").changed() && editor.active() {
                *requests = self.update(editor, true);
            }
        });
        // The selected step's form, the command and the file: the rest of the panel.
        let message_height = if self.message.is_empty() { 0.0 } else { 2.0 * (crate::theme::BASE_FONT * 1.4 + 4.0) };
        let height = (ui.available_height() - message_height - 30.0).max(120.0);
        dialog_layout::child_frame(ui, |ui| {
            egui::ScrollArea::vertical().id_salt("step editor").max_height(height).auto_shrink([false, false]).show(ui, |ui| {
                self.editor(ui, editor);
                self.text_section(ui, editor);
                dialog_layout::separator_text(ui, "Pipeline file");
                ui.horizontal_wrapped(|ui| {
                    requests.save_file |= ui.button("Save Pipeline").clicked();
                    requests.load_file |= ui.button("Load Pipeline").on_hover_text("JSON or the text syntax").clicked();
                });
            });
        });
        if !self.message.is_empty() {
            let red = dialog_layout::status_error(ui);
            ui.add(egui::Label::new(egui::RichText::new(&self.message).color(red)).wrap());
        }
    }

    fn list(&mut self, ui: &mut egui::Ui, editor: &mut PipelineEditor, structure: &Structure) {
        let stages = editor.stages().to_vec();
        let p = crate::theme::palette(ui.visuals().dark_mode);
        let frame = crate::theme::FRAME_HEIGHT;
        let spacing = ui.spacing().item_spacing;
        let indent = 20.0;
        let rail = ui.cursor().min.x + 7.0;
        let font = egui::TextStyle::Body.resolve(ui.style());
        // Screen y of each node of the pipe: the input, every step, the output.
        let mut nodes: Vec<f32> = Vec::new();
        let text_row = |ui: &mut egui::Ui, text: &str, colour: egui::Color32, nodes: &mut Vec<f32>| {
            ui.horizontal(|ui| {
                ui.add_space(indent);
                let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), frame), egui::Sense::hover());
                nodes.push(rect.center().y);
                ui.painter().text(rect.left_center(), egui::Align2::LEFT_CENTER, text, egui::TextStyle::Body.resolve(ui.style()), colour);
            });
        };
        match editor.input() {
            Some(input) => text_row(ui, &format!("Input: {} atoms{}", input.atoms.len(), if input.has_unit_cell { ", periodic cell" } else { "" }), p.text, &mut nodes),
            None => text_row(ui, "Input: not started (the active structure is used when it starts)", p.text_muted, &mut nodes),
        }
        let (mut move_from, mut move_to, mut remove, mut duplicate) = (None, None, None, None);
        let (mut replace, mut insert): (Option<(usize, String)>, Option<(usize, String)>) = (None, None);
        for i in 0..editor.size() {
            let m = editor.modifier(i).clone();
            let mut failed = None;
            let mut status = String::new();
            if let Some(stage) = stages.get(i) {
                if !stage.error.is_empty() {
                    failed = Some(stage.error.clone());
                } else if stage.skipped {
                    status = (if m.enabled { "not run (an earlier step failed)" } else { "bypassed" }).to_string();
                } else {
                    status = format!("{} atoms, {} selected, {:.1} ms", stage.atoms, stage.selected, stage.milliseconds);
                }
            }
            let mut status_below = false;
            ui.horizontal(|ui| {
                ui.add_space(indent);
                // The grip: six dots; dragging it (or the name) moves the step.
                let (grip, grip_response) = ui.allocate_exact_size(egui::vec2(12.0, frame), egui::Sense::drag());
                nodes.push(grip.center().y);
                let hot = grip_response.hovered() || grip_response.dragged();
                for r in 0..3 {
                    for c in 0..2 {
                        ui.painter().circle_filled(egui::pos2(grip.left() + 3.5 + 5.0 * c as f32, grip.top() + frame * (0.3 + 0.2 * r as f32)), 1.3, if hot { p.accent } else { p.text_muted });
                    }
                }
                if hot {
                    ui.ctx().set_cursor_icon(egui::CursorIcon::Move);
                }
                grip_response.dnd_set_drag_payload(i);
                let grip_response = grip_response.on_hover_text("Drag to move this step");
                if let Some(from) = grip_response.dnd_release_payload::<usize>() {
                    move_from = Some(*from);
                    move_to = Some(i);
                }
                let mut enabled = m.enabled;
                if dialog_layout::check_box(ui, &mut enabled).on_hover_text(if enabled { "Applied; untick to bypass this step" } else { "Bypassed" }).changed() {
                    let _ = editor.set_enabled(i, enabled);
                }
                // The name, with the result right-aligned in its row when it fits.
                let buttons = 3.0 * frame + 3.0 * spacing.x;
                let name_width = (ui.available_width() - buttons).max(80.0);
                let title = format!("{}  {}", i + 1, title_of(&m.type_id, &m.label));
                let title_colour = if m.enabled { p.text } else { p.text_muted };
                let title_galley = ui.painter().layout_no_wrap(title.clone(), font.clone(), title_colour);
                let title_width = title_galley.size().x + 20.0;
                let fits = |ui: &egui::Ui, text: &str| ui.painter().layout_no_wrap(text.to_string(), font.clone(), p.text_muted).size().x + 2.0 * spacing.x + title_width <= name_width;
                if !status.is_empty() && !fits(ui, &status) {
                    if let Some(stage) = stages.get(i).filter(|s| s.error.is_empty() && !s.skipped) {
                        let shorter = format!("{} atoms, {} sel.", stage.atoms, stage.selected);
                        if fits(ui, &shorter) {
                            status = shorter;
                        }
                    }
                }
                status_below = !status.is_empty() && !fits(ui, &status);
                let (rect, name) = ui.allocate_exact_size(egui::vec2(name_width, frame), egui::Sense::click_and_drag());
                let selected = self.selected == Some(i);
                if selected || name.hovered() {
                    let fill = if selected { ui.visuals().selection.bg_fill } else { egui::Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), if ui.visuals().dark_mode { 61 } else { 36 }) };
                    ui.painter().rect_filled(rect, 4.0, fill);
                }
                ui.painter().galley(egui::pos2(rect.left() + 4.0, rect.center().y - title_galley.size().y / 2.0), title_galley, title_colour);
                if !status_below && !status.is_empty() {
                    ui.painter().text(egui::pos2(rect.right() - 4.0, rect.center().y), egui::Align2::RIGHT_CENTER, &status, font.clone(), p.text_muted);
                }
                if name.clicked() {
                    self.select(Some(i));
                }
                if name.double_clicked() {
                    self.edit_request = Some(i);
                }
                // Drag a step onto another to move it there.
                name.dnd_set_drag_payload(i);
                if let Some(from) = name.dnd_release_payload::<usize>() {
                    move_from = Some(*from);
                    move_to = Some(i);
                }
                let path = crate::menu_table::menu_path_of(&m.type_id);
                let name = name.on_hover_text(format!("Drag to move; right-click to replace, insert or remove; double-click to edit in the {path} dialog"));
                name.context_menu(|ui| {
                    self.select(Some(i));
                    ui.weak(title_of(&m.type_id, &m.label));
                    ui.separator();
                    if ui.add_enabled(!path.is_empty(), egui::Button::new("Edit in Dialog")).on_hover_text(format!("Open {path} on this step")).clicked() {
                        self.edit_request = Some(i);
                        ui.close();
                    }
                    ui.menu_button("Replace With", |ui| {
                        if let Some(t) = modifier_menu(ui, &m.type_id) {
                            replace = Some((i, t));
                        }
                    });
                    ui.menu_button("Insert Above", |ui| {
                        if let Some(t) = modifier_menu(ui, "") {
                            insert = Some((i, t));
                        }
                    });
                    ui.menu_button("Insert Below", |ui| {
                        if let Some(t) = modifier_menu(ui, "") {
                            insert = Some((i + 1, t));
                        }
                    });
                    ui.separator();
                    if ui.add_enabled(i > 0, egui::Button::new("Move Up").shortcut_text("Alt+Up")).clicked() {
                        (move_from, move_to) = (Some(i), Some(i - 1));
                        ui.close();
                    }
                    if ui.add_enabled(i + 1 < editor.size(), egui::Button::new("Move Down").shortcut_text("Alt+Down")).clicked() {
                        (move_from, move_to) = (Some(i), Some(i + 1));
                        ui.close();
                    }
                    if ui.button("Duplicate").clicked() {
                        duplicate = Some(i);
                        ui.close();
                    }
                    if ui.button(if m.enabled { "Disable" } else { "Enable" }).clicked() {
                        let _ = editor.set_enabled(i, !m.enabled);
                        ui.close();
                    }
                    ui.separator();
                    if ui.add(egui::Button::new("Delete").shortcut_text("Del")).clicked() {
                        remove = Some(i);
                        ui.close();
                    }
                });
                if dialog_layout::arrow_button(ui, true, i > 0).on_hover_text("Move up").clicked() {
                    (move_from, move_to) = (Some(i), Some(i - 1));
                }
                if dialog_layout::arrow_button(ui, false, i + 1 < editor.size()).on_hover_text("Move down").clicked() {
                    (move_from, move_to) = (Some(i), Some(i + 1));
                }
                if dialog_layout::cross_button(ui, true).on_hover_text("Remove this step").clicked() {
                    remove = Some(i);
                }
            });
            // An error, or a status too long for the name's row, below it.
            if failed.is_some() || status_below {
                ui.horizontal(|ui| {
                    ui.add_space(indent + 12.0 + frame + 2.0 * spacing.x);
                    match &failed {
                        Some(e) => {
                            let red = dialog_layout::status_error(ui);
                            ui.add(egui::Label::new(egui::RichText::new(e).color(red)).wrap());
                        }
                        None => {
                            ui.label(egui::RichText::new(&status).color(p.text_muted));
                        }
                    }
                });
            }
        }
        // The hint sits beside the pipe, not over it.
        ui.horizontal(|ui| {
            ui.add_space(indent);
            let hint = if editor.size() == 0 { "No modifiers yet: use Add Modifier, or right-click here. The output equals the input." } else { "Drag a step by its grip or name to reorder; right-click it for more." };
            ui.add(egui::Label::new(egui::RichText::new(hint).color(p.text_muted)).wrap());
        });
        if !stages.is_empty() && stages.len() == editor.size() {
            let last = stages.last().unwrap();
            text_row(ui, &format!("Output: {} atoms, {} selected", last.atoms, last.selected), p.text, &mut nodes);
        } else {
            text_row(ui, "Output", p.text, &mut nodes);
        }
        // The pipe: one line through every node, a dot per step.
        if nodes.len() >= 2 {
            ui.painter().vline(rail, nodes[0]..=nodes[nodes.len() - 1], egui::Stroke::new(2.0, p.text_muted));
        }
        for (k, &y) in nodes.iter().enumerate() {
            let ends = k == 0 || k + 1 == nodes.len();
            let off = !ends && !editor.modifier(k - 1).enabled;
            let fill = if ends { p.accent } else if off { p.frame } else { p.button_active };
            let radius = if ends { 5.5 } else { 4.5 };
            ui.painter().circle(egui::pos2(rail, y), radius, fill, egui::Stroke::new(1.2, p.text_muted));
        }
        // Keys for the selected step.
        if let Some(s) = self.selected.filter(|&s| s < editor.size()) {
            if !ui.ctx().egui_wants_keyboard_input() && ui.ui_contains_pointer() {
                let (del, up, down) = ui.input(|i| (i.key_pressed(egui::Key::Delete), i.modifiers.alt && i.key_pressed(egui::Key::ArrowUp), i.modifiers.alt && i.key_pressed(egui::Key::ArrowDown)));
                if del {
                    remove = Some(s);
                } else if up && s > 0 {
                    (move_from, move_to) = (Some(s), Some(s - 1));
                } else if down && s + 1 < editor.size() {
                    (move_from, move_to) = (Some(s), Some(s + 1));
                }
            }
        }
        if let (Some(from), Some(to)) = (move_from, move_to) {
            self.move_step(editor, from, to);
        }
        if let Some(i) = duplicate {
            if let Ok(j) = editor.duplicate(i) {
                self.select(Some(j));
            }
        }
        if let Some(i) = remove {
            self.remove_step(editor, i);
        }
        if let Some((i, t)) = replace {
            match editor.replace(i, &t) {
                Ok(()) => self.select(Some(i)),
                Err(e) => self.message = e,
            }
        }
        if let Some((i, t)) = insert {
            self.add_step(editor, structure, &t, Some(i));
        }
    }

    fn editor(&mut self, ui: &mut egui::Ui, editor: &mut PipelineEditor) {
        let muted = crate::theme::palette(ui.visuals().dark_mode).text_muted;
        let Some(index) = self.selected.filter(|&s| s < editor.size()) else {
            ui.label(egui::RichText::new("Select a step to edit its parameters.").color(muted));
            return;
        };
        let m = editor.modifier(index).clone();
        let Some(t) = find_modifier_type(&m.type_id) else { return };
        if self.texts_for.as_ref() != Some(&(index, m.type_id.clone())) {
            self.texts.clear();
            self.texts_for = Some((index, m.type_id.clone()));
        }
        dialog_layout::separator_text(ui, &format!("Step {}: {}", index + 1, t.title));
        ui.add(egui::Label::new(egui::RichText::new(t.help.trim_end()).color(muted)).wrap());
        // Parameters first: they are what a step is selected for.
        dialog_layout::form_row(ui, "Display name", 0.42, |ui| {
            let mut label = m.label.clone();
            if ui.add(egui::TextEdit::singleline(&mut label).hint_text("optional").desired_width(ui.available_width())).changed() {
                let _ = editor.set_label(index, &label);
            }
        });
        for (name, value) in parameter_form(ui, "pipeline step parameters", &t, &m.parameters, &mut self.texts) {
            if let Err(e) = editor.set_parameter(index, &name, value) {
                self.message = e;
            }
        }
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            // The same dialog as the Build or Edit menu, editing this step.
            let path = crate::menu_table::menu_path_of(&m.type_id);
            if ui.add_enabled(!path.is_empty(), egui::Button::new("Edit in Dialog")).on_hover_text(format!("Open {path} on this step")).clicked() {
                self.edit_request = Some(index);
            }
            let mut enabled = m.enabled;
            if dialog_layout::labelled_check(ui, &mut enabled, "Enabled").on_hover_text("Untick to bypass this step").changed() {
                let _ = editor.set_enabled(index, enabled);
            }
            if ui.add_enabled(index > 0, egui::Button::new("Move Up")).clicked() {
                self.move_step(editor, index, index - 1);
            }
            if ui.add_enabled(index + 1 < editor.size(), egui::Button::new("Move Down")).clicked() {
                self.move_step(editor, index, index + 1);
            }
            if ui.button("Duplicate").clicked() {
                if let Ok(j) = editor.duplicate(index) {
                    self.select(Some(j));
                }
            }
            if dialog_layout::danger_button(ui, "Remove").clicked() {
                self.remove_step(editor, index);
            }
        });
    }

    fn text_section(&mut self, ui: &mut egui::Ui, editor: &mut PipelineEditor) {
        // The pipeline in the shell-pipe syntax, usable with AtomForge --pipe.
        dialog_layout::separator_text(ui, "As a command");
        let text = editor.pipeline().to_text().unwrap_or_else(|e| e);
        if text.is_empty() {
            let muted = crate::theme::palette(ui.visuals().dark_mode).text_muted;
            ui.label(egui::RichText::new("(empty pipeline)").color(muted));
        } else {
            ui.add(egui::Label::new(&text).wrap());
        }
        ui.horizontal_wrapped(|ui| {
            if ui.button("Copy as Command").on_hover_text("Copy an AtomForge --pipe command line that runs these steps").clicked() {
                ui.ctx().copy_text(format!("AtomForge --pipe \"{text}\" --input INPUT --output OUTPUT"));
            }
            if ui.button(if self.show_text { "Cancel Text Edit" } else { "Edit as Text" }).clicked() {
                self.show_text = !self.show_text;
                self.text = text.clone();
            }
        });
        if self.show_text {
            ui.add(egui::TextEdit::multiline(&mut self.text).desired_rows(4).desired_width(f32::INFINITY));
            if ui.button("Apply Text").clicked() {
                match self.apply_text(editor, &self.text.clone()) {
                    Ok(()) => self.show_text = false,
                    Err(e) => self.message = e,
                }
            }
        }
    }

    /// Replaces the steps by a pipeline in the text syntax (Edit as Text, UI scripts).
    pub fn apply_text(&mut self, editor: &mut PipelineEditor, text: &str) -> Result<(), String> {
        let pipeline = Pipeline::parse(text)?;
        editor.set_pipeline(pipeline);
        self.select(if editor.size() > 0 { Some(0) } else { None });
        self.message.clear();
        Ok(())
    }

    /// Writes the pipeline (JSON) to a file.
    pub fn save_file(&mut self, editor: &PipelineEditor, path: &str) {
        let text = editor.pipeline().to_json().dump(2) + "\n";
        if std::fs::write(path, text).is_err() {
            self.message = format!("Cannot write {path}");
        } else {
            self.message.clear();
        }
    }

    /// Reads a pipeline file (JSON or the text syntax).
    pub fn load_file(&mut self, editor: &mut PipelineEditor, path: &str) {
        let result = std::fs::read_to_string(path).map_err(|_| format!("Cannot read {path}")).and_then(|body| {
            if body.trim_start().starts_with('{') {
                Pipeline::from_json(&Json::parse(&body)?)
            } else {
                Pipeline::parse(&body)
            }
        });
        match result {
            Ok(p) => {
                editor.set_pipeline(p);
                self.select(if editor.size() > 0 { Some(0) } else { None });
                self.message.clear();
            }
            Err(e) => self.message = e,
        }
    }

    /// The panel's options and the pipeline for a project (PipelineDialog::snapshot).
    pub fn snapshot(&self, editor: &PipelineEditor) -> String {
        let mut state = editor.to_json();
        state.set("open", Json::from(self.open));
        state.set("auto_update", Json::from(self.auto_update));
        state.set("colour_selection", Json::from(self.colour_selection));
        state.set("selected", Json::from(self.selected.map(|s| s as f64).unwrap_or(-1.0)));
        state.dump(-1)
    }

    /// Restores a project's pipeline (and its input while it was running).
    pub fn restore(&mut self, editor: &mut PipelineEditor, text: &str, input: Option<&Structure>) {
        editor.reset();
        self.selected = None;
        if text.is_empty() {
            return;
        }
        let result = Json::parse(text).and_then(|state| {
            editor.from_json(&state)?;
            if let (Some(input), Some(true)) = (input, state.find("active").and_then(|v| v.boolean().ok())) {
                editor.set_input(input);
            }
            if let Some(v) = state.find("open").and_then(|v| v.boolean().ok()) {
                self.open = v;
            }
            if let Some(v) = state.find("auto_update").and_then(|v| v.boolean().ok()) {
                self.auto_update = v;
            }
            if let Some(v) = state.find("colour_selection").and_then(|v| v.boolean().ok()) {
                self.colour_selection = v;
            }
            if let Some(v) = state.find("selected").and_then(|v| v.number().ok()) {
                self.selected = (v >= 0.0 && (v as usize) < editor.size()).then_some(v as usize);
            }
            Ok(())
        });
        if let Err(e) = result {
            self.message = format!("The saved pipeline could not be restored: {e}");
        }
    }
}
