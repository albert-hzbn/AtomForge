//! One structure operation from the menus, as the pipeline step it also is
//! (the C++ ui/OperationDialog and ui/ModifierEditor): the step's settings,
//! the atoms it acts on, and Apply (or Select for selection steps).

use af_core::json::Json;
use af_core::structure::Structure;
use af_pipeline::pipeline::{find_modifier_type, make_modifier, ModifierType, PipelineData};
use eframe::egui;

use crate::menu_table::menu_path_of;
use crate::dialog_layout::{self as dl, Note};
use crate::params_form::{description, parameter_form, FormTexts};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Target {
    AllAtoms,
    ViewSelection,
    Condition,
}

/// What an Apply produced.
pub enum Outcome {
    /// The new structure.
    Structure(Structure),
    /// The atoms (indices in Structure::atoms) to select.
    Selection(Vec<usize>),
    /// "Update Step" on a pipeline step opened in this dialog.
    UpdateStep(crate::structure_dialog::StepEdit),
}

pub struct OperationDialog {
    pub open: bool,
    step: String,
    parameters: Json,
    /// Text of each parameter while it is edited.
    texts: FormTexts,
    target: Target,
    condition: String,
    message: String,
    error: bool,
    /// Editing a pipeline step instead of the structure.
    step_edit: Option<crate::structure_dialog::StepEdit>,
    /// Add to Pipeline: steps for the app to append to the active tab's pipeline.
    pub added: Option<Vec<af_pipeline::pipeline::Modifier>>,
}

impl Default for OperationDialog {
    fn default() -> Self {
        Self {
            open: false,
            step: String::new(),
            parameters: Json::object(),
            texts: FormTexts::default(),
            target: Target::AllAtoms,
            condition: String::new(),
            message: String::new(),
            error: false,
            step_edit: None,
            added: None,
        }
    }
}

/// Builders that make a structure need nothing loaded; the rest work on one.
pub fn needs_structure(step: &str) -> bool {
    !matches!(
        step,
        "build-bulk" | "build-nanotube" | "build-monolayer" | "build-nanoribbon" | "build-twisted-bilayer" | "build-cluster" | "build-molecule" | "build-packing" | "build-amorphous"
    )
}

fn acts_on_atoms(step: &str) -> bool {
    step == "delete-selected" || find_modifier_type(step).is_some_and(|t| t.parameters.iter().any(|p| p.name == "target"))
}

fn is_selection_step(step: &str) -> bool {
    find_modifier_type(step).is_some_and(|t| t.category == "Selection")
}

impl OperationDialog {
    pub fn step(&self) -> &str {
        &self.step
    }

    pub fn open_step(&mut self, step: &str) -> bool {
        let Ok(modifier) = make_modifier(step) else { return false };
        if self.step_edit.take().is_some() {
            self.step.clear();
        }
        if self.step != step {
            self.step = step.to_string();
            self.parameters = modifier.parameters;
            self.texts.clear();
            self.message.clear();
        }
        self.open = true;
        true
    }

    /// Opens the dialog on a pipeline step: its parameters, and Update Step.
    pub fn edit_step(&mut self, edit: crate::structure_dialog::StepEdit) -> bool {
        if find_modifier_type(&edit.step).is_none() {
            return false;
        }
        self.step = edit.step.clone();
        self.parameters = edit.parameters.clone();
        self.texts.clear();
        self.message.clear();
        self.error = false;
        self.step_edit = Some(edit);
        self.open = true;
        true
    }

    /// Update Step: the step with the parameters now in the form; ends the edit.
    pub fn update_step(&mut self) -> Option<crate::structure_dialog::StepEdit> {
        let mut edit = self.step_edit.take()?;
        edit.parameters = self.parameters.clone();
        self.open = false;
        Some(edit)
    }

    fn apply(&mut self, t: &ModifierType, structure: &Structure, view_selection: &[usize]) -> Result<Outcome, String> {
        let mut data = PipelineData::new(structure.clone());
        // Selection steps start from the view's selection, so add/subtract/intersect combine with it.
        if is_selection_step(&self.step) || (acts_on_atoms(&self.step) && self.target == Target::ViewSelection) {
            for &i in view_selection {
                if i < data.selected.len() {
                    data.selected[i] = 1;
                }
            }
        }
        let mut parameters = self.parameters.clone();
        if acts_on_atoms(&self.step) {
            match self.target {
                Target::AllAtoms => data.selected.iter_mut().for_each(|s| *s = 1),
                Target::Condition => {
                    let mut condition = make_modifier("select-expression")?.parameters;
                    condition.set("expression", self.condition.as_str());
                    let select = find_modifier_type("select-expression").ok_or("select-expression is missing")?;
                    (select.apply)(&condition, &mut data)?;
                }
                Target::ViewSelection => {}
            }
            if self.target != Target::AllAtoms && data.selected_count() == 0 {
                return Err("No atoms match; nothing was changed".into());
            }
            if parameters.contains("target") {
                parameters.set("target", "selected");
            }
        }
        data.notes.clear();
        (t.apply)(&parameters, &mut data)?;
        self.message = data.notes.last().cloned().unwrap_or_else(|| format!("{} applied", t.title));
        if is_selection_step(&self.step) {
            return Ok(Outcome::Selection(data.selected.iter().enumerate().filter(|(_, &s)| s != 0).map(|(i, _)| i).collect()));
        }
        Ok(Outcome::Structure(data.structure))
    }

    /// Sets a parameter from its text (as typed in the window); UI scripts use this.
    pub fn set_parameter_text(&mut self, name: &str, text: &str) -> Result<(), String> {
        let t = find_modifier_type(&self.step).ok_or("No operation is open")?;
        let value = self.texts.set(&t, name, text)?;
        self.parameters.set(name, value);
        Ok(())
    }

    pub fn set_target(&mut self, target: Target, condition: &str) {
        self.target = target;
        self.condition = condition.to_string();
    }

    /// Presses Apply (or Select); UI scripts use this.
    pub fn apply_now(&mut self, structure: &Structure, view_selection: &[usize]) -> Result<Outcome, String> {
        let t = find_modifier_type(&self.step).ok_or("No operation is open")?;
        let result = self.apply(&t, structure, view_selection);
        if let Err(e) = &result {
            self.message = e.clone();
            self.error = true;
        }
        result
    }

    /// Draws the window; returns the outcome of an Apply.
    pub fn show(&mut self, ctx: &egui::Context, structure: &Structure, view_selection: &[usize]) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        let Some(t) = find_modifier_type(&self.step) else {
            self.open = false;
            return None;
        };
        let path = menu_path_of(&self.step);
        let editing = self.step_edit.is_some();
        // Titled like the menu item that opens it (the step's own title otherwise).
        let name = path.rsplit(" > ").next().filter(|s| !s.is_empty()).unwrap_or(&t.title).to_string();
        let title = if editing { format!("{name} (pipeline step)") } else { name };
        let selecting = is_selection_step(&self.step);
        let mut outcome = None;
        let mut open = self.open;
        // A fixed width and a height that follows the content.
        const WIDTH: f32 = 500.0 - 32.0;
        crate::theme::window(title).id(egui::Id::unique("structure operation")).open(&mut open).resizable(false).min_width(WIDTH).max_width(WIDTH).show(ctx, |ui| {
            ui.set_width(WIDTH);
            if !path.is_empty() {
                dl::muted(ui, &path);
            }
            if let Some(edit) = &self.step_edit {
                dl::note(ui, Note::Info, &format!("Editing the pipeline step {}. Update step writes these settings into the pipeline; the structure changes through the pipeline.", edit.step));
            }
            ui.add(egui::Label::new(description(&t)).wrap());
            if !t.parameters.is_empty() {
                dl::section(ui, "Settings");
                for (name, value) in parameter_form(ui, "operation parameters", &t, &self.parameters, &mut self.texts) {
                    self.parameters.set(&name, value);
                }
            }
            if !editing && acts_on_atoms(&self.step) {
                dl::section(ui, "Apply to");
                ui.horizontal(|ui| {
                    crate::dialog_widgets::radio_value(ui, &mut self.target, Target::AllAtoms, "All atoms");
                    crate::dialog_widgets::radio_value(ui, &mut self.target, Target::ViewSelection, format!("Selected in the view ({})", view_selection.len()));
                    crate::dialog_widgets::radio_value(ui, &mut self.target, Target::Condition, "Atoms where");
                });
                if self.target == Target::Condition {
                    dl::text_field(ui, &mut self.condition, ui.available_width(), "e.g. fz > 0.5 && element == O");
                }
            } else if !editing && selecting {
                ui.add_space(4.0);
                let combines = t.parameters.iter().any(|p| p.name == "mode");
                let text = if combines {
                    format!("Sets the selection in the view ({} atoms selected now), combined with it by the mode above.", view_selection.len())
                } else {
                    format!("Changes the selection in the view ({} atoms selected now).", view_selection.len())
                };
                dl::note(ui, Note::Info, &text);
            }
            if !self.message.is_empty() {
                ui.add_space(4.0);
                dl::note(ui, if self.error { Note::Error } else { Note::Success }, &self.message);
            }
            if editing {
                // A pipeline step acts on the selection made by the steps above it.
                match dl::footer(ui, "Update Step", Some("Close"), true) {
                    1 => outcome = self.update_step().map(Outcome::UpdateStep),
                    2 => self.open = false,
                    _ => {}
                }
                return;
            }
            let action = if selecting { "Select" } else { "Apply" };
            ui.add_space(2.0);
            ui.separator();
            ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT), egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let w = dl::footer_button_width();
                let hint = if selecting { "Change the selection in the view" } else { "Change the structure (Edit > Undo reverts it)" };
                if crate::theme::primary_button(ui, action, Some(egui::vec2(w, crate::theme::FRAME_HEIGHT))).on_hover_text(hint).clicked() {
                    match self.apply(&t, structure, view_selection) {
                        Ok(o) => {
                            self.error = false;
                            outcome = Some(o);
                        }
                        Err(e) => {
                            self.message = e;
                            self.error = true;
                        }
                    }
                }
                if crate::theme::secondary_button(ui, "Add to Pipeline", None).on_hover_text("Append this step to Edit > Structure pipeline instead of changing the structure now").clicked() {
                    self.add_to_pipeline();
                }
            });
        });
        if !open {
            self.step_edit = None;
        }
        self.open = open && self.open;
        outcome
    }

    /// Add to Pipeline: the step (after a select-expression step for a
    /// condition), left in `added` for the app to append to the pipeline.
    fn add_to_pipeline(&mut self) {
        let mut steps = Vec::new();
        if acts_on_atoms(&self.step) && self.target == Target::Condition {
            if let Ok(mut condition) = make_modifier("select-expression") {
                condition.parameters.set("expression", self.condition.as_str());
                steps.push(condition);
            }
        }
        let Ok(mut step) = make_modifier(&self.step) else { return };
        step.parameters = self.parameters.clone();
        if acts_on_atoms(&self.step) && step.parameters.contains("target") {
            step.parameters.set("target", if self.target == Target::AllAtoms { "all" } else { "selected" });
        }
        if self.step == "delete-selected" && self.target == Target::AllAtoms {
            self.message = "Choose a condition: deleting every atom is not a useful pipeline step".into();
            self.error = true;
            return;
        }
        steps.push(step);
        self.added = Some(steps);
        self.message = format!(
            "Added to the structure pipeline{}",
            if self.target == Target::ViewSelection && acts_on_atoms(&self.step) { " (add a selection step there: the view's selection is not part of a pipeline)" } else { "" }
        );
        self.error = false;
    }

}
