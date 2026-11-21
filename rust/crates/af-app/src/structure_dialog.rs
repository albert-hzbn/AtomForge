//! The interface every dedicated builder, edit and analysis dialog implements,
//! so the app owns them as a list: the menus open them by step id, the frame
//! draws them, and UI scripts drive them with `dialog STEP COMMAND ARGS`.

use af_core::structure::Structure;
use eframe::egui;

/// What the dialog sees of the app: the structure in the view (empty when no
/// tab is open) and its element colours.
pub struct Input<'a> {
    pub source: &'a Structure,
    pub colours: &'a [[f32; 3]],
}

/// A pipeline step opened in its Build or Edit dialog (the C++ ui/StepEditing.h
/// StepEdit): the dialog is filled from the step's parameters (and, for steps
/// that take a structure, from the structure entering the step), and its main
/// button ("Update Step") hands the new parameters back as Outcome::UpdateStep.
#[derive(Clone, Debug, Default)]
pub struct StepEdit {
    /// The tab (Document::id) whose pipeline holds the step.
    pub tab: u64,
    /// The step's position in the pipeline (0-based) when the dialog opened.
    pub index: usize,
    /// The modifier type, e.g. "build-gb".
    pub step: String,
    pub parameters: af_core::json::Json,
    /// The structure entering the step (empty when the pipeline is not running).
    pub input: Structure,
}

impl StepEdit {
    /// Build steps keep their command-line flags in "options".
    pub fn options(&self) -> String {
        self.parameters.find("options").and_then(|o| o.string().ok()).unwrap_or("").to_string()
    }

    /// The options' flags (an unterminated quote reads as no flags).
    pub fn reader(&self) -> af_pipeline::options::Reader {
        af_pipeline::options::Reader::new(&self.options()).unwrap_or_default()
    }

    /// The step with new options (its other parameters kept).
    pub fn with_options(&self, text: &str) -> StepEdit {
        let mut out = self.clone();
        out.parameters.set("options", text);
        out
    }

    /// The step with one parameter changed.
    pub fn with_parameter(&self, name: &str, value: impl Into<af_core::json::Json>) -> StepEdit {
        let mut out = self.clone();
        out.parameters.set(name, value);
        out
    }

    /// Whether the dialog should start from the step's input structure
    /// ("use_input", true when absent).
    pub fn use_input(&self) -> bool {
        self.parameters.find("use_input").and_then(|v| v.boolean().ok()).unwrap_or(true)
    }
}

/// The note at the top of a dialog editing a pipeline step (dialogLayout::stepBanner).
pub fn step_banner(ui: &mut egui::Ui, step: &str) {
    egui::Frame::new().fill(ui.visuals().faint_bg_color).inner_margin(6.0).corner_radius(4.0).show(ui, |ui| {
        ui.colored_label(
            egui::Color32::from_rgb(90, 150, 230),
            format!("Editing the pipeline step {step}. Update step writes these settings into the pipeline; the structure changes through the pipeline."),
        );
    });
}

/// What a dialog hands back to the app.
pub enum Outcome {
    /// A new structure from a builder: replaces the open one as one undo step
    /// (or opens a tab), framed anew; the text is the status message.
    Built(Structure, String),
    /// An edit of the structure in the view: one undo step, view kept.
    Edited(Structure, String),
    /// Atoms to select in the view (0-based indices into the structure).
    Select(Vec<usize>, String),
    /// Just a status message.
    Status(String),
    /// "Update Step": the pipeline step's new parameters.
    UpdateStep(StepEdit),
    /// "Add to Pipeline": steps appended to the active tab's pipeline.
    AddSteps(Vec<af_pipeline::pipeline::Modifier>, String),
}

pub trait StructureDialog {
    /// The menu step ids (MENU_LOCATIONS) or analysis ids this dialog serves.
    fn steps(&self) -> &'static [&'static str];
    /// The window title (also the palette label when there is no menu entry).
    fn title(&self) -> &'static str;
    fn open(&mut self, step: &str);
    fn is_open(&self) -> bool;
    /// Closes the window (UI tour steps).
    fn close(&mut self);
    /// Draws the window when open; called every frame (also when closed, so
    /// background work can finish).
    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome>;
    /// A UI-script command: `dialog STEP COMMAND ARGS`. Must run synchronously
    /// (no background threads), never open native file pickers or external
    /// programs, and return Err with a readable message on failure. Commands
    /// that check results are named `expect-...`.
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String>;
    /// Opens the dialog on a pipeline step (editStep): fill the controls from
    /// `edit.reader()` / `edit.parameters` (missing flags take the builder's
    /// defaults) and `edit.input`, show "Update Step" instead of the build or
    /// apply button, and return Outcome::UpdateStep with the new options when
    /// it is pressed (also the `update-step` script command). Closing the
    /// dialog, or opening it from the menu, ends the step edit.
    fn edit_step(&mut self, _edit: StepEdit) -> bool {
        false
    }
    /// A file dropped on the window while this dialog is open (feedDroppedFile);
    /// false leaves it to the app, which opens it in a tab.
    fn drop_file(&mut self, _path: &str) -> bool {
        false
    }
    /// Whether this dialog takes every dropped file at once (Custom Structure
    /// and Merge Structures do, ahead of all others, as in the C++).
    fn takes_all_drops(&self) -> bool {
        false
    }
    /// The dialog's state as part of a project tab (the electronic workspace).
    fn workspace(&self) -> Option<af_electronic::project::Workspace> {
        None
    }
    /// Restores the state a project saved; Err when it cannot.
    fn restore_workspace(&mut self, _saved: &af_electronic::project::Workspace) -> Result<(), String> {
        Ok(())
    }
    /// A calculation is running (projects are not saved meanwhile).
    fn busy(&self) -> bool {
        false
    }
}

/// Every dedicated dialog, in menu order.
pub fn all() -> Vec<Box<dyn StructureDialog>> {
    vec![
        Box::new(crate::stacking_fault_dialog::StackingFaultDialog::default()),
        Box::new(crate::csl_dialog::CslDialog::default()),
        Box::new(crate::interface_dialog::InterfaceDialog::default()),
        Box::new(crate::nano_dialog::NanoDialog::default()),
        Box::new(crate::poly_dialog::PolyDialog::default()),
        Box::new(crate::amorphous_dialog::AmorphousDialog::default()),
        Box::new(crate::sss_dialog::SssDialog::default()),
        Box::new(crate::custom_dialog::CustomDialog::default()),
        Box::new(crate::shape_dialog::ShapeDialog::default()),
        Box::new(crate::edit_structure_dialog::EditStructureDialog::default()),
        Box::new(crate::transform_dialog::TransformDialog::default()),
        Box::new(crate::merge_dialog::MergeDialog::default()),
        Box::new(crate::interstitial_dialog::InterstitialDialog::default()),
        Box::new(crate::sculptor_dialog::SculptorDialog::default()),
        Box::new(crate::electronic_dialog::ElectronicDialog::default()),
        Box::new(crate::wannier_dialog::WannierDialog::default()),
        Box::new(crate::lobster_dialog::LobsterDialog::default()),
        Box::new(crate::dislocation_analysis_dialog::DislocationAnalysisDialog::default()),
    ]
}
