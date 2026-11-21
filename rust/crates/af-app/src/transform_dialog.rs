//! Edit > Transform Structure (the C++ ui/TransformAtomsDialog): three rows
//! of integer coefficients give the new cell vectors as combinations of a, b
//! and c; the determinant says how many current cells the new one holds, and
//! Apply rebuilds the structure as that supercell (one undo step). On a
//! pipeline `supercell` step the matrix comes from and goes back to its
//! "matrix" parameter (Update Step).

use af_core::json::Json;
use af_core::structure::Structure;
use af_core::supercell::build_supercell;
use eframe::egui;

use crate::dialog_layout::{self as dl, Note};

use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

const IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

pub struct TransformDialog {
    open: bool,
    matrix: [[i32; 3]; 3],
    step: Option<StepEdit>,
}

impl Default for TransformDialog {
    fn default() -> Self {
        Self { open: false, matrix: IDENTITY, step: None }
    }
}

/// The determinant of an integer matrix (the number of cells in the new one).
pub fn determinant(m: &[[i32; 3]; 3]) -> i64 {
    let m = m.map(|r| r.map(|v| v as i64));
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// The note under the matrix: (text, is a warning).
pub fn determinant_note(det: i64) -> (String, bool) {
    match det {
        0 => ("The determinant is 0: these rows would collapse the cell.".into(), true),
        d if d < 0 => (format!("The determinant is {d}: the new cell would be left-handed. Negate a row to keep it right-handed."), true),
        1 => ("Determinant 1: the same volume and number of atoms.".into(), false),
        d => (format!("Determinant {d}: the new cell holds {d} times the atoms."), false),
    }
}

/// The matrix of a supercell step (editStep): nine integers row by row from
/// its "matrix" text; a missing or short matrix gives the identity.
pub fn matrix_from_parameters(parameters: &Json) -> [[i32; 3]; 3] {
    let mut m = IDENTITY;
    if let Some(text) = parameters.find("matrix").and_then(|v| v.string().ok()) {
        // istream >> int: stops at the first word that is not an integer.
        let values: Vec<i32> = text.split_whitespace().map_while(|w| w.parse::<i32>().ok()).take(9).collect();
        if values.len() == 9 {
            for (k, v) in values.into_iter().enumerate() {
                m[k / 3][k % 3] = v;
            }
        }
    }
    m
}

/// The step's parameters with the matrix written back (stepParameters).
pub fn parameters_with_matrix(parameters: &Json, m: &[[i32; 3]; 3]) -> Json {
    let mut p = parameters.clone();
    let text = m.iter().flatten().map(|v| v.to_string()).collect::<Vec<_>>().join(" ");
    p.set("matrix", text);
    p
}

impl TransformDialog {
    fn apply(&self, source: &Structure) -> Result<Outcome, String> {
        if !source.has_unit_cell {
            return Err("Transform Structure needs a structure with a unit cell.".into());
        }
        let det = determinant(&self.matrix);
        if det == 0 {
            return Err(determinant_note(det).0);
        }
        let s = build_supercell(source, &self.matrix);
        let n = s.atoms.len();
        Ok(Outcome::Edited(s, format!("Transformed the structure: {n} atoms")))
    }

    /// Update Step: the new parameters, ending the step edit.
    fn update_step(&mut self) -> Result<Outcome, String> {
        let edit = self.step.take().ok_or("Transform Structure is not editing a pipeline step")?;
        if determinant(&self.matrix) == 0 {
            self.step = Some(edit);
            return Err(determinant_note(0).0);
        }
        let parameters = parameters_with_matrix(&edit.parameters, &self.matrix);
        self.open = false;
        Ok(Outcome::UpdateStep(StepEdit { parameters, ..edit }))
    }
}

impl StructureDialog for TransformDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn steps(&self) -> &'static [&'static str] {
        &["supercell"]
    }

    fn title(&self) -> &'static str {
        "Transform Structure"
    }

    fn open(&mut self, _step: &str) {
        self.step = None;
        self.open = true;
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        if edit.step != "supercell" {
            return false;
        }
        self.matrix = matrix_from_parameters(&edit.parameters);
        self.step = Some(edit);
        self.open = true;
        true
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        let mut result = None;
        let mut open = true;
        let mut close = false;
        let editing = self.step.as_ref().map(|s| s.step.clone());
        crate::theme::modal_window("Transform Structure").id(egui::Id::unique("transform structure")).open(&mut open).resizable(false).collapsible(false).show(ctx, |ui| {
            if let Some(step) = &editing {
                ui.set_max_width(460.0);
                crate::structure_dialog::step_banner(ui, step);
            }
            dl::section(ui, "Transformation matrix");
            ui.label("Enter three rows of integer coefficients.");
            for (r, row) in self.matrix.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    dl::int_fields(ui, ("transform row", r), row, 300.0, -99..=99);
                });
            }
            ui.add_space(4.0);
            // The determinant is the number of current cells in the new one.
            let det = determinant(&self.matrix);
            let (note, warn) = determinant_note(det);
            if warn {
                dl::note(ui, Note::Warning, &note);
            } else {
                dl::muted(ui, note);
            }
            if editing.is_none() && !input.source.has_unit_cell {
                dl::note(ui, Note::Warning, "This structure has no unit cell to transform.");
            }
            let primary = if editing.is_some() { "Update Step" } else { "Apply" };
            let enabled = det != 0 && (editing.is_some() || input.source.has_unit_cell);
            match dl::footer(ui, primary, Some("Cancel"), enabled) {
                1 if editing.is_some() => {
                    result = self.update_step().ok();
                    close = true;
                }
                1 => {
                    result = self.apply(input.source).ok();
                    close = true;
                }
                2 => close = true,
                _ => {}
            }
        });
        self.open = open && !close;
        if !self.open {
            self.step = None;
        }
        result
    }

    /// matrix M11 M12 ... M33 | apply | update-step | expect-det N
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        match command {
            "matrix" => {
                let v: Vec<i32> = args.split_whitespace().map(|w| w.parse::<i32>()).collect::<Result<_, _>>().map_err(|_| "matrix wants nine integers, row by row".to_string())?;
                if v.len() != 9 {
                    return Err("matrix wants nine integers, row by row".into());
                }
                for k in 0..9 {
                    self.matrix[k / 3][k % 3] = v[k];
                }
                Ok(None)
            }
            "apply" => self.apply(input.source).map(Some),
            "update-step" => self.update_step().map(Some),
            "expect-det" => {
                let want: i64 = args.trim().parse().map_err(|_| "expect-det N")?;
                let det = determinant(&self.matrix);
                if det != want {
                    return Err(format!("determinant {det}"));
                }
                Ok(None)
            }
            other => Err(format!("Transform Structure has no command {other}")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn determinants_and_notes() {
        assert_eq!(determinant(&[[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 1);
        assert_eq!(determinant(&[[2, 0, 0], [0, 3, 0], [0, 0, 1]]), 6);
        assert_eq!(determinant(&[[0, 1, 1], [1, 0, 1], [1, 1, 0]]), 2);
        assert_eq!(determinant(&[[1, 0, 0], [0, 0, 1], [0, 1, 0]]), -1);
        assert!(determinant_note(0).1);
        assert!(determinant_note(-1).1);
        assert_eq!(determinant_note(4).0, "Determinant 4: the new cell holds 4 times the atoms.");
    }

    #[test]
    fn matrix_parameter_round_trip() {
        // As the C++ editStep: missing or short matrices keep the identity.
        assert_eq!(matrix_from_parameters(&Json::object()), IDENTITY);
        let mut p = Json::object();
        p.set("matrix", "2 0 0 0 2");
        assert_eq!(matrix_from_parameters(&p), IDENTITY);
        p.set("matrix", "0 1 1 1 0 1 1 1 0");
        let m = matrix_from_parameters(&p);
        assert_eq!(m, [[0, 1, 1], [1, 0, 1], [1, 1, 0]]);
        p.set("other", 3);
        let back = parameters_with_matrix(&p, &[[2, 0, 0], [0, -1, 0], [0, 0, 3]]);
        assert_eq!(back.find("matrix").unwrap().string().unwrap(), "2 0 0 0 -1 0 0 0 3");
        assert!(back.find("other").is_some());
        assert_eq!(matrix_from_parameters(&back), [[2, 0, 0], [0, -1, 0], [0, 0, 3]]);
    }
}
