//! Edit > Edit Structure (the C++ ui/EditMenuDialogs Edit Structure window,
//! which also serves the add-atom and set-cell menu entries): the cell
//! vectors, and one row per atom with its element (click to substitute),
//! its position in Cartesian or direct coordinates, the selective-dynamics
//! flags and Delete. Each finished edit is one undo step: a drag commits when
//! the mouse is released, a typed value when the field loses focus.

use af_core::cell::{cartesian_to_fractional, fractional_to_cartesian};
use af_core::elements;
use af_core::structure::{AtomSite, Structure};
use eframe::egui;

use crate::dialog_layout as dl;
use glam::DVec3;

use crate::display_settings::periodic_table;
use af_core::json::Json;
use af_pipeline::pipeline::{find_modifier_type, parameter};

use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

#[derive(Default)]
pub struct EditStructureDialog {
    open: bool,
    /// The structure being edited, and the key of the structure it was taken from.
    work: Structure,
    key: Option<u64>,
    /// Edits not yet handed to the app (while a drag or typing is in progress).
    pending: Option<String>,
    direct: bool,
    /// The row whose element is being picked, and the element chosen.
    picking: Option<usize>,
    pick_element: i32,
    scroll_to_bottom: bool,
    /// A pipeline add-atom or set-cell step being edited (the step editor),
    /// and the editor's controls (kept between edits, as the C++ members).
    step: Option<StepEdit>,
    step_state: StepState,
}

/// The controls of the Edit Structure step editor (the C++ m_step* members),
/// with the C++ defaults.
#[derive(Clone, Debug, PartialEq)]
pub struct StepState {
    pub cell: [[f64; 3]; 3],
    pub scale_atoms: bool,
    pub element: i32,
    pub position: [f64; 3],
    pub fractional: bool,
    pub select: bool,
}

impl Default for StepState {
    fn default() -> Self {
        Self { cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], scale_atoms: true, element: 29, position: [0.0; 3], fractional: false, select: true }
    }
}

fn read_vector(v: &Json, out: &mut [f64; 3]) {
    if let Ok(items) = v.items() {
        for (k, item) in items.iter().take(3).enumerate() {
            out[k] = item.number().unwrap_or(0.0);
        }
    }
}

fn vector_json(v: &[f64; 3]) -> Json {
    Json::array(v.iter().map(|&x| Json::from(x)).collect())
}

impl StepState {
    /// EditMenuDialogs::editStep: the controls from the step's parameters
    /// (missing ones take the step's defaults). Default cell vectors stand for
    /// "not set yet", so the cell entering the step is shown. Starts from
    /// `self` (the previous step editor's values) as the C++ members do; the
    /// error is the C++ "Step parameters not read" message.
    pub fn read(&mut self, edit: &StepEdit) -> Result<(), String> {
        let t = find_modifier_type(&edit.step).ok_or(format!("unknown step {}", edit.step))?;
        let get = |name: &str| parameter(&edit.parameters, &t, name);
        if edit.step == "set-cell" {
            let mut defaults = true;
            for (r, name) in ["a", "b", "c"].iter().enumerate() {
                let value = get(name)?;
                defaults = defaults && value.dump(-1) == parameter(&Json::object(), &t, name)?.dump(-1);
                read_vector(&value, &mut self.cell[r]);
            }
            if defaults && edit.input.has_unit_cell {
                self.cell = edit.input.cell_vectors;
            }
            self.scale_atoms = get("scale_atoms")?.boolean().map_err(|e| e.to_string())?;
        } else {
            let z = elements::atomic_number(get("element")?.string().map_err(|e| e.to_string())?);
            if (1..=118).contains(&z) {
                self.element = z;
            }
            read_vector(&get("position")?, &mut self.position);
            self.fractional = get("fractional")?.boolean().map_err(|e| e.to_string())?;
            self.select = get("select")?.boolean().map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    /// EditMenuDialogs::stepParameters: the step's parameters with these values.
    pub fn parameters(&self, step: &str, base: &Json) -> Json {
        let mut p = base.clone();
        if step == "set-cell" {
            p.set("a", vector_json(&self.cell[0]));
            p.set("b", vector_json(&self.cell[1]));
            p.set("c", vector_json(&self.cell[2]));
            p.set("scale_atoms", self.scale_atoms);
        } else {
            p.set("element", elements::symbol(self.element));
            p.set("position", vector_json(&self.position));
            p.set("fractional", self.fractional);
            p.set("select", self.select);
        }
        p
    }

    /// Position Mode: switches Cartesian/Direct, keeping the point where it is
    /// when the cell entering the step allows converting.
    pub fn set_fractional(&mut self, fractional: bool, input: &Structure) {
        if fractional == self.fractional {
            return;
        }
        let from = DVec3::from(self.position);
        let to = if fractional { cartesian_to_fractional(input, from) } else { fractional_to_cartesian(input, from) };
        if let Some(to) = to {
            self.position = [to.x, to.y, to.z];
        }
        self.fractional = fractional;
    }

    /// drawStepEditor.
    fn ui(&mut self, ui: &mut egui::Ui, step: &str, input: &Structure, colours: &[[f32; 3]]) {
        let warn = egui::Color32::from_rgb(230, 170, 60);
        if step == "set-cell" {
            ui.label("Lattice Vectors (Cartesian)");
            for (r, name) in ["a", "b", "c"].iter().enumerate() {
                ui.horizontal(|ui| {
                    for v in self.cell[r].iter_mut() {
                        ui.add(egui::DragValue::new(v).speed(0.01).fixed_decimals(6));
                    }
                    ui.label(*name);
                });
            }
            if input.has_unit_cell && ui.button("Use Cell Entering The Step").clicked() {
                self.cell = input.cell_vectors;
            }
            crate::dialog_widgets::checkbox(ui, &mut self.scale_atoms, "Scale atoms with the cell (keep fractional coordinates)");
            if self.scale_atoms && !input.has_unit_cell {
                ui.colored_label(warn, "Scaling atoms needs a structure with a unit cell.");
            }
        } else {
            ui.label("Add Atom");
            periodic_table(ui, &mut self.element, Some(colours));
            if !(1..=118).contains(&self.element) {
                self.element = 29;
            }
            ui.label(format!("Element: {}", elements::symbol(self.element)));
            let mut fractional = self.fractional;
            egui::ComboBox::from_label("Position Mode").selected_text(if fractional { "Direct" } else { "Cartesian" }).show_ui(ui, |ui| {
                ui.selectable_value(&mut fractional, false, "Cartesian");
                ui.selectable_value(&mut fractional, true, "Direct");
            });
            self.set_fractional(fractional, input);
            ui.horizontal(|ui| {
                for v in self.position.iter_mut() {
                    ui.add(egui::DragValue::new(v).speed(0.005).fixed_decimals(6));
                }
                ui.label(if self.fractional { "Position (fractional)" } else { "Position (Å)" });
            });
            if self.fractional && !input.has_unit_cell {
                ui.colored_label(warn, "A fractional position needs a structure with a unit cell.");
            }
            crate::dialog_widgets::checkbox(ui, &mut self.select, "Select the new atom");
        }
    }
}

fn key(s: &Structure) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    for a in &s.atoms {
        a.atomic_number.hash(&mut h);
        a.fixed.hash(&mut h);
        for v in [a.x, a.y, a.z] {
            v.to_bits().hash(&mut h);
        }
    }
    for row in &s.cell_vectors {
        for v in row {
            v.to_bits().hash(&mut h);
        }
    }
    s.has_unit_cell.hash(&mut h);
    h.finish()
}

fn colour_of(z: i32, colours: &[[f32; 3]]) -> [f32; 3] {
    if z >= 0 && (z as usize) < colours.len() {
        colours[z as usize]
    } else {
        elements::default_colour(z)
    }
}

/// Edit Structure's Add Atom: hydrogen at the origin.
pub fn add_hydrogen(s: &mut Structure, colours: &[[f32; 3]]) {
    let mut atom = AtomSite { symbol: "H".into(), atomic_number: 1, ..Default::default() };
    [atom.r, atom.g, atom.b] = colour_of(1, colours);
    s.atoms.push(atom);
}

/// Gives atom `i` element `z` and its colour.
pub fn set_element(s: &mut Structure, i: usize, z: i32, colours: &[[f32; 3]]) -> Result<(), String> {
    if !(1..=118).contains(&z) {
        return Err(format!("no element with Z = {z}"));
    }
    let atom = s.atoms.get_mut(i).ok_or(format!("no atom {i}"))?;
    atom.atomic_number = z;
    atom.symbol = elements::symbol(z).to_string();
    [atom.r, atom.g, atom.b] = colour_of(z, colours);
    Ok(())
}

/// Sets atom `i`'s direct coordinates, changing only the components that
/// differ from its current ones (so the others stay exact).
pub fn set_fractional(s: &mut Structure, i: usize, frac: [f64; 3]) -> Result<(), String> {
    let atom = s.atoms.get(i).ok_or(format!("no atom {i}"))?;
    let current = cartesian_to_fractional(s, DVec3::new(atom.x, atom.y, atom.z)).ok_or("direct coordinates need a unit cell")?;
    let mut f = current;
    for k in 0..3 {
        if frac[k] != clean(current[k]) {
            f[k] = frac[k];
        }
    }
    let p = fractional_to_cartesian(s, f).ok_or("direct coordinates need a unit cell")?;
    s.atoms[i].set_position([p.x, p.y, p.z]);
    Ok(())
}

/// No "-0.00000": tiny values show as zero.
fn clean(v: f64) -> f64 {
    if v.abs() < 5e-11 {
        0.0
    } else {
        v
    }
}

impl EditStructureDialog {
    /// The atoms table and cell; returns a description of a finished edit.
    fn contents(&mut self, ui: &mut egui::Ui, colours: &[[f32; 3]]) -> (bool, Option<String>) {
        let mut changed = false;
        let mut done: Option<String> = None;
        let p = crate::theme::palette(ui.visuals().dark_mode);
        let h = crate::theme::FRAME_HEIGHT;
        // ImGui drag fields are framed like text fields.
        let framed = |ui: &mut egui::Ui| {
            let w = &mut ui.visuals_mut().widgets;
            w.inactive.weak_bg_fill = p.frame;
            w.hovered.weak_bg_fill = p.frame_hover;
            w.active.weak_bg_fill = p.frame_active;
        };
        dl::muted(ui, "Edit the cell and each atom; click an element to substitute it. Edit > Undo reverts changes.");
        dl::section(ui, "Cell Vectors (Å)");
        if !self.work.has_unit_cell {
            dl::muted(ui, "This structure has no unit cell.");
        } else {
            for (r, name) in ["a", "b", "c"].iter().enumerate() {
                ui.horizontal(|ui| {
                    framed(ui);
                    let each = (520.0 - 2.0 * ui.spacing().item_spacing.x) / 3.0;
                    for c in 0..3 {
                        changed |= ui.add_sized([each, h], egui::DragValue::new(&mut self.work.cell_vectors[r][c]).speed(0.01).fixed_decimals(5)).changed();
                    }
                    ui.label(*name);
                });
            }
        }
        dl::section(ui, &format!("Atoms ({})", self.work.atoms.len()));
        let cell = self.work.has_unit_cell;
        if !cell {
            self.direct = false;
        }
        ui.horizontal(|ui| {
            ui.add_enabled_ui(cell, |ui| {
                egui::ComboBox::from_id_salt("edit coordinates").width(212.0).selected_text(if self.direct { "Direct (fractional)" } else { "Cartesian (Å)" }).show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.direct, false, "Cartesian (Å)");
                    ui.selectable_value(&mut self.direct, true, "Direct (fractional)");
                });
                ui.label("Coordinates");
            });
            if ui.button("Add Atom").on_hover_text("Adds a hydrogen atom at the origin; set its element and position in its row").clicked() {
                add_hydrogen(&mut self.work, colours);
                self.scroll_to_bottom = true;
                done = Some("Added an atom".into());
            }
        });
        if self.work.atoms.is_empty() {
            dl::muted(ui, "No atoms. Add Atom places one at the origin.");
            return (changed, done);
        }
        let header = if self.direct { "Position (u, v, w)" } else { "Position (x, y, z in Å)" };
        // The table: bordered, a header row, rows of fields, filling the
        // height above the footer.
        let footer_height = h + 3.0 * ui.spacing().item_spacing.y + 4.0;
        let table_height = (ui.available_height() - footer_height).max(160.0);
        let spacing = ui.spacing().item_spacing.x;
        let mut delete = None;
        egui::Frame::new().stroke(egui::Stroke::new(1.0, p.border)).show(ui, |ui| {
            ui.set_width(ui.available_width());
            ui.set_height(table_height);
            let width = ui.available_width();
            let (number_w, element_w, fixed_w, delete_w) = (56.0, 76.0, 3.0 * (h - 8.0 + 7.0 + 10.0) + 2.0 * spacing, 64.0);
            let position_w = width - number_w - element_w - fixed_w - delete_w - 4.0 * spacing - 8.0;
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 0.0;
                dl::header_cell(ui, "#", number_w + spacing);
                dl::header_cell(ui, "Element", element_w + spacing);
                dl::header_cell(ui, header, position_w + spacing);
                dl::header_cell(ui, "Fixed", fixed_w + spacing);
                dl::header_cell(ui, "", ui.available_width());
            });
            let row_height = h + ui.spacing().item_spacing.y;
            let mut area = egui::ScrollArea::vertical().id_salt("edit atom rows").auto_shrink([false, false]);
            if std::mem::take(&mut self.scroll_to_bottom) {
                area = area.vertical_scroll_offset(row_height * self.work.atoms.len() as f32);
            }
            let n = self.work.atoms.len();
            area.show_rows(ui, row_height, n, |ui, rows| {
                for i in rows {
                    ui.horizontal(|ui| {
                        ui.add_space(4.0);
                        let muted = p.text_muted;
                        let (r, _) = ui.allocate_exact_size(egui::vec2(number_w - 4.0, h), egui::Sense::hover());
                        ui.painter().text(egui::pos2(r.left(), r.center().y), egui::Align2::LEFT_CENTER, i.to_string(), egui::TextStyle::Body.resolve(ui.style()), muted);
                        let symbol = self.work.atoms[i].symbol.clone();
                        if ui.add_sized([element_w, h], egui::Button::new(symbol)).on_hover_text("Substitute this atom's element").clicked() {
                            self.picking = Some(i);
                            self.pick_element = self.work.atoms[i].atomic_number.clamp(1, 118);
                        }
                        let each = (position_w - 2.0 * spacing) / 3.0;
                        ui.scope(|ui| {
                            framed(ui);
                            if self.direct {
                                let a = &self.work.atoms[i];
                                let frac = cartesian_to_fractional(&self.work, DVec3::new(a.x, a.y, a.z)).unwrap_or_default();
                                let mut edit = [clean(frac.x), clean(frac.y), clean(frac.z)];
                                let mut moved = false;
                                for v in edit.iter_mut() {
                                    moved |= ui.add_sized([each, h], egui::DragValue::new(v).speed(0.001).fixed_decimals(5)).changed();
                                }
                                if moved && set_fractional(&mut self.work, i, edit).is_ok() {
                                    changed = true;
                                }
                            } else {
                                let a = &mut self.work.atoms[i];
                                for v in [&mut a.x, &mut a.y, &mut a.z] {
                                    changed |= ui.add_sized([each, h], egui::DragValue::new(v).speed(0.005).fixed_decimals(5)).changed();
                                }
                            }
                        });
                        let atom = &mut self.work.atoms[i];
                        ui.allocate_ui_with_layout(egui::vec2(fixed_w, h), egui::Layout::left_to_right(egui::Align::Center), |ui| {
                            for (k, axis) in ["x", "y", "z"].iter().enumerate() {
                                let mut fixed = (atom.fixed >> k) & 1 != 0;
                                if crate::dialog_widgets::checkbox(ui, &mut fixed, *axis).on_hover_text("Held fixed in a relaxation (POSCAR selective dynamics)").changed() {
                                    atom.fixed = if fixed { atom.fixed | (1 << k) } else { atom.fixed & !(1 << k) };
                                    changed = true;
                                }
                            }
                        });
                        let red = ui.visuals().error_fg_color;
                        if ui.add(egui::Button::new(egui::RichText::new("Delete").color(red)).min_size(egui::vec2(delete_w, h))).clicked() {
                            delete = Some(i);
                        }
                    });
                }
            });
        });
        if let Some(i) = delete {
            let atom = self.work.atoms[i].clone();
            self.work.erase_atom(i);
            done = Some(format!("Deleted atom {i} ({})", atom.symbol));
        }
        (changed, done)
    }

    fn picker(&mut self, ctx: &egui::Context, colours: &[[f32; 3]]) -> Option<String> {
        let row = self.picking?;
        let mut open = true;
        let mut apply = false;
        let mut cancel = false;
        crate::theme::modal_window("Edit Atom Element").id(egui::Id::unique("edit atom element")).open(&mut open).collapsible(false).resizable(false).show(ctx, |ui| {
            ui.weak("Click an element, then Apply.");
            periodic_table(ui, &mut self.pick_element, Some(colours));
            ui.label(format!("Element: {}", elements::symbol(self.pick_element)));
            ui.horizontal(|ui| {
                apply = ui.button(egui::RichText::new("Apply").strong()).clicked() || ui.input(|i| i.key_pressed(egui::Key::Enter));
                cancel = ui.button("Cancel").clicked();
            });
        });
        if !open || cancel {
            self.picking = None;
            return None;
        }
        if apply {
            self.picking = None;
            let old = self.work.atoms.get(row).map(|a| a.symbol.clone()).unwrap_or_default();
            if set_element(&mut self.work, row, self.pick_element, colours).is_ok() {
                return Some(format!("Substituted atom {row}: {old} → {}", elements::symbol(self.pick_element)));
            }
        }
        None
    }
}

impl EditStructureDialog {
    /// Update Step: the new parameters, ending the step edit.
    fn update_step(&mut self) -> Result<Outcome, String> {
        let edit = self.step.take().ok_or("Edit Structure is not editing a pipeline step")?;
        let parameters = self.step_state.parameters(&edit.step, &edit.parameters);
        self.open = false;
        Ok(Outcome::UpdateStep(StepEdit { parameters, ..edit }))
    }

    /// The Edit Structure window on a pipeline step: only what the step can express.
    fn show_step(&mut self, ctx: &egui::Context, edit: &StepEdit, colours: &[[f32; 3]]) -> Option<Outcome> {
        let mut open = true;
        let (mut update, mut close) = (false, false);
        crate::theme::modal_window("Edit Structure").id(egui::Id::unique("edit structure")).open(&mut open).resizable(false).show(ctx, |ui| {
            ui.set_max_width(560.0);
            crate::structure_dialog::step_banner(ui, &edit.step);
            ui.separator();
            self.step_state.ui(ui, &edit.step, &edit.input, colours);
            ui.separator();
            ui.horizontal(|ui| {
                update = ui.button(egui::RichText::new("Update Step").strong()).clicked();
                close = ui.button("Close").clicked();
            });
        });
        if update {
            return self.update_step().ok();
        }
        if !open || close {
            self.step = None;
            self.open = false;
        }
        None
    }

    /// Step-editor script commands: cell a|b|c X Y Z | scale-atoms on|off |
    /// element SYMBOL | mode cartesian|direct | position X Y Z | select on|off |
    /// update-step | expect-step NAME VALUE (the parameter's JSON).
    fn step_script(&mut self, edit: &StepEdit, command: &str, args: &str) -> Result<Option<Outcome>, String> {
        let words: Vec<&str> = args.split_whitespace().collect();
        let numbers = |from: usize| -> Result<[f64; 3], String> {
            let v: Vec<f64> = words.iter().skip(from).map(|w| w.parse::<f64>()).collect::<Result<_, _>>().map_err(|_| format!("{command} wants three numbers"))?;
            (v.len() == 3).then(|| [v[0], v[1], v[2]]).ok_or(format!("{command} wants three numbers"))
        };
        let on = args.trim() == "on";
        let st = &mut self.step_state;
        match command {
            "cell" => {
                let r = ["a", "b", "c"].iter().position(|n| words.first() == Some(n)).ok_or("cell a|b|c X Y Z")?;
                st.cell[r] = numbers(1)?;
            }
            "scale-atoms" => st.scale_atoms = on,
            "element" => {
                let z = elements::atomic_number(args.trim());
                if !(1..=118).contains(&z) {
                    return Err(format!("no element {}", args.trim()));
                }
                st.element = z;
            }
            "mode" => st.set_fractional(args.trim() == "direct", &edit.input),
            "position" => st.position = numbers(0)?,
            "select" => st.select = on,
            "update-step" => return self.update_step().map(Some),
            "expect-step" => {
                let (name, want) = args.trim().split_once(' ').ok_or("expect-step NAME VALUE")?;
                let got = st.parameters(&edit.step, &edit.parameters).find(name).map(|v| v.dump(-1)).unwrap_or_default();
                if got != want.trim() {
                    return Err(format!("{name} is {got}"));
                }
            }
            other => return Err(format!("the Edit Structure step editor has no command {other}")),
        }
        Ok(None)
    }
}

impl StructureDialog for EditStructureDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn steps(&self) -> &'static [&'static str] {
        &["add-atom", "set-cell"]
    }

    fn title(&self) -> &'static str {
        "Edit Structure"
    }

    fn open(&mut self, _step: &str) {
        self.step = None;
        self.open = true;
        self.key = None;
        self.pending = None;
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        if edit.step != "set-cell" && edit.step != "add-atom" {
            return false;
        }
        if let Err(e) = self.step_state.read(&edit) {
            eprintln!("[Edit Structure] Step parameters not read: {e}");
        }
        self.step = Some(edit);
        self.picking = None;
        self.open = true;
        true
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        if let Some(edit) = self.step.clone() {
            return self.show_step(ctx, &edit, input.colours);
        }
        // Follow the structure in the view (undo, other edits) unless an edit is in progress.
        if self.pending.is_none() {
            let k = key(input.source);
            if self.key != Some(k) {
                self.work = input.source.clone();
                self.key = Some(k);
            }
        }
        let mut open = true;
        let mut close = false;
        let mut finished = None;
        crate::theme::modal_window("Edit Structure").id(egui::Id::unique("edit structure")).open(&mut open).default_size([980.0, 640.0]).resizable(false).show(ctx, |ui| {
            ui.set_min_size(egui::vec2(948.0, 560.0));
            let (changed, done) = self.contents(ui, input.colours);
            if changed && self.pending.is_none() {
                self.pending = Some("Edited the structure".into());
            }
            if done.is_some() {
                finished = done;
            }
            if dl::footer(ui, "Close", None, true) == 1 {
                close = true;
            }
        });
        if let Some(text) = self.picker(ctx, input.colours) {
            finished = Some(text);
        }
        self.open = open && !close;
        // Hand the edit over once nothing is being dragged or typed.
        let settled = !ctx.egui_is_using_pointer() && ctx.memory(|m| m.focused().is_none());
        let text = finished.or_else(|| if settled || !self.open { self.pending.take() } else { None })?;
        self.pending = None;
        self.key = Some(key(&self.work));
        Some(Outcome::Edited(self.work.clone(), text))
    }

    /// cell a|b|c X Y Z | add-atom | element I SYMBOL | position I X Y Z |
    /// fractional I U V W | fixed I FX FY FZ | delete I (rows from 0)
    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        if let Some(edit) = self.step.clone() {
            return self.step_script(&edit, command, args);
        }
        let mut s = input.source.clone();
        let words: Vec<&str> = args.split_whitespace().collect();
        let index = || -> Result<usize, String> {
            let i: usize = words.first().and_then(|w| w.parse().ok()).ok_or(format!("{command} wants an atom row first"))?;
            if i >= s.atoms.len() {
                return Err(format!("no atom {i}"));
            }
            Ok(i)
        };
        let triple = |from: usize| -> Result<[f64; 3], String> {
            let v: Vec<f64> = words.iter().skip(from).map(|w| w.parse::<f64>()).collect::<Result<_, _>>().map_err(|_| format!("{command} wants three numbers"))?;
            (v.len() == 3).then(|| [v[0], v[1], v[2]]).ok_or(format!("{command} wants three numbers"))
        };
        let text = match command {
            "cell" => {
                if !s.has_unit_cell {
                    return Err("This structure has no unit cell.".into());
                }
                let r = ["a", "b", "c"].iter().position(|n| words.first() == Some(n)).ok_or("cell a|b|c X Y Z")?;
                s.cell_vectors[r] = triple(1)?;
                "Edited the cell".to_string()
            }
            "add-atom" => {
                add_hydrogen(&mut s, input.colours);
                "Added an atom".to_string()
            }
            "element" => {
                let i = index()?;
                let z = words.get(1).map(|w| elements::atomic_number(w)).unwrap_or(0);
                set_element(&mut s, i, z, input.colours)?;
                format!("Substituted atom {i}")
            }
            "position" => {
                let i = index()?;
                s.atoms[i].set_position(triple(1)?);
                "Edited the structure".to_string()
            }
            "fractional" => {
                let i = index()?;
                set_fractional(&mut s, i, triple(1)?)?;
                "Edited the structure".to_string()
            }
            "fixed" => {
                let i = index()?;
                let f = triple(1)?;
                s.atoms[i].fixed = (0..3).filter(|&k| f[k] != 0.0).fold(0u8, |m, k| m | (1 << k));
                "Edited the structure".to_string()
            }
            "delete" => {
                let i = index()?;
                s.erase_atom(i);
                format!("Deleted atom {i}")
            }
            "expect-fixed" => {
                let i = index()?;
                let want: u8 = words.get(1).and_then(|w| w.parse().ok()).ok_or("expect-fixed I MASK")?;
                if s.atoms[i].fixed != want {
                    return Err(format!("atom {i} has fixed mask {}", s.atoms[i].fixed));
                }
                return Ok(None);
            }
            other => return Err(format!("Edit Structure has no command {other}")),
        };
        self.key = None;
        self.pending = None;
        Ok(Some(Outcome::Edited(s, text)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic() -> Structure {
        let mut s = Structure { has_unit_cell: true, cell_vectors: [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]], ..Default::default() };
        s.atoms.push(AtomSite::new("Cu", [1.0, 2.0, 3.0]));
        s
    }

    #[test]
    fn direct_edit_keeps_other_components() {
        let mut s = cubic();
        set_fractional(&mut s, 0, [0.5, 0.5, 0.75]).unwrap();
        assert_eq!(s.atoms[0].position(), [2.0, 2.0, 3.0]);
    }

    fn edit(step: &str, parameters: Json, input: Structure) -> StepEdit {
        StepEdit { index: 0, step: step.into(), parameters, input, ..Default::default() }
    }

    #[test]
    fn set_cell_step_round_trip() {
        // Default vectors stand for "not set yet": the cell entering the step is shown.
        let mut st = StepState::default();
        st.read(&edit("set-cell", Json::object(), cubic())).unwrap();
        assert_eq!(st.cell, cubic().cell_vectors);
        assert!(st.scale_atoms);
        // Without a cell entering the step, the defaults stay.
        let mut st = StepState::default();
        st.read(&edit("set-cell", Json::object(), Structure::default())).unwrap();
        assert_eq!(st.cell, StepState::default().cell);
        // Values written back read the same.
        st.cell[0] = [5.0, 0.5, 0.0];
        st.scale_atoms = false;
        let p = st.parameters("set-cell", &Json::object());
        assert_eq!(p.find("a").unwrap().dump(-1), "[5,0.5,0]");
        let mut again = StepState::default();
        again.read(&edit("set-cell", p, cubic())).unwrap();
        assert_eq!(again.cell, st.cell);
        assert!(!again.scale_atoms);
    }

    #[test]
    fn add_atom_step_round_trip() {
        let mut st = StepState::default();
        st.read(&edit("add-atom", Json::object(), cubic())).unwrap();
        assert_eq!((st.element, st.position, st.fractional, st.select), (29, [0.0; 3], false, true));
        st.element = 8;
        st.position = [1.0, 2.0, 3.0];
        st.set_fractional(true, &cubic());
        assert_eq!(st.position, [0.25, 0.5, 0.75]);
        st.select = false;
        let p = st.parameters("add-atom", &Json::object());
        assert_eq!(p.find("element").unwrap().string().unwrap(), "O");
        let mut again = StepState::default();
        again.read(&edit("add-atom", p, cubic())).unwrap();
        assert_eq!(again, st);
        // An unknown element keeps the previous one.
        let mut p = Json::object();
        p.set("element", "Qq");
        again.read(&edit("add-atom", p, cubic())).unwrap();
        assert_eq!(again.element, 8);
    }

    #[test]
    fn add_and_substitute() {
        let mut s = cubic();
        let colours = elements::default_colours();
        add_hydrogen(&mut s, &colours);
        assert_eq!(s.atoms[1].symbol, "H");
        assert_eq!(s.atoms[1].position(), [0.0, 0.0, 0.0]);
        set_element(&mut s, 1, 8, &colours).unwrap();
        assert_eq!((s.atoms[1].symbol.as_str(), s.atoms[1].atomic_number), ("O", 8));
        assert!(set_element(&mut s, 1, 0, &colours).is_err());
    }
}
