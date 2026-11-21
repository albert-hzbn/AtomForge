//! Build > Bulk Crystal (the C++ ui/BulkCrystalBuilderDialog): crystal
//! system and space group (origin choice where the International Tables give
//! two), the lattice parameters the system leaves free, and the asymmetric
//! unit's sites, expanded by the group's operations.

use af_algorithms::bulk_crystal::{
    add_default_asymmetric_atom, apply_element_to_atom, apply_system_constraints, build_bulk_crystal, current_range, has_two_origins, space_group_label, BulkBuildResult, CrystalSystem,
    LatticeParameters,
};
use af_core::structure::{AtomSite, Structure};
use eframe::egui;

use crate::dialog_layout::{self as dl, Note};
use crate::display_settings::periodic_table;
use crate::structure_dialog::StepEdit;
use af_pipeline::options::Writer;

const SYSTEMS: [CrystalSystem; 7] =
    [CrystalSystem::Triclinic, CrystalSystem::Monoclinic, CrystalSystem::Orthorhombic, CrystalSystem::Tetragonal, CrystalSystem::Trigonal, CrystalSystem::Hexagonal, CrystalSystem::Cubic];
const SYSTEM_LABELS: [&str; 7] = ["Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal", "Trigonal", "Hexagonal", "Cubic"];
/// --system names, in CrystalSystem order.
const SYSTEM_NAMES: [&str; 7] = ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"];

pub struct BulkDialog {
    pub open: bool,
    /// Holds the crystal-system list open (manual screenshots; cleared with the dialog).
    pub show_systems: bool,
    system: usize,
    space_group: i32,
    origin: i32,
    lattice: LatticeParameters,
    atoms: Vec<AtomSite>,
    result: BulkBuildResult,
    /// The site whose element is being picked.
    picking: Option<usize>,
    pick_element: i32,
    /// Editing a pipeline step (Edit in Dialog).
    step_edit: Option<StepEdit>,
    /// Sites read from a step get their colours when next drawn.
    recolour: bool,
    /// Update Step was pressed: the step with its new options.
    pub updated: Option<StepEdit>,
}

impl Default for BulkDialog {
    fn default() -> Self {
        Self {
            open: false,
            show_systems: false,
            system: 6,
            space_group: 225,
            origin: 1,
            lattice: LatticeParameters::default(),
            atoms: Vec::new(),
            result: BulkBuildResult::default(),
            picking: None,
            pick_element: 6,
            step_edit: None,
            recolour: false,
            updated: None,
        }
    }
}

impl BulkDialog {
    /// The options of the matching pipeline step (`build-bulk`): --system,
    /// --spacegroup, the lattice, --origin 2 and each --atom.
    pub fn step_options(&self) -> String {
        let mut lp = self.lattice;
        apply_system_constraints(SYSTEMS[self.system], &mut lp);
        let mut out = Writer::new();
        out.add("--system", SYSTEM_NAMES[self.system.min(6)])
            .add_integer("--spacegroup", self.space_group as i64)
            .add_number("--a", lp.a)
            .add_number("--b", lp.b)
            .add_number("--c", lp.c)
            .add_number("--alpha", lp.alpha)
            .add_number("--beta", lp.beta)
            .add_number("--gamma", lp.gamma);
        if has_two_origins(self.space_group) && self.origin == 2 {
            out.add_integer("--origin", 2);
        }
        for a in &self.atoms {
            out.add("--atom", &format!("{} {} {} {}", a.symbol, Writer::number(a.x), Writer::number(a.y), Writer::number(a.z)));
        }
        out.str().to_string()
    }

    /// Opens the dialog on a build-bulk step (editStep); missing flags take
    /// the defaults of AtomForge --build bulk.
    pub fn edit_step(&mut self, edit: StepEdit) -> bool {
        let reader = edit.reader();
        let system = reader.text("--system", "cubic").to_ascii_lowercase();
        self.system = SYSTEM_NAMES.iter().position(|s| *s == system).unwrap_or(6);
        self.space_group = reader.integer("--spacegroup", 225);
        self.origin = if reader.integer("--origin", 1) == 2 { 2 } else { 1 };
        let a = reader.number("--a", 4.0);
        self.lattice = LatticeParameters {
            a,
            b: reader.number("--b", a),
            c: reader.number("--c", a),
            alpha: reader.number("--alpha", 90.0),
            beta: reader.number("--beta", 90.0),
            gamma: reader.number("--gamma", 90.0),
        };
        apply_system_constraints(SYSTEMS[self.system], &mut self.lattice);
        // Each --atom is "SYMBOL fx fy fz"; colours are applied when drawn.
        self.atoms.clear();
        for value in reader.values("--atom") {
            let w: Vec<&str> = value.split_whitespace().collect();
            if w.len() < 4 {
                continue;
            }
            let (Ok(x), Ok(y), Ok(z)) = (w[1].parse::<f64>(), w[2].parse::<f64>(), w[3].parse::<f64>()) else { continue };
            let number = af_core::elements::atomic_number(w[0]);
            if number <= 0 {
                continue;
            }
            self.atoms.push(AtomSite { x, y, z, atomic_number: number, symbol: af_core::elements::symbol(number).to_string(), ..Default::default() });
        }
        self.recolour = true;
        self.step_edit = Some(edit);
        self.open = true;
        true
    }

    /// Bulk Crystal from the menu: builds again (ends a step edit).
    pub fn open_from_menu(&mut self) {
        self.step_edit = None;
        self.open = true;
    }

    /// Update Step: hands back the step with these options and ends the edit.
    pub fn update_step(&mut self) -> Option<StepEdit> {
        let edit = self.step_edit.take()?;
        self.open = false;
        Some(edit.with_options(&self.step_options()))
    }

    /// Sets the dialog from text (UI scripts): system, space group, a, b, c, angles and sites.
    pub fn set(&mut self, what: &str, value: &str, colours: &[[f32; 3]]) -> Result<(), String> {
        let number = |v: &str| v.trim().parse::<f64>().map_err(|e| format!("{what}: {e}"));
        match what {
            "system" => {
                self.system = SYSTEM_NAMES.iter().position(|s| *s == value.trim()).ok_or(format!("no crystal system {value}"))?;
                self.space_group = current_range(SYSTEMS[self.system]).first;
                apply_system_constraints(SYSTEMS[self.system], &mut self.lattice);
            }
            "spacegroup" => {
                let sg = number(value)? as i32;
                let range = current_range(SYSTEMS[self.system]);
                if sg < range.first || sg > range.last {
                    return Err(format!("space group {sg} is not {}", SYSTEM_LABELS[self.system]));
                }
                self.space_group = sg;
            }
            "origin" => self.origin = number(value)? as i32,
            "a" => self.lattice.a = number(value)?,
            "b" => self.lattice.b = number(value)?,
            "c" => self.lattice.c = number(value)?,
            "alpha" => self.lattice.alpha = number(value)?,
            "beta" => self.lattice.beta = number(value)?,
            "gamma" => self.lattice.gamma = number(value)?,
            "sites" => {
                // "Na 0 0 0; Cl 0.5 0.5 0.5"
                let mut atoms = Vec::new();
                for site in value.split(';').filter(|s| !s.trim().is_empty()) {
                    let w: Vec<&str> = site.split_whitespace().collect();
                    if w.len() != 4 {
                        return Err(format!("a site is SYMBOL fx fy fz, not {site}"));
                    }
                    let z = af_core::elements::atomic_number(w[0]);
                    if z <= 0 {
                        return Err(format!("no element {}", w[0]));
                    }
                    let mut a = AtomSite { x: number(w[1])?, y: number(w[2])?, z: number(w[3])?, ..Default::default() };
                    apply_element_to_atom(&mut a, z, colours);
                    atoms.push(a);
                }
                self.atoms = atoms;
            }
            _ => return Err(format!("the bulk dialog has no {what}")),
        }
        Ok(())
    }

    /// Build: the crystal, or why not.
    pub fn build(&mut self, colours: &[[f32; 3]]) -> Option<Structure> {
        let system = SYSTEMS[self.system];
        apply_system_constraints(system, &mut self.lattice);
        let mut s = Structure::default();
        let origin = if has_two_origins(self.space_group) { self.origin } else { 1 };
        self.result = build_bulk_crystal(&mut s, system, self.space_group, &self.lattice, &self.atoms, colours, origin);
        self.result.success.then_some(s)
    }

    pub fn result(&self) -> &BulkBuildResult {
        &self.result
    }

    /// Draws the dialog; returns a built crystal.
    pub fn show(&mut self, ctx: &egui::Context, colours: &[[f32; 3]]) -> Option<Structure> {
        if !self.open {
            return None;
        }
        if self.recolour {
            for a in self.atoms.iter_mut() {
                let z = a.atomic_number;
                apply_element_to_atom(a, z, colours);
            }
            self.recolour = false;
        }
        if self.atoms.is_empty() {
            add_default_asymmetric_atom(&mut self.atoms, colours);
        }
        let mut built = None;
        let mut open = self.open;
        crate::theme::modal_window("Build Bulk Crystal").open(&mut open).default_size([656.0, 520.0]).min_width(640.0).min_height(440.0).show(ctx, |ui| {
            if let Some(edit) = &self.step_edit {
                crate::structure_dialog::step_banner(ui, &edit.step);
            }
            let spacing = ui.spacing().item_spacing.x;
            let half = ((ui.available_width() - spacing) * 0.5).floor();
            dl::section(ui, "Crystal symmetry");
            ui.horizontal(|ui| {
                let r = ui.label("Crystal system");
                ui.add_space((half - r.rect.width()).max(0.0));
                ui.label("Space group");
            });
            ui.horizontal(|ui| {
                let previous = self.system;
                let combo = egui::ComboBox::from_id_salt("bulk system").width(half - 8.0).selected_text(SYSTEM_LABELS[self.system]).show_ui(ui, |ui| {
                    for (k, label) in SYSTEM_LABELS.iter().enumerate() {
                        ui.selectable_value(&mut self.system, k, *label);
                    }
                });
                if self.show_systems {
                    // Scripts have no pointer to open the list with: the same list, held open.
                    let rect = combo.response.rect;
                    egui::Area::new(egui::Id::unique("bulk systems list")).order(egui::Order::Tooltip).fixed_pos(rect.left_bottom() + egui::vec2(0.0, 2.0)).show(ui.ctx(), |ui| {
                        egui::Frame::menu(ui.style()).show(ui, |ui| {
                            ui.set_width(rect.width() - 12.0);
                            ui.with_layout(egui::Layout::top_down_justified(egui::Align::Min), |ui| {
                                for (k, label) in SYSTEM_LABELS.iter().enumerate() {
                                    ui.add(egui::Button::selectable(k == self.system, *label).frame_when_inactive(false).min_size(egui::vec2(rect.width() - 12.0, crate::theme::FRAME_HEIGHT - 6.0)));
                                }
                            });
                        });
                    });
                }
                if self.system != previous {
                    self.space_group = current_range(SYSTEMS[self.system]).first;
                    apply_system_constraints(SYSTEMS[self.system], &mut self.lattice);
                }
                let range = current_range(SYSTEMS[self.system]);
                egui::ComboBox::from_id_salt("bulk space group").width(half - 8.0).selected_text(space_group_label(self.space_group)).height(400.0).show_ui(ui, |ui| {
                    for sg in range.first..=range.last {
                        ui.selectable_value(&mut self.space_group, sg, space_group_label(sg));
                    }
                });
            });
            if has_two_origins(self.space_group) {
                ui.horizontal(|ui| {
                    ui.label("Origin choice");
                    crate::dialog_widgets::radio_value(ui, &mut self.origin, 1, "1");
                    crate::dialog_widgets::radio_value(ui, &mut self.origin, 2, "2");
                })
                .response
                .on_hover_text("The International Tables give this group two origins; site coordinates differ between them\n(Si in Fd-3m: 0 0 0 in choice 1, 1/8 1/8 1/8 in choice 2). Match the source of your coordinates.");
            }
            dl::section(ui, "Lattice parameters (Å, °)");
            self.lattice_inputs(ui);
            dl::section(ui, "Asymmetric unit");
            ui.horizontal(|ui| {
                dl::muted(ui, "Sites in fractional coordinates; the space group generates the rest.");
                ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT), egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Reset").on_hover_text("Replaces every site with a single hydrogen site at the origin").clicked() {
                        self.atoms.clear();
                        add_default_asymmetric_atom(&mut self.atoms, colours);
                    }
                    if ui.button("Add Site").on_hover_text("Adds a hydrogen site at the origin; set its element and position in its row").clicked() {
                        add_default_asymmetric_atom(&mut self.atoms, colours);
                    }
                });
            });
            // The sites: a bordered table with a header row, filling the height
            // left above the result and the footer.
            let result_height = if self.result.message.is_empty() && !self.result.success { 0.0 } else { 52.0 };
            let footer_height = crate::theme::FRAME_HEIGHT + 3.0 * ui.spacing().item_spacing.y + 4.0;
            let list_height = (ui.available_height() - footer_height - result_height).max(110.0);
            let mut delete = None;
            let p = crate::theme::palette(ui.visuals().dark_mode);
            egui::Frame::new().stroke(egui::Stroke::new(1.0, p.border)).show(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.set_height(list_height);
                let width = ui.available_width();
                let (number_w, element_w, delete_w) = (40.0, 76.0, 64.0);
                let position_w = width - number_w - element_w - delete_w - 3.0 * spacing - 8.0;
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 0.0;
                    dl::header_cell(ui, "#", number_w + spacing);
                    dl::header_cell(ui, "Element", element_w + spacing);
                    dl::header_cell(ui, "Position (u, v, w)", position_w + spacing);
                    dl::header_cell(ui, "", ui.available_width());
                });
                egui::ScrollArea::vertical().id_salt("bulk sites").auto_shrink([false, false]).show(ui, |ui| {
                    let single = self.atoms.len() == 1;
                    for (i, atom) in self.atoms.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.add_space(4.0);
                            ui.add_sized([number_w - 4.0, crate::theme::FRAME_HEIGHT], egui::Label::new(egui::RichText::new(i.to_string()).color(p.text_muted)).halign(egui::Align::Min));
                            if ui.add_sized([element_w, crate::theme::FRAME_HEIGHT], egui::Button::new(&atom.symbol)).on_hover_text("Substitute this site's element").clicked() {
                                self.picking = Some(i);
                                self.pick_element = atom.atomic_number;
                            }
                            let each = ((position_w - 2.0 * spacing) / 3.0).max(40.0);
                            ui.scope(|ui| {
                                // ImGui drag fields are framed like text fields.
                                let w = &mut ui.visuals_mut().widgets;
                                w.inactive.weak_bg_fill = p.frame;
                                w.hovered.weak_bg_fill = p.frame_hover;
                                w.active.weak_bg_fill = p.frame_active;
                                for v in [&mut atom.x, &mut atom.y, &mut atom.z] {
                                    ui.add_sized([each, crate::theme::FRAME_HEIGHT], egui::DragValue::new(v).speed(0.005).max_decimals(4).min_decimals(4));
                                }
                            });
                            let red = ui.visuals().error_fg_color;
                            if ui.add_enabled(!single, egui::Button::new(egui::RichText::new("Delete").color(red)).min_size(egui::vec2(delete_w, crate::theme::FRAME_HEIGHT))).clicked() {
                                delete = Some(i);
                            }
                        });
                    }
                });
            });
            if let Some(i) = delete {
                self.atoms.remove(i);
                if self.atoms.is_empty() {
                    add_default_asymmetric_atom(&mut self.atoms, colours);
                }
            }
            if !self.result.message.is_empty() || self.result.success {
                ui.add_space(4.0);
                if self.result.success {
                    dl::note(ui, Note::Success, &format!("Built {} atoms in space group {} ({}).", self.result.generated_atoms, self.result.space_group_number, self.result.space_group_symbol));
                } else {
                    dl::note(ui, Note::Error, &self.result.message);
                }
            }
            let primary = if self.step_edit.is_some() { "Update Step" } else { "Build" };
            match dl::footer(ui, primary, Some("Close"), true) {
                1 if self.step_edit.is_some() => self.updated = self.update_step(),
                1 => built = self.build(colours),
                2 => {
                    self.step_edit = None;
                    self.open = false;
                }
                _ => {}
            }
        });
        if !open {
            self.step_edit = None;
        }
        self.open = open && self.open;
        // The element picker for a site.
        if let Some(i) = self.picking {
            let mut keep = true;
            let mut chosen = false;
            crate::theme::modal_window("Site Element").id(egui::Id::unique("bulk element")).collapsible(false).open(&mut keep).show(ctx, |ui| {
                ui.weak("Click an element.");
                chosen = periodic_table(ui, &mut self.pick_element, Some(colours));
            });
            if chosen {
                if let Some(atom) = self.atoms.get_mut(i) {
                    apply_element_to_atom(atom, self.pick_element, colours);
                }
            }
            if chosen || !keep {
                self.picking = None;
            }
        }
        built
    }

    fn lattice_inputs(&mut self, ui: &mut egui::Ui) {
        // Fields 100 px wide with their names after them, and the constraints
        // the system imposes as a muted hint (drawLatticeParameterInputs).
        let system = SYSTEMS[self.system];
        let lp = &mut self.lattice;
        let field = |ui: &mut egui::Ui, salt: &str, label: &str, v: &mut f64| {
            dl::number_field(ui, ("bulk lattice", salt.to_string()), v, 4, 100.0);
            ui.label(label);
            ui.add_space(4.0);
        };
        let hint = |ui: &mut egui::Ui, text: &str| {
            ui.add_space(12.0);
            dl::muted(ui, text);
        };
        match system {
            CrystalSystem::Triclinic => {
                ui.horizontal(|ui| {
                    field(ui, "a", "a", &mut lp.a);
                    field(ui, "b", "b", &mut lp.b);
                    field(ui, "c", "c", &mut lp.c);
                });
                ui.horizontal(|ui| {
                    field(ui, "alpha", "α", &mut lp.alpha);
                    field(ui, "beta", "β", &mut lp.beta);
                    field(ui, "gamma", "γ", &mut lp.gamma);
                });
            }
            CrystalSystem::Monoclinic => {
                ui.horizontal(|ui| {
                    field(ui, "a", "a", &mut lp.a);
                    field(ui, "b", "b", &mut lp.b);
                    field(ui, "c", "c", &mut lp.c);
                });
                ui.horizontal(|ui| {
                    field(ui, "beta", "β", &mut lp.beta);
                    hint(ui, "α = γ = 90°");
                });
            }
            CrystalSystem::Orthorhombic => {
                ui.horizontal(|ui| {
                    field(ui, "a", "a", &mut lp.a);
                    field(ui, "b", "b", &mut lp.b);
                    field(ui, "c", "c", &mut lp.c);
                });
                dl::muted(ui, "α = β = γ = 90°");
            }
            CrystalSystem::Tetragonal | CrystalSystem::Trigonal | CrystalSystem::Hexagonal => {
                ui.horizontal(|ui| {
                    field(ui, "a", "a", &mut lp.a);
                    field(ui, "c", "c", &mut lp.c);
                    hint(
                        ui,
                        match system {
                            CrystalSystem::Tetragonal => "b = a, all angles 90°",
                            CrystalSystem::Trigonal => "hexagonal axes: b = a, γ = 120°",
                            _ => "b = a, γ = 120°",
                        },
                    );
                });
            }
            CrystalSystem::Cubic => {
                ui.horizontal(|ui| {
                    field(ui, "a", "a", &mut lp.a);
                    hint(ui, "b = c = a, all angles 90°");
                });
            }
        }
        apply_system_constraints(system, lp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rock_salt_from_the_dialog() {
        let colours = af_core::elements::default_colours();
        let mut d = BulkDialog::default();
        d.set("system", "cubic", &colours).unwrap();
        d.set("spacegroup", "225", &colours).unwrap();
        d.set("a", "5.64", &colours).unwrap();
        d.set("sites", "Na 0 0 0; Cl 0.5 0.5 0.5", &colours).unwrap();
        let s = d.build(&colours).expect("built");
        assert_eq!(s.atoms.len(), 8);
        assert_eq!(d.result().space_group_number, 225);
        assert_eq!(s.cell_vectors[1][1], 5.64);
        assert!(d.set("spacegroup", "12", &colours).is_err(), "not cubic");
        assert_eq!(d.step_options(), "--system cubic --spacegroup 225 --a 5.64 --b 5.64 --c 5.64 --alpha 90 --beta 90 --gamma 90 --atom \"Na 0 0 0\" --atom \"Cl 0.5 0.5 0.5\"");
        // Edit in Dialog reads the same options back.
        let options = d.step_options();
        let mut e = BulkDialog::default();
        let mut parameters = af_core::json::Json::object();
        parameters.set("options", options.as_str());
        assert!(e.edit_step(StepEdit { step: "build-bulk".into(), parameters, ..Default::default() }));
        assert_eq!(e.step_options(), options);
        assert_eq!(e.update_step().unwrap().options(), options);
        d.set("system", "hexagonal", &colours).unwrap();
        assert_eq!(d.space_group, 168);
        assert_eq!(d.lattice.gamma, 120.0);
    }
}
