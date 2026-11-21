//! Build > Interface Builder (the C++ ui/InterfaceBuilderDialog): two crystals
//! cut along their planes, a search for coincident in-plane supercells (Stradi
//! et al., J. Phys.: Condens. Matter 29, 185901 (2017)), the matches as a
//! sortable table and a strain-against-atoms plot with the Pareto front, and
//! how the selected interface is assembled (strained side, gap, vacuum,
//! registry shift, in-plane repeats).

use std::sync::mpsc::{channel, Receiver};

use af_algorithms::interface::{build_interface, find_matches, prepare_slab, repeat_interface_xy, suggested_match, Assembly, Match, SearchParams, Side, StrainOn};
use af_core::structure::{AtomSite, Structure};
use eframe::egui;

use crate::structure_dialog::{Input, Outcome, StepEdit, StructureDialog};

/// The structure files a crystal can be loaded from (as File > Open).
const STRUCTURE_EXTENSIONS: &[&str] = &["cif", "xyz", "extxyz", "vasp", "poscar", "contcar", "pdb", "mol", "sdf", "mol2", "gjf", "xsf", "pwi", "in", "lmp", "data", "cell", "gen"];

/// One crystal of the interface.
struct Crystal {
    crystal: Structure,
    path: String,
    status: String,
    cleave: bool,
    hkl: [i32; 3],
    layers: i32,
}

impl Default for Crystal {
    fn default() -> Self {
        Self { crystal: Structure::default(), path: String::new(), status: String::new(), cleave: true, hkl: [1, 1, 1], layers: 4 }
    }
}

impl Crystal {
    fn side(&self) -> Side {
        Side { crystal: self.crystal.clone(), cleave: self.cleave, h: self.hkl[0], k: self.hkl[1], l: self.hkl[2], layers: self.layers }
    }
}

type SearchOutput = Result<(Structure, Structure, Vec<Match>), String>;

#[derive(Clone, Copy, PartialEq)]
enum SortColumn {
    Index,
    Atoms,
    Strain,
    Area,
}

pub struct InterfaceDialog {
    open: bool,
    a: Crystal,
    b: Crystal,
    search: SearchParams,
    assembly: Assembly,
    repeat: [i32; 2],
    slab_a: Structure,
    slab_b: Structure,
    matches: Vec<Match>,
    searched: bool,
    /// Settings changed since the search.
    stale: bool,
    search_message: String,
    selected: Option<usize>,
    suggested: Option<usize>,
    running: Option<Receiver<SearchOutput>>,
    preview: Structure,
    preview_dirty: bool,
    sort: (SortColumn, bool),
    scroll_to_selected: bool,
    /// Editing a pipeline step: Update Step instead of Build.
    step: Option<StepEdit>,
    /// Search once the step's crystals are loaded, then select its --pick.
    step_search: bool,
    step_pick: Option<usize>,
    /// The preview of the selected match (drag to turn, scroll to zoom).
    orbit: crate::builder_preview::Orbit,
    /// Load File...: the picker, and which crystal it loads.
    browser: crate::file_browser::FileBrowser,
    browse_slot: usize,
}

impl Default for InterfaceDialog {
    fn default() -> Self {
        Self {
            open: false,
            a: Crystal::default(),
            b: Crystal::default(),
            search: SearchParams::default(),
            assembly: Assembly::default(),
            repeat: [1, 1],
            slab_a: Structure::default(),
            slab_b: Structure::default(),
            matches: Vec::new(),
            searched: false,
            stale: false,
            search_message: String::new(),
            selected: None,
            suggested: None,
            running: None,
            preview: Structure::default(),
            preview_dirty: true,
            sort: (SortColumn::Index, false),
            scroll_to_selected: false,
            step: None,
            step_search: false,
            step_pick: None,
            orbit: crate::builder_preview::Orbit::default(),
            browser: crate::file_browser::FileBrowser::default(),
            browse_slot: 1,
        }
    }
}

fn file_name(path: &str) -> String {
    std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.into())
}

fn strain_on_text(on: StrainOn) -> &'static str {
    match on {
        StrainOn::B => "B (A keeps its lattice)",
        StrainOn::A => "A (B keeps its lattice)",
        StrainOn::Both => "Both (to their average)",
    }
}

/// An fcc crystal (the C++ openExample).
fn fcc(symbol: &str, z: i32, a: f64) -> Structure {
    let mut s = Structure { has_unit_cell: true, cell_vectors: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], ..Default::default() };
    let c = af_core::elements::default_colour(z);
    for f in [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]] {
        s.atoms.push(AtomSite { symbol: symbol.into(), atomic_number: z, x: f[0] * a, y: f[1] * a, z: f[2] * a, r: c[0], g: c[1], b: c[2], ..Default::default() });
    }
    s
}

fn search(a: Side, b: Side, params: SearchParams) -> SearchOutput {
    let slab_a = prepare_slab(&a).map_err(|e| format!("Crystal A: {e}"))?;
    let slab_b = prepare_slab(&b).map_err(|e| format!("Crystal B: {e}"))?;
    let matches = find_matches(&slab_a, &slab_b, &params);
    Ok((slab_a, slab_b, matches))
}

impl InterfaceDialog {
    fn load(&mut self, slot: usize, path: &str) -> Result<(), String> {
        let c = if slot == 1 { &mut self.a } else { &mut self.b };
        match af_io::loader::load_structure_from_file(path) {
            Err(e) => {
                c.status = format!("Could not load {}: {e}", file_name(path));
                Err(c.status.clone())
            }
            Ok(s) if !s.has_unit_cell => {
                c.status = format!("{} has no unit cell; the builder needs a crystal.", file_name(path));
                Err(c.status.clone())
            }
            Ok(s) => {
                c.status = format!("{}  ({} atoms)", file_name(path), s.atoms.len());
                c.crystal = s;
                c.path = path.into();
                self.stale = true;
                Ok(())
            }
        }
    }

    fn use_structure(&mut self, slot: usize, s: &Structure, what: &str) {
        let c = if slot == 1 { &mut self.a } else { &mut self.b };
        c.crystal = s.clone();
        c.path.clear();
        c.status = format!("{what}  ({} atoms)", s.atoms.len());
        self.stale = true;
    }

    fn example(&mut self, cube_on_hexagonal: bool) {
        self.use_structure(1, &fcc("Cu", 29, 3.615), "Example: fcc Cu");
        self.use_structure(2, &fcc("Ni", 28, 3.524), "Example: fcc Ni");
        self.a.cleave = true;
        self.b.cleave = true;
        self.a.hkl = if cube_on_hexagonal { [1, 0, 0] } else { [1, 1, 1] };
        self.b.hkl = [1, 1, 1];
        self.search.max_strain = 0.03;
        self.reset_matches();
        self.searched = false;
    }

    fn reset_matches(&mut self) {
        self.matches.clear();
        self.selected = None;
        self.suggested = None;
        self.preview = Structure::default();
        self.preview_dirty = true;
    }

    fn ready(&self) -> bool {
        !self.a.crystal.atoms.is_empty() && !self.b.crystal.atoms.is_empty()
    }

    fn start_search(&mut self) {
        self.reset_matches();
        self.searched = true;
        self.stale = false;
        self.search_message = "Searching...".into();
    }

    fn finish_search(&mut self, out: SearchOutput) {
        match out {
            Err(e) => self.search_message = e,
            Ok((slab_a, slab_b, matches)) => {
                self.slab_a = slab_a;
                self.slab_b = slab_b;
                self.matches = matches;
                if self.matches.is_empty() {
                    self.search_message = format!(
                        "No matches within {}% strain and {:.0} Å². Raise the strain or area limit, or the repeats.",
                        af_core::format::format_g(self.search.max_strain * 100.0, 3),
                        self.search.max_area
                    );
                    return;
                }
                let s = suggested_match(&self.matches);
                self.suggested = (s >= 0).then_some(s as usize);
                self.search_message = format!("{} {}", self.matches.len(), if self.matches.len() == 1 { "match" } else { "matches" });
                if let Some(s) = self.suggested {
                    self.select(s);
                }
            }
        }
    }

    fn run_search_now(&mut self) {
        self.start_search();
        let out = search(self.a.side(), self.b.side(), self.search);
        self.finish_search(out);
    }

    /// The dialog's settings as the step's options (the flags of AtomForge
    /// --build interface; the C++ stepOptions).
    pub fn step_options(&self) -> String {
        use af_pipeline::options::Writer;
        let mut out = Writer::new();
        if !self.a.path.is_empty() {
            out.add("--layerA", &self.a.path);
        }
        if !self.b.path.is_empty() {
            out.add("--layerB", &self.b.path);
        }
        let hkl = |c: &Crystal| [c.hkl[0] as f64, c.hkl[1] as f64, c.hkl[2] as f64];
        if self.a.cleave {
            out.add_numbers("--hklA", &hkl(&self.a));
        }
        if self.b.cleave {
            out.add_numbers("--hklB", &hkl(&self.b));
        }
        out.add_integer("--layersA", self.a.layers as i64).add_integer("--layersB", self.b.layers as i64);
        out.add_number("--max-strain", self.search.max_strain).add_number("--max-area", self.search.max_area).add_integer("--repeats", self.search.max_repeats as i64);
        out.add_number("--min-angle", self.search.min_cell_angle);
        out.add_number("--angle-min", self.search.angle_min).add_number("--angle-max", self.search.angle_max).add_number("--angle-step", self.search.angle_step);
        let on = match self.assembly.strain_on {
            StrainOn::A => "A",
            StrainOn::Both => "both",
            StrainOn::B => "B",
        };
        out.add("--strain-on", on).add_number("--gap", self.assembly.gap).add_number("--vacuum", self.assembly.vacuum);
        if self.assembly.shift[0] != 0.0 || self.assembly.shift[1] != 0.0 {
            out.add_numbers("--shift", &self.assembly.shift);
        }
        out.add_integer("--repx", self.repeat[0] as i64).add_integer("--repy", self.repeat[1] as i64);
        if let Some(pick) = self.selected.or(self.step_pick) {
            out.add_integer("--pick", pick as i64);
        }
        out.str().to_string()
    }

    /// Fills the dialog from a pipeline step (the C++ editStep): missing flags
    /// take the defaults of AtomForge --build interface.
    fn fill_from_step(&mut self, edit: &StepEdit) {
        let o = edit.reader();
        self.a = Crystal::default();
        self.b = Crystal::default();
        for (flag, c) in [("--hklA", &mut self.a), ("--hklB", &mut self.b)] {
            let values = o.numbers(flag);
            c.cleave = values.len() == 3;
            if c.cleave {
                c.hkl = [values[0].round() as i32, values[1].round() as i32, values[2].round() as i32];
            }
        }
        self.a.layers = o.integer("--layersA", 4).max(1);
        self.b.layers = o.integer("--layersB", 4).max(1);
        self.search = SearchParams::default();
        self.search.max_strain = o.number("--max-strain", 0.03).max(1e-5);
        self.search.max_area = o.number("--max-area", 150.0).max(1.0);
        self.search.max_repeats = o.integer("--repeats", o.integer("--nmax", 6)).clamp(1, 30);
        self.search.min_cell_angle = o.number("--min-angle", 20.0).clamp(0.0, 89.0);
        self.search.angle_min = o.number("--angle-min", 0.0);
        self.search.angle_max = o.number("--angle-max", 180.0);
        self.search.angle_step = o.number("--angle-step", 1.0).max(0.01);
        self.assembly = Assembly::default();
        let on = o.text("--strain-on", "");
        self.assembly.strain_on = if on == "A" || on == "a" { StrainOn::A } else if on == "both" { StrainOn::Both } else { StrainOn::B };
        self.assembly.gap = o.number("--gap", 2.0).max(0.0);
        self.assembly.vacuum = o.number("--vacuum", 0.0).max(0.0);
        let shift = o.numbers("--shift");
        if shift.len() == 2 {
            self.assembly.shift = [shift[0], shift[1]];
        }
        self.repeat = [o.integer("--repx", 1).max(1), o.integer("--repy", 1).max(1)];
        // Crystal A: a file named in the options, else the structure entering the step.
        let path_a = o.text("--layerA", "");
        if !path_a.is_empty() {
            if self.load(1, &path_a).is_err() {
                self.a.path = path_a;
            }
        } else if !edit.input.atoms.is_empty() {
            self.use_structure(1, &edit.input, "Step input");
        }
        let path_b = o.text("--layerB", "");
        if !path_b.is_empty() && self.load(2, &path_b).is_err() {
            self.b.path = path_b;
        }
        self.reset_matches();
        self.searched = false;
        self.step_pick = if o.has("--pick") { usize::try_from(o.integer("--pick", -1)).ok() } else { None };
        self.step_search = true;
    }

    /// Editing a step: search with its settings once both crystals are set,
    /// and select its match.
    fn step_search_if_ready(&mut self) {
        if self.step_search && self.ready() {
            self.step_search = false;
            let pick = self.step_pick;
            self.run_search_now();
            if let Some(pick) = pick {
                self.select(pick);
            }
        }
    }

    fn update_step(&mut self) -> Option<Outcome> {
        let edit = self.step.take()?;
        self.open = false;
        self.step_search = false;
        Some(Outcome::UpdateStep(edit.with_options(&self.step_options())))
    }

    fn select(&mut self, index: usize) {
        if index >= self.matches.len() {
            return;
        }
        if Some(index) != self.selected {
            self.scroll_to_selected = true;
        }
        self.selected = Some(index);
        self.preview_dirty = true;
    }

    fn rebuild_preview(&mut self) {
        if !self.preview_dirty {
            return;
        }
        self.preview_dirty = false;
        self.preview = Structure::default();
        let Some(m) = self.selected.and_then(|i| self.matches.get(i)) else { return };
        let mut p = build_interface(&self.slab_a, &self.slab_b, m, &self.assembly);
        if self.repeat[0] > 1 || self.repeat[1] > 1 {
            p = repeat_interface_xy(&p, self.repeat[0], self.repeat[1]);
        }
        self.preview = p;
    }

    fn footer_status(&self) -> String {
        let Some(i) = self.selected else { return String::new() };
        let m = &self.matches[i];
        let on = match self.assembly.strain_on {
            StrainOn::Both => "both crystals",
            StrainOn::A => "A",
            StrainOn::B => "B",
        };
        format!("Match {i}: {} atoms, mean strain {:.3}% on {on}, rotation {:.1}°", self.preview.atoms.len(), m.mean_strain * 100.0, m.theta_deg)
    }

    fn crystal_panel(&mut self, ui: &mut egui::Ui, title: &str, slot: usize, active: &Structure) {
        use crate::dialog_layout as dl;
        dl::section(ui, title);
        let mut use_active = false;
        let width = ui.available_width();
        let label_width = (width * 0.47).min(170.0);
        let field_width = width - label_width - 8.0;
        {
            let c = if slot == 1 { &mut self.a } else { &mut self.b };
            if c.crystal.atoms.is_empty() && c.status.is_empty() {
                dl::muted(ui, "Load a crystal, or use the structure in the view.");
            } else {
                ui.add(egui::Label::new(&c.status).wrap());
            }
        }
        ui.horizontal(|ui| {
            if ui.button("Load File...").on_hover_text("Load this crystal from a structure file").clicked() {
                self.browse_slot = slot;
                self.browser.show_for(crate::file_browser::Purpose::Open, if slot == 1 { "Load crystal A" } else { "Load crystal B" }, "", STRUCTURE_EXTENSIONS);
            }
            let ok = !active.atoms.is_empty() && active.has_unit_cell;
            if ui
                .add_enabled(ok, egui::Button::new("Use Active Structure"))
                .on_hover_text("Use the structure in the view as this crystal")
                .on_disabled_hover_text("The structure in the view has no unit cell")
                .clicked()
            {
                use_active = true;
            }
        });
        let c = if slot == 1 { &mut self.a } else { &mut self.b };
        let mut stale = false;
        dl::row_hint(ui, "Interface plane", "Cut the crystal along a plane with the surface builder. Untick for a structure that is already a slab: its a and b vectors span the interface.", label_width, |ui| {
            stale |= crate::dialog_widgets::checkbox(ui, &mut c.cleave, "Cut along (h k l)").changed();
        });
        if c.cleave {
            dl::row(ui, "(h k l)", label_width, |ui| {
                stale |= dl::int_fields(ui, ("interface hkl", slot), &mut c.hkl, field_width, -20..=20);
            });
        }
        dl::row_hint(ui, "Layers", "Repeats of the crystal's translation along the interface normal (atomic layers for most planes).", label_width, |ui| {
            stale |= dl::int_input(ui, ("interface layers", slot), &mut c.layers, field_width, 1..=200);
        });
        self.stale |= stale;
        if use_active {
            self.use_structure(slot, active, "Active structure");
        }
    }

    fn search_settings(&mut self, ui: &mut egui::Ui) {
        use crate::dialog_layout as dl;
        dl::section(ui, "Matching");
        let width = ui.available_width();
        let label_width = (width * 0.47).min(170.0);
        let field_width = width - label_width - 8.0;
        let s = &mut self.search;
        let mut stale = false;
        let mut strain = s.max_strain * 100.0;
        dl::row_hint(ui, "Max strain (%)", "Largest mean absolute strain (|e11| + |e22| + |e12|) / 3 of the strained crystal.", label_width, |ui| {
            stale |= dl::float_input(ui, "interface max strain", &mut strain, 0.1, 2, field_width, 0.01..=50.0);
        });
        s.max_strain = strain / 100.0;
        dl::row_hint(ui, "Max area (Å²)", "Largest interface cell area: the main limit on the number of atoms.", label_width, |ui| {
            stale |= dl::float_input(ui, "interface max area", &mut s.max_area, 10.0, 0, field_width, 1.0..=100000.0);
        });
        let mut repeats = s.max_repeats as i32;
        dl::row_hint(ui, "Repeats", "Largest multiple i, j of B's surface vectors (i b1 + j b2) tried for the interface cell.", label_width, |ui| {
            stale |= dl::int_input(ui, "interface repeats", &mut repeats, field_width, 1..=30);
        });
        s.max_repeats = repeats as _;
        dl::row_hint(ui, "Min cell angle (°)", "Smallest angle between the interface cell vectors (excludes very slanted cells).", label_width, |ui| {
            stale |= dl::float_input(ui, "interface min angle", &mut s.min_cell_angle, 1.0, 0, field_width, 0.0..=89.0);
        });
        let mut angles = [s.angle_min, s.angle_max, s.angle_step];
        dl::row_hint(ui, "Rotation (°)", "Rotations of B about the interface normal that are searched: from, to and step.", label_width, |ui| {
            if dl::float_fields(ui, "interface rotation", &mut angles, 1, field_width) {
                stale = true;
            }
        });
        s.angle_min = angles[0].clamp(-180.0, 360.0);
        s.angle_max = angles[1].clamp(s.angle_min, 360.0);
        s.angle_step = angles[2].clamp(0.05, 90.0);
        self.stale |= stale;
    }

    fn results_table(&mut self, ui: &mut egui::Ui, height: f32) {
        let mut order: Vec<usize> = (0..self.matches.len()).collect();
        let (column, descending) = self.sort;
        let key = |m: &Match, i: usize| match column {
            SortColumn::Atoms => m.atoms as f64,
            SortColumn::Strain => m.mean_strain,
            SortColumn::Area => m.area,
            SortColumn::Index => i as f64,
        };
        order.sort_by(|&x, &y| {
            let (kx, ky) = (key(&self.matches[x], x), key(&self.matches[y], y));
            let o = kx.partial_cmp(&ky).unwrap_or(std::cmp::Ordering::Equal);
            if descending { o.reverse() } else { o }
        });
        let mut clicked = None;
        let mut sort = self.sort;
        let p = crate::theme::palette(ui.visuals().dark_mode);
        let font = egui::TextStyle::Body.resolve(ui.style());
        let widths = [44.0, 70.0, 76.0, 58.0, 58.0, 58.0, 74.0, 78.0, 72.0, 84.0];
        let row_h = crate::theme::FRAME_HEIGHT - 2.0;
        let total: f32 = widths.iter().sum();
        let frame = egui::Frame::new().stroke(egui::Stroke::new(1.0, p.border)).corner_radius(2).inner_margin(egui::Margin::same(1));
        frame.show(ui, |ui| {
            ui.set_width(ui.available_width() - 2.0);
            ui.spacing_mut().item_spacing = egui::vec2(0.0, 0.0);
            // The header row (TableHeadersRow), sortable columns clickable.
            let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), row_h), egui::Sense::hover());
            ui.painter().rect_filled(rect, 0.0, p.button);
            let headers: [(&str, Option<SortColumn>); 10] = [
                ("#", Some(SortColumn::Index)),
                ("Atoms", Some(SortColumn::Atoms)),
                ("Strain %", Some(SortColumn::Strain)),
                ("e11 %", None),
                ("e22 %", None),
                ("e12 %", None),
                ("Area Å²", Some(SortColumn::Area)),
                ("Cell angle", None),
                ("Rotation", None),
                ("Cells A / B", None),
            ];
            let mut x = rect.left();
            for (k, (text, col)) in headers.iter().enumerate() {
                let cell = egui::Rect::from_min_size(egui::pos2(x, rect.top()), egui::vec2(widths[k], row_h));
                let arrow = match col {
                    Some(c) if sort.0 == *c => if sort.1 { " ▼" } else { " ▲" },
                    _ => "",
                };
                ui.painter().text(egui::pos2(cell.left() + 8.0, cell.center().y), egui::Align2::LEFT_CENTER, format!("{text}{arrow}"), font.clone(), p.text);
                if let Some(c) = col {
                    let r = ui.interact(cell, ui.scope_id().with(("interface header", k)), egui::Sense::click());
                    if r.clicked() {
                        sort = if sort.0 == *c { (*c, !sort.1) } else { (*c, false) };
                    }
                }
                if k > 0 {
                    ui.painter().vline(cell.left(), cell.y_range(), egui::Stroke::new(1.0, p.separator));
                }
                x += widths[k];
            }
            egui::ScrollArea::vertical().id_salt("interface matches").max_height(height - row_h).auto_shrink([false, true]).show(ui, |ui| {
                ui.spacing_mut().item_spacing = egui::vec2(0.0, 0.0);
                for &i in &order {
                    let m = &self.matches[i];
                    let width = ui.available_width().max(total);
                    let (row, response) = ui.allocate_exact_size(egui::vec2(width, row_h), egui::Sense::click());
                    let selected = Some(i) == self.selected;
                    if selected {
                        ui.painter().rect_filled(row, 0.0, egui::Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), if ui.visuals().dark_mode { 87 } else { 56 }));
                    } else if response.hovered() {
                        ui.painter().rect_filled(row, 0.0, egui::Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), 30));
                    }
                    ui.painter().hline(row.x_range(), row.bottom(), egui::Stroke::new(1.0, p.separator));
                    let text = if m.pareto { p.text } else { p.text_muted };
                    let cells = [
                        format!("{i}{}", if Some(i) == self.suggested { "  •" } else { "" }),
                        m.atoms.to_string(),
                        format!("{:.3}", m.mean_strain * 100.0),
                        format!("{:.3}", m.strain[0] * 100.0),
                        format!("{:.3}", m.strain[1] * 100.0),
                        format!("{:.3}", m.strain[2] * 100.0),
                        format!("{:.1}", m.area),
                        format!("{:.1}°", m.cell_angle),
                        format!("{:.1}°", m.theta_deg),
                        format!("{} / {}", m.cells_a, m.cells_b),
                    ];
                    let mut x = row.left();
                    for (k, t) in cells.iter().enumerate() {
                        let colour = if k == 1 { text } else { p.text };
                        ui.painter().text(egui::pos2(x + 8.0, row.center().y), egui::Align2::LEFT_CENTER, t, font.clone(), colour);
                        if k > 0 {
                            ui.painter().vline(x, row.y_range(), egui::Stroke::new(1.0, p.separator));
                        }
                        x += widths[k];
                    }
                    let mut tip = if m.pareto { "Pareto-optimal: no other match has as few atoms and less strain." } else { "Another match has as few atoms and less strain." }.to_string();
                    if Some(i) == self.suggested {
                        tip += "\nSuggested (•): lowest strain (%) + area / 25 Å² of the Pareto-optimal matches.";
                    }
                    let response = response.on_hover_text(tip);
                    if response.clicked() {
                        clicked = Some(i);
                    }
                    if selected && self.scroll_to_selected {
                        response.scroll_to_me(Some(egui::Align::Center));
                        self.scroll_to_selected = false;
                    }
                }
            });
        });
        self.sort = sort;
        if let Some(i) = clicked {
            self.select(i);
        }
    }

    /// Strain against atoms, the Pareto front joined; click a point to select it.
    fn plot(&mut self, ui: &mut egui::Ui, size: egui::Vec2) {
        let (response, painter) = ui.allocate_painter(size, egui::Sense::click());
        let rect = response.rect;
        let visuals = ui.visuals().clone();
        let font = egui::FontId::proportional(crate::theme::BASE_FONT * 0.92);
        let text = visuals.weak_text_color();
        painter.text(rect.min + egui::vec2(6.0, 2.0), egui::Align2::LEFT_TOP, "Strain against atoms (click a point)", egui::FontId::proportional(crate::theme::BASE_FONT), text);
        if self.matches.is_empty() {
            return;
        }
        let lo = egui::pos2(rect.min.x + 50.0, rect.min.y + 28.0);
        let hi = egui::pos2(rect.max.x - 8.0, rect.max.y - 34.0);
        if hi.x - lo.x < 40.0 || hi.y - lo.y < 40.0 {
            return;
        }
        let (mut x_max, mut y_min, mut y_max) = (0.0f64, f64::MAX, 0.0f64);
        for m in &self.matches {
            x_max = x_max.max(m.mean_strain * 100.0);
            y_min = y_min.min(m.atoms as f64);
            y_max = y_max.max(m.atoms as f64);
        }
        x_max = (x_max * 1.08).max(1e-3);
        let span = y_max - y_min;
        let y_min = (y_min - 0.05 * span - 1.0).max(0.0);
        let y_max = y_max + 0.05 * span + 1.0;
        let px = |x: f64| lo.x + (x / x_max) as f32 * (hi.x - lo.x);
        let py = |y: f64| hi.y - ((y - y_min) / (y_max - y_min)) as f32 * (hi.y - lo.y);
        let grid = egui::Stroke::new(1.0, visuals.widgets.noninteractive.bg_stroke.color.gamma_multiply(0.45));
        painter.rect_stroke(egui::Rect::from_min_max(lo, hi), 0.0, visuals.widgets.noninteractive.bg_stroke, egui::StrokeKind::Inside);
        for t in 0..=4 {
            let x = x_max * t as f64 / 4.0;
            let y = y_min + (y_max - y_min) * t as f64 / 4.0;
            painter.line_segment([egui::pos2(px(x), lo.y), egui::pos2(px(x), hi.y)], grid);
            painter.line_segment([egui::pos2(lo.x, py(y)), egui::pos2(hi.x, py(y))], grid);
            painter.text(egui::pos2(px(x), hi.y + 2.0), egui::Align2::CENTER_TOP, format!("{x:.2}"), font.clone(), text);
            painter.text(egui::pos2(lo.x - 4.0, py(y)), if t == 0 { egui::Align2::RIGHT_BOTTOM } else { egui::Align2::RIGHT_CENTER }, format!("{y:.0}"), font.clone(), text);
        }
        painter.text(egui::pos2((lo.x + hi.x) * 0.5, hi.y + 16.0), egui::Align2::CENTER_TOP, "mean strain (%)", font, text);
        let accent = crate::theme::palette(visuals.dark_mode).accent;
        let muted = text.gamma_multiply(0.75);
        let front: Vec<egui::Pos2> = self.matches.iter().filter(|m| m.pareto).map(|m| egui::pos2(px(m.mean_strain * 100.0), py(m.atoms as f64))).collect();
        if front.len() >= 2 {
            painter.add(egui::Shape::line(front, egui::Stroke::new(1.5, accent)));
        }
        let mouse = response.hover_pos();
        let mut hovered = None;
        let mut nearest = 49.0f32;
        for (i, m) in self.matches.iter().enumerate() {
            let p = egui::pos2(px(m.mean_strain * 100.0), py(m.atoms as f64));
            painter.circle_filled(p, if m.pareto { 4.0 } else { 3.0 }, if m.pareto { accent } else { muted });
            if let Some(mp) = mouse {
                let d = (mp - p).length_sq();
                if d < nearest {
                    nearest = d;
                    hovered = Some(i);
                }
            }
        }
        if let Some(m) = self.selected.and_then(|i| self.matches.get(i)) {
            painter.circle_stroke(egui::pos2(px(m.mean_strain * 100.0), py(m.atoms as f64)), 7.0, egui::Stroke::new(2.0, visuals.text_color()));
        }
        if let Some(i) = hovered {
            let m = &self.matches[i];
            let tip = format!("Match {i}: {} atoms, strain {:.3}%, area {:.1} Å², rotation {:.1}°", m.atoms, m.mean_strain * 100.0, m.area, m.theta_deg);
            if response.clicked() {
                self.select(i);
            }
            response.on_hover_text(tip);
        }
    }

    fn assembly_panel(&mut self, ui: &mut egui::Ui) {
        use crate::dialog_layout as dl;
        dl::section(ui, "Assembly");
        let a = &mut self.assembly;
        let mut dirty = false;
        let half = ((ui.available_width() - 16.0 - 2.0 * ui.spacing().item_spacing.x) / 2.0).floor();
        let label_width = 180.0f32.min(half * 0.45);
        let field = half - label_width - 8.0;
        ui.horizontal_top(|ui| {
            ui.vertical(|ui| {
                ui.set_width(half);
                dl::row(ui, "Strained crystal", label_width, |ui| {
                    egui::ComboBox::from_id_salt("interface strain on").width(field - 8.0).selected_text(strain_on_text(a.strain_on)).show_ui(ui, |ui| {
                        for on in [StrainOn::B, StrainOn::A, StrainOn::Both] {
                            dirty |= ui.selectable_value(&mut a.strain_on, on, strain_on_text(on)).changed();
                        }
                    });
                });
                dl::row_hint(ui, "Vacuum (Å)", "0 keeps the cell periodic along the normal (two interfaces, each at the gap); more adds vacuum above and below.", label_width, |ui| {
                    dirty |= dl::float_input(ui, "interface vacuum", &mut a.vacuum, 1.0, 1, field, 0.0..=500.0);
                });
                dl::row(ui, "Repeat in plane", label_width, |ui| {
                    dirty |= dl::int_fields(ui, "interface repeat", &mut self.repeat, field * 0.68, 1..=50);
                });
            });
            ui.add_space(16.0);
            ui.vertical(|ui| {
                ui.set_width(half);
                dl::row_hint(ui, "Gap (Å)", "Distance between the facing atomic layers of A and B.", label_width, |ui| {
                    dirty |= dl::float_input(ui, "interface gap", &mut a.gap, 0.1, 2, field, 0.0..=50.0);
                });
                dl::row_hint(ui, "Shift of B (a, b)", "In-plane shift of B as fractions of the interface cell vectors: sets the registry of the two crystals.", label_width, |ui| {
                    dirty |= dl::float_fields(ui, "interface shift", &mut a.shift, 3, field);
                });
            });
        });
        self.preview_dirty |= dirty;
    }

    fn body(&mut self, ui: &mut egui::Ui, input: &Input) -> bool {
        use crate::dialog_layout as dl;
        let mut build = false;
        // The dialog keeps its size: the content fills it and never widens it.
        ui.set_max_width(ui.available_width());
        let footer = 56.0;
        let height = (ui.available_height() - footer).max(300.0);
        ui.horizontal_top(|ui| {
            ui.allocate_ui_with_layout(egui::vec2(360.0, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                egui::ScrollArea::vertical().id_salt("interface settings").max_height(height).show(ui, |ui| {
                    ui.set_width(340.0);
                    self.crystal_panel(ui, "Crystal A (bottom)", 1, input.source);
                    ui.add_space(10.0);
                    self.crystal_panel(ui, "Crystal B (top)", 2, input.source);
                    ui.add_space(10.0);
                    self.search_settings(ui);
                    ui.add_space(6.0);
                    let ready = self.ready() && self.running.is_none();
                    let w = ui.available_width();
                    if ui.add_enabled_ui(ready, |ui| crate::theme::primary_button(ui, "Find Matches", Some(egui::vec2(w, crate::theme::FRAME_HEIGHT)))).inner.clicked() {
                        self.start_search();
                        let (a, b, p) = (self.a.side(), self.b.side(), self.search);
                        let (tx, rx) = channel();
                        std::thread::spawn(move || {
                            let _ = tx.send(search(a, b, p));
                        });
                        self.running = Some(rx);
                    }
                    if !self.ready() {
                        dl::muted(ui, "Set both crystals first.");
                    } else if self.running.is_some() {
                        dl::muted(ui, "Searching...");
                    } else if self.searched && self.stale {
                        dl::note(ui, dl::Note::Warning, "Settings changed: find matches again.");
                    }
                });
            });
            ui.add_sized([9.0, height], egui::Separator::default().vertical());
            ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), height), egui::Layout::top_down(egui::Align::Min), |ui| {
                let top = ui.cursor().top();
                let left = |ui: &egui::Ui| height - (ui.cursor().top() - top);
                if self.matches.is_empty() {
                    if self.searched && self.running.is_none() {
                        dl::note(ui, dl::Note::Warning, &self.search_message.clone());
                    }
                    ui.add_space(20.0);
                    ui.vertical_centered(|ui| {
                        dl::muted(ui, if self.searched { "No matches yet." } else { "Set both crystals and their planes, then Find Matches.\nThe matches are listed with the fewest atoms first." });
                    });
                    return;
                }
                ui.horizontal_wrapped(|ui| {
                    ui.label(&self.search_message);
                    dl::muted(ui, "(fewest atoms first; dimmed: another match has as few atoms and less strain; • suggested)");
                });
                let assembly_height = 170.0;
                let rest = (left(ui) - assembly_height).max(320.0);
                self.results_table(ui, (rest * 0.42).max(140.0));
                ui.add_space(8.0);
                let panel_height = (left(ui) - assembly_height).max(180.0);
                ui.horizontal_top(|ui| {
                    let w = ui.available_width() - ui.spacing().item_spacing.x - 2.0;
                    let plot_w = (w * 0.42).floor();
                    dl::panel(ui, egui::vec2(plot_w, panel_height), |ui| {
                        let size = ui.available_size();
                        self.plot(ui, egui::vec2(size.x, panel_height - 22.0));
                    });
                    let colours = input.colours.to_vec();
                    dl::panel(ui, egui::vec2(w - plot_w, panel_height), |ui| {
                        let caption = if self.preview.atoms.is_empty() { "The selected match is shown here.".to_string() } else { format!("{} atoms; drag to turn, scroll to zoom", self.preview.atoms.len()) };
                        dl::muted(ui, caption);
                        let size = egui::vec2(ui.available_width(), (panel_height - 50.0).max(60.0));
                        let preview = self.preview.clone();
                        crate::builder_preview::structure(ui, size, &preview, &colours, &mut self.orbit, &[]);
                    });
                });
                self.assembly_panel(ui);
            });
        });
        self.rebuild_preview();
        let label = if self.step.is_some() { "Update Step" } else { "Build" };
        let can_build = self.step.is_some() || !self.preview.atoms.is_empty();
        ui.add_space(2.0);
        ui.separator();
        match dl::action_bar(ui, &self.footer_status(), dl::Note::Info, label, Some("Close"), can_build) {
            1 => build = true,
            2 => self.open = false,
            _ => {}
        }
        build
    }

    fn outcome(&self) -> Outcome {
        Outcome::Built(self.preview.clone(), "Built an interface".into())
    }
}

impl StructureDialog for InterfaceDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        let slot = if self.a.crystal.atoms.is_empty() { 1 } else { 2 };
        let _ = self.load(slot, path);
        true
    }
    fn steps(&self) -> &'static [&'static str] {
        &["build-interface"]
    }

    fn title(&self) -> &'static str {
        "Interface Builder"
    }

    fn edit_step(&mut self, edit: StepEdit) -> bool {
        self.fill_from_step(&edit);
        self.step = Some(edit);
        self.open = true;
        self.step_search_if_ready();
        true
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
        // From the menu: build normally.
        self.step = None;
        self.step_search = false;
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, input: &Input) -> Option<Outcome> {
        if let Some(rx) = &self.running {
            match rx.try_recv() {
                Ok(out) => {
                    self.running = None;
                    self.finish_search(out);
                }
                Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                    self.running = None;
                    self.search_message = "The search stopped unexpectedly.".into();
                }
                Err(_) => ctx.request_repaint_after(std::time::Duration::from_millis(100)),
            }
        }
        if !self.open {
            return None;
        }
        self.step_search_if_ready();
        self.rebuild_preview();
        let mut build = false;
        let mut open = self.open;
        crate::theme::modal_window(self.title()).open(&mut open).default_size([1240.0, 820.0]).min_size([900.0, 600.0]).resizable(true).show(ctx, |ui| {
            if let Some(edit) = &self.step {
                crate::structure_dialog::step_banner(ui, &edit.step);
            }
            build = self.body(ui, input);
        });
        if let Some(path) = self.browser.show(ctx) {
            let _ = self.load(self.browse_slot, &path.to_string_lossy());
        }
        self.open = self.open && open;
        if build && self.step.is_some() {
            return self.update_step();
        }
        if !self.open {
            self.step = None;
            self.step_search = false;
        }
        if build && !self.preview.atoms.is_empty() {
            self.open = false;
            return Some(self.outcome());
        }
        None
    }

    fn script(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let out = self.script_command(command, args, input);
        self.step_search_if_ready();
        out
    }
}

impl InterfaceDialog {
    fn script_command(&mut self, command: &str, args: &str, input: &Input) -> Result<Option<Outcome>, String> {
        let a = args.trim();
        let words: Vec<&str> = a.split_whitespace().collect();
        let num = |i: usize| words.get(i).and_then(|w| w.parse::<f64>().ok()).ok_or(format!("{command}: expected a number in {a:?}"));
        // Commands on one crystal take A or B first.
        let slot = || match words.first().copied() {
            Some("A") | Some("a") => Ok(1),
            Some("B") | Some("b") => Ok(2),
            _ => Err(format!("{command} A|B ...")),
        };
        match command {
            "close" => {
                self.open = false;
                self.step = None;
                self.step_search = false;
            }
            "update-step" => return self.update_step().map(Some).ok_or("not editing a pipeline step".into()),
            "expect-options" => {
                if self.step_options() != a {
                    return Err(format!("options: {}", self.step_options()));
                }
            }
            "example" => {
                // The example pair, shown in the open dialog.
                self.example(a == "cube");
                self.open = true;
            }
            "load" => {
                let s = slot()?;
                let path = a[1..].trim().to_string();
                self.load(s, &path)?;
            }
            "use-view" => {
                let s = slot()?;
                if input.source.atoms.is_empty() || !input.source.has_unit_cell {
                    return Err("the structure in the view has no unit cell".into());
                }
                self.use_structure(s, input.source, "Active structure");
            }
            "hkl" => {
                let s = slot()?;
                let c = if s == 1 { &mut self.a } else { &mut self.b };
                c.hkl = [num(1)? as i32, num(2)? as i32, num(3)? as i32];
                c.cleave = true;
                self.stale = true;
            }
            "slab" => {
                // slab A|B: the crystal is already a slab (no cut).
                let s = slot()?;
                if s == 1 { &mut self.a } else { &mut self.b }.cleave = false;
                self.stale = true;
            }
            "layers" => {
                let s = slot()?;
                if s == 1 { &mut self.a } else { &mut self.b }.layers = (num(1)? as i32).clamp(1, 200);
                self.stale = true;
            }
            "max-strain" => self.search.max_strain = (num(0)? / 100.0).clamp(0.0001, 0.5),
            "max-area" => self.search.max_area = num(0)?.clamp(1.0, 100000.0),
            "repeats" => self.search.max_repeats = (num(0)? as i32).clamp(1, 30),
            "min-angle" => self.search.min_cell_angle = num(0)?.clamp(0.0, 89.0),
            "rotation" => {
                self.search.angle_min = num(0)?;
                self.search.angle_max = num(1)?;
                self.search.angle_step = num(2)?.max(0.05);
            }
            "find" => {
                if !self.ready() {
                    return Err("Set both crystals first.".into());
                }
                self.run_search_now();
                if self.matches.is_empty() {
                    return Err(self.search_message.clone());
                }
            }
            "pick" => {
                let i = num(0)? as usize;
                if i >= self.matches.len() {
                    return Err(format!("no match {i} of {}", self.matches.len()));
                }
                self.select(i);
            }
            "strain-on" => {
                self.assembly.strain_on = match a {
                    "A" | "a" => StrainOn::A,
                    "both" => StrainOn::Both,
                    _ => StrainOn::B,
                };
                self.preview_dirty = true;
            }
            "gap" => {
                self.assembly.gap = num(0)?.max(0.0);
                self.preview_dirty = true;
            }
            "vacuum" => {
                self.assembly.vacuum = num(0)?.max(0.0);
                self.preview_dirty = true;
            }
            "shift" => {
                self.assembly.shift = [num(0)?, num(1)?];
                self.preview_dirty = true;
            }
            "repeat" => {
                self.repeat = [(num(0)? as i32).clamp(1, 50), (num(1)? as i32).clamp(1, 50)];
                self.preview_dirty = true;
            }
            "build" if self.step.is_some() => return Err("editing a pipeline step: use update-step".into()),
            "build" => {
                self.rebuild_preview();
                if self.preview.atoms.is_empty() {
                    return Err("select a match first".into());
                }
                self.open = false;
                return Ok(Some(self.outcome()));
            }
            "expect-matches" => {
                // expect-matches N: at least N matches.
                if (self.matches.len() as f64) < num(0)? {
                    return Err(format!("{} matches ({})", self.matches.len(), self.search_message));
                }
            }
            "expect-selected" => {
                // expect-selected ATOMS: the selected match's interface has this many atoms.
                self.rebuild_preview();
                if self.preview.atoms.len() as f64 != num(0)? {
                    return Err(format!("selected {:?}: {} atoms", self.selected, self.preview.atoms.len()));
                }
            }
            _ => return Err(format!("unknown interface dialog command {command}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edit(options: &str) -> StepEdit {
        let mut parameters = af_core::json::Json::object();
        parameters.set("options", options);
        StepEdit { index: 0, step: "build-interface".into(), parameters, input: Structure::default(), ..Default::default() }
    }

    #[test]
    fn step_options_round_trip() {
        // The C++ stepOptions for the same state (the crystals' files missing:
        // no search, the step's --pick kept).
        for text in [
            "--layersA 4 --layersB 4 --max-strain 0.03 --max-area 150 --repeats 6 --min-angle 20 --angle-min 0 --angle-max 180 --angle-step 1 --strain-on B --gap 2 --vacuum 0 --repx 1 --repy 1",
            "--layerA a.cif --layerB \"b c.cif\" --hklA \"1 1 1\" --hklB \"1 0 0\" --layersA 3 --layersB 5 --max-strain 0.05 --max-area 200 --repeats 8 --min-angle 30 --angle-min -10 --angle-max 90 --angle-step 0.5 --strain-on both --gap 2.5 --vacuum 10 --shift \"0.25 0.5\" --repx 2 --repy 3 --pick 4",
        ] {
            let mut d = InterfaceDialog::default();
            assert!(d.edit_step(edit(text)));
            assert_eq!(d.step_options(), text);
        }
    }

    #[test]
    fn step_search_selects_the_pick() {
        let mut d = InterfaceDialog::default();
        let mut e = edit("--hklA \"1 1 1\" --hklB \"1 1 1\" --pick 1 --repx 2 --repy 2");
        e.input = fcc("Cu", 29, 3.615);
        d.edit_step(e);
        assert!(d.matches.is_empty(), "B is not set yet");
        d.use_structure(2, &fcc("Ni", 28, 3.524), "Ni");
        d.step_search_if_ready();
        assert_eq!(d.selected, Some(1));
        d.rebuild_preview();
        assert_eq!(d.preview.atoms.len(), 400);
    }
}

/// The CSL, Interface and Stacking Fault dialogs' editStep then stepOptions
/// against the C++ (rust/tests/golden/step_grain.tsv, written by
/// tools/golden/golden_step_grain.cpp).
#[cfg(test)]
mod golden_step_options {
    use crate::structure_dialog::{StepEdit, StructureDialog};

    #[test]
    fn matches_cpp_step_options() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/step_grain.tsv");
        let text = std::fs::read_to_string(path).expect("golden step_grain.tsv");
        let mut n = 0;
        for line in text.lines().filter(|l| !l.is_empty()) {
            let mut f = line.splitn(3, '\t');
            let (step, options, expected) = (f.next().unwrap(), f.next().unwrap_or(""), f.next().unwrap_or(""));
            let mut parameters = af_core::json::Json::object();
            parameters.set("options", options);
            let edit = StepEdit { index: 0, step: step.into(), parameters, input: Default::default(), ..Default::default() };
            let got = match step {
                "build-gb" => {
                    let mut d = crate::csl_dialog::CslDialog::default();
                    d.edit_step(edit);
                    d.step_options()
                }
                "build-interface" => {
                    let mut d = super::InterfaceDialog::default();
                    d.edit_step(edit);
                    d.step_options()
                }
                "build-stacking-fault" => {
                    let mut d = crate::stacking_fault_dialog::StackingFaultDialog::default();
                    d.edit_step(edit);
                    d.step_options()
                }
                other => panic!("unknown step {other}"),
            };
            assert_eq!(got, expected, "{step} {options}");
            n += 1;
        }
        assert_eq!(n, 16);
    }
}
