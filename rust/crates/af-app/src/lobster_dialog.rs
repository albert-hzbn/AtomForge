//! Analysis > Electronic Structure > LOBSTER Bonding Analysis (the C++
//! ui/LobsterAnalysisDialog): COHP/COOP/COBI curves per bond (or their
//! average) with the running integral, and the per-bond Fermi-level
//! ICOHP/ICOOP/ICOBI table, from LOBSTER output files.

use af_electronic::lobster::{read_cohpcar, read_icohplist, CohpData, IcohpList};
use eframe::egui;

use crate::analysis_panels::{empty_state, file_name, line_plot, section, two_panels};
use crate::file_browser::{FileBrowser, Purpose};
use crate::structure_dialog::{Input, Outcome, StructureDialog};

#[derive(Default)]
pub struct LobsterDialog {
    open: bool,
    cohpcar_path: String,
    cohpcar_error: String,
    cohpcar: CohpData,
    /// None: the "average" entry.
    selected_bond: Option<usize>,
    spin_channel: usize,
    icohplist_path: String,
    icohplist_error: String,
    icohplist: IcohpList,
    picker: FileBrowser,
    /// The picker is choosing an ICOHPLIST (else a COHPCAR).
    picking_list: bool,
}

impl LobsterDialog {
    pub fn load_cohpcar(&mut self, path: &str) {
        self.cohpcar_error.clear();
        match read_cohpcar(path) {
            Ok(data) => {
                self.cohpcar = data;
                self.cohpcar_path = path.to_string();
                self.selected_bond = None;
                self.spin_channel = 0;
            }
            Err(e) => {
                self.cohpcar_error = e;
                self.cohpcar = CohpData::default();
                self.cohpcar_path.clear();
            }
        }
    }

    pub fn load_icohplist(&mut self, path: &str) {
        self.icohplist_error.clear();
        match read_icohplist(path) {
            Ok(list) => {
                self.icohplist = list;
                self.icohplist_path = path.to_string();
            }
            Err(e) => {
                self.icohplist_error = e;
                self.icohplist = IcohpList::default();
                self.icohplist_path.clear();
            }
        }
    }

    /// The plotted pCOHP and IpCOHP of the selected bond (or the average).
    fn curves(&self) -> (&[f64], &[f64]) {
        let spin = self.spin_channel;
        match self.selected_bond {
            None => (&self.cohpcar.average_cohp[spin], &self.cohpcar.average_icohp[spin]),
            Some(b) => (&self.cohpcar.cohp[spin][b], &self.cohpcar.icohp[spin][b]),
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        section(ui, "COHP / COOP / COBI curves");
        if ui.button("Open COHPCAR").clicked() {
            self.picking_list = false;
            let name = if self.cohpcar_path.is_empty() { "COHPCAR.lobster".to_string() } else { file_name(&self.cohpcar_path) };
            self.picker.show_for(Purpose::Open, "Open LOBSTER COHPCAR", &name, &[]);
        }
        if !self.cohpcar_path.is_empty() {
            ui.label(format!(
                "Loaded: {} ({} bonds, {})",
                file_name(&self.cohpcar_path),
                self.cohpcar.bonds.len(),
                if self.cohpcar.spin_polarized { "spin-polarized" } else { "non-spin-polarized" }
            ));
        }
        if !self.cohpcar_error.is_empty() {
            ui.label(&self.cohpcar_error);
        }
        if !self.cohpcar.energies.is_empty() {
            ui.label(format!("Fermi energy: {:.4} eV", self.cohpcar.fermi_energy));
            if self.cohpcar.spin_polarized {
                egui::ComboBox::from_label("Spin channel").selected_text(["Up", "Down"][self.spin_channel.min(1)]).show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.spin_channel, 0, "Up");
                    ui.selectable_value(&mut self.spin_channel, 1, "Down");
                });
            }
            self.spin_channel = self.spin_channel.min(if self.cohpcar.spin_polarized { 1 } else { 0 });
            egui::Frame::group(ui.style()).show(ui, |ui| {
                egui::ScrollArea::vertical().id_salt("bond list").max_height(160.0).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    if ui.selectable_label(self.selected_bond.is_none(), "Average (all bonds)").clicked() {
                        self.selected_bond = None;
                    }
                    for (i, bond) in self.cohpcar.bonds.iter().enumerate() {
                        let label = format!("Bond {i}: atom {} - atom {} ({:.3} A)", bond.atom1, bond.atom2, bond.length);
                        if ui.selectable_label(self.selected_bond == Some(i), label).clicked() {
                            self.selected_bond = Some(i);
                        }
                    }
                });
            });
        }
        section(ui, "ICOHPLIST / ICOOPLIST / ICOBILIST summary");
        if ui.button("Open ICOHPLIST").clicked() {
            self.picking_list = true;
            let name = if self.icohplist_path.is_empty() { "ICOHPLIST.lobster".to_string() } else { file_name(&self.icohplist_path) };
            self.picker.show_for(Purpose::Open, "Open LOBSTER ICOHPLIST", &name, &[]);
        }
        if !self.icohplist_path.is_empty() {
            ui.label(format!(
                "Loaded: {} ({} entries, {})",
                file_name(&self.icohplist_path),
                self.icohplist.entries.len(),
                if self.icohplist.spin_polarized { "spin-polarized" } else { "non-spin-polarized" }
            ));
        }
        if !self.icohplist_error.is_empty() {
            ui.label(&self.icohplist_error);
        }
    }

    fn results(&self, ui: &mut egui::Ui) {
        if self.cohpcar.energies.is_empty() && self.icohplist.entries.is_empty() {
            empty_state(ui, "Open a LOBSTER COHPCAR.lobster (bonding curves) or ICOHPLIST.lobster (integrated bond strengths); they appear here.");
        }
        if !self.cohpcar.energies.is_empty() {
            let (cohp, icohp) = self.curves();
            section(ui, if self.selected_bond.is_none() { "Average pCOHP / IpCOHP" } else { "Selected-bond pCOHP / IpCOHP" });
            let e = &self.cohpcar.energies;
            ui.label(format!(
                "Series 1: pCOHP. Series 2: IpCOHP (running integral). Fermi level at 0 eV, range {:.3} to {:.3} eV.",
                e.first().copied().unwrap_or(0.0),
                e.last().copied().unwrap_or(0.0)
            ));
            line_plot(ui, "pCOHP and IpCOHP", "Energy - Fermi (eV)", "eV", e, vec![("pCOHP".into(), cohp.to_vec()), ("IpCOHP".into(), icohp.to_vec())], 260.0);
        }
        if !self.icohplist.entries.is_empty() {
            section(ui, "ICOHPLIST entries");
            let spin = self.icohplist.spin_polarized;
            egui::Frame::group(ui.style()).show(ui, |ui| {
                egui::ScrollArea::vertical().id_salt("icohp table").max_height(220.0).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    let mono = |ui: &mut egui::Ui, text: String| ui.label(egui::RichText::new(text).monospace());
                    mono(ui, if spin { "atom1  atom2  length_A  numBonds  icohp_up  icohp_down" } else { "atom1  atom2  length_A  numBonds  icohp" }.to_string());
                    for e in &self.icohplist.entries {
                        if spin {
                            mono(ui, format!("{:5}  {:5}  {:8.4}  {:8}  {:8.5}  {:8.5}", e.atom1, e.atom2, e.length, e.num_bonds, e.icohp[0], e.icohp[1]));
                        } else {
                            mono(ui, format!("{:5}  {:5}  {:8.4}  {:8}  {:8.5}", e.atom1, e.atom2, e.length, e.num_bonds, e.icohp[0]));
                        }
                    }
                });
            });
        }
    }
}

fn expect_close(what: &str, got: f64, want: &str) -> Result<(), String> {
    let want: f64 = want.trim().parse().map_err(|_| format!("{what} wants a number"))?;
    if (got - want).abs() > 1e-6 * want.abs().max(1.0) {
        return Err(format!("{what} is {got}, expected {want}"));
    }
    Ok(())
}

impl StructureDialog for LobsterDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn steps(&self) -> &'static [&'static str] {
        &["lobster"]
    }

    fn title(&self) -> &'static str {
        "LOBSTER Bonding Analysis"
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, _input: &Input) -> Option<Outcome> {
        if !self.open {
            return None;
        }
        let mut open = self.open;
        crate::theme::window("LOBSTER Bonding Analysis").id(egui::Id::unique("lobster analysis")).open(&mut open).collapsible(false).default_size([1180.0, 820.0]).show(ctx, |ui| {
            two_panels(ui, "lobster", self, |s, ui| s.controls(ui), |s, ui| s.results(ui));
        });
        self.open = open;
        if let Some(path) = self.picker.show(ctx) {
            let path = path.to_string_lossy().into_owned();
            if self.picking_list {
                self.load_icohplist(&path);
            } else {
                self.load_cohpcar(&path);
            }
        }
        None
    }

    /// cohpcar PATH | icohplist PATH | bond N|average | spin up|down |
    /// expect-bonds N | expect-entries N | expect-fermi E | expect-icohp VALUE
    /// (the plotted IpCOHP at the Fermi level) | expect-error TEXT
    fn script(&mut self, command: &str, args: &str, _input: &Input) -> Result<Option<Outcome>, String> {
        let args = args.trim();
        match command {
            "cohpcar" => self.load_cohpcar(args),
            "icohplist" => self.load_icohplist(args),
            "bond" => {
                self.selected_bond = if args == "average" {
                    None
                } else {
                    let b: usize = args.parse().map_err(|_| "bond N|average")?;
                    if b >= self.cohpcar.bonds.len() {
                        return Err(format!("no bond {b}: {} bonds", self.cohpcar.bonds.len()));
                    }
                    Some(b)
                }
            }
            "spin" => {
                let s = if args == "down" { 1 } else { 0 };
                if s == 1 && !self.cohpcar.spin_polarized {
                    return Err("the COHPCAR is not spin-polarized".into());
                }
                self.spin_channel = s;
            }
            "expect-bonds" => {
                let n: usize = args.parse().map_err(|_| "expect-bonds N")?;
                if self.cohpcar.bonds.len() != n {
                    return Err(format!("{} bonds (error: {})", self.cohpcar.bonds.len(), self.cohpcar_error));
                }
            }
            "expect-entries" => {
                let n: usize = args.parse().map_err(|_| "expect-entries N")?;
                if self.icohplist.entries.len() != n {
                    return Err(format!("{} ICOHPLIST entries (error: {})", self.icohplist.entries.len(), self.icohplist_error));
                }
            }
            "expect-fermi" => expect_close("Fermi energy", self.cohpcar.fermi_energy, args)?,
            "expect-icohp" => {
                if self.cohpcar.energies.is_empty() {
                    return Err("no COHPCAR loaded".into());
                }
                // The sample closest to 0 eV.
                let e = &self.cohpcar.energies;
                let i = (0..e.len()).min_by(|&a, &b| e[a].abs().total_cmp(&e[b].abs())).unwrap_or(0);
                expect_close("IpCOHP at the Fermi level", self.curves().1[i], args)?;
            }
            "expect-error" => {
                let error = format!("{}{}", self.cohpcar_error, self.icohplist_error);
                if error.is_empty() || !error.contains(args) {
                    return Err(format!("error is {error:?}"));
                }
            }
            other => return Err(format!("unknown lobster command {other}")),
        }
        Ok(None)
    }
}
