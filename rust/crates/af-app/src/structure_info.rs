//! View > Structure Info (the C++ ui/StructureInfoDialog): summary (atoms,
//! formula, cell, density), composition, symmetry from spglib (computed in
//! the background) and the atomic positions.

use std::collections::BTreeMap;
use std::sync::mpsc::{channel, Receiver};

use af_core::elements;
use af_core::glm::{self, cross, dot, length};
use af_core::structure::Structure;
use eframe::egui;

use crate::dialog_layout::{self as dl, Note};
use glam::{Mat3, Vec3};

pub type Rows = Vec<(String, String)>;

/// What the background worker computes (computeStructureInfo).
#[derive(Clone, Debug)]
pub struct Computed {
    pub element_counts: BTreeMap<i32, usize>,
    pub formula: String,
    pub symmetry: Result<(i32, String, String, String), String>,
}

pub fn element_counts(s: &Structure) -> BTreeMap<i32, usize> {
    let mut counts = BTreeMap::new();
    for a in &s.atoms {
        *counts.entry(a.atomic_number).or_insert(0) += 1;
    }
    counts
}

/// "Cu4 O2"-style formula over Z order (buildFormula).
pub fn formula(counts: &BTreeMap<i32, usize>) -> String {
    let mut parts = Vec::new();
    for (&z, &n) in counts {
        if z <= 0 || n == 0 {
            continue;
        }
        parts.push(if n > 1 { format!("{}{n}", elements::symbol(z)) } else { elements::symbol(z).to_string() });
    }
    parts.join(" ")
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// The formula reduced by the common divisor of the counts, e.g. Al2Cu3.
pub fn reduced_formula(counts: &BTreeMap<i32, usize>) -> String {
    let divisor = counts.iter().filter(|(&z, &n)| z > 0 && n > 0).fold(0, |d, (_, &n)| gcd(d, n));
    let mut text = String::new();
    for (&z, &n) in counts {
        if z <= 0 || n == 0 {
            continue;
        }
        text += elements::symbol(z);
        if n / divisor > 1 {
            text += &(n / divisor).to_string();
        }
    }
    text
}

fn lattice(s: &Structure) -> Option<(Mat3, Vec3)> {
    if !s.has_unit_cell {
        return None;
    }
    let cell = glm::cell_matrix(s);
    let origin = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
    (glm::determinant(&cell).abs() > 1e-10).then_some((cell, origin))
}

/// spglib's space group (analyzeSymmetryWithSpglib, symprec 1e-5): number,
/// international symbol, Hall symbol and point group.
pub fn symmetry(s: &Structure) -> Result<(i32, String, String, String), String> {
    if !s.has_unit_cell {
        return Err("No unit cell available.".into());
    }
    if s.atoms.is_empty() {
        return Err("No atoms available.".into());
    }
    let (cell, origin) = lattice(s).ok_or("Invalid unit cell (zero volume).")?;
    let inverse = glm::inverse(&cell);
    let mut positions = Vec::with_capacity(s.atoms.len());
    let mut types = Vec::with_capacity(s.atoms.len());
    for (i, a) in s.atoms.iter().enumerate() {
        let f = glm::mul(&inverse, Vec3::new(a.x as f32, a.y as f32, a.z as f32) - origin);
        let f = f - glm::floor(f);
        positions.push([f.x as f64, f.y as f64, f.z as f64]);
        types.push(if a.atomic_number > 0 { a.atomic_number } else { i as i32 + 1000 });
    }
    let d = af_science::spg::dataset(&s.cell_vectors, &positions, &types, 1e-5).ok_or("spglib could not identify symmetry (try larger symprec).")?;
    Ok((d.spacegroup_number, d.international_symbol, d.hall_symbol, d.pointgroup_symbol))
}

/// Identifies a structure's content (structureFingerprint): atom count,
/// cell and the first 64 atoms' elements and positions.
pub fn fingerprint(s: &Structure) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    s.atoms.len().hash(&mut h);
    s.has_unit_cell.hash(&mut h);
    if s.has_unit_cell {
        for row in &s.cell_vectors {
            for v in row {
                v.to_bits().hash(&mut h);
            }
        }
    }
    for a in s.atoms.iter().take(64) {
        a.atomic_number.hash(&mut h);
        for v in [a.x, a.y, a.z] {
            v.to_bits().hash(&mut h);
        }
    }
    h.finish()
}

pub fn compute(s: &Structure) -> Computed {
    let element_counts = element_counts(s);
    let formula = formula(&element_counts);
    Computed { element_counts, formula, symmetry: symmetry(s) }
}

fn angle(u: Vec3, v: Vec3) -> f32 {
    let (lu, lv) = (length(u), length(v));
    if lu <= 1e-10 || lv <= 1e-10 {
        return 0.0;
    }
    (dot(u, v) / (lu * lv)).clamp(-1.0, 1.0).acos() * (180.0 / std::f32::consts::PI)
}

pub fn summary_rows(s: &Structure, data: Option<&Computed>) -> Rows {
    let mut rows: Rows = vec![("Atoms".into(), s.atoms.len().to_string())];
    let fixed = s.atoms.iter().filter(|a| a.fixed != 0).count();
    if fixed > 0 {
        rows.push(("Fixed atoms".into(), format!("{fixed} (selective dynamics)")));
    }
    match data {
        Some(d) => {
            let reduced = reduced_formula(&d.element_counts);
            let full = d.formula.replace(' ', "");
            rows.push(("Formula".into(), if reduced == full { reduced } else { format!("{reduced}   (cell: {full})") }));
        }
        None => rows.push(("Formula".into(), "computing...".into())),
    }
    if !s.has_unit_cell {
        rows.push(("Cell".into(), "none (an isolated molecule or cluster)".into()));
        return rows;
    }
    let v = |k: usize| Vec3::new(s.cell_vectors[k][0] as f32, s.cell_vectors[k][1] as f32, s.cell_vectors[k][2] as f32);
    let (a, b, c) = (v(0), v(1), v(2));
    let volume = dot(a, cross(b, c)).abs();
    rows.push(("Cell lengths".into(), format!("a = {:.4} Å    b = {:.4} Å    c = {:.4} Å", length(a), length(b), length(c))));
    rows.push(("Cell angles".into(), format!("α = {:.2}°    β = {:.2}°    γ = {:.2}°", angle(b, c), angle(a, c), angle(a, b))));
    rows.push(("Volume".into(), format!("{volume:.3} Å³")));
    if let Some(d) = data {
        if volume > 1e-6 {
            let mass: f64 = d.element_counts.iter().filter(|(&z, _)| z > 0).map(|(&z, &n)| n as f64 * elements::atomic_mass(z)).sum();
            let volume = volume as f64;
            rows.push(("Density".into(), format!("{:.4} g/cm³    ({:.5} atoms/Å³)", mass * 1.66053906660 / volume, s.atoms.len() as f64 / volume)));
        }
    }
    rows
}

pub fn symmetry_rows(data: Option<&Computed>) -> Rows {
    match data.map(|d| &d.symmetry) {
        None => vec![("Space group".into(), "computing...".into())],
        Some(Err(e)) => vec![("Space group".into(), format!("unavailable ({e})"))],
        Some(Ok((number, international, hall, point))) => vec![
            ("Space group".into(), format!("{international} (No. {number})")),
            ("Point group".into(), point.clone()),
            ("Hall symbol".into(), hall.clone()),
        ],
    }
}

/// The text of Copy Summary.
pub fn summary_text(summary: &Rows, symmetry: &Rows, data: Option<&Computed>) -> String {
    let mut text = String::new();
    for (label, value) in summary.iter().chain(symmetry.iter()) {
        text += &format!("{label}: {value}\n");
    }
    if let Some(d) = data {
        for (&z, &n) in &d.element_counts {
            text += &format!("{}: {n} atoms\n", elements::symbol(z));
        }
    }
    text
}

#[derive(Default)]
pub struct StructureInfoWindow {
    pub open: bool,
    /// The structure revision the data is for and the data.
    cached: Option<(u64, Computed)>,
    pending: Option<(u64, Receiver<Computed>)>,
}

const MAX_POSITION_ROWS: usize = 500;

impl StructureInfoWindow {
    /// The data for this revision, once computed (UI scripts wait for it).
    pub fn data(&self, revision: u64) -> Option<&Computed> {
        self.cached.as_ref().filter(|(r, _)| *r == revision).map(|(_, d)| d)
    }

    /// Starts (or finishes) the background computation for a revision.
    pub fn update(&mut self, s: &Structure, revision: u64) {
        if let Some((r, rx)) = &self.pending {
            if let Ok(data) = rx.try_recv() {
                self.cached = Some((*r, data));
                self.pending = None;
            }
        }
        let current = self.cached.as_ref().is_some_and(|(r, _)| *r == revision);
        if self.open && !current && self.pending.as_ref().is_none_or(|(r, _)| *r != revision) {
            let (tx, rx) = channel();
            let snapshot = s.clone();
            std::thread::spawn(move || {
                let _ = tx.send(compute(&snapshot));
            });
            self.pending = Some((revision, rx));
        }
    }

    pub fn show(&mut self, ctx: &egui::Context, s: &Structure, revision: u64) {
        self.update(s, revision);
        if !self.open {
            return;
        }
        if self.pending.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }
        let data = self.data(revision).cloned();
        let summary = summary_rows(s, data.as_ref());
        let symmetry = symmetry_rows(data.as_ref());
        let mut open = self.open;
        // A fixed width and a height that follows the content.
        const WIDTH: f32 = 760.0 - 32.0;
        crate::theme::modal_window("Structure Info").open(&mut open).resizable(false).min_width(WIDTH).max_width(WIDTH).show(ctx, |ui| {
            ui.set_width(WIDTH);
            let p = crate::theme::palette(ui.visuals().dark_mode);
            // Labels muted in a column of 7.5 font sizes, values after them (drawRows).
            let rows = |ui: &mut egui::Ui, rows: &Rows| {
                ui.scope(|ui| {
                    ui.spacing_mut().item_spacing.y = 6.0;
                    for (label, value) in rows {
                        ui.horizontal(|ui| {
                            let (r, _) = ui.allocate_exact_size(egui::vec2(18.0 * 7.5, 22.0), egui::Sense::hover());
                            ui.painter().text(egui::pos2(r.left(), r.center().y), egui::Align2::LEFT_CENTER, label, egui::TextStyle::Body.resolve(ui.style()), p.text_muted);
                            ui.label(value);
                        });
                    }
                });
            };
            dl::section(ui, "Summary");
            rows(ui, &summary);
            dl::section(ui, "Composition");
            match &data {
                Some(d) => {
                    let counts: Vec<(i32, usize)> = d.element_counts.iter().map(|(&z, &n)| (z, n)).collect();
                    let total = s.atoms.len();
                    dl::table(ui, &["Element", "Z", "Count", "Fraction"], &[0.25, 0.25, 0.25], counts.len(), |ui, r, k| {
                        let (z, n) = counts[r];
                        let text = match k {
                            0 => elements::symbol(z).to_string(),
                            1 => z.to_string(),
                            2 => n.to_string(),
                            _ => format!("{:.2} %", if total == 0 { 0.0 } else { 100.0 * n as f64 / total as f64 }),
                        };
                        ui.label(text);
                    });
                }
                None => {
                    dl::muted(ui, "Counting elements...");
                }
            }
            dl::section(ui, "Symmetry");
            rows(ui, &symmetry);
            dl::section(ui, "Atomic Positions");
            if s.atoms.is_empty() {
                dl::muted(ui, "No atoms.");
            } else if s.atoms.len() > MAX_POSITION_ROWS {
                dl::note(
                    ui,
                    Note::Info,
                    &format!("The position table lists up to {MAX_POSITION_ROWS} atoms; this structure has {}. Save it (File > Save As) as an XYZ or CIF file to read every position.", s.atoms.len()),
                );
            } else {
                let cell = glm::cell_matrices(s);
                let origin = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
                let has_lattice = lattice(s).is_some();
                let shown = |v: f64| if v.abs() < 5e-6 { 0.0 } else { v };
                let values: Vec<[f64; 6]> = s
                    .atoms
                    .iter()
                    .map(|a| {
                        let f = match (&cell, has_lattice) {
                            (Some((_, inv)), true) => glm::mul(inv, Vec3::new(a.x as f32, a.y as f32, a.z as f32) - origin),
                            _ => Vec3::ZERO,
                        };
                        [a.x, a.y, a.z, f.x as f64, f.y as f64, f.z as f64]
                    })
                    .collect();
                // Stretch weights 0.6, 0.6 and 1 (as the C++ table): fractions of the width.
                let w = [0.6 / 7.2, 0.6 / 7.2, 1.0 / 7.2, 1.0 / 7.2, 1.0 / 7.2, 1.0 / 7.2, 1.0 / 7.2];
                let headers = ["#", "Element", "x (Å)", "y (Å)", "z (Å)", "u", "v", "w"];
                // The header, then up to ten rows that scroll under it.
                dl::table(ui, &headers, &w, 0, |_, _, _| {});
                let row_h = crate::theme::FRAME_HEIGHT - 6.0;
                egui::ScrollArea::vertical().id_salt("positions").max_height(row_h * s.atoms.len().min(10) as f32 + 2.0).show(ui, |ui| {
                    dl::table(ui, &[""; 8], &w, s.atoms.len(), |ui, r, k| {
                        let text = match k {
                            0 => r.to_string(),
                            1 => s.atoms[r].symbol.clone(),
                            _ => format!("{:.5}", shown(values[r][k - 2])),
                        };
                        ui.label(text);
                    });
                });
            }
            match dl::footer(ui, "Close", Some("Copy Summary"), true) {
                2 => ui.ctx().copy_text(summary_text(&summary, &symmetry, data.as_ref())),
                1 => self.open = false,
                _ => {}
            }
        });
        self.open &= open;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::structure::AtomSite;

    fn atom(symbol: &str, p: [f64; 3]) -> AtomSite {
        AtomSite { symbol: symbol.into(), atomic_number: elements::atomic_number(symbol), x: p[0], y: p[1], z: p[2], ..Default::default() }
    }

    #[test]
    fn rock_salt_summary() {
        let a = 5.64;
        let mut s = Structure { cell_vectors: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], has_unit_cell: true, ..Default::default() };
        for f in [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]] {
            s.atoms.push(atom("Na", [f[0] * a, f[1] * a, f[2] * a]));
            s.atoms.push(atom("Cl", [(f[0] + 0.5) * a, f[1] * a, f[2] * a]));
        }
        let data = compute(&s);
        assert_eq!(data.formula, "Na4 Cl4");
        assert_eq!(reduced_formula(&data.element_counts), "NaCl");
        let rows = summary_rows(&s, Some(&data));
        let get = |k: &str| rows.iter().find(|(l, _)| l == k).map(|(_, v)| v.clone()).unwrap();
        assert_eq!(get("Formula"), "NaCl   (cell: Na4Cl4)");
        assert_eq!(get("Cell lengths"), "a = 5.6400 Å    b = 5.6400 Å    c = 5.6400 Å");
        assert_eq!(get("Cell angles"), "α = 90.00°    β = 90.00°    γ = 90.00°");
        assert_eq!(get("Volume"), "179.406 Å³");
        // 4 NaCl (58.44 u) in 179.4 Å^3: 2.164 g/cm3.
        assert!(get("Density").starts_with("2.16"), "{}", get("Density"));
        let sym = symmetry_rows(Some(&data));
        assert_eq!(sym[0].1, "Fm-3m (No. 225)");
        assert_eq!(sym[1].1, "m-3m");
        assert_eq!(sym[2].1, "-F 4 2 3");
    }

    #[test]
    fn molecule_summary() {
        let s = Structure { atoms: vec![atom("O", [0.0; 3]), atom("H", [0.96, 0.0, 0.0]), atom("H", [-0.24, 0.93, 0.0])], ..Default::default() };
        let data = compute(&s);
        assert_eq!(data.formula, "H2 O");
        let rows = summary_rows(&s, Some(&data));
        assert_eq!(rows[1].1, "H2O");
        assert_eq!(rows[2], ("Cell".to_string(), "none (an isolated molecule or cluster)".to_string()));
        assert_eq!(symmetry_rows(Some(&data))[0].1, "unavailable (No unit cell available.)");
    }
}
