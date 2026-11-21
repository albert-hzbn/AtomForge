//! LAMMPS data files (read_data; the C++ readLammpsData): box with tilts,
//! Masses (with element labels after '#'), Atom Type Labels, Atoms in the
//! atomic, charge, molecular or full style, Velocities.

use std::collections::BTreeMap;

use af_core::{elements, AtomSite, Structure};

use crate::text::{lines, stod, stoi, stoll, words};

/// The element of an atom type: a label (trusted when its mass agrees) or the
/// element whose standard mass is nearest (within 0.5 u).
pub fn type_element(label: &str, mass: f64) -> i32 {
    if !label.is_empty() {
        let mut symbol: Vec<char> = label.chars().collect();
        symbol[0] = symbol[0].to_ascii_uppercase();
        for c in symbol.iter_mut().skip(1) {
            *c = c.to_ascii_lowercase();
        }
        let symbol: String = symbol.into_iter().collect();
        // Open Babel's GetAtomicNum: symbols, and D/T for hydrogen isotopes.
        let z = match symbol.as_str() {
            "D" | "T" => 1,
            s => elements::atomic_number(s),
        };
        if z > 0 && (mass <= 0.0 || (elements::atomic_mass(z) - mass).abs() <= 1.0f64.max(0.03 * mass)) {
            return z;
        }
    }
    let mut best = 0;
    let mut best_difference = 0.5;
    if mass > 0.0 {
        for z in 1..=118 {
            let difference = (elements::atomic_mass(z) - mass).abs();
            if difference < best_difference {
                best = z;
                best_difference = difference;
            }
        }
    }
    best
}

fn trim(s: &str) -> &str {
    s.trim_matches(|c| c == ' ' || c == '\t' || c == '\r' || c == '\n')
}

pub fn read(text: &str, structure: &mut Structure) -> Result<(), String> {
    let all = lines(text, false);
    let mut atom_count: i64 = -1;
    let mut type_count: i32 = 0;
    let (mut lo, mut hi) = ([0.0f64; 3], [0.0f64; 3]);
    let (mut xy, mut xz, mut yz) = (0.0, 0.0, 0.0);
    let mut have = [false; 3];
    let mut masses: BTreeMap<i32, f64> = BTreeMap::new();
    let mut labels: BTreeMap<i32, String> = BTreeMap::new();
    let mut atom_rows: Vec<Vec<String>> = Vec::new();
    let mut velocity_rows: BTreeMap<i64, [f64; 3]> = BTreeMap::new();
    let mut style = String::new();
    let mut section = String::new();
    for line in all.iter().skip(1) {
        let (body_raw, comment) = match line.find('#') {
            Some(h) => (&line[..h], &line[h + 1..]),
            None => (*line, ""),
        };
        let body = trim(body_raw);
        if body.is_empty() {
            continue;
        }
        let w = words(body);
        let first = w[0].as_bytes()[0];
        if !first.is_ascii_digit() && first != b'-' && first != b'+' && first != b'.' {
            section = body.to_string();
            if section == "Atoms" {
                style = trim(comment).to_string();
            }
            continue;
        }
        let result: Result<(), String> = (|| {
            if section.is_empty() {
                if w.len() >= 2 && w[1] == "atoms" {
                    atom_count = stoll(w[0])?;
                } else if w.len() >= 3 && w[1] == "atom" && w[2] == "types" {
                    type_count = stoi(w[0])?;
                } else if w.len() >= 4 && w[2] == "xlo" && w[3] == "xhi" {
                    lo[0] = stod(w[0])?;
                    hi[0] = stod(w[1])?;
                    have[0] = true;
                } else if w.len() >= 4 && w[2] == "ylo" && w[3] == "yhi" {
                    lo[1] = stod(w[0])?;
                    hi[1] = stod(w[1])?;
                    have[1] = true;
                } else if w.len() >= 4 && w[2] == "zlo" && w[3] == "zhi" {
                    lo[2] = stod(w[0])?;
                    hi[2] = stod(w[1])?;
                    have[2] = true;
                } else if w.len() >= 6 && w[3] == "xy" && w[4] == "xz" && w[5] == "yz" {
                    xy = stod(w[0])?;
                    xz = stod(w[1])?;
                    yz = stod(w[2])?;
                }
            } else if section == "Masses" && w.len() >= 2 {
                let t = stoi(w[0])?;
                masses.insert(t, stod(w[1])?);
                if let Some(label) = words(comment).first() {
                    labels.entry(stoi(w[0])?).or_insert_with(|| label.to_string());
                }
            } else if section == "Atom Type Labels" && w.len() >= 2 {
                labels.insert(stoi(w[0])?, w[1].to_string());
            } else if section == "Atoms" {
                atom_rows.push(w.iter().map(|s| s.to_string()).collect());
            } else if section == "Velocities" && w.len() >= 4 {
                velocity_rows.insert(stoll(w[0])?, [stod(w[1])?, stod(w[2])?, stod(w[3])?]);
            }
            Ok(())
        })();
        if result.is_err() {
            return Err(format!("Cannot read this line of the LAMMPS data file: {}", trim(line)));
        }
    }
    if atom_count < 0 || !have.iter().all(|&h| h) {
        return Err("Not a LAMMPS data file: it needs the atom count and the xlo xhi, ylo yhi and zlo zhi box lines.".into());
    }
    if atom_rows.len() as i64 != atom_count {
        return Err(format!(
            "The LAMMPS data file declares {atom_count} atoms but its Atoms section lists {}.",
            atom_rows.len()
        ));
    }
    let (mut type_column, mut x_column) = (1usize, 2usize);
    let style_words = words(&style);
    let kind = style_words.first().copied().unwrap_or("");
    let columns = atom_rows.first().map_or(0, |r| r.len());
    match kind {
        "atomic" => {}
        "charge" => x_column = 3,
        "full" => {
            type_column = 2;
            x_column = 4;
        }
        "molecular" | "bond" | "angle" => {
            type_column = 2;
            x_column = 3;
        }
        "" => match columns {
            5 | 8 => {}
            6 | 9 => {
                let charge = atom_rows.iter().any(|row| {
                    let (v, used) = af_core::format::parse_prefix_f64(&row[2]).unwrap_or((0.0, 0));
                    used != row[2].len() || v != v.floor() || v < 1.0 || v > type_count.max(1) as f64
                });
                if !charge {
                    return Err("Cannot tell whether the LAMMPS atoms are in the charge or the molecular style; \
                                add it after the section name, e.g. \"Atoms # charge\"."
                        .into());
                }
                x_column = 3;
            }
            7 | 10 => {
                type_column = 2;
                x_column = 4;
            }
            _ => return Err("Cannot tell the LAMMPS atom style; add it after the section name, e.g. \"Atoms # atomic\".".into()),
        },
        other => return Err(format!("LAMMPS atom style '{other}' is not supported; use atomic, charge, molecular or full.")),
    }
    let mut element_of_type: BTreeMap<i32, i32> = BTreeMap::new();
    let mut unknown = String::new();
    for t in 1..=type_count.max(1) {
        let mass = masses.get(&t).copied().unwrap_or(0.0);
        let label = labels.get(&t).map(String::as_str).unwrap_or("");
        let z = type_element(label, mass);
        if z > 0 {
            element_of_type.insert(t, z);
        } else {
            if !unknown.is_empty() {
                unknown.push_str(", ");
            }
            unknown.push_str(&t.to_string());
        }
    }
    if !unknown.is_empty() {
        return Err(format!(
            "The element of LAMMPS atom type {unknown} is unknown: give its mass in a Masses section, or its symbol after the mass (\"1 63.546 # Cu\")."
        ));
    }
    let mut sites: Vec<(i64, AtomSite)> = Vec::with_capacity(atom_rows.len());
    for row in &atom_rows {
        if row.len() < x_column + 3 {
            return Err(format!(
                "A LAMMPS atom line has too few columns for the {} style.",
                if kind.is_empty() { "detected" } else { kind }
            ));
        }
        let parsed: Result<(i64, AtomSite), String> = (|| {
            let t = stoi(&row[type_column])?;
            let Some(&z) = element_of_type.get(&t) else {
                return Err(format!(
                    "LAMMPS atom {} has type {}, which is not among the declared types.",
                    row[0], row[type_column]
                ));
            };
            let mut site = AtomSite::new(elements::symbol(z), [stod(&row[x_column])?, stod(&row[x_column + 1])?, stod(&row[x_column + 2])?]);
            site.atomic_number = z;
            if kind == "charge" || (kind.is_empty() && x_column == 3 && type_column == 1) {
                site.charge = stod(&row[2])?;
            }
            if kind == "full" || (kind.is_empty() && x_column == 4) {
                site.charge = stod(&row[3])?;
            }
            let id = stoll(&row[0])?;
            if let Some(v) = velocity_rows.get(&id) {
                for k in 0..3 {
                    site.velocity[k] = v[k] / 1000.0;
                }
            }
            Ok((id, site))
        })();
        match parsed {
            Ok(entry) => sites.push(entry),
            Err(message) if message.starts_with("LAMMPS atom") => return Err(message),
            Err(_) => return Err(format!("Cannot read LAMMPS atom line starting {}.", row[0])),
        }
    }
    sites.sort_by_key(|(id, _)| *id);
    structure.atoms.extend(sites.into_iter().map(|(_, s)| s));
    structure.has_unit_cell = true;
    structure.cell_vectors = [[hi[0] - lo[0], 0.0, 0.0], [xy, hi[1] - lo[1], 0.0], [xz, yz, hi[2] - lo[2]]];
    structure.cell_offset = lo;
    Ok(())
}
