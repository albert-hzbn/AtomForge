//! Molecular formats the C++ read and wrote through Open Babel, now native:
//! PDB (with the CRYST1 cell), MDL MOL/SDF (V2000 and V3000; the first
//! record), Tripos MOL2 and Gaussian input (GJF, Cartesian with Tv lattice
//! vectors). Writers follow Open Babel's layouts; bonds come from covalent
//! radii (1.2 x the sum) for molecules of up to 500 atoms without a cell.

use af_core::{elements, AtomSite, Structure};

use crate::cif::cell_from_parameters;
use crate::text::{lines, words};

fn normal_symbol(text: &str) -> String {
    let mut chars = text.chars();
    match chars.next() {
        Some(first) => std::iter::once(first.to_ascii_uppercase()).chain(chars.map(|c| c.to_ascii_lowercase())).collect(),
        None => String::new(),
    }
}

fn element_number(text: &str) -> i32 {
    match normal_symbol(text).as_str() {
        "D" | "T" => 1,
        s => elements::atomic_number(s),
    }
}

fn field(line: &str, start: usize, width: usize) -> &str {
    let bytes = line.as_bytes();
    if start >= bytes.len() {
        return "";
    }
    let end = (start + width).min(bytes.len());
    std::str::from_utf8(&bytes[start..end]).unwrap_or("")
}

// ---------------------------------------------------------------- PDB

/// The element of a PDB atom: columns 77-78, else from the atom name.
fn pdb_element(line: &str) -> i32 {
    let element: String = field(line, 76, 2).chars().filter(|c| c.is_ascii_alphabetic()).collect();
    if !element.is_empty() {
        let z = element_number(&element);
        if z > 0 {
            return z;
        }
    }
    let name = field(line, 12, 4);
    let b = name.as_bytes();
    if b.is_empty() {
        return 0;
    }
    // A name starting with a space or digit has a one-letter element in column 14.
    if b[0] == b' ' || b[0].is_ascii_digit() {
        let letters: String = name[1..].chars().take_while(|c| c.is_ascii_alphabetic()).take(1).collect();
        return element_number(&letters);
    }
    let two: String = name.chars().take(2).filter(|c| c.is_ascii_alphabetic()).collect();
    if two.len() == 2 {
        let z = element_number(&two);
        if z > 0 {
            return z;
        }
    }
    element_number(&name[..1])
}

/// The CRYST1 cell at the precision written.
pub fn pdb_cell(text: &str) -> Option<[[f64; 3]; 3]> {
    for line in lines(text, true) {
        if !line.starts_with("CRYST1") {
            continue;
        }
        let widths = [(6, 9), (15, 9), (24, 9), (33, 7), (40, 7), (47, 7)];
        let mut values = [0.0; 6];
        for (k, &(start, width)) in widths.iter().enumerate() {
            if line.len() < start + 1 {
                return None;
            }
            values[k] = af_core::format::parse_prefix_f64(field(line, start, width))?.0;
        }
        if values[0] <= 0.0 || values[1] <= 0.0 || values[2] <= 0.0 {
            return None;
        }
        let cell = cell_from_parameters([values[0], values[1], values[2]], [values[3], values[4], values[5]])?;
        return Some([cell.col(0).to_array(), cell.col(1).to_array(), cell.col(2).to_array()]);
    }
    None
}

pub fn read_pdb(text: &str, structure: &mut Structure) -> Result<(), String> {
    for line in lines(text, true) {
        if line.starts_with("ENDMDL") || line.starts_with("END ") || line == "END" {
            if !structure.atoms.is_empty() {
                break;
            }
            continue;
        }
        if !line.starts_with("ATOM") && !line.starts_with("HETATM") {
            continue;
        }
        let coordinate = |start: usize| af_core::format::parse_prefix_f64(field(line, start, 8)).map(|v| v.0);
        let (Some(x), Some(y), Some(z)) = (coordinate(30), coordinate(38), coordinate(46)) else {
            return Err("Failed to load file. The file may be corrupted or unreadable.".into());
        };
        let number = pdb_element(line);
        if number == 0 {
            return Err("Failed to load file. The file may be corrupted or unreadable.".into());
        }
        let mut site = AtomSite::new(elements::symbol(number), [x, y, z]);
        site.atomic_number = number;
        structure.atoms.push(site);
    }
    if let Some(cell) = pdb_cell(text) {
        structure.has_unit_cell = true;
        structure.cell_vectors = cell;
    }
    Ok(())
}

// ---------------------------------------------------------------- MOL / SDF

pub fn read_mol(text: &str, structure: &mut Structure) -> Result<(), String> {
    let all = lines(text, true);
    let bad = || "Failed to load file. The file may be corrupted or unreadable.".to_string();
    let counts = all.get(3).ok_or_else(bad)?;
    if counts.contains("V3000") {
        let mut inside = false;
        for line in &all[4..] {
            let t = line.trim();
            if t.starts_with("M  V30 BEGIN ATOM") {
                inside = true;
                continue;
            }
            if t.starts_with("M  V30 END ATOM") {
                break;
            }
            if inside && t.starts_with("M  V30") {
                let w = words(t);
                // M V30 index symbol x y z ...
                if w.len() < 7 {
                    return Err(bad());
                }
                let number = element_number(w[3]);
                let coordinate = |s: &str| s.parse::<f64>().map_err(|_| bad());
                let mut site = AtomSite::new(elements::symbol(number.max(0)), [coordinate(w[4])?, coordinate(w[5])?, coordinate(w[6])?]);
                site.atomic_number = number;
                structure.atoms.push(site);
            }
        }
        return Ok(());
    }
    let count = crate::text::strtol(field(counts, 0, 3)).0;
    if count < 0 {
        return Err(bad());
    }
    for k in 0..count as usize {
        let line = all.get(4 + k).ok_or_else(bad)?;
        let coordinate = |start: usize| af_core::format::parse_prefix_f64(field(line, start, 10)).map(|v| v.0).ok_or_else(bad);
        let symbol = field(line, 31, 3).trim();
        let number = element_number(symbol);
        let mut site = AtomSite::new(elements::symbol(number.max(0)), [coordinate(0)?, coordinate(10)?, coordinate(20)?]);
        site.atomic_number = number;
        // Charge column (0, 1-7: +3, +2, +1, doublet, -1, -2, -3) is not kept.
        structure.atoms.push(site);
    }
    Ok(())
}

// ---------------------------------------------------------------- MOL2

pub fn read_mol2(text: &str, structure: &mut Structure) -> Result<(), String> {
    let mut inside = false;
    for line in lines(text, true) {
        if line.starts_with("@<TRIPOS>") {
            if inside {
                break;
            }
            inside = line.starts_with("@<TRIPOS>ATOM");
            continue;
        }
        if !inside {
            continue;
        }
        let w = words(line);
        if w.len() < 6 {
            continue;
        }
        let kind = w[5];
        let element = &kind[..kind.find('.').unwrap_or(kind.len())];
        let mut number = element_number(element);
        if number == 0 {
            // Types such as "Du" (dummy) or "LP": fall back to the atom name.
            number = element_number(&w[1].chars().take_while(|c| c.is_ascii_alphabetic()).collect::<String>());
        }
        let coordinate = |s: &str| s.parse::<f64>().map_err(|_| "Failed to load file. The file may be corrupted or unreadable.".to_string());
        let mut site = AtomSite::new(elements::symbol(number), [coordinate(w[2])?, coordinate(w[3])?, coordinate(w[4])?]);
        site.atomic_number = number;
        // Partial charges (column 9) are not kept, as in the C++ build.
        structure.atoms.push(site);
    }
    Ok(())
}

// ---------------------------------------------------------------- Gaussian input

pub fn read_gjf(text: &str, structure: &mut Structure) -> Result<(), String> {
    let all = lines(text, true);
    let mut i = 0;
    // Link 0 commands, then the route section up to a blank line.
    while i < all.len() && all[i].trim_start().starts_with('%') {
        i += 1;
    }
    while i < all.len() && !all[i].trim().is_empty() {
        i += 1;
    }
    // After the title, the charge and multiplicity: a line of integers only
    // (titles may be blank, as Open Babel writes them).
    let integers = |line: &str| {
        let w = words(line);
        w.len() >= 2 && w.len() % 2 == 0 && w.iter().all(|x| x.parse::<i64>().is_ok())
    };
    while i < all.len() && !integers(all[i]) {
        i += 1;
    }
    if i >= all.len() {
        return Err("The Gaussian input has no charge and multiplicity line.".into());
    }
    i += 1;
    let mut lattice = Vec::new();
    while i < all.len() && !all[i].trim().is_empty() {
        let line = all[i].replace(',', " ");
        let w = words(&line);
        i += 1;
        if w.is_empty() {
            continue;
        }
        let label = w[0];
        let symbol: String = label.chars().take_while(|c| c.is_ascii_alphabetic()).collect();
        let numbers: Vec<f64> = w[1..].iter().filter_map(|s| crate::text::number_full(s)).collect();
        // "El x y z" or "El freeze x y z" (Cartesian); Z-matrix lines are not read here.
        let xyz = match (w.len(), numbers.len()) {
            (4, 3) => [numbers[0], numbers[1], numbers[2]],
            (5, 4) => [numbers[1], numbers[2], numbers[3]],
            _ => return Err("Only Cartesian Gaussian input is read (an element and x y z per line).".into()),
        };
        if symbol.eq_ignore_ascii_case("tv") {
            lattice.push(xyz);
            continue;
        }
        let number = if symbol.is_empty() { crate::text::atoi(label) } else { element_number(&symbol) };
        if !(1..=118).contains(&number) {
            return Err(format!("Unknown element '{label}' in the Gaussian input."));
        }
        let mut site = AtomSite::new(elements::symbol(number), xyz);
        site.atomic_number = number;
        if w.len() == 5 && numbers[0] == -1.0 {
            site.fixed = 7;
        }
        structure.atoms.push(site);
    }
    if !lattice.is_empty() {
        // One, two or three translation vectors: periodic along them.
        structure.has_unit_cell = true;
        let mut cell = [[0.0; 3]; 3];
        let mut pbc = [false; 3];
        for (k, v) in lattice.iter().take(3).enumerate() {
            cell[k] = *v;
            pbc[k] = true;
        }
        // Missing vectors (a slab or a chain): along the remaining axes, long
        // enough to hold the atoms with 10 A to spare, the cell starting below them.
        let mut offset = [0.0; 3];
        for k in lattice.len().min(3)..3 {
            let (low, high) = structure.atoms.iter().fold((f64::MAX, f64::MIN), |(l, h), a| {
                let v = [a.x, a.y, a.z][k];
                (l.min(v), h.max(v))
            });
            let (low, high) = if structure.atoms.is_empty() { (0.0, 0.0) } else { (low, high) };
            cell[k] = [0.0; 3];
            cell[k][k] = high - low + 20.0;
            offset[k] = low - 10.0;
        }
        structure.cell_vectors = cell;
        structure.cell_offset = offset;
        structure.pbc = pbc;
    }
    Ok(())
}

// ---------------------------------------------------------------- writers

/// Bonds by covalent radii (1.2 x the sum, at least 0.4 A), i < j.
pub fn guess_bonds(s: &Structure) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    if s.has_unit_cell || s.atoms.len() > 500 {
        return out;
    }
    for i in 0..s.atoms.len() {
        for j in i + 1..s.atoms.len() {
            let (a, b) = (&s.atoms[i], &s.atoms[j]);
            let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
            let limit = 1.2 * (elements::covalent_radius(a.atomic_number) + elements::covalent_radius(b.atomic_number)) as f64;
            if d > 0.4 && d <= limit {
                out.push((i, j));
            }
        }
    }
    out
}

fn symbol_of(a: &AtomSite) -> String {
    if a.symbol.is_empty() {
        elements::symbol(a.atomic_number).to_string()
    } else {
        a.symbol.clone()
    }
}

fn clean(v: f64) -> f64 {
    if v.abs() < 5e-6 {
        0.0
    } else {
        v
    }
}

pub fn write_pdb(s: &Structure) -> String {
    let mut out = String::from("COMPND    UNNAMED\nAUTHOR    GENERATED BY ATOMFORGE\n");
    if s.has_unit_cell {
        let [a, b, c] = s.cell_vectors.map(glam::DVec3::from);
        let angle = |u: glam::DVec3, v: glam::DVec3| (u.dot(v) / (u.length() * v.length())).clamp(-1.0, 1.0).acos().to_degrees();
        out.push_str(&format!(
            "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} P1          1\n",
            a.length(),
            b.length(),
            c.length(),
            angle(b, c),
            angle(a, c),
            angle(a, b)
        ));
    }
    for (i, a) in s.atoms.iter().enumerate() {
        let symbol = symbol_of(a);
        let name = if symbol.len() == 1 { format!(" {symbol:<3}") } else { format!("{:<4}", symbol.to_ascii_uppercase()) };
        out.push_str(&format!(
            "HETATM{:>5} {} UNL     1    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00          {:>2}  \n",
            (i + 1) % 100000,
            name,
            clean(a.x),
            clean(a.y),
            clean(a.z),
            symbol
        ));
    }
    let bonds = guess_bonds(s);
    let mut partners: Vec<Vec<usize>> = vec![Vec::new(); s.atoms.len()];
    for &(i, j) in &bonds {
        partners[i].push(j);
        partners[j].push(i);
    }
    for (i, list) in partners.iter().enumerate() {
        for chunk in list.chunks(4) {
            let mut line = format!("CONECT{:>5}", i + 1);
            for &j in chunk {
                line.push_str(&format!("{:>5}", j + 1));
            }
            out.push_str(&format!("{line:<70}\n"));
        }
    }
    out.push_str(&format!("MASTER        0    0    0    0    0    0    0    0{:>5}    0{:>5}    0\nEND\n", s.atoms.len(), s.atoms.len()));
    out
}

pub fn write_mol(s: &Structure, sdf: bool) -> String {
    let bonds = guess_bonds(s);
    let mut out = String::from("\n AtomForge          3D\n\n");
    out.push_str(&format!("{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n", s.atoms.len().min(999), bonds.len().min(999)));
    for a in &s.atoms {
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            clean(a.x),
            clean(a.y),
            clean(a.z),
            symbol_of(a)
        ));
    }
    for (i, j) in bonds {
        out.push_str(&format!("{:>3}{:>3}  1  0  0  0  0\n", i + 1, j + 1));
    }
    out.push_str("M  END\n");
    if sdf {
        out.push_str("$$$$\n");
    }
    out
}

pub fn write_mol2(s: &Structure) -> String {
    let bonds = guess_bonds(s);
    let mut out = format!("@<TRIPOS>MOLECULE\n*****\n {} {} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n", s.atoms.len(), bonds.len());
    for (i, a) in s.atoms.iter().enumerate() {
        let symbol = symbol_of(a);
        out.push_str(&format!(
            "{:>7} {:<8}{:>10.4}{:>10.4}{:>10.4} {:<7}{:>2}  UNL1  {:>12.4}\n",
            i + 1,
            symbol.to_ascii_uppercase(),
            clean(a.x),
            clean(a.y),
            clean(a.z),
            symbol,
            1,
            a.charge
        ));
    }
    out.push_str("@<TRIPOS>BOND\n");
    for (k, (i, j)) in bonds.iter().enumerate() {
        out.push_str(&format!("{:>6}{:>6}{:>6}    1\n", k + 1, i + 1, j + 1));
    }
    if s.has_unit_cell {
        let [a, b, c] = s.cell_vectors.map(glam::DVec3::from);
        let angle = |u: glam::DVec3, v: glam::DVec3| (u.dot(v) / (u.length() * v.length())).clamp(-1.0, 1.0).acos().to_degrees();
        out.push_str(&format!(
            "@<TRIPOS>CRYSIN\n{:.4} {:.4} {:.4} {:.4} {:.4} {:.4} 1 1\n",
            a.length(),
            b.length(),
            c.length(),
            angle(b, c),
            angle(a, c),
            angle(a, b)
        ));
    }
    out
}

pub fn write_gjf(s: &Structure) -> String {
    let mut out = String::from("!Put Keywords Here, check Charge and Multiplicity.\n#\n\n AtomForge structure\n\n0  1\n");
    for a in &s.atoms {
        out.push_str(&format!("{:<3}{:>17.5}{:>16.5}{:>16.5}\n", symbol_of(a), clean(a.x), clean(a.y), clean(a.z)));
    }
    if s.has_unit_cell {
        for k in 0..3 {
            if s.pbc[k] {
                let v = s.cell_vectors[k];
                out.push_str(&format!("{:<3}{:>17.5}{:>16.5}{:>16.5}\n", "Tv", v[0], v[1], v[2]));
            }
        }
    }
    out.push_str("\n\n");
    out
}
