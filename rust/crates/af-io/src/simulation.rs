//! Inputs of simulation codes (the C++ io/SimulationFormats): Quantum
//! ESPRESSO pw.x (pwi), CASTEP .cell, DFTB+ .gen, XCrySDen .xsf (read and
//! write) and a CP2K &SUBSYS section (written).

use std::collections::BTreeMap;

use af_core::format::format_g;
use af_core::{elements, AtomSite, Structure};

use crate::text::{lines, stod_used, strtol, words};

const BOHR: f64 = 0.529177210903;

fn canonical(format: &str) -> String {
    let f = format.to_ascii_lowercase();
    if f == "pwscf" || f == "pwi" || f == "in" {
        "pwi".into()
    } else {
        f
    }
}

pub fn is_simulation_format(format: &str) -> bool {
    matches!(canonical(format).as_str(), "pwi" | "cell" | "gen" | "xsf" | "cp2k")
}

pub fn can_read_simulation_format(format: &str) -> bool {
    matches!(canonical(format).as_str(), "pwi" | "cell" | "gen" | "xsf")
}

/// An element symbol from a label such as "Fe1", "fe" or "Fe_up".
fn element_from(label: &str) -> String {
    let mut letters = String::new();
    for c in label.chars() {
        if !c.is_ascii_alphabetic() {
            break;
        }
        letters.push(c);
        if letters.len() == 2 {
            break;
        }
    }
    if letters.is_empty() {
        return String::new();
    }
    let mut chars: Vec<char> = letters.chars().collect();
    chars[0] = chars[0].to_ascii_uppercase();
    if chars.len() == 2 {
        chars[1] = chars[1].to_ascii_lowercase();
        let two: String = chars.iter().collect();
        if elements::atomic_number(&two) > 0 {
            return two;
        }
        chars.truncate(1);
    }
    let one: String = chars.iter().collect();
    if elements::atomic_number(&one) > 0 {
        one
    } else {
        String::new()
    }
}

fn atom(label: &str, position: [f64; 3]) -> Result<AtomSite, String> {
    let mut symbol = element_from(label);
    if symbol.is_empty() {
        let (z, used) = strtol(label);
        if used == label.len() && used > 0 && (1..=118).contains(&z) {
            symbol = elements::symbol(z as i32).to_string();
        }
    }
    if symbol.is_empty() {
        return Err(format!("Unknown element '{label}'"));
    }
    Ok(AtomSite::new(&symbol, position))
}

fn species(s: &Structure) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for a in &s.atoms {
        if !out.contains(&a.symbol) {
            out.push(a.symbol.clone());
        }
    }
    out
}

fn from_fractional(s: &Structure, f: [f64; 3]) -> [f64; 3] {
    let mut x = [0.0; 3];
    for r in 0..3 {
        for c in 0..3 {
            x[c] += f[r] * s.cell_vectors[r][c];
        }
    }
    x
}

fn length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// `std::setw(18) << std::fixed << std::setprecision(10)`.
fn number(out: &mut String, value: f64) {
    out.push_str(&format!("{:>18.10}", value));
}

fn vector_line(out: &mut String, v: [f64; 3]) {
    out.push_str("  ");
    for x in v {
        number(out, x);
    }
    out.push('\n');
}

fn left(text: &str, width: usize) -> String {
    format!("{text:<width$}")
}

fn write_pwi(s: &Structure) -> String {
    struct Kind {
        label: String,
        element: String,
        magmom: f64,
    }
    let magnetic = s.atoms.iter().any(|a| a.magmom != 0.0);
    let mut kinds: Vec<Kind> = Vec::new();
    let mut kind_of = Vec::new();
    for a in &s.atoms {
        let k = kinds.iter().position(|k| k.element == a.symbol && (!magnetic || k.magmom == a.magmom)).unwrap_or_else(|| {
            kinds.push(Kind { label: a.symbol.clone(), element: a.symbol.clone(), magmom: if magnetic { a.magmom } else { 0.0 } });
            kinds.len() - 1
        });
        kind_of.push(k);
    }
    let mut per_element: BTreeMap<String, i32> = BTreeMap::new();
    for k in &kinds {
        *per_element.entry(k.element.clone()).or_default() += 1;
    }
    let mut numbered: BTreeMap<String, i32> = BTreeMap::new();
    for k in kinds.iter_mut() {
        if per_element[&k.element] > 1 {
            let n = numbered.entry(k.element.clone()).or_default();
            *n += 1;
            k.label = format!("{}{}", k.element, n);
        }
    }
    let mut out = String::new();
    out.push_str("&CONTROL\n  calculation = 'scf'\n  prefix = 'atomforge'\n  pseudo_dir = './'\n  outdir = './out'\n/\n");
    out.push_str(&format!("&SYSTEM\n  ibrav = 0\n  nat = {}\n  ntyp = {}\n  ecutwfc = 40\n  ecutrho = 320\n", s.atoms.len(), kinds.len()));
    if magnetic {
        out.push_str("  nspin = 2\n");
        for (k, kind) in kinds.iter().enumerate() {
            // The stream is still in its default (general) format here: precision 4 is %.4g.
            out.push_str(&format!("  starting_magnetization({}) = {}\n", k + 1, format_g((kind.magmom / 5.0).clamp(-1.0, 1.0), 4)));
        }
    }
    out.push_str("/\n&ELECTRONS\n/\n");
    let fixed = s.atoms.iter().any(|a| a.fixed != 0);
    if fixed {
        out.push_str("&IONS\n/\n");
    }
    out.push_str("ATOMIC_SPECIES\n");
    for k in &kinds {
        out.push_str(&format!("  {} {:.4} {}.UPF\n", left(&k.label, 4), elements::atomic_mass(elements::atomic_number(&k.element)), k.element));
    }
    out.push_str("CELL_PARAMETERS angstrom\n");
    for row in s.cell_vectors {
        vector_line(&mut out, row);
    }
    out.push_str("ATOMIC_POSITIONS angstrom\n");
    for (i, a) in s.atoms.iter().enumerate() {
        out.push_str("  ");
        out.push_str(&left(&kinds[kind_of[i]].label, 4));
        number(&mut out, a.x);
        number(&mut out, a.y);
        number(&mut out, a.z);
        if fixed {
            let free = |bit: u8| if a.fixed & bit != 0 { 0 } else { 1 };
            out.push_str(&format!(" {} {} {}", free(1), free(2), free(4)));
        }
        out.push('\n');
    }
    out.push_str("K_POINTS automatic\n ");
    for k in 0..3 {
        let len = length(s.cell_vectors[k]);
        let n = if s.pbc[k] { 1.max((30.0 / len).ceil() as i32) } else { 1 };
        out.push_str(&format!(" {n}"));
    }
    out.push_str(" 0 0 0\n");
    out
}

fn write_cell(s: &Structure) -> String {
    let mut out = String::from("%BLOCK LATTICE_CART\nANG\n");
    for row in s.cell_vectors {
        vector_line(&mut out, row);
    }
    out.push_str("%ENDBLOCK LATTICE_CART\n\n%BLOCK POSITIONS_ABS\nANG\n");
    let magnetic = s.atoms.iter().any(|a| a.magmom != 0.0);
    for a in &s.atoms {
        out.push_str("  ");
        out.push_str(&left(&a.symbol, 3));
        number(&mut out, a.x);
        number(&mut out, a.y);
        number(&mut out, a.z);
        if magnetic {
            // std::fixed is still set by number(): precision 6 is %.6f.
            out.push_str(&format!(" SPIN={:.6}", a.magmom));
        }
        out.push('\n');
    }
    out.push_str("%ENDBLOCK POSITIONS_ABS\n");
    let mut lines_out = Vec::new();
    let mut counter: BTreeMap<String, i32> = BTreeMap::new();
    let mut constraint = 0;
    for a in &s.atoms {
        let index = {
            let c = counter.entry(a.symbol.clone()).or_default();
            *c += 1;
            *c
        };
        for k in 0..3 {
            if a.fixed & (1 << k) != 0 {
                constraint += 1;
                lines_out.push(format!("  {constraint} {} {index} {} {} {}", a.symbol, (k == 0) as i32, (k == 1) as i32, (k == 2) as i32));
            }
        }
    }
    if !lines_out.is_empty() {
        out.push_str("\n%BLOCK IONIC_CONSTRAINTS\n");
        for l in lines_out {
            out.push_str(&l);
            out.push('\n');
        }
        out.push_str("%ENDBLOCK IONIC_CONSTRAINTS\n");
    }
    out
}

fn write_gen(s: &Structure) -> String {
    let kinds = species(s);
    let mut out = format!("{} {}\n", s.atoms.len(), if s.has_unit_cell { 'S' } else { 'C' });
    for k in &kinds {
        out.push_str(k);
        out.push(' ');
    }
    out.push('\n');
    for (i, a) in s.atoms.iter().enumerate() {
        let kind = kinds.iter().position(|k| *k == a.symbol).unwrap() + 1;
        out.push_str(&format!("{:>6} {:>3}", i + 1, kind));
        number(&mut out, a.x);
        number(&mut out, a.y);
        number(&mut out, a.z);
        out.push('\n');
    }
    if s.has_unit_cell {
        vector_line(&mut out, s.cell_offset);
        for row in s.cell_vectors {
            vector_line(&mut out, row);
        }
    }
    out
}

fn write_xsf(s: &Structure) -> String {
    let periodic = (0..3).filter(|&k| s.has_unit_cell && s.pbc[k]).count();
    let mut out = String::new();
    if !s.has_unit_cell {
        out.push_str("ATOMS\n");
    } else {
        let kind = if periodic == 3 {
            "CRYSTAL"
        } else if periodic == 2 && s.pbc[0] && s.pbc[1] {
            "SLAB"
        } else if periodic == 1 && s.pbc[0] {
            "POLYMER"
        } else if periodic == 0 {
            "MOLECULE"
        } else {
            "CRYSTAL"
        };
        out.push_str(kind);
        out.push_str("\nPRIMVEC\n");
        for row in s.cell_vectors {
            vector_line(&mut out, row);
        }
        out.push_str(&format!("PRIMCOORD\n{} 1\n", s.atoms.len()));
    }
    for a in &s.atoms {
        out.push_str(&format!("{:>4}", a.atomic_number));
        number(&mut out, a.x);
        number(&mut out, a.y);
        number(&mut out, a.z);
        out.push('\n');
    }
    out
}

fn write_cp2k(s: &Structure) -> String {
    let mut periodic: String = (0..3).filter(|&k| s.pbc[k]).map(|k| ['X', 'Y', 'Z'][k]).collect();
    if periodic.is_empty() {
        periodic = "NONE".into();
    }
    let mut out = String::from("&SUBSYS\n  &CELL\n");
    for k in 0..3 {
        out.push_str("    ");
        out.push(['A', 'B', 'C'][k]);
        for x in s.cell_vectors[k] {
            number(&mut out, x);
        }
        out.push('\n');
    }
    out.push_str(&format!("    PERIODIC {periodic}\n  &END CELL\n  &COORD\n"));
    for a in &s.atoms {
        out.push_str("    ");
        out.push_str(&left(&a.symbol, 3));
        number(&mut out, a.x);
        number(&mut out, a.y);
        number(&mut out, a.z);
        out.push('\n');
    }
    out.push_str("  &END COORD\n");
    for k in species(s) {
        out.push_str(&format!("  &KIND {k}\n    ELEMENT {k}\n    BASIS_SET DZVP-MOLOPT-SR-GTH\n    POTENTIAL GTH-PBE\n  &END KIND\n"));
    }
    out.push_str("&END SUBSYS\n");
    out
}

pub fn write_simulation_format(s: &Structure, format: &str) -> Result<String, String> {
    let f = canonical(format);
    if !s.has_unit_cell && (f == "pwi" || f == "cell" || f == "cp2k") {
        return Err("This format needs a unit cell; give the structure one with Edit > Structure Operations > Add Vacuum or Shrink-Wrap Cell.".into());
    }
    Ok(match f.as_str() {
        "pwi" => write_pwi(s),
        "cell" => write_cell(s),
        "gen" => write_gen(s),
        "xsf" => write_xsf(s),
        "cp2k" => write_cp2k(s),
        _ => return Err(format!("Unknown simulation format {format}")),
    })
}

// ---------------------------------------------------------------- readers

/// std::stod after replacing Fortran exponents (d, D); the whole word.
fn to_number(w: &str) -> Result<f64, String> {
    let t: String = w.chars().map(|c| if c == 'd' || c == 'D' { 'e' } else { c }).collect();
    let (v, used) = stod_used(&t)?;
    if used != t.len() {
        return Err(format!("'{w}' is not a number"));
    }
    Ok(v)
}

fn vector3(t: &[&str], from: usize) -> Result<[f64; 3], String> {
    if t.len() < from + 3 {
        return Err("Expected three numbers".into());
    }
    Ok([to_number(t[from])?, to_number(t[from + 1])?, to_number(t[from + 2])?])
}

fn line_at<'a>(lines: &[&'a str], i: usize) -> Result<&'a str, String> {
    lines.get(i).copied().ok_or_else(|| "The file ends early".to_string())
}

fn stol(word: Option<&&str>) -> Result<i64, String> {
    let w = word.ok_or("Missing number")?;
    let (v, used) = strtol(w);
    if used == 0 {
        return Err(format!("'{w}' is not a number"));
    }
    Ok(v)
}

fn set_cell(s: &mut Structure, rows: [[f64; 3]; 3]) {
    s.has_unit_cell = true;
    s.cell_vectors = rows;
}

fn read_xsf(text: &str, s: &mut Structure) -> Result<(), String> {
    let lines = lines(text, true);
    let mut rows = [[0.0; 3]; 3];
    let mut have_cell = false;
    let mut pbc = [true; 3];
    let mut i = 0;
    while i < lines.len() {
        let t = words(lines[i]);
        if t.is_empty() || t[0].starts_with('#') {
            i += 1;
            continue;
        }
        match t[0] {
            "SLAB" => pbc = [true, true, false],
            "POLYMER" => pbc = [true, false, false],
            "MOLECULE" => pbc = [false; 3],
            "PRIMVEC" => {
                for r in 0..3 {
                    rows[r] = vector3(&words(line_at(&lines, i + 1 + r)?), 0)?;
                }
                have_cell = true;
                i += 3;
            }
            "PRIMCOORD" => {
                let header = words(line_at(&lines, i + 1)?);
                let n = stol(header.first())?;
                for k in 0..n.max(0) as usize {
                    let a = words(line_at(&lines, i + 2 + k)?);
                    let label = a.first().ok_or("Missing atom")?;
                    s.atoms.push(atom(label, vector3(&a, 1)?)?);
                }
                i += 1 + n.max(0) as usize;
            }
            "ATOMS" => {
                let mut k = i + 1;
                while k < lines.len() {
                    let a = words(lines[k]);
                    if a.len() < 4 || (a[0].as_bytes()[0].is_ascii_alphabetic() && element_from(a[0]).is_empty()) {
                        break;
                    }
                    match vector3(&a, 1).and_then(|p| atom(a[0], p)) {
                        Ok(site) => s.atoms.push(site),
                        Err(_) => break,
                    }
                    i = k;
                    k += 1;
                }
            }
            _ => {}
        }
        i += 1;
    }
    if have_cell {
        set_cell(s, rows);
        s.pbc = pbc;
    }
    Ok(())
}

fn read_gen(text: &str, s: &mut Structure) -> Result<(), String> {
    let lines: Vec<&str> = lines(text, true)
        .into_iter()
        .filter(|line| {
            let t = words(line);
            !t.is_empty() && !t[0].starts_with('#')
        })
        .collect();
    let header = words(line_at(&lines, 0)?);
    let n = stol(header.first())?;
    let kind_char = header.get(1).ok_or("Missing geometry type")?.as_bytes()[0].to_ascii_uppercase();
    let kinds = words(line_at(&lines, 1)?);
    let mut raw = Vec::new();
    for k in 0..n.max(0) as usize {
        let a = words(line_at(&lines, 2 + k)?);
        let kind = stol(a.get(1))?;
        if kind < 1 || kind > kinds.len() as i64 {
            return Err("Species index out of range".into());
        }
        raw.push((kinds[(kind - 1) as usize].to_string(), vector3(&a, 2)?));
    }
    let n = n.max(0) as usize;
    if kind_char == b'S' || kind_char == b'F' {
        let origin = vector3(&words(line_at(&lines, 2 + n)?), 0)?;
        let mut rows = [[0.0; 3]; 3];
        for r in 0..3 {
            rows[r] = vector3(&words(line_at(&lines, 3 + n + r)?), 0)?;
        }
        set_cell(s, rows);
        s.cell_offset = origin;
    }
    for (symbol, position) in raw {
        let p = if kind_char == b'F' { from_fractional(s, position) } else { position };
        s.atoms.push(atom(&symbol, p)?);
    }
    Ok(())
}

fn block<'a>(lines: &[&'a str], name: &str) -> Option<Vec<Vec<&'a str>>> {
    for (i, line) in lines.iter().enumerate() {
        let t = words(line);
        if t.len() >= 2 && t[0].eq_ignore_ascii_case("%block") && t[1].to_ascii_lowercase() == name {
            let mut out = Vec::new();
            for row_line in &lines[i + 1..] {
                let row = words(row_line);
                if !row.is_empty() && row[0].eq_ignore_ascii_case("%endblock") {
                    return Some(out);
                }
                if !row.is_empty() && !row[0].starts_with('!') && !row[0].starts_with('#') {
                    out.push(row);
                }
            }
            return Some(out);
        }
    }
    None
}

fn unit_scale(rows: &mut Vec<Vec<&str>>) -> f64 {
    if !rows.is_empty() && rows[0].len() == 1 && rows[0][0].as_bytes()[0].is_ascii_alphabetic() {
        let unit = rows[0][0].to_ascii_lowercase();
        rows.remove(0);
        if unit == "bohr" || unit == "a0" {
            return BOHR;
        }
        if unit == "nm" {
            return 10.0;
        }
    }
    1.0
}

fn spin_of(row: &[&str]) -> Result<f64, String> {
    for k in 4..row.len() {
        let word = row[k].to_ascii_lowercase();
        if !word.starts_with("spin") {
            continue;
        }
        let mut value = word[4..].to_string();
        let mut next = k + 1;
        if value.starts_with('=') || value.starts_with(':') {
            value = value[1..].to_string();
        } else if next < row.len() && (row[next] == "=" || row[next] == ":") {
            next += 1;
        }
        if value.is_empty() && next < row.len() {
            value = row[next].to_string();
        }
        if !value.is_empty() {
            return to_number(&value);
        }
    }
    Ok(0.0)
}

fn cell_entry(rows: &[Vec<&str>], r: usize, c: usize) -> Result<f64, String> {
    to_number(rows.get(r).and_then(|row| row.get(c)).ok_or("The block ends early")?)
}

fn read_cell(text: &str, s: &mut Structure) -> Result<(), String> {
    let lines = lines(text, true);
    if let Some(mut lattice) = block(&lines, "lattice_cart") {
        let scale = unit_scale(&mut lattice);
        let mut rows = [[0.0; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                rows[r][c] = cell_entry(&lattice, r, c)? * scale;
            }
        }
        set_cell(s, rows);
    } else {
        let mut abc = block(&lines, "lattice_abc").ok_or("No LATTICE_CART or LATTICE_ABC block")?;
        let scale = unit_scale(&mut abc);
        let (a, b, c) = (cell_entry(&abc, 0, 0)? * scale, cell_entry(&abc, 0, 1)? * scale, cell_entry(&abc, 0, 2)? * scale);
        let rad = std::f64::consts::PI / 180.0;
        let (al, be, ga) = (cell_entry(&abc, 1, 0)? * rad, cell_entry(&abc, 1, 1)? * rad, cell_entry(&abc, 1, 2)? * rad);
        let cx = c * be.cos();
        let cy = c * (al.cos() - be.cos() * ga.cos()) / ga.sin();
        set_cell(s, [[a, 0.0, 0.0], [b * ga.cos(), b * ga.sin(), 0.0], [cx, cy, (c * c - cx * cx - cy * cy).max(0.0).sqrt()]]);
    }
    if let Some(mut positions) = block(&lines, "positions_abs") {
        let scale = unit_scale(&mut positions);
        for row in &positions {
            let x = vector3(row, 1)?;
            let mut site = atom(row.first().ok_or("Missing atom")?, [x[0] * scale, x[1] * scale, x[2] * scale])?;
            site.magmom = spin_of(row)?;
            s.atoms.push(site);
        }
    } else {
        let positions = block(&lines, "positions_frac").ok_or("No POSITIONS_ABS or POSITIONS_FRAC block")?;
        for row in &positions {
            let p = from_fractional(s, vector3(row, 1)?);
            let mut site = atom(row.first().ok_or("Missing atom")?, p)?;
            site.magmom = spin_of(row)?;
            s.atoms.push(site);
        }
    }
    let constraints = block(&lines, "ionic_constraints").unwrap_or_default();
    let mut by_species: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, a) in s.atoms.iter().enumerate() {
        by_species.entry(a.symbol.clone()).or_default().push(i);
    }
    for row in &constraints {
        if row.len() < 6 {
            continue;
        }
        let symbol = element_from(row[1]);
        let index = stol(row.get(2))?;
        let Some(list) = by_species.get(&symbol) else { continue };
        if index < 1 || index > list.len() as i64 {
            continue;
        }
        let d = vector3(row, 3)?;
        for k in 0..3 {
            if d[k].abs() > 0.5 && d[(k + 1) % 3].abs() < 1e-9 && d[(k + 2) % 3].abs() < 1e-9 {
                s.atoms[list[(index - 1) as usize]].fixed |= 1 << k;
            }
        }
    }
    Ok(())
}

/// The &SYSTEM assignments of a pw.x input (lower-case keys, values as written).
fn namelist_values(text: &str) -> std::collections::HashMap<String, String> {
    let mut out = std::collections::HashMap::new();
    for line in text.lines() {
        let line = line.split('!').next().unwrap_or("");
        for part in line.split(',') {
            let Some((key, value)) = part.split_once('=') else { continue };
            let key: String = key.trim().to_ascii_lowercase().chars().filter(|c| !c.is_whitespace()).collect();
            let value = value.trim().trim_end_matches('/').trim().trim_matches(|c| c == '\'' || c == '"').to_string();
            if !key.is_empty() && !value.is_empty() {
                out.insert(key, value);
            }
        }
    }
    out
}

/// The cell of a Bravais-lattice index (Quantum ESPRESSO's latgen), rows in
/// Angstrom, from celldm(1..6) (celldm(1) in bohr) or A, B, C (Angstrom)
/// with cosAB, cosAC, cosBC. Returns the cell and alat in Angstrom.
pub fn ibrav_cell(ibrav: i64, values: &std::collections::HashMap<String, String>) -> Result<([[f64; 3]; 3], f64), String> {
    let number = |key: &str| -> Result<Option<f64>, String> {
        match values.get(key) {
            None => Ok(None),
            Some(v) => v.replace(['d', 'D'], "e").parse::<f64>().map(Some).map_err(|_| format!("{key} = {v} is not a number")),
        }
    };
    let mut celldm = [0.0f64; 6];
    for (k, c) in celldm.iter_mut().enumerate() {
        *c = number(&format!("celldm({})", k + 1))?.unwrap_or(0.0);
    }
    if let Some(a) = number("a")? {
        // A, B, C in Angstrom and the cosines (cosBC for ibrav 14, cosAB for 12 and 13, cosAC for -12 and -13).
        let (b, c) = (number("b")?.unwrap_or(0.0), number("c")?.unwrap_or(0.0));
        let (cos_ab, cos_ac, cos_bc) = (number("cosab")?.unwrap_or(0.0), number("cosac")?.unwrap_or(0.0), number("cosbc")?.unwrap_or(0.0));
        celldm[0] = a / BOHR;
        celldm[1] = b / a;
        celldm[2] = c / a;
        match ibrav {
            14 => {
                celldm[3] = cos_bc;
                celldm[4] = cos_ac;
                celldm[5] = cos_ab;
            }
            -12 | -13 => celldm[4] = cos_ac,
            5 | -5 | 12 | 13 => celldm[3] = cos_ab,
            _ => {}
        }
    }
    if celldm[0] <= 0.0 {
        return Err(format!("ibrav = {ibrav} needs celldm(1) or A"));
    }
    let a = celldm[0] * BOHR;
    let (ba, ca) = (celldm[1], celldm[2]);
    let need = |ok: bool, what: &str| if ok { Ok(()) } else { Err(format!("ibrav = {ibrav} needs {what}")) };
    let v = |x: f64, y: f64, z: f64| [a * x, a * y, a * z];
    let rows = match ibrav {
        1 => [v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)],
        2 => [v(-0.5, 0.0, 0.5), v(0.0, 0.5, 0.5), v(-0.5, 0.5, 0.0)],
        3 => [v(0.5, 0.5, 0.5), v(-0.5, 0.5, 0.5), v(-0.5, -0.5, 0.5)],
        -3 => [v(-0.5, 0.5, 0.5), v(0.5, -0.5, 0.5), v(0.5, 0.5, -0.5)],
        4 => {
            need(ca > 0.0, "celldm(3) or C")?;
            [v(1.0, 0.0, 0.0), v(-0.5, 3f64.sqrt() / 2.0, 0.0), v(0.0, 0.0, ca)]
        }
        5 | -5 => {
            let c = celldm[3];
            need(c > -0.5 && c < 1.0, "celldm(4) (cos gamma) between -1/2 and 1")?;
            let tx = ((1.0 - c) / 2.0).sqrt();
            let ty = ((1.0 - c) / 6.0).sqrt();
            let tz = ((1.0 + 2.0 * c) / 3.0).sqrt();
            if ibrav == 5 {
                [v(tx, -ty, tz), v(0.0, 2.0 * ty, tz), v(-tx, -ty, tz)]
            } else {
                let u = tz - 2.0 * 2f64.sqrt() * ty;
                let w = tz + 2f64.sqrt() * ty;
                let s3 = 1.0 / 3f64.sqrt();
                [v(u * s3, w * s3, w * s3), v(w * s3, u * s3, w * s3), v(w * s3, w * s3, u * s3)]
            }
        }
        6 => {
            need(ca > 0.0, "celldm(3) or C")?;
            [v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, ca)]
        }
        7 => {
            need(ca > 0.0, "celldm(3) or C")?;
            [v(0.5, -0.5, ca / 2.0), v(0.5, 0.5, ca / 2.0), v(-0.5, -0.5, ca / 2.0)]
        }
        8 | 9 | -9 | 91 | 10 | 11 => {
            need(ba > 0.0 && ca > 0.0, "celldm(2) and celldm(3), or B and C")?;
            match ibrav {
                8 => [v(1.0, 0.0, 0.0), v(0.0, ba, 0.0), v(0.0, 0.0, ca)],
                9 => [v(0.5, ba / 2.0, 0.0), v(-0.5, ba / 2.0, 0.0), v(0.0, 0.0, ca)],
                -9 => [v(0.5, -ba / 2.0, 0.0), v(0.5, ba / 2.0, 0.0), v(0.0, 0.0, ca)],
                91 => [v(1.0, 0.0, 0.0), v(0.0, ba / 2.0, -ca / 2.0), v(0.0, ba / 2.0, ca / 2.0)],
                10 => [v(0.5, 0.0, ca / 2.0), v(0.5, ba / 2.0, 0.0), v(0.0, ba / 2.0, ca / 2.0)],
                _ => [v(0.5, ba / 2.0, ca / 2.0), v(-0.5, ba / 2.0, ca / 2.0), v(-0.5, -ba / 2.0, ca / 2.0)],
            }
        }
        12 | 13 => {
            need(ba > 0.0 && ca > 0.0, "celldm(2) and celldm(3), or B and C")?;
            let cg = celldm[3];
            need(cg.abs() < 1.0, "celldm(4) (cos gamma) or cosAB")?;
            let sg = (1.0 - cg * cg).sqrt();
            if ibrav == 12 {
                [v(1.0, 0.0, 0.0), v(ba * cg, ba * sg, 0.0), v(0.0, 0.0, ca)]
            } else {
                [v(0.5, 0.0, -ca / 2.0), v(ba * cg, ba * sg, 0.0), v(0.5, 0.0, ca / 2.0)]
            }
        }
        -12 | -13 => {
            need(ba > 0.0 && ca > 0.0, "celldm(2) and celldm(3), or B and C")?;
            let cb = celldm[4];
            need(cb.abs() < 1.0, "celldm(5) (cos beta) or cosAC")?;
            let sb = (1.0 - cb * cb).sqrt();
            if ibrav == -12 {
                [v(1.0, 0.0, 0.0), v(0.0, ba, 0.0), v(ca * cb, 0.0, ca * sb)]
            } else {
                [v(0.5, ba / 2.0, 0.0), v(-0.5, ba / 2.0, 0.0), v(ca * cb, 0.0, ca * sb)]
            }
        }
        14 => {
            need(ba > 0.0 && ca > 0.0, "celldm(2) and celldm(3), or B and C")?;
            let (cos_a, cos_b, cos_g) = (celldm[3], celldm[4], celldm[5]);
            let sg = (1.0 - cos_g * cos_g).sqrt();
            let volume = 1.0 + 2.0 * cos_a * cos_b * cos_g - cos_a * cos_a - cos_b * cos_b - cos_g * cos_g;
            need(sg > 0.0 && volume > 0.0, "cosines that make a cell")?;
            [v(1.0, 0.0, 0.0), v(ba * cos_g, ba * sg, 0.0), v(ca * cos_b, ca * (cos_a - cos_b * cos_g) / sg, ca * volume.sqrt() / sg)]
        }
        _ => return Err(format!("ibrav = {ibrav} is not a Quantum ESPRESSO lattice")),
    };
    Ok((rows, a))
}

fn read_pwi(text: &str, s: &mut Structure) -> Result<(), String> {
    let namelist = namelist_values(text);
    let ibrav = namelist.get("ibrav").map(|v| strtol(v).0).unwrap_or(0);
    let lattice = if ibrav != 0 { Some(ibrav_cell(ibrav, &namelist)?) } else { None };
    let lines = lines(text, true);
    let mut alat = 0.0;
    let mut nat: i64 = -1;
    for line in &lines {
        let mut l: String = line.to_ascii_lowercase().chars().filter(|&c| c != ' ').collect();
        if let Some(bang) = l.find('!') {
            l.truncate(bang);
        }
        let value = |key: &str| l.find(&format!("{key}=")).map(|at| l[at + key.len() + 1..].to_string());
        if let Some(v) = value("celldm(1)") {
            let end = v.find([',', '/']).unwrap_or(v.len());
            alat = to_number(&v[..end])? * BOHR;
        }
        if let Some(v) = value("nat") {
            nat = strtol(&v).0;
        }
        if l.starts_with("a=") || l.contains(",a=") {
            let at = if l.starts_with("a=") { 0 } else { l.find(",a=").unwrap() + 1 };
            let rest = &l[at + 2..];
            let end = rest.find([',', '/']).unwrap_or(rest.len());
            alat = to_number(&rest[..end])?;
        }
    }
    let strip = |u: &str| -> String { u.to_ascii_lowercase().chars().filter(|c| !"(){}".contains(*c)).collect() };
    let mut rows = [[0.0; 3]; 3];
    let mut have_cell = false;
    for (i, line) in lines.iter().enumerate() {
        let t = words(line);
        if t.is_empty() {
            continue;
        }
        if t[0].to_ascii_lowercase() == "cell_parameters" {
            let unit = t.get(1).map(|u| strip(u)).unwrap_or_default();
            let scale = if unit == "bohr" {
                BOHR
            } else if unit == "angstrom" {
                1.0
            } else if alat > 0.0 {
                alat
            } else {
                1.0
            };
            for r in 0..3 {
                rows[r] = vector3(&words(line_at(&lines, i + 1 + r)?), 0)?;
                for x in rows[r].iter_mut() {
                    *x *= scale;
                }
            }
            have_cell = true;
        }
    }
    if let Some((cell, a)) = lattice {
        // ibrav != 0: the lattice comes from the namelist, and alat is its a.
        rows = cell;
        alat = a;
        have_cell = true;
    }
    if !have_cell {
        return Err("No CELL_PARAMETERS card".into());
    }
    set_cell(s, rows);
    for (i, line) in lines.iter().enumerate() {
        let t = words(line);
        if t.is_empty() || t[0].to_ascii_lowercase() != "atomic_positions" {
            continue;
        }
        let unit = t.get(1).map(|u| strip(u)).unwrap_or_else(|| "alat".into());
        let mut k = i + 1;
        while k < lines.len() && (nat < 0 || (s.atoms.len() as i64) < nat) {
            let a = words(lines[k]);
            if a.len() < 4 || element_from(a[0]).is_empty() {
                break;
            }
            let mut x = vector3(&a, 1)?;
            if unit == "crystal" {
                x = from_fractional(s, x);
            } else {
                let scale = if unit == "bohr" {
                    BOHR
                } else if unit == "angstrom" {
                    1.0
                } else if alat > 0.0 {
                    alat
                } else {
                    1.0
                };
                for v in x.iter_mut() {
                    *v *= scale;
                }
            }
            let mut site = atom(a[0], x)?;
            if a.len() >= 7 {
                for c in 0..3 {
                    if a[4 + c] == "0" {
                        site.fixed |= 1 << c;
                    }
                }
            }
            s.atoms.push(site);
            k += 1;
        }
        break;
    }
    if s.atoms.is_empty() {
        return Err("No ATOMIC_POSITIONS card".into());
    }
    Ok(())
}

pub fn read_simulation_format(text: &str, format: &str) -> Result<Structure, String> {
    let mut s = Structure::default();
    let result = match canonical(format).as_str() {
        "xsf" => read_xsf(text, &mut s),
        "gen" => read_gen(text, &mut s),
        "cell" => read_cell(text, &mut s),
        "pwi" => read_pwi(text, &mut s),
        _ => Err(format!("Unknown simulation format {format}")),
    }
    .and_then(|_| if s.atoms.is_empty() { Err("No atoms found".to_string()) } else { Ok(()) });
    match result {
        Ok(()) => Ok(s),
        Err(e) => Err(format!("Could not read the {format} file: {e}")),
    }
}

#[cfg(test)]
mod ibrav_tests {
    use super::*;

    fn cell(text: &str) -> [[f64; 3]; 3] {
        let values = namelist_values(text);
        let ibrav = strtol(values.get("ibrav").unwrap()).0;
        ibrav_cell(ibrav, &values).unwrap().0
    }

    fn volume(r: [[f64; 3]; 3]) -> f64 {
        let (a, b, c) = (glam::DVec3::from(r[0]), glam::DVec3::from(r[1]), glam::DVec3::from(r[2]));
        a.dot(b.cross(c)).abs()
    }

    #[test]
    fn lattices() {
        // fcc: a quarter of the cube; A in Angstrom equals celldm(1) in bohr.
        let fcc = cell("&system\n ibrav = 2, celldm(1) = 6.82, nat = 1\n/");
        let a = 6.82 * BOHR;
        assert!((volume(fcc) - a * a * a / 4.0).abs() < 1e-9);
        let same = cell(&format!("&system\n ibrav=2, A={a}\n/"));
        for r in 0..3 {
            for k in 0..3 {
                assert!((fcc[r][k] - same[r][k]).abs() < 1e-12);
            }
        }
        // bcc half, hexagonal, rhombohedral (both settings) and triclinic volumes.
        assert!((volume(cell("ibrav=3, A=2.87")) - 2.87f64.powi(3) / 2.0).abs() < 1e-9);
        assert!((volume(cell("ibrav=4, A=3.21, C=5.21")) - 3.21 * 3.21 * 5.21 * 3f64.sqrt() / 2.0).abs() < 1e-9);
        let c: f64 = 0.3;
        let rhombohedral = 1.0 - 3.0 * c * c + 2.0 * c * c * c;
        assert!((volume(cell("ibrav=5, A=4, cosAB=0.3")) - 64.0 * rhombohedral.sqrt()).abs() < 1e-9);
        assert!((volume(cell("ibrav=-5, A=4, cosAB=0.3")) - 64.0 * rhombohedral.sqrt()).abs() < 1e-9);
        let t = cell("ibrav=14, A=3, B=4, C=5, cosBC=0.1, cosAC=0.2, cosAB=0.3");
        let g = 1.0 + 2.0 * 0.1 * 0.2 * 0.3 - 0.01 - 0.04 - 0.09;
        assert!((volume(t) - 60.0 * f64::sqrt(g)).abs() < 1e-9);
        assert!((glam::DVec3::from(t[2]).length() - 5.0).abs() < 1e-12);
        assert!(ibrav_cell(4, &namelist_values("ibrav=4, A=3")).is_err());
    }

    /// Every lattice: the primitive volume is the conventional one over the
    /// lattice points per cell, and the oblique angles come out as given.
    #[test]
    fn every_lattice_has_its_volume() {
        let (a, b, c): (f64, f64, f64) = (3.0, 4.0, 5.0);
        let (sg, sb) = ((1.0f64 - 0.04).sqrt(), (1.0f64 - 0.04).sqrt());
        let cases: [(i64, &str, f64); 18] = [
            (1, "A=3", a * a * a),
            (2, "A=3", a * a * a / 4.0),
            (3, "A=3", a * a * a / 2.0),
            (-3, "A=3", a * a * a / 2.0),
            (6, "A=3, C=5", a * a * c),
            (7, "A=3, C=5", a * a * c / 2.0),
            (8, "A=3, B=4, C=5", a * b * c),
            (9, "A=3, B=4, C=5", a * b * c / 2.0),
            (-9, "A=3, B=4, C=5", a * b * c / 2.0),
            (91, "A=3, B=4, C=5", a * b * c / 2.0),
            (10, "A=3, B=4, C=5", a * b * c / 4.0),
            (11, "A=3, B=4, C=5", a * b * c / 2.0),
            (12, "A=3, B=4, C=5, cosAB=0.2", a * b * c * sg),
            (13, "A=3, B=4, C=5, cosAB=0.2", a * b * c * sg / 2.0),
            (-12, "A=3, B=4, C=5, cosAC=0.2", a * b * c * sb),
            (-13, "A=3, B=4, C=5, cosAC=0.2", a * b * c * sb / 2.0),
            (4, "A=3, C=5", a * a * c * 3f64.sqrt() / 2.0),
            (14, "A=3, B=4, C=5, cosBC=0, cosAC=0, cosAB=0.2", a * b * c * sg),
        ];
        for (ibrav, text, expected) in cases {
            let r = cell(&format!("ibrav={ibrav}, {text}"));
            assert!((volume(r) - expected).abs() < 1e-9, "ibrav {ibrav}: {} != {expected}", volume(r));
        }
        // Oblique cells keep their lengths and angles.
        let m = cell("ibrav=12, A=3, B=4, C=5, cosAB=0.2");
        let (va, vb) = (glam::DVec3::from(m[0]), glam::DVec3::from(m[1]));
        assert!((va.length() - 3.0).abs() < 1e-12 && (vb.length() - 4.0).abs() < 1e-12);
        assert!((va.dot(vb) / 12.0 - 0.2).abs() < 1e-12);
        let t = cell("ibrav=14, A=3, B=4, C=5, cosBC=0.1, cosAC=0.2, cosAB=0.3");
        let (ta, tb, tc) = (glam::DVec3::from(t[0]), glam::DVec3::from(t[1]), glam::DVec3::from(t[2]));
        assert!((tb.dot(tc) / 20.0 - 0.1).abs() < 1e-12 && (ta.dot(tc) / 15.0 - 0.2).abs() < 1e-12 && (ta.dot(tb) / 12.0 - 0.3).abs() < 1e-12);
        // Rhombohedral: three equal vectors at the given angle.
        for ibrav in [5, -5] {
            let h = cell(&format!("ibrav={ibrav}, A=4, cosAB=0.3"));
            for i in 0..3 {
                let (p, q) = (glam::DVec3::from(h[i]), glam::DVec3::from(h[(i + 1) % 3]));
                assert!((p.length() - 4.0).abs() < 1e-12 && (p.dot(q) / 16.0 - 0.3).abs() < 1e-12, "ibrav {ibrav}");
            }
        }
    }

    #[test]
    fn reads_ibrav_input() {
        let text = "&control\n/\n&system\n  ibrav = 2, celldm(1) = 6.82, nat = 2, ntyp = 1\n/\nATOMIC_SPECIES\nCu 63.5 Cu.upf\nATOMIC_POSITIONS crystal\nCu 0 0 0\nCu 0.5 0.5 0.5\n";
        let s = read_simulation_format(text, "pwi").unwrap();
        assert_eq!(s.atoms.len(), 2);
        assert!(s.has_unit_cell);
    }
}
