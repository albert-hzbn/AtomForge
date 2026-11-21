//! VASP POSCAR/CONTCAR (the C++ readPoscar): scale (negative = volume, three
//! = per axis), lattice, VASP 5 element line or VASP 4 title species, counts,
//! selective dynamics and Direct or Cartesian coordinates.

use af_core::{elements, AtomSite, Structure};
use glam::{DMat3, DVec3};

use crate::text::{lines, number_full, words};

fn symbol_number(word: &str) -> i32 {
    let letters = word.chars().take_while(|c| c.is_ascii_alphabetic()).count();
    if letters == 0 {
        return 0;
    }
    let mut symbol: String = word.chars().take(letters.min(2)).collect();
    let first = symbol.remove(0).to_ascii_uppercase();
    let rest = symbol.to_ascii_lowercase();
    elements::atomic_number(&format!("{first}{rest}"))
}

pub fn read(text: &str, structure: &mut Structure) -> Result<(), String> {
    let lines = lines(text, true);
    let fail = |message: String| Err(format!("Invalid POSCAR: {message}"));
    if lines.len() < 7 {
        return fail("the file is incomplete.".into());
    }
    let scale_words = words(lines[1]);
    let mut scale = [1.0f64; 3];
    let mut volume = 0.0;
    let per_axis = scale_words.len() >= 3 && scale_words[..3].iter().all(|w| number_full(w).is_some());
    if !scale_words.is_empty() && !per_axis {
        match number_full(scale_words[0]) {
            Some(v) if v != 0.0 => scale[0] = v,
            _ => return fail("the scale factor on line 2 is not a nonzero number.".into()),
        }
        if scale[0] < 0.0 {
            volume = -scale[0];
            scale[0] = 1.0;
        }
        scale[1] = scale[0];
        scale[2] = scale[0];
    } else if per_axis {
        for k in 0..3 {
            match number_full(scale_words[k]) {
                Some(v) if v > 0.0 => scale[k] = v,
                _ => return fail("the scale factors on line 2 must be positive.".into()),
            }
        }
    } else {
        return fail("line 2 needs the scale factor.".into());
    }
    let mut columns = [DVec3::ZERO; 3];
    for r in 0..3 {
        let row = words(lines[2 + r]);
        let values: Option<Vec<f64>> = if row.len() >= 3 { row[..3].iter().map(|w| number_full(w)).collect() } else { None };
        let Some(v) = values else { return fail(format!("lattice vector {} needs three numbers.", r + 1)) };
        columns[r] = DVec3::new(v[0] * scale[0], v[1] * scale[1], v[2] * scale[2]);
    }
    let mut cell = DMat3::from_cols(columns[0], columns[1], columns[2]);
    let determinant = cell.determinant();
    if determinant.abs() < 1e-12 {
        return fail("the lattice vectors have no volume.".into());
    }
    if volume > 0.0 {
        cell *= (volume / determinant.abs()).cbrt();
    }

    let mut line = 5usize;
    let species_words = words(lines[line]);
    let mut elements_list: Vec<i32> = Vec::new();
    let has_element_line = !species_words.is_empty() && !species_words[0].as_bytes()[0].is_ascii_digit();
    if has_element_line {
        for word in &species_words {
            let z = symbol_number(word);
            if z == 0 {
                return fail(format!("unknown element '{word}' on line 6."));
            }
            elements_list.push(z);
        }
        line += 1;
    }
    if line >= lines.len() {
        return fail("the atom counts are missing.".into());
    }
    let mut counts: Vec<i64> = Vec::new();
    for word in words(lines[line]) {
        let Some(value) = number_full(word) else { break };
        if value < 0.0 || value != value.floor() {
            return fail(format!("the atom counts on line {} must be whole numbers.", line + 1));
        }
        counts.push(value as i64);
    }
    if counts.is_empty() {
        return fail(format!("the atom counts on line {} are missing.", line + 1));
    }
    line += 1;
    if !has_element_line {
        let title = words(lines[0]);
        for word in title.iter().take(counts.len()) {
            let z = symbol_number(word);
            if z == 0 {
                elements_list.clear();
                break;
            }
            elements_list.push(z);
        }
        if elements_list.len() != counts.len() {
            return fail("it names no elements (VASP 4): put them on the title line, e.g. \"Na Cl\", or add an element line.".into());
        }
    }
    if elements_list.len() != counts.len() {
        return fail("the element line and the counts do not match.".into());
    }
    let total: i64 = counts.iter().sum();
    if total <= 0 {
        return fail("it has no atoms.".into());
    }
    if total > 5_000_000 {
        return fail("it has more than 5000000 atoms.".into());
    }
    let mut selective = false;
    if line < lines.len() {
        let w = words(lines[line]);
        if !w.is_empty() && (w[0].starts_with('S') || w[0].starts_with('s')) {
            selective = true;
            line += 1;
        }
    }
    if line >= lines.len() {
        return fail("the coordinate type (Direct or Cartesian) is missing.".into());
    }
    let mode = words(lines[line]);
    if mode.is_empty() {
        return fail("the coordinate type (Direct or Cartesian) is missing.".into());
    }
    let kind = mode[0].as_bytes()[0].to_ascii_lowercase();
    let cartesian = kind == b'c' || kind == b'k';
    line += 1;
    let cartesian_scale = if volume > 0.0 { (volume / determinant.abs()).cbrt() } else { 1.0 };
    let mut fixed_count = 0usize;
    structure.atoms.reserve(total as usize);
    for (species, &z) in elements_list.iter().enumerate() {
        for _ in 0..counts[species] {
            if line >= lines.len() {
                return fail(format!("it declares {total} atoms but lists only {}.", structure.atoms.len()));
            }
            let w = words(lines[line]);
            let values: Option<Vec<f64>> = if w.len() >= 3 { w[..3].iter().map(|x| number_full(x)).collect() } else { None };
            let Some(v) = values else { return fail(format!("atom {} needs three coordinates.", structure.atoms.len() + 1)) };
            let position = if cartesian {
                DVec3::new(v[0] * scale[0], v[1] * scale[1], v[2] * scale[2]) * cartesian_scale
            } else {
                cell * DVec3::new(v[0], v[1], v[2])
            };
            let mut site = AtomSite::new(elements::symbol(z), position.to_array());
            site.atomic_number = z;
            if selective && w.len() >= 6 {
                for k in 0..3 {
                    let flag = w[3 + k].as_bytes()[0];
                    if flag == b'F' || flag == b'f' {
                        site.fixed |= 1 << k;
                    }
                }
            }
            if site.fixed != 0 {
                fixed_count += 1;
            }
            structure.atoms.push(site);
            line += 1;
        }
    }
    structure.has_unit_cell = true;
    for r in 0..3 {
        structure.cell_vectors[r] = cell.col(r).to_array();
    }
    if selective {
        structure.load_note = format!(
            "{fixed_count}{} fixed coordinates (selective dynamics).",
            if fixed_count == 1 { " atom has" } else { " atoms have" }
        );
    }
    Ok(())
}
