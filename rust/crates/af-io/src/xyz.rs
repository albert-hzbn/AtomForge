//! XYZ and extended XYZ (the first frame): count, comment with Lattice,
//! Properties and pbc keys, one line per atom (the C++ readXyz).

use af_core::{elements, AtomSite, Structure};

use crate::text::{is_space, strtol};

/// Lattice="ax ay az bx by bz cx cy cz" (the key in either case of 'L').
pub fn extended_lattice(comment: &str) -> Option<[[f64; 3]; 3]> {
    let bytes = comment.as_bytes();
    let mut key = None;
    let mut from = 0;
    while let Some(found) = comment[from..].find("attice=") {
        let at = from + found;
        if at > 0 && (bytes[at - 1] == b'L' || bytes[at - 1] == b'l') {
            key = Some(at + 7);
            break;
        }
        from = at + 1;
    }
    let key = key?;
    if key >= comment.len() {
        return None;
    }
    let quote = bytes[key];
    if quote != b'"' && quote != b'\'' {
        return None;
    }
    let close = key + 1 + comment[key + 1..].find(quote as char)?;
    let inner = &comment[key + 1..close];
    // istringstream >> double: each word must start with a number.
    let mut cell = [[0.0; 3]; 3];
    let mut rest = inner;
    for row in cell.iter_mut() {
        for value in row.iter_mut() {
            let (v, used) = af_core::format::parse_prefix_f64(rest)?;
            if !v.is_finite() {
                return None;
            }
            *value = v;
            rest = &rest[used..];
        }
    }
    Some(cell)
}

#[derive(Default)]
struct Extra {
    magmom: i32,
    charge: i32,
    velocity: i32,
}

/// Columns of species and positions from Properties=name:type:count:...; 1
/// when found, 0 when the key names neither, -1 without the key.
fn extended_columns(comment: &str, species: &mut i32, position: &mut i32, extra: &mut Extra) -> i32 {
    let bytes = comment.as_bytes();
    let mut key = None;
    let mut from = 0;
    while let Some(found) = comment[from..].find("roperties=") {
        let at = from + found;
        if at > 0 && (bytes[at - 1] == b'P' || bytes[at - 1] == b'p') {
            key = Some(at + 10);
            break;
        }
        from = at + 1;
    }
    let Some(key) = key else { return -1 };
    let end = comment[key..].find([' ', '\t', '\r', '\n']).map_or(comment.len(), |e| key + e);
    let spec = &comment[key..end];
    *species = -1;
    *position = -1;
    let mut number_column = -1;
    let mut column = 0;
    // getline(spec, piece, ':') three times per field; std::getline gives no
    // empty piece after a final ':'.
    let mut parts: Vec<&str> = spec.split(':').collect();
    if spec.ends_with(':') {
        parts.pop();
    }
    let mut k = 0;
    while k + 2 < parts.len() {
        let (name, _type, count) = (parts[k], parts[k + 1], parts[k + 2]);
        let (n, used) = strtol(count);
        if used == 0 {
            return 0;
        }
        if n < 1 {
            return 0;
        }
        let n = n as i32;
        let lower = name.to_ascii_lowercase();
        if lower == "species" && n == 1 {
            *species = column;
        }
        if name == "Z" && n == 1 {
            number_column = column;
        }
        if (lower == "pos" || lower == "positions") && n == 3 {
            *position = column;
        }
        if (lower == "initial_magmoms" || lower == "magmoms" || lower == "magmom") && n == 1 {
            extra.magmom = column;
        }
        if (lower == "initial_charges" || lower == "charges" || lower == "charge") && n == 1 {
            extra.charge = column;
        }
        if (lower == "velo" || lower == "vel" || lower == "velocities") && n == 3 {
            extra.velocity = column;
        }
        column += n;
        k += 3;
    }
    if *species < 0 {
        *species = number_column;
    }
    if *species >= 0 && *position >= 0 {
        1
    } else {
        0
    }
}

fn element_of_label(label: &str) -> i32 {
    if label.as_bytes()[0].is_ascii_digit() {
        return crate::text::atoi(label);
    }
    let mut symbol = String::new();
    for c in label.chars() {
        if !c.is_ascii_alphabetic() || symbol.len() == 2 {
            break;
        }
        symbol.push(if symbol.is_empty() { c.to_ascii_uppercase() } else { c.to_ascii_lowercase() });
    }
    let mut z = elements::atomic_number(&symbol);
    // "CA1" style labels: a one-letter element followed by more letters.
    if z == 0 && symbol.len() == 2 {
        z = elements::atomic_number(&symbol[..1]);
    }
    z
}

fn real(word: &str) -> Option<f64> {
    crate::text::number_full(word)
}

pub fn read(text: &str, structure: &mut Structure) -> Result<(), String> {
    let mut pos = 0usize;
    let mut next_line = || -> Option<&str> {
        if pos >= text.len() {
            return None;
        }
        let end = text[pos..].find('\n').map_or(text.len(), |e| pos + e);
        let mut line = &text[pos..end];
        if let Some(stripped) = line.strip_suffix('\r') {
            line = stripped;
        }
        pos = end + 1;
        Some(line)
    };
    let mut line;
    loop {
        line = next_line().ok_or_else(|| "The XYZ file is empty.".to_string())?;
        if line.chars().any(|c| c != ' ' && c != '\t') {
            break;
        }
    }
    let (count, used) = strtol(line);
    if used == 0 || count < 0 || line[used..].chars().any(|c| c != ' ' && c != '\t') {
        return Err("Not an XYZ file: the first line must be the number of atoms.".into());
    }
    if count > 5_000_000 {
        return Err(format!(
            "The file has {count} atoms, more than the 5000000 that can be edited; open it to view it as an atom cloud."
        ));
    }
    let comment = next_line().unwrap_or("").to_string();
    let (mut species, mut position) = (0, 1);
    let mut extra = Extra { magmom: -1, charge: -1, velocity: -1 };
    let columns = extended_columns(&comment, &mut species, &mut position, &mut extra);
    if columns == 0 {
        return Err("The extended XYZ Properties key names no species (or Z) and pos columns.".into());
    }
    if columns < 0 {
        species = 0;
        position = 1;
        extra = Extra { magmom: -1, charge: -1, velocity: -1 };
    }
    let needed = [species, position + 2, extra.magmom, extra.charge, if extra.velocity >= 0 { extra.velocity + 2 } else { -1 }]
        .into_iter()
        .max()
        .unwrap()
        + 1;
    structure.atoms.reserve(count as usize);
    for n in 0..count {
        let Some(line) = next_line() else {
            return Err(format!("The XYZ file declares {count} atoms but has only {n}."));
        };
        let words: Vec<&str> =
            line.split(|c| c == ' ' || c == '\t').filter(|w| !w.is_empty()).take(needed as usize).collect();
        if (words.len() as i32) < needed {
            return Err(format!("XYZ atom {} needs an element and three coordinates: {line}", n + 1));
        }
        let label = words[species as usize];
        let z = element_of_label(label);
        if !(1..=118).contains(&z) {
            return Err(format!("Unknown element '{label}' on XYZ atom {}.", n + 1));
        }
        let mut xyz = [0.0; 3];
        for k in 0..3 {
            let word = words[(position + k as i32) as usize];
            xyz[k] = real(word).ok_or_else(|| format!("XYZ atom {} has an invalid coordinate '{word}'.", n + 1))?;
        }
        let mut site = AtomSite::new(elements::symbol(z), xyz);
        site.atomic_number = z;
        let mut ok = true;
        if extra.magmom >= 0 {
            match real(words[extra.magmom as usize]) {
                Some(v) => site.magmom = v,
                None => ok = false,
            }
        }
        if extra.charge >= 0 {
            match real(words[extra.charge as usize]) {
                Some(v) => site.charge = v,
                None => ok = false,
            }
        }
        if extra.velocity >= 0 {
            for k in 0..3 {
                match real(words[(extra.velocity + k) as usize]) {
                    Some(v) => site.velocity[k as usize] = v,
                    None => ok = false,
                }
            }
        }
        if !ok {
            return Err(format!("XYZ atom {} has an invalid magnetic moment, charge or velocity.", n + 1));
        }
        structure.atoms.push(site);
    }
    if let Some(cell) = extended_lattice(&comment) {
        let mut periodic = true;
        let mut axes = [true; 3];
        let bytes = comment.as_bytes();
        let mut from = 0;
        while let Some(found) = comment[from..].find("pbc=") {
            let at = from + found;
            from = at + 1;
            if at > 0 && !is_space(bytes[at - 1] as char) {
                continue;
            }
            let start = at + 4;
            let quote = bytes.get(start).copied().unwrap_or(b' ');
            let value = if quote == b'"' || quote == b'\'' {
                let close = comment[start + 1..].find(quote as char).map_or(comment.len(), |c| start + 1 + c);
                &comment[start + 1..close]
            } else {
                let stop = comment[start..].find([' ', '\t']).map_or(comment.len(), |c| start + c);
                &comment[start..stop]
            };
            let flag = |w: &str| w.contains(['T', 't', 'Y', 'y', '1']);
            periodic = flag(value);
            for (axis, word) in crate::text::words(value).into_iter().take(3).enumerate() {
                axes[axis] = flag(word);
            }
            break;
        }
        if periodic {
            structure.has_unit_cell = true;
            structure.cell_vectors = cell;
            structure.pbc = axes;
        } else {
            structure.load_note = "The file's cell is not periodic (pbc=\"F F F\"), so the structure opens without a cell.".into();
        }
    }
    Ok(())
}
