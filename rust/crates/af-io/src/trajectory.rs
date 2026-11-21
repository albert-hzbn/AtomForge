//! Trajectories (the C++ io/Trajectory): XYZ/extXYZ frames and LAMMPS dump
//! frames (element or type columns, Cartesian or scaled positions,
//! orthogonal or triclinic boxes, atoms sorted by id).

use std::io::BufRead;
use std::sync::Mutex;

use af_core::format::parse_prefix_f64;
use af_core::task::checkpoint;
use af_core::{elements, AtomSite, Structure};

fn standard_symbol(symbol: &str) -> String {
    if symbol.is_empty() || !symbol.chars().all(|c| c.is_ascii_alphabetic()) {
        return symbol.to_string();
    }
    let mut chars = symbol.chars();
    let first = chars.next().unwrap().to_ascii_uppercase();
    std::iter::once(first).chain(chars.map(|c| c.to_ascii_lowercase())).collect()
}

fn blank(line: &str) -> bool {
    line.chars().all(|c| c == ' ' || c == '\t' || c == '\r')
}

/// `istream >> double` on one word.
fn read_double(word: &str) -> Option<f64> {
    parse_prefix_f64(word).map(|v| v.0)
}

fn read_lines(path: &str) -> Result<Vec<String>, String> {
    let file = std::fs::File::open(path).map_err(|_| "Cannot open trajectory".to_string())?;
    std::io::BufReader::new(file).lines().collect::<Result<_, _>>().map_err(|_| "Cannot open trajectory".to_string())
}

pub fn load_xyz_trajectory(path: &str) -> Result<Vec<Structure>, String> {
    let lines = read_lines(path)?;
    let mut frames = Vec::new();
    let mut total = 0usize;
    let mut i = 0;
    while i < lines.len() {
        checkpoint()?;
        let line = &lines[i];
        i += 1;
        if blank(line) {
            continue;
        }
        let w: Vec<&str> = line.split_whitespace().collect();
        let count = w.first().and_then(|c| crate::text::stoll(c).ok().filter(|_| crate::text::strtol(c).1 == c.len()));
        let count = match (count, w.len()) {
            (Some(n), 1) if (1..=1_000_000).contains(&n) => n as usize,
            _ => return Err("Invalid XYZ frame atom count".into()),
        };
        total += count;
        if total > 10_000_000 || frames.len() >= 100_000 {
            return Err("Trajectory exceeds memory limit".into());
        }
        let comment = lines.get(i).ok_or("Missing XYZ frame comment")?.clone();
        i += 1;
        let mut frame = Structure::default();
        if let Some(at) = comment.find("Lattice=\"") {
            let mut rest = &comment[at + 9..];
            for r in 0..3 {
                for c in 0..3 {
                    let (v, used) = parse_prefix_f64(rest).ok_or("Invalid XYZ lattice")?;
                    if !v.is_finite() {
                        return Err("Invalid XYZ lattice".into());
                    }
                    frame.cell_vectors[r][c] = v;
                    rest = &rest[used..];
                }
            }
            frame.has_unit_cell = true;
        }
        let (mut species_column, mut position_column, mut column_count) = (0usize, 1usize, 4usize);
        if let Some(at) = comment.find("Properties=") {
            let spec = &comment[at + 11..];
            let spec = &spec[..spec.find([' ', '\t', '\r']).unwrap_or(spec.len())];
            let mut parts: Vec<&str> = spec.split(':').collect();
            if spec.ends_with(':') || spec.is_empty() {
                parts.pop();
            }
            if parts.len() % 3 != 0 {
                return Err("Invalid extXYZ Properties".into());
            }
            let (mut has_species, mut has_positions) = (false, false);
            column_count = 0;
            for p in (0..parts.len()).step_by(3) {
                let width = crate::text::atoi(parts[p + 2]);
                if width < 1 {
                    return Err("Invalid extXYZ Properties".into());
                }
                if parts[p] == "species" && width == 1 {
                    species_column = column_count;
                    has_species = true;
                }
                if (parts[p] == "pos" || parts[p] == "positions") && width == 3 {
                    position_column = column_count;
                    has_positions = true;
                }
                column_count += width as usize;
            }
            if !has_species || !has_positions {
                return Err("extXYZ frames need species and pos columns".into());
            }
        }
        for _ in 0..count {
            checkpoint()?;
            let row = lines.get(i).ok_or("Truncated XYZ trajectory")?;
            i += 1;
            let tokens: Vec<&str> = row.split_whitespace().collect();
            if tokens.len() < column_count {
                return Err("Invalid XYZ atom row".into());
            }
            let symbol = standard_symbol(tokens[species_column]);
            let coordinate = |k: usize| crate::text::stod(tokens[position_column + k]).map_err(|_| "Invalid XYZ atom row".to_string());
            let p = [coordinate(0)?, coordinate(1)?, coordinate(2)?];
            if p.iter().any(|v| !v.is_finite()) {
                return Err("Invalid XYZ atom row".into());
            }
            let z = elements::atomic_number(&symbol);
            if z == 0 {
                return Err(format!("Unknown XYZ element: {symbol}"));
            }
            let mut atom = AtomSite::new(&symbol, p);
            atom.atomic_number = z;
            frame.atoms.push(atom);
        }
        frames.push(frame);
    }
    if frames.is_empty() {
        return Err("Trajectory contains no frames".into());
    }
    Ok(frames)
}

static TYPE_ELEMENTS: Mutex<(Vec<String>, String)> = Mutex::new((Vec::new(), String::new()));

/// The elements of LAMMPS atom types 1, 2, ... ("Cu Ni" or "Cu,Ni").
pub fn set_lammps_type_elements(symbols: &str) -> Result<(), String> {
    let text: String = symbols.chars().map(|c| if c == ',' || c == ';' { ' ' } else { c }).collect();
    let mut list = Vec::new();
    for word in text.split_whitespace() {
        let symbol = standard_symbol(word);
        if elements::atomic_number(&symbol) == 0 {
            return Err(format!("Unknown element in the LAMMPS types: {word}"));
        }
        list.push(symbol);
    }
    *TYPE_ELEMENTS.lock().unwrap() = (list, symbols.to_string());
    Ok(())
}

pub fn lammps_type_elements() -> String {
    TYPE_ELEMENTS.lock().unwrap().1.clone()
}

pub fn load_lammps_trajectory(path: &str) -> Result<Vec<Structure>, String> {
    let lines = read_lines(path)?;
    let mut iter = lines.iter().map(String::as_str);
    load_lammps_frames(&mut iter, usize::MAX)
}

/// Up to `max_frames` dump frames from the lines.
pub fn load_lammps_frames<'a>(lines: &mut dyn Iterator<Item = &'a str>, max_frames: usize) -> Result<Vec<Structure>, String> {
    let mut frames = Vec::new();
    let mut total = 0usize;
    let strip = |l: &'a str| l.strip_suffix('\r').unwrap_or(l);
    while frames.len() < max_frames {
        let Some(line) = lines.next() else { break };
        checkpoint()?;
        if blank(line) {
            continue;
        }
        if !line.starts_with("ITEM: TIMESTEP") {
            return Err(format!("Expected 'ITEM: TIMESTEP', got: {line}"));
        }
        lines.next().ok_or("Truncated LAMMPS dump: missing timestep value")?;
        if !lines.next().is_some_and(|l| l.starts_with("ITEM: NUMBER OF ATOMS")) {
            return Err("Expected 'ITEM: NUMBER OF ATOMS'".into());
        }
        let count_line = lines.next().ok_or("Truncated LAMMPS dump: missing atom count")?;
        let count = count_line.split_whitespace().next().map(crate::text::strtol).filter(|&(_, used)| used > 0).map(|(v, _)| v);
        let count = match count {
            Some(n) if (1..=1_000_000).contains(&n) => n as usize,
            _ => return Err("Invalid LAMMPS atom count".into()),
        };
        total += count;
        if total > 10_000_000 || frames.len() >= 100_000 {
            return Err("Trajectory exceeds memory limit".into());
        }
        let box_line = lines.next().filter(|l| l.starts_with("ITEM: BOX BOUNDS")).ok_or("Expected 'ITEM: BOX BOUNDS'")?;
        let triclinic = box_line.contains("xy xz yz");
        let (mut lo, mut hi, mut tilt) = ([0.0f64; 3], [0.0f64; 3], [0.0f64; 3]);
        for axis in 0..3 {
            let row = lines.next().ok_or("Truncated LAMMPS box bounds")?;
            let w: Vec<&str> = row.split_whitespace().collect();
            let get = |k: usize| w.get(k).and_then(|x| read_double(x));
            let ok = if triclinic {
                match (get(0), get(1), get(2)) {
                    (Some(a), Some(b), Some(c)) => {
                        lo[axis] = a;
                        hi[axis] = b;
                        tilt[axis] = c;
                        true
                    }
                    _ => false,
                }
            } else {
                match (get(0), get(1)) {
                    (Some(a), Some(b)) => {
                        lo[axis] = a;
                        hi[axis] = b;
                        true
                    }
                    _ => false,
                }
            };
            if !ok || !lo[axis].is_finite() || !hi[axis].is_finite() || !tilt[axis].is_finite() {
                return Err("Invalid LAMMPS box bounds".into());
            }
        }
        let header = lines.next().map(strip).filter(|l| l.starts_with("ITEM: ATOMS ")).ok_or("Expected 'ITEM: ATOMS <columns>'")?;
        let columns: Vec<&str> = header["ITEM: ATOMS ".len()..].split_whitespace().collect();
        let (mut element_col, mut x_col, mut y_col, mut z_col, mut id_col, mut type_col) = (-1i32, -1i32, -1i32, -1i32, -1i32, -1i32);
        for (i, &c) in columns.iter().enumerate() {
            let i = i as i32;
            match c {
                "element" => element_col = i,
                "id" => id_col = i,
                "type" => type_col = i,
                "x" | "xu" => x_col = i,
                "y" | "yu" => y_col = i,
                "z" | "zu" => z_col = i,
                _ => {}
            }
        }
        let mut scaled = false;
        if x_col < 0 || y_col < 0 || z_col < 0 {
            x_col = -1;
            y_col = -1;
            z_col = -1;
            for (i, &c) in columns.iter().enumerate() {
                let i = i as i32;
                match c {
                    "xs" | "xsu" => x_col = i,
                    "ys" | "ysu" => y_col = i,
                    "zs" | "zsu" => z_col = i,
                    _ => {}
                }
            }
            scaled = true;
        }
        let type_elements = TYPE_ELEMENTS.lock().unwrap().0.clone();
        if element_col < 0 && (type_col < 0 || type_elements.is_empty()) {
            return Err("LAMMPS dump has no 'element' column: name the elements of types 1, 2, ... \
                        (LAMMPS types in the Trajectory window, --types on the command line), \
                        or add 'dump_modify <id> element <symbols...>' to the LAMMPS input"
                .into());
        }
        if x_col < 0 || y_col < 0 || z_col < 0 {
            return Err("LAMMPS dump is missing x/y/z (or xu, xs, xsu) position columns".into());
        }
        let mut frame = Structure::default();
        let (xy, xz, yz) = if triclinic { (tilt[0], tilt[1], tilt[2]) } else { (0.0, 0.0, 0.0) };
        let min4 = 0.0f64.min(xy).min(xz).min(xy + xz);
        let max4 = 0.0f64.max(xy).max(xz).max(xy + xz);
        let (xlo, xhi) = (lo[0] - min4, hi[0] - max4);
        let (ylo, yhi) = (lo[1] - 0.0f64.min(yz), hi[1] - 0.0f64.max(yz));
        frame.cell_offset = [xlo, ylo, lo[2]];
        frame.cell_vectors = [[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, hi[2] - lo[2]]];
        frame.has_unit_cell = true;
        let max_column = [element_col, x_col, y_col, z_col, id_col, type_col].into_iter().max().unwrap();
        let mut rows: Vec<(i64, AtomSite)> = Vec::with_capacity(count);
        for i in 0..count {
            checkpoint()?;
            let line = lines.next().ok_or("Truncated LAMMPS trajectory")?;
            let values: Vec<&str> = line.split_whitespace().collect();
            if values.len() as i32 <= max_column {
                return Err("LAMMPS atom row has too few columns".into());
            }
            let symbol = if element_col >= 0 {
                standard_symbol(values[element_col as usize])
            } else {
                let word = values[type_col as usize];
                let (t, used) = crate::text::strtol(word);
                if used == 0 || t < 1 || t > type_elements.len() as i64 {
                    return Err(format!("LAMMPS atom type {word} has no element in the LAMMPS types ({} given)", type_elements.len()));
                }
                type_elements[(t - 1) as usize].clone()
            };
            let coordinate = |col: i32, axis: &str| {
                read_double(values[col as usize]).filter(|v| v.is_finite()).ok_or_else(|| format!("Invalid LAMMPS {axis} coordinate"))
            };
            let mut p = [coordinate(x_col, "x")?, coordinate(y_col, "y")?, coordinate(z_col, "z")?];
            let z = elements::atomic_number(&symbol);
            if z == 0 {
                return Err(format!("Unknown LAMMPS element: {symbol}"));
            }
            if scaled {
                let cv = frame.cell_vectors;
                let (fa, fb, fc) = (p[0], p[1], p[2]);
                for k in 0..3 {
                    p[k] = frame.cell_offset[k] + fa * cv[0][k] + fb * cv[1][k] + fc * cv[2][k];
                }
            }
            let mut atom = AtomSite::new(&symbol, p);
            atom.atomic_number = z;
            let id = if id_col >= 0 {
                let (v, used) = crate::text::strtol(values[id_col as usize]);
                if used == 0 {
                    return Err("Invalid LAMMPS atom id".into());
                }
                v
            } else {
                i as i64
            };
            rows.push((id, atom));
        }
        rows.sort_by_key(|r| r.0);
        frame.atoms = rows.into_iter().map(|r| r.1).collect();
        frames.push(frame);
    }
    if frames.is_empty() {
        return Err("Trajectory contains no frames".into());
    }
    Ok(frames)
}

/// XYZ/extXYZ or LAMMPS dump, from the first non-blank line.
pub fn load_trajectory(path: &str) -> Result<Vec<Structure>, String> {
    let lines = read_lines(path)?;
    let first = lines.iter().find(|l| !blank(l)).map(String::as_str).unwrap_or("");
    if first.starts_with("ITEM:") {
        load_lammps_trajectory(path)
    } else {
        load_xyz_trajectory(path)
    }
}
