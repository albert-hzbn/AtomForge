//! Volumetric files (the C++ electronic/VolumeIO): Gaussian Cube, VASP
//! CHGCAR/LOCPOT/ELFCAR-style grids and XSF DATAGRID_3D, read and written.

use std::path::Path;

use af_core::format::{format_g, stod};
use af_core::glm::{determinant_d, inverse_d, mul_d};
use glam::{DMat3, DVec3};

use crate::volume::{require_aligned, Grid, Site, Volume, BOHR, HARTREE};

fn atomic_number_from_symbol(symbol: &str) -> i32 {
    (1..=118).find(|&z| af_core::elements::symbol(z) == symbol).unwrap_or(0)
}

/// An istream over a text-mode file: std::getline and operator>> sharing one position.
struct Input {
    text: Vec<u8>,
    position: usize,
}

impl Input {
    fn new(bytes: &[u8]) -> Self {
        // Text mode: "\r\n" reads as "\n".
        let mut text = Vec::with_capacity(bytes.len());
        let mut i = 0;
        while i < bytes.len() {
            if bytes[i] == b'\r' && bytes.get(i + 1) == Some(&b'\n') {
                i += 1;
                continue;
            }
            text.push(bytes[i]);
            i += 1;
        }
        Self { text, position: 0 }
    }

    /// std::getline: None at the end with nothing read.
    fn getline(&mut self) -> Option<String> {
        if self.position >= self.text.len() {
            return None;
        }
        let rest = &self.text[self.position..];
        let end = rest.iter().position(|&c| c == b'\n');
        let line = String::from_utf8_lossy(&rest[..end.unwrap_or(rest.len())]).into_owned();
        self.position += end.map_or(rest.len(), |e| e + 1);
        Some(line)
    }

    fn line(&mut self) -> Result<String, String> {
        self.getline().ok_or_else(|| "Unexpected end of volumetric file".to_string())
    }

    /// operator>> into a string.
    fn next_token(&mut self) -> Option<String> {
        while self.position < self.text.len() && is_space(self.text[self.position]) {
            self.position += 1;
        }
        let start = self.position;
        while self.position < self.text.len() && !is_space(self.text[self.position]) {
            self.position += 1;
        }
        (self.position > start).then(|| String::from_utf8_lossy(&self.text[start..self.position]).into_owned())
    }

    fn token(&mut self) -> Result<String, String> {
        self.next_token().ok_or_else(|| "Truncated volumetric data".to_string())
    }
}

fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\x0B' | b'\x0C' | b'\r')
}

fn words(s: &str) -> Vec<String> {
    s.split(|c: char| c.is_ascii() && is_space(c as u8)).filter(|w| !w.is_empty()).map(str::to_string).collect()
}

fn number(s: &str) -> Result<f64, String> {
    let s: String = s.chars().map(|c| match c {
        'D' => 'E',
        'd' => 'e',
        c => c,
    }).collect();
    let (value, end) = stod(&s).map_err(|_| "stod".to_string())?;
    if end != s.len() || !value.is_finite() {
        return Err(format!("Invalid finite number: {s}"));
    }
    Ok(value)
}

/// std::stoi with the whole text consumed.
fn integer(s: &str) -> Result<i32, String> {
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() && is_space(bytes[i]) {
        i += 1;
    }
    let sign_start = i;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    let digits_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == digits_start {
        return Err("stoi".into());
    }
    let value: i64 = s[sign_start..i].parse().unwrap_or(if bytes[sign_start] == b'-' { i64::MIN } else { i64::MAX });
    if value < i64::from(i32::MIN) || value > i64::from(i32::MAX) {
        return Err("stoi".into());
    }
    if i != s.len() {
        return Err(format!("Invalid integer: {s}"));
    }
    Ok(value as i32)
}

fn vector(w: &[String], offset: usize) -> Result<DVec3, String> {
    if w.len() < offset + 3 {
        return Err("Missing vector components".into());
    }
    Ok(DVec3::new(number(&w[offset])?, number(&w[offset + 1])?, number(&w[offset + 2])?))
}

fn allocate(g: &mut Grid) -> Result<(), String> {
    let mut n: usize = 1;
    for &d in &g.shape {
        if !(2..=100000).contains(&d) || n > 100_000_000 / d as usize {
            return Err("Invalid or excessive grid dimensions".into());
        }
        n *= d as usize;
    }
    g.values.resize(n, 0.0);
    Ok(())
}

fn set_quantity(g: &mut Grid, quantity: &str, density_scale: f64, potential_scale: f64) -> Result<(), String> {
    let mut factor = 1.0;
    match quantity {
        "density" => {
            g.unit = "e/A^3".into();
            factor = density_scale;
        }
        "potential" => {
            g.unit = "eV".into();
            factor = potential_scale;
        }
        "elf" => g.unit = "1".into(),
        "raw" | "auto" => {}
        _ => return Err("Quantity must be auto, raw, density, potential or elf".into()),
    }
    for v in g.values.iter_mut() {
        *v *= factor;
    }
    g.validate()
}

fn set_column(m: &mut DMat3, i: usize, v: DVec3) {
    *m.col_mut(i) = v;
}

fn read_cube(input: &mut Input, quantity: &str, coordinate_unit: &str) -> Result<Volume, String> {
    let mut out = Volume::default();
    let mut g = Grid { name: input.line()?, ..Default::default() };
    input.line()?;
    let header = words(&input.line()?);
    if header.len() < 4 {
        return Err("Invalid Cube header".into());
    }
    let atoms = integer(&header[0])?;
    if !(-10_000_000..=10_000_000).contains(&atoms) {
        return Err("Excessive Cube atom count".into());
    }
    let mut channels = if header.len() > 4 { integer(&header[4])? } else { 1 };
    if atoms < 0 && channels != 1 {
        return Err("Orbital Cube NVAL must be 1".into());
    }
    if coordinate_unit != "bohr" && coordinate_unit != "angstrom" {
        return Err("Cube coordinate unit must be bohr or angstrom".into());
    }
    let length = if coordinate_unit == "bohr" { BOHR } else { 1.0 };
    g.origin = vector(&header, 1)? * length;
    for i in 0..3 {
        let w = words(&input.line()?);
        if w.len() < 4 {
            return Err("Invalid Cube axis".into());
        }
        // Negative dimension conventions vary among writers; units are explicit.
        let n = integer(&w[0])?;
        if n == 0 || !(-100000..=100000).contains(&n) {
            return Err("Invalid Cube dimension".into());
        }
        g.shape[i] = n.abs();
        set_column(&mut g.cell, i, vector(&w, 1)? * (length * (g.shape[i] - 1) as f64));
    }
    allocate(&mut g)?;
    for _ in 0..atoms.unsigned_abs() {
        let w = words(&input.line()?);
        if w.len() < 5 {
            return Err("Invalid Cube atom".into());
        }
        out.sites.push(Site { number: integer(&w[0])?, position: vector(&w, 2)? * length });
    }
    let mut names = Vec::new();
    if atoms < 0 {
        channels = integer(&input.token()?)?;
        if !(1..=1024).contains(&channels) {
            return Err("Invalid Cube channel count".into());
        }
        for _ in 0..channels {
            names.push(format!("orbital {}", integer(&input.token()?)?));
        }
    }
    if !(1..=1024).contains(&channels) || g.values.len() > 100_000_000 / channels as usize {
        return Err("Excessive Cube channels or samples".into());
    }
    out.fields = vec![g.clone(); channels as usize];
    for (c, field) in out.fields.iter_mut().enumerate() {
        field.name = if names.is_empty() { format!("field {c}") } else { names[c].clone() };
    }
    for x in 0..g.shape[0] {
        for y in 0..g.shape[1] {
            for z in 0..g.shape[2] {
                let index = g.index(x, y, z);
                for field in out.fields.iter_mut() {
                    field.values[index] = number(&input.token()?)?;
                }
            }
        }
    }
    if input.next_token().is_some() {
        return Err("Extra Cube values after expected data".into());
    }
    for field in out.fields.iter_mut() {
        set_quantity(field, quantity, 1.0 / (BOHR * BOHR * BOHR), HARTREE)?;
    }
    Ok(out)
}

fn read_vasp(input: &mut Input, quantity: &str) -> Result<Volume, String> {
    let mut out = Volume::default();
    let mut g = Grid { periodic: true, name: input.line()?, ..Default::default() };
    let scales = words(&input.line()?);
    if scales.len() != 1 && scales.len() != 3 {
        return Err("VASP requires one or three scale factors".into());
    }
    for i in 0..3 {
        let v = vector(&words(&input.line()?), 0)?;
        set_column(&mut g.cell, i, v);
    }
    let scaling;
    if scales.len() == 1 {
        let mut s = number(&scales[0])?;
        if s == 0.0 || determinant_d(&g.cell).abs() < 1e-12 {
            return Err("Invalid VASP lattice or scale".into());
        }
        if s < 0.0 {
            s = (-s / determinant_d(&g.cell).abs()).cbrt();
        }
        scaling = DVec3::splat(s);
    } else {
        scaling = vector(&scales, 0)?;
        if scaling.x <= 0.0 || scaling.y <= 0.0 || scaling.z <= 0.0 {
            return Err("VASP component scales must be positive".into());
        }
    }
    for i in 0..3 {
        let v = g.cell.col(i) * scaling;
        set_column(&mut g.cell, i, v);
    }
    let mut symbols = words(&input.line()?);
    if symbols.is_empty() {
        return Err("Missing VASP species".into());
    }
    let counts;
    if symbols[0].as_bytes()[0].is_ascii_digit() {
        counts = symbols.clone();
        symbols = vec!["X".to_string(); counts.len()];
    } else {
        counts = words(&input.line()?);
    }
    if counts.len() != symbols.len() {
        return Err("VASP species and counts differ".into());
    }
    let mut mode = words(&input.line()?);
    if !mode.is_empty() && (mode[0].starts_with('S') || mode[0].starts_with('s')) {
        mode = words(&input.line()?);
    }
    if mode.is_empty() {
        return Err("Missing VASP coordinate mode".into());
    }
    let m = mode[0].as_bytes()[0].to_ascii_lowercase();
    if m != b'c' && m != b'k' && m != b'd' {
        return Err("Unknown VASP coordinate mode".into());
    }
    for (species, count) in counts.iter().enumerate() {
        let count = integer(count)?;
        if !(0..=10_000_000).contains(&count) || out.sites.len() + count as usize > 10_000_000 {
            return Err("Invalid VASP atom count".into());
        }
        for _ in 0..count {
            let p = vector(&words(&input.line()?), 0)?;
            let p = if m == b'd' { mul_d(&g.cell, p) } else { scaling * p };
            out.sites.push(Site { number: atomic_number_from_symbol(&symbols[species]), position: p });
        }
    }
    while let Some(s) = input.getline() {
        let w = words(&s);
        if w.is_empty() {
            continue;
        }
        if w[0] == "augmentation" {
            if w.len() < 4 || w[1] != "occupancies" {
                return Err("Invalid VASP augmentation header".into());
            }
            let count = integer(&w[3])?;
            if !(0..=100_000_000).contains(&count) {
                return Err("Invalid augmentation count".into());
            }
            for _ in 0..count {
                number(&input.token()?)?;
            }
            if count != 0 {
                input.getline();
            }
            continue;
        }
        if w.len() != 3 {
            return Err("Unexpected VASP content after grid/atoms; trajectories are not supported".into());
        }
        for i in 0..3 {
            g.shape[i] = integer(&w[i])?;
        }
        allocate(&mut g)?;
        if !out.fields.is_empty() && g.shape != out.fields[0].shape {
            return Err("VASP channel dimensions differ".into());
        }
        for k in 0..g.values.len() {
            g.values[k] = number(&input.token()?)?;
        }
        let rest = input.getline().unwrap_or_default();
        if !words(&rest).is_empty() {
            return Err("Extra values on VASP final grid line".into());
        }
        let mut field = g.clone();
        let density_scale = 1.0 / determinant_d(&field.cell).abs();
        set_quantity(&mut field, quantity, density_scale, 1.0)?;
        g = field.clone();
        out.fields.push(field);
        if out.fields.len() > 4 {
            return Err("VASP supports at most four scalar channels".into());
        }
    }
    let n = out.fields.len();
    if n != 1 && n != 2 && n != 4 {
        return Err("Expected one, two or four VASP grids".into());
    }
    let names: &[&str] = if n == 4 { &["total", "mx", "my", "mz"] } else { &["total", "magnetization"] };
    for (field, name) in out.fields.iter_mut().zip(names) {
        field.name = name.to_string();
    }
    Ok(out)
}

fn read_xsf(input: &mut Input, quantity: &str) -> Result<Volume, String> {
    let mut out = Volume::default();
    let mut atoms_section = false;
    let mut total_samples: usize = 0;
    let z_of = |w: &str| -> Result<i32, String> { if w.as_bytes()[0].is_ascii_digit() { integer(w) } else { Ok(atomic_number_from_symbol(w)) } };
    while let Some(s) = input.getline() {
        let w = words(&s);
        if w.is_empty() || w[0].starts_with('#') {
            continue;
        }
        if atoms_section {
            let z = z_of(&w[0])?;
            if w.len() >= 4 && (z > 0 || w[0] == "0") {
                if out.sites.len() >= 10_000_000 {
                    return Err("Excessive XSF atom count".into());
                }
                out.sites.push(Site { number: z, position: vector(&w, 1)? });
                continue;
            }
            atoms_section = false;
        }
        if w[0] == "ATOMS" {
            atoms_section = true;
            continue;
        }
        if w[0] == "ANIMSTEPS" || w[0].contains("BANDGRID") {
            return Err("XSF animations and reciprocal-space BANDGRID are not real-space scalar grids".into());
        }
        if w[0] == "PRIMCOORD" {
            let counts = words(&input.line()?);
            if counts.is_empty() {
                return Err("Missing XSF atom count".into());
            }
            let count = integer(&counts[0])?;
            if !(0..=10_000_000).contains(&count) {
                return Err("Invalid XSF atom count".into());
            }
            for _ in 0..count {
                let atom = words(&input.line()?);
                if atom.len() < 4 {
                    return Err("Invalid XSF atom".into());
                }
                out.sites.push(Site { number: z_of(&atom[0])?, position: vector(&atom, 1)? });
            }
        } else if w[0].starts_with("BEGIN_DATAGRID_3D") {
            let mut g = Grid { name: if w[0].len() > 17 { w[0][17..].to_string() } else { "field".into() }, ..Default::default() };
            let dimensions = words(&input.line()?);
            if dimensions.len() != 3 {
                return Err("Invalid XSF dimensions".into());
            }
            for i in 0..3 {
                g.shape[i] = integer(&dimensions[i])?;
            }
            g.origin = vector(&words(&input.line()?), 0)?;
            for i in 0..3 {
                let v = vector(&words(&input.line()?), 0)?;
                set_column(&mut g.cell, i, v);
            }
            allocate(&mut g)?;
            total_samples += g.values.len();
            if total_samples > 100_000_000 {
                return Err("Excessive total XSF samples".into());
            }
            for k in 0..g.values.len() {
                g.values[k] = number(&input.token()?)?;
            }
            if input.token()? != "END_DATAGRID_3D" {
                return Err("Missing XSF grid terminator or extra values".into());
            }
            set_quantity(&mut g, quantity, 1.0, 1.0)?;
            out.fields.push(g);
        }
    }
    if out.fields.is_empty() {
        return Err("No XSF DATAGRID_3D found".into());
    }
    Ok(out)
}

/// Loads a volumetric file: .cube/.cub, .xsf, else VASP (CHGCAR, LOCPOT, ELFCAR, ...).
/// `quantity`: auto, raw, density, potential or elf; `coordinate_unit` (Cube): bohr or angstrom.
pub fn load_volume(path: &str, quantity: &str, coordinate_unit: &str) -> Result<Volume, String> {
    let bytes = std::fs::read(path).map_err(|_| format!("Cannot open volumetric file: {path}"))?;
    let p = Path::new(path);
    let ext = p.extension().map(|e| format!(".{}", e.to_string_lossy().to_ascii_lowercase())).unwrap_or_default();
    let mut input = Input::new(&bytes);
    if ext == ".cube" || ext == ".cub" {
        return read_cube(&mut input, quantity, coordinate_unit);
    }
    if ext == ".xsf" {
        return read_xsf(&mut input, quantity);
    }
    let mut q = quantity.to_string();
    if q == "auto" {
        let name = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        q = if name.starts_with("CHGCAR") || name == "CHG" || name.starts_with("PARCHG") || name.starts_with("AECCAR") {
            "density"
        } else if name.starts_with("LOCPOT") {
            "potential"
        } else if name.starts_with("ELFCAR") {
            "elf"
        } else {
            "raw"
        }
        .to_string();
    }
    read_vasp(&mut input, &q)
}

fn g17(x: f64) -> String {
    format_g(x, 17)
}

fn write_vector(out: &mut String, p: DVec3) {
    out.push_str(&format!("{} {} {}\n", g17(p.x), g17(p.y), g17(p.z)));
}

/// Saves the fields (aligned) and sites as cube, xsf or vasp.
pub fn save_volume(volume: &Volume, path: &str, format: &str) -> Result<(), String> {
    if volume.fields.is_empty() {
        return Err("No fields to save".into());
    }
    for field in &volume.fields {
        require_aligned(&volume.fields[0], field)?;
    }
    if format != "cube" && format != "xsf" && format != "vasp" {
        return Err("Output format must be cube, xsf or vasp".into());
    }
    let first = &volume.fields[0];
    if format == "vasp" && (!first.periodic || first.origin.length() > 1e-10) {
        return Err("VASP output requires a periodic grid with zero origin".into());
    }
    if format == "vasp" && ![1, 2, 4].contains(&volume.fields.len()) {
        return Err("VASP requires one, two or four fields".into());
    }
    if format == "vasp" && volume.sites.is_empty() {
        return Err("VASP output requires atomic sites".into());
    }
    let mut out = String::new();
    if format == "cube" {
        // Add periodic endpoint planes so finite-grid interchange preserves extent.
        let mut shape = first.shape;
        if first.periodic {
            for n in shape.iter_mut() {
                *n += 1;
            }
        }
        out += &format!("AtomForge volume\nUnits: atomic units for density/potential; raw otherwise\n{} ", volume.sites.len());
        let p = first.origin / BOHR;
        out += &format!("{} {} {} {}\n", g17(p.x), g17(p.y), g17(p.z), volume.fields.len());
        let steps = first.steps();
        for i in 0..3 {
            out += &format!("{} ", shape[i]);
            write_vector(&mut out, steps.col(i) / BOHR);
        }
        for site in &volume.sites {
            out += &format!("{} 0 ", site.number);
            write_vector(&mut out, site.position / BOHR);
        }
        for x in 0..shape[0] {
            for y in 0..shape[1] {
                for z in 0..shape[2] {
                    for f in &volume.fields {
                        let factor = if f.unit == "e/A^3" {
                            BOHR * BOHR * BOHR
                        } else if f.unit == "eV" {
                            1.0 / HARTREE
                        } else {
                            1.0
                        };
                        out += &g17(f.values[f.index(x % f.shape[0], y % f.shape[1], z % f.shape[2])] * factor);
                        out.push(' ');
                    }
                    out.push('\n');
                }
            }
        }
    } else if format == "xsf" {
        out += "CRYSTAL\nPRIMVEC\n";
        for i in 0..3 {
            write_vector(&mut out, first.cell.col(i));
        }
        out += &format!("PRIMCOORD\n{} 1\n", volume.sites.len());
        for site in &volume.sites {
            out += &format!("{} ", site.number);
            write_vector(&mut out, site.position);
        }
        out += "BEGIN_BLOCK_DATAGRID_3D\nAtomForge\n";
        for (c, f) in volume.fields.iter().enumerate() {
            let mut shape = f.shape;
            if f.periodic {
                for n in shape.iter_mut() {
                    *n += 1;
                }
            }
            out += &format!("BEGIN_DATAGRID_3D_field{c}\n{} {} {}\n", shape[0], shape[1], shape[2]);
            write_vector(&mut out, f.origin);
            for i in 0..3 {
                write_vector(&mut out, f.cell.col(i));
            }
            for z in 0..shape[2] {
                for y in 0..shape[1] {
                    for x in 0..shape[0] {
                        out += &g17(f.values[f.index(x % f.shape[0], y % f.shape[1], z % f.shape[2])]);
                        out.push('\n');
                    }
                }
            }
            out += "END_DATAGRID_3D\n";
        }
        out += "END_BLOCK_DATAGRID_3D\n";
    } else {
        out += "AtomForge post-processing grid (no PAW augmentation)\n1.0\n";
        for i in 0..3 {
            write_vector(&mut out, first.cell.col(i));
        }
        // Each site is its own species group, preserving site order without loss.
        for site in &volume.sites {
            out += &format!("{} ", af_core::elements::symbol(site.number));
        }
        out.push('\n');
        for _ in &volume.sites {
            out += "1 ";
        }
        out += "\nDirect\n";
        let inverse = inverse_d(&first.cell);
        for site in &volume.sites {
            write_vector(&mut out, mul_d(&inverse, site.position));
        }
        for f in &volume.fields {
            out += &format!("\n{} {} {}\n", f.shape[0], f.shape[1], f.shape[2]);
            let factor = if f.unit == "e/A^3" { determinant_d(&f.cell).abs() } else { 1.0 };
            for (count, v) in f.values.iter().enumerate() {
                out += &g17(v * factor);
                out.push(if (count + 1) % 5 != 0 { ' ' } else { '\n' });
            }
            out.push('\n');
        }
    }
    // The C++ writes through a text-mode stream: CRLF line ends on Windows.
    if cfg!(windows) {
        out = out.replace('\n', "\r\n");
    }
    std::fs::write(path, out).map_err(|_| format!("Cannot write volume: {path}"))
}
