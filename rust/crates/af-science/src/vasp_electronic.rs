//! VASP electronic structure (the C++ science/VaspElectronic): EIGENVAL,
//! DOSCAR, line-mode KPOINTS and PROCAR readers, and the "vasp-electronic"
//! tool (bands, gap, DOS, fat-band weights).

use std::collections::BTreeMap;
use std::f64::consts::TAU;
use std::path::Path;

use af_core::format::stod;
use af_core::json::Json;

use crate::analysis::{band_gap, jv};
use crate::core::*;
use crate::data::{read_text_file, Parameters};
use crate::tools::ToolOutput;

/// std::getline over a text: lines without the '\n' (a '\r' stays).
struct Lines<'a> {
    text: &'a str,
    position: usize,
}

impl<'a> Lines<'a> {
    fn new(text: &'a str) -> Self {
        Self { text, position: 0 }
    }

    fn next_line(&mut self) -> Option<&'a str> {
        if self.position >= self.text.len() {
            return None;
        }
        let rest = &self.text[self.position..];
        match rest.find('\n') {
            Some(end) => {
                self.position += end + 1;
                Some(&rest[..end])
            }
            None => {
                self.position = self.text.len();
                Some(rest)
            }
        }
    }

    /// An unchecked getline: an empty line at the end.
    fn line(&mut self) -> &'a str {
        self.next_line().unwrap_or("")
    }
}

/// Whitespace-separated words (istream >> std::string).
fn tokens(line: &str) -> Vec<&str> {
    line.split(|c: char| matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')).filter(|t| !t.is_empty()).collect()
}

fn number(token: &str, what: &str) -> Result<f64, String> {
    stod(token).map(|(v, _)| v).map_err(|_| format!("Invalid number in {what}: {token}"))
}

fn is_blank(line: &str) -> bool {
    line.chars().all(|c| matches!(c, ' ' | '\t' | '\r'))
}

pub struct Eigenvalues {
    /// Fractional.
    pub kpoints: Vec<Vec3>,
    pub weights: Vec<f64>,
    /// energies[spin][k][band] (eV).
    pub energies: Vec<Vec<Vec<f64>>>,
}

pub fn read_eigenval(path: &Path) -> Result<Eigenvalues, String> {
    let text = read_text_file(path)?;
    let mut input = Lines::new(&text);
    let first = tokens(input.line());
    if first.len() < 4 {
        return Err("Invalid EIGENVAL header".into());
    }
    let spins = number(first[3], "EIGENVAL")? as i32;
    if spins != 1 && spins != 2 {
        return Err("EIGENVAL ISPIN must be 1 or 2".into());
    }
    for _ in 0..4 {
        input.line();
    }
    let sizes = tokens(input.line());
    if sizes.len() < 3 {
        return Err("Invalid EIGENVAL dimensions".into());
    }
    let kpoints = number(sizes[1], "EIGENVAL")? as i32;
    let bands = number(sizes[2], "EIGENVAL")? as i32;
    let (spins, kn, bn) = (spins as usize, kpoints.max(0) as usize, bands.max(0) as usize);
    let mut result = Eigenvalues { kpoints: Vec::new(), weights: Vec::new(), energies: vec![vec![vec![0.0; bn]; kn]; spins] };
    for k in 0..kn {
        let line = loop {
            let line = input.next_line().ok_or("Truncated EIGENVAL")?;
            if !is_blank(line) {
                break line;
            }
        };
        let kline = tokens(line);
        if kline.len() < 4 {
            return Err("Invalid EIGENVAL k point".into());
        }
        result.kpoints.push([number(kline[0], "EIGENVAL")?, number(kline[1], "EIGENVAL")?, number(kline[2], "EIGENVAL")?]);
        result.weights.push(number(kline[3], "EIGENVAL")?);
        for b in 0..bn {
            let row = tokens(input.next_line().ok_or("Truncated EIGENVAL")?);
            if row.len() < 1 + spins {
                return Err("Invalid EIGENVAL band row".into());
            }
            for s in 0..spins {
                result.energies[s][k][b] = number(row[1 + s], "EIGENVAL")?;
            }
        }
    }
    Ok(result)
}

#[derive(Default)]
pub struct DensityOfStates {
    pub fermi: f64,
    /// eV.
    pub energy: Vec<f64>,
    /// [spin][point].
    pub total: Vec<Vec<f64>>,
    /// Site-projected DOS summed over orbitals: [atom][spin][point]; empty when absent.
    pub projected: Vec<Vec<Vec<f64>>>,
}

pub fn read_doscar(path: &Path) -> Result<DensityOfStates, String> {
    let text = read_text_file(path)?;
    let mut input = Lines::new(&text);
    let first = tokens(input.line());
    if first.len() < 3 {
        return Err("Invalid DOSCAR header".into());
    }
    let atoms = number(first[0], "DOSCAR")? as i32;
    let partial = number(first[2], "DOSCAR")? != 0.0;
    for _ in 0..4 {
        input.line();
    }
    let header = tokens(input.line());
    if header.len() < 4 {
        return Err("Invalid DOSCAR energy header".into());
    }
    let points = number(header[2], "DOSCAR")? as i32;
    let mut dos = DensityOfStates { fermi: number(header[3], "DOSCAR")?, ..Default::default() };
    let mut spins = 0usize;
    for i in 0..points {
        let row = tokens(input.next_line().ok_or("Truncated DOSCAR")?);
        if i == 0 {
            spins = if row.len() >= 5 { 2 } else { 1 };
            dos.total = vec![Vec::new(); spins];
        }
        dos.energy.push(number(row.first().copied().unwrap_or(""), "DOSCAR")?);
        for s in 0..spins {
            dos.total[s].push(number(row.get(1 + s).copied().unwrap_or(""), "DOSCAR")?);
        }
    }
    if partial {
        for _ in 0..atoms {
            match input.next_line() {
                Some(line) if !is_blank(line) => {} // header of the site block
                _ => break,
            }
            let mut site = vec![vec![0.0; points.max(0) as usize]; spins];
            for i in 0..points.max(0) as usize {
                let row = tokens(input.next_line().ok_or("Truncated DOSCAR site block")?);
                // Spin-polarized orbital columns alternate up/down.
                for (c, token) in row.iter().enumerate().skip(1) {
                    site[if spins == 2 { (c - 1) % 2 } else { 0 }][i] += number(token, "DOSCAR")?;
                }
            }
            dos.projected.push(site);
        }
    }
    Ok(dos)
}

/// The labels of a line-mode KPOINTS file and its points per segment.
pub fn read_line_mode_labels(path: &Path) -> Result<(Vec<String>, i32), String> {
    let text = read_text_file(path)?;
    let mut input = Lines::new(&text);
    input.line();
    let count = tokens(input.line());
    if count.is_empty() {
        return Err("Invalid KPOINTS file".into());
    }
    let points_per_segment = number(count[0], "KPOINTS")? as i32;
    let mode = input.line();
    match mode.find(|c| c != ' ' && c != '\t') {
        Some(i) if mode.as_bytes()[i].to_ascii_lowercase() == b'l' => {}
        _ => return Err("KPOINTS is not in line mode".into()),
    }
    input.line(); // Reciprocal / Cartesian
    let mut labels = Vec::new();
    while let Some(full) = input.next_line() {
        if is_blank(full) {
            continue;
        }
        let mut line = full;
        let mut label = String::new();
        if let Some(bang) = line.find('!') {
            if let Some(first) = tokens(&line[bang + 1..]).first() {
                label = first.to_string();
            }
            line = &line[..bang];
        }
        let row = tokens(line);
        if row.len() < 3 {
            continue;
        }
        if label.is_empty() && row.len() >= 4 {
            label = row[3].to_string();
        }
        if label == "\\Gamma" || label == "Gamma" || label == "G" || label == "GAMMA" {
            label = "G".into();
        }
        labels.push(label);
    }
    Ok((labels, points_per_segment))
}

pub struct Projections {
    pub orbitals: Vec<String>,
    /// weight[spin][k][band][ion][orbital].
    pub weight: Vec<Vec<Vec<Vec<Vec<f64>>>>>,
}

pub fn read_procar(path: &Path) -> Result<Projections, String> {
    let text = read_text_file(path)?;
    let mut input = Lines::new(&text);
    let mut result = Projections { orbitals: Vec::new(), weight: Vec::new() };
    let (mut spin, mut k, mut band) = (-1i64, -1i64, -1i64);
    let (mut in_block, mut block_done) = (false, false);
    while let Some(line) = input.next_line() {
        let row = tokens(line);
        if row.is_empty() {
            continue;
        }
        if line.contains("# of k-points") {
            spin += 1;
            result.weight.push(Vec::new());
            k = -1;
            continue;
        }
        if row[0] == "k-point" {
            k += 1;
            if spin >= 0 {
                result.weight[spin as usize].push(Vec::new());
            }
            band = -1;
            continue;
        }
        if row[0] == "band" {
            band += 1;
            if spin >= 0 && k >= 0 {
                result.weight[spin as usize][k as usize].push(Vec::new());
            }
            in_block = false;
            block_done = false;
            continue;
        }
        if row[0] == "ion" {
            // The first ion block of a band holds the projections (a later
            // one, with LORBIT = 12, holds phases).
            if !block_done {
                in_block = true;
                if result.orbitals.is_empty() && row.len() >= 2 {
                    result.orbitals = row[1..row.len() - 1].iter().map(|s| s.to_string()).collect();
                }
            }
            continue;
        }
        if row[0] == "tot" {
            if in_block {
                block_done = true;
            }
            in_block = false;
            continue;
        }
        if in_block && spin >= 0 && k >= 0 && band >= 0 && row[0].as_bytes()[0].is_ascii_digit() {
            let mut values = Vec::new();
            let mut c = 1;
            while c < row.len() - 1 && values.len() < result.orbitals.len() {
                values.push(number(row[c], "PROCAR")?);
                c += 1;
            }
            result.weight[spin as usize][k as usize][band as usize].push(values);
        }
    }
    if result.weight.is_empty() || result.orbitals.is_empty() {
        return Err("No projections found in PROCAR".into());
    }
    Ok(result)
}

pub fn vasp_electronic(p: &Parameters) -> Result<ToolOutput, String> {
    let eigen = read_eigenval(Path::new(p.json("eigenval_file")?.string()?))?;
    let (spins, kpoints) = (eigen.energies.len(), eigen.kpoints.len());
    let bands = eigen.energies.first().and_then(|s| s.first()).map_or(0, |k| k.len());
    let mut result = Json::object();
    let mut fermi_source = "none";
    let mut fermi = 0.0;
    let have_dos = p.has("doscar_file");
    let dos = if have_dos { read_doscar(Path::new(p.json("doscar_file")?.string()?))? } else { DensityOfStates::default() };
    if p.has("fermi_eV") {
        fermi = p.number("fermi_eV")?;
        fermi_source = "input";
    } else if have_dos {
        fermi = dos.fermi;
        fermi_source = "DOSCAR";
    }
    // Reciprocal lattice (2 pi convention) from the structure, else fractional distances.
    let mut reciprocal = identity();
    let mut symbols: Vec<String> = Vec::new();
    if p.has("structure") {
        let structure = &p.structure("structure")?.structure;
        if !structure.has_unit_cell {
            return Err("The structure needs a cell for k-space distances".into());
        }
        let inv = inverse(&structure.cell_vectors)?;
        for i in 0..3 {
            for j in 0..3 {
                reciprocal[i][j] = TAU * inv[j][i];
            }
        }
        symbols = structure.atoms.iter().map(|a| a.symbol.clone()).collect();
        result.set("distance_unit", "1/Angstrom (2 pi included)");
    } else {
        result.set("distance_unit", "fractional reciprocal coordinates");
    }
    // Path distances, with line-mode labels and path breaks.
    let mut point_labels = vec![String::new(); kpoints];
    let mut break_before = vec![false; kpoints];
    if p.has("kpoints_file") {
        let (labels, per_segment) = read_line_mode_labels(Path::new(p.json("kpoints_file")?.string()?))?;
        let segments = labels.len() / 2;
        if per_segment > 1 && segments * per_segment as usize == kpoints {
            for s in 0..segments {
                let start = s * per_segment as usize;
                let end = start + per_segment as usize - 1;
                point_labels[start] = labels[2 * s].clone();
                point_labels[end] = labels[2 * s + 1].clone();
                if s > 0 && labels[2 * s] != labels[2 * s - 1] {
                    break_before[start] = true;
                }
            }
        } else {
            result.set("kpoints_note", "KPOINTS segments do not match the EIGENVAL k points; labels omitted");
        }
    }
    let mut distance = vec![0.0; kpoints];
    for k in 1..kpoints {
        distance[k] = distance[k - 1] + if break_before[k] { 0.0 } else { norm(row_times(sub(eigen.kpoints[k], eigen.kpoints[k - 1]), &reciprocal)) };
    }
    let mut markers: Vec<(f64, String)> = Vec::new();
    for k in 0..kpoints {
        if point_labels[k].is_empty() {
            continue;
        }
        let mut label = point_labels[k].clone();
        if k + 1 < kpoints && break_before[k + 1] && !point_labels[k + 1].is_empty() {
            label += &format!("|{}", point_labels[k + 1]);
        }
        if break_before[k] {
            continue; // merged into the previous point's label
        }
        // A continuous junction repeats its point at the end and start of two segments.
        if let Some(last) = markers.last() {
            if (last.0 - distance[k]).abs() < 1e-12 && last.1 == label {
                continue;
            }
        }
        markers.push((distance[k], label));
    }
    let markers_json: Vec<Json> = markers
        .into_iter()
        .map(|(d, l)| {
            let mut marker = Json::object();
            marker.set("distance", d);
            marker.set("label", l);
            marker
        })
        .collect();
    let energies = Json::Array(
        eigen.energies.iter().map(|spin| Json::Array(spin.iter().map(|row| jv(&row.iter().map(|e| e - fermi).collect::<Vec<_>>())).collect())).collect(),
    );
    result.set("fermi_eV", fermi);
    result.set("fermi_source", fermi_source);
    result.set("spins", spins);
    result.set("kpoints", kpoints);
    result.set("bands", bands);
    result.set("distance", jv(&distance));
    result.set("labels", Json::Array(markers_json));
    result.set("energies_minus_fermi_eV", energies.clone());
    if fermi_source != "none" {
        let gap = (|| {
            let mut request = Json::object();
            request.set("energies_eV", energies.clone());
            request.set("fermi_eV", 0.0);
            band_gap(&Parameters::new("band-gap", &request, Path::new("."), None)?)
        })();
        match gap {
            Ok(g) => result.set("gap", g),
            Err(error) => {
                result.set("gap", Json::Null);
                result.set("gap_note", error);
            }
        }
    }
    if have_dos {
        let mut dos_json = Json::object();
        dos_json.set("energy_minus_fermi_eV", jv(&dos.energy.iter().map(|e| e - fermi).collect::<Vec<_>>()));
        dos_json.set("total", Json::Array(dos.total.iter().map(|s| jv(s)).collect()));
        if !dos.projected.is_empty() && dos.projected.len() == symbols.len() {
            let mut by_element: BTreeMap<&str, Vec<Vec<f64>>> = BTreeMap::new();
            for (a, symbol) in symbols.iter().enumerate() {
                let sum = by_element.entry(symbol).or_default();
                if sum.is_empty() {
                    *sum = dos.projected[a].clone();
                } else {
                    for s in 0..sum.len() {
                        for i in 0..sum[s].len() {
                            sum[s][i] += dos.projected[a][s][i];
                        }
                    }
                }
            }
            let mut projected = Json::object();
            for (element, spins_data) in &by_element {
                projected.set(element, Json::Array(spins_data.iter().map(|s| jv(s)).collect()));
            }
            dos_json.set("projected_by_element", projected);
        }
        result.set("dos", dos_json);
    }
    if p.has("procar_file") {
        let projections = read_procar(Path::new(p.json("procar_file")?.string()?))?;
        if projections.weight.len() != spins || projections.weight[0].len() != kpoints {
            return Err("PROCAR k points or spins do not match EIGENVAL".into());
        }
        let orbital = if p.has("projection_orbital") { p.json("projection_orbital")?.string()?.to_string() } else { "all".to_string() };
        let elements: Vec<String> = if p.has("projection_element") {
            vec![p.json("projection_element")?.string()?.to_string()]
        } else {
            let mut e = symbols.clone();
            e.sort();
            e.dedup();
            e
        };
        if elements.is_empty() {
            return Err("Fat bands need the structure (for element names) or projection_element".into());
        }
        let selected_orbital = |name: &str| {
            orbital == "all" || name == orbital || (orbital.len() == 1 && !name.is_empty() && name.as_bytes()[0] == orbital.as_bytes()[0] && name != "tot")
        };
        let mut fat = Json::object();
        for element in &elements {
            let mut by_spin = Vec::new();
            for s in 0..spins {
                let mut by_k = Vec::new();
                for k in 0..kpoints {
                    let mut row = Vec::new();
                    for ions in &projections.weight[s][k] {
                        let (mut selected, mut total) = (0.0, 0.0);
                        for (ion, values) in ions.iter().enumerate() {
                            for (o, &v) in values.iter().enumerate() {
                                total += v;
                                if ion < symbols.len() && &symbols[ion] == element && selected_orbital(&projections.orbitals[o]) {
                                    selected += v;
                                }
                            }
                        }
                        row.push(if total > 0.0 { selected / total } else { 0.0 });
                    }
                    by_k.push(jv(&row));
                }
                by_spin.push(Json::Array(by_k));
            }
            fat.set(element, Json::Array(by_spin));
        }
        result.set("fat_band_orbital", orbital);
        result.set("fat_band_weights", fat);
        result.set("procar_orbitals", Json::Array(projections.orbitals.iter().map(|o| Json::from(o.as_str())).collect()));
    }
    Ok(ToolOutput { result, ..Default::default() })
}
