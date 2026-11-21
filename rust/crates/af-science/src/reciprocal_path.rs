//! HPKOT band paths through the Brillouin zone and DFT band-structure inputs
//! (the C++ science/ReciprocalPath): spglib standardisation, the extended
//! Bravais lattice, the HPKOT primitive cell and SeeK-path's special points.

use std::collections::BTreeMap;
use std::f64::consts::TAU;
use std::fmt::Write as _;

use af_core::elements;
use af_core::json::Json;
use af_core::structure::{AtomSite, Structure};

use crate::analysis::jv;
use crate::core::*;
use crate::data::{atomic_mass, structure_json, StructureInput};
use crate::hpkot_data::LATTICES;
use crate::tools::ToolOutput;

/// Arithmetic over the named cell and k-vector parameters of the HPKOT tables.
struct Expression<'a> {
    text: &'a [u8],
    source: &'a str,
    variables: &'a BTreeMap<String, f64>,
    position: usize,
}

impl<'a> Expression<'a> {
    fn evaluate(source: &'a str, variables: &'a BTreeMap<String, f64>) -> Result<f64, String> {
        let mut e = Expression { text: source.as_bytes(), source, variables, position: 0 };
        let value = e.sum()?;
        e.skip();
        if e.position != e.text.len() {
            return Err(e.invalid());
        }
        Ok(value)
    }

    fn invalid(&self) -> String {
        format!("Invalid k-path expression: {}", self.source)
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.position).copied()
    }

    fn skip(&mut self) {
        while self.peek() == Some(b' ') {
            self.position += 1;
        }
    }

    fn sum(&mut self) -> Result<f64, String> {
        let mut value = self.product()?;
        loop {
            self.skip();
            match self.peek() {
                Some(op @ (b'+' | b'-')) => {
                    self.position += 1;
                    let rhs = self.product()?;
                    value = if op == b'+' { value + rhs } else { value - rhs };
                }
                _ => return Ok(value),
            }
        }
    }

    fn product(&mut self) -> Result<f64, String> {
        let mut value = self.factor()?;
        loop {
            self.skip();
            match self.peek() {
                Some(op @ (b'*' | b'/')) => {
                    self.position += 1;
                    let rhs = self.factor()?;
                    value = if op == b'*' { value * rhs } else { value / rhs };
                }
                _ => return Ok(value),
            }
        }
    }

    /// The longest decimal prefix (strtod without hex, inf or nan).
    fn number(&mut self) -> Result<f64, String> {
        let t = self.text;
        let start = self.position;
        let mut i = start;
        let digits = |i: &mut usize| {
            let from = *i;
            while *i < t.len() && t[*i].is_ascii_digit() {
                *i += 1;
            }
            *i - from
        };
        let mut count = digits(&mut i);
        if i < t.len() && t[i] == b'.' {
            i += 1;
            count += digits(&mut i);
        }
        if count == 0 {
            return Err("stod".into());
        }
        if i < t.len() && (t[i] == b'e' || t[i] == b'E') {
            let mut j = i + 1;
            if j < t.len() && (t[j] == b'+' || t[j] == b'-') {
                j += 1;
            }
            if digits(&mut j) > 0 {
                i = j;
            }
        }
        self.position = i;
        std::str::from_utf8(&t[start..i]).unwrap().parse::<f64>().map_err(|_| "stod".to_string())
    }

    fn factor(&mut self) -> Result<f64, String> {
        self.skip();
        let Some(c) = self.peek() else {
            return Err(self.invalid());
        };
        if c == b'-' || c == b'+' {
            self.position += 1;
            let value = self.factor()?;
            return Ok(if c == b'-' { -value } else { value });
        }
        if c == b'(' {
            self.position += 1;
            let value = self.sum()?;
            self.skip();
            if self.peek() != Some(b')') {
                return Err(self.invalid());
            }
            self.position += 1;
            return Ok(value);
        }
        if c.is_ascii_digit() || c == b'.' {
            return self.number();
        }
        let start = self.position;
        while self.peek().is_some_and(|c| c.is_ascii_alphanumeric() || c == b'_') {
            self.position += 1;
        }
        let name = &self.source[start..self.position];
        match self.variables.get(name) {
            Some(v) if !name.is_empty() => Ok(*v),
            _ => Err(format!("Unknown k-path symbol '{name}' in {}", self.source)),
        }
    }
}

#[derive(Clone, Copy)]
struct CellParameters {
    a: f64,
    b: f64,
    c: f64,
    cos_alpha: f64,
    cos_beta: f64,
    cos_gamma: f64,
}

fn parameters(cell: &Mat3) -> CellParameters {
    let (a, b, c) = (norm(cell[0]), norm(cell[1]), norm(cell[2]));
    CellParameters {
        a,
        b,
        c,
        cos_alpha: dot(cell[1], cell[2]) / b / c,
        cos_beta: dot(cell[0], cell[2]) / a / c,
        cos_gamma: dot(cell[0], cell[1]) / a / b,
    }
}

/// Rows are reciprocal vectors with real . reciprocal^T = 2 pi I.
fn reciprocal_rows(real: &Mat3) -> Result<Mat3, String> {
    let inv = inverse(real)?;
    Ok(std::array::from_fn(|i| std::array::from_fn(|j| TAU * inv[j][i])))
}

fn integer_matrix(v: [i32; 9]) -> Mat3 {
    std::array::from_fn(|r| std::array::from_fn(|c| v[r * 3 + c] as f64))
}

/// HPKOT Table 3: conventional (a, b, c) -> primitive (a, b, c) P, and P^-1.
fn primitive_matrices(bravais: &str) -> (Mat3, Mat3) {
    let scaled = |v: [i32; 9], s: f64| {
        let mut m = integer_matrix(v);
        for row in m.iter_mut() {
            *row = scale(*row, s);
        }
        m
    };
    match bravais {
        "cF" | "oF" => (scaled([0, 1, 1, 1, 0, 1, 1, 1, 0], 0.5), integer_matrix([-1, 1, 1, 1, -1, 1, 1, 1, -1])),
        "cI" | "tI" | "oI" => (scaled([-1, 1, 1, 1, -1, 1, 1, 1, -1], 0.5), integer_matrix([0, 1, 1, 1, 0, 1, 1, 1, 0])),
        "hR" => (scaled([2, -1, -1, 1, 1, -2, 1, 1, 1], 1.0 / 3.0), integer_matrix([1, 0, 1, -1, 1, 1, 0, -1, 1])),
        "oC" => (scaled([1, 1, 0, -1, 1, 0, 0, 0, 2], 0.5), integer_matrix([1, -1, 0, 1, 1, 0, 0, 0, 1])),
        "oA" => (scaled([0, 0, 2, 1, 1, 0, -1, 1, 0], 0.5), integer_matrix([0, 1, -1, 0, 1, 1, 1, 0, 0])),
        "mC" => (scaled([1, -1, 0, 1, 1, 0, 0, 0, 2], 0.5), integer_matrix([1, 1, 0, -1, 1, 0, 0, 0, 1])),
        _ => (identity(), identity()),
    }
}

fn centrosymmetric(point_group: &str) -> bool {
    ["-1", "2/m", "mmm", "4/m", "4/mmm", "-3", "-3m", "6/m", "6/mmm", "m-3", "m-3m"].contains(&point_group)
}

pub fn reciprocal_path(input: &StructureInput, spacing: f64, symprec: f64, time_reversal: bool) -> Result<ToolOutput, String> {
    positive(spacing, "spacing_inv_A")?;
    positive(symprec, "symprec_A")?;
    let structure = &input.structure;
    if !structure.has_unit_cell || !(input.pbc[0] && input.pbc[1] && input.pbc[2]) || structure.atoms.is_empty() {
        return Err("Reciprocal paths require a nonempty periodic structure".into());
    }
    let cell = structure.cell_vectors;
    cell_matrix(NdArray { shape: vec![3, 3], values: cell.iter().flatten().copied().collect() })?;
    let inv = inverse(&cell)?;
    let fractions: Vec<Vec3> = structure.atoms.iter().map(|a| row_times([a.x, a.y, a.z], &inv)).collect();
    let types: Vec<i32> = structure.atoms.iter().map(|a| a.atomic_number).collect();
    let d = af_spglib::get_dataset(&cell, &fractions, &types, symprec)
        .ok_or_else(|| "Spglib could not detect the symmetry of the system".to_string())?;
    let mut conventional = d.std_lattice;
    let mut conventional_positions: Vec<Vec3> = d.std_positions.clone();
    let conventional_types = d.std_types.clone();
    let number = d.spacegroup_number;
    let international = d.international_symbol.clone();
    let point_group = d.pointgroup_symbol.clone();

    let family = match number {
        ..=2 => 'a',
        ..=15 => 'm',
        ..=74 => 'o',
        ..=142 => 't',
        ..=194 => 'h',
        _ => 'c',
    };
    let bravais = format!("{family}{}", international.chars().next().map(String::from).unwrap_or_default());
    let has_inversion = centrosymmetric(&point_group);
    let mut cp = parameters(&conventional);
    let extended: String = match bravais.as_str() {
        "cP" => if number <= 206 { "cP1" } else { "cP2" }.into(),
        "cF" => if number <= 206 { "cF1" } else { "cF2" }.into(),
        "cI" => "cI1".into(),
        "tP" => "tP1".into(),
        "tI" => if cp.c <= cp.a { "tI1" } else { "tI2" }.into(),
        "oP" => "oP1".into(),
        "oF" => {
            let (a, b, c) = (cp.a, cp.b, cp.c);
            if 1.0 / (a * a) > 1.0 / (b * b) + 1.0 / (c * c) {
                "oF1"
            } else if 1.0 / (c * c) > 1.0 / (a * a) + 1.0 / (b * b) {
                "oF2"
            } else {
                "oF3"
            }
            .into()
        }
        "oI" => {
            // The longest conventional vector selects oI1 (c), oI2 (a) or oI3 (b).
            let mut largest = (cp.c, 1);
            for candidate in [(cp.b, 3), (cp.a, 2)] {
                if largest.0 < candidate.0 || (largest.0 == candidate.0 && largest.1 < candidate.1) {
                    largest = candidate;
                }
            }
            format!("oI{}", largest.1)
        }
        "oC" => if cp.a <= cp.b { "oC1" } else { "oC2" }.into(),
        "oA" => if cp.b <= cp.c { "oA1" } else { "oA2" }.into(),
        "hP" => {
            const HP1: [i32; 15] = [143, 144, 145, 146, 147, 148, 149, 151, 153, 157, 159, 160, 161, 162, 163];
            if HP1.contains(&number) { "hP1" } else { "hP2" }.into()
        }
        "hR" => if 3f64.sqrt() * cp.a <= 2f64.sqrt() * cp.c { "hR1" } else { "hR2" }.into(),
        "mP" => "mP1".into(),
        "mC" => {
            let (a, b, c) = (cp.a, cp.b, cp.c);
            let sin_beta2 = 1.0 - cp.cos_beta * cp.cos_beta;
            if b < a * sin_beta2.sqrt() {
                "mC1"
            } else if -a * cp.cos_beta / c + a * a * sin_beta2 / (b * b) <= 1.0 {
                "mC2"
            } else {
                "mC3"
            }
            .into()
        }
        "aP" => {
            // Reciprocal-space Niggli reduction, then the HPKOT all-acute/all-obtuse choice.
            let reciprocal = af_spglib::niggli_reduce(&reciprocal_rows(&conventional)?, 1e-5)
                .ok_or_else(|| "Niggli reduction failed for the triclinic cell".to_string())?;
            let real2 = reciprocal_rows(&reciprocal)?;
            let k2 = parameters(&reciprocal);
            let conditions = [(k2.b * k2.c * k2.cos_alpha).abs(), (k2.c * k2.a * k2.cos_beta).abs(), (k2.a * k2.b * k2.cos_gamma).abs()];
            let m2_choices = [integer_matrix([0, 0, 1, 1, 0, 0, 0, 1, 0]), integer_matrix([0, 1, 0, 0, 0, 1, 1, 0, 0]), identity()];
            let mut smallest = 0;
            for k in 1..3 {
                if conditions[k] < conditions[smallest] {
                    smallest = k;
                }
            }
            let real3 = multiply(&transpose(&m2_choices[smallest]), &real2);
            let k3 = parameters(&reciprocal_rows(&real3)?);
            let (pa, pb, pg) = (k3.cos_alpha > 0.0, k3.cos_beta > 0.0, k3.cos_gamma > 0.0);
            let m3 = if (pa && !pb && !pg) || (!pa && pb && pg) {
                integer_matrix([1, 0, 0, 0, -1, 0, 0, 0, -1])
            } else if (!pa && pb && !pg) || (pa && !pb && pg) {
                integer_matrix([-1, 0, 0, 0, 1, 0, 0, 0, -1])
            } else if (!pa && !pb && pg) || (pa && pb && !pg) {
                integer_matrix([-1, 0, 0, 0, -1, 0, 0, 0, 1])
            } else {
                identity()
            };
            let last = multiply(&transpose(&m3), &real3);
            let kf = parameters(&reciprocal_rows(&last)?);
            let name = if kf.cos_alpha <= 0.0 && kf.cos_beta <= 0.0 && kf.cos_gamma <= 0.0 {
                "aP2"
            } else if kf.cos_alpha >= 0.0 && kf.cos_beta >= 0.0 && kf.cos_gamma >= 0.0 {
                "aP3"
            } else {
                return Err("Unexpected aP triclinic lattice: neither all-obtuse nor all-acute".into());
            };
            let last_inverse = inverse(&last)?;
            for position in conventional_positions.iter_mut() {
                *position = row_times(row_times(*position, &conventional), &last_inverse);
            }
            conventional = last;
            cp = parameters(&last);
            name.into()
        }
        _ => return Err(format!("Unknown Bravais lattice {bravais} for space group {number}")),
    };
    let CellParameters { a, b, c, cos_alpha, cos_beta, cos_gamma } = cp;

    // Primitive cell in the HPKOT convention (positions not wrapped).
    let (p, inverse_p) = primitive_matrices(&bravais);
    let primitive = multiply(&transpose(&p), &conventional);
    let inverse_pt = transpose(&inverse_p);
    let mapped: Vec<Vec3> = conventional_positions.iter().map(|x| row_times(*x, &inverse_pt)).collect();
    let ratio = determinant(&inverse_p).round() as i32;
    let mut group = vec![-1i32; mapped.len()];
    let mut representatives = Vec::new();
    for i in 0..mapped.len() {
        if group[i] >= 0 {
            continue;
        }
        let mut members = 0;
        for j in i..mapped.len() {
            let mut same = true;
            for k in 0..3 {
                let d = mapped[j][k] - mapped[i][k];
                same = same && (d - d.round()).abs() < 1e-6;
            }
            if !same {
                continue;
            }
            if group[j] >= 0 || conventional_types[j] != conventional_types[i] {
                return Err("Problem creating the primitive cell: overlapping atoms of different type".into());
            }
            group[j] = representatives.len() as i32;
            members += 1;
        }
        if members != ratio {
            return Err("Problem creating the primitive cell from the conventional cell".into());
        }
        representatives.push(i);
    }
    let mut primitive_structure = Structure { has_unit_cell: true, cell_vectors: primitive, ..Default::default() };
    for &index in &representatives {
        let z = conventional_types[index];
        let position = row_times(mapped[index], &primitive);
        let [r, g, bl] = elements::default_colour(z);
        primitive_structure.atoms.push(AtomSite {
            atomic_number: z,
            symbol: elements::symbol(z).to_string(),
            x: position[0],
            y: position[1],
            z: position[2],
            r,
            g,
            b: bl,
            ..Default::default()
        });
    }

    let table = LATTICES.iter().find(|l| l.name == extended).ok_or_else(|| format!("Missing HPKOT data for {extended}"))?;
    let mut variables: BTreeMap<String, f64> = BTreeMap::from([
        ("a".to_string(), a),
        ("b".to_string(), b),
        ("c".to_string(), c),
        ("cosalpha".to_string(), cos_alpha),
        ("cosbeta".to_string(), cos_beta),
        ("cosgamma".to_string(), cos_gamma),
        ("sinbeta".to_string(), (1.0 - cos_beta * cos_beta).sqrt()),
    ]);
    for (name, expression) in table.parameters {
        let value = Expression::evaluate(expression, &variables)?;
        variables.insert(name.to_string(), value);
    }
    let mut points: Vec<(String, Vec3)> = Vec::new();
    for (label, x, y, z) in table.points {
        points.push((
            label.to_string(),
            [Expression::evaluate(x, &variables)?, Expression::evaluate(y, &variables)?, Expression::evaluate(z, &variables)?],
        ));
    }
    let mut path: Vec<(String, String)> = table.path.iter().map(|(f, t)| (f.to_string(), t.to_string())).collect();
    // Without inversion or time reversal, k and -k are inequivalent: add the -k path.
    let augmented = !has_inversion && !time_reversal;
    if augmented {
        let original = points.clone();
        for (label, coordinates) in original {
            if label != "GAMMA" {
                points.push((format!("{label}'"), scale(coordinates, -1.0)));
            }
        }
        let original_path = path.clone();
        let prime = |s: &String| if s == "GAMMA" { s.clone() } else { format!("{s}'") };
        for (from, to) in &original_path {
            path.push((prime(from), prime(to)));
        }
    }
    let coordinate = |label: &str| -> Result<Vec3, String> {
        points.iter().find(|p| p.0 == label).map(|p| p.1).ok_or_else(|| format!("Undefined k-point {label}"))
    };
    let reciprocal = reciprocal_rows(&primitive)?;
    let mut k_fractional: Vec<Vec3> = Vec::new();
    let mut labels: Vec<String> = Vec::new();
    let mut distance = Vec::new();
    let mut segments = Vec::new();
    let mut previous = 0.0;
    for (from, to) in &path {
        let (start, stop) = (coordinate(from)?, coordinate(to)?);
        let length = norm(sub(row_times(stop, &reciprocal), row_times(start, &reciprocal)));
        let samples = 2.max((length / spacing) as i32);
        let mut segment_start = labels.len() as i64;
        for i in 0..samples {
            if i == 0 && labels.last().is_some_and(|l| l == from) {
                segment_start -= 1;
                continue;
            }
            let fraction = i as f64 / (samples - 1) as f64;
            k_fractional.push(add(start, scale(sub(stop, start), fraction)));
            labels.push(if i == 0 {
                from.clone()
            } else if i == samples - 1 {
                to.clone()
            } else {
                String::new()
            });
            distance.push(previous + length * fraction);
        }
        previous += length;
        segments.push(Json::Array(vec![Json::from(segment_start), Json::from(labels.len() as i64)]));
    }
    let mut output = ToolOutput::default();
    let k_fractional_json: Vec<Json> = k_fractional.iter().map(|k| jv(k)).collect();
    let k_cartesian: Vec<Json> = k_fractional.iter().map(|k| jv(&row_times(*k, &reciprocal))).collect();
    let mut special = Json::object();
    for (label, coordinates) in &points {
        special.set(label, jv(coordinates));
    }
    let mut result = Json::object();
    result.set("primitive_structure", structure_json(&primitive_structure));
    output.frames.push(primitive_structure);
    result.set("spacegroup_number", number);
    result.set("spacegroup_symbol", international);
    result.set("bravais_lattice_extended", extended);
    result.set("has_inversion_symmetry", has_inversion);
    result.set("augmented_path", augmented);
    result.set("path", Json::Array(path.iter().map(|(f, t)| Json::Array(vec![Json::from(f.as_str()), Json::from(t.as_str())])).collect()));
    result.set("special_points_fractional", special);
    result.set("kpoints_fractional", Json::Array(k_fractional_json));
    result.set("kpoints_inv_A", Json::Array(k_cartesian));
    result.set("distance_inv_A", jv(&distance));
    result.set("labels", Json::Array(labels.into_iter().map(Json::from).collect()));
    result.set("segments", Json::Array(segments));
    result.set("reciprocal_lattice_inv_A", mat3_json(&reciprocal));
    output.result = result;
    Ok(output)
}

/// VASP KPOINTS/POSCAR and Quantum ESPRESSO cards for the HPKOT path.
pub fn dft_inputs(input: &StructureInput, points_per_segment: i64, symprec: f64, time_reversal: bool) -> Result<ToolOutput, String> {
    if !(2..=1000).contains(&points_per_segment) {
        return Err("points_per_segment must be between 2 and 1000".into());
    }
    let mut path = reciprocal_path(input, 0.05, symprec, time_reversal)?;
    let primitive = path.frames[0].clone();
    let special = path.result.at("special_points_fractional")?.clone();
    let mut segments: Vec<(String, String)> = Vec::new();
    for segment in path.result.at("path")?.items()? {
        let pair = segment.items()?;
        segments.push((pair[0].string()?.to_string(), pair[1].string()?.to_string()));
    }
    let coordinate = |label: &str| -> Result<Vec<f64>, String> { special.at(label)?.items()?.iter().map(|v| Ok(v.number()?)).collect() };
    let qe_label = |label: &str| if label == "GAMMA" { "G".to_string() } else { label.to_string() };
    let (mut kpoints, mut poscar, mut qe) = (String::new(), String::new(), String::new());
    let title = format!(
        "HPKOT path, space group {} ({})",
        path.result.at("spacegroup_symbol")?.string()?,
        path.result.at("spacegroup_number")?.number()? as i32
    );
    let _ = write!(kpoints, "{title}\n{points_per_segment}\nLine-mode\nReciprocal\n");
    for (s, (from, to)) in segments.iter().enumerate() {
        if s > 0 {
            kpoints.push('\n');
        }
        for label in [from, to] {
            let k = coordinate(label)?;
            let shown = if label == "GAMMA" { "\\Gamma".to_string() } else { label.clone() };
            let _ = writeln!(kpoints, "{:.10} {:.10} {:.10} ! {shown}", k[0], k[1], k[2]);
        }
    }
    // POSCAR (VASP 5), species grouped in order of first appearance.
    let mut species: Vec<String> = Vec::new();
    for atom in &primitive.atoms {
        if !species.contains(&atom.symbol) {
            species.push(atom.symbol.clone());
        }
    }
    let _ = writeln!(poscar, "{title}\n1.0");
    for row in &primitive.cell_vectors {
        let _ = writeln!(poscar, "  {:.12} {:.12} {:.12}", row[0], row[1], row[2]);
    }
    for s in &species {
        let _ = write!(poscar, " {s}");
    }
    poscar.push('\n');
    for s in &species {
        let _ = write!(poscar, " {}", primitive.atoms.iter().filter(|a| &a.symbol == s).count());
    }
    poscar.push_str("\nDirect\n");
    let cell = primitive.cell_vectors;
    for s in &species {
        for atom in primitive.atoms.iter().filter(|a| &a.symbol == s) {
            let f = fractional([atom.x, atom.y, atom.z], &cell)?;
            let _ = writeln!(poscar, "  {:.12} {:.12} {:.12}", f[0], f[1], f[2]);
        }
    }
    // Quantum ESPRESSO: a point weight is the number of steps to the next point; 1 jumps.
    let mut sequence: Vec<(String, i64)> = Vec::new();
    for (from, to) in &segments {
        if sequence.last().is_none_or(|last| &last.0 != from) {
            if let Some(last) = sequence.last_mut() {
                last.1 = 1;
            }
            sequence.push((from.clone(), points_per_segment));
        } else if let Some(last) = sequence.last_mut() {
            last.1 = points_per_segment;
        }
        sequence.push((to.clone(), 1));
    }
    qe.push_str("CELL_PARAMETERS angstrom\n");
    for row in &primitive.cell_vectors {
        let _ = writeln!(qe, "  {:.10} {:.10} {:.10}", row[0], row[1], row[2]);
    }
    qe.push_str("ATOMIC_SPECIES\n");
    for s in &species {
        let _ = writeln!(qe, "  {s} {:.10} {s}.UPF", atomic_mass(s)?);
    }
    qe.push_str("ATOMIC_POSITIONS crystal\n");
    for atom in &primitive.atoms {
        let f = fractional([atom.x, atom.y, atom.z], &cell)?;
        let _ = writeln!(qe, "  {} {:.10} {:.10} {:.10}", atom.symbol, f[0], f[1], f[2]);
    }
    let _ = writeln!(qe, "K_POINTS crystal_b\n{}", sequence.len());
    for (label, weight) in &sequence {
        let k = coordinate(label)?;
        let _ = writeln!(qe, "  {:.10} {:.10} {:.10} {weight} ! {}", k[0], k[1], k[2], qe_label(label));
    }
    let mut files = Json::object();
    files.set("KPOINTS", kpoints);
    files.set("POSCAR", poscar);
    files.set("qe_band_cards.in", qe);
    path.result.set("points_per_segment", points_per_segment);
    path.result.set("files", files);
    Ok(path)
}
