//! Native interatomic potentials (the C++ science/Potentials and
//! ManyBodyPotentials): EMT, Lennard-Jones, EAM (LAMMPS tables), Tersoff,
//! Stillinger-Weber and Buckingham with Ewald-summed Coulomb interactions.
//! Energies in eV, forces in eV/Angstrom, stress in eV/Angstrom^3 (tension
//! positive, ASE's sign).

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use af_core::elements;
use af_core::format::stod;
use af_core::json::Json;
use af_core::structure::{AtomSite, Structure};

use crate::core::*;
use crate::data::{atomic_mass, display, join_normal};

extern "C" {
    #[link_name = "erfc"]
    fn c_erfc(x: f64) -> f64;
}

fn erfc(x: f64) -> f64 {
    unsafe { c_erfc(x) }
}

/// Atoms for force evaluation. Units: Angstrom, amu.
#[derive(Clone, Debug, Default)]
pub struct Configuration {
    pub symbols: Vec<String>,
    pub numbers: Vec<i32>,
    pub positions: Vec<Vec3>,
    pub masses: Vec<f64>,
    pub cell: Mat3,
    pub pbc: Pbc,
}

impl Configuration {
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    pub fn fully_periodic(&self) -> bool {
        self.pbc[0] && self.pbc[1] && self.pbc[2]
    }
}

pub fn configuration_from(structure: &Structure, pbc: &Pbc) -> Result<Configuration, String> {
    let mut result = Configuration::default();
    for atom in &structure.atoms {
        result.symbols.push(atom.symbol.clone());
        result.numbers.push(if atom.atomic_number != 0 { atom.atomic_number } else { elements::atomic_number(&atom.symbol) });
        result.positions.push([atom.x, atom.y, atom.z]);
        result.masses.push(atomic_mass(&atom.symbol)?);
    }
    if structure.has_unit_cell {
        result.cell = structure.cell_vectors;
    }
    result.pbc = if structure.has_unit_cell { *pbc } else { [false; 3] };
    Ok(result)
}

pub fn structure_from(configuration: &Configuration) -> Structure {
    let mut structure = Structure::default();
    for i in 0..configuration.len() {
        let z = configuration.numbers[i];
        let [r, g, b] = elements::default_colour(z);
        let p = configuration.positions[i];
        structure.atoms.push(AtomSite {
            symbol: configuration.symbols[i].clone(),
            atomic_number: z,
            x: p[0],
            y: p[1],
            z: p[2],
            r,
            g,
            b,
            ..Default::default()
        });
    }
    if determinant(&configuration.cell).abs() > 1e-12 {
        structure.has_unit_cell = true;
        structure.cell_vectors = configuration.cell;
    }
    structure
}

#[derive(Clone, Debug, Default)]
pub struct PotentialResult {
    /// eV.
    pub energy: f64,
    /// eV/Angstrom.
    pub forces: Vec<Vec3>,
    /// eV/Angstrom^3, tension positive (ASE sign convention).
    pub stress: Mat3,
    pub has_stress: bool,
}

/// A native interatomic potential; stateless between calls.
pub trait Potential: Send + Sync {
    fn compute(&self, configuration: &Configuration, stress: bool) -> Result<PotentialResult, String>;
    fn description(&self) -> String;
}

/// Potentials are stateless, so one parsed instance serves every image or run.
pub type PotentialHandle = Arc<dyn Potential>;

fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

fn min(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

/// Each unordered pair once, including periodic images (half neighbour list).
fn half_list(c: &Configuration, cutoff: f64) -> Result<Vec<Neighbor>, String> {
    let periodic = c.pbc.iter().any(|&p| p);
    let cell = if periodic { c.cell } else { identity() };
    let mut result = Vec::new();
    for neighbor in neighbor_list(&c.positions, &cell, c.pbc, cutoff)? {
        let mut keep = neighbor.i < neighbor.j;
        if neighbor.i == neighbor.j {
            for k in 0..3 {
                if neighbor.shift[k] != 0 {
                    keep = neighbor.shift[k] > 0;
                    break;
                }
            }
        }
        if keep {
            result.push(neighbor);
        }
    }
    Ok(result)
}

fn add_virial(virial: &mut Mat3, force_on_i: Vec3, d: Vec3) {
    for a in 0..3 {
        for b in 0..3 {
            virial[a][b] += 0.5 * (force_on_i[a] * d[b] + force_on_i[b] * d[a]);
        }
    }
}

fn finish_stress(result: &mut PotentialResult, c: &Configuration, virial: &Mat3, stress: bool) -> Result<(), String> {
    if !stress {
        return Ok(());
    }
    if !c.fully_periodic() {
        return Err("Stress requires a full periodic cell".into());
    }
    let volume = cell_volume(&c.cell);
    for a in 0..3 {
        for b in 0..3 {
            result.stress[a][b] = virial[a][b] / volume;
        }
    }
    result.has_stress = true;
    Ok(())
}

// ---------------------------------------------------------------- EMT

const BOHR: f64 = 0.52917721067;
/// ASE's historical rounding of (16 pi / 3)^(1/3) / sqrt(2).
const BETA: f64 = 1.809;

/// E0 (eV), s0 (bohr), V0 (eV), eta2, kappa, lambda (1/bohr), n0 (1/bohr^3).
const EMT_PARAMETERS: [(&str, [f64; 7]); 11] = [
    ("Al", [-3.28, 3.00, 1.493, 1.240, 2.000, 1.169, 0.00700]),
    ("Cu", [-3.51, 2.67, 2.476, 1.652, 2.740, 1.906, 0.00910]),
    ("Ag", [-2.96, 3.01, 2.132, 1.652, 2.790, 1.892, 0.00547]),
    ("Au", [-3.80, 3.00, 2.321, 1.674, 2.873, 2.182, 0.00703]),
    ("Ni", [-4.44, 2.60, 3.673, 1.669, 2.757, 1.948, 0.01030]),
    ("Pd", [-3.90, 2.87, 2.773, 1.818, 3.107, 2.155, 0.00688]),
    ("Pt", [-5.85, 2.90, 4.067, 1.812, 3.145, 2.192, 0.00802]),
    ("H", [-3.21, 1.31, 0.132, 2.652, 2.790, 3.892, 0.00547]),
    ("C", [-3.50, 1.81, 0.332, 1.652, 2.790, 1.892, 0.01322]),
    ("N", [-5.10, 1.88, 0.132, 1.652, 2.790, 1.892, 0.01222]),
    ("O", [-4.60, 1.95, 0.332, 1.652, 2.790, 1.892, 0.00850]),
];

#[derive(Clone, Copy)]
struct EmtElement {
    e0: f64,
    s0: f64,
    v0: f64,
    eta2: f64,
    kappa: f64,
    lambda: f64,
    n0: f64,
    gamma1: f64,
    gamma2: f64,
}

pub struct Emt {
    rc: f64,
    acut: f64,
    list_cutoff: f64,
    elements: BTreeMap<&'static str, EmtElement>,
}

impl Emt {
    pub fn new() -> Self {
        let mut max_s0 = 0.0;
        for (_, raw) in &EMT_PARAMETERS {
            max_s0 = max(max_s0, raw[1] * BOHR);
        }
        let rc = BETA * max_s0 * 0.5 * (3f64.sqrt() + 4f64.sqrt());
        let rr = rc * 2.0 * 4f64.sqrt() / (3f64.sqrt() + 4f64.sqrt());
        let acut = 9999f64.ln() / (rr - rc);
        let mut elements = BTreeMap::new();
        for (symbol, raw) in &EMT_PARAMETERS {
            let mut p = EmtElement {
                e0: raw[0],
                s0: raw[1] * BOHR,
                v0: raw[2],
                eta2: raw[3] / BOHR,
                kappa: raw[4] / BOHR,
                lambda: raw[5] / BOHR,
                n0: raw[6] / (BOHR * BOHR * BOHR),
                gamma1: 0.0,
                gamma2: 0.0,
            };
            let shells = [12.0, 6.0, 24.0];
            for i in 0..3 {
                let r = p.s0 * BETA * (i as f64 + 1.0).sqrt();
                let x = shells[i] / (12.0 * (1.0 + (acut * (r - rc)).exp()));
                p.gamma1 += x * (-p.eta2 * (r - BETA * p.s0)).exp();
                p.gamma2 += x * (-p.kappa / BETA * (r - BETA * p.s0)).exp();
            }
            elements.insert(*symbol, p);
        }
        Self { rc, acut, list_cutoff: rc + 0.5, elements }
    }
}

impl Default for Emt {
    fn default() -> Self {
        Self::new()
    }
}

impl Potential for Emt {
    fn description(&self) -> String {
        "EMT (effective-medium theory)".into()
    }

    fn compute(&self, c: &Configuration, stress: bool) -> Result<PotentialResult, String> {
        let n = c.len();
        let mut parameters = Vec::with_capacity(n);
        for symbol in &c.symbols {
            let Some(p) = self.elements.get(symbol.as_str()) else {
                return Err(format!("No EMT potential for {symbol}; EMT supports Al, Cu, Ag, Au, Ni, Pd, Pt, H, C, N and O"));
            };
            parameters.push(*p);
        }
        let pairs = half_list(c, self.list_cutoff)?;
        let mut result = PotentialResult { forces: vec![[0.0; 3]; n], ..Default::default() };
        let mut virial = [[0.0; 3]; 3];
        let mut sigma1 = vec![0.0; n];
        let mut deds = vec![0.0; n];
        for pair in &pairs {
            let (a1, a2) = (pair.i, pair.j);
            let (p1, p2) = (&parameters[a1], &parameters[a2]);
            let ksi = p2.n0 / p1.n0;
            let d = pair.vector;
            let r = norm(d);
            let x = (self.acut * (r - self.rc)).exp();
            let theta = 1.0 / (1.0 + x);
            let y1 = 0.5 * p1.v0 * (-p2.kappa * (r / BETA - p2.s0)).exp() * ksi / p1.gamma2 * theta;
            let y2 = 0.5 * p2.v0 * (-p1.kappa * (r / BETA - p1.s0)).exp() / ksi / p2.gamma2 * theta;
            result.energy -= y1 + y2;
            let f = scale(d, ((y1 * p2.kappa + y2 * p1.kappa) / BETA + (y1 + y2) * self.acut * theta * x) / r);
            result.forces[a1] = add(result.forces[a1], f);
            result.forces[a2] = sub(result.forces[a2], f);
            add_virial(&mut virial, f, d);
            sigma1[a1] += (-p2.eta2 * (r - BETA * p2.s0)).exp() * ksi * theta / p1.gamma1;
            sigma1[a2] += (-p1.eta2 * (r - BETA * p1.s0)).exp() / ksi * theta / p2.gamma1;
        }
        for a in 0..n {
            let p = &parameters[a];
            if !(sigma1[a] > 0.0) {
                deds[a] = 0.0;
                result.energy -= p.e0;
                continue;
            }
            let ds = -(sigma1[a] / 12.0).ln() / (BETA * p.eta2);
            let x = p.lambda * ds;
            let y = (-x).exp();
            let z = 6.0 * p.v0 * (-p.kappa * ds).exp();
            deds[a] = (x * y * p.e0 * p.lambda + p.kappa * z) / (sigma1[a] * BETA * p.eta2);
            result.energy += p.e0 * ((1.0 + x) * y - 1.0) + z;
        }
        for pair in &pairs {
            let (a1, a2) = (pair.i, pair.j);
            let (p1, p2) = (&parameters[a1], &parameters[a2]);
            let ksi = p2.n0 / p1.n0;
            let d = pair.vector;
            let r = norm(d);
            let x = (self.acut * (r - self.rc)).exp();
            let theta = 1.0 / (1.0 + x);
            let y1 = (-p2.eta2 * (r - BETA * p2.s0)).exp() * ksi / p1.gamma1 * theta * deds[a1];
            let y2 = (-p1.eta2 * (r - BETA * p1.s0)).exp() / ksi / p2.gamma1 * theta * deds[a2];
            let f = scale(d, ((y1 * p2.eta2 + y2 * p1.eta2) + (y1 + y2) * self.acut * theta * x) / r);
            result.forces[a1] = sub(result.forces[a1], f);
            result.forces[a2] = add(result.forces[a2], f);
            add_virial(&mut virial, scale(f, -1.0), d);
        }
        finish_stress(&mut result, c, &virial, stress)?;
        Ok(result)
    }
}

// --------------------------------------------------------- Lennard-Jones

pub struct LennardJones {
    epsilon: f64,
    sigma: f64,
    cutoff: f64,
    shift: f64,
}

impl LennardJones {
    pub fn new(epsilon: f64, sigma: f64, cutoff: f64) -> Result<Self, String> {
        positive(epsilon, "epsilon")?;
        positive(sigma, "sigma")?;
        positive(cutoff, "cutoff")?;
        let ratio6 = (sigma / cutoff).powf(6.0);
        Ok(Self { epsilon, sigma, cutoff, shift: 4.0 * epsilon * (ratio6 * ratio6 - ratio6) })
    }
}

impl Potential for LennardJones {
    fn description(&self) -> String {
        "Lennard-Jones".into()
    }

    fn compute(&self, c: &Configuration, stress: bool) -> Result<PotentialResult, String> {
        let mut result = PotentialResult { forces: vec![[0.0; 3]; c.len()], ..Default::default() };
        let mut virial = [[0.0; 3]; 3];
        for pair in half_list(c, self.cutoff)? {
            let r2 = dot(pair.vector, pair.vector);
            let c6 = (self.sigma * self.sigma / r2).powf(3.0);
            let c12 = c6 * c6;
            result.energy += 4.0 * self.epsilon * (c12 - c6) - self.shift;
            // Force on j along d = r_j - r_i; the reaction acts on i.
            let on_j = scale(pair.vector, 24.0 * self.epsilon * (2.0 * c12 - c6) / r2);
            result.forces[pair.j] = add(result.forces[pair.j], on_j);
            result.forces[pair.i] = sub(result.forces[pair.i], on_j);
            add_virial(&mut virial, scale(on_j, -1.0), pair.vector);
        }
        finish_stress(&mut result, c, &virial, stress)?;
        Ok(result)
    }
}

// ---------------------------------------------------------------- EAM

/// LAMMPS pair_eam cubic interpolation of a uniformly tabulated function
/// (PairEAM::interpolate), with zero-based indices.
#[derive(Clone, Default)]
struct Spline {
    delta: f64,
    c: Vec<[f64; 7]>,
}

impl Spline {
    fn build(values: &[f64], delta: f64) -> Result<Self, String> {
        let n = values.len();
        if n < 5 || !(delta > 0.0) {
            return Err("EAM tables need at least five samples and a positive spacing".into());
        }
        let mut c = vec![[0.0; 7]; n];
        for m in 0..n {
            c[m][6] = values[m];
        }
        c[0][5] = c[1][6] - c[0][6];
        c[1][5] = 0.5 * (c[2][6] - c[0][6]);
        c[n - 2][5] = 0.5 * (c[n - 1][6] - c[n - 3][6]);
        c[n - 1][5] = c[n - 1][6] - c[n - 2][6];
        for m in 2..n - 2 {
            c[m][5] = ((c[m - 2][6] - c[m + 2][6]) + 8.0 * (c[m + 1][6] - c[m - 1][6])) / 12.0;
        }
        for m in 0..n - 1 {
            c[m][4] = 3.0 * (c[m + 1][6] - c[m][6]) - 2.0 * c[m][5] - c[m + 1][5];
            c[m][3] = c[m][5] + c[m + 1][5] - 2.0 * (c[m + 1][6] - c[m][6]);
        }
        c[n - 1][4] = 0.0;
        c[n - 1][3] = 0.0;
        for m in 0..n {
            c[m][2] = c[m][5] / delta;
            c[m][1] = 2.0 * c[m][4] / delta;
            c[m][0] = 3.0 * c[m][3] / delta;
        }
        Ok(Self { delta, c })
    }

    /// (value, derivative).
    fn evaluate(&self, x: f64) -> (f64, f64) {
        let mut p = max(0.0, x) / self.delta;
        let mut m = p as i32;
        m = m.min(self.c.len() as i32 - 2);
        p -= m as f64;
        p = min(p, 1.0);
        let k = &self.c[m as usize];
        (((k[3] * p + k[4]) * p + k[5]) * p + k[6], (k[0] * p + k[1]) * p + k[2])
    }
}

/// Reads a text file as std::ifstream in text mode: lines without "\r\n".
fn text_lines(path: &Path) -> Option<Vec<String>> {
    let bytes = std::fs::read(path).ok()?;
    let text = String::from_utf8_lossy(&bytes).replace("\r\n", "\n");
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    if text.ends_with('\n') {
        lines.pop();
    }
    Some(lines)
}

fn words(text: &str) -> impl Iterator<Item = &str> {
    text.split(|c: char| matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')).filter(|t| !t.is_empty())
}

/// istream >> int: the leading integer of a word, 0 when there is none.
fn leading_int(word: Option<&str>) -> i64 {
    let Some(word) = word else { return 0 };
    let bytes = word.as_bytes();
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    let start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == start {
        return 0;
    }
    word[..i].parse::<i64>().unwrap_or(if bytes[0] == b'-' { i64::from(i32::MIN) } else { i64::from(i32::MAX) }).clamp(i64::from(i32::MIN), i64::from(i32::MAX))
}

fn file_name_of(path: &Path) -> String {
    path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()
}

pub struct Eam {
    elements: Vec<String>,
    embedding: Vec<Spline>,
    density: Vec<Vec<Spline>>,
    pair: Vec<Spline>,
    cutoff: f64,
    rho_max: f64,
    finnis_sinclair: bool,
    description: String,
}

impl Eam {
    pub fn new(file: &Path, format: &str) -> Result<Self, String> {
        let lines = text_lines(file).ok_or_else(|| format!("Cannot open EAM potential file {}", display(file)))?;
        let mut format = format.to_string();
        if format == "auto" {
            let lower = file_name_of(file).to_ascii_lowercase();
            format = if lower.ends_with(".fs") || lower.ends_with(".eam.fs") {
                "fs"
            } else if lower.ends_with(".alloy") || lower.ends_with(".setfl") {
                "setfl"
            } else {
                "funcfl"
            }
            .to_string();
        }
        if format != "setfl" && format != "fs" && format != "funcfl" {
            return Err("EAM format must be setfl, fs or funcfl".into());
        }
        let line = |k: usize| lines.get(k).map(String::as_str).unwrap_or("");
        let truncated = || format!("Truncated EAM potential file {}", display(file));
        let mut eam = Eam {
            elements: Vec::new(),
            embedding: Vec::new(),
            density: Vec::new(),
            pair: Vec::new(),
            cutoff: 0.0,
            rho_max: 0.0,
            finnis_sinclair: false,
            description: String::new(),
        };
        let skip = if format == "funcfl" { 1 } else { 4 };
        let mut t = Tokens { tokens: lines.iter().skip(skip).flat_map(|l| words(l)).map(str::to_string).collect(), cursor: 0, truncated: truncated() };
        if format == "funcfl" {
            let z = t.number()? as i32;
            t.number()?; // mass
            t.number()?; // lattice constant
            t.skip(); // lattice name
            let nrho = t.number()? as i32;
            let drho = t.number()?;
            let nr = t.number()? as i32;
            let dr = t.number()?;
            eam.cutoff = t.number()?;
            eam.elements = vec![elements::symbol(z).to_string()];
            eam.embedding = vec![Spline::build(&t.table(nrho)?, drho)?];
            eam.rho_max = (nrho - 1) as f64 * drho;
            let charge = t.table(nr)?;
            // funcfl Z(r) in e: r*phi = 27.2 * 0.529 * Z(r)^2 (Hartree-Bohr in eV-Angstrom), as LAMMPS.
            let rphi: Vec<f64> = charge.iter().map(|q| 27.2 * 0.529 * q * q).collect();
            eam.density = vec![vec![Spline::build(&t.table(nr)?, dr)?]];
            eam.pair = vec![Spline::build(&rphi, dr)?];
        } else {
            let mut header = words(line(3));
            let count = leading_int(header.next());
            if !(1..=20).contains(&count) {
                return Err(format!("Invalid element count in EAM file {}", display(file)));
            }
            let count = count as usize;
            eam.elements = (0..count).map(|_| header.next().unwrap_or("").to_string()).collect();
            let nrho = t.number()? as i32;
            let drho = t.number()?;
            let nr = t.number()? as i32;
            let dr = t.number()?;
            eam.cutoff = t.number()?;
            eam.rho_max = (nrho - 1) as f64 * drho;
            let per_element = if format == "fs" { count } else { 1 };
            for _ in 0..count {
                t.number()?; // atomic number
                t.number()?; // mass
                t.number()?; // lattice constant
                t.skip(); // lattice name
                eam.embedding.push(Spline::build(&t.table(nrho)?, drho)?);
                let mut densities = Vec::new();
                for _ in 0..per_element {
                    densities.push(Spline::build(&t.table(nr)?, dr)?);
                }
                eam.density.push(densities);
            }
            for i in 0..count {
                for _ in 0..=i {
                    eam.pair.push(Spline::build(&t.table(nr)?, dr)?);
                }
            }
        }
        eam.finnis_sinclair = format == "fs";
        if !(eam.cutoff > 0.0) {
            return Err("EAM cutoff must be positive".into());
        }
        eam.description = format!("EAM ({format}, {})", file_name_of(file));
        Ok(eam)
    }

    /// Density at i contributed by a neighbour of type `from`: rho_from (setfl) or rho_{from,at} (fs).
    fn density(&self, from: usize, at: usize) -> &Spline {
        if self.finnis_sinclair {
            &self.density[from][at]
        } else {
            &self.density[from][0]
        }
    }
}

/// The whitespace-separated words of an EAM file, read in order.
struct Tokens {
    tokens: Vec<String>,
    cursor: usize,
    truncated: String,
}

impl Tokens {
    fn number(&mut self) -> Result<f64, String> {
        let token = self.tokens.get(self.cursor).ok_or_else(|| self.truncated.clone())?;
        self.cursor += 1;
        stod(token).map(|(v, _)| v).map_err(|_| "stod".to_string())
    }

    /// `++cursor`: a word skipped without parsing.
    fn skip(&mut self) {
        self.cursor += 1;
    }

    fn table(&mut self, count: i32) -> Result<Vec<f64>, String> {
        (0..count.max(0)).map(|_| self.number()).collect()
    }
}

impl Potential for Eam {
    fn description(&self) -> String {
        self.description.clone()
    }

    fn compute(&self, c: &Configuration, stress: bool) -> Result<PotentialResult, String> {
        let n = c.len();
        let mut type_ = Vec::with_capacity(n);
        for symbol in &c.symbols {
            let Some(t) = self.elements.iter().position(|e| e == symbol) else {
                return Err(format!("The EAM potential has no parameters for {symbol}"));
            };
            type_.push(t);
        }
        let pairs = half_list(c, self.cutoff)?;
        let mut rho = vec![0.0; n];
        let mut embedding_slope = vec![0.0; n];
        for pair in &pairs {
            let r = norm(pair.vector);
            let (i, j) = (pair.i, pair.j);
            rho[i] += self.density(type_[j], type_[i]).evaluate(r).0;
            rho[j] += self.density(type_[i], type_[j]).evaluate(r).0;
        }
        let mut result = PotentialResult { forces: vec![[0.0; 3]; n], ..Default::default() };
        for i in 0..n {
            let (value, slope) = self.embedding[type_[i]].evaluate(rho[i]);
            result.energy += value;
            // Linear continuation beyond the tabulated density, as LAMMPS.
            if rho[i] > self.rho_max {
                result.energy += slope * (rho[i] - self.rho_max);
            }
            embedding_slope[i] = slope;
        }
        let mut virial = [[0.0; 3]; 3];
        for pair in &pairs {
            let r = norm(pair.vector);
            let (i, j) = (pair.i, pair.j);
            let (_, d_rho_ji) = self.density(type_[j], type_[i]).evaluate(r);
            let (_, d_rho_ij) = self.density(type_[i], type_[j]).evaluate(r);
            let (a, b) = (type_[i].max(type_[j]), type_[i].min(type_[j]));
            let (rphi, d_rphi) = self.pair[a * (a + 1) / 2 + b].evaluate(r);
            let phi = rphi / r;
            let d_phi = (d_rphi - phi) / r;
            result.energy += phi;
            let de_dr = embedding_slope[i] * d_rho_ji + embedding_slope[j] * d_rho_ij + d_phi;
            let on_i = scale(pair.vector, de_dr / r);
            result.forces[i] = add(result.forces[i], on_i);
            result.forces[j] = sub(result.forces[j], on_i);
            add_virial(&mut virial, on_i, pair.vector);
        }
        finish_stress(&mut result, c, &virial, stress)?;
        Ok(result)
    }
}

// ------------------------------------------------- many-body helpers

const COULOMB: f64 = 14.3996454784255; // e^2 / (4 pi eps0) in eV Angstrom

/// Accumulates forces and the virial from gradients with respect to bond vectors.
struct Accumulator {
    result: PotentialResult,
    virial: Mat3,
}

impl Accumulator {
    fn new(atoms: usize) -> Self {
        Self { result: PotentialResult { forces: vec![[0.0; 3]; atoms], ..Default::default() }, virial: [[0.0; 3]; 3] }
    }

    fn bond(&mut self, i: usize, j: usize, d: Vec3, gradient: Vec3) {
        self.result.forces[i] = add(self.result.forces[i], gradient);
        self.result.forces[j] = sub(self.result.forces[j], gradient);
        for a in 0..3 {
            for b in 0..3 {
                self.virial[a][b] += 0.5 * (gradient[a] * d[b] + gradient[b] * d[a]);
            }
        }
    }

    fn finish(mut self, c: &Configuration, stress: bool) -> Result<PotentialResult, String> {
        let virial = self.virial;
        finish_stress(&mut self.result, c, &virial, stress)?;
        Ok(self.result)
    }
}

/// Full neighbour lists (both directions, periodic images included), grouped by atom.
fn neighbours_by_atom(c: &Configuration, cutoff: f64) -> Result<Vec<Vec<Neighbor>>, String> {
    let periodic = c.pbc.iter().any(|&p| p);
    let mut lists = vec![Vec::new(); c.len()];
    for n in neighbor_list(&c.positions, &if periodic { c.cell } else { identity() }, c.pbc, cutoff)? {
        lists[n.i].push(n);
    }
    Ok(lists)
}

fn tokens_of(file: &Path) -> Result<Vec<String>, String> {
    let lines = text_lines(file).ok_or_else(|| format!("Cannot read potential file {}", display(file)))?;
    let mut tokens = Vec::new();
    for line in &lines {
        let line = &line[..line.find('#').unwrap_or(line.len())];
        tokens.extend(words(line).map(str::to_string));
    }
    Ok(tokens)
}

fn file_number(text: &str, file: &Path) -> Result<f64, String> {
    match stod(text) {
        Ok((value, used)) if used == text.len() && value.is_finite() => Ok(value),
        _ => Err(format!("Invalid number '{text}' in {}", display(file))),
    }
}

type Triplet = (String, String, String);

fn lookup<'a, P>(table: &'a BTreeMap<Triplet, P>, a: &str, b: &str, c: &str, kind: &str) -> Result<&'a P, String> {
    table.get(&(a.to_string(), b.to_string(), c.to_string())).ok_or_else(|| format!("{kind} has no parameters for {a}-{b}-{c}"))
}

fn require_elements<P>(table: &BTreeMap<Triplet, P>, c: &Configuration, kind: &str) -> Result<(), String> {
    for a in &c.symbols {
        for b in &c.symbols {
            lookup(table, a, b, b, kind)?;
        }
    }
    Ok(())
}

// ------------------------------------------------------------ Tersoff

#[derive(Clone, Copy, Debug)]
struct TersoffParameters {
    m: f64,
    gamma: f64,
    lambda3: f64,
    c: f64,
    d: f64,
    h: f64,
    n: f64,
    beta: f64,
    lambda2: f64,
    b: f64,
    r: f64,
    dd: f64,
    lambda1: f64,
    a: f64,
}

/// (value, derivative) of the Tersoff cutoff function.
fn cutoff_function(r: f64, big_r: f64, d: f64) -> (f64, f64) {
    if r < big_r - d {
        return (1.0, 0.0);
    }
    if r > big_r + d {
        return (0.0, 0.0);
    }
    let x = 0.5 * PI * (r - big_r) / d;
    (0.5 - 0.5 * x.sin(), -0.25 * PI / d * x.cos())
}

pub struct Tersoff {
    table: BTreeMap<Triplet, TersoffParameters>,
    label: String,
    cutoff: f64,
}

impl Tersoff {
    fn with(table: BTreeMap<Triplet, TersoffParameters>, label: String) -> Self {
        let mut cutoff = 0.0;
        for p in table.values() {
            cutoff = max(cutoff, p.r + p.dd);
        }
        Self { table, label, cutoff }
    }

    /// From a LAMMPS tersoff file; None: built-in Si (Tersoff 1988).
    pub fn new(file: Option<&Path>) -> Result<Self, String> {
        let Some(file) = file else {
            let p = TersoffParameters {
                m: 3.0,
                gamma: 1.0,
                lambda3: 0.0,
                c: 1.0039e5,
                d: 16.217,
                h: -0.59825,
                n: 0.78734,
                beta: 1.1e-6,
                lambda2: 1.7322,
                b: 471.18,
                r: 2.85,
                dd: 0.15,
                lambda1: 2.4799,
                a: 1830.8,
            };
            return Ok(Self::with(BTreeMap::from([(("Si".into(), "Si".into(), "Si".into()), p)]), "Si, Tersoff 1988".into()));
        };
        let tokens = tokens_of(file)?;
        if tokens.is_empty() || tokens.len() % 17 != 0 {
            return Err(format!("A Tersoff file needs 17 entries per line: {}", display(file)));
        }
        let mut table = BTreeMap::new();
        for t in (0..tokens.len()).step_by(17) {
            let mut v = [0.0; 14];
            for k in 0..14 {
                v[k] = file_number(&tokens[t + 3 + k], file)?;
            }
            let p = TersoffParameters {
                m: v[0],
                gamma: v[1],
                lambda3: v[2],
                c: v[3],
                d: v[4],
                h: v[5],
                n: v[6],
                beta: v[7],
                lambda2: v[8],
                b: v[9],
                r: v[10],
                dd: v[11],
                lambda1: v[12],
                a: v[13],
            };
            if p.m != 1.0 && p.m != 3.0 {
                return Err("Tersoff m must be 1 or 3".into());
            }
            table.insert((tokens[t].clone(), tokens[t + 1].clone(), tokens[t + 2].clone()), p);
        }
        Ok(Self::with(table, file_name_of(file)))
    }
}

impl Potential for Tersoff {
    fn description(&self) -> String {
        format!("Tersoff ({})", self.label)
    }

    fn compute(&self, c: &Configuration, stress: bool) -> Result<PotentialResult, String> {
        require_elements(&self.table, c, "Tersoff")?;
        let mut acc = Accumulator::new(c.len());
        let lists = neighbours_by_atom(c, self.cutoff)?;
        for i in 0..c.len() {
            let list = &lists[i];
            let si = &c.symbols[i];
            for a in 0..list.len() {
                let nj = &list[a];
                let sj = &c.symbols[nj.j];
                let pij = lookup(&self.table, si, sj, sj, "Tersoff")?;
                let dij = nj.vector;
                let rij = norm(dij);
                if rij >= pij.r + pij.dd {
                    continue;
                }
                let (fc, dfc) = cutoff_function(rij, pij.r, pij.dd);
                // Bond order: zeta over the other neighbours k of i.
                let mut zeta = 0.0;
                let mut d_zeta_dij = [0.0; 3];
                let mut d_zeta_dik: Vec<(usize, Vec3)> = Vec::new();
                for b in 0..list.len() {
                    if b == a {
                        continue;
                    }
                    let nk = &list[b];
                    let sk = &c.symbols[nk.j];
                    let pijk = lookup(&self.table, si, sj, sk, "Tersoff")?;
                    let dik = nk.vector;
                    let rik = norm(dik);
                    if rik >= pijk.r + pijk.dd {
                        continue;
                    }
                    let (fck, dfck) = cutoff_function(rik, pijk.r, pijk.dd);
                    let cosine = dot(dij, dik) / (rij * rik);
                    let h = pijk.h - cosine;
                    let denominator = pijk.d * pijk.d + h * h;
                    let g = pijk.gamma * (1.0 + pijk.c * pijk.c / (pijk.d * pijk.d) - pijk.c * pijk.c / denominator);
                    let dg = -pijk.gamma * 2.0 * pijk.c * pijk.c * h / (denominator * denominator); // dg/dcos
                    let l3m = pijk.lambda3.powf(pijk.m);
                    let delta = rij - rik;
                    let ex = (l3m * delta.powf(pijk.m)).exp();
                    let dex = ex * l3m * pijk.m * delta.powf(pijk.m - 1.0); // d ex / d(rij - rik)
                    zeta += fck * g * ex;
                    let uij = scale(dij, 1.0 / rij);
                    let uik = scale(dik, 1.0 / rik);
                    let dcos_dij = sub(scale(dik, 1.0 / (rij * rik)), scale(dij, cosine / (rij * rij)));
                    let dcos_dik = sub(scale(dij, 1.0 / (rij * rik)), scale(dik, cosine / (rik * rik)));
                    d_zeta_dij = add(d_zeta_dij, add(scale(uij, fck * g * dex), scale(dcos_dij, fck * dg * ex)));
                    d_zeta_dik.push((b, add(add(scale(uik, dfck * g * ex - fck * g * dex), scale(dcos_dik, fck * dg * ex)), [0.0; 3])));
                }
                let (mut bij, mut db_dzeta) = (1.0, 0.0);
                if zeta > 0.0 && pij.beta > 0.0 {
                    let t = (pij.beta * zeta).powf(pij.n);
                    bij = (1.0 + t).powf(-0.5 / pij.n);
                    db_dzeta = -0.5 * (1.0 + t).powf(-0.5 / pij.n - 1.0) * t / zeta;
                }
                let fr = pij.a * (-pij.lambda1 * rij).exp();
                let dfr = -pij.lambda1 * fr;
                let fa = -pij.b * (-pij.lambda2 * rij).exp();
                let dfa = -pij.lambda2 * fa;
                acc.result.energy += 0.5 * fc * (fr + bij * fa);
                let radial = 0.5 * (dfc * (fr + bij * fa) + fc * (dfr + bij * dfa));
                let bond_order = 0.5 * fc * fa * db_dzeta;
                acc.bond(i, nj.j, dij, add(scale(dij, radial / rij), scale(d_zeta_dij, bond_order)));
                for (b, gradient) in d_zeta_dik {
                    acc.bond(i, list[b].j, list[b].vector, scale(gradient, bond_order));
                }
            }
        }
        acc.finish(c, stress)
    }
}

// --------------------------------------------------- Stillinger-Weber

#[derive(Clone, Copy, Debug)]
struct SwParameters {
    epsilon: f64,
    sigma: f64,
    a: f64,
    lambda: f64,
    gamma: f64,
    cos0: f64,
    big_a: f64,
    big_b: f64,
    p: f64,
    q: f64,
    _tol: f64,
}

pub struct StillingerWeber {
    table: BTreeMap<Triplet, SwParameters>,
    label: String,
    cutoff: f64,
}

impl StillingerWeber {
    fn with(table: BTreeMap<Triplet, SwParameters>, label: String) -> Self {
        let mut cutoff = 0.0;
        for p in table.values() {
            cutoff = max(cutoff, p.a * p.sigma);
        }
        Self { table, label, cutoff }
    }

    /// From a LAMMPS sw file; None: built-in Si (Stillinger and Weber 1985).
    pub fn new(file: Option<&Path>) -> Result<Self, String> {
        let Some(file) = file else {
            let p = SwParameters {
                epsilon: 2.1683,
                sigma: 2.0951,
                a: 1.80,
                lambda: 21.0,
                gamma: 1.20,
                cos0: -1.0 / 3.0,
                big_a: 7.049556277,
                big_b: 0.6022245584,
                p: 4.0,
                q: 0.0,
                _tol: 0.0,
            };
            return Ok(Self::with(BTreeMap::from([(("Si".into(), "Si".into(), "Si".into()), p)]), "Si, Stillinger-Weber 1985".into()));
        };
        let tokens = tokens_of(file)?;
        if tokens.is_empty() || tokens.len() % 14 != 0 {
            return Err(format!("A Stillinger-Weber file needs 14 entries per line: {}", display(file)));
        }
        let mut table = BTreeMap::new();
        for t in (0..tokens.len()).step_by(14) {
            let mut v = [0.0; 11];
            for k in 0..11 {
                v[k] = file_number(&tokens[t + 3 + k], file)?;
            }
            let p = SwParameters {
                epsilon: v[0],
                sigma: v[1],
                a: v[2],
                lambda: v[3],
                gamma: v[4],
                cos0: v[5],
                big_a: v[6],
                big_b: v[7],
                p: v[8],
                q: v[9],
                _tol: v[10],
            };
            table.insert((tokens[t].clone(), tokens[t + 1].clone(), tokens[t + 2].clone()), p);
        }
        Ok(Self::with(table, file_name_of(file)))
    }
}

impl Potential for StillingerWeber {
    fn description(&self) -> String {
        format!("Stillinger-Weber ({})", self.label)
    }

    fn compute(&self, c: &Configuration, stress: bool) -> Result<PotentialResult, String> {
        let kind = "Stillinger-Weber";
        require_elements(&self.table, c, kind)?;
        let mut acc = Accumulator::new(c.len());
        let lists = neighbours_by_atom(c, self.cutoff)?;
        for i in 0..c.len() {
            let list = &lists[i];
            let si = &c.symbols[i];
            // Two-body term, half per direction.
            for nj in list {
                let sj = &c.symbols[nj.j];
                let p = lookup(&self.table, si, sj, sj, kind)?;
                let r = norm(nj.vector);
                let cut = p.a * p.sigma;
                if r >= cut {
                    continue;
                }
                let s = p.sigma / r;
                let sp = s.powf(p.p);
                let sq = s.powf(p.q);
                let ex = (p.sigma / (r - cut)).exp();
                let phi = p.big_a * p.epsilon * (p.big_b * sp - sq) * ex;
                let dphi = p.big_a * p.epsilon * ex * ((-p.p * p.big_b * sp + p.q * sq) / r - (p.big_b * sp - sq) * p.sigma / ((r - cut) * (r - cut)));
                acc.result.energy += 0.5 * phi;
                acc.bond(i, nj.j, nj.vector, scale(nj.vector, 0.5 * dphi / r));
            }
            // Three-body term for each pair of neighbours j < k around i.
            for a in 0..list.len() {
                for b in a + 1..list.len() {
                    let sj = &c.symbols[list[a].j];
                    let sk = &c.symbols[list[b].j];
                    let pij = lookup(&self.table, si, sj, sj, kind)?;
                    let pik = lookup(&self.table, si, sk, sk, kind)?;
                    let pijk = lookup(&self.table, si, sj, sk, kind)?;
                    let (dij, dik) = (list[a].vector, list[b].vector);
                    let (rij, rik) = (norm(dij), norm(dik));
                    let (cutj, cutk) = (pij.a * pij.sigma, pik.a * pik.sigma);
                    if rij >= cutj || rik >= cutk {
                        continue;
                    }
                    let ej = (pij.gamma * pij.sigma / (rij - cutj)).exp();
                    let ek = (pik.gamma * pik.sigma / (rik - cutk)).exp();
                    let cosine = dot(dij, dik) / (rij * rik);
                    let delta = cosine - pijk.cos0;
                    let le = pijk.lambda * pijk.epsilon;
                    acc.result.energy += le * delta * delta * ej * ek;
                    let de_j = -le * delta * delta * ej * ek * pij.gamma * pij.sigma / ((rij - cutj) * (rij - cutj));
                    let de_k = -le * delta * delta * ej * ek * pik.gamma * pik.sigma / ((rik - cutk) * (rik - cutk));
                    let de_cos = 2.0 * le * delta * ej * ek;
                    let dcos_dij = sub(scale(dik, 1.0 / (rij * rik)), scale(dij, cosine / (rij * rij)));
                    let dcos_dik = sub(scale(dij, 1.0 / (rij * rik)), scale(dik, cosine / (rik * rik)));
                    acc.bond(i, list[a].j, dij, add(scale(dij, de_j / rij), scale(dcos_dij, de_cos)));
                    acc.bond(i, list[b].j, dik, add(scale(dik, de_k / rik), scale(dcos_dik, de_cos)));
                }
            }
        }
        acc.finish(c, stress)
    }
}

// ------------------------------------------------ Buckingham + Coulomb

#[derive(Clone, Copy, Debug)]
struct BuckinghamPair {
    a: f64,
    rho: f64,
    c: f64,
}

pub struct Buckingham {
    pairs: BTreeMap<(String, String), BuckinghamPair>,
    charges: BTreeMap<String, f64>,
    cutoff: f64,
    accuracy: f64,
}

impl Buckingham {
    /// {"pairs": [{"elements": [a, b], "A", "rho", "C"}], "charges": {element: q},
    /// "cutoff": 10, "ewald_accuracy": 1e-6}.
    pub fn new(options: &Json) -> Result<Self, String> {
        let mut pairs = BTreeMap::new();
        if let Some(list) = options.find("pairs").filter(|l| l.is_array()) {
            for entry in list.items()? {
                let elements = entry.at("elements")?.items()?;
                if elements.len() != 2 {
                    return Err("Each Buckingham pair needs two elements".into());
                }
                let pair = BuckinghamPair {
                    a: entry.at("A")?.number()?,
                    rho: entry.at("rho")?.number()?,
                    c: if entry.contains("C") { entry.at("C")?.number()? } else { 0.0 },
                };
                if !(pair.rho > 0.0) {
                    return Err("Buckingham rho must be positive".into());
                }
                pairs.insert((elements[0].string()?.to_string(), elements[1].string()?.to_string()), pair);
            }
        }
        let mut charges = BTreeMap::new();
        if let Some(table) = options.find("charges").filter(|t| t.is_object()) {
            for (element, charge) in table.members()? {
                charges.insert(element.clone(), charge.number()?);
            }
        }
        if pairs.is_empty() && charges.is_empty() {
            return Err("The Buckingham potential needs \"pairs\" and/or \"charges\"".into());
        }
        let cutoff = if options.contains("cutoff") { options.at("cutoff")?.number()? } else { 10.0 };
        let accuracy = if options.contains("ewald_accuracy") { options.at("ewald_accuracy")?.number()? } else { 1e-6 };
        if !(cutoff > 0.0) || !(accuracy > 0.0 && accuracy < 1.0) {
            return Err("Buckingham cutoff must be positive and ewald_accuracy in (0, 1)".into());
        }
        Ok(Self { pairs, charges, cutoff, accuracy })
    }

    fn pair_for(&self, a: &str, b: &str) -> Option<&BuckinghamPair> {
        self.pairs.get(&(a.to_string(), b.to_string())).or_else(|| self.pairs.get(&(b.to_string(), a.to_string())))
    }

    /// Reciprocal-space and self terms of the Ewald sum, with forces and virial.
    fn reciprocal(&self, c: &Configuration, q: &[f64], alpha: f64, acc: &mut Accumulator) -> Result<(), String> {
        let volume = cell_volume(&c.cell);
        let inv = inverse(&c.cell)?;
        let mut b = [[0.0; 3]; 3]; // reciprocal rows (2 pi included)
        for r in 0..3 {
            for k in 0..3 {
                b[r][k] = 2.0 * PI * inv[k][r];
            }
        }
        let k_max = 2.0 * alpha * (-self.accuracy.ln()).sqrt();
        let limits: [i32; 3] = std::array::from_fn(|r| (k_max * norm(c.cell[r]) / (2.0 * PI)).ceil() as i32 + 1);
        let mut self_term = 0.0;
        for qi in q {
            self_term += qi * qi;
        }
        acc.result.energy -= COULOMB * alpha / PI.sqrt() * self_term;
        for h in -limits[0]..=limits[0] {
            for k in -limits[1]..=limits[1] {
                for l in -limits[2]..=limits[2] {
                    if h == 0 && k == 0 && l == 0 {
                        continue;
                    }
                    let g = add(add(scale(b[0], h as f64), scale(b[1], k as f64)), scale(b[2], l as f64));
                    let g2 = dot(g, g);
                    if g2 > k_max * k_max {
                        continue;
                    }
                    let (mut sc, mut ss) = (0.0, 0.0);
                    let mut phase = vec![0.0; c.len()];
                    for i in 0..c.len() {
                        phase[i] = dot(g, c.positions[i]);
                        sc += q[i] * phase[i].cos();
                        ss += q[i] * phase[i].sin();
                    }
                    let prefactor = COULOMB * 2.0 * PI / volume * (-g2 / (4.0 * alpha * alpha)).exp() / g2;
                    let structure2 = sc * sc + ss * ss;
                    // E_g = (2 pi / V) exp(-g^2 / 4 alpha^2) / g^2 |S(g)|^2, summed over all g != 0.
                    let eg = prefactor * structure2;
                    acc.result.energy += eg;
                    for i in 0..c.len() {
                        // F_i = -dE/dr_i = 2 prefactor q_i g (sin(g.r_i) S_c - cos(g.r_i) S_s).
                        let weight = 2.0 * prefactor * q[i] * (phase[i].sin() * sc - phase[i].cos() * ss);
                        acc.result.forces[i] = add(acc.result.forces[i], scale(g, weight));
                    }
                    // dE_g/d strain = E_g (-delta_ab + 2 (1 + g^2 / 4 alpha^2) g_a g_b / g^2).
                    let factor = 2.0 * (1.0 + g2 / (4.0 * alpha * alpha)) / g2;
                    for a in 0..3 {
                        for d in 0..3 {
                            acc.virial[a][d] += eg * (factor * g[a] * g[d] - if a == d { 1.0 } else { 0.0 });
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

impl Potential for Buckingham {
    fn description(&self) -> String {
        if self.charges.is_empty() { "Buckingham" } else { "Buckingham + Coulomb (Ewald)" }.into()
    }

    fn compute(&self, c: &Configuration, stress: bool) -> Result<PotentialResult, String> {
        let mut acc = Accumulator::new(c.len());
        let mut q = vec![0.0; c.len()];
        let mut total = 0.0;
        for i in 0..c.len() {
            if !self.charges.is_empty() {
                let Some(charge) = self.charges.get(&c.symbols[i]) else {
                    return Err(format!("No charge given for {}", c.symbols[i]));
                };
                q[i] = *charge;
            }
            total += q[i];
        }
        let coulomb = !self.charges.is_empty();
        let periodic = c.fully_periodic();
        if coulomb && periodic && f64::abs(total) > 1e-8 {
            return Err("Ewald summation needs a charge-neutral cell".into());
        }
        if coulomb && !periodic && c.pbc.iter().any(|&p| p) {
            return Err("Coulomb interactions need a fully periodic or an open structure".into());
        }
        // Ewald splitting: real-space cutoff r_c with erfc(alpha r_c) ~ accuracy.
        let alpha = if coulomb && periodic { (-self.accuracy.ln()).sqrt() / self.cutoff } else { 0.0 };
        let real_cutoff = if coulomb && !periodic { 1e9 } else { self.cutoff };
        let open = !c.pbc.iter().any(|&p| p);
        // Real-space pairs (each once).
        let mut pairs = Vec::new();
        if open && coulomb {
            // Direct Coulomb sum over all pairs of an open structure.
            for i in 0..c.len() {
                for j in i + 1..c.len() {
                    pairs.push(Neighbor { i, j, shift: [0; 3], vector: sub(c.positions[j], c.positions[i]) });
                }
            }
        } else {
            for n in neighbor_list(&c.positions, &if open { identity() } else { c.cell }, c.pbc, real_cutoff)? {
                let mut keep = n.i < n.j;
                if n.i == n.j {
                    for k in 0..3 {
                        if n.shift[k] != 0 {
                            keep = n.shift[k] > 0;
                            break;
                        }
                    }
                }
                if keep {
                    pairs.push(n);
                }
            }
        }
        for n in &pairs {
            let r = norm(n.vector);
            let (mut energy, mut derivative) = (0.0, 0.0);
            if r < self.cutoff {
                if let Some(b) = self.pair_for(&c.symbols[n.i], &c.symbols[n.j]) {
                    let ex = b.a * (-r / b.rho).exp();
                    let r6 = r.powf(6.0);
                    energy += ex - b.c / r6;
                    derivative += -ex / b.rho + 6.0 * b.c / (r6 * r);
                }
            }
            if coulomb {
                let qq = COULOMB * q[n.i] * q[n.j];
                if periodic {
                    let e = erfc(alpha * r);
                    energy += qq * e / r;
                    derivative += -qq * (e / (r * r) + 2.0 * alpha / PI.sqrt() * (-alpha * alpha * r * r).exp() / r);
                } else {
                    energy += qq / r;
                    derivative += -qq / (r * r);
                }
            }
            acc.result.energy += energy;
            acc.bond(n.i, n.j, n.vector, scale(n.vector, derivative / r));
        }
        if coulomb && periodic {
            self.reciprocal(c, &q, alpha, &mut acc)?;
        }
        acc.finish(c, stress)
    }
}

// ------------------------------------------------------------- factory

fn kwarg(object: &Json, name: &str, fallback: f64) -> Result<f64, String> {
    match object.find(name) {
        None => Ok(fallback),
        Some(v) if v.is_null() => Ok(fallback),
        Some(v) => positive(v.number()?, name),
    }
}

/// Parses {"potential": "EMT"}, {"potential": "LennardJones", "epsilon",
/// "sigma", "cutoff"}, {"potential": "EAM", "file", "format"},
/// {"potential": "Tersoff" | "StillingerWeber", "file" (optional)} or
/// {"potential": "Buckingham", ...}; relative files resolve against base.
/// ASE-style {"module", "attribute", "kwargs"} descriptions of EMT and
/// LennardJones are accepted.
pub fn potential_factory(specification: &Json, base: &Path) -> Result<PotentialHandle, String> {
    if !specification.is_object() {
        return Err("The interatomic potential must be an object such as {\"potential\": \"EMT\"}".into());
    }
    let name: String;
    let mut options = Json::object();
    if let Some(potential) = specification.find("potential") {
        name = potential.string()?.to_string();
        options = specification.clone();
    } else if specification.contains("module") && specification.contains("attribute") {
        // Accept the earlier ASE-style description for the native equivalents.
        name = specification.at("attribute")?.string()?.to_string();
        if let Some(kwargs) = specification.find("kwargs").filter(|k| !k.is_null()) {
            options = kwargs.clone();
        }
        if options.contains("rc") {
            let rc = options.at("rc")?.clone();
            options.set("cutoff", rc);
        }
    } else {
        return Err("The interatomic potential must name \"potential\": EMT or LennardJones".into());
    }
    let key: String = name.chars().filter(|c| c.is_ascii_alphanumeric()).map(|c| c.to_ascii_lowercase()).collect();
    match key.as_str() {
        "emt" => Ok(Arc::new(Emt::new())),
        "lennardjones" | "lj" => {
            let sigma = kwarg(&options, "sigma", 1.0)?;
            let epsilon = kwarg(&options, "epsilon", 1.0)?;
            let cutoff = kwarg(&options, "cutoff", 3.0 * sigma)?;
            Ok(Arc::new(LennardJones::new(epsilon, sigma, cutoff)?))
        }
        "eam" | "eamalloy" | "eamfs" | "finnissinclair" => {
            let file = options.find("file").filter(|f| f.is_string() && !f.string().unwrap_or("").is_empty());
            let Some(file) = file else {
                return Err("The EAM potential needs a \"file\"".into());
            };
            let path = join_normal(base, file.string()?);
            let format = if options.contains("format") {
                options.at("format")?.string()?.to_string()
            } else if key == "eamfs" || key == "finnissinclair" {
                "fs".to_string()
            } else {
                "auto".to_string()
            };
            Ok(Arc::new(Eam::new(&path, &format)?))
        }
        "tersoff" | "stillingerweber" | "sw" => {
            let path: Option<PathBuf> = match options.find("file") {
                Some(f) if f.is_string() && !f.string()?.is_empty() => Some(join_normal(base, f.string()?)),
                _ => None,
            };
            if key == "tersoff" {
                Ok(Arc::new(Tersoff::new(path.as_deref())?))
            } else {
                Ok(Arc::new(StillingerWeber::new(path.as_deref())?))
            }
        }
        "buckingham" | "buck" => Ok(Arc::new(Buckingham::new(&options)?)),
        _ => Err(format!(
            "Unsupported interatomic potential '{name}'; AtomForge provides EMT, LennardJones, EAM, Tersoff, StillingerWeber and Buckingham natively"
        )),
    }
}
