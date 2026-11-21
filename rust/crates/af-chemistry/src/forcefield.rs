//! The Universal Force Field (the C++ chemistry/ForceField): atom typing,
//! bond, angle, torsion and inversion terms, Lennard-Jones and Coulomb pairs
//! (periodic images included), restraints, and an L-BFGS minimiser.

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::OnceLock;

use af_core::structure::Structure;

use crate::molecule::*;
use crate::uff_parameters::{UffParameter, UFF_TABLE};

const PI: f64 = 3.14159265358979323846;
/// kcal/mol A / e^2.
const COULOMB: f64 = 332.0637;

fn parameter(name: &str) -> Option<&'static UffParameter> {
    static INDEX: OnceLock<HashMap<&'static str, &'static UffParameter>> = OnceLock::new();
    INDEX.get_or_init(|| UFF_TABLE.iter().map(|p| (p.name, p)).collect()).get(name).copied()
}

fn accumulate(g: &mut Option<&mut Vec<Vec3>>, i: i32, v: Vec3, s: f64) {
    if let Some(g) = g.as_deref_mut() {
        let gi = &mut g[i as usize];
        gi[0] += v[0] * s;
        gi[1] += v[1] * s;
        gi[2] += v[2] * s;
    }
}

fn group16(z: i32) -> bool {
    matches!(z, 8 | 16 | 34 | 52 | 84)
}

fn geometry_of(t: &str) -> u8 {
    t.as_bytes().get(2).copied().unwrap_or(b'_')
}

fn is_sp2(t: &str) -> bool {
    matches!(geometry_of(t), b'2' | b'R')
}

fn is_sp3(t: &str) -> bool {
    geometry_of(t) == b'3'
}

/// Bond order for UFF: aromatic 1.5, amide C-N 1.41, else the integer order.
fn uff_order(m: &Molecule, k: i32) -> f64 {
    let b = &m.bonds[k as usize];
    if b.dative {
        return 1.0;
    }
    if b.aromatic {
        return 1.5;
    }
    if b.order == 1 {
        let amide = |c: i32, n: i32| {
            m.atoms[c as usize].z == 6
                && m.atoms[n as usize].z == 7
                && m.bonds_of(c).iter().any(|&e| m.bonds[e as usize].order == 2 && m.atoms[m.other(e, c) as usize].z == 8)
        };
        if amide(b.a, b.b) || amide(b.b, b.a) {
            return 1.41;
        }
    }
    f64::from(b.order)
}

/// The image vector of bond k seen from atom `from`.
fn image_from(m: &Molecule, k: i32, from: i32) -> Vec3 {
    let b = &m.bonds[k as usize];
    if b.a == from {
        b.image
    } else {
        scale(b.image, -1.0)
    }
}

/// Dihedral of p1-p2-p3-p4 from b1 = p2-p1, b2 = p3-p2, b3 = p4-p3, with the
/// gradients for the four points (Blondel and Karplus 1996).
pub(crate) fn dihedral(b1: Vec3, b2: Vec3, b3: Vec3, g: Option<&mut [Vec3; 4]>) -> f64 {
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let length = norm(b2);
    let y = length * dot(b1, n2);
    let x = dot(n1, n2);
    let phi = y.atan2(x);
    if let Some(g) = g {
        let a2 = fmax(dot(n1, n1), 1e-12);
        let b2n = fmax(dot(n2, n2), 1e-12);
        let l = fmax(length, 1e-12);
        let f12 = dot(b1, b2);
        let f32_ = dot(b3, b2);
        g[0] = scale(n1, -l / a2);
        g[3] = scale(n2, l / b2n);
        g[1] = add(scale(n1, l / a2 + f12 / (a2 * l)), scale(n2, f32_ / (b2n * l)));
        g[2] = add(scale(n2, -l / b2n - f32_ / (b2n * l)), scale(n1, -f12 / (a2 * l)));
    }
    phi
}

fn element_prefix(z: i32) -> String {
    let mut s = symbol_of(z);
    if s.len() == 1 {
        s.push('_');
    }
    s
}

/// Whether UFF parameters exist for an element.
pub fn uff_has_element(z: i32) -> bool {
    let prefix = element_prefix(z);
    UFF_TABLE.iter().any(|p| p.name.as_bytes().get(..2) == Some(prefix.as_bytes()))
}

/// Angle theta0 (degrees) of a UFF type.
pub fn uff_angle(t: &str) -> f64 {
    parameter(t).map_or(109.47, |p| p.theta0)
}

/// vdW radius x1/2 (Angstrom) of a UFF type.
pub fn uff_vdw_radius(t: &str) -> f64 {
    parameter(t).map_or(1.7, |p| p.x1 / 2.0)
}

/// The UFF type of an atom (e.g. C_R, N_3, O_2, Fe6+2).
pub fn uff_type(m: &Molecule, i: i32) -> Result<String, String> {
    let a = &m.atoms[i as usize];
    let z = a.z;
    let hyb = hybridisation(m, i);
    let mut sum = 0;
    for &k in m.bonds_of(i) {
        if !m.bonds[k as usize].dative {
            sum += m.bonds[k as usize].order;
        }
    }
    sum += a.hydrogens;
    let degree = m.degree(i) + a.hydrogens;
    let conjugated = || {
        // A neighbour (through a single bond) with a pi bond, or a lone-pair N/O next to a double bond.
        for &k in m.bonds_of(i) {
            if m.bonds[k as usize].dative {
                continue;
            }
            let b = m.other(k, i);
            for &e in m.bonds_of(b) {
                let be = &m.bonds[e as usize];
                if e == k || be.dative {
                    continue;
                }
                if be.order >= 2 || be.aromatic {
                    return true;
                }
            }
            if m.bonds[k as usize].order >= 2 && matches!(m.atoms[b as usize].z, 7 | 8) {
                for &e in m.bonds_of(i) {
                    if e != k && matches!(m.atoms[m.other(e, i) as usize].z, 7 | 8) {
                        return true;
                    }
                }
            }
        }
        false
    };
    let named = match z {
        1 => Some("H_"),
        6 => Some(if hyb == 1 {
            "C_1"
        } else if hyb == 2 {
            if a.aromatic || conjugated() {
                "C_R"
            } else {
                "C_2"
            }
        } else {
            "C_3"
        }),
        7 => Some(if hyb == 1 {
            "N_1"
        } else if hyb == 2 {
            if a.aromatic || conjugated() || sum == degree {
                "N_R"
            } else {
                "N_2"
            }
        } else {
            "N_3"
        }),
        8 => Some(if hyb == 1 {
            "O_1"
        } else if a.aromatic {
            "O_R"
        } else if hyb == 2 {
            "O_2"
        } else if degree == 2 && conjugated() {
            "O_R"
        } else {
            "O_3"
        }),
        16 => Some(if a.aromatic {
            "S_R"
        } else if degree == 1 && sum == 2 {
            "S_2"
        } else if sum >= 6 {
            "S_3+6"
        } else if sum >= 3 {
            "S_3+4"
        } else {
            "S_3+2"
        }),
        15 => {
            let metal = m.bonds_of(i).iter().any(|&k| m.bonds[k as usize].dative);
            Some(if metal {
                "P_3+q"
            } else if sum >= 5 || degree >= 4 {
                "P_3+5"
            } else {
                "P_3+3"
            })
        }
        5 => Some(if degree >= 4 { "B_3" } else { "B_2" }),
        9 => Some("F_"),
        17 => Some("Cl"),
        35 => Some("Br"),
        53 => Some("I_"),
        _ => None,
    };
    if let Some(name) = named {
        return Ok(name.to_string());
    }
    // Other elements: the type of that element whose coordination is closest.
    let prefix = element_prefix(z);
    let coordination = m.degree(i);
    let wanted = if coordination <= 2 {
        b'1'
    } else if coordination == 3 {
        b'2'
    } else if coordination == 4 {
        // Square planar when the four bonds lie in a plane.
        let mut sum4 = 0.0;
        let nb = m.bonds_of(i);
        for x in 0..nb.len() {
            for y in x + 1..nb.len() {
                let u = add(sub(m.atoms[m.other(nb[x], i) as usize].position, a.position), image_from(m, nb[x], i));
                let v = add(sub(m.atoms[m.other(nb[y], i) as usize].position, a.position), image_from(m, nb[y], i));
                let c = dot(u, v) / fmax(norm(u) * norm(v), 1e-12);
                sum4 += c.clamp(-1.0, 1.0).acos() * 180.0 / PI;
            }
        }
        if sum4 > 680.0 {
            b'4'
        } else {
            b'3'
        }
    } else if coordination == 5 {
        b'5'
    } else {
        b'6'
    };
    let code = |g: u8| -> i32 {
        if g == b'R' {
            2
        } else if g.is_ascii_digit() {
            i32::from(g - b'0')
        } else {
            3
        }
    };
    let mut best = "";
    let mut best_score = 1 << 20;
    for p in UFF_TABLE {
        let name = p.name;
        if name.as_bytes().get(..2) != Some(prefix.as_bytes()) {
            continue;
        }
        if name.len() == 2 && prefix.len() == 2 && prefix.as_bytes()[1] != b'_' && name != prefix {
            continue;
        }
        let score = (code(geometry_of(name)) - code(wanted)).abs() * 4 + if name.contains('+') { 0 } else { 1 };
        if score < best_score {
            best_score = score;
            best = name;
        }
    }
    if best.is_empty() {
        return Err(format!("UFF has no parameters for {}", symbol_of(z)));
    }
    Ok(best.to_string())
}

/// Natural UFF bond length between two types for a bond order (1.5 aromatic).
pub fn uff_bond_length(ti: &str, tj: &str, order: f64) -> Result<f64, String> {
    let (Some(pi), Some(pj)) = (parameter(ti), parameter(tj)) else {
        return Err(format!("unknown UFF type {}", if parameter(ti).is_some() { tj } else { ti }));
    };
    let (ri, rj) = (pi.r1, pj.r1);
    let rbo = if order > 0.0 { -0.1332 * (ri + rj) * order.ln() } else { 0.0 };
    let (ci, cj) = (pi.chi, pj.chi);
    let difference = ci.sqrt() - cj.sqrt();
    let ren = ri * rj * (difference * difference) / (ci * ri + cj * rj);
    Ok(ri + rj + rbo - ren)
}

#[derive(Clone, Copy, Debug)]
pub struct ForceFieldOptions {
    /// Coulomb with the atom charges (formal by default).
    pub electrostatics: bool,
    pub dielectric: f64,
    /// Non-bonded cutoff (Angstrom).
    pub cutoff: f64,
}

impl Default for ForceFieldOptions {
    fn default() -> Self {
        Self { electrostatics: false, dielectric: 1.0, cutoff: 12.0 }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Restraint {
    /// 2 distance (A), 3 angle (deg), 4 dihedral (deg).
    pub kind: i32,
    pub atoms: [i32; 4],
    pub target: f64,
    /// kcal/mol per A^2 or per rad^2.
    pub k: f64,
}

impl Default for Restraint {
    fn default() -> Self {
        Self { kind: 2, atoms: [0; 4], target: 0.0, k: 500.0 }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnergyTerms {
    pub bond: f64,
    pub angle: f64,
    pub torsion: f64,
    pub inversion: f64,
    pub vdw: f64,
    pub electrostatic: f64,
    pub restraint: f64,
}

impl EnergyTerms {
    pub fn total(&self) -> f64 {
        self.bond + self.angle + self.torsion + self.inversion + self.vdw + self.electrostatic + self.restraint
    }
}

#[derive(Clone, Debug)]
struct BondTerm {
    i: i32,
    j: i32,
    image: Vec3,
    k: f64,
    r0: f64,
}

#[derive(Clone, Debug)]
struct AngleTerm {
    i: i32,
    j: i32,
    k: i32,
    image_ji: Vec3,
    image_jk: Vec3,
    force: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    form: i32,
}

#[derive(Clone, Debug)]
struct TorsionTerm {
    i: i32,
    j: i32,
    k: i32,
    l: i32,
    image_ij: Vec3,
    image_jk: Vec3,
    image_kl: Vec3,
    v: f64,
    n: i32,
    cos_phase: f64,
}

#[derive(Clone, Debug)]
struct InversionTerm {
    centre: i32,
    j: i32,
    k: i32,
    l: i32,
    image_j: Vec3,
    image_k: Vec3,
    image_l: Vec3,
    force: f64,
    c0: f64,
    c1: f64,
    c2: f64,
}

#[derive(Clone, Debug)]
struct Pair {
    i: i32,
    j: i32,
    image: Vec3,
    d: f64,
    x: f64,
    q: f64,
}

/// UFF for a molecule (explicit hydrogens), optionally in a periodic cell.
#[derive(Clone, Debug)]
pub struct UffForceField {
    atoms: i32,
    types: Vec<String>,
    bonds: Vec<BondTerm>,
    angles: Vec<AngleTerm>,
    torsions: Vec<TorsionTerm>,
    inversions: Vec<InversionTerm>,
    vdw_distance: Vec<f64>,
    vdw_depth: Vec<f64>,
    /// Excluded (1-2, 1-3) pairs by i * n + j (i < j), with their image vectors.
    excluded: Vec<(i64, Vec3)>,
    options: ForceFieldOptions,
    periodic: bool,
    cell: Structure,
    pairs: RefCell<Vec<Pair>>,
    pair_positions: RefCell<Vec<Vec3>>,
    pub restraints: Vec<Restraint>,
    /// Per-atom charges for the Coulomb term (default: formal charges).
    pub charges: Vec<f64>,
}

impl UffForceField {
    /// `cell`: a periodic structure's cell and periodicity (non-bonded images);
    /// None or no cell for a molecule in vacuum.
    pub fn new(m: &Molecule, cell: Option<&Structure>, options: ForceFieldOptions) -> Result<Self, String> {
        let n = m.atoms.len() as i32;
        if m.atoms.iter().any(|a| a.hydrogens > 0) {
            return Err("UFF needs explicit hydrogens".into());
        }
        let types = (0..n).map(|i| uff_type(m, i)).collect::<Result<Vec<_>, _>>()?;
        let mut p: Vec<&UffParameter> = Vec::with_capacity(n as usize);
        let (mut vdw_distance, mut vdw_depth) = (Vec::new(), Vec::new());
        for t in &types {
            let q = parameter(t).ok_or_else(|| format!("UFF has no parameters for type {t}"))?;
            vdw_distance.push(q.x1);
            vdw_depth.push(q.d1);
            p.push(q);
        }
        let charges = m.atoms.iter().map(|a| f64::from(a.charge)).collect();
        let mut field = UffForceField {
            atoms: n,
            types,
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
            inversions: Vec::new(),
            vdw_distance,
            vdw_depth,
            excluded: Vec::new(),
            options,
            periodic: false,
            cell: Structure::default(),
            pairs: RefCell::new(Vec::new()),
            pair_positions: RefCell::new(Vec::new()),
            restraints: Vec::new(),
            charges,
        };
        if let Some(c) = cell.filter(|c| c.has_unit_cell && c.pbc.iter().any(|&b| b)) {
            field.periodic = true;
            field.cell.has_unit_cell = true;
            field.cell.cell_vectors = c.cell_vectors;
            field.cell.pbc = c.pbc;
            field.cell.cell_offset = c.cell_offset;
        }
        let pair_key = |i: i32, j: i32| i64::from(i.min(j)) * i64::from(n) + i64::from(i.max(j));
        let mut bond_length = vec![0.0; m.bonds.len()];
        for (k, b) in m.bonds.iter().enumerate() {
            let r0 = uff_bond_length(&field.types[b.a as usize], &field.types[b.b as usize], uff_order(m, k as i32))?;
            bond_length[k] = r0;
            let force = 664.12 * p[b.a as usize].z1 * p[b.b as usize].z1 / (r0 * r0 * r0);
            field.bonds.push(BondTerm { i: b.a, j: b.b, image: b.image, k: force, r0 });
            field.excluded.push((pair_key(b.a, b.b), if b.a < b.b { b.image } else { scale(b.image, -1.0) }));
        }
        // Angles.
        for j in 0..n {
            let nb = m.bonds_of(j);
            if nb.len() < 2 || m.atoms[j as usize].z == 1 {
                continue;
            }
            let octahedral = nb.len() == 6 && !is_metal(m.atoms[j as usize].z);
            let theta0 = (if octahedral { 90.0 } else { p[j as usize].theta0 }) * PI / 180.0;
            let g = geometry_of(&field.types[j as usize]);
            for x in 0..nb.len() {
                for y in x + 1..nb.len() {
                    let (i, k) = (m.other(nb[x], j), m.other(nb[y], j));
                    let (rij, rjk) = (bond_length[nb[x] as usize], bond_length[nb[y] as usize]);
                    let c0 = theta0.cos();
                    let rik = (rij * rij + rjk * rjk - 2.0 * rij * rjk * c0).sqrt();
                    let force = 664.12 * p[i as usize].z1 * p[k as usize].z1 / rik.powf(5.0) * (3.0 * rij * rjk * (1.0 - c0 * c0) - rik * rik * c0);
                    let mut t = AngleTerm { i, j, k, image_ji: image_from(m, nb[x], j), image_jk: image_from(m, nb[y], j), force, c0: 0.0, c1: 0.0, c2: 0.0, form: 0 };
                    // Hypervalent main-group centres with six ligands (SF6, PF6-) are octahedral.
                    if octahedral {
                        t.form = 4;
                    } else if g == b'1' || (p[j as usize].theta0 - 180.0).abs() < 1e-6 {
                        t.form = 1;
                    } else if g == b'4' || g == b'6' {
                        t.form = 4;
                    } else if g == b'2' || g == b'R' {
                        t.form = 3;
                    } else {
                        let s2 = fmax(1.0 - c0 * c0, 1e-8);
                        t.c2 = 1.0 / (4.0 * s2);
                        t.c1 = -4.0 * t.c2 * c0;
                        t.c0 = t.c2 * (2.0 * c0 * c0 + 1.0);
                    }
                    let image_ik = sub(t.image_jk, t.image_ji);
                    field.excluded.push((pair_key(i, k), if i < k { image_ik } else { scale(image_ik, -1.0) }));
                    field.angles.push(t);
                }
            }
        }
        // Torsions about bonds between sp2/sp3 main-group atoms.
        for (kb, b) in m.bonds.iter().enumerate() {
            if b.dative {
                continue;
            }
            let (j, k) = (b.a, b.b);
            let (tj, tk) = (field.types[j as usize].clone(), field.types[k as usize].clone());
            if !(is_sp2(&tj) || is_sp3(&tj)) || !(is_sp2(&tk) || is_sp3(&tk)) {
                continue;
            }
            if is_metal(m.atoms[j as usize].z) || is_metal(m.atoms[k as usize].z) {
                continue;
            }
            let (dj, dk) = (m.degree(j), m.degree(k));
            if dj < 2 || dk < 2 {
                continue;
            }
            let order = uff_order(m, kb as i32);
            let (mut v, mut cos_phase, mut n_fold): (f64, f64, i32);
            let torsion_v = |a: i32| {
                // sp3 torsional barrier V1 (group 16: O 2.0, others 6.8).
                if group16(m.atoms[a as usize].z) {
                    if m.atoms[a as usize].z == 8 {
                        2.0
                    } else {
                        6.8
                    }
                } else {
                    p[a as usize].v1
                }
            };
            if is_sp3(&tj) && is_sp3(&tk) {
                if group16(m.atoms[j as usize].z) && group16(m.atoms[k as usize].z) {
                    n_fold = 2;
                    cos_phase = -1.0;
                    v = (torsion_v(j) * torsion_v(k)).sqrt();
                } else {
                    n_fold = 3;
                    cos_phase = -1.0;
                    v = (p[j as usize].v1 * p[k as usize].v1).sqrt();
                }
            } else if is_sp2(&tj) && is_sp2(&tk) {
                n_fold = 2;
                cos_phase = 1.0;
                v = 5.0 * (p[j as usize].u1 * p[k as usize].u1).sqrt() * (1.0 + 4.18 * order.ln());
            } else {
                let sp2 = if is_sp2(&tj) { j } else { k };
                let sp3 = if sp2 == j { k } else { j };
                n_fold = 6;
                cos_phase = 1.0;
                v = 1.0;
                if group16(m.atoms[sp3 as usize].z) && !group16(m.atoms[sp2 as usize].z) {
                    n_fold = 2;
                    cos_phase = -1.0;
                    v = 5.0 * (p[j as usize].u1 * p[k as usize].u1).sqrt() * (1.0 + 4.18 * order.ln());
                } else {
                    for &e in m.bonds_of(sp2) {
                        let x = m.other(e, sp2);
                        if x != sp3 && is_sp2(&field.types[x as usize]) {
                            n_fold = 3;
                            cos_phase = -1.0;
                            v = 2.0;
                            break;
                        }
                    }
                }
            }
            if v == 0.0 {
                continue;
            }
            v /= f64::from((dj - 1) * (dk - 1));
            for &e1 in m.bonds_of(j) {
                let i = m.other(e1, j);
                if i == k {
                    continue;
                }
                for &e3 in m.bonds_of(k) {
                    let l = m.other(e3, k);
                    if l == j || l == i {
                        continue;
                    }
                    field.torsions.push(TorsionTerm {
                        i,
                        j,
                        k,
                        l,
                        image_ij: image_from(m, e1, i),
                        image_jk: b.image,
                        image_kl: image_from(m, e3, k),
                        v,
                        n: n_fold,
                        cos_phase,
                    });
                }
            }
        }
        // Inversions at sp2 carbon and nitrogen.
        for c in 0..n {
            let t = field.types[c as usize].as_str();
            if m.degree(c) != 3 {
                continue;
            }
            let force = if t == "C_2" || t == "C_R" {
                if m.bonds_of(c).iter().any(|&e| field.types[m.other(e, c) as usize] == "O_2") {
                    50.0
                } else {
                    6.0
                }
            } else if t == "N_2" || t == "N_R" {
                6.0
            } else {
                continue;
            };
            let nb = m.bonds_of(c);
            let (a0, a1, a2) = (m.other(nb[0], c), m.other(nb[1], c), m.other(nb[2], c));
            let (i0, i1, i2) = (image_from(m, nb[0], c), image_from(m, nb[1], c), image_from(m, nb[2], c));
            // The out-of-plane atom L in turn; E = K/3 (1 - cos w).
            for (j, k, l, ij, ik, il) in [(a0, a1, a2, i0, i1, i2), (a1, a2, a0, i1, i2, i0), (a2, a0, a1, i2, i0, i1)] {
                field.inversions.push(InversionTerm { centre: c, j, k, l, image_j: ij, image_k: ik, image_l: il, force: force / 3.0, c0: 1.0, c1: -1.0, c2: 0.0 });
            }
        }
        field.excluded.sort_by_key(|e| e.0);
        Ok(field)
    }

    pub fn types(&self) -> &[String] {
        &self.types
    }

    fn build_pairs(&self, x: &[Vec3]) -> Result<(), String> {
        let mut pairs = self.pairs.borrow_mut();
        pairs.clear();
        *self.pair_positions.borrow_mut() = x.to_vec();
        let n = self.atoms;
        let excluded = |i: i32, j: i32, image: Vec3| {
            let key = i64::from(i.min(j)) * i64::from(n) + i64::from(i.max(j));
            let im = if i < j { image } else { scale(image, -1.0) };
            let start = self.excluded.partition_point(|e| e.0 < key);
            self.excluded[start..].iter().take_while(|e| e.0 == key).any(|e| norm(sub(e.1, im)) < 1e-3)
        };
        let mut add_pair = |i: i32, j: i32, image: Vec3| {
            if excluded(i, j, image) {
                return;
            }
            let d = (self.vdw_depth[i as usize] * self.vdw_depth[j as usize]).sqrt();
            let xij = (self.vdw_distance[i as usize] * self.vdw_distance[j as usize]).sqrt();
            let q = if self.options.electrostatics { COULOMB * self.charges[i as usize] * self.charges[j as usize] / self.options.dielectric } else { 0.0 };
            pairs.push(Pair { i, j, image, d, x: xij, q });
        };
        if !self.periodic && n <= 1500 {
            for i in 0..n {
                for j in i + 1..n {
                    add_pair(i, j, [0.0; 3]);
                }
            }
            return Ok(());
        }
        let mut s = self.cell.clone();
        if !self.periodic {
            s.has_unit_cell = false;
            s.pbc = [false; 3];
        }
        s.atoms = x.iter().map(|p| af_core::structure::AtomSite { x: p[0], y: p[1], z: p[2], ..Default::default() }).collect();
        let reach = self.options.cutoff + 2.0;
        let within = af_core::neighbours::within(&s, reach, self.periodic).map_err(|e| e.to_string())?;
        for pr in within {
            if pr.j < pr.i {
                continue;
            }
            let d = [pr.d.x, pr.d.y, pr.d.z];
            if pr.i == pr.j {
                // A self image: count each lattice vector once.
                if d[0] < 0.0 || (d[0] == 0.0 && (d[1] < 0.0 || (d[1] == 0.0 && d[2] < 0.0))) {
                    continue;
                }
            }
            let image = sub(d, sub(x[pr.j], x[pr.i]));
            add_pair(pr.i as i32, pr.j as i32, image);
        }
        Ok(())
    }

    /// Energy (kcal/mol) and, with `gradient`, dE/dx (kcal/mol/A).
    pub fn energy(&self, x: &[Vec3], mut g: Option<&mut Vec<Vec3>>, terms: Option<&mut EnergyTerms>) -> Result<f64, String> {
        if let Some(g) = g.as_deref_mut() {
            *g = vec![[0.0; 3]; self.atoms as usize];
        }
        let want_gradient = g.is_some();
        let mut e = EnergyTerms::default();
        for t in &self.bonds {
            let d = add(sub(x[t.j as usize], x[t.i as usize]), t.image);
            let r = fmax(norm(d), 1e-9);
            let dr = r - t.r0;
            e.bond += 0.5 * t.k * dr * dr;
            let f = scale(d, t.k * dr / r);
            accumulate(&mut g, t.j, f, 1.0);
            accumulate(&mut g, t.i, f, -1.0);
        }
        for t in &self.angles {
            let u = add(sub(x[t.i as usize], x[t.j as usize]), t.image_ji);
            let v = add(sub(x[t.k as usize], x[t.j as usize]), t.image_jk);
            let lu = fmax(norm(u), 1e-9);
            let lv = fmax(norm(v), 1e-9);
            let c = (dot(u, v) / (lu * lv)).clamp(-1.0, 1.0);
            let (energy, de_dc) = match t.form {
                1 => (t.force * (1.0 + c), t.force),
                4 => (t.force * (c * c - c * c * c * c) / 2.0, t.force * (c - 2.0 * c * c * c)),
                3 => (t.force * (1.0 - (4.0 * c * c * c - 3.0 * c)) / 9.0, -t.force * (12.0 * c * c - 3.0) / 9.0),
                _ => (t.force * (t.c0 + t.c1 * c + t.c2 * (2.0 * c * c - 1.0)), t.force * (t.c1 + 4.0 * t.c2 * c)),
            };
            e.angle += energy;
            if want_gradient {
                let du = sub(scale(v, 1.0 / (lu * lv)), scale(u, c / (lu * lu)));
                let dv = sub(scale(u, 1.0 / (lu * lv)), scale(v, c / (lv * lv)));
                accumulate(&mut g, t.i, du, de_dc);
                accumulate(&mut g, t.k, dv, de_dc);
                accumulate(&mut g, t.j, add(du, dv), -de_dc);
            }
        }
        for t in &self.torsions {
            let b1 = add(sub(x[t.j as usize], x[t.i as usize]), t.image_ij);
            let b2 = add(sub(x[t.k as usize], x[t.j as usize]), t.image_jk);
            let b3 = add(sub(x[t.l as usize], x[t.k as usize]), t.image_kl);
            let mut grad = [[0.0; 3]; 4];
            let phi = dihedral(b1, b2, b3, if want_gradient { Some(&mut grad) } else { None });
            let nphi = f64::from(t.n) * phi;
            e.torsion += 0.5 * t.v * (1.0 - t.cos_phase * nphi.cos());
            if want_gradient {
                let de_dphi = 0.5 * t.v * t.cos_phase * f64::from(t.n) * nphi.sin();
                accumulate(&mut g, t.i, grad[0], de_dphi);
                accumulate(&mut g, t.j, grad[1], de_dphi);
                accumulate(&mut g, t.k, grad[2], de_dphi);
                accumulate(&mut g, t.l, grad[3], de_dphi);
            }
        }
        for t in &self.inversions {
            let c = x[t.centre as usize];
            let a = add(sub(x[t.j as usize], c), t.image_j);
            let b = add(sub(x[t.k as usize], c), t.image_k);
            let rl = add(sub(x[t.l as usize], c), t.image_l);
            let nrm = cross(a, b);
            let ln = fmax(norm(nrm), 1e-12);
            let ll = fmax(norm(rl), 1e-12);
            let nh = scale(nrm, 1.0 / ln);
            let el = scale(rl, 1.0 / ll);
            let s = dot(nh, el).clamp(-1.0, 1.0);
            let cw = fmax(1.0 - s * s, 1e-12).sqrt();
            e.inversion += t.force * (t.c0 + t.c1 * cw + t.c2 * (1.0 - 2.0 * s * s));
            if want_gradient {
                let de_ds = t.force * (t.c1 * (-s / cw) + t.c2 * (-4.0 * s));
                let gn = scale(sub(el, scale(nh, s)), 1.0 / ln);
                let ds_j = cross(b, gn);
                let ds_k = cross(gn, a);
                let ds_l = scale(sub(nh, scale(el, s)), 1.0 / ll);
                accumulate(&mut g, t.j, ds_j, de_ds);
                accumulate(&mut g, t.k, ds_k, de_ds);
                accumulate(&mut g, t.l, ds_l, de_ds);
                accumulate(&mut g, t.centre, add(add(ds_j, ds_k), ds_l), -de_ds);
            }
        }
        // Non-bonded pairs: rebuilt when an atom moved by more than a skin.
        let large = self.periodic || self.atoms > 1500;
        let mut rebuild = self.pair_positions.borrow().len() != x.len();
        if !rebuild && large {
            let previous = self.pair_positions.borrow();
            rebuild = (0..self.atoms as usize).any(|i| norm(sub(x[i], previous[i])) > 0.9);
        }
        if rebuild {
            self.build_pairs(x)?;
        }
        let cutoff = if large { self.options.cutoff } else { 1e300 };
        for pr in self.pairs.borrow().iter() {
            let d = add(sub(x[pr.j as usize], x[pr.i as usize]), pr.image);
            let r2 = dot(d, d);
            let r = fmax(r2, 1e-12).sqrt();
            if r > cutoff {
                continue;
            }
            let s2 = pr.x * pr.x / fmax(r2, 1e-12);
            let s6 = s2 * s2 * s2;
            let s12 = s6 * s6;
            e.vdw += pr.d * (s12 - 2.0 * s6);
            let mut de_dr = 12.0 * pr.d / r * (s6 - s12);
            if pr.q != 0.0 {
                let ec = pr.q / r;
                e.electrostatic += ec;
                de_dr -= ec / r;
            }
            if want_gradient {
                let f = scale(d, de_dr / r);
                accumulate(&mut g, pr.j, f, 1.0);
                accumulate(&mut g, pr.i, f, -1.0);
            }
        }
        for rs in &self.restraints {
            let a = rs.atoms;
            let p = |q: usize| x[a[q] as usize];
            match rs.kind {
                2 => {
                    let d = sub(p(1), p(0));
                    let r = fmax(norm(d), 1e-9);
                    let dr = r - rs.target;
                    e.restraint += rs.k * dr * dr;
                    let f = scale(d, 2.0 * rs.k * dr / r);
                    accumulate(&mut g, a[1], f, 1.0);
                    accumulate(&mut g, a[0], f, -1.0);
                }
                3 => {
                    let u = sub(p(0), p(1));
                    let v = sub(p(2), p(1));
                    let lu = fmax(norm(u), 1e-9);
                    let lv = fmax(norm(v), 1e-9);
                    let c = (dot(u, v) / (lu * lv)).clamp(-1.0, 1.0);
                    let theta = c.acos();
                    let d = theta - rs.target * PI / 180.0;
                    e.restraint += rs.k * d * d;
                    if want_gradient {
                        let sin_t = fmax((1.0 - c * c).sqrt(), 1e-6);
                        let de_dc = 2.0 * rs.k * d * (-1.0 / sin_t);
                        let du = sub(scale(v, 1.0 / (lu * lv)), scale(u, c / (lu * lu)));
                        let dv = sub(scale(u, 1.0 / (lu * lv)), scale(v, c / (lv * lv)));
                        accumulate(&mut g, a[0], du, de_dc);
                        accumulate(&mut g, a[2], dv, de_dc);
                        accumulate(&mut g, a[1], add(du, dv), -de_dc);
                    }
                }
                4 => {
                    let mut grad = [[0.0; 3]; 4];
                    let phi = dihedral(sub(p(1), p(0)), sub(p(2), p(1)), sub(p(3), p(2)), if want_gradient { Some(&mut grad) } else { None });
                    let mut d = phi - rs.target * PI / 180.0;
                    d = d.sin().atan2(d.cos());
                    e.restraint += rs.k * d * d;
                    if want_gradient {
                        for q in 0..4 {
                            accumulate(&mut g, a[q], grad[q], 2.0 * rs.k * d);
                        }
                    }
                }
                _ => {}
            }
        }
        if let Some(terms) = terms {
            *terms = e;
        }
        Ok(e.total())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MinimizeResult {
    pub initial_energy: f64,
    pub energy: f64,
    pub steps: i32,
    pub converged: bool,
    /// kcal/mol/A.
    pub max_force: f64,
}

/// L-BFGS minimisation over x. fixed[i] bits 0-2 hold x, y, z of atom i
/// (empty: nothing fixed). Converged when the largest force component is
/// below `tolerance`.
pub fn minimize(field: &UffForceField, x: &mut Vec<Vec3>, fixed: &[u8], max_steps: i32, tolerance: f64) -> Result<MinimizeResult, String> {
    let n = x.len();
    let dim = 3 * n;
    let mask = |g: &mut Vec<Vec3>| {
        if fixed.is_empty() {
            return;
        }
        for i in 0..n {
            for c in 0..3 {
                if (fixed[i] >> c) & 1 != 0 {
                    g[i][c] = 0.0;
                }
            }
        }
    };
    let flat = |v: &[Vec3]| -> Vec<f64> { v.iter().flat_map(|p| p.iter().copied()).collect() };
    let max_component = |g: &[Vec3]| {
        let mut m = 0.0f64;
        for v in g {
            for &c in v {
                m = fmax(m, c.abs());
            }
        }
        m
    };
    let mut g = Vec::new();
    let mut e = field.energy(x, Some(&mut g), None)?;
    mask(&mut g);
    let mut result = MinimizeResult { initial_energy: e, ..Default::default() };
    let memory = 8;
    let (mut s_hist, mut y_hist, mut rho_hist): (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>) = (Vec::new(), Vec::new(), Vec::new());
    let mut step = 0;
    while step < max_steps {
        result.max_force = max_component(&g);
        if result.max_force < tolerance {
            result.converged = true;
            break;
        }
        // Two-loop recursion for the search direction d = -H g.
        let mut q = flat(&g);
        let mut alpha = vec![0.0; s_hist.len()];
        for k in (0..s_hist.len()).rev() {
            let mut a = 0.0;
            for c in 0..dim {
                a += s_hist[k][c] * q[c];
            }
            a *= rho_hist[k];
            alpha[k] = a;
            for c in 0..dim {
                q[c] -= a * y_hist[k][c];
            }
        }
        let mut gamma = 0.01;
        if let (Some(s), Some(y)) = (s_hist.last(), y_hist.last()) {
            let (mut sy, mut yy) = (0.0, 0.0);
            for c in 0..dim {
                sy += s[c] * y[c];
                yy += y[c] * y[c];
            }
            gamma = if yy > 0.0 { sy / yy } else { 0.01 };
        }
        for v in q.iter_mut() {
            *v *= gamma;
        }
        for k in 0..s_hist.len() {
            let mut b = 0.0;
            for c in 0..dim {
                b += y_hist[k][c] * q[c];
            }
            b *= rho_hist[k];
            for c in 0..dim {
                q[c] += s_hist[k][c] * (alpha[k] - b);
            }
        }
        let mut d = vec![0.0; dim];
        let mut slope = 0.0;
        let gf = flat(&g);
        for c in 0..dim {
            d[c] = -q[c];
            slope += d[c] * gf[c];
        }
        if slope >= 0.0 {
            // Not a descent direction: restart along the gradient.
            s_hist.clear();
            y_hist.clear();
            rho_hist.clear();
            slope = 0.0;
            for c in 0..dim {
                d[c] = -gf[c] * 0.01;
                slope += d[c] * gf[c];
            }
        }
        if !fixed.is_empty() {
            for i in 0..n {
                for c in 0..3 {
                    if (fixed[i] >> c) & 1 != 0 {
                        d[3 * i + c] = 0.0;
                    }
                }
            }
        }
        // At most 0.3 A per atom per step.
        let mut longest = 0.0f64;
        for i in 0..n {
            longest = fmax(longest, (d[3 * i] * d[3 * i] + d[3 * i + 1] * d[3 * i + 1] + d[3 * i + 2] * d[3 * i + 2]).sqrt());
        }
        let mut t = if longest > 0.3 { 0.3 / longest } else { 1.0 };
        let mut trial = vec![[0.0; 3]; n];
        let mut g_trial = Vec::new();
        let mut e_trial = e;
        let mut accepted = false;
        for _ in 0..30 {
            for i in 0..n {
                for c in 0..3 {
                    trial[i][c] = x[i][c] + t * d[3 * i + c];
                }
            }
            e_trial = field.energy(&trial, Some(&mut g_trial), None)?;
            if e_trial <= e + 1e-4 * t * slope {
                accepted = true;
                break;
            }
            t *= 0.5;
        }
        if !accepted {
            if s_hist.is_empty() {
                break; // no progress along the gradient either
            }
            s_hist.clear();
            y_hist.clear();
            rho_hist.clear();
            step += 1;
            continue;
        }
        mask(&mut g_trial);
        let mut s = vec![0.0; dim];
        let mut y = vec![0.0; dim];
        let mut sy = 0.0;
        for i in 0..n {
            for c in 0..3 {
                s[3 * i + c] = trial[i][c] - x[i][c];
                y[3 * i + c] = g_trial[i][c] - g[i][c];
                sy += s[3 * i + c] * y[3 * i + c];
            }
        }
        *x = trial;
        let previous = e;
        e = e_trial;
        g = g_trial;
        if sy > 1e-10 {
            s_hist.push(s);
            y_hist.push(y);
            rho_hist.push(1.0 / sy);
            if s_hist.len() > memory {
                s_hist.remove(0);
                y_hist.remove(0);
                rho_hist.remove(0);
            }
        }
        if (previous - e).abs() < 1e-10 && max_component(&g) < 10.0 * tolerance {
            result.converged = true;
            step += 1;
            break;
        }
        step += 1;
    }
    result.max_force = max_component(&g);
    if result.max_force < tolerance {
        result.converged = true;
    }
    result.energy = e;
    result.steps = step;
    Ok(result)
}
