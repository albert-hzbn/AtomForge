//! Small molecules with experimental gas-phase geometries (NIST CCCBDB bond
//! lengths and angles), centred on the centre of mass (the C++
//! pipeline/MoleculeLibrary). Where GCC folds a sin or cos of a constant at
//! compile time (correctly rounded: pyramid angles, methanol, ethylene, H2CO)
//! af_core::crmath is used; bent, ethane, benzene and pyramid's cos(angle)
//! call the runtime UCRT functions in the C++ build, as f64::sin/cos do.

use std::collections::BTreeMap;
use std::path::Path;

use af_core::crmath;
use af_core::cxx::sort_by;
use af_core::elements::atomic_mass;
use af_science::core::{add, cross, dot, norm, scale, sub, Vec3};
use af_science::data::{atomic_number, read_frames};

const PI: f64 = 3.14159265358979323846;
const DEGREE: f64 = PI / 180.0;

#[derive(Clone, Debug, PartialEq)]
pub struct MoleculeAtom {
    pub symbol: String,
    pub position: Vec3,
}

type Atoms = Vec<MoleculeAtom>;

fn a(symbol: &str, position: Vec3) -> MoleculeAtom {
    MoleculeAtom { symbol: symbol.into(), position }
}

fn tetrahedron() -> [Vec3; 4] {
    let s = 1.0 / 3.0f64.sqrt();
    [[s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]]
}

fn diatomic(x: &str, y: &str, length: f64) -> Atoms {
    vec![a(x, [0.0, 0.0, 0.0]), a(y, [0.0, 0.0, length])]
}

fn bent(apex: &str, arm: &str, length: f64, angle_degrees: f64) -> Atoms {
    let half = 0.5 * angle_degrees * DEGREE;
    vec![
        a(apex, [0.0, 0.0, 0.0]),
        a(arm, [length * half.sin(), 0.0, length * half.cos()]),
        a(arm, [-length * half.sin(), 0.0, length * half.cos()]),
    ]
}

fn pyramid(apex: &str, arm: &str, length: f64, angle_degrees: f64) -> Atoms {
    let c = (angle_degrees * DEGREE).cos();
    let sin_beta = (2.0 * (1.0 - c) / 3.0).sqrt();
    let cos_beta = (1.0 - sin_beta * sin_beta).sqrt();
    let mut out = vec![a(apex, [0.0, 0.0, 0.0])];
    for k in 0..3 {
        let phi = 2.0 * PI * k as f64 / 3.0;
        out.push(a(arm, [length * sin_beta * crmath::cos(phi), length * sin_beta * crmath::sin(phi), -length * cos_beta]));
    }
    out
}

fn tetrahedral_hydride(centre: &str, arm: &str, length: f64) -> Atoms {
    let mut out = vec![a(centre, [0.0, 0.0, 0.0])];
    for d in tetrahedron() {
        out.push(a(arm, scale(d, length)));
    }
    out
}

fn ethane() -> Atoms {
    let (cc, ch, hcc) = (1.535, 1.094, 111.2 * DEGREE);
    let mut out = vec![a("C", [0.0, 0.0, -cc / 2.0]), a("C", [0.0, 0.0, cc / 2.0])];
    for side in 0..2 {
        let sign = if side == 0 { -1.0 } else { 1.0 };
        for k in 0..3 {
            let phi = 2.0 * PI * k as f64 / 3.0 + if side != 0 { PI / 3.0 } else { 0.0 };
            let z = sign * (cc / 2.0 - ch * hcc.cos());
            out.push(a("H", [ch * hcc.sin() * phi.cos(), ch * hcc.sin() * phi.sin(), z]));
        }
    }
    out
}

fn ethylene() -> Atoms {
    let (cc, ch, hcc) = (1.339, 1.086, 121.3 * DEGREE);
    let mut out = vec![a("C", [0.0, 0.0, -cc / 2.0]), a("C", [0.0, 0.0, cc / 2.0])];
    for side in 0..2 {
        let sign = if side == 0 { -1.0 } else { 1.0 };
        for k in 0..2 {
            let x = (if k != 0 { 1.0 } else { -1.0 }) * ch * crmath::sin(hcc);
            out.push(a("H", [x, 0.0, sign * (cc / 2.0 - ch * crmath::cos(hcc))]));
        }
    }
    out
}

fn acetylene() -> Atoms {
    let (cc, ch) = (1.203, 1.063);
    vec![a("C", [0.0, 0.0, -cc / 2.0]), a("C", [0.0, 0.0, cc / 2.0]), a("H", [0.0, 0.0, -cc / 2.0 - ch]), a("H", [0.0, 0.0, cc / 2.0 + ch])]
}

fn benzene() -> Atoms {
    let (cc, ch) = (1.397, 1.084);
    let mut out = Vec::new();
    for k in 0..6 {
        let phi = k as f64 * 60.0 * DEGREE;
        out.push(a("C", [cc * phi.cos(), cc * phi.sin(), 0.0]));
    }
    for k in 0..6 {
        let phi = k as f64 * 60.0 * DEGREE;
        out.push(a("H", [(cc + ch) * phi.cos(), (cc + ch) * phi.sin(), 0.0]));
    }
    out
}

fn methanol() -> Atoms {
    let (co, oh, ch, coh, hco) = (1.427, 0.956, 1.096, 108.5 * DEGREE, 109.5 * DEGREE);
    let mut out = vec![a("C", [0.0, 0.0, 0.0]), a("O", [0.0, 0.0, co]), a("H", [oh * crmath::sin(coh), 0.0, co - oh * crmath::cos(coh)])];
    for k in 0..3 {
        let phi = PI + 2.0 * PI * k as f64 / 3.0;
        out.push(a("H", [ch * crmath::sin(hco) * crmath::cos(phi), ch * crmath::sin(hco) * crmath::sin(phi), ch * crmath::cos(hco)]));
    }
    out
}

fn fullerene() -> Atoms {
    let p = (1.0 + 5.0f64.sqrt()) / 2.0;
    let seeds: [Vec3; 3] = [[0.0, 1.0, 3.0 * p], [1.0, 2.0 + p, 2.0 * p], [p, 2.0, p * p * p]];
    let mut points: Vec<Vec3> = Vec::new();
    for seed in seeds {
        for sx in [-1.0, 1.0] {
            for sy in [-1.0, 1.0] {
                for sz in [-1.0, 1.0] {
                    let v = [seed[0] * sx, seed[1] * sy, seed[2] * sz];
                    for q in [v, [v[1], v[2], v[0]], [v[2], v[0], v[1]]] {
                        if !points.iter().any(|r| norm(sub(q, *r)) < 1e-6) {
                            points.push(q);
                        }
                    }
                }
            }
        }
    }
    points.into_iter().map(|q| a("C", scale(q, 1.43 / 2.0))).collect()
}

fn h2co() -> Atoms {
    let (co, ch, half) = (1.205, 1.111, 0.5 * 116.1 * DEGREE);
    vec![
        a("C", [0.0, 0.0, 0.0]),
        a("O", [0.0, 0.0, co]),
        a("H", [ch * crmath::sin(half), 0.0, -ch * crmath::cos(half)]),
        a("H", [-ch * crmath::sin(half), 0.0, -ch * crmath::cos(half)]),
    ]
}

/// The library in std::map (byte-wise) order.
fn table() -> BTreeMap<&'static str, fn() -> Atoms> {
    let entries: [(&str, fn() -> Atoms); 28] = [
        ("H2", || diatomic("H", "H", 0.741)),
        ("N2", || diatomic("N", "N", 1.098)),
        ("O2", || diatomic("O", "O", 1.208)),
        ("F2", || diatomic("F", "F", 1.412)),
        ("Cl2", || diatomic("Cl", "Cl", 1.988)),
        ("CO", || diatomic("C", "O", 1.128)),
        ("NO", || diatomic("N", "O", 1.151)),
        ("HF", || diatomic("H", "F", 0.917)),
        ("HCl", || diatomic("H", "Cl", 1.275)),
        ("H2O", || bent("O", "H", 0.9572, 104.52)),
        ("H2S", || bent("S", "H", 1.336, 92.11)),
        ("CO2", || vec![a("C", [0.0, 0.0, 0.0]), a("O", [0.0, 0.0, 1.160]), a("O", [0.0, 0.0, -1.160])]),
        ("SO2", || bent("S", "O", 1.431, 119.3)),
        ("NO2", || bent("N", "O", 1.193, 134.1)),
        ("O3", || bent("O", "O", 1.278, 116.8)),
        ("NH3", || pyramid("N", "H", 1.012, 106.67)),
        ("PH3", || pyramid("P", "H", 1.420, 93.3)),
        ("CH4", || tetrahedral_hydride("C", "H", 1.087)),
        ("SiH4", || tetrahedral_hydride("Si", "H", 1.480)),
        ("CF4", || tetrahedral_hydride("C", "F", 1.315)),
        ("C2H2", acetylene),
        ("C2H4", ethylene),
        ("C2H6", ethane),
        ("CH3OH", methanol),
        ("C6H6", benzene),
        ("HCN", || vec![a("H", [0.0, 0.0, -1.064]), a("C", [0.0, 0.0, 0.0]), a("N", [0.0, 0.0, 1.156])]),
        ("H2CO", h2co),
        ("C60", fullerene),
    ];
    entries.into_iter().collect()
}

fn aliases() -> BTreeMap<&'static str, &'static str> {
    [
        ("water", "H2O"),
        ("ammonia", "NH3"),
        ("methane", "CH4"),
        ("silane", "SiH4"),
        ("ethane", "C2H6"),
        ("ethylene", "C2H4"),
        ("ethene", "C2H4"),
        ("acetylene", "C2H2"),
        ("ethyne", "C2H2"),
        ("benzene", "C6H6"),
        ("methanol", "CH3OH"),
        ("fullerene", "C60"),
        ("formaldehyde", "H2CO"),
        ("hydrogen cyanide", "HCN"),
        ("buckminsterfullerene", "C60"),
        ("ozone", "O3"),
        ("hydrogen", "H2"),
        ("nitrogen", "N2"),
        ("oxygen", "O2"),
        ("carbon dioxide", "CO2"),
        ("carbon monoxide", "CO"),
    ]
    .into_iter()
    .collect()
}

fn centre(atoms: &mut Atoms) {
    let mut total = 0.0;
    let mut c = [0.0; 3];
    for x in atoms.iter() {
        let m = atomic_mass(atomic_number(&x.symbol));
        total += m;
        c = add(c, scale(x.position, m));
    }
    c = scale(c, 1.0 / total);
    for x in atoms.iter_mut() {
        x.position = sub(x.position, c);
    }
}

/// The names of the library's molecules.
pub fn molecule_names() -> Vec<String> {
    table().keys().map(|k| k.to_string()).collect()
}

/// A molecule by name (case-insensitive), a structure file or a single element.
pub fn molecule(name: &str) -> Result<Vec<MoleculeAtom>, String> {
    let key = name.to_ascii_lowercase();
    for (entry, make) in table() {
        if entry.to_ascii_lowercase() == key {
            let mut atoms = make();
            centre(&mut atoms);
            return Ok(atoms);
        }
    }
    if let Some(alias) = aliases().get(key.as_str()) {
        return molecule(alias);
    }
    if name.contains(['.', '/', '\\']) && Path::new(name).exists() {
        let frames = read_frames(Path::new(name), None)?;
        let mut atoms: Atoms = frames.first().map(|f| f.structure.atoms.iter().map(|x| a(&x.symbol, [x.x, x.y, x.z])).collect()).unwrap_or_default();
        if atoms.is_empty() {
            return Err(format!("The molecule file {name} has no atoms"));
        }
        centre(&mut atoms);
        return Ok(atoms);
    }
    if !name.is_empty() {
        let mut symbol: String = name.chars().take(1).map(|c| c.to_ascii_uppercase()).collect();
        symbol.extend(name.chars().skip(1).map(|c| c.to_ascii_lowercase()));
        if atomic_number(&symbol) > 0 {
            return Ok(vec![a(&symbol, [0.0; 3])]);
        }
    }
    Err(format!("Unknown molecule '{name}'; give an element symbol or one of {}", molecule_names().join(", ")))
}

/// The mass of a molecule (u).
pub fn molecule_mass(atoms: &[MoleculeAtom]) -> f64 {
    let mut mass = 0.0;
    for x in atoms {
        mass += atomic_mass(atomic_number(&x.symbol));
    }
    mass
}

/// Principal axes of inertia (unit vectors, smallest moment first, right-handed).
pub fn principal_axes(atoms: &[MoleculeAtom]) -> [Vec3; 3] {
    let mut total = 0.0;
    let mut c = [0.0; 3];
    for x in atoms {
        let m = atomic_mass(atomic_number(&x.symbol));
        total += m;
        c = add(c, scale(x.position, m));
    }
    c = scale(c, 1.0 / total);
    let mut t = [[0.0f64; 3]; 3];
    for x in atoms {
        let m = atomic_mass(atomic_number(&x.symbol));
        let r = sub(x.position, c);
        let r2 = dot(r, r);
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] += m * ((if i == j { r2 } else { 0.0 }) - r[i] * r[j]);
            }
        }
    }
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let off = t[0][1].abs() + t[0][2].abs() + t[1][2].abs();
        if off < 1e-14 {
            break;
        }
        for p in 0..2 {
            for q in p + 1..3 {
                if t[p][q].abs() < 1e-18 {
                    continue;
                }
                let theta = 0.5 * (2.0 * t[p][q]).atan2(t[q][q] - t[p][p]);
                let (cs, sn) = (theta.cos(), theta.sin());
                for k in 0..3 {
                    let (tkp, tkq) = (t[k][p], t[k][q]);
                    t[k][p] = cs * tkp - sn * tkq;
                    t[k][q] = sn * tkp + cs * tkq;
                }
                for k in 0..3 {
                    let (tpk, tqk) = (t[p][k], t[q][k]);
                    t[p][k] = cs * tpk - sn * tqk;
                    t[q][k] = sn * tpk + cs * tqk;
                }
                for k in 0..3 {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = cs * vkp - sn * vkq;
                    v[k][q] = sn * vkp + cs * vkq;
                }
            }
        }
    }
    let mut order = [0usize, 1, 2];
    sort_by(&mut order, |&x, &y| t[x][x] < t[y][y]);
    let mut axes = [[0.0; 3]; 3];
    for k in 0..3 {
        axes[k] = [v[0][order[k]], v[1][order[k]], v[2][order[k]]];
    }
    if dot(cross(axes[0], axes[1]), axes[2]) < 0.0 {
        axes[2] = scale(axes[2], -1.0);
    }
    axes
}
