//! Building blocks shared by the builder and edit steps (the C++
//! pipeline/StructureEdits): covalent bonds, rotations, appending atoms and
//! saturating dangling bonds.

use std::collections::{BTreeMap, BTreeSet};

use af_core::{AtomSite, Structure};
use af_science::core::{add, cross, dot, identity, neighbor_list, norm, scale, sub, Mat3, Vec3};

const TETRAHEDRAL: f64 = 109.4712206;
const DEGREE: f64 = 3.14159265358979323846 / 180.0;
const PI: f64 = 3.14159265358979323846;

/// A covalent bond i-j: the vector from i to the bonded image of j.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub vector: Vec3,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Neighbours {
    pub atoms: Vec<Vec<usize>>,
    pub vectors: Vec<Vec<Vec3>>,
}

fn unit(v: Vec3) -> Vec3 {
    let n = norm(v);
    if n > 1e-12 {
        scale(v, 1.0 / n)
    } else {
        [0.0; 3]
    }
}

fn perpendicular(v: Vec3) -> Vec3 {
    let trial = if v[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    unit(cross(v, trial))
}

/// The rotation taking direction `from` to `to` (shortest arc); v' = R v.
pub fn rotation_between(from: Vec3, to: Vec3) -> Mat3 {
    let (a, b) = (unit(from), unit(to));
    let v = cross(a, b);
    let c = dot(a, b);
    if c < -1.0 + 1e-12 {
        return rotation_about(perpendicular(a), PI);
    }
    if norm(v) < 1e-14 {
        return identity();
    }
    rotation_about(v, norm(v).atan2(c))
}

/// The rotation by `angle` (radians, right-handed) about `axis`.
pub fn rotation_about(axis: Vec3, angle: f64) -> Mat3 {
    let k = unit(axis);
    let (c, s) = (angle.cos(), angle.sin());
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = (if i == j { c } else { 0.0 }) + (1.0 - c) * k[i] * k[j];
        }
    }
    r[0][1] -= s * k[2];
    r[0][2] += s * k[1];
    r[1][0] += s * k[2];
    r[1][2] -= s * k[0];
    r[2][0] -= s * k[1];
    r[2][1] += s * k[0];
    r
}

pub fn rotate(r: &Mat3, v: Vec3) -> Vec3 {
    [dot(r[0], v), dot(r[1], v), dot(r[2], v)]
}

/// Covalent radius (A, Cordero et al. 2008); 1.5 when unknown.
pub fn covalent_radius(z: i32) -> f64 {
    let radii = af_core::elements::literature_covalent_radii();
    if z > 0 && (z as usize) < radii.len() && radii[z as usize] > 0.0 {
        return radii[z as usize] as f64;
    }
    1.5
}

/// Pairs closer than scale x (sum of covalent radii), across periodic cell
/// vectors; each pair once (i < j, or i == j for a bond to an own image).
pub fn covalent_bonds(s: &Structure, bond_scale: f64) -> Result<Vec<Bond>, String> {
    let mut bonds = Vec::new();
    let n = s.atoms.len();
    if n < 2 && !(n == 1 && s.has_unit_cell) {
        return Ok(bonds);
    }
    let mut largest = 0.0f64;
    let positions: Vec<Vec3> = s
        .atoms
        .iter()
        .map(|a| {
            let r = covalent_radius(a.atomic_number);
            if largest < r {
                largest = r;
            }
            [a.x, a.y, a.z]
        })
        .collect();
    let (mut pbc, mut cell) = ([false; 3], identity());
    if s.has_unit_cell {
        cell = s.cell_vectors;
        pbc = s.pbc;
    }
    for link in neighbor_list(&positions, &cell, pbc, bond_scale * 2.0 * largest)? {
        let (i, j) = (link.i as usize, link.j as usize);
        if j < i {
            continue;
        }
        if i == j {
            let t = link.shift;
            if t[0] < 0 || (t[0] == 0 && (t[1] < 0 || (t[1] == 0 && t[2] < 0))) {
                continue;
            }
        }
        let limit = bond_scale * (covalent_radius(s.atoms[i].atomic_number) + covalent_radius(s.atoms[j].atomic_number));
        if norm(link.vector) <= limit && norm(link.vector) > 1e-6 {
            bonds.push(Bond { i, j, vector: link.vector });
        }
    }
    Ok(bonds)
}

pub fn bond_neighbours(s: &Structure, bond_scale: f64) -> Result<Neighbours, String> {
    let n = s.atoms.len();
    let mut out = Neighbours { atoms: vec![Vec::new(); n], vectors: vec![Vec::new(); n] };
    for b in covalent_bonds(s, bond_scale)? {
        out.atoms[b.i].push(b.j);
        out.vectors[b.i].push(b.vector);
        out.atoms[b.j].push(b.i);
        out.vectors[b.j].push(scale(b.vector, -1.0));
    }
    Ok(out)
}

/// A new atom (default colour) at the end, extending the aligned per-atom
/// metadata and the selection.
pub fn append_atom(s: &mut Structure, selected: &mut Vec<u8>, symbol: &str, position: Vec3, select: bool) {
    let z = af_science::data::atomic_number(symbol);
    let [r, g, b] = af_core::elements::default_colour(z);
    let atom = AtomSite { symbol: symbol.to_string(), atomic_number: z, x: position[0], y: position[1], z: position[2], r, g, b, ..Default::default() };
    let n = s.atoms.len();
    if s.grain_colors.len() == n && n > 0 {
        s.grain_colors.push([0.6, 0.6, 0.6]);
    }
    if s.grain_region_ids.len() == n && n > 0 {
        s.grain_region_ids.push(-1);
    }
    if s.atom_property.len() == n && n > 0 {
        s.atom_property.push(f64::NAN);
    }
    if s.frozen_property.len() == n && n > 0 {
        s.frozen_property.push(f64::NAN);
    }
    s.atoms.push(atom);
    selected.resize(n, 0);
    selected.push(select as u8);
}

/// Saturation with a terminating atom: an explicit bond length (0: H from
/// its table, other atoms from covalent radii), at most `max_per_atom` new
/// atoms per atom, and only atoms with mask[i] (when given).
#[derive(Clone, Debug, PartialEq)]
pub struct SaturationOptions<'a> {
    pub terminator: String,
    pub bond_length: f64,
    pub max_per_atom: i32,
    pub mask: Option<&'a [u8]>,
}

impl Default for SaturationOptions<'_> {
    fn default() -> Self {
        Self { terminator: "H".into(), bond_length: 0.0, max_per_atom: 4, mask: None }
    }
}

pub fn saturate_with_hydrogen(s: &mut Structure, selected: &mut Vec<u8>, elements: &BTreeSet<String>, force_sp3: bool, bond_scale: f64, select_new: bool) -> Result<usize, String> {
    saturate(s, selected, elements, force_sp3, bond_scale, select_new, &SaturationOptions::default())
}

pub fn saturate(
    s: &mut Structure,
    selected: &mut Vec<u8>,
    elements: &BTreeSet<String>,
    force_sp3: bool,
    bond_scale: f64,
    select_new: bool,
    options: &SaturationOptions,
) -> Result<usize, String> {
    let valence: BTreeMap<&str, i32> =
        [("C", 4), ("Si", 4), ("Ge", 4), ("Sn", 4), ("N", 3), ("P", 3), ("B", 3), ("Al", 3), ("O", 2), ("S", 2), ("Se", 2)].into_iter().collect();
    let hydrogen_bond: BTreeMap<&str, f64> = [
        ("C", 1.09),
        ("Si", 1.48),
        ("Ge", 1.53),
        ("Sn", 1.71),
        ("N", 1.01),
        ("P", 1.42),
        ("B", 1.19),
        ("Al", 1.58),
        ("O", 0.97),
        ("S", 1.34),
        ("Se", 1.46),
    ]
    .into_iter()
    .collect();
    let nb = bond_neighbours(s, bond_scale)?;
    let mut additions: Vec<Vec3> = Vec::new();
    for i in 0..s.atoms.len() {
        let symbol = s.atoms[i].symbol.as_str();
        if !elements.contains(symbol) {
            continue;
        }
        if let Some(mask) = options.mask {
            if i >= mask.len() || mask[i] == 0 {
                continue;
            }
        }
        let Some(&val) = valence.get(symbol) else { continue };
        let u: Vec<Vec3> = nb.vectors[i].iter().map(|&d| unit(d)).collect();
        let k = u.len() as i32;
        let mut missing = val - k;
        let planar_capable = !force_sp3 && (symbol == "C" || symbol == "B" || symbol == "N");
        if planar_capable && k == 3 && norm(add(add(u[0], u[1]), u[2])) < 0.5 {
            missing = 0;
        }
        let sp2_pair = planar_capable && k == 2 && dot(u[0], u[1]) < (115.0 * DEGREE).cos();
        if sp2_pair {
            missing = missing.min(1);
        }
        missing = missing.min(options.max_per_atom);
        if missing <= 0 {
            continue;
        }
        let terminator_z = af_science::data::atomic_number(&options.terminator);
        let bond = if options.bond_length > 0.0 {
            options.bond_length
        } else if options.terminator == "H" && hydrogen_bond.contains_key(symbol) {
            hydrogen_bond[symbol]
        } else {
            covalent_radius(s.atoms[i].atomic_number) + covalent_radius(terminator_z)
        };
        let origin = [s.atoms[i].x, s.atoms[i].y, s.atoms[i].z];
        let mut directions: Vec<Vec3> = Vec::new();
        if k == 0 {
            let c = 1.0 / 3.0f64.sqrt();
            let tetra = [[c, c, c], [c, -c, -c], [-c, c, -c], [-c, -c, c]];
            if val == 2 {
                let half = 0.5 * 104.5 * DEGREE;
                directions = vec![[half.sin(), 0.0, half.cos()], [-half.sin(), 0.0, half.cos()]];
            } else {
                directions.extend_from_slice(&tetra[..4.min(missing) as usize]);
            }
        } else if k == 1 {
            let alpha = (if val == 2 { 104.5 } else { TETRAHEDRAL }) * DEGREE;
            let e1 = perpendicular(u[0]);
            let e2 = unit(cross(u[0], e1));
            for m in 0..missing {
                let phi = 2.0 * PI * m as f64 / 3.0;
                directions.push(add(scale(u[0], alpha.cos()), scale(add(scale(e1, phi.cos()), scale(e2, phi.sin())), alpha.sin())));
            }
        } else if k == 2 {
            let b = unit(scale(add(u[0], u[1]), -1.0));
            let p = unit(cross(u[0], u[1]));
            let half = 0.5 * TETRAHEDRAL * DEGREE;
            if sp2_pair || missing == 1 {
                if sp2_pair {
                    directions.push(b);
                } else {
                    directions.push(unit(add(scale(b, half.cos()), scale(p, half.sin()))));
                }
            } else {
                directions.push(unit(add(scale(b, half.cos()), scale(p, half.sin()))));
                directions.push(unit(sub(scale(b, half.cos()), scale(p, half.sin()))));
            }
        } else {
            let mut d = scale(add(add(u[0], u[1]), u[2]), -1.0);
            if norm(d) < 1e-3 {
                d = cross(sub(u[1], u[0]), sub(u[2], u[0]));
            }
            directions.push(unit(d));
        }
        for d in directions {
            additions.push(add(origin, scale(unit(d), bond)));
        }
    }
    let count = additions.len();
    for position in additions {
        append_atom(s, selected, &options.terminator, position, select_new);
    }
    Ok(count)
}
