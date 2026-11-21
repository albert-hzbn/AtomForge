//! Interfaces between two crystals by geometric lattice matching (the C++
//! algorithms/InterfaceBuilder and InterfaceMatching), after Stradi et al.,
//! J. Phys.: Condens. Matter 29, 185901 (2017): each crystal is cut along a
//! plane, B's surface lattice is rotated through a range of angles, the
//! supercells of the two surface lattices that coincide within a strain limit
//! are listed, and a chosen match is assembled into one cell.

use std::collections::BTreeMap;

use af_core::cxx::sort_by;
use af_core::{AtomSite, Structure};

use crate::surface::{build_surface, SurfaceParams};

const PI: f64 = 3.14159265358979323846;

// ---------------------------------------------------------------- glm dvec2 / dmat2

type V2 = [f64; 2];
/// glm::dmat2: m[column][row].
type M2 = [[f64; 2]; 2];

fn add(a: V2, b: V2) -> V2 {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub(a: V2, b: V2) -> V2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn scale(s: f64, v: V2) -> V2 {
    [s * v[0], s * v[1]]
}
fn dot2(a: V2, b: V2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn length2(v: V2) -> f64 {
    dot2(v, v).sqrt()
}
fn cross2(a: V2, b: V2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn det(m: &M2) -> f64 {
    m[0][0] * m[1][1] - m[1][0] * m[0][1]
}
fn inverse(m: &M2) -> M2 {
    let o = 1.0 / (m[0][0] * m[1][1] - m[1][0] * m[0][1]);
    [[m[1][1] * o, -m[0][1] * o], [-m[1][0] * o, m[0][0] * o]]
}
fn transpose(m: &M2) -> M2 {
    [[m[0][0], m[1][0]], [m[0][1], m[1][1]]]
}
fn mul(a: &M2, b: &M2) -> M2 {
    [
        [a[0][0] * b[0][0] + a[1][0] * b[0][1], a[0][1] * b[0][0] + a[1][1] * b[0][1]],
        [a[0][0] * b[1][0] + a[1][0] * b[1][1], a[0][1] * b[1][0] + a[1][1] * b[1][1]],
    ]
}
fn mulv(m: &M2, v: V2) -> V2 {
    [m[0][0] * v[0] + m[1][0] * v[1], m[0][1] * v[0] + m[1][1] * v[1]]
}
fn madd(a: &M2, b: &M2) -> M2 {
    [[a[0][0] + b[0][0], a[0][1] + b[0][1]], [a[1][0] + b[1][0], a[1][1] + b[1][1]]]
}
fn msub(a: &M2, b: &M2) -> M2 {
    [[a[0][0] - b[0][0], a[0][1] - b[0][1]], [a[1][0] - b[1][0], a[1][1] - b[1][1]]]
}
fn mscale(s: f64, m: &M2) -> M2 {
    [[s * m[0][0], s * m[0][1]], [s * m[1][0], s * m[1][1]]]
}
fn mdiv(m: &M2, s: f64) -> M2 {
    [[m[0][0] / s, m[0][1] / s], [m[1][0] / s, m[1][1] / s]]
}
const I2: M2 = [[1.0, 0.0], [0.0, 1.0]];

// ---------------------------------------------------------------- InterfaceBuilder.cpp

/// The structure stacked `layers` times along its third cell vector.
pub fn repeat_layers_z(base: &Structure, layers: i32) -> Structure {
    if layers <= 1 {
        return base.clone();
    }
    let mut result = base.clone();
    let cz = base.cell_vectors[2];
    for l in 1..layers {
        for a in &base.atoms {
            let mut a = a.clone();
            a.x += l as f64 * cz[0];
            a.y += l as f64 * cz[1];
            a.z += l as f64 * cz[2];
            result.atoms.push(a);
        }
    }
    for v in result.cell_vectors[2].iter_mut() {
        *v *= layers as f64;
    }
    result
}

/// Atom positions and the in-plane cell vectors mapped by the 2x2 matrix F
/// (rows: x' = F[0][0] x + F[0][1] y, y' = F[1][0] x + F[1][1] y).
pub fn apply_transform_2d(b: &Structure, f: &[[f64; 2]; 2]) -> Structure {
    let mut result = b.clone();
    for atom in result.atoms.iter_mut() {
        let (ox, oy) = (atom.x, atom.y);
        atom.x = f[0][0] * ox + f[0][1] * oy;
        atom.y = f[1][0] * ox + f[1][1] * oy;
    }
    if result.has_unit_cell {
        for i in 0..2 {
            let (ox, oy) = (result.cell_vectors[i][0], result.cell_vectors[i][1]);
            result.cell_vectors[i][0] = f[0][0] * ox + f[0][1] * oy;
            result.cell_vectors[i][1] = f[1][0] * ox + f[1][1] * oy;
        }
    }
    result
}

/// The interface repeated rx times along a and ry times along b.
pub fn repeat_interface_xy(iface: &Structure, rx: i32, ry: i32) -> Structure {
    if rx <= 1 && ry <= 1 {
        return iface.clone();
    }
    let mut result = Structure { has_unit_cell: iface.has_unit_cell, cell_vectors: iface.cell_vectors, cell_offset: iface.cell_offset, ..Default::default() };
    let (ax, ay) = (iface.cell_vectors[0][0], iface.cell_vectors[0][1]);
    let (bx, by) = (iface.cell_vectors[1][0], iface.cell_vectors[1][1]);
    for t in 0..rx * ry {
        let (ix, iy) = (t / ry, t % ry);
        let dx = ix as f64 * ax + iy as f64 * bx;
        let dy = ix as f64 * ay + iy as f64 * by;
        for a in &iface.atoms {
            let mut a = a.clone();
            a.x += dx;
            a.y += dy;
            result.atoms.push(a);
        }
    }
    result.cell_vectors[0][0] *= rx as f64;
    result.cell_vectors[0][1] *= rx as f64;
    result.cell_vectors[1][0] *= ry as f64;
    result.cell_vectors[1][1] *= ry as f64;
    result
}

// ---------------------------------------------------------------- InterfaceMatching.cpp

/// One side of the interface: a crystal cut along (h k l), or a structure
/// that is already a slab (its a and b vectors span the interface plane).
#[derive(Clone, Debug, PartialEq)]
pub struct Side {
    pub crystal: Structure,
    pub cleave: bool,
    pub h: i32,
    pub k: i32,
    pub l: i32,
    /// Surface layers (repeats of the stacking translation) in the slab.
    pub layers: i32,
}

impl Default for Side {
    fn default() -> Self {
        Self { crystal: Structure::default(), cleave: true, h: 0, k: 0, l: 1, layers: 4 }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SearchParams {
    pub max_repeats: i32,
    pub max_strain: f64,
    pub max_area: f64,
    pub min_cell_angle: f64,
    pub angle_min: f64,
    pub angle_max: f64,
    pub angle_step: f64,
    pub max_matches: usize,
}

impl Default for SearchParams {
    fn default() -> Self {
        Self { max_repeats: 6, max_strain: 0.03, max_area: 150.0, min_cell_angle: 20.0, angle_min: 0.0, angle_max: 180.0, angle_step: 1.0, max_matches: 2000 }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Match {
    /// Rows: the interface cell in A's surface vectors.
    pub a: [[i32; 2]; 2],
    /// Rows: the interface cell in B's surface vectors.
    pub b: [[i32; 2]; 2],
    pub theta_deg: f64,
    /// e11, e22, e12 taking B's cell onto A's.
    pub strain: [f64; 3],
    pub mean_strain: f64,
    pub area: f64,
    pub cell_angle: f64,
    pub cells_a: i32,
    pub cells_b: i32,
    pub atoms: i32,
    pub pareto: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StrainOn {
    #[default]
    B,
    A,
    Both,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Assembly {
    pub gap: f64,
    pub vacuum: f64,
    pub shift: [f64; 2],
    pub strain_on: StrainOn,
}

impl Default for Assembly {
    fn default() -> Self {
        Self { gap: 2.0, vacuum: 0.0, shift: [0.0, 0.0], strain_on: StrainOn::B }
    }
}

type D3 = [f64; 3];
fn dot3(a: D3, b: D3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross3(a: D3, b: D3) -> D3 {
    [a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0], a[0] * b[1] - b[0] * a[1]]
}
fn normalize3(v: D3) -> D3 {
    let s = 1.0 / dot3(v, v).sqrt();
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Rotates the structure so a lies along x and a x b along z.
fn orient(s: &mut Structure) {
    let a = s.cell_vectors[0];
    let b = s.cell_vectors[1];
    let x = normalize3(a);
    let z = normalize3(cross3(a, b));
    let y = cross3(z, x);
    let rotate = |v: D3| [dot3(x, v), dot3(y, v), dot3(z, v)];
    for atom in s.atoms.iter_mut() {
        let p = rotate([atom.x, atom.y, atom.z]);
        atom.x = p[0];
        atom.y = p[1];
        atom.z = p[2];
    }
    for row in s.cell_vectors.iter_mut() {
        *row = rotate(*row);
    }
    s.cell_offset = rotate(s.cell_offset);
    if s.cell_vectors[2][2] < 0.0 {
        for v in s.cell_vectors[2].iter_mut() {
            *v = -*v;
        }
    }
    s.cell_vectors[0][1] = 0.0;
    s.cell_vectors[0][2] = 0.0;
    s.cell_vectors[1][2] = 0.0;
}

/// The strain taking B's cell (columns v) onto A's (u) from the polar
/// decomposition, and the rotation that remains.
fn polar_strain(u1: V2, u2: V2, v1: V2, v2: V2, strain: &mut [f64; 3], rotation_deg: &mut f64) -> bool {
    let u: M2 = [u1, u2];
    let v: M2 = [v1, v2];
    if det(&v).abs() < 1e-12 {
        return false;
    }
    let f = mul(&u, &inverse(&v));
    let b = mul(&f, &transpose(&f));
    let d = det(&b);
    let root = (if 0.0 < d { d } else { 0.0 }).sqrt(); // std::max(0.0, d)
    let norm = (b[0][0] + b[1][1] + 2.0 * root).sqrt();
    if norm < 1e-12 {
        return false;
    }
    let s = mdiv(&madd(&b, &mscale(root, &I2)), norm);
    let r = mul(&inverse(&s), &f);
    *rotation_deg = r[0][1].atan2(r[0][0]) * 180.0 / PI;
    let phi = u1[1].atan2(u1[0]);
    let (c, sn) = (phi.cos(), phi.sin());
    let q: M2 = [[c, -sn], [sn, c]];
    let e = mul(&mul(&q, &msub(&s, &I2)), &transpose(&q));
    strain[0] = e[0][0];
    strain[1] = e[1][1];
    strain[2] = e[1][0];
    true
}

fn gcdl(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// The canonical basis {(p, 0), (q, s)}, 0 <= q < p, of the sublattice of
/// Z^2 spanned by two integer rows.
fn canonical(m: &[[i32; 2]; 2]) -> (i64, i64, i64) {
    let det = (m[0][0] as i64 * m[1][1] as i64 - m[0][1] as i64 * m[1][0] as i64).abs();
    let s = gcdl(m[0][1] as i64, m[1][1] as i64);
    if s == 0 || det == 0 {
        return (0, 0, 0);
    }
    let p = det / s;
    let (mut old_r, mut r, mut old_x, mut xx, mut old_y, mut yy) = (m[0][1] as i64, m[1][1] as i64, 1i64, 0i64, 0i64, 1i64);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_x, xx) = (xx, old_x - q * xx);
        (old_y, yy) = (yy, old_y - q * yy);
    }
    if old_r < 0 {
        old_x = -old_x;
        old_y = -old_y;
    }
    let mut q = (old_x * m[0][0] as i64 + old_y * m[1][0] as i64) % p;
    if q < 0 {
        q += p;
    }
    (p, q, s)
}

#[derive(Clone, Copy)]
struct Candidate {
    sa: [i32; 2],
    ij: [i32; 2],
    u: V2,
    v: V2,
}

/// The slab of one side with its layers, a along x, b in the xy plane, the
/// stacking direction along +z.
pub fn prepare_slab(side: &Side) -> Result<Structure, String> {
    if !side.crystal.has_unit_cell || side.crystal.atoms.is_empty() {
        return Err("needs a structure with a unit cell".into());
    }
    let mut slab = if side.cleave {
        if side.h == 0 && side.k == 0 && side.l == 0 {
            return Err("(0 0 0) is not a plane".into());
        }
        let params = SurfaceParams { h: side.h, k: side.k, l: side.l, layers: side.layers.max(1), vacuum: 0.0, ..Default::default() };
        let result = build_surface(&side.crystal, &params);
        if !result.success {
            return Err(result.message);
        }
        result.structure
    } else {
        repeat_layers_z(&side.crystal, side.layers.max(1))
    };
    let c = cross3(slab.cell_vectors[0], slab.cell_vectors[1]);
    if dot3(c, c).sqrt() < 1e-8 {
        return Err("its a and b cell vectors are parallel".into());
    }
    orient(&mut slab);
    Ok(slab)
}

fn llround(x: f64) -> i64 {
    x.round() as i64
}

/// Every match of the two surface cells within the limits, sorted by atoms
/// and then strain, with the Pareto-optimal ones marked.
pub fn find_matches(slab_a: &Structure, slab_b: &Structure, params: &SearchParams) -> Vec<Match> {
    let mut matches = Vec::new();
    let a1 = [slab_a.cell_vectors[0][0], slab_a.cell_vectors[0][1]];
    let a2 = [slab_a.cell_vectors[1][0], slab_a.cell_vectors[1][1]];
    let b1 = [slab_b.cell_vectors[0][0], slab_b.cell_vectors[0][1]];
    let b2 = [slab_b.cell_vectors[1][0], slab_b.cell_vectors[1][1]];
    let det_a = cross2(a1, a2);
    if det_a.abs() < 1e-10 || cross2(b1, b2).abs() < 1e-10 || params.angle_step <= 0.0 || params.max_repeats < 1 {
        return matches;
    }
    let (atoms_a, atoms_b) = (slab_a.atoms.len() as i32, slab_b.atoms.len() as i32);
    let min_angle = params.min_cell_angle.clamp(0.0, 89.0);
    let vector_tolerance = 3.0 * params.max_strain + 1e-3;
    let mut best: BTreeMap<(i64, i64, i64, i64, i64, i64), Match> = BTreeMap::new();
    let span = params.angle_max - params.angle_min;
    let steps = 1.max((span / params.angle_step + 1e-9).floor() as i32 + if span < 180.0 - 1e-9 { 1 } else { 0 });
    let mut kept: Vec<Candidate> = Vec::new();
    for step in 0..steps {
        let theta = params.angle_min + step as f64 * params.angle_step;
        let (c, s) = ((theta * PI / 180.0).cos(), (theta * PI / 180.0).sin());
        let c1 = [c * b1[0] - s * b1[1], s * b1[0] + c * b1[1]];
        let c2 = [c * b2[0] - s * b2[1], s * b2[0] + c * b2[1]];
        kept.clear();
        for i in -params.max_repeats..=params.max_repeats {
            for j in -params.max_repeats..=params.max_repeats {
                if i == 0 && j == 0 {
                    continue;
                }
                let v = add(scale(i as f64, c1), scale(j as f64, c2));
                let fa = cross2(v, a2) / det_a;
                let fb = cross2(a1, v) / det_a;
                let (sa0, sa1) = (llround(fa) as i32, llround(fb) as i32);
                if sa0 == 0 && sa1 == 0 {
                    continue;
                }
                let u = add(scale(sa0 as f64, a1), scale(sa1 as f64, a2));
                if length2(sub(v, u)) > vector_tolerance * length2(u) {
                    continue;
                }
                kept.push(Candidate { sa: [sa0, sa1], ij: [i, j], u, v });
            }
        }
        for p in 0..kept.len() {
            for q in p + 1..kept.len() {
                let (mut first, mut second) = (kept[p], kept[q]);
                let cross_u = cross2(first.u, second.u);
                let cross_v = cross2(first.v, second.v);
                if cross_u * cross_v <= 0.0 {
                    continue;
                }
                if cross_u < 0.0 {
                    std::mem::swap(&mut first, &mut second);
                }
                let area = cross_u.abs();
                if area > params.max_area || area < 1e-8 {
                    continue;
                }
                for _ in 0..64 {
                    if dot2(second.u, second.u) < dot2(first.u, first.u) {
                        std::mem::swap(&mut first, &mut second);
                    }
                    let mu = (dot2(first.u, second.u) / dot2(first.u, first.u)).round();
                    if mu == 0.0 {
                        break;
                    }
                    let k = mu as i32;
                    second.u = sub(second.u, scale(mu, first.u));
                    second.v = sub(second.v, scale(mu, first.v));
                    second.sa[0] -= k * first.sa[0];
                    second.sa[1] -= k * first.sa[1];
                    second.ij[0] -= k * first.ij[0];
                    second.ij[1] -= k * first.ij[1];
                }
                if cross2(first.u, second.u) < 0.0 {
                    second.u = [-second.u[0], -second.u[1]];
                    second.v = [-second.v[0], -second.v[1]];
                    second.sa = [-second.sa[0], -second.sa[1]];
                    second.ij = [-second.ij[0], -second.ij[1]];
                }
                let cosine = dot2(first.u, second.u) / (length2(first.u) * length2(second.u));
                let angle = cosine.clamp(-1.0, 1.0).acos() * 180.0 / PI;
                if angle < min_angle || angle > 180.0 - min_angle {
                    continue;
                }
                let (mut e, mut residual) = ([0.0; 3], 0.0);
                if !polar_strain(first.u, second.u, first.v, second.v, &mut e, &mut residual) {
                    continue;
                }
                let mean = (e[0].abs() + e[1].abs() + e[2].abs()) / 3.0;
                if mean > params.max_strain {
                    continue;
                }
                let mut m = Match { a: [first.sa, second.sa], b: [first.ij, second.ij], ..Default::default() };
                m.cells_a = m.a[0][0] * m.a[1][1] - m.a[0][1] * m.a[1][0];
                m.cells_b = m.b[0][0] * m.b[1][1] - m.b[0][1] * m.b[1][0];
                if m.cells_a <= 0 || m.cells_b <= 0 {
                    continue;
                }
                m.theta_deg = (theta + residual + 540.0) % 360.0 - 180.0;
                m.strain = e;
                m.mean_strain = mean;
                m.area = area;
                m.cell_angle = angle;
                m.atoms = atoms_a * m.cells_a + atoms_b * m.cells_b;
                let (pa, qa, sa) = canonical(&m.a);
                let (pb, qb, sb) = canonical(&m.b);
                let key = (pa, qa, sa, pb, qb, sb);
                let replace = match best.get(&key) {
                    None => true,
                    Some(f) => mean < f.mean_strain - 1e-12 || (mean < f.mean_strain + 1e-12 && m.theta_deg.abs() < f.theta_deg.abs() - 1e-9),
                };
                if replace {
                    best.insert(key, m);
                }
            }
        }
    }
    let mut candidates: Vec<Match> = best.into_values().collect();
    candidates.sort_by(|x, y| x.theta_deg.abs().partial_cmp(&y.theta_deg.abs()).unwrap_or(std::cmp::Ordering::Equal));
    let mut seen = std::collections::BTreeSet::new();
    for m in candidates {
        let l1 = (m.a[0][0] as f64 * a1[0] + m.a[0][1] as f64 * a2[0]).hypot(m.a[0][0] as f64 * a1[1] + m.a[0][1] as f64 * a2[1]);
        let l2 = (m.a[1][0] as f64 * a1[0] + m.a[1][1] as f64 * a2[0]).hypot(m.a[1][0] as f64 * a1[1] + m.a[1][1] as f64 * a2[1]);
        let lo = if l2 < l1 { l2 } else { l1 };
        let hi = if l1 < l2 { l2 } else { l1 };
        let other = 180.0 - m.cell_angle;
        let ang = if other < m.cell_angle { other } else { m.cell_angle };
        let shape = (m.atoms, llround(lo * 1e3), llround(hi * 1e3), llround(ang * 1e2), llround(m.mean_strain * 1e6));
        if seen.insert(shape) {
            matches.push(m);
        }
    }
    sort_by(&mut matches, |x, y| {
        if x.atoms != y.atoms {
            return x.atoms < y.atoms;
        }
        if x.mean_strain != y.mean_strain {
            return x.mean_strain < y.mean_strain;
        }
        x.theta_deg < y.theta_deg
    });
    let within = |big: &[[i32; 2]; 2], small: &[[i32; 2]; 2]| {
        let det = small[0][0] as i64 * small[1][1] as i64 - small[0][1] as i64 * small[1][0] as i64;
        if det == 0 {
            return false;
        }
        for r in 0..2 {
            let x0 = big[r][0] as i64 * small[1][1] as i64 - big[r][1] as i64 * small[1][0] as i64;
            let x1 = big[r][1] as i64 * small[0][0] as i64 - big[r][0] as i64 * small[0][1] as i64;
            if x0 % det != 0 || x1 % det != 0 {
                return false;
            }
        }
        true
    };
    let mut distinct: Vec<Match> = Vec::new();
    for m in &matches {
        let mut redundant = false;
        for smaller in &distinct {
            if smaller.cells_a >= m.cells_a
                || m.cells_a % smaller.cells_a != 0
                || m.cells_b % smaller.cells_b != 0
                || m.cells_a / smaller.cells_a != m.cells_b / smaller.cells_b
            {
                continue;
            }
            if m.mean_strain < smaller.mean_strain - 1e-9 {
                continue;
            }
            let same_rotation = ((m.theta_deg - smaller.theta_deg + 540.0) % 360.0 - 180.0).abs() <= 0.5;
            if (same_rotation && within(&m.a, &smaller.a) && within(&m.b, &smaller.b)) || (m.mean_strain - smaller.mean_strain).abs() < 1e-6 {
                redundant = true;
                break;
            }
        }
        if !redundant {
            distinct.push(*m);
        }
    }
    let mut matches = distinct;
    for m in matches.iter_mut() {
        if m.theta_deg.abs() < 0.005 {
            m.theta_deg = 0.0;
        }
    }
    matches.truncate(params.max_matches);
    let mut lowest = 1e300;
    for m in matches.iter_mut() {
        if m.mean_strain < lowest - 1e-12 {
            m.pareto = true;
            lowest = m.mean_strain;
        }
    }
    matches
}

/// Of the Pareto-optimal matches, the one with the lowest mean strain (%)
/// plus area / 25 A^2; -1 when none.
pub fn suggested_match(matches: &[Match]) -> i32 {
    let (mut chosen, mut score) = (-1, 1e300);
    for (i, m) in matches.iter().enumerate() {
        if !m.pareto {
            continue;
        }
        let value = m.mean_strain * 100.0 + m.area / 25.0;
        if value < score {
            score = value;
            chosen = i as i32;
        }
    }
    chosen
}

fn floor2(v: V2) -> V2 {
    [v[0].floor(), v[1].floor()]
}

/// The slab repeated in its plane by the integer rows of m.
pub fn in_plane_supercell(slab: &Structure, m: &[[i32; 2]; 2]) -> Structure {
    let a = [slab.cell_vectors[0][0], slab.cell_vectors[0][1]];
    let b = [slab.cell_vectors[1][0], slab.cell_vectors[1][1]];
    let big_a = add(scale(m[0][0] as f64, a), scale(m[0][1] as f64, b));
    let big_b = add(scale(m[1][0] as f64, a), scale(m[1][1] as f64, b));
    let mut result = slab.clone();
    result.atoms.clear();
    result.cell_vectors[0] = [big_a[0], big_a[1], 0.0];
    result.cell_vectors[1] = [big_b[0], big_b[1], 0.0];
    let base: M2 = [a, b];
    let cell: M2 = [big_a, big_b];
    if det(&cell).abs() < 1e-12 {
        return result;
    }
    let (inverse_base, inverse_cell) = (inverse(&base), inverse(&cell));
    let min4 = |w: i32, x: i32, y: i32, z: i32| w.min(x).min(y).min(z);
    let max4 = |w: i32, x: i32, y: i32, z: i32| w.max(x).max(y).max(z);
    let lo0 = min4(0, m[0][0], m[1][0], m[0][0] + m[1][0]) - 1;
    let hi0 = max4(0, m[0][0], m[1][0], m[0][0] + m[1][0]) + 1;
    let lo1 = min4(0, m[0][1], m[1][1], m[0][1] + m[1][1]) - 1;
    let hi1 = max4(0, m[0][1], m[1][1], m[0][1] + m[1][1]) + 1;
    const EPS: f64 = 1e-8;
    for atom in &slab.atoms {
        let f = mulv(&inverse_base, [atom.x, atom.y]);
        let f = sub(f, floor2([f[0] + EPS, f[1] + EPS]));
        let home = mulv(&base, f);
        for n0 in lo0..=hi0 {
            for n1 in lo1..=hi1 {
                let p = add(add(home, scale(n0 as f64, a)), scale(n1 as f64, b));
                let g = mulv(&inverse_cell, p);
                if g[0] < -EPS || g[0] >= 1.0 - EPS || g[1] < -EPS || g[1] >= 1.0 - EPS {
                    continue;
                }
                let mut copy = atom.clone();
                copy.x = p[0];
                copy.y = p[1];
                result.atoms.push(copy);
            }
        }
    }
    result
}

/// The interface of a match: slab A below, slab B above, B rotated by the
/// match's angle and strained as asked.
pub fn build_interface(slab_a: &Structure, slab_b: &Structure, m: &Match, assembly: &Assembly) -> Structure {
    let mut a = in_plane_supercell(slab_a, &m.a);
    let mut b = in_plane_supercell(slab_b, &m.b);
    let u: M2 = [[a.cell_vectors[0][0], a.cell_vectors[0][1]], [a.cell_vectors[1][0], a.cell_vectors[1][1]]];
    let v: M2 = [[b.cell_vectors[0][0], b.cell_vectors[0][1]], [b.cell_vectors[1][0], b.cell_vectors[1][1]]];
    let (c, s) = ((m.theta_deg * PI / 180.0).cos(), (m.theta_deg * PI / 180.0).sin());
    let r: M2 = [[c, s], [-s, c]];
    let rotated_v = mul(&r, &v);
    let w = match assembly.strain_on {
        StrainOn::B => u,
        StrainOn::A => rotated_v,
        StrainOn::Both => mscale(0.5, &madd(&u, &rotated_v)),
    };
    let map_a = mul(&w, &inverse(&u));
    let map_b = mul(&w, &inverse(&v));
    let transform = |st: &Structure, m: &M2| apply_transform_2d(st, &[[m[0][0], m[1][0]], [m[0][1], m[1][1]]]);
    a = transform(&a, &map_a);
    b = transform(&b, &map_b);
    let z_range = |st: &Structure| {
        let (mut low, mut high) = (1e300f64, -1e300f64);
        for atom in &st.atoms {
            low = if atom.z < low { atom.z } else { low };
            high = if high < atom.z { atom.z } else { high };
        }
        (low, high)
    };
    let (low_a, high_a) = z_range(&a);
    let (low_b, high_b) = z_range(&b);
    let (height_a, height_b) = (high_a - low_a, high_b - low_b);
    let gap = if 0.0 < assembly.gap { assembly.gap } else { 0.0 };
    let vacuum = if 0.0 < assembly.vacuum { assembly.vacuum } else { 0.0 };
    let base = if vacuum > 0.0 { 0.5 * vacuum } else { 0.5 * gap };
    let height = height_a + gap + height_b + if vacuum > 0.0 { vacuum } else { gap };
    let shift = add(scale(assembly.shift[0], w[0]), scale(assembly.shift[1], w[1]));
    let mut result = Structure {
        has_unit_cell: true,
        cell_vectors: [[w[0][0], w[0][1], 0.0], [w[1][0], w[1][1], 0.0], [0.0, 0.0, height]],
        cell_offset: [0.0; 3],
        ..Default::default()
    };
    for atom in &a.atoms {
        let mut atom: AtomSite = atom.clone();
        atom.z = atom.z - low_a + base;
        result.atoms.push(atom);
    }
    for atom in &b.atoms {
        let mut atom: AtomSite = atom.clone();
        atom.x += shift[0];
        atom.y += shift[1];
        atom.z = atom.z - low_b + base + height_a + gap;
        result.atoms.push(atom);
    }
    let inverse_w = inverse(&w);
    for atom in result.atoms.iter_mut() {
        let f = mulv(&inverse_w, [atom.x, atom.y]);
        let f = sub(f, floor2([f[0] + 1e-9, f[1] + 1e-9]));
        let p = mulv(&w, f);
        atom.x = p[0];
        atom.y = p[1];
    }
    result
}
