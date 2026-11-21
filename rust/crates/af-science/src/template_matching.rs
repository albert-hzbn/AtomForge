//! Polyhedral template matching and what builds on it (the C++
//! science/TemplateMatching): structure types, orientations and lattice
//! constants, chemical ordering, grain segmentation, elastic strain and fcc
//! planar faults. Kabsch fits use af_core::eigen::jacobi_svd3, Eigen's own
//! algorithm, and 3x3 products use Eigen's summation order (eigen_matmul3,
//! eigen_mul3), so every fit and tie-break matches the C++ bit for bit.

use std::collections::BTreeMap;
use std::sync::OnceLock;

use af_core::cxx::{partial_sort_by, sort_by, stable_sort_by};
use af_core::eigen::{determinant3, eigen_matmul3 as matmul3, inverse3, jacobi_svd3, eigen_mul3 as mul3, transpose3, Matrix3 as M3, IDENTITY3};
use af_core::json::Json;

use crate::analysis::{ji, jv, R};
use crate::clusters::union_root;
use crate::core::{finite_array, integer, neighbor_list, positive};
use crate::data::Parameters;
use crate::local_environment::{add_classification, read_atom_set, AtomSet};
use crate::voronoi::mean_atomic_spacing;

type V3 = [f64; 3];
const PI: f64 = std::f64::consts::PI;

fn ideal_ca() -> f64 {
    (8.0f64 / 3.0).sqrt()
}

const OTHER: i32 = 0;
const FCC: i32 = 1;
const HCP: i32 = 2;
const BCC: i32 = 3;
const ICO: i32 = 4;
const SC: i32 = 5;
const CUBIC_DIAMOND: i32 = 6;
const HEX_DIAMOND: i32 = 7;
const KIND_COUNT: usize = 8;
const NAMES: [&str; KIND_COUNT] = ["other", "fcc", "hcp", "bcc", "icosahedral", "simple cubic", "cubic diamond", "hexagonal diamond"];
const KEYS: [&str; KIND_COUNT] = ["other", "fcc", "hcp", "bcc", "ico", "sc", "cubic-diamond", "hex-diamond"];

// ---------------------------------------------------------------- vector helpers
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn squared(a: V3) -> f64 {
    dot(a, a)
}
fn norm(a: V3) -> f64 {
    squared(a).sqrt()
}
/// Eigen's normalized(): each component divided by the norm.
fn normalized(a: V3) -> V3 {
    let n = squared(a);
    if n > 0.0 {
        let n = n.sqrt();
        [a[0] / n, a[1] / n, a[2] / n]
    } else {
        a
    }
}
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn outer(a: V3, b: V3) -> M3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i] * b[j]))
}
fn madd(a: &M3, b: &M3) -> M3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}

// ---------------------------------------------------------------- quaternions (Eigen)
#[derive(Clone, Copy, Debug)]
struct Quat {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

impl Quat {
    const IDENTITY: Quat = Quat { w: 1.0, x: 0.0, y: 0.0, z: 0.0 };

    /// Eigen's quaternion of a rotation matrix.
    fn from_matrix(m: &M3) -> Quat {
        let t = m[0][0] + m[1][1] + m[2][2];
        if t > 0.0 {
            let t = (t + 1.0).sqrt();
            let w = 0.5 * t;
            let t = 0.5 / t;
            Quat { w, x: (m[2][1] - m[1][2]) * t, y: (m[0][2] - m[2][0]) * t, z: (m[1][0] - m[0][1]) * t }
        } else {
            let mut i = 0;
            if m[1][1] > m[0][0] {
                i = 1;
            }
            if m[2][2] > m[i][i] {
                i = 2;
            }
            let j = (i + 1) % 3;
            let k = (j + 1) % 3;
            let t = (m[i][i] - m[j][j] - m[k][k] + 1.0).sqrt();
            let mut c = [0.0; 3];
            c[i] = 0.5 * t;
            let t = 0.5 / t;
            let w = (m[k][j] - m[j][k]) * t;
            c[j] = (m[j][i] + m[i][j]) * t;
            c[k] = (m[k][i] + m[i][k]) * t;
            Quat { w, x: c[0], y: c[1], z: c[2] }
        }
    }

    fn from_angle_axis(angle: f64, axis: V3) -> Quat {
        let ha = 0.5 * angle;
        let s = ha.sin();
        Quat { w: ha.cos(), x: s * axis[0], y: s * axis[1], z: s * axis[2] }
    }

    fn mul(self, b: Quat) -> Quat {
        let a = self;
        Quat {
            w: a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            x: a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            y: a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
            z: a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
        }
    }

    fn conjugate(self) -> Quat {
        Quat { w: self.w, x: -self.x, y: -self.y, z: -self.z }
    }

    /// (x, y, z, w), Eigen's coefficient order.
    fn coeffs(self) -> [f64; 4] {
        [self.x, self.y, self.z, self.w]
    }

    fn from_coeffs(c: [f64; 4]) -> Quat {
        Quat { x: c[0], y: c[1], z: c[2], w: c[3] }
    }

    fn normalized(self) -> Quat {
        let c = self.coeffs();
        let n = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2] + c[3] * c[3]).sqrt();
        Quat::from_coeffs(c.map(|v| v / n))
    }

    fn to_matrix(self) -> M3 {
        let (tx, ty, tz) = (2.0 * self.x, 2.0 * self.y, 2.0 * self.z);
        let (twx, twy, twz) = (tx * self.w, ty * self.w, tz * self.w);
        let (txx, txy, txz) = (tx * self.x, ty * self.x, tz * self.x);
        let (tyy, tyz, tzz) = (ty * self.y, tz * self.y, tz * self.z);
        [[1.0 - (tyy + tzz), txy - twz, txz + twy], [txy + twz, 1.0 - (txx + tzz), tyz - twx], [txz - twy, tyz + twx, 1.0 - (txx + tyy)]]
    }
}

fn quaternion_json(q: &Quat) -> Json {
    jv(&[q.w, q.x, q.y, q.z])
}

// ---------------------------------------------------------------- templates
type Cell = [V3; 3];

fn cartesian(cell: &Cell, f: V3) -> V3 {
    add(add(scale(cell[0], f[0]), scale(cell[1], f[1])), scale(cell[2], f[2]))
}

fn shell(ideal: &Cell, actual: &Cell, basis: &[V3], count: usize) -> Vec<V3> {
    let mut candidates: Vec<(V3, V3, i64)> = Vec::new();
    for i in -3..=3 {
        for j in -3..=3 {
            for k in -3..=3 {
                for b in basis {
                    let f = sub(add(*b, [i as f64, j as f64, k as f64]), basis[0]);
                    let v = cartesian(ideal, f);
                    if norm(v) < 1e-9 {
                        continue;
                    }
                    candidates.push((v, cartesian(actual, f), (norm(v) * 1e6).round() as i64));
                }
            }
        }
    }
    stable_sort_by(&mut candidates, |a, b| a.2 < b.2);
    assert!(candidates.len() > count && candidates[count].2 != candidates[count - 1].2, "template shell is not closed");
    candidates[..count].iter().map(|c| c.1).collect()
}

fn hexagonal_cell(ca: f64) -> Cell {
    [[1.0, 0.0, 0.0], [-0.5, 3.0f64.sqrt() / 2.0, 0.0], [0.0, 0.0, ca]]
}

const CUBIC: Cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn template_vectors(kind: i32, ca: f64) -> Vec<V3> {
    let fcc_basis = vec![[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]];
    match kind {
        FCC => shell(&CUBIC, &CUBIC, &fcc_basis, 12),
        HCP => shell(&hexagonal_cell(ideal_ca()), &hexagonal_cell(ca), &[[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.5]], 12),
        BCC => shell(&CUBIC, &CUBIC, &[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], 14),
        SC => shell(&CUBIC, &CUBIC, &[[0.0, 0.0, 0.0]], 6),
        CUBIC_DIAMOND => {
            let mut basis = fcc_basis.clone();
            for b in &fcc_basis {
                basis.push(add(*b, [0.25, 0.25, 0.25]));
            }
            shell(&CUBIC, &CUBIC, &basis, 16)
        }
        HEX_DIAMOND => {
            let z = 1.0 / 16.0;
            let basis = [[1.0 / 3.0, 2.0 / 3.0, z], [2.0 / 3.0, 1.0 / 3.0, 0.5 + z], [2.0 / 3.0, 1.0 / 3.0, 1.0 - z], [1.0 / 3.0, 2.0 / 3.0, 0.5 - z]];
            shell(&hexagonal_cell(ideal_ca()), &hexagonal_cell(ca), &basis, 16)
        }
        ICO => {
            let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
            let mut out = Vec::new();
            for s in [-1.0, 1.0] {
                for t in [-phi, phi] {
                    out.push([0.0, s, t]);
                    out.push([s, t, 0.0]);
                    out.push([t, 0.0, s]);
                }
            }
            let d = (1.0 + phi * phi).sqrt();
            out.into_iter().map(|v| [v[0] / d, v[1] / d, v[2] / d]).collect()
        }
        _ => panic!("unknown template"),
    }
}

/// Orthonormal frame (columns) with e1 along u and e2 in the (u, w) plane.
fn frame(u: V3, w: V3) -> M3 {
    let e1 = normalized(u);
    let e2 = normalized(sub(w, scale(e1, dot(w, e1))));
    let e3 = cross(e1, e2);
    [[e1[0], e2[0], e3[0]], [e1[1], e2[1], e3[1]], [e1[2], e2[2], e3[2]]]
}

fn parallel(a: V3, b: V3) -> bool {
    norm(cross(a, b)) < 0.05 * norm(a) * norm(b)
}

struct Template {
    unit: Vec<V3>,
    mean_length: f64,
    symmetry_quats: Vec<Quat>,
    anchors: Vec<(usize, usize)>,
    anchor_ratio: Vec<f64>,
    anchor_cos: Vec<f64>,
}

fn make_template(kind: i32) -> Template {
    let vectors = template_vectors(kind, ideal_ca());
    let mut mean_length = 0.0;
    for v in &vectors {
        mean_length += norm(*v);
    }
    mean_length /= vectors.len() as f64;
    let u: Vec<V3> = vectors.iter().map(|v| [v[0] / mean_length, v[1] / mean_length, v[2] / mean_length]).collect();
    let n = u.len();
    let mut b = 1;
    while parallel(u[0], u[b]) {
        b += 1;
    }
    let reference = frame(u[0], u[b]);
    let mut symmetry_quats = Vec::new();
    let mut permutations: Vec<Vec<i32>> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let (a, c) = (u[i], u[j]);
            if i == j || (norm(a) - norm(u[0])).abs() > 1e-6 || (norm(c) - norm(u[b])).abs() > 1e-6 || (dot(a, c) - dot(u[0], u[b])).abs() > 1e-6 {
                continue;
            }
            let rotation = matmul3(&frame(a, c), &transpose3(&reference));
            let mut permutation = vec![-1i32; n];
            let mut ok = true;
            let mut k = 0;
            while k < n && ok {
                let image = mul3(&rotation, u[k]);
                for m in 0..n {
                    if norm(sub(image, u[m])) < 1e-6 {
                        permutation[k] = m as i32;
                    }
                }
                ok = permutation[k] >= 0;
                k += 1;
            }
            if !ok {
                continue;
            }
            symmetry_quats.push(Quat::from_matrix(&rotation).normalized());
            permutations.push(permutation);
        }
    }
    let (mut anchors, mut anchor_ratio, mut anchor_cos) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        for j in 0..n {
            if i == j || parallel(u[i], u[j]) {
                continue;
            }
            let smallest = permutations.iter().all(|perm| !((perm[i], perm[j]) < (i as i32, j as i32)));
            if !smallest {
                continue;
            }
            let (a, c) = (u[i], u[j]);
            anchors.push((i, j));
            anchor_ratio.push(norm(c) / norm(a));
            anchor_cos.push(dot(a, c) / (norm(a) * norm(c)));
        }
    }
    Template { unit: u, mean_length, symmetry_quats, anchors, anchor_ratio, anchor_cos }
}

fn template_of(kind: i32) -> &'static Template {
    static TEMPLATES: OnceLock<Vec<Option<Template>>> = OnceLock::new();
    TEMPLATES.get_or_init(|| (0..KIND_COUNT as i32).map(|k| if k == 0 { None } else { Some(make_template(k)) }).collect())[kind as usize].as_ref().unwrap()
}

fn shell_size(kind: i32) -> usize {
    template_of(kind).unit.len()
}

// ---------------------------------------------------------------- fitting
fn hungarian(cost: &[Vec<f64>]) -> Vec<i32> {
    let n = cost.len();
    let inf = f64::INFINITY;
    let (mut u, mut v, mut minv) = (vec![0.0; n + 1], vec![0.0; n + 1], vec![0.0; n + 1]);
    let (mut p, mut way) = (vec![0usize; n + 1], vec![0usize; n + 1]);
    let mut used = vec![false; n + 1];
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        minv.iter_mut().for_each(|x| *x = inf);
        used.iter_mut().for_each(|x| *x = false);
        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = inf;
            let mut j1 = 0;
            for j in 1..=n {
                if used[j] {
                    continue;
                }
                let current = cost[i0 - 1][j - 1] - u[i0] - v[j];
                if current < minv[j] {
                    minv[j] = current;
                    way[j] = j0;
                }
                if minv[j] < delta {
                    delta = minv[j];
                    j1 = j;
                }
            }
            for j in 0..=n {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }
    let mut row_to_column = vec![0i32; n];
    for j in 1..=n {
        row_to_column[p[j] - 1] = (j - 1) as i32;
    }
    row_to_column
}

fn assign(t: &[V3], b: &[V3], rotation: &M3) -> Vec<i32> {
    let n = b.len();
    let rotated: Vec<V3> = (0..n).map(|j| mul3(rotation, t[j])).collect();
    let mut map = vec![0i32; n];
    let mut used = vec![0i32; n];
    let mut bijective = true;
    for k in 0..n {
        let mut best = f64::INFINITY;
        for j in 0..n {
            let d = squared(sub(rotated[j], b[k]));
            if d < best {
                best = d;
                map[k] = j as i32;
            }
        }
        if used[map[k] as usize] > 0 {
            bijective = false;
        }
        used[map[k] as usize] += 1;
    }
    if bijective {
        return map;
    }
    let cost: Vec<Vec<f64>> = (0..n).map(|k| (0..n).map(|j| squared(sub(rotated[j], b[k]))).collect()).collect();
    hungarian(&cost)
}

/// The rotation minimising sum |R t_map(k) - b_k|^2 (Kabsch, Eigen's JacobiSVD).
fn kabsch(t: &[V3], b: &[V3], map: &[i32]) -> M3 {
    let mut a = [[0.0; 3]; 3];
    for k in 0..b.len() {
        a = madd(&a, &outer(b[k], t[map[k] as usize]));
    }
    let (u, _, v) = jacobi_svd3(&a);
    let vt = transpose3(&v);
    let mut d = IDENTITY3;
    if determinant3(&matmul3(&u, &vt)) < 0.0 {
        d[2][2] = -1.0;
    }
    matmul3(&matmul3(&u, &d), &vt)
}

#[derive(Clone)]
struct Fit {
    kind: i32,
    rmsd: f64,
    scale: f64,
    rotation: M3,
    map: Vec<i32>,
}

impl Default for Fit {
    fn default() -> Self {
        Self { kind: OTHER, rmsd: f64::INFINITY, scale: 0.0, rotation: IDENTITY3, map: Vec::new() }
    }
}

fn fit_template(kind: i32, data: &[V3]) -> Fit {
    let tpl = template_of(kind);
    let n = data.len();
    let mut mean = 0.0;
    for v in data {
        mean += norm(*v);
    }
    mean /= n as f64;
    let b: Vec<V3> = data.iter().map(|v| [v[0] / mean, v[1] / mean, v[2] / mean]).collect();
    let mut first = 0;
    for k in 1..n {
        if norm(b[k]) < norm(b[first]) {
            first = k;
        }
    }
    let mut by_angle: Vec<(f64, usize)> =
        (0..n).filter(|&k| k != first && !parallel(b[first], b[k])).map(|k| (-dot(b[first], b[k]) / (norm(b[first]) * norm(b[k])), k)).collect();
    sort_by(&mut by_angle, |x, y| x < y);
    let mut best = Fit { kind, ..Default::default() };
    for c in 0..by_angle.len().min(2) {
        let second = by_angle[c].1;
        let ratio = norm(b[second]) / norm(b[first]);
        let cosine = -by_angle[c].0;
        let data_frame = frame(b[first], b[second]);
        for a in 0..tpl.anchors.len() {
            if (tpl.anchor_ratio[a] - ratio).abs() > 0.12 * tpl.anchor_ratio[a] {
                continue;
            }
            if (tpl.anchor_cos[a].clamp(-1.0, 1.0).acos() - cosine.clamp(-1.0, 1.0).acos()).abs() > 15.0 * PI / 180.0 {
                continue;
            }
            let (i, j) = tpl.anchors[a];
            let mut rotation = matmul3(&data_frame, &transpose3(&frame(tpl.unit[i], tpl.unit[j])));
            let mut map = assign(&tpl.unit, &b, &rotation);
            rotation = kabsch(&tpl.unit, &b, &map);
            map = assign(&tpl.unit, &b, &rotation);
            rotation = kabsch(&tpl.unit, &b, &map);
            let (mut crossing, mut template_norm) = (0.0, 0.0);
            for k in 0..n {
                let rt = mul3(&rotation, tpl.unit[map[k] as usize]);
                crossing += dot(b[k], rt);
                template_norm += squared(rt);
            }
            let s = crossing / template_norm;
            let mut residual = 0.0;
            for k in 0..n {
                // Eigen folds s * (R * u) into (s * R) * u.
                let scaled: M3 = std::array::from_fn(|r| std::array::from_fn(|c| s * rotation[r][c]));
                residual += squared(sub(b[k], mul3(&scaled, tpl.unit[map[k] as usize])));
            }
            let rmsd = (residual / n as f64).sqrt();
            if rmsd < best.rmsd {
                best.rmsd = rmsd;
                best.rotation = rotation;
                best.map = map;
                best.scale = mean * s / tpl.mean_length;
            }
        }
        if best.rmsd < 0.02 {
            break;
        }
    }
    best
}

fn lattice_symmetry(kind: i32) -> &'static [Quat] {
    static HEXAGONAL: OnceLock<Vec<Quat>> = OnceLock::new();
    if kind == HCP || kind == HEX_DIAMOND {
        return HEXAGONAL.get_or_init(|| {
            let mut out = Vec::new();
            for k in 0..6 {
                out.push(Quat::from_angle_axis(k as f64 * PI / 3.0, [0.0, 0.0, 1.0]));
                let a = k as f64 * PI / 6.0;
                out.push(Quat::from_angle_axis(PI, [a.cos(), a.sin(), 0.0]));
            }
            out
        });
    }
    if kind == ICO {
        return &template_of(ICO).symmetry_quats;
    }
    &template_of(FCC).symmetry_quats
}

fn reduced_orientation(kind: i32, rotation: &M3) -> Quat {
    let q = Quat::from_matrix(rotation);
    let mut best = q;
    let mut best_w = -1.0;
    for s in lattice_symmetry(kind) {
        let candidate = q.mul(*s).normalized();
        if candidate.w.abs() > best_w {
            best_w = candidate.w.abs();
            best = candidate;
        }
    }
    if best.w < 0.0 {
        best = Quat::from_coeffs(best.coeffs().map(|c| c * -1.0));
    }
    best
}

fn misorientation(kind: i32, a: &Quat, b: &Quat) -> f64 {
    let delta = a.conjugate().mul(*b);
    let mut best_w = 0.0f64;
    for s in lattice_symmetry(kind) {
        let w = delta.mul(*s).w.abs();
        if best_w < w {
            best_w = w;
        }
    }
    2.0 * (if 1.0 < best_w { 1.0 } else { best_w }).acos()
}

// ---------------------------------------------------------------- neighbours and matching
struct Neighbors {
    index: Vec<Vec<i32>>,
    vectors: Vec<Vec<V3>>,
}

fn kind_from_key(key: &str) -> Result<i32, String> {
    for k in 1..KIND_COUNT {
        if key == KEYS[k] {
            return Ok(k as i32);
        }
    }
    Err(format!("Unknown structure '{key}'; choose from {}", KEYS[1..].join(", ")))
}

fn enabled_kinds(p: &Parameters) -> Result<Vec<i32>, String> {
    if !p.has("structures") {
        return Ok(vec![FCC, HCP, BCC, ICO]);
    }
    let list = p.json("structures")?;
    let kinds = match list {
        Json::String(s) => vec![kind_from_key(s)?],
        _ => list.items()?.iter().map(|item| kind_from_key(item.string()?)).collect::<Result<Vec<_>, _>>()?,
    };
    if kinds.is_empty() {
        return Err("structures must name at least one structure".into());
    }
    Ok(kinds)
}

struct Geometry {
    atoms: AtomSet,
    cutoff: f64,
    automatic: bool,
}

fn geometry(p: &Parameters) -> Result<Geometry, String> {
    let atoms = read_atom_set(p, "positions")?;
    let mut cutoff = p.number_or("cutoff_A", 0.0)?;
    if cutoff < 0.0 {
        return Err("cutoff_A must be nonnegative".into());
    }
    let mut automatic = false;
    if cutoff == 0.0 && !atoms.positions.is_empty() {
        cutoff = 1.8 * mean_atomic_spacing(&atoms)?;
        automatic = true;
    }
    Ok(Geometry { atoms, cutoff, automatic })
}

fn nearest_neighbors(g: &mut Geometry, needed: usize) -> Result<Neighbors, String> {
    let n = g.atoms.positions.len();
    let mut attempt = 0;
    loop {
        let mut result = Neighbors { index: vec![Vec::new(); n], vectors: vec![Vec::new(); n] };
        let mut lists: Vec<Vec<(f64, i32, V3)>> = vec![Vec::new(); n];
        for link in neighbor_list(&g.atoms.positions, &g.atoms.cell, g.atoms.pbc, g.cutoff)? {
            let v = link.vector;
            if norm(v) < 1e-10 {
                return Err("Coincident atoms make local environment analysis undefined".into());
            }
            lists[link.i].push((norm(v), link.j as i32, v));
        }
        let mut enough = 0;
        for i in 0..n {
            let list = &mut lists[i];
            let keep = list.len().min(16);
            partial_sort_by(list, keep, |a, b| a.0 < b.0);
            for item in &list[..keep] {
                result.index[i].push(item.1);
                result.vectors[i].push(item.2);
            }
            if list.len() >= needed {
                enough += 1;
            }
        }
        if !g.automatic || attempt >= 3 || n == 0 || 2 * enough >= n {
            return Ok(result);
        }
        g.cutoff *= 1.3;
        attempt += 1;
    }
}

fn environment(kind: i32, i: usize, nb: &Neighbors, out: &mut Vec<V3>, atoms: &mut Vec<i32>) -> bool {
    out.clear();
    atoms.clear();
    let vectors = &nb.vectors[i];
    if kind != CUBIC_DIAMOND && kind != HEX_DIAMOND {
        let count = shell_size(kind);
        if vectors.len() < count {
            return false;
        }
        out.extend_from_slice(&vectors[..count]);
        atoms.extend_from_slice(&nb.index[i][..count]);
        return true;
    }
    if vectors.len() < 4 {
        return false;
    }
    for k in 0..4 {
        out.push(vectors[k]);
        atoms.push(nb.index[i][k]);
    }
    for k in 0..4 {
        let j = nb.index[i][k] as usize;
        let second = &nb.vectors[j];
        if second.len() < 4 {
            return false;
        }
        let mut back = 0;
        for m in 1..4 {
            if norm(add(second[m], vectors[k])) < norm(add(second[back], vectors[k])) {
                back = m;
            }
        }
        for m in 0..4 {
            if m != back {
                out.push(add(vectors[k], second[m]));
                atoms.push(nb.index[j][m]);
            }
        }
    }
    true
}

struct Matching {
    cutoff: f64,
    neighbors: Neighbors,
    fits: Vec<Fit>,
}

fn match_all(p: &Parameters, kinds: &[i32], rmsd_cutoff: f64) -> Result<Matching, String> {
    let mut needed = 4;
    for &k in kinds {
        if k != CUBIC_DIAMOND && k != HEX_DIAMOND {
            needed = needed.max(shell_size(k));
        }
    }
    let mut g = geometry(p)?;
    let neighbors = nearest_neighbors(&mut g, needed)?;
    let n = g.atoms.positions.len();
    let mut fits = vec![Fit::default(); n];
    let (mut env, mut atoms) = (Vec::new(), Vec::new());
    for (i, slot) in fits.iter_mut().enumerate() {
        let mut best = Fit::default();
        for &kind in kinds {
            if !environment(kind, i, &neighbors, &mut env, &mut atoms) {
                continue;
            }
            let fit = fit_template(kind, &env);
            if fit.rmsd < best.rmsd {
                best = fit;
            }
        }
        if best.kind == OTHER || !(best.rmsd <= rmsd_cutoff || rmsd_cutoff <= 0.0) {
            best = Fit { rmsd: best.rmsd, ..Default::default() };
        }
        *slot = best;
    }
    Ok(Matching { cutoff: g.cutoff, neighbors, fits })
}

fn rmsd_cutoff_of(p: &Parameters) -> Result<f64, String> {
    let cutoff = p.number_or("rmsd_cutoff", 0.1)?;
    if cutoff < 0.0 {
        return Err("rmsd_cutoff must be nonnegative (0 = no cutoff)".into());
    }
    Ok(cutoff)
}

fn coplanar(v: &[V3]) -> bool {
    if v.len() < 3 {
        return true;
    }
    let mut normal = [0.0; 3];
    'outer: for a in 0..v.len() {
        for b in a + 1..v.len() {
            if norm(normal) >= 1e-6 {
                break 'outer;
            }
            normal = cross(normalized(v[a]), normalized(v[b]));
        }
        if norm(normal) >= 1e-6 {
            break;
        }
    }
    if norm(normal) < 1e-6 {
        return true;
    }
    let normal = normalized(normal);
    v.iter().all(|x| dot(normal, normalized(*x)).abs() <= 0.25)
}

pub fn polyhedral_template_matching(p: &Parameters) -> R {
    let m = match_all(p, &enabled_kinds(p)?, rmsd_cutoff_of(p)?)?;
    let n = m.fits.len();
    let mut types = vec![0i32; n];
    let (mut rmsd, mut lattice, mut angle) = (vec![f64::NAN; n], vec![f64::NAN; n], vec![f64::NAN; n]);
    let mut orientations = Vec::new();
    for i in 0..n {
        let fit = &m.fits[i];
        types[i] = fit.kind;
        if fit.rmsd.is_finite() {
            rmsd[i] = fit.rmsd;
        }
        if fit.kind == OTHER {
            orientations.push(jv(&[f64::NAN; 4]));
            continue;
        }
        lattice[i] = fit.scale;
        let q = reduced_orientation(fit.kind, &fit.rotation);
        angle[i] = 2.0 * (if 1.0 < q.w { 1.0 } else { q.w }).acos() * 180.0 / PI;
        orientations.push(quaternion_json(&q));
    }
    let mut result = Json::object();
    add_classification(&mut result, "structure_type", &types, &NAMES, "");
    result.set("rmsd", jv(&rmsd));
    result.set("orientation", Json::Array(orientations));
    result.set("orientation_angle_deg", jv(&angle));
    result.set("lattice_constant_A", jv(&lattice));
    result.set("cutoff_A", m.cutoff);
    if p.has("types") {
        let species = finite_array(p.array("types")?, "types", 1)?.values;
        if species.len() != n {
            return Err("types needs one value per atom".into());
        }
        let order_names = ["none", "pure", "L1_0", "L1_2 (minority site)", "L1_2 (majority site)", "B2", "zincblende or wurtzite", "other"];
        let mut ordering = vec![0i32; n];
        let (mut env, mut atoms) = (Vec::new(), Vec::new());
        for i in 0..n {
            let kind = m.fits[i].kind;
            if kind == OTHER || !environment(kind, i, &m.neighbors, &mut env, &mut atoms) {
                continue;
            }
            let first = if kind == CUBIC_DIAMOND || kind == HEX_DIAMOND {
                4
            } else if kind == BCC {
                8
            } else {
                env.len()
            };
            let (mut same, mut unlike) = (Vec::new(), Vec::new());
            let mut outer_same = true;
            for k in 0..env.len() {
                let equal = species[atoms[k] as usize] == species[i];
                if k < first {
                    if equal {
                        same.push(env[k]);
                    } else {
                        unlike.push(env[k]);
                    }
                } else {
                    outer_same = outer_same && equal;
                }
            }
            ordering[i] = if unlike.is_empty() && outer_same {
                1
            } else if kind == FCC && unlike.len() == 12 {
                3
            } else if kind == FCC && unlike.len() == 4 && coplanar(&unlike) {
                4
            } else if kind == FCC && same.len() == 4 && coplanar(&same) {
                2
            } else if kind == BCC && same.is_empty() && outer_same {
                5
            } else if (kind == CUBIC_DIAMOND || kind == HEX_DIAMOND) && same.is_empty() && outer_same {
                6
            } else {
                7
            };
        }
        add_classification(&mut result, "ordering", &ordering, &order_names, "ordering_");
    }
    Ok(result)
}

pub fn grain_segmentation(p: &Parameters) -> R {
    let threshold = positive(p.number_or("merge_threshold_deg", 4.0)?, "merge_threshold_deg")? * PI / 180.0;
    let minimum = integer(p.number_or("min_grain_size", 100.0)?, "min_grain_size", 1)?;
    let adopt = p.boolean("adopt_orphans", true)?;
    let m = match_all(p, &enabled_kinds(p)?, rmsd_cutoff_of(p)?)?;
    let n = m.fits.len();
    let orientation: Vec<Quat> = (0..n).map(|i| if m.fits[i].kind != OTHER { reduced_orientation(m.fits[i].kind, &m.fits[i].rotation) } else { Quat::IDENTITY }).collect();
    let mut edges: Vec<(f64, usize, usize)> = Vec::new();
    let (mut env, mut atoms) = (Vec::new(), Vec::new());
    for i in 0..n {
        let kind = m.fits[i].kind;
        if kind == OTHER || !environment(kind, i, &m.neighbors, &mut env, &mut atoms) {
            continue;
        }
        let first = if kind == CUBIC_DIAMOND || kind == HEX_DIAMOND { 4 } else { atoms.len() };
        for &a in &atoms[..first] {
            let j = a as usize;
            if j <= i || m.fits[j].kind != kind {
                continue;
            }
            let angle = misorientation(kind, &orientation[i], &orientation[j]);
            if angle < threshold {
                edges.push((angle, i, j));
            }
        }
    }
    sort_by(&mut edges, |a, b| a.0 < b.0);
    let mut parent: Vec<usize> = (0..n).collect();
    let mut size = vec![1usize; n];
    let mut sum: Vec<[f64; 4]> = orientation.iter().map(|q| q.coeffs()).collect();
    let mean_of = |sum: &[[f64; 4]], r: usize| {
        let c = sum[r];
        let norm = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2] + c[3] * c[3]).sqrt();
        Quat::from_coeffs(c.map(|v| v / norm))
    };
    for &(_, ea, eb) in &edges {
        let mut a = union_root(&mut parent, ea);
        let mut b = union_root(&mut parent, eb);
        if a == b {
            continue;
        }
        let kind = m.fits[ea].kind;
        let (qa, qb) = (mean_of(&sum, a), mean_of(&sum, b));
        if misorientation(kind, &qa, &qb) >= threshold {
            continue;
        }
        if size[a] < size[b] {
            std::mem::swap(&mut a, &mut b);
        }
        // qa is the mean before the swap, as in the C++.
        let delta = qa.conjugate().mul(mean_of(&sum, b));
        let (mut best_s, mut best_w) = (Quat::IDENTITY, -1.0);
        for s in lattice_symmetry(kind) {
            let w = delta.mul(*s).w.abs();
            if w > best_w {
                best_w = w;
                best_s = *s;
            }
        }
        let aligned = Quat::from_coeffs(sum[b]).mul(best_s);
        let mut addition = aligned.coeffs();
        let d: f64 = (0..4).map(|k| addition[k] * sum[a][k]).sum();
        if d < 0.0 {
            addition = addition.map(|v| -v);
        }
        for k in 0..4 {
            sum[a][k] += addition[k];
        }
        size[a] += size[b];
        parent[b] = a;
    }
    let mut clusters: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        if m.fits[i].kind != OTHER {
            let root = union_root(&mut parent, i);
            clusters.entry(root).or_default().push(i);
        }
    }
    let mut grains: Vec<(usize, Vec<usize>)> = clusters.into_iter().filter(|(_, members)| members.len() as i64 >= minimum).collect();
    sort_by(&mut grains, |a, b| a.1.len() > b.1.len());
    let mut grain = vec![0i32; n];
    for (g, (_, members)) in grains.iter().enumerate() {
        for &i in members {
            grain[i] = g as i32 + 1;
        }
    }
    let mut adopted = 0usize;
    if adopt && !grains.is_empty() {
        loop {
            let mut claims = Vec::new();
            for i in 0..n {
                if grain[i] != 0 {
                    continue;
                }
                if let Some(&j) = m.neighbors.index[i].iter().find(|&&j| grain[j as usize] != 0) {
                    claims.push((i, grain[j as usize]));
                }
            }
            if claims.is_empty() {
                break;
            }
            for &(i, g) in &claims {
                grain[i] = g;
            }
            adopted += claims.len();
        }
    }
    let mut atoms_per_grain = vec![0i64; grains.len() + 1];
    for &g in &grain {
        atoms_per_grain[g as usize] += 1;
    }
    let table: Vec<Json> = grains
        .iter()
        .enumerate()
        .map(|(g, (r, members))| {
            let kind = m.fits[members[0]].kind;
            let mut item = Json::object();
            item.set("id", (g + 1) as f64);
            item.set("atoms", atoms_per_grain[g + 1] as f64);
            item.set("crystalline_atoms", members.len() as f64);
            item.set("structure", NAMES[kind as usize]);
            item.set("orientation", quaternion_json(&reduced_orientation(kind, &mean_of(&sum, *r).to_matrix())));
            item
        })
        .collect();
    let mut contacts: BTreeMap<(i32, i32), i64> = BTreeMap::new();
    for i in 0..n {
        for &j in &m.neighbors.index[i] {
            let gj = grain[j as usize];
            if grain[i] != 0 && gj != 0 && grain[i] < gj {
                *contacts.entry((grain[i], gj)).or_insert(0) += 1;
            }
        }
    }
    let boundaries: Vec<Json> = contacts
        .iter()
        .map(|(&(ga, gb), &count)| {
            let (gra, grb) = (&grains[(ga - 1) as usize], &grains[(gb - 1) as usize]);
            let (ka, kb) = (m.fits[gra.1[0]].kind, m.fits[grb.1[0]].kind);
            let mut item = Json::object();
            item.set("grains", jv(&[ga as f64, gb as f64]));
            item.set("neighbor_pairs", count as f64);
            item.set("misorientation_deg", if ka == kb { misorientation(ka, &mean_of(&sum, gra.0), &mean_of(&sum, grb.0)) * 180.0 / PI } else { f64::NAN });
            item
        })
        .collect();
    let types: Vec<i32> = m.fits.iter().map(|f| f.kind).collect();
    let mut result = Json::object();
    result.set("grain_id", ji(&grain));
    result.set("structure_type", ji(&types));
    result.set("grain_count", grains.len() as f64);
    result.set("grains", Json::Array(table));
    result.set("boundaries", Json::Array(boundaries));
    result.set("adopted_atoms", adopted as f64);
    result.set("unassigned_atoms", atoms_per_grain[0] as f64);
    result.set("cutoff_A", m.cutoff);
    Ok(result)
}

pub fn elastic_strain(p: &Parameters) -> R {
    let lattice = positive(p.number("lattice_constant_A")?, "lattice_constant_A")?;
    let ca = positive(p.number_or("c_over_a", ideal_ca())?, "c_over_a")?;
    let kind = kind_from_key(if p.has("lattice") { p.json("lattice")?.string()? } else { "fcc" })?;
    if kind == ICO {
        return Err("Elastic strain needs a crystal lattice, not icosahedral order".into());
    }
    let m = match_all(p, &[kind], rmsd_cutoff_of(p)?)?;
    let reference = template_vectors(kind, ca);
    let n = m.fits.len();
    let (mut volumetric, mut shear) = (vec![f64::NAN; n], vec![f64::NAN; n]);
    let (mut strain, mut gradient) = (Vec::new(), Vec::new());
    let (mut env, mut atoms) = (Vec::new(), Vec::new());
    let mut identified = 0usize;
    let mut volumetric_sum = 0.0;
    for i in 0..n {
        let fit = &m.fits[i];
        if fit.kind != kind || !environment(kind, i, &m.neighbors, &mut env, &mut atoms) {
            strain.push(jv(&[f64::NAN; 6]));
            gradient.push(jv(&[f64::NAN; 9]));
            continue;
        }
        let (mut bt, mut tt) = ([[0.0; 3]; 3], [[0.0; 3]; 3]);
        for k in 0..env.len() {
            let t = scale(reference[fit.map[k] as usize], lattice);
            bt = madd(&bt, &outer(env[k], t));
            tt = madd(&tt, &outer(t, t));
        }
        let f = matmul3(&bt, &inverse3(&tt));
        let ff = matmul3(&f, &transpose3(&f));
        let e: M3 = std::array::from_fn(|r| std::array::from_fn(|c| 0.5 * (ff[r][c] - IDENTITY3[r][c])));
        strain.push(jv(&[e[0][0], e[1][1], e[2][2], e[1][2], e[0][2], e[0][1]]));
        let local = matmul3(&transpose3(&fit.rotation), &f);
        gradient.push(jv(&[local[0][0], local[0][1], local[0][2], local[1][0], local[1][1], local[1][2], local[2][0], local[2][1], local[2][2]]));
        volumetric[i] = (e[0][0] + e[1][1] + e[2][2]) / 3.0;
        shear[i] = (e[0][1] * e[0][1]
            + e[1][2] * e[1][2]
            + e[0][2] * e[0][2]
            + ((e[0][0] - e[1][1]) * (e[0][0] - e[1][1]) + (e[1][1] - e[2][2]) * (e[1][1] - e[2][2]) + (e[0][0] - e[2][2]) * (e[0][0] - e[2][2])) / 6.0)
            .sqrt();
        identified += 1;
        volumetric_sum += volumetric[i];
    }
    let mut result = Json::object();
    result.set("strain_voigt", Json::Array(strain));
    result.set("strain_components", Json::Array(["xx", "yy", "zz", "yz", "xz", "xy"].iter().map(|s| Json::String(s.to_string())).collect()));
    result.set("volumetric_strain", jv(&volumetric));
    result.set("shear_strain", jv(&shear));
    result.set("deformation_gradient_lattice", Json::Array(gradient));
    result.set("structure", NAMES[kind as usize]);
    result.set("identified_atoms", identified as f64);
    result.set("mean_volumetric_strain", if identified > 0 { volumetric_sum / identified as f64 } else { f64::NAN });
    result.set("cutoff_A", m.cutoff);
    Ok(result)
}

pub fn fcc_planar_faults(p: &Parameters) -> R {
    let m = match_all(p, &[FCC, HCP], rmsd_cutoff_of(p)?)?;
    let n = m.fits.len();
    let names = ["non-hcp", "intrinsic stacking fault", "coherent twin boundary", "multilayer stacking fault", "other hcp"];
    let mut fault = vec![0i32; n];
    let (mut env, mut atoms) = (Vec::new(), Vec::new());
    let (mut spacing, mut spacing_count) = (0.0, 0usize);
    for i in 0..n {
        if m.fits[i].kind == FCC || m.fits[i].kind == HCP {
            spacing += if m.fits[i].kind == FCC { m.fits[i].scale / 2.0f64.sqrt() } else { m.fits[i].scale };
            spacing_count += 1;
        }
        if m.fits[i].kind != HCP || !environment(HCP, i, &m.neighbors, &mut env, &mut atoms) {
            continue;
        }
        let normal = mul3(&m.fits[i].rotation, [0.0, 0.0, 1.0]);
        let (mut up, mut down) = (0, 0);
        for k in 0..env.len() {
            let c = dot(normalized(env[k]), normal);
            if m.fits[atoms[k] as usize].kind != HCP {
                continue;
            }
            if c > 0.5 {
                up += 1;
            } else if c < -0.5 {
                down += 1;
            }
        }
        fault[i] = match (up, down) {
            (0, 0) => 2,
            (3, 0) | (0, 3) => 1,
            (3, 3) => 3,
            _ => 4,
        };
    }
    let mut count = vec![0i64; names.len()];
    for &f in &fault {
        count[f as usize] += 1;
    }
    let d = if spacing_count > 0 { spacing / spacing_count as f64 } else { f64::NAN };
    let area_per_atom = 3.0f64.sqrt() / 2.0 * d * d;
    let mut result = Json::object();
    add_classification(&mut result, "fault_type", &fault, &names, "");
    result.set("nearest_neighbor_distance_A", d);
    result.set("stacking_fault_area_A2", count[1] as f64 / 2.0 * area_per_atom);
    result.set("twin_boundary_area_A2", count[2] as f64 * area_per_atom);
    result.set("cutoff_A", m.cutoff);
    Ok(result)
}
