//! Coincidence-site-lattice grain boundaries (the C++
//! algorithms/CSLComputation), following the aimsgb workflow (Cheng, Luo and
//! Yang, Comput. Mater. Sci. 155 (2018) 92): Σ candidates for a rotation
//! axis, CSL matrices, the two grains and the stacked bicrystal.

use std::collections::BTreeMap;

use af_core::cxx::sort_by;
use af_core::{AtomSite, Structure};

const PI: f64 = 3.14159265358979323846;

pub type M3 = [[f64; 3]; 3];

fn gcd2(a: i32, b: i32) -> i32 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    if a == 0 {
        1
    } else {
        a
    }
}

fn gcd_n(v: &[i32]) -> i32 {
    let mut g = 0;
    for &x in v {
        g = gcd2(g, x.abs());
    }
    if g == 0 {
        1
    } else {
        g
    }
}

fn is_integer(x: f64, tol: f64) -> bool {
    (x - x.round()).abs() < tol
}

fn is_integer_vec(v: &[f64], tol: f64) -> bool {
    let mut sum = 0.0;
    for &x in v {
        sum += (x - x.round()) * (x - x.round());
    }
    sum.sqrt() < tol
}

fn smallest_multiplier(v: &[f64]) -> i32 {
    for m in 1..=10000 {
        if v.iter().all(|&x| is_integer(m as f64 * x, 1e-5)) {
            return m;
        }
    }
    1
}

pub fn eye3() -> M3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

pub fn transpose3(a: &M3) -> M3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}

pub fn mul3(a: &M3, b: &M3) -> M3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

pub fn det3(m: &M3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

pub fn inv3(m: &M3) -> M3 {
    let d = det3(m);
    if d.abs() < 1e-15 {
        return eye3();
    }
    let id = 1.0 / d;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * id,
            -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * id,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id,
        ],
        [
            -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * id,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id,
            -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * id,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * id,
            -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * id,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id,
        ],
    ]
}

pub fn round_m3(m: &M3, decimals: i32) -> M3 {
    let mut scale = 1.0;
    for _ in 0..decimals {
        scale *= 10.0;
    }
    m.map(|row| row.map(|x| (x * scale).round() / scale))
}

fn round_to_int(m: &M3) -> [[i32; 3]; 3] {
    m.map(|row| row.map(|x| x.round() as i32))
}

fn from_int(a: &[[i32; 3]; 3]) -> M3 {
    a.map(|row| row.map(|x| x as f64))
}

fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm3(v: &[f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Cramer's rule (zero for a singular matrix).
pub fn solve3x3(a: &M3, b: &[f64; 3]) -> [f64; 3] {
    let d = det3(a);
    if d.abs() < 1e-15 {
        return [0.0; 3];
    }
    std::array::from_fn(|col| {
        let mut ai = *a;
        for row in 0..3 {
            ai[row][col] = b[row];
        }
        det3(&ai) / d
    })
}

pub fn reduce_int_vec(v: &mut [i32]) {
    let g = gcd_n(v);
    if g > 1 {
        for x in v.iter_mut() {
            *x /= g;
        }
    }
}

fn to_integer_vec3(v: &[f64; 3]) -> [i32; 3] {
    let mult = smallest_multiplier(v) as f64;
    let mut out = v.map(|x| (x * mult).round() as i32);
    reduce_int_vec(&mut out);
    out
}

const UNIMODULAR: [[[i32; 3]; 3]; 5] = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 1], [0, 1, 0], [0, 1, 1]],
    [[1, 0, 1], [0, 1, 0], [0, 1, -1]],
    [[1, 0, 1], [0, 1, 0], [-1, 1, 0]],
    [[1, 0, 1], [1, 1, 0], [1, 1, 1]],
];

/// Rodrigues rotation about an integer axis by degrees.
pub fn rotate_matrix(axis: &[i32; 3], angle_deg: f64) -> M3 {
    let mut ax = axis.map(|x| x as f64);
    let n = norm3(&ax);
    if n < 1e-12 {
        return eye3();
    }
    for x in ax.iter_mut() {
        *x /= n;
    }
    let angle = angle_deg * PI / 180.0;
    let (c, s) = (angle.cos(), angle.sin());
    let mut k = [[0.0; 3]; 3];
    k[0][1] = -ax[2];
    k[0][2] = ax[1];
    k[1][0] = ax[2];
    k[1][2] = -ax[0];
    k[2][0] = -ax[1];
    k[2][1] = ax[0];
    let kk = mul3(&k, &k);
    let mut r = eye3();
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] += k[i][j] * s + kk[i][j] * (1.0 - c);
        }
    }
    r
}

/// (|sum| max, |orig| max, |sum| sum, |orig| sum) for row i plus sign * row j.
fn shrink_measure(csl: &M3, i: usize, j: usize, sign: f64) -> bool {
    let (mut sum_max, mut orig_max, mut sum_sum, mut orig_sum) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for c in 0..3 {
        let s = (csl[i][c] + sign * csl[j][c]).abs();
        let o = csl[i][c].abs();
        if s > sum_max {
            sum_max = s;
        }
        if o > orig_max {
            orig_max = o;
        }
        sum_sum += s;
        orig_sum += o;
    }
    sum_max < orig_max || (sum_max == orig_max && sum_sum < orig_sum)
}

pub fn o_lattice_to_csl(o_lattice: &M3, mut n: f64) -> M3 {
    let mut csl = transpose3(o_lattice);
    if n < 0.0 {
        for j in 0..3 {
            csl[0][j] *= -1.0;
        }
        n *= -1.0;
    }
    let n_int = n.round() as i32;
    for _ in 0..200 {
        let m: [i32; 3] = std::array::from_fn(|i| smallest_multiplier(&csl[i]));
        let m_prod = m[0] * m[1] * m[2];
        if m_prod <= n_int {
            for i in 0..3 {
                for j in 0..3 {
                    csl[i][j] *= m[i] as f64;
                }
            }
            if m_prod < n_int {
                let factor = n_int / m_prod;
                for j in 0..3 {
                    csl[0][j] *= factor as f64;
                }
            }
            break;
        }
        let mut changed = false;
        'search: for i in 0..3 {
            for j in 0..3 {
                if i == j || m[i] == 1 || m[j] == 1 {
                    continue;
                }
                let (a, b) = if m[i] <= m[j] { (i, j) } else { (j, i) };
                for ki in 1..m[b] {
                    for kk in [ki, -ki] {
                        let handle: [f64; 3] = std::array::from_fn(|c| csl[a][c] + kk as f64 * csl[b][c]);
                        if smallest_multiplier(&handle) < m[a] {
                            for c in 0..3 {
                                csl[a][c] += kk as f64 * csl[b][c];
                            }
                            changed = true;
                            break 'search;
                        }
                    }
                }
            }
        }
        if !changed {
            for i in 0..3 {
                for j in 0..3 {
                    csl[i][j] *= m[i] as f64;
                }
            }
            break;
        }
    }
    csl = csl.map(|row| row.map(|x| x.round()));
    for _ in 0..200 {
        let mut changed = false;
        'outer: for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    continue;
                }
                while shrink_measure(&csl, i, j, 1.0) {
                    for c in 0..3 {
                        csl[i][c] += csl[j][c];
                    }
                    changed = true;
                }
                if !changed {
                    while shrink_measure(&csl, i, j, -1.0) {
                        for c in 0..3 {
                            csl[i][c] -= csl[j][c];
                        }
                        changed = true;
                    }
                }
                if changed {
                    break 'outer;
                }
            }
        }
        if !changed {
            break;
        }
    }
    transpose3(&csl)
}

pub fn reduce_csl(csl: &M3) -> M3 {
    let mut t = transpose3(csl);
    for row in t.iter_mut() {
        let mut iv = row.map(|x| x.round() as i32);
        reduce_int_vec(&mut iv);
        *row = iv.map(|x| x as f64);
    }
    transpose3(&t)
}

pub fn orthogonalize_csl(csl_in: &M3, axis: &[i32; 3]) -> M3 {
    let mut csl = transpose3(csl_in);
    let axis_d = axis.map(|x| x as f64);
    let mut c = solve3x3(&transpose3(&csl), &axis_d);
    if !is_integer_vec(&c, 1e-5) {
        let mult = smallest_multiplier(&c) as f64;
        for x in c.iter_mut() {
            *x *= mult;
        }
    }
    c = c.map(|x| x.round());
    let mut ind: i32 = -1;
    let mut min_abs = 1e30;
    for i in 0..3 {
        if c[i].abs() > 1e-8 && c[i].abs() < min_abs {
            min_abs = c[i].abs();
            ind = i as i32;
        }
    }
    let ind = if ind < 0 { 2 } else { ind as usize };
    if ind != 2 {
        let tmp = csl[2];
        csl[2] = csl[ind].map(|x| -x);
        csl[ind] = tmp;
        let ct = c[2];
        c[2] = -c[ind];
        c[ind] = ct;
    }
    let mut new_row2 = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            new_row2[j] += c[i] * csl[i][j];
        }
    }
    csl[2] = new_row2;
    if c[2] < 0.0 {
        csl[1] = csl[1].map(|x| x * -1.0);
    }
    let n1 = norm3(&csl[2]);
    if n1 < 1e-12 {
        return transpose3(&csl);
    }
    let u1 = csl[2].map(|x| x / n1);
    let d1 = dot3(&csl[1], &u1);
    let d0 = dot3(&csl[0], &u1);
    let y2_1: [f64; 3] = std::array::from_fn(|i| csl[1][i] - d1 * u1[i]);
    let y2_2: [f64; 3] = std::array::from_fn(|i| csl[0][i] - d0 * u1[i]);
    let c0_1 = to_integer_vec3(&y2_1);
    let c0_2 = to_integer_vec3(&y2_2);
    let sum1 = c0_1[0].abs() + c0_1[1].abs() + c0_1[2].abs();
    let sum2 = c0_2[0].abs() + c0_2[1].abs() + c0_2[2].abs();
    if sum1 > sum2 {
        let n2 = norm3(&y2_2);
        if n2 < 1e-12 {
            return transpose3(&csl);
        }
        let u2 = y2_2.map(|x| x / n2);
        let (a, b) = (dot3(&csl[1], &u1), dot3(&csl[1], &u2));
        let y3: [f64; 3] = std::array::from_fn(|i| csl[1][i] - a * u1[i] - b * u2[i]);
        csl[1] = to_integer_vec3(&y3).map(|x| x as f64);
        csl[0] = c0_2.map(|x| x as f64);
    } else {
        let n2 = norm3(&y2_1);
        if n2 < 1e-12 {
            return transpose3(&csl);
        }
        let u2 = y2_1.map(|x| x / n2);
        let (a, b) = (dot3(&csl[0], &u1), dot3(&csl[0], &u2));
        let y3: [f64; 3] = std::array::from_fn(|i| csl[0][i] - a * u1[i] - b * u2[i]);
        let y3i = to_integer_vec3(&y3);
        csl[1] = c0_1.map(|x| x as f64);
        csl[0] = y3i.map(|x| x as f64);
    }
    csl = csl.map(|row| row.map(|x| x.round()));
    transpose3(&csl)
}

pub fn csl_matrix(sigma: i32, rotation: &M3) -> M3 {
    let r_inv = inv3(rotation);
    let mut t = eye3();
    for u in &UNIMODULAR {
        let ur = mul3(&from_int(u), &r_inv);
        t = eye3();
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] -= ur[i][j];
            }
        }
        if det3(&t).abs() > 1e-6 {
            break;
        }
    }
    let o_lattice = round_m3(&inv3(&t), 12);
    let n = (sigma as f64 / det3(&o_lattice) * 1e6).round() / 1e6;
    o_lattice_to_csl(&o_lattice, n)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SigmaCandidate {
    pub sigma: i32,
    pub m: i32,
    pub n: i32,
    pub theta_deg: f32,
    pub csl: [[i32; 3]; 3],
    /// The three boundary planes: the columns of the CSL.
    pub plane: [[i32; 3]; 3],
}

pub fn classify_boundary_type(axis: &[i32; 3], plane: &[i32; 3]) -> &'static str {
    let dot = axis[0] * plane[0] + axis[1] * plane[1] + axis[2] * plane[2];
    let cx = axis[1] * plane[2] - axis[2] * plane[1];
    let cy = axis[2] * plane[0] - axis[0] * plane[2];
    let cz = axis[0] * plane[1] - axis[1] * plane[0];
    if cx == 0 && cy == 0 && cz == 0 {
        return "Twist";
    }
    if dot == 0 {
        return "Tilt";
    }
    "Mixed"
}

fn fill_candidate(sigma: i32, m: i32, n: i32, theta: f32, csl: &[[i32; 3]; 3]) -> SigmaCandidate {
    SigmaCandidate { sigma, m, n, theta_deg: theta, csl: *csl, plane: std::array::from_fn(|col| [csl[0][col], csl[1][col], csl[2][col]]) }
}

/// Σ candidates up to `max_sigma` for a rotation axis (two per Σ for principal axes).
pub fn compute_gb_info(axis_in: &[i32; 3], max_sigma: i32) -> Vec<SigmaCandidate> {
    let mut axis = *axis_in;
    reduce_int_vec(&mut axis);
    let axis_norm2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    if axis_norm2 <= 0 {
        return Vec::new();
    }
    let axis_norm = (axis_norm2 as f32).sqrt();
    let mut sorted = axis.map(|x| x.abs());
    sorted.sort();
    let principal = sorted == [0, 0, 1];
    #[derive(Clone, Copy)]
    struct Tmn {
        theta: f64,
        m: i32,
        n: i32,
    }
    let mut sigma_tmn: BTreeMap<i32, Vec<Tmn>> = BTreeMap::new();
    let m_max = (4.0 * max_sigma as f64).sqrt().ceil() as i32;
    for m in 0..=m_max {
        for n in 0..=m_max {
            if gcd2(m, n) > 1 {
                continue;
            }
            let mut sigma = m * m + n * n * axis_norm2;
            while sigma != 0 && sigma % 2 == 0 {
                sigma /= 2;
            }
            if sigma <= 1 || sigma > max_sigma {
                continue;
            }
            let theta = if m == 0 { 180.0 } else { 2.0 * (n as f64 * axis_norm as f64 / m as f64).atan() * 180.0 / PI };
            sigma_tmn.entry(sigma).or_default().push(Tmn { theta, m, n });
        }
    }
    let is_std_basis_row = |v: &[i32; 3]| v[0] * v[0] + v[1] * v[1] + v[2] * v[2] == 1 && v[0] + v[1] + v[2] == 1;
    let mut result = Vec::new();
    for (&sigma, entries) in sigma_tmn.iter_mut() {
        sort_by(entries, |a, b| a.theta < b.theta);
        let min_theta = entries[0].theta;
        let rotation = rotate_matrix(&axis, min_theta);
        let csl = round_to_int(&reduce_csl(&orthogonalize_csl(&csl_matrix(sigma, &rotation), &axis)));
        if principal && entries.len() > 1 {
            // A second CSL from the rows with m and n exchanged (aimsgb ext_csl).
            let mut found: Option<(usize, [[i32; 3]; 3])> = None;
            for (ei, e) in entries.iter().enumerate() {
                let mut tried = csl;
                for row in 0..3 {
                    if is_std_basis_row(&csl[row]) {
                        continue;
                    }
                    let (a, b) = (csl[row][0], csl[row][1]);
                    let (sa, sb) = (if a >= 0 { 1 } else { -1 }, if b >= 0 { 1 } else { -1 });
                    tried[row] = if a.abs() > b.abs() { [sa * e.m, sb * e.n, 0] } else { [sa * e.n, sb * e.m, 0] };
                }
                if tried != csl {
                    found = Some((ei, tried));
                    break;
                }
            }
            if let Some((ind, ext)) = found {
                let theta2 = 90.0f32 - min_theta as f32;
                let (first, second) = if ind > 0 {
                    (
                        fill_candidate(sigma, entries[0].m, entries[0].n, min_theta as f32, &csl),
                        fill_candidate(sigma, entries[ind].m, entries[ind].n, theta2, &ext),
                    )
                } else {
                    (
                        fill_candidate(sigma, entries[0].m, entries[0].n, theta2, &csl),
                        fill_candidate(sigma, entries[ind].m, entries[ind].n, min_theta as f32, &ext),
                    )
                };
                result.push(first);
                result.push(second);
                continue;
            }
        }
        result.push(fill_candidate(sigma, entries[0].m, entries[0].n, min_theta as f32, &csl));
    }
    result
}

/// Atoms and cell vectors (rows) of one grain.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Grain {
    pub atoms: Vec<AtomSite>,
    pub cell: M3,
}

pub fn structure_to_grain(s: &Structure) -> Grain {
    Grain { atoms: s.atoms.clone(), cell: s.cell_vectors }
}

pub fn grain_to_structure(g: &Grain) -> Structure {
    Structure { atoms: g.atoms.clone(), has_unit_cell: true, cell_vectors: g.cell, cell_offset: [0.0; 3], ..Default::default() }
}

pub fn cart_to_frac(cart: &[f64; 3], cell_inv: &M3) -> [f64; 3] {
    let mut frac = [0.0; 3];
    for j in 0..3 {
        for i in 0..3 {
            frac[j] += cart[i] * cell_inv[i][j];
        }
    }
    frac
}

pub fn frac_to_cart(frac: &[f64; 3], cell: &M3) -> [f64; 3] {
    let mut cart = [0.0; 3];
    for j in 0..3 {
        for i in 0..3 {
            cart[j] += frac[i] * cell[i][j];
        }
    }
    cart
}

pub fn wrap_frac(mut f: f64) -> f64 {
    f -= f.floor();
    if f >= 1.0 - 1e-8 {
        f = 0.0;
    }
    if f < -1e-10 {
        f += 1.0;
    }
    f
}

fn clamp_unit(x: f64) -> f64 {
    // std::max(-1.0, std::min(1.0, x))
    let m = if x < 1.0 { x } else { 1.0 };
    if -1.0 < m {
        m
    } else {
        -1.0
    }
}

pub fn cell_to_parameters(cell: &M3) -> [f64; 6] {
    let (a, b, c) = (norm3(&cell[0]), norm3(&cell[1]), norm3(&cell[2]));
    let cos_alpha = if a > 1e-12 && b > 1e-12 && c > 1e-12 { dot3(&cell[1], &cell[2]) / (b * c) } else { 0.0 };
    let cos_beta = if a > 1e-12 && c > 1e-12 { dot3(&cell[0], &cell[2]) / (a * c) } else { 0.0 };
    let cos_gamma = if a > 1e-12 && b > 1e-12 { dot3(&cell[0], &cell[1]) / (a * b) } else { 0.0 };
    let deg = |x: f64| clamp_unit(x).acos() * 180.0 / PI;
    [a, b, c, deg(cos_alpha), deg(cos_beta), deg(cos_gamma)]
}

pub fn parameters_to_cell(p: &[f64; 6]) -> M3 {
    let [a, b, c, alpha, beta, gamma] = *p;
    let (ar, br, gr) = (alpha * PI / 180.0, beta * PI / 180.0, gamma * PI / 180.0);
    let (cos_alpha, cos_beta, cos_gamma, sin_gamma) = (ar.cos(), br.cos(), gr.cos(), gr.sin());
    let mut cell = [[0.0; 3]; 3];
    cell[0][0] = a;
    cell[1][0] = b * cos_gamma;
    cell[1][1] = b * sin_gamma;
    cell[2][0] = c * cos_beta;
    let cy = if sin_gamma.abs() > 1e-10 { c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma } else { 0.0 };
    cell[2][1] = cy;
    let cz2 = c * c - cell[2][0] * cell[2][0] - cy * cy;
    cell[2][2] = if cz2 > 0.0 { cz2.sqrt() } else { 0.0 };
    cell
}

fn remap(atoms: &[AtomSite], from_inv: &M3, to_cell: &M3) -> Vec<AtomSite> {
    atoms
        .iter()
        .map(|a| {
            let frac = cart_to_frac(&[a.x, a.y, a.z], from_inv).map(wrap_frac);
            let p = frac_to_cart(&frac, to_cell);
            let mut b = a.clone();
            b.x = p[0];
            b.y = p[1];
            b.z = p[2];
            b
        })
        .collect()
}

/// The cell in standard orientation (a along x, b in xy) with the atoms re-placed.
pub fn canonicalize_grain(g: &mut Grain) {
    let old_inv = inv3(&g.cell);
    let new_cell = parameters_to_cell(&cell_to_parameters(&g.cell));
    g.atoms = remap(&g.atoms, &old_inv, &new_cell);
    g.cell = new_cell;
}

/// The supercell with rows of M (integers) in terms of the grain's vectors.
pub fn make_supercell(grain: &Grain, m: &[[i32; 3]; 3]) -> Grain {
    let mut result = Grain::default();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                result.cell[i][j] += m[i][k] as f64 * grain.cell[k][j];
            }
        }
    }
    let m_inv = inv3(&from_int(m));
    let cell_inv = inv3(&grain.cell);
    let range = m.iter().map(|row| row.iter().map(|x| x.abs()).sum::<i32>()).max().unwrap_or(0) + 1;
    for atom in &grain.atoms {
        let old = cart_to_frac(&[atom.x, atom.y, atom.z], &cell_inv);
        for ti in -range..=range {
            for tj in -range..=range {
                for tk in -range..=range {
                    let shifted = [old[0] + ti as f64, old[1] + tj as f64, old[2] + tk as f64];
                    let mut f = [0.0; 3];
                    for j in 0..3 {
                        for k in 0..3 {
                            f[j] += shifted[k] * m_inv[k][j];
                        }
                    }
                    if f.iter().any(|&x| x < -1e-8 || x >= 1.0 - 1e-8) {
                        continue;
                    }
                    let p = frac_to_cart(&f.map(wrap_frac), &result.cell);
                    let mut a = atom.clone();
                    a.x = p[0];
                    a.y = p[1];
                    a.z = p[2];
                    result.atoms.push(a);
                }
            }
        }
    }
    canonicalize_grain(&mut result);
    result
}

pub fn make_supercell_diag(grain: &Grain, na: i32, nb: i32, nc: i32) -> Grain {
    make_supercell(grain, &[[na, 0, 0], [0, nb, 0], [0, 0, nc]])
}

pub fn is_orthogonal(cell: &M3, tol: f64) -> bool {
    let d = |a: usize, b: usize| dot3(&cell[a], &cell[b]);
    let (l0, l1, l2) = (norm3(&cell[0]), norm3(&cell[1]), norm3(&cell[2]));
    if l0 < 1e-10 || l1 < 1e-10 || l2 < 1e-10 {
        return true;
    }
    d(0, 1).abs() / (l0 * l1) < tol && d(0, 2).abs() / (l0 * l2) < tol && d(1, 2).abs() / (l1 * l2) < tol
}

/// The `direction` vector replaced by its component normal to the other two.
pub fn set_orthogonal_grain(g: &Grain, direction: usize) -> Grain {
    let mut result = Grain { cell: g.cell, atoms: Vec::new() };
    let (v1, v2) = (result.cell[(direction + 1) % 3], result.cell[(direction + 2) % 3]);
    let mut cross = [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]];
    let cn = norm3(&cross);
    if cn < 1e-12 {
        result.atoms = g.atoms.clone();
        return result;
    }
    for x in cross.iter_mut() {
        *x /= cn;
    }
    let proj = dot3(&g.cell[direction], &cross);
    result.cell[direction] = cross.map(|x| proj * x);
    result.atoms = remap(&g.atoms, &inv3(&result.cell), &result.cell);
    result
}

/// Grain B: grain A rotated 180 degrees about the first axis that changes it.
pub fn b_from_a(a: &Grain) -> Grain {
    let cell_inv = inv3(&a.cell);
    for try_axis in 0..3 {
        let mut b = Grain { cell: a.cell, atoms: Vec::new() };
        let mut different = false;
        for atom in &a.atoms {
            let mut cart = [atom.x, atom.y, atom.z];
            for (c, x) in cart.iter_mut().enumerate() {
                if c != try_axis {
                    *x = -*x;
                }
            }
            let frac = cart_to_frac(&cart, &cell_inv).map(wrap_frac);
            let p = frac_to_cart(&frac, &a.cell);
            let orig = cart_to_frac(&[atom.x, atom.y, atom.z], &cell_inv).map(wrap_frac);
            let mut df = 0.0;
            for c in 0..3 {
                let mut dd = frac[c] - orig[c];
                if dd > 0.5 {
                    dd -= 1.0;
                }
                if dd < -0.5 {
                    dd += 1.0;
                }
                df += dd * dd;
            }
            if df > 1e-6 {
                different = true;
            }
            let mut moved = atom.clone();
            moved.x = p[0];
            moved.y = p[1];
            moved.z = p[2];
            b.atoms.push(moved);
        }
        if different {
            return b;
        }
    }
    let mut b = Grain { cell: a.cell, atoms: Vec::new() };
    for atom in &a.atoms {
        let cart = [atom.x, -atom.y, -atom.z];
        let frac = cart_to_frac(&cart, &cell_inv).map(wrap_frac);
        let p = frac_to_cart(&frac, &a.cell);
        let mut moved = atom.clone();
        moved.x = p[0];
        moved.y = p[1];
        moved.z = p[2];
        b.atoms.push(moved);
    }
    b
}

/// Grain B on top of grain A along `dir`, with a gap and vacuum.
pub fn stack_grains(a: &Grain, b: &Grain, dir: usize, vacuum: f64, gap: f64) -> Grain {
    let pa = cell_to_parameters(&a.cell);
    let pb = cell_to_parameters(&b.cell);
    let l = if dir == 1 { (pa[dir] + gap) * (pb[5] * PI / 180.0).sin() } else { pa[dir] + gap };
    let mut abc = [pa[0], pa[1], pa[2]];
    abc[dir] += pb[dir] + 2.0 * gap + vacuum;
    let cell = parameters_to_cell(&[abc[0], abc[1], abc[2], pb[3], pb[4], pb[5]]);
    let inv = inv3(&cell);
    let mut result = Grain { cell, atoms: remap(&a.atoms, &inv, &cell) };
    let mut l_vec = [0.0; 3];
    l_vec[dir] = l;
    for atom in &b.atoms {
        let cart = [atom.x + l_vec[0], atom.y + l_vec[1], atom.z + l_vec[2]];
        let p = frac_to_cart(&cart_to_frac(&cart, &inv).map(wrap_frac), &cell);
        let mut moved = atom.clone();
        moved.x = p[0];
        moved.y = p[1];
        moved.z = p[2];
        result.atoms.push(moved);
    }
    result
}

/// Removes the later of every pair nearer than the threshold; the count removed.
pub fn remove_overlaps(atoms: &mut Vec<AtomSite>, threshold: f64) -> i32 {
    if threshold <= 0.0 || atoms.len() < 2 {
        return 0;
    }
    let dist_sq = threshold * threshold;
    let mut keep = vec![true; atoms.len()];
    let mut removed = 0;
    for i in 0..atoms.len() {
        if !keep[i] {
            continue;
        }
        for j in i + 1..atoms.len() {
            if !keep[j] {
                continue;
            }
            let (dx, dy, dz) = (atoms[i].x - atoms[j].x, atoms[i].y - atoms[j].y, atoms[i].z - atoms[j].z);
            if dx * dx + dy * dy + dz * dz < dist_sq {
                keep[j] = false;
                removed += 1;
            }
        }
    }
    if removed > 0 {
        let mut k = keep.iter();
        atoms.retain(|_| *k.next().unwrap());
    }
    removed
}

/// In-plane primitive reduction keeping the stacking vector (the boundary
/// stays at fractions 0 and 0.5 along it).
pub fn reduce_to_primitive_gb(s: &mut Structure, stack_dir: usize, tol: f64) -> bool {
    if !s.has_unit_cell || s.atoms.len() < 2 {
        return false;
    }
    let natom = s.atoms.len();
    let cell = s.cell_vectors;
    let cell_inv = inv3(&cell);
    let fc: Vec<[f64; 3]> = s.atoms.iter().map(|a| cart_to_frac(&[a.x, a.y, a.z], &cell_inv).map(wrap_frac)).collect();
    let tol_sq = tol * tol;
    let to_cart = |f: &[f64; 3], cell: &M3| -> [f64; 3] {
        let mut c = [0.0; 3];
        for d in 0..3 {
            for k in 0..3 {
                c[k] += f[d] * cell[d][k];
            }
        }
        c
    };
    let len2 = |c: &[f64; 3]| c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    let maps_all = |t: &[f64; 3]| -> bool {
        (0..natom).all(|i| {
            let sf = [wrap_frac(fc[i][0] + t[0]), wrap_frac(fc[i][1] + t[1]), wrap_frac(fc[i][2] + t[2])];
            (0..natom).any(|j| {
                if s.atoms[j].atomic_number != s.atoms[i].atomic_number {
                    return false;
                }
                let mut diff = [0.0; 3];
                for d in 0..3 {
                    diff[d] = sf[d] - fc[j][d];
                    diff[d] -= diff[d].round();
                }
                len2(&to_cart(&diff, &cell)) < tol_sq
            })
        })
    };
    // The most common species (the first of equal counts in Z order).
    let mut counts: BTreeMap<i32, i32> = BTreeMap::new();
    for a in &s.atoms {
        *counts.entry(a.atomic_number).or_default() += 1;
    }
    let mut ref_z = s.atoms[0].atomic_number;
    let mut max_c = 0;
    for (&z, &c) in &counts {
        if c > max_c {
            max_c = c;
            ref_z = z;
        }
    }
    let ref_idx = (0..natom).find(|&i| s.atoms[i].atomic_number == ref_z).unwrap_or(0);
    let mut valid: Vec<[f64; 3]> = Vec::new();
    for j in 0..natom {
        if j == ref_idx || s.atoms[j].atomic_number != ref_z {
            continue;
        }
        let mut t = [0.0; 3];
        for d in 0..3 {
            t[d] = fc[j][d] - fc[ref_idx][d];
            if t[d] > 0.5 {
                t[d] -= 1.0;
            }
            if t[d] < -0.5 {
                t[d] += 1.0;
            }
        }
        if t[stack_dir].abs() > 0.01 {
            continue;
        }
        t[stack_dir] = 0.0;
        if len2(&to_cart(&t, &cell)) < tol_sq {
            continue;
        }
        if maps_all(&t) {
            valid.push(t);
        }
    }
    if valid.is_empty() {
        return false;
    }
    let (d1, d2) = ((stack_dir + 1) % 3, (stack_dir + 2) % 3);
    let mut cand: Vec<[f64; 3]> = Vec::new();
    let mut add = |f: [f64; 3]| {
        let c = to_cart(&f, &cell);
        if len2(&c) > tol_sq {
            cand.push(c);
        }
    };
    let (mut e1, mut e2) = ([0.0; 3], [0.0; 3]);
    e1[d1] = 1.0;
    e2[d2] = 1.0;
    let neg = |v: [f64; 3]| v.map(|x| -x);
    let plus = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] { std::array::from_fn(|k| a[k] + b[k]) };
    let minus = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] { std::array::from_fn(|k| a[k] - b[k]) };
    add(e1);
    add(e2);
    add(neg(e1));
    add(neg(e2));
    for &t in &valid {
        add(t);
        add(neg(t));
        add(plus(e1, t));
        add(minus(e1, t));
        add(plus(e2, t));
        add(minus(e2, t));
    }
    for i in 0..valid.len() {
        for j in i + 1..valid.len() {
            add(plus(valid[i], valid[j]));
            add(minus(valid[i], valid[j]));
        }
    }
    let vlen = |a: &[f64; 3]| len2(a).sqrt();
    sort_by(&mut cand, |a, b| vlen(a) < vlen(b));
    let v1 = cand[0];
    let mut v2 = None;
    for c in &cand[1..] {
        let cr = [v1[1] * c[2] - v1[2] * c[1], v1[2] * c[0] - v1[0] * c[2], v1[0] * c[1] - v1[1] * c[0]];
        if len2(&cr).sqrt() > tol * vlen(&v1) {
            v2 = Some(*c);
            break;
        }
    }
    let Some(mut v2) = v2 else { return false };
    // Keep the handedness relative to the stacking direction.
    let cross = |a: &[f64; 3], b: &[f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let st = cell[stack_dir];
    let dot_stack = dot3(&cross(&v1, &v2), &st);
    let orig_dot = dot3(&cross(&cell[d1], &cell[d2]), &st);
    if dot_stack * orig_dot < 0.0 {
        v2 = neg(v2);
    }
    let mut new_cell = [[0.0; 3]; 3];
    new_cell[d1] = v1;
    new_cell[d2] = v2;
    new_cell[stack_dir] = cell[stack_dir];
    if det3(&new_cell).abs() >= det3(&cell).abs() - 1e-4 {
        return false;
    }
    let nc_inv = inv3(&new_cell);
    let mut new_atoms: Vec<AtomSite> = Vec::new();
    for atom in &s.atoms {
        let f = cart_to_frac(&[atom.x, atom.y, atom.z], &nc_inv).map(wrap_frac);
        let nc = frac_to_cart(&f, &new_cell);
        let duplicate = new_atoms.iter().any(|a| {
            let df = cart_to_frac(&[a.x, a.y, a.z], &nc_inv).map(wrap_frac);
            let mut dd = [0.0; 3];
            for d in 0..3 {
                dd[d] = f[d] - df[d];
                dd[d] -= dd[d].round();
            }
            len2(&to_cart(&dd, &new_cell)) < tol_sq
        });
        if !duplicate {
            let mut a = atom.clone();
            a.x = nc[0];
            a.y = nc[1];
            a.z = nc[2];
            new_atoms.push(a);
        }
    }
    s.cell_vectors = new_cell;
    s.atoms = new_atoms;
    true
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CslBicrystal {
    pub candidate: SigmaCandidate,
    pub plane: [i32; 3],
    pub removed: i32,
    pub structure: Structure,
}

/// The CLI/dialog pipeline: grain A from the CSL, grain B rotated, stacked
/// along the chosen plane, overlaps removed, optionally reduced in-plane,
/// with grain region ids (0 below the boundary at 0.5, 1 above).
#[allow(clippy::too_many_arguments)]
pub fn build_csl_bicrystal(
    input: &Structure,
    axis: &[i32; 3],
    sigma_target: i32,
    sigma_max: i32,
    plane_index: usize,
    uc_a: i32,
    uc_b: i32,
    vacuum: f32,
    gap: f32,
    overlap: f32,
    conventional: bool,
) -> Result<CslBicrystal, String> {
    let candidates = compute_gb_info(axis, sigma_max);
    if candidates.is_empty() {
        return Err(format!("no \u{3a3} candidates found for axis [{} {} {}] with max \u{3a3}={sigma_max}", axis[0], axis[1], axis[2]));
    }
    let selected = if sigma_target > 0 {
        candidates.iter().position(|c| c.sigma == sigma_target).ok_or_else(|| format!("\u{3a3}={sigma_target} not found for axis [{} {} {}]", axis[0], axis[1], axis[2]))?
    } else {
        0
    };
    if plane_index > 2 {
        return Err("--plane must be 0, 1, or 2".into());
    }
    Ok(build_bicrystal_from_candidate(input, &candidates[selected], plane_index, uc_a, uc_b, vacuum, gap, overlap, conventional))
}

/// The bicrystal of one Σ candidate (the dialog's selected entry): `plane_index`
/// (0..=2) picks the CSL column used as boundary plane and stacking direction.
#[allow(clippy::too_many_arguments)]
pub fn build_bicrystal_from_candidate(
    input: &Structure,
    sel: &SigmaCandidate,
    plane_index: usize,
    uc_a: i32,
    uc_b: i32,
    vacuum: f32,
    gap: f32,
    overlap: f32,
    conventional: bool,
) -> CslBicrystal {
    let sel = *sel;
    let direction = plane_index;
    let input_grain = structure_to_grain(input);
    let csl_t: [[i32; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| sel.csl[j][i]));
    let mut grain_a = make_supercell(&input_grain, &csl_t);
    if !is_orthogonal(&grain_a.cell, 0.01) {
        grain_a = set_orthogonal_grain(&grain_a, direction);
    }
    let mut scale_b = [1; 3];
    scale_b[direction] = uc_b;
    let temp_a = make_supercell_diag(&grain_a, scale_b[0], scale_b[1], scale_b[2]);
    let grain_b = make_supercell(&b_from_a(&temp_a), &[[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    let mut scale_a = [1; 3];
    scale_a[direction] = uc_a;
    grain_a = make_supercell_diag(&grain_a, scale_a[0], scale_a[1], scale_a[2]);
    let mut gb = stack_grains(&grain_a, &grain_b, direction, vacuum as f64, gap as f64);
    let removed = remove_overlaps(&mut gb.atoms, overlap as f64);
    let mut structure = grain_to_structure(&gb);
    // The periodic-boundary tolerance from the smallest layer spacing.
    let inv = inv3(&gb.cell);
    let mut layers: Vec<f64> = structure.atoms.iter().map(|a| wrap_frac(cart_to_frac(&[a.x, a.y, a.z], &inv)[direction])).collect();
    sort_by(&mut layers, |a, b| a < b);
    let mut min_spacing = 1.0;
    for w in layers.windows(2) {
        let d = w[1] - w[0];
        if d > 1e-8 && d < min_spacing {
            min_spacing = d;
        }
    }
    structure.pbc_boundary_tol = (min_spacing * 0.5 + 1e-6) as f32;
    if !conventional {
        reduce_to_primitive_gb(&mut structure, direction, 0.25);
    }
    let inv = inv3(&structure.cell_vectors);
    structure.grain_region_ids =
        structure.atoms.iter().map(|a| if wrap_frac(cart_to_frac(&[a.x, a.y, a.z], &inv)[direction]) >= 0.5 { 1 } else { 0 }).collect();
    CslBicrystal { candidate: sel, plane: sel.plane[plane_index], removed, structure }
}

/// The primitive cell from spglib's standardization (the C++ reduceToPrimitive):
/// atoms keep the symbol and colour of their element.
pub fn reduce_to_primitive(s: &mut Structure, symprec: f64) -> bool {
    if !s.has_unit_cell || s.atoms.is_empty() {
        return false;
    }
    let natom = s.atoms.len();
    let cell_inv = inv3(&s.cell_vectors);
    let positions: Vec<[f64; 3]> = s
        .atoms
        .iter()
        .map(|a| {
            let mut f = [0.0; 3];
            for j in 0..3 {
                for i in 0..3 {
                    f[j] += [a.x, a.y, a.z][i] * cell_inv[i][j];
                }
            }
            f.map(|x| {
                let w = x - x.floor();
                if w >= 1.0 - 1e-8 {
                    0.0
                } else {
                    w
                }
            })
        })
        .collect();
    let types: Vec<i32> = s.atoms.iter().map(|a| a.atomic_number).collect();
    let Some(prim) = af_science::spg::standardize_primitive(&s.cell_vectors, &positions, &types, symprec) else { return false };
    let n = prim.positions.len();
    if n == 0 || n > natom {
        return false;
    }
    s.cell_vectors = prim.lattice;
    let mut info: BTreeMap<i32, (String, f32, f32, f32)> = BTreeMap::new();
    for a in &s.atoms {
        info.insert(a.atomic_number, (a.symbol.clone(), a.r, a.g, a.b));
    }
    s.atoms = prim
        .positions
        .iter()
        .zip(&prim.types)
        .map(|(f, &z)| {
            let mut cart = [0.0; 3];
            for j in 0..3 {
                for i in 0..3 {
                    cart[j] += f[i] * prim.lattice[i][j];
                }
            }
            let (symbol, r, g, b) = info.get(&z).cloned().unwrap_or((af_core::elements::symbol(z).to_string(), 0.5, 0.5, 0.5));
            AtomSite { x: cart[0], y: cart[1], z: cart[2], atomic_number: z, symbol, r, g, b, ..Default::default() }
        })
        .collect();
    true
}
