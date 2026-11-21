//! Three-dimensional coordinates for a molecular graph by distance geometry
//! (the C++ chemistry/Embedding): bounds from UFF bond lengths and
//! hybridisation angles, cis/trans limits for 1-4 pairs, van der Waals lower
//! limits, triangle smoothing, random distances embedded through the metric
//! matrix in four dimensions, refinement with chirality and planarity terms,
//! then a UFF minimisation. Stereocentres and stereo double bonds are enforced.
//!
//! Parity: the metric-matrix eigenvectors come from a different eigensolver
//! than Eigen's, so embeddings are not bit-identical to the C++ ones; the
//! tests check the stereochemistry, the geometry and the energy instead
//! (docs/rust-port-plan.md section 6a).

use std::collections::{BTreeSet, HashMap};

use af_core::cxx_random::{uniform_real, Mt19937};
use nalgebra::DMatrix;

use crate::forcefield::{minimize, uff_bond_length, uff_type, ForceFieldOptions, UffForceField};
use crate::molecule::*;
use crate::smiles::{parse_smiles, SmilesError};

const PI: f64 = 3.14159265358979323846;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EmbedOptions {
    /// Embeddings tried (the lowest UFF energy is kept).
    pub attempts: i32,
    pub seed: u32,
    /// UFF minimisation of each embedding.
    pub optimize: bool,
    pub optimize_steps: i32,
    pub force_tolerance: f64,
}

impl Default for EmbedOptions {
    fn default() -> Self {
        Self { attempts: 5, seed: 1, optimize: true, optimize_steps: 3000, force_tolerance: 0.01 }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EmbedResult {
    /// UFF energy (kcal/mol) of the result.
    pub energy: f64,
    /// Embeddings that kept the stereochemistry.
    pub succeeded: i32,
    /// Stereocentres without a descriptor (given a random configuration).
    pub unspecified_stereo: i32,
}

fn vdw_lower(z: i32) -> f64 {
    match z {
        1 => 1.10,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        15 | 16 => 1.80,
        17 => 1.75,
        _ => 1.80,
    }
}

/// Limited-memory BFGS on a flat vector; f returns the value and fills the gradient.
fn minimize_flat(x: &mut Vec<f64>, mut f: impl FnMut(&[f64], &mut [f64]) -> f64, max_iter: i32, gtol: f64) {
    let dim = x.len();
    let mut g = vec![0.0; dim];
    let mut e = f(x, &mut g);
    let (mut s_h, mut y_h, mut rho): (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>) = (Vec::new(), Vec::new(), Vec::new());
    for _ in 0..max_iter {
        let mut gmax = 0.0f64;
        for &v in &g {
            gmax = fmax(gmax, v.abs());
        }
        if gmax < gtol {
            break;
        }
        let mut q = g.clone();
        let mut alpha = vec![0.0; s_h.len()];
        for k in (0..s_h.len()).rev() {
            let mut a = 0.0;
            for c in 0..dim {
                a += s_h[k][c] * q[c];
            }
            a *= rho[k];
            alpha[k] = a;
            for c in 0..dim {
                q[c] -= a * y_h[k][c];
            }
        }
        let mut gamma = 1e-3;
        if let (Some(s), Some(y)) = (s_h.last(), y_h.last()) {
            let (mut sy, mut yy) = (0.0, 0.0);
            for c in 0..dim {
                sy += s[c] * y[c];
                yy += y[c] * y[c];
            }
            if yy > 0.0 {
                gamma = sy / yy;
            }
        }
        for v in q.iter_mut() {
            *v *= gamma;
        }
        for k in 0..s_h.len() {
            let mut b = 0.0;
            for c in 0..dim {
                b += y_h[k][c] * q[c];
            }
            b *= rho[k];
            for c in 0..dim {
                q[c] += s_h[k][c] * (alpha[k] - b);
            }
        }
        let mut slope = 0.0;
        for c in 0..dim {
            q[c] = -q[c];
            slope += q[c] * g[c];
        }
        if slope >= 0.0 {
            s_h.clear();
            y_h.clear();
            rho.clear();
            slope = 0.0;
            for c in 0..dim {
                q[c] = -g[c] * 1e-3;
                slope += q[c] * g[c];
            }
        }
        let mut longest = 0.0f64;
        for &v in &q {
            longest = fmax(longest, v.abs());
        }
        let mut t = if longest > 0.5 { 0.5 / longest } else { 1.0 };
        let mut trial = vec![0.0; dim];
        let mut gt = vec![0.0; dim];
        let mut et = e;
        let mut ok = false;
        for _ in 0..30 {
            for c in 0..dim {
                trial[c] = x[c] + t * q[c];
            }
            et = f(&trial, &mut gt);
            if et <= e + 1e-4 * t * slope {
                ok = true;
                break;
            }
            t *= 0.5;
        }
        if !ok {
            if s_h.is_empty() {
                break;
            }
            s_h.clear();
            y_h.clear();
            rho.clear();
            continue;
        }
        let mut s = vec![0.0; dim];
        let mut y = vec![0.0; dim];
        let mut sy = 0.0;
        for c in 0..dim {
            s[c] = trial[c] - x[c];
            y[c] = gt[c] - g[c];
            sy += s[c] * y[c];
        }
        let previous = e;
        x.copy_from_slice(&trial);
        g.copy_from_slice(&gt);
        e = et;
        if sy > 1e-12 {
            s_h.push(s);
            y_h.push(y);
            rho.push(1.0 / sy);
            if s_h.len() > 8 {
                s_h.remove(0);
                y_h.remove(0);
                rho.remove(0);
            }
        }
        if (previous - e).abs() < 1e-12 * fmax(1.0, e.abs()) {
            break;
        }
    }
}

struct Bounds {
    n: usize,
    lower: Vec<f64>,
    upper: Vec<f64>,
}

impl Bounds {
    fn index(&self, i: i32, j: i32) -> usize {
        (i.min(j) as usize) * self.n + i.max(j) as usize
    }
    fn lo(&mut self, i: i32, j: i32) -> &mut f64 {
        let k = self.index(i, j);
        &mut self.lower[k]
    }
    fn up(&mut self, i: i32, j: i32) -> &mut f64 {
        let k = self.index(i, j);
        &mut self.upper[k]
    }
    fn get(&self, i: i32, j: i32) -> (f64, f64) {
        let k = self.index(i, j);
        (self.lower[k], self.upper[k])
    }
}

fn four_body_distance(rij: f64, rjk: f64, rkl: f64, t1: f64, t2: f64, phi: f64) -> f64 {
    let (ix, iy) = (rij * t1.cos(), rij * t1.sin());
    let (lx, ly, lz) = (rjk - rkl * t2.cos(), rkl * t2.sin() * phi.cos(), rkl * t2.sin() * phi.sin());
    ((ix - lx) * (ix - lx) + (iy - ly) * (iy - ly) + lz * lz).sqrt()
}

fn tetrahedral_volume(x: &[f64], dim: usize, n: &[i32; 4], g: &mut [[f64; 3]; 4]) -> f64 {
    let p = |a: i32, c: usize| x[a as usize * dim + c];
    let mut u = [0.0; 3];
    let mut v = [0.0; 3];
    let mut w = [0.0; 3];
    for c in 0..3 {
        u[c] = p(n[1], c) - p(n[0], c);
        v[c] = p(n[2], c) - p(n[0], c);
        w[c] = p(n[3], c) - p(n[0], c);
    }
    let vxw = [v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2], v[0] * w[1] - v[1] * w[0]];
    let wxu = [w[1] * u[2] - w[2] * u[1], w[2] * u[0] - w[0] * u[2], w[0] * u[1] - w[1] * u[0]];
    let uxv = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    for c in 0..3 {
        g[1][c] = vxw[c];
        g[2][c] = wxu[c];
        g[3][c] = uxv[c];
        g[0][c] = -(vxw[c] + wxu[c] + uxv[c]);
    }
    u[0] * vxw[0] + u[1] * vxw[1] + u[2] * vxw[2]
}

/// Checks the stereo descriptors against the coordinates.
pub fn stereo_matches(m: &Molecule) -> bool {
    for t in &m.tetrahedral {
        if t.n.iter().any(|&v| v < 0) {
            continue;
        }
        let p = |i: i32| m.atoms[i as usize].position;
        let v = signed_volume(p(t.centre), p(t.n[1]), p(t.n[2]), p(t.n[3]));
        if v * f64::from(t.sign) <= 0.0 {
            return false;
        }
    }
    for d in &m.double_bonds {
        let (pa, px, py, pb) = (m.atoms[d.a as usize].position, m.atoms[d.x as usize].position, m.atoms[d.y as usize].position, m.atoms[d.b as usize].position);
        let mut axis = [py[0] - px[0], py[1] - px[1], py[2] - px[2]];
        let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        if len < 1e-9 {
            return false;
        }
        for v in axis.iter_mut() {
            *v /= len;
        }
        let (mut u, mut w, mut du, mut dw) = ([0.0; 3], [0.0; 3], 0.0, 0.0);
        for c in 0..3 {
            u[c] = pa[c] - px[c];
            w[c] = pb[c] - py[c];
            du += u[c] * axis[c];
            dw += w[c] * axis[c];
        }
        let mut dot_uw = 0.0;
        for c in 0..3 {
            dot_uw += (u[c] - du * axis[c]) * (w[c] - dw * axis[c]);
        }
        if (dot_uw > 0.0) != d.cis {
            return false;
        }
    }
    true
}

/// Coordinates for a molecule whose hydrogens are explicit atoms. Fails when
/// no attempt reproduces the requested stereochemistry.
pub fn embed(m: &mut Molecule, options: &EmbedOptions) -> Result<EmbedResult, String> {
    let n = m.atoms.len();
    if m.atoms.iter().any(|a| a.hydrogens > 0) {
        return Err("embedding needs explicit hydrogens".into());
    }
    if n == 0 {
        return Err("no atoms to embed".into());
    }
    if n > 1500 {
        return Err("distance geometry is limited to 1500 atoms; build large molecules from fragments".into());
    }
    let mut result = EmbedResult::default();
    if n == 1 {
        m.atoms[0].position = [0.0; 3];
        return Ok(result);
    }
    let ni = n as i32;
    perceive_aromaticity(m);
    let rings = smallest_rings(m);
    let types: Vec<String> = (0..ni).map(|i| uff_type(m, i)).collect::<Result<_, _>>()?;
    // Ideal bond lengths.
    let mut length = Vec::with_capacity(m.bonds.len());
    for b in &m.bonds {
        let mut order = if b.aromatic { 1.5 } else { f64::from(b.order) };
        if b.dative {
            order = 1.0;
        }
        length.push(uff_bond_length(&types[b.a as usize], &types[b.b as usize], order)?);
    }
    // Ring angles: smallest ring through consecutive i-j-k.
    let mut ring_of: HashMap<(i32, i32, i32), usize> = HashMap::new();
    for r in &rings {
        let s = r.len();
        for q in 0..s {
            let (i, j, k) = (r[(q + s - 1) % s], r[q], r[(q + 1) % s]);
            let key = (i.min(k), j, i.max(k));
            match ring_of.get(&key) {
                Some(&v) if s >= v => {}
                _ => {
                    ring_of.insert(key, s);
                }
            }
        }
    }
    let ideal_angle = |i: i32, j: i32, k: i32| -> f64 {
        if is_metal(m.atoms[j as usize].z) {
            return -1.0;
        }
        let hyb = hybridisation(m, j);
        if let Some(&s) = ring_of.get(&(i.min(k), j, i.max(k))) {
            match s {
                3 => return 60.0,
                4 => return 89.0,
                5 => return if hyb == 2 { 108.0 } else { 104.5 },
                6 if hyb == 2 => return 120.0,
                _ => {}
            }
        }
        if hyb == 1 {
            return 180.0;
        }
        if hyb == 2 {
            return 120.0;
        }
        if m.atoms[j as usize].z == 8 {
            return 104.5;
        }
        if m.atoms[j as usize].z == 16 && m.degree(j) == 2 {
            return 97.0;
        }
        if m.degree(j) > 4 {
            return -1.0;
        }
        109.47
    };
    let mut bounds = Bounds { n, lower: vec![0.0; n * n], upper: vec![1000.0; n * n] };
    for i in 0..ni {
        for j in i + 1..ni {
            *bounds.lo(i, j) = 0.8 * (vdw_lower(m.atoms[i as usize].z) + vdw_lower(m.atoms[j as usize].z));
        }
    }
    let mut bond_to: Vec<HashMap<i32, f64>> = vec![HashMap::new(); n];
    for (k, b) in m.bonds.iter().enumerate() {
        bond_to[b.a as usize].insert(b.b, length[k]);
        bond_to[b.b as usize].insert(b.a, length[k]);
    }
    let bt = |a: i32, b: i32| bond_to[a as usize].get(&b).copied().unwrap_or(0.0);
    // 1-4 pairs: between the cis and trans distances (fixed for stereo double
    // bonds and ring bonds of aromatic rings).
    for (kb, b) in m.bonds.iter().enumerate() {
        let (j, k) = (b.a, b.b);
        for i in m.neighbours(j) {
            if i == k {
                continue;
            }
            let t1 = ideal_angle(i, j, k);
            if t1 < 0.0 {
                continue;
            }
            for l in m.neighbours(k) {
                if l == j || l == i {
                    continue;
                }
                let t2 = ideal_angle(j, k, l);
                if t2 < 0.0 {
                    continue;
                }
                if m.bond_between(i, l) >= 0 {
                    continue; // 3-ring
                }
                let cis = four_body_distance(bt(i, j), length[kb], bt(k, l), t1 * PI / 180.0, t2 * PI / 180.0, 0.0);
                let trans = four_body_distance(bt(i, j), length[kb], bt(k, l), t1 * PI / 180.0, t2 * PI / 180.0, PI);
                let (mut lo, mut hi) = (fmin(cis, trans), fmax(cis, trans));
                let planar = (b.order == 2 && !b.dative) || b.aromatic;
                if planar {
                    let mut fixed = false;
                    let mut want_cis = false;
                    // The neighbours on the same end of a stereo double bond.
                    if b.order == 2 && !b.aromatic {
                        for d in &m.double_bonds {
                            if !((d.x == j && d.y == k) || (d.x == k && d.y == j)) {
                                continue;
                            }
                            let (a, c) = if d.x == j { (d.a, d.b) } else { (d.b, d.a) };
                            let (same_a, same_c) = (i == a, l == c);
                            want_cis = d.cis == (same_a == same_c);
                            fixed = true;
                        }
                    }
                    // Aromatic or ring double bonds: i and l both in or both out of a common ring -> cis.
                    if !fixed && bond_in_ring(m, kb as i32) {
                        let (mut i_in, mut l_in) = (false, false);
                        for r in &rings {
                            if !(r.contains(&j) && r.contains(&k)) {
                                continue;
                            }
                            for &v in r {
                                i_in |= v == i;
                                l_in |= v == l;
                            }
                        }
                        want_cis = i_in == l_in;
                        fixed = true;
                    }
                    if fixed {
                        lo = if want_cis { cis } else { trans };
                        hi = lo;
                    }
                }
                *bounds.lo(i, l) = lo - 0.05;
                let u = bounds.up(i, l);
                *u = fmin(*u, hi + 0.05);
            }
        }
    }
    // 1-3 pairs from the ideal angles.
    for j in 0..ni {
        let nb = m.neighbours(j);
        for x in 0..nb.len() {
            for y in x + 1..nb.len() {
                let (i, k) = (nb[x], nb[y]);
                let t = ideal_angle(i, j, k);
                let (a, c) = (bt(j, i), bt(j, k));
                if t < 0.0 {
                    let l = bounds.lo(i, k);
                    *l = fmin(*l, (a * a + c * c).sqrt());
                    continue;
                }
                let d = (a * a + c * c - 2.0 * a * c * (t * PI / 180.0).cos()).sqrt();
                *bounds.lo(i, k) = d - 0.03;
                *bounds.up(i, k) = d + 0.03;
            }
        }
    }
    // 1-2 pairs.
    for (k, b) in m.bonds.iter().enumerate() {
        *bounds.lo(b.a, b.b) = length[k] - 0.01;
        *bounds.up(b.a, b.b) = length[k] + 0.01;
    }
    // Triangle smoothing.
    if n <= 600 {
        for k in 0..ni {
            for i in 0..ni {
                if i == k {
                    continue;
                }
                let (lik, uik) = bounds.get(i, k);
                for j in i + 1..ni {
                    if j == k {
                        continue;
                    }
                    let (lkj, ukj) = bounds.get(k, j);
                    let idx = bounds.index(i, j);
                    let mut u = bounds.upper[idx];
                    let mut l = bounds.lower[idx];
                    if uik + ukj < u {
                        u = uik + ukj;
                    }
                    let candidate = fmax(lik - ukj, lkj - uik);
                    if candidate > l {
                        l = candidate;
                    }
                    if l > u {
                        l = u;
                    }
                    bounds.upper[idx] = u;
                    bounds.lower[idx] = l;
                }
            }
        }
    }
    let cap = 4.0 + 3.0 * (n as f64).cbrt() * 1.6;
    for u in bounds.upper.iter_mut() {
        *u = fmin(*u, fmax(cap, 3.0));
    }

    // Stereo and planarity targets.
    let chiral: Vec<([i32; 4], i32)> = m.tetrahedral.iter().filter(|t| t.n.iter().all(|&v| v >= 0)).map(|t| (t.n, t.sign)).collect();
    let mut planar: Vec<[i32; 4]> = Vec::new();
    for j in 0..ni {
        if m.degree(j) == 3 && hybridisation(m, j) == 2 {
            let nb = m.neighbours(j);
            planar.push([j, nb[0], nb[1], nb[2]]);
        }
    }
    // Unspecified stereocentres (reported).
    {
        let rank = canonical_ranks(m);
        for i in 0..ni {
            if m.degree(i) != 4 || hybridisation(m, i) != 3 || is_metal(m.atoms[i as usize].z) {
                continue;
            }
            let distinct: BTreeSet<i32> = m.neighbours(i).iter().map(|&x| rank[x as usize]).collect();
            if distinct.len() != 4 {
                continue;
            }
            if !m.tetrahedral.iter().any(|t| t.centre == i) {
                result.unspecified_stereo += 1;
            }
        }
    }

    let mut random = Mt19937::new(options.seed);
    let mut best_energy = 1e300;
    let mut best: Vec<Vec3> = Vec::new();
    let mut failures = 0;
    let dim = 4usize;
    let mut attempt = 0;
    while attempt < options.attempts.max(1) + failures && failures < 40 {
        attempt += 1;
        // Random distances within the bounds and the metric matrix.
        let mut d2 = DMatrix::<f64>::zeros(n, n);
        for i in 0..n {
            for j in i + 1..n {
                let (l, u0) = bounds.get(i as i32, j as i32);
                let u = fmax(u0, l);
                let d = l + uniform_real(&mut random, 0.0, 1.0) * (u - l);
                d2[(i, j)] = d * d;
                d2[(j, i)] = d * d;
            }
        }
        let total = d2.sum() / (2.0 * (n * n) as f64);
        let d0: Vec<f64> = (0..n).map(|i| d2.row(i).sum() / n as f64 - total).collect();
        let metric = DMatrix::<f64>::from_fn(n, n, |i, j| 0.5 * (d0[i] + d0[j] - d2[(i, j)]));
        let solver = metric.symmetric_eigen();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| solver.eigenvalues[a].total_cmp(&solver.eigenvalues[b]));
        let mut x = vec![0.0; n * dim];
        for c in 0..dim.min(n) {
            let col = order[n - 1 - c];
            let lambda = solver.eigenvalues[col];
            let s = if lambda > 0.0 { lambda.sqrt() } else { 0.0 };
            for i in 0..n {
                x[i * dim + c] = s * solver.eigenvectors[(i, col)];
            }
        }
        // Small noise breaks exact degeneracies (planar or linear starts).
        for v in x.iter_mut() {
            *v += uniform_real(&mut random, -0.05, 0.05);
        }
        let error = |p: &[f64], g: &mut [f64], active: usize, fourth: f64| -> f64 {
            g.fill(0.0);
            let mut e = 0.0;
            for i in 0..n {
                for j in i + 1..n {
                    let mut diff = [0.0; 4];
                    let mut r2 = 0.0;
                    for c in 0..active {
                        diff[c] = p[i * dim + c] - p[j * dim + c];
                        r2 += diff[c] * diff[c];
                    }
                    let (l, u) = bounds.get(i as i32, j as i32);
                    let coef;
                    if r2 > u * u {
                        let t = r2 / (u * u) - 1.0;
                        e += t * t;
                        coef = 2.0 * t / (u * u);
                    } else if r2 < l * l {
                        let t = 2.0 * l * l / (l * l + r2) - 1.0;
                        e += t * t;
                        coef = 2.0 * t * (-2.0 * l * l / ((l * l + r2) * (l * l + r2)));
                    } else {
                        continue;
                    }
                    for c in 0..active {
                        g[i * dim + c] += coef * 2.0 * diff[c];
                        g[j * dim + c] -= coef * 2.0 * diff[c];
                    }
                }
            }
            let mut gv = [[0.0; 3]; 4];
            for (atoms, sign) in &chiral {
                let v = tetrahedral_volume(p, dim, atoms, &mut gv);
                let target = 2.0;
                let s = f64::from(*sign) * v;
                if s < target {
                    let t = target - s;
                    e += t * t;
                    for q in 0..4 {
                        for c in 0..3 {
                            g[atoms[q] as usize * dim + c] += 2.0 * t * f64::from(-*sign) * gv[q][c];
                        }
                    }
                }
            }
            for pl in &planar {
                let v = tetrahedral_volume(p, dim, pl, &mut gv);
                let w = 0.1;
                e += w * v * v;
                for q in 0..4 {
                    for c in 0..3 {
                        g[pl[q] as usize * dim + c] += w * 2.0 * v * gv[q][c];
                    }
                }
            }
            if fourth > 0.0 {
                for i in 0..n {
                    let w4 = p[i * dim + 3];
                    e += fourth * w4 * w4;
                    g[i * dim + 3] += 2.0 * fourth * w4;
                }
            }
            if active < 4 {
                for i in 0..n {
                    g[i * dim + 3] = 0.0;
                }
            }
            e
        };
        minimize_flat(&mut x, |p, g| error(p, g, 4, 0.0), 400, 1e-4);
        minimize_flat(&mut x, |p, g| error(p, g, 4, 1.0), 400, 1e-4);
        for i in 0..n {
            x[i * dim + 3] = 0.0;
        }
        minimize_flat(&mut x, |p, g| error(p, g, 3, 0.0), 600, 1e-5);
        for i in 0..n {
            m.atoms[i].position = [x[i * dim], x[i * dim + 1], x[i * dim + 2]];
        }
        if !stereo_matches(m) {
            failures += 1;
            continue;
        }
        let mut positions: Vec<Vec3> = m.atoms.iter().map(|a| a.position).collect();
        let energy = (|| -> Result<f64, String> {
            let field = UffForceField::new(m, None, ForceFieldOptions::default())?;
            if options.optimize {
                Ok(minimize(&field, &mut positions, &[], options.optimize_steps, options.force_tolerance)?.energy)
            } else {
                field.energy(&positions, None, None)
            }
        })()
        .unwrap_or(0.0); // elements without parameters: keep the geometric embedding
        let saved: Vec<Vec3> = m.atoms.iter().map(|a| a.position).collect();
        for i in 0..n {
            m.atoms[i].position = positions[i];
        }
        if !stereo_matches(m) {
            for i in 0..n {
                m.atoms[i].position = saved[i];
            }
            failures += 1;
            continue;
        }
        result.succeeded += 1;
        if energy < best_energy {
            best_energy = energy;
            best = positions;
        }
        if result.succeeded >= options.attempts.max(1) {
            break;
        }
    }
    if best.is_empty() {
        return Err("could not embed the molecule with the requested stereochemistry".into());
    }
    for i in 0..n {
        m.atoms[i].position = best[i];
    }
    // Centre on the centroid.
    let mut centre = [0.0; 3];
    for a in &m.atoms {
        for c in 0..3 {
            centre[c] += a.position[c] / n as f64;
        }
    }
    for a in m.atoms.iter_mut() {
        for c in 0..3 {
            a.position[c] -= centre[c];
        }
    }
    result.energy = best_energy;
    Ok(result)
}

/// Parses a SMILES, adds explicit hydrogens and embeds it.
pub fn molecule_from_smiles(smiles: &str, options: &EmbedOptions) -> Result<(Molecule, EmbedResult), SmilesError> {
    let mut m = parse_smiles(smiles)?;
    make_hydrogens_explicit(&mut m, false);
    // Dummy atoms ('*') become hydrogens for the geometry.
    for a in m.atoms.iter_mut() {
        if a.dummy {
            a.z = 1;
        }
    }
    let r = embed(&mut m, options).map_err(|message| SmilesError { message, position: None })?;
    Ok((m, r))
}
