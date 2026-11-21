//! Field analysis on grids (the C++ electronic/FieldAnalysis): finite
//! differences, gradient, Laplacian, energy densities, real-space bonding
//! descriptors (RDG, sign(lambda_2) rho, DORI), smoothing, profiles,
//! sections and contours.

use af_core::glm::{inverse_d, mul_d, normalize_d};
use glam::{DMat3, DVec2, DVec3, IVec2, IVec3};

use crate::volume::{Grid, HARTREE, PI};

// Constants GCC folds with MPFR (correctly rounded), as in the C++:
/// std::pow(bohr, 3).
const BOHR3: f64 = f64::from_bits(0x3fc2f7b774c52326);
/// std::pow(bohr, 5).
const BOHR5: f64 = f64::from_bits(0x3fa53ef4e1d9bb05);
/// std::pow(bohr, 8).
const BOHR8: f64 = f64::from_bits(0x3f792fc2fd929001);
/// std::pow(3 pi^2, 2.0 / 3).
const PI2_TWO_THIRDS: f64 = f64::from_bits(0x4023243d46b7cf6d);
/// std::pow(3 pi^2, 1.0 / 3).
const PI2_THIRD: f64 = f64::from_bits(0x4008bfd4dd6882ca);

fn derivative(g: &Grid, axis: usize, second: bool) -> Result<Grid, String> {
    let mut out = g.clone();
    let n = g.shape[axis];
    if !g.periodic && n < if second { 4 } else { 3 } {
        return Err("Finite differences require >=3 samples per axis (>=4 for Laplacian)".into());
    }
    for z in 0..g.shape[2] {
        for y in 0..g.shape[1] {
            for x in 0..g.shape[0] {
                let p = IVec3::new(x, y, z);
                let i = p[axis];
                let at = |j: i32| {
                    let mut q = p;
                    q[axis] = if g.periodic { (j % n + n) % n } else { j };
                    g.values[g.index(q.x, q.y, q.z)]
                };
                let value = if g.periodic || (i > 0 && i < n - 1) {
                    if second {
                        at(i + 1) - 2.0 * at(i) + at(i - 1)
                    } else {
                        (at(i + 1) - at(i - 1)) * 0.5
                    }
                } else if i == 0 {
                    if second {
                        2.0 * at(0) - 5.0 * at(1) + 4.0 * at(2) - at(3)
                    } else {
                        (-3.0 * at(0) + 4.0 * at(1) - at(2)) * 0.5
                    }
                } else if second {
                    2.0 * at(i) - 5.0 * at(i - 1) + 4.0 * at(i - 2) - at(i - 3)
                } else {
                    (3.0 * at(i) - 4.0 * at(i - 1) + at(i - 2)) * 0.5
                };
                out.values[g.index(x, y, z)] = value;
            }
        }
    }
    Ok(out)
}

/// Cartesian gradient components (1/Angstrom), central differences.
pub fn gradient(grid: &Grid) -> Result<[Grid; 3], String> {
    grid.validate()?;
    let mut result = [grid.clone(), grid.clone(), grid.clone()];
    let d = [derivative(grid, 0, false)?, derivative(grid, 1, false)?, derivative(grid, 2, false)?];
    let transform = inverse_d(&grid.steps()).transpose();
    for i in 0..grid.values.len() {
        let v = mul_d(&transform, DVec3::new(d[0].values[i], d[1].values[i], d[2].values[i]));
        for axis in 0..3 {
            result[axis].values[i] = v[axis];
        }
    }
    for (axis, r) in result.iter_mut().enumerate() {
        r.name = format!("gradient {axis}");
        r.unit = format!("({})/A", grid.unit);
        r.validate()?;
    }
    Ok(result)
}

pub fn laplacian(grid: &Grid) -> Result<Grid, String> {
    grid.validate()?;
    let mut result = grid.clone();
    result.values.iter_mut().for_each(|v| *v = 0.0);
    let inverse = inverse_d(&grid.steps());
    let metric: DMat3 = inverse * inverse.transpose();
    for a in 0..3 {
        for b in a..3 {
            let m = metric.col(a)[b];
            if m.abs() < 1e-16 {
                continue;
            }
            let d = if a == b { derivative(grid, a, true)? } else { derivative(&derivative(grid, a, false)?, b, false)? };
            for i in 0..result.values.len() {
                result.values[i] += (if a == b { 1.0 } else { 2.0 }) * m * d.values[i];
            }
        }
    }
    result.name = "Laplacian".into();
    result.unit = format!("({})/A^2", grid.unit);
    result.validate()?;
    Ok(result)
}

/// Kinetic (gradient expansion), potential (local virial) and total energy densities (eV/A^3).
pub fn energy_density(density: &Grid, floor: f64) -> Result<[Grid; 3], String> {
    if density.unit != "e/A^3" {
        return Err("Energy conversion requires electron density in e/A^3".into());
    }
    if !floor.is_finite() || floor <= 0.0 {
        return Err("Density floor must be positive and finite".into());
    }
    let grad = gradient(density)?;
    let lap = laplacian(density)?;
    let mut result = [density.clone(), density.clone(), density.clone()];
    for i in 0..density.values.len() {
        if density.values[i] < 0.0 {
            return Err("Energy conversion requires nonnegative density, not spin/difference density".into());
        }
        let rho = density.values[i] * BOHR3;
        if density.values[i] <= floor {
            for f in result.iter_mut() {
                f.values[i] = 0.0;
            }
            continue;
        }
        let mut norm = 0.0;
        for f in &grad {
            norm += f.values[i] * f.values[i] * BOHR8;
        }
        let l = lap.values[i] * BOHR5;
        let g = 0.3 * PI2_TWO_THIRDS * rho.powf(5.0 / 3.0) + norm / (72.0 * rho) + l / 12.0;
        let v = l / 4.0 - 2.0 * g;
        let factor = HARTREE / BOHR3;
        result[0].values[i] = g * factor;
        result[1].values[i] = v * factor;
        result[2].values[i] = (g + v) * factor;
    }
    let names = ["kinetic energy density (gradient expansion)", "potential energy density (local virial)", "total energy density (local virial)"];
    for (r, name) in result.iter_mut().zip(names) {
        r.name = name.into();
        r.unit = "eV/A^3".into();
        r.validate()?;
    }
    Ok(result)
}

fn validate_density(density: &Grid, floor: f64) -> Result<(), String> {
    if density.unit != "e/A^3" {
        return Err("This descriptor requires electron density in e/A^3".into());
    }
    if !floor.is_finite() || floor <= 0.0 {
        return Err("Density floor must be positive and finite".into());
    }
    if density.values.iter().any(|&v| v < 0.0) {
        return Err("This descriptor requires nonnegative density, not spin/difference density".into());
    }
    Ok(())
}

/// Middle eigenvalue (lambda_2) of a real symmetric 3x3 matrix via the closed-form
/// trigonometric solution (Smith, Comm. ACM 4, 168, 1961).
fn middle_eigenvalue(a11: f64, a22: f64, a33: f64, a12: f64, a13: f64, a23: f64) -> f64 {
    let p1 = a12 * a12 + a13 * a13 + a23 * a23;
    if p1 < 1e-300 {
        let mut lo = a11;
        let mut hi = a11;
        for v in [a22, a33] {
            if v < lo {
                lo = v;
            }
            if hi < v {
                hi = v;
            }
        }
        return a11 + a22 + a33 - lo - hi;
    }
    let q = (a11 + a22 + a33) / 3.0;
    let p = (((a11 - q) * (a11 - q) + (a22 - q) * (a22 - q) + (a33 - q) * (a33 - q) + 2.0 * p1) / 6.0).sqrt();
    let (b11, b22, b33, b12, b13, b23) = ((a11 - q) / p, (a22 - q) / p, (a33 - q) / p, a12 / p, a13 / p, a23 / p);
    let mut r = (b11 * (b22 * b33 - b23 * b23) - b12 * (b12 * b33 - b23 * b13) + b13 * (b12 * b23 - b22 * b13)) / 2.0;
    r = r.clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let eig1 = q + 2.0 * p * phi.cos();
    let eig3 = q + 2.0 * p * (phi + 2.0 * PI / 3.0).cos();
    3.0 * q - eig1 - eig3
}

/// d(rho)/dx,dy,dz and the Hessian xx, yy, zz, xy, xz, yz (mixed partials symmetrized).
fn density_hessian(density: &Grid) -> Result<([Grid; 3], [Grid; 6]), String> {
    let grad = gradient(density)?;
    let gx = gradient(&grad[0])?;
    let gy = gradient(&grad[1])?;
    let gz = gradient(&grad[2])?;
    let (mut hxy, mut hxz, mut hyz) = (gx[1].clone(), gx[2].clone(), gy[2].clone());
    for i in 0..hxy.values.len() {
        hxy.values[i] = (gx[1].values[i] + gy[0].values[i]) / 2.0;
        hxz.values[i] = (gx[2].values[i] + gz[0].values[i]) / 2.0;
        hyz.values[i] = (gy[2].values[i] + gz[1].values[i]) / 2.0;
    }
    let [gx0, _, _] = gx;
    let [_, gy1, _] = gy;
    let [_, _, gz2] = gz;
    Ok((grad, [gx0, gy1, gz2, hxy, hxz, hyz]))
}

pub fn reduced_density_gradient(density: &Grid, floor: f64) -> Result<Grid, String> {
    validate_density(density, floor)?;
    let grad = gradient(density)?;
    let mut result = density.clone();
    let constant = 2.0 * PI2_THIRD;
    for i in 0..density.values.len() {
        let rho = density.values[i];
        if rho <= floor {
            result.values[i] = 0.0;
            continue;
        }
        let grad_norm = (grad[0].values[i] * grad[0].values[i] + grad[1].values[i] * grad[1].values[i] + grad[2].values[i] * grad[2].values[i]).sqrt();
        result.values[i] = grad_norm / (constant * rho.powf(4.0 / 3.0));
    }
    result.name = "reduced density gradient".into();
    result.unit = "dimensionless".into();
    result.validate()?;
    Ok(result)
}

pub fn signed_density(density: &Grid, floor: f64) -> Result<Grid, String> {
    validate_density(density, floor)?;
    let (_, h) = density_hessian(density)?;
    let mut result = density.clone();
    for i in 0..density.values.len() {
        let rho = density.values[i];
        if rho <= floor {
            result.values[i] = 0.0;
            continue;
        }
        let lambda2 = middle_eigenvalue(h[0].values[i], h[1].values[i], h[2].values[i], h[3].values[i], h[4].values[i], h[5].values[i]);
        result.values[i] = if lambda2 == 0.0 { 0.0 } else { rho.copysign(lambda2) };
    }
    result.name = "sign(lambda_2) * density".into();
    result.unit = "e/A^3".into();
    result.validate()?;
    Ok(result)
}

/// Density overlap regions indicator (de Silva and Corminboeuf).
pub fn dori(density: &Grid, floor: f64) -> Result<Grid, String> {
    validate_density(density, floor)?;
    let (grad, h6) = density_hessian(density)?;
    let mut result = density.clone();
    for i in 0..density.values.len() {
        let rho = density.values[i];
        if rho <= floor {
            result.values[i] = 0.0;
            continue;
        }
        let g = [grad[0].values[i], grad[1].values[i], grad[2].values[i]];
        let h = [
            [h6[0].values[i], h6[3].values[i], h6[4].values[i]],
            [h6[3].values[i], h6[1].values[i], h6[5].values[i]],
            [h6[4].values[i], h6[5].values[i], h6[2].values[i]],
        ];
        let k = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) / (rho * rho);
        let mut grad_gamma_squared = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                let jacobian = (h[a][b] * rho - g[a] * g[b]) / (rho * rho);
                grad_gamma_squared += jacobian * jacobian;
            }
        }
        let theta = grad_gamma_squared / (k * k * k);
        result.values[i] = theta / (1.0 + theta);
    }
    result.name = "DORI".into();
    result.unit = "dimensionless".into();
    result.validate()?;
    Ok(result)
}

/// Gaussian smoothing (sigma in Angstrom, kernel radius in grid steps).
pub fn smooth(grid: &Grid, sigma: f64, radius: i32) -> Result<Grid, String> {
    grid.validate()?;
    if !sigma.is_finite() || sigma <= 0.0 || !(1..=50).contains(&radius) {
        return Err("Smoothing needs sigma > 0 Angstrom and radius 1..50 grid steps".into());
    }
    let mut out = grid.clone();
    let steps = grid.steps();
    let mut kernel: Vec<(IVec3, f64)> = Vec::new();
    for z in -radius..=radius {
        for y in -radius..=radius {
            for x in -radius..=radius {
                let delta = mul_d(&steps, DVec3::new(x as f64, y as f64, z as f64));
                kernel.push((IVec3::new(x, y, z), (-delta.dot(delta) / (2.0 * sigma * sigma)).exp()));
            }
        }
    }
    if kernel.len() as f64 * grid.values.len() as f64 > 2e9 {
        return Err("Smoothing workload too large; reduce radius or grid resolution".into());
    }
    for z in 0..grid.shape[2] {
        for y in 0..grid.shape[1] {
            for x in 0..grid.shape[0] {
                let (mut sum, mut weights) = (0.0, 0.0);
                for (offset, w) in &kernel {
                    let mut p = IVec3::new(x, y, z) + *offset;
                    let mut valid = true;
                    for a in 0..3 {
                        if grid.periodic {
                            p[a] = (p[a] % grid.shape[a] + grid.shape[a]) % grid.shape[a];
                        } else if p[a] < 0 || p[a] >= grid.shape[a] {
                            valid = false;
                        }
                    }
                    if valid {
                        sum += w * grid.values[grid.index(p.x, p.y, p.z)];
                        weights += w;
                    }
                }
                out.values[grid.index(x, y, z)] = sum / weights;
            }
        }
    }
    out.name = format!("Gaussian smoothed {}", grid.name);
    Ok(out)
}

/// (distance, value) samples along a segment.
pub fn line_profile(grid: &Grid, start: DVec3, end: DVec3, count: i32) -> Result<Vec<DVec2>, String> {
    grid.validate()?;
    if !(2..=1_000_000).contains(&count) {
        return Err("Line profile count must be 2..1000000".into());
    }
    let length = (end - start).length();
    let mut out = Vec::with_capacity(count as usize);
    for i in 0..count {
        let t = i as f64 / (count - 1) as f64;
        out.push(DVec2::new(length * t, grid.sample(start + t * (end - start))?));
    }
    Ok(out)
}

/// (height, plane average) along a cell axis.
pub fn planar_average(grid: &Grid, axis: i32) -> Result<Vec<DVec2>, String> {
    grid.validate()?;
    if !(0..=2).contains(&axis) {
        return Err("Axis must be 0, 1 or 2".into());
    }
    let axis = axis as usize;
    let n = grid.shape[axis] as usize;
    let mut out = vec![DVec2::ZERO; n];
    let mut weights = vec![0.0; n];
    let inverse = inverse_d(&grid.cell);
    let height = 1.0 / DVec3::new(inverse.col(0)[axis], inverse.col(1)[axis], inverse.col(2)[axis]).length();
    for z in 0..grid.shape[2] {
        for y in 0..grid.shape[1] {
            for x in 0..grid.shape[0] {
                let i = [x, y, z][axis] as usize;
                let w = grid.weight(x, y, z);
                out[i].y += grid.values[grid.index(x, y, z)] * w;
                weights[i] += w;
            }
        }
    }
    let divisions = (grid.shape[axis] - if grid.periodic { 0 } else { 1 }) as f64;
    for i in 0..n {
        out[i].x = height * i as f64 / divisions;
        out[i].y /= weights[i];
    }
    Ok(out)
}

/// The planar average smoothed by a window of `window` (odd) planes.
pub fn macroscopic_average(grid: &Grid, axis: i32, window: i32) -> Result<Vec<DVec2>, String> {
    let mut out = planar_average(grid, axis)?;
    let raw = out.clone();
    let n = out.len() as i32;
    if window < 1 || window % 2 == 0 || window > n {
        return Err("Macroscopic window must be odd and no longer than the profile".into());
    }
    for i in 0..n {
        let (mut sum, mut count) = (0.0, 0);
        for d in -window / 2..=window / 2 {
            let mut j = i + d;
            if grid.periodic {
                j = (j + n) % n;
            }
            if j >= 0 && j < n {
                sum += raw[j as usize].y;
                count += 1;
            }
        }
        out[i as usize].y = sum / count as f64;
    }
    Ok(out)
}

/// A planar section spanned by u and v (nu x nv samples), as a two-plane grid.
pub fn section(grid: &Grid, origin: DVec3, u: DVec3, v: DVec3, nu: i32, nv: i32) -> Result<Grid, String> {
    grid.validate()?;
    if nu < 2 || nv < 2 || nu as f64 * nv as f64 > 10_000_000.0 || u.cross(v).length() < 1e-12 {
        return Err("Section requires independent spanning vectors and valid dimensions".into());
    }
    let mut out = Grid {
        shape: [nu, nv, 2],
        origin,
        cell: DMat3::from_cols(u, v, normalize_d(u.cross(v))),
        name: "section".into(),
        unit: grid.unit.clone(),
        values: vec![0.0; nu as usize * nv as usize * 2],
        ..Default::default()
    };
    // Both planes hold the same samples; z = 0 is the actual section surface.
    for y in 0..nv {
        for x in 0..nu {
            let value = grid.sample(origin + u * (x as f64 / (nu - 1) as f64) + v * (y as f64 / (nv - 1) as f64))?;
            let (i0, i1) = (out.index(x, y, 0), out.index(x, y, 1));
            out.values[i0] = value;
            out.values[i1] = value;
        }
    }
    out.validate()?;
    Ok(out)
}

/// Contour segments (point pairs) of a section at `level`.
pub fn contours(slice: &Grid, level: f64) -> Result<Vec<DVec3>, String> {
    slice.validate()?;
    if !level.is_finite() {
        return Err("Contour level must be finite".into());
    }
    let mut out = Vec::new();
    // Consistent triangulation removes the ambiguous marching-square saddle.
    let triangles = [[0usize, 1, 2], [0, 2, 3]];
    for y in 0..slice.shape[1] - 1 {
        for x in 0..slice.shape[0] - 1 {
            let p = [IVec2::new(x, y), IVec2::new(x + 1, y), IVec2::new(x + 1, y + 1), IVec2::new(x, y + 1)];
            for tri in &triangles {
                let mut hits = Vec::new();
                for e in 0..3 {
                    let (a, b) = (p[tri[e]], p[tri[(e + 1) % 3]]);
                    let (va, vb) = (slice.values[slice.index(a.x, a.y, 0)], slice.values[slice.index(b.x, b.y, 0)]);
                    if (va < level) != (vb < level) {
                        let t = (level - va) / (vb - va);
                        // glm::mix(x, y, a) = x * (1 - a) + y * a.
                        hits.push(slice.position(a.x, a.y, 0) * (1.0 - t) + slice.position(b.x, b.y, 0) * t);
                    }
                }
                if hits.len() == 2 {
                    out.extend(hits);
                }
            }
        }
    }
    Ok(out)
}
