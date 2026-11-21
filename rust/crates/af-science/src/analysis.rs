//! Trajectory, local-environment, electronic and fitting analyses (the C++
//! science/Analysis): MSD, diffusion, VACF, vibrational spectra, local
//! strain, Ackland-Jones and diamond structure types, centrosymmetry,
//! Steinhardt bond order, Wigner-Seitz occupancy, structure factors, band
//! gaps, effective masses, work functions, Birch-Murnaghan EOS, elastic
//! tensors, phonon DOS and harmonic thermodynamics.
//!
//! Eigen's SVD least squares, symmetric eigensolver and LDLT go through
//! nalgebra (equal to about 1e-12 relative, docs/rust-port-plan.md 6a);
//! everything else follows the C++ operation order.

use af_core::cxx::{sort_by, stable_sort_by};
use af_core::json::Json;
use nalgebra::{DMatrix, DVector, Matrix3, Matrix4, Vector4};
use num_complex::Complex64;

use crate::core::*;
use crate::data::Parameters;

const PI: f64 = std::f64::consts::PI;
const BOLTZMANN: f64 = 8.617333262145e-5;
const EV_PER_A3_TO_GPA: f64 = 160.2176634;

pub(crate) type R = Result<Json, String>;

// ---------------------------------------------------------------- JSON helpers
pub(crate) fn jv(values: &[f64]) -> Json {
    Json::Array(values.iter().map(|&v| Json::Number(v)).collect())
}
pub(crate) fn ji(values: &[i32]) -> Json {
    Json::Array(values.iter().map(|&v| Json::Number(v as f64)).collect())
}
pub(crate) fn jb(values: &[bool]) -> Json {
    Json::Array(values.iter().map(|&v| Json::Bool(v)).collect())
}
fn matrix_json(m: &DMatrix<f64>) -> Json {
    Json::Array((0..m.nrows()).map(|i| Json::Array((0..m.ncols()).map(|j| Json::Number(m[(i, j)])).collect())).collect())
}
fn vector_json(v: &[f64]) -> Json {
    jv(v)
}

fn fmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}
fn fmin(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

// ---------------------------------------------------------------- helpers
fn trajectory(values: NdArray, name: &str) -> Result<NdArray, String> {
    let result = finite_array(values, name, 3)?;
    if result.shape[0] < 2 || result.shape[2] != 3 {
        return Err(format!("{name} must have shape (at least two frames, atoms, 3)"));
    }
    Ok(result)
}

pub(crate) fn vector1d(values: NdArray, name: &str) -> Result<Vec<f64>, String> {
    Ok(finite_array(values, name, 1)?.values)
}

fn optional_integer(p: &Parameters, name: &str, minimum: i64) -> Result<Option<i64>, String> {
    if !p.has(name) {
        return Ok(None);
    }
    Ok(Some(integer(p.number(name)?, name, minimum)?))
}

fn lag_times(limit: i64, dt: f64) -> Vec<f64> {
    (0..=limit as usize).map(|i| i as f64 * dt).collect()
}

fn origins(count: usize, limit: i64) -> Vec<i32> {
    (0..=limit as usize).map(|i| (count - i) as i32).collect()
}

fn trapezoid(y: &[f64], x: &[f64]) -> f64 {
    let mut area = 0.0;
    for i in 1..y.len() {
        area += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }
    area
}

/// NumPy-compatible numerical rank (max(M,N) eps sigma_max tolerance).
fn matrix_rank(m: &DMatrix<f64>) -> usize {
    if m.is_empty() {
        return 0;
    }
    let sv = m.clone().svd(false, false).singular_values;
    let top = sv.iter().cloned().fold(0.0f64, f64::max);
    let tolerance = top * m.nrows().max(m.ncols()) as f64 * f64::EPSILON;
    sv.iter().filter(|&&s| s > tolerance).count()
}

/// Minimum-norm least squares (Eigen's complete orthogonal decomposition).
fn least_squares(a: &DMatrix<f64>, b: &DMatrix<f64>) -> DMatrix<f64> {
    let svd = a.clone().svd(true, true);
    let top = svd.singular_values.iter().cloned().fold(0.0f64, f64::max);
    let tolerance = top * a.nrows().max(a.ncols()) as f64 * f64::EPSILON;
    svd.solve(b, tolerance).unwrap_or_else(|_| DMatrix::zeros(a.ncols(), b.ncols()))
}

/// Symmetric eigenpairs, eigenvalues ascending (Eigen's order).
fn symmetric_eigen(m: &DMatrix<f64>) -> (Vec<f64>, DMatrix<f64>) {
    let e = m.clone().symmetric_eigen();
    let mut order: Vec<usize> = (0..e.eigenvalues.len()).collect();
    order.sort_by(|&a, &b| e.eigenvalues[a].partial_cmp(&e.eigenvalues[b]).unwrap_or(std::cmp::Ordering::Equal));
    let values = order.iter().map(|&k| e.eigenvalues[k]).collect();
    let vectors = DMatrix::from_fn(m.nrows(), m.ncols(), |i, j| e.eigenvectors[(i, order[j])]);
    (values, vectors)
}

struct Linear {
    slope: f64,
    intercept: f64,
    r: f64,
    standard_error: f64,
}

fn linear_regression(x: &[f64], y: &[f64]) -> Linear {
    let n = x.len() as f64;
    let xm = x.iter().fold(0.0, |s, v| s + v) / n;
    let ym = y.iter().fold(0.0, |s, v| s + v) / n;
    let (mut ssxm, mut ssym, mut ssxym) = (0.0, 0.0, 0.0);
    for i in 0..x.len() {
        ssxm += (x[i] - xm) * (x[i] - xm);
        ssym += (y[i] - ym) * (y[i] - ym);
        ssxym += (x[i] - xm) * (y[i] - ym);
    }
    ssxm /= n;
    ssym /= n;
    ssxym /= n;
    let mut r = if ssxm == 0.0 || ssym == 0.0 { 0.0 } else { ssxym / (ssxm * ssym).sqrt() };
    r = r.clamp(-1.0, 1.0);
    let slope = ssxym / ssxm;
    let df = n - 2.0;
    let standard_error = fmax(0.0, (1.0 - r * r) * ssym / ssxm / df).sqrt();
    Linear { slope, intercept: ym - slope * xm, r, standard_error }
}

pub(crate) fn neighbors_by_center(neighbors: &[Neighbor], atoms: usize) -> Vec<Vec<usize>> {
    let mut result = vec![Vec::new(); atoms];
    for (n, nb) in neighbors.iter().enumerate() {
        result[nb.i].push(n);
    }
    result
}

/// The fully normalised associated Legendre function (LocalGeometry).
pub(crate) fn normalized_legendre(l: i32, m: i32, x: f64) -> f64 {
    let mut pmm = 1.0;
    if m > 0 {
        let omx2 = (1.0 - x) * (1.0 + x);
        let mut fact = 1.0;
        for _ in 1..=m {
            pmm *= omx2 * fact / (fact + 1.0);
            fact += 2.0;
        }
    }
    pmm = ((2 * m + 1) as f64 * pmm / (4.0 * PI)).sqrt();
    if m & 1 != 0 {
        pmm = -pmm;
    }
    if l == m {
        return pmm;
    }
    let mut pmmp1 = x * (2.0 * m as f64 + 3.0).sqrt() * pmm;
    if l == m + 1 {
        return pmmp1;
    }
    let mut old_factor = (2.0 * m as f64 + 3.0).sqrt();
    let mut pll = 0.0;
    for ll in m + 2..=l {
        let llf = ll as f64;
        let factor = ((4.0 * llf * llf - 1.0) / (llf * llf - m as f64 * m as f64)).sqrt();
        pll = (x * pmmp1 - pmm / old_factor) * factor;
        old_factor = factor;
        pmm = pmmp1;
        pmmp1 = pll;
    }
    pll
}

fn local_neighbors(positions: &[Vec3], cutoff: f64, p: &Parameters, cell_name: &str, pbc: Pbc) -> Result<(Vec<Neighbor>, Mat3), String> {
    let cell = if p.has(cell_name) { cell_matrix(p.array(cell_name)?)? } else { identity() };
    if (pbc[0] || pbc[1] || pbc[2]) && !p.has(cell_name) {
        return Err("Periodic analysis requires a cell".into());
    }
    let neighbors = neighbor_list(positions, &cell, pbc, cutoff)?;
    if neighbors.iter().any(|n| norm(n.vector) < 1e-10) {
        return Err("Coincident atoms make local environment analysis undefined".into());
    }
    Ok((neighbors, cell))
}

fn remove_frame_drift(a: &mut NdArray) {
    let (count, atoms) = (a.shape[0], a.shape[1]);
    for f in 0..count {
        for k in 0..3 {
            let mut mean = 0.0;
            for i in 0..atoms {
                mean += a.values[(f * atoms + i) * 3 + k];
            }
            mean /= atoms as f64;
            for i in 0..atoms {
                a.values[(f * atoms + i) * 3 + k] -= mean;
            }
        }
    }
}

// ---------------------------------------------------------------- dynamics
pub fn mean_square_displacement(p: &Parameters) -> R {
    let mut c = trajectory(p.array("positions")?, "positions")?;
    let wrapped = p.boolean("wrapped", false)?;
    let remove_drift = p.boolean("remove_drift", false)?;
    let dt = positive(p.number("timestep_fs")?, "timestep_fs")?;
    let (count, atoms) = (c.shape[0], c.shape[1]);
    let at = |f: usize, a: usize, k: usize| (f * atoms + a) * 3 + k;
    if wrapped {
        if !p.has("cell") {
            return Err("cell must be a nonempty finite array of the required dimension".into());
        }
        let mic = MinimumImage::new(&cell_matrix(p.array("cell")?)?, p.pbc("pbc", [true; 3])?)?;
        let raw = c.values.clone();
        for f in 1..count {
            for a in 0..atoms {
                let step = mic.apply([
                    raw[at(f, a, 0)] - raw[at(f - 1, a, 0)],
                    raw[at(f, a, 1)] - raw[at(f - 1, a, 1)],
                    raw[at(f, a, 2)] - raw[at(f - 1, a, 2)],
                ]);
                for k in 0..3 {
                    c.values[at(f, a, k)] = c.values[at(f - 1, a, k)] + step[k];
                }
            }
        }
    }
    if remove_drift {
        remove_frame_drift(&mut c);
    }
    let limit = optional_integer(p, "max_lag", 1)?.unwrap_or(count as i64 - 1);
    if limit >= count as i64 {
        return Err("max_lag must be smaller than frame count".into());
    }
    let mut covariance = NdArray::new(vec![limit as usize + 1, 3, 3], 0.0);
    for lag in 1..=limit as usize {
        let mut sum = [[0.0f64; 3]; 3];
        let samples = count - lag;
        for t in 0..samples {
            for a in 0..atoms {
                let d: [f64; 3] = std::array::from_fn(|k| c.values[at(t + lag, a, k)] - c.values[at(t, a, k)]);
                for i in 0..3 {
                    for j in 0..3 {
                        sum[i][j] += d[i] * d[j];
                    }
                }
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                covariance.values[(lag * 3 + i) * 3 + j] = sum[i][j] / (samples * atoms) as f64;
            }
        }
    }
    let n = limit as usize + 1;
    let mut msd = vec![0.0; n];
    let mut components = NdArray::new(vec![n, 3], 0.0);
    for lag in 0..n {
        for k in 0..3 {
            components.values[lag * 3 + k] = covariance.values[(lag * 3 + k) * 3 + k];
        }
        msd[lag] = components.values[lag * 3] + components.values[lag * 3 + 1] + components.values[lag * 3 + 2];
    }
    let mut result = Json::object();
    result.set("lag_fs", jv(&lag_times(limit, dt)));
    result.set("msd_A2", jv(&msd));
    result.set("components_A2", array_json(&components));
    result.set("displacement_tensor_A2", array_json(&covariance));
    result.set("origins", ji(&origins(count, limit)));
    Ok(result)
}

pub fn diffusion_coefficient(p: &Parameters) -> R {
    let times = vector1d(p.array("lag_fs")?, "lag_fs")?;
    let values = vector1d(p.array("msd_A2")?, "msd_A2")?;
    let dimension_value = p.number_or("dimensions", 3.0)?;
    if dimension_value != 1.0 && dimension_value != 2.0 && dimension_value != 3.0 {
        return Err("dimensions must be 1, 2 or 3".into());
    }
    let dimensions = dimension_value as i32;
    let mut valid = times.len() == values.len() && times[0] >= 0.0;
    let mut i = 0;
    while valid && i < times.len() {
        valid = values[i] >= 0.0 && (i == 0 || times[i] > times[i - 1]);
        i += 1;
    }
    if !valid {
        return Err("MSD must be nonnegative with matching, increasing nonnegative times".into());
    }
    let interval = vector1d(p.array("fit_range_fs")?, "fit_range_fs")?;
    if interval.len() != 2 || interval[0] >= interval[1] || interval[0] < times[0] || interval[1] > *times.last().unwrap() {
        return Err("fit_range_fs must lie inside the sampled time interval".into());
    }
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for i in 0..times.len() {
        if times[i] >= interval[0] && times[i] <= interval[1] {
            x.push(times[i]);
            y.push(values[i]);
        }
    }
    if x.len() < 3 {
        return Err("At least three fit points are required".into());
    }
    let fit = linear_regression(&x, &y);
    if fit.slope < 0.0 {
        return Err("Negative MSD slope is not a physical diffusion estimate".into());
    }
    let diffusion = fit.slope / (2 * dimensions) as f64;
    let mut result = Json::object();
    result.set("D_A2_per_fs", diffusion);
    result.set("D_m2_per_s", diffusion * 1e-5);
    result.set("slope_A2_per_fs", fit.slope);
    result.set("intercept_A2", fit.intercept);
    result.set("r_squared", if fit.r.is_finite() { fit.r * fit.r } else { 0.0 });
    result.set("ols_slope_stderr", fit.standard_error);
    result.set("fit_points", x.len() as f64);
    result.set("fit_range_fs", jv(&interval));
    result.set("lag_fs", jv(&times));
    result.set("msd_A2", jv(&values));
    Ok(result)
}

pub fn velocity_autocorrelation(p: &Parameters) -> R {
    let mut v = trajectory(p.array("velocities")?, "velocities")?;
    let remove_drift = p.boolean("remove_drift", false)?;
    let normalize = p.boolean("normalize", false)?;
    let dt = positive(p.number("timestep_fs")?, "timestep_fs")?;
    let (count, atoms) = (v.shape[0], v.shape[1]);
    if remove_drift {
        remove_frame_drift(&mut v);
    }
    let limit = optional_integer(p, "max_lag", 1)?.unwrap_or(count as i64 - 1);
    if limit >= count as i64 {
        return Err("max_lag must be smaller than frame count".into());
    }
    let mut size = 1usize;
    while size < 2 * count - 1 {
        size <<= 1;
    }
    let mut correlation = vec![0.0; limit as usize + 1];
    for a in 0..atoms {
        for k in [0usize, 2] {
            let mut series = vec![Complex64::new(0.0, 0.0); size];
            for f in 0..count {
                series[f] = Complex64::new(v.values[(f * atoms + a) * 3 + k], if k + 1 < 3 { v.values[(f * atoms + a) * 3 + k + 1] } else { 0.0 });
            }
            let spectrum = fft(series, false);
            let mut power = vec![Complex64::new(0.0, 0.0); size];
            for q in 0..size {
                let other = spectrum[(size - q) % size].conj();
                let s = spectrum[q] + other;
                let x = Complex64::new(0.5 * s.re, 0.5 * s.im);
                let y = Complex64::new(0.0, -0.5) * (spectrum[q] - other);
                power[q] = Complex64::new(x.norm_sqr() + y.norm_sqr(), 0.0);
            }
            let back = fft(power, true);
            for lag in 0..correlation.len() {
                correlation[lag] += back[lag].re;
            }
        }
    }
    for lag in 0..correlation.len() {
        correlation[lag] /= atoms as f64 * (count - lag) as f64;
    }
    if normalize {
        if correlation[0] <= 0.0 {
            return Err("Cannot normalize zero velocity autocorrelation".into());
        }
        let zero = correlation[0];
        for value in correlation.iter_mut() {
            *value /= zero;
        }
    }
    let mut result = Json::object();
    result.set("lag_fs", jv(&lag_times(limit, dt)));
    result.set("vacf", jv(&correlation));
    result.set("unit", if normalize { "dimensionless" } else { "A^2/fs^2" });
    result.set("origins", ji(&origins(count, limit)));
    Ok(result)
}

pub fn vibrational_spectrum(p: &Parameters) -> R {
    let remove_mean = p.boolean("remove_mean", true)?;
    let mut v = trajectory(p.array("velocities")?, "velocities")?;
    let dt = positive(p.number("timestep_fs")?, "timestep_fs")?;
    let (count, atoms) = (v.shape[0], v.shape[1]);
    if count < 4 {
        return Err("At least four velocity frames are required".into());
    }
    let at = |f: usize, a: usize, k: usize| (f * atoms + a) * 3 + k;
    if remove_mean {
        for a in 0..atoms {
            for k in 0..3 {
                let mut mean = 0.0;
                for f in 0..count {
                    mean += v.values[at(f, a, k)];
                }
                mean /= count as f64;
                for f in 0..count {
                    v.values[at(f, a, k)] -= mean;
                }
            }
        }
    }
    if p.has("masses") {
        let weights = vector1d(p.array("masses")?, "masses")?;
        if weights.len() != atoms || !weights.iter().all(|&w| w > 0.0) {
            return Err("One positive mass is required per atom".into());
        }
        for f in 0..count {
            for a in 0..atoms {
                for k in 0..3 {
                    v.values[at(f, a, k)] *= weights[a].sqrt();
                }
            }
        }
    }
    let window = if p.has("window") { p.json("window")?.string()?.to_string() } else { "hann".to_string() };
    if window != "hann" && window != "none" {
        return Err("window must be hann or none".into());
    }
    let mut taper = vec![1.0; count];
    if window == "hann" {
        for f in 0..count {
            taper[f] = 0.5 - 0.5 * (2.0 * PI * f as f64 / (count - 1) as f64).cos();
        }
    }
    let bins = count / 2 + 1;
    let mut power = vec![0.0; bins];
    for a in 0..atoms {
        for k in 0..3 {
            let series: Vec<Complex64> = (0..count).map(|f| Complex64::new(v.values[at(f, a, k)] * taper[f], 0.0)).collect();
            let spectrum = fft(series, false);
            for q in 0..bins {
                power[q] += spectrum[q].norm_sqr();
            }
        }
    }
    let end = if count % 2 == 0 { bins - 1 } else { bins };
    for q in 1..end {
        power[q] *= 2.0;
    }
    let frequency: Vec<f64> = (0..bins).map(|q| q as f64 / (count as f64 * dt) * 1000.0).collect();
    let area = trapezoid(&power, &frequency);
    if !(area > 0.0) {
        return Err("No nonzero vibrational spectral weight".into());
    }
    for value in power.iter_mut() {
        *value /= area;
    }
    let mut result = Json::object();
    result.set("frequency_THz", jv(&frequency));
    result.set("density_per_THz", jv(&power));
    result.set("resolution_THz", 1000.0 / (count as f64 * dt));
    result.set("nyquist_THz", 500.0 / dt);
    Ok(result)
}

// ---------------------------------------------------------------- local environments
pub fn local_strain(p: &Parameters) -> R {
    let reference = points(p.array("reference")?, "positions")?;
    let cutoff = positive(p.number("cutoff_A")?, "cutoff_A")?;
    let pbc = p.pbc("pbc", [false; 3])?;
    let periodic = pbc.iter().any(|&b| b);
    let (neighbors, reference_cell) = local_neighbors(&reference, cutoff, p, "reference_cell", pbc)?;
    let now = points(p.array("current")?, "current")?;
    if now.len() != reference.len() {
        return Err("Reference/current atom ordering and counts must correspond".into());
    }
    let cell = if p.has("current_cell") { cell_matrix(p.array("current_cell")?)? } else { reference_cell };
    let mic = MinimumImage::new(&cell, pbc)?;
    let mapping = multiply(&inverse(&reference_cell)?, &cell);
    let count = reference.len();
    let mut deformation = NdArray::new(vec![count, 3, 3], f64::NAN);
    let mut strain = NdArray::new(vec![count, 3, 3], f64::NAN);
    let mut residual = vec![f64::NAN; count];
    let mut valid = vec![false; count];
    let mut coordination = vec![0i32; count];
    let groups = neighbors_by_center(&neighbors, count);
    for index in 0..count {
        let group = &groups[index];
        coordination[index] = group.len() as i32;
        if group.len() < 3 {
            continue;
        }
        let rows = group.len();
        let mut original = DMatrix::zeros(rows, 3);
        let mut displaced = DMatrix::zeros(rows, 3);
        for (n, &g) in group.iter().enumerate() {
            let nb = &neighbors[g];
            let mut d = add(sub(now[nb.j], now[index]), row_times([nb.shift[0] as f64, nb.shift[1] as f64, nb.shift[2] as f64], &cell));
            if periodic {
                let prediction = row_times(nb.vector, &mapping);
                d = add(prediction, mic.apply(sub(d, prediction)));
            }
            for k in 0..3 {
                original[(n, k)] = nb.vector[k];
                displaced[(n, k)] = d[k];
            }
        }
        if matrix_rank(&original) < 3 {
            continue;
        }
        let transposed = least_squares(&original, &displaced);
        let green = (&transposed * transposed.transpose() - DMatrix::<f64>::identity(3, 3)) / 2.0;
        for i in 0..3 {
            for j in 0..3 {
                deformation.values[(index * 3 + i) * 3 + j] = transposed[(j, i)];
                strain.values[(index * 3 + i) * 3 + j] = green[(i, j)];
            }
        }
        residual[index] = (&displaced - &original * &transposed).norm_squared();
        valid[index] = true;
    }
    let mut result = Json::object();
    result.set("deformation_gradient", array_json(&deformation));
    result.set("green_lagrange_strain", array_json(&strain));
    result.set("d2min_A2", jv(&residual));
    result.set("valid", jb(&valid));
    result.set("coordination", ji(&coordination));
    Ok(result)
}

/// Ackland and Jones (2006): 0 other, 1 fcc, 2 hcp, 3 bcc, 4 icosahedral.
pub(crate) fn ackland_jones_type(mut bonds: Vec<Vec3>) -> i32 {
    if bonds.len() < 6 {
        return 0;
    }
    sort_by(&mut bonds, |a, b| dot(*a, *a) < dot(*b, *b));
    let mut r0 = 0.0;
    for b in bonds.iter().take(6) {
        r0 += dot(*b, *b);
    }
    r0 /= 6.0;
    let mut shell = Vec::new();
    let mut n1 = 0;
    for b in &bonds {
        let r2 = dot(*b, *b);
        if r2 < 1.45 * r0 {
            shell.push(*b);
        }
        if r2 < 1.55 * r0 {
            n1 += 1;
        }
    }
    let mut chi = [0i32; 8];
    for j in 0..shell.len() {
        for k in j + 1..shell.len() {
            let c = dot(shell[j], shell[k]) / (norm(shell[j]) * norm(shell[k]));
            let bin = if c < -0.945 {
                0
            } else if c < -0.915 {
                1
            } else if c < -0.755 {
                2
            } else if c < -0.195 {
                3
            } else if c < 0.195 {
                4
            } else if c < 0.245 {
                5
            } else if c < 0.795 {
                6
            } else {
                7
            };
            chi[bin] += 1;
        }
    }
    let cp = (1.0 - chi[6] as f64 / 24.0).abs();
    let bcc = if chi[5] + chi[6] - chi[4] != 0 { 0.35 * chi[4] as f64 / (chi[5] + chi[6] - chi[4]) as f64 } else { f64::INFINITY };
    let fcc = 0.61 * (((chi[0] + chi[1]) as f64 - 6.0).abs() + chi[2] as f64) / 6.0;
    let hcp = ((chi[0] as f64 - 3.0).abs() + ((chi[0] + chi[1] + chi[2] + chi[3]) as f64 - 9.0).abs()) / 12.0;
    let (mut delta_bcc, mut delta_fcc, mut delta_hcp) = (bcc, fcc, hcp);
    if chi[0] == 7 {
        delta_bcc = 0.0;
    } else if chi[0] == 6 {
        delta_fcc = 0.0;
    } else if chi[0] <= 3 {
        delta_hcp = 0.0;
    }
    if chi[7] > 0 {
        return 0;
    }
    if chi[4] < 3 {
        return if !(11..=13).contains(&n1) { 0 } else { 4 };
    }
    if delta_bcc <= cp {
        return if n1 < 11 { 0 } else { 3 };
    }
    if !(11..=12).contains(&n1) {
        return 0;
    }
    if delta_fcc < delta_hcp {
        1
    } else {
        2
    }
}

fn automatic_cutoff(p: &Parameters, positions: usize, factor: f64) -> Result<f64, String> {
    let mut cutoff = p.number_or("cutoff_A", 0.0)?;
    if cutoff <= 0.0 {
        if !p.has("cell") {
            return Err("Set cutoff_A for structures without a cell".into());
        }
        let volume = cell_volume(&cell_matrix(p.array("cell")?)?) / positions as f64;
        cutoff = factor * (3.0 * volume / (4.0 * PI)).cbrt();
    }
    Ok(cutoff)
}

fn classification(result: &mut Json, key: &str, types: &[i32], names: &[&str]) {
    let mut counts = vec![0i32; names.len()];
    for &t in types {
        counts[t as usize] += 1;
    }
    result.set(key, ji(types));
    result.set("type_names", Json::Array(names.iter().map(|n| Json::String(n.to_string())).collect()));
    let (mut count_json, mut fractions) = (Json::object(), Json::object());
    for (t, name) in names.iter().enumerate() {
        count_json.set(name, counts[t] as f64);
        fractions.set(name, if types.is_empty() { 0.0 } else { counts[t] as f64 / types.len() as f64 });
    }
    result.set("counts", count_json);
    result.set("fractions", fractions);
}

pub fn structure_type(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    let pbc = p.pbc("pbc", [false; 3])?;
    let cutoff = automatic_cutoff(p, positions.len(), 3.0)?;
    let (neighbors, _) = local_neighbors(&positions, cutoff, p, "cell", pbc)?;
    let groups = neighbors_by_center(&neighbors, positions.len());
    let types: Vec<i32> = (0..positions.len()).map(|i| ackland_jones_type(groups[i].iter().map(|&g| neighbors[g].vector).collect())).collect();
    let mut result = Json::object();
    classification(&mut result, "structure_type", &types, &["other", "fcc", "hcp", "bcc", "icosahedral"]);
    result.set("cutoff_A", cutoff);
    Ok(result)
}

fn diamond_shell_type(shell: &[Vec3]) -> i32 {
    if shell.len() != 12 {
        return 0;
    }
    let mut mean = 0.0;
    for v in shell {
        mean += norm(*v);
    }
    let cutoff = mean / 12.0 * (1.0 + 2.0f64.sqrt()) / 2.0;
    let near = |a: usize, b: usize| norm(sub(shell[a], shell[b])) < cutoff;
    let (mut n421, mut n422) = (0, 0);
    for k in 0..12 {
        let common: Vec<usize> = (0..12).filter(|&m| m != k && near(k, m)).collect();
        if common.len() != 4 {
            return 0;
        }
        let mut bonds = Vec::new();
        for a in 0..4 {
            for b in a + 1..4 {
                if near(common[a], common[b]) {
                    bonds.push((a, b));
                }
            }
        }
        if bonds.len() != 2 {
            return 0;
        }
        let chained = bonds[0].0 == bonds[1].0 || bonds[0].0 == bonds[1].1 || bonds[0].1 == bonds[1].0 || bonds[0].1 == bonds[1].1;
        if chained {
            n422 += 1;
        } else {
            n421 += 1;
        }
    }
    if n421 == 12 {
        1
    } else if n421 == 6 && n422 == 6 {
        2
    } else {
        0
    }
}

pub fn diamond_structure(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    let pbc = p.pbc("pbc", [false; 3])?;
    let cutoff = automatic_cutoff(p, positions.len(), 2.0)?;
    let (neighbors, _) = local_neighbors(&positions, cutoff, p, "cell", pbc)?;
    let mut groups = neighbors_by_center(&neighbors, positions.len());
    let count = positions.len();
    let mut fourfold = vec![false; count];
    for i in 0..count {
        let group = &mut groups[i];
        stable_sort_by(group, |&a, &b| norm(neighbors[a].vector) < norm(neighbors[b].vector));
        if group.len() >= 4 {
            group.truncate(4);
            fourfold[i] = true;
        }
    }
    let mut types = vec![0i32; count];
    let mut second: Vec<Vec<usize>> = vec![Vec::new(); count];
    for i in 0..count {
        if !fourfold[i] {
            continue;
        }
        let mut shell = Vec::new();
        let mut complete = true;
        for &g in &groups[i] {
            let first = &neighbors[g];
            if !fourfold[first.j] {
                complete = false;
                break;
            }
            let mut back = 0;
            for &h in &groups[first.j] {
                let v = add(first.vector, neighbors[h].vector);
                if norm(v) < 0.25 * norm(first.vector) {
                    back += 1;
                    continue;
                }
                shell.push(v);
                second[i].push(neighbors[h].j);
            }
            if back != 1 {
                complete = false;
                break;
            }
        }
        if !complete {
            continue;
        }
        types[i] = match diamond_shell_type(&shell) {
            1 => 1,
            2 => 4,
            _ => 0,
        };
    }
    let mut labelled = types.clone();
    for i in 0..count {
        if types[i] == 1 || types[i] == 4 {
            for &g in &groups[i] {
                let j = neighbors[g].j;
                if labelled[j] == 0 {
                    labelled[j] = types[i] + 1;
                }
            }
        }
    }
    for i in 0..count {
        if types[i] == 1 || types[i] == 4 {
            for &j in &second[i] {
                if labelled[j] == 0 {
                    labelled[j] = types[i] + 2;
                }
            }
        }
    }
    let names = [
        "other",
        "cubic diamond",
        "cubic diamond (1st neighbor)",
        "cubic diamond (2nd neighbor)",
        "hexagonal diamond",
        "hexagonal diamond (1st neighbor)",
        "hexagonal diamond (2nd neighbor)",
    ];
    let mut result = Json::object();
    classification(&mut result, "diamond_type", &labelled, &names);
    result.set("cutoff_A", cutoff);
    Ok(result)
}

pub fn centrosymmetry(p: &Parameters) -> R {
    let neighbor_count = integer(p.number_or("neighbors", 12.0)?, "neighbors", 2)?;
    if neighbor_count % 2 != 0 {
        return Err("neighbors must be even".into());
    }
    let positions = points(p.array("positions")?, "positions")?;
    let cutoff = positive(p.number("cutoff_A")?, "cutoff_A")?;
    let (neighbors, _) = local_neighbors(&positions, cutoff, p, "cell", p.pbc("pbc", [false; 3])?)?;
    let groups = neighbors_by_center(&neighbors, positions.len());
    let n = neighbor_count as usize;
    let mut values = vec![f64::NAN; positions.len()];
    let mut valid = vec![false; positions.len()];
    for index in 0..positions.len() {
        if groups[index].len() < n {
            continue;
        }
        let mut bonds: Vec<Vec3> = groups[index].iter().map(|&g| neighbors[g].vector).collect();
        stable_sort_by(&mut bonds, |a, b| norm(*a) < norm(*b));
        bonds.truncate(n);
        let mut costs = Vec::new();
        for a in 0..n {
            for b in a + 1..n {
                let s = add(bonds[a], bonds[b]);
                costs.push(dot(s, s));
            }
        }
        sort_by(&mut costs, |a, b| a < b);
        values[index] = costs[..n / 2].iter().fold(0.0, |s, v| s + v);
        valid[index] = true;
    }
    let mut result = Json::object();
    result.set("centrosymmetry_A2", jv(&values));
    result.set("valid", jb(&valid));
    Ok(result)
}

pub fn bond_order(p: &Parameters) -> R {
    let mut orders: Vec<i32> = Vec::new();
    if p.has("degrees") {
        for degree in p.array("degrees")?.values {
            orders.push(integer(degree, "degree", 0)? as i32);
        }
    } else {
        orders = vec![4, 6];
    }
    let mut sorted = orders.clone();
    sorted.sort();
    if orders.is_empty() || sorted.windows(2).any(|w| w[0] == w[1]) {
        return Err("Supply distinct degrees".into());
    }
    if orders.iter().any(|&o| o > 32) {
        return Err("degree must be <= 32".into());
    }
    let positions = points(p.array("positions")?, "positions")?;
    let cutoff = positive(p.number("cutoff_A")?, "cutoff_A")?;
    let (neighbors, _) = local_neighbors(&positions, cutoff, p, "cell", p.pbc("pbc", [false; 3])?)?;
    let groups = neighbors_by_center(&neighbors, positions.len());
    let coordination: Vec<i32> = groups.iter().map(|g| g.len() as i32).collect();
    let valid: Vec<bool> = groups.iter().map(|g| !g.is_empty()).collect();
    let mut order_result = Json::object();
    for &l in &orders {
        let mut q = vec![f64::NAN; positions.len()];
        for index in 0..positions.len() {
            let group = &groups[index];
            if group.is_empty() {
                continue;
            }
            let mut total = 0.0;
            for m in 0..=l {
                let mut mean = Complex64::new(0.0, 0.0);
                for &g in group {
                    let bond = neighbors[g].vector;
                    let cosine = (bond[2] / norm(bond)).clamp(-1.0, 1.0);
                    let azimuth = bond[1].atan2(bond[0]);
                    let legendre = normalized_legendre(l, m, cosine);
                    let angle = m as f64 * azimuth;
                    let polar = Complex64::new(angle.cos(), angle.sin());
                    mean += Complex64::new(legendre * polar.re, legendre * polar.im);
                }
                mean = Complex64::new(mean.re / group.len() as f64, mean.im / group.len() as f64);
                total += (if m != 0 { 2.0 } else { 1.0 }) * mean.norm_sqr();
            }
            q[index] = (4.0 * PI / (2 * l + 1) as f64 * total).sqrt();
        }
        order_result.set(&format!("q{l}"), jv(&q));
    }
    let mut result = Json::object();
    result.set("order", order_result);
    result.set("coordination", ji(&coordination));
    result.set("valid", jb(&valid));
    Ok(result)
}

pub fn wigner_seitz(p: &Parameters) -> R {
    let sites = points(p.array("reference_sites")?, "reference_sites")?;
    let raw = p.array("positions")?;
    let mut atoms = Vec::new();
    if !(raw.size() == 0 && (raw.shape == vec![0] || raw.shape == vec![0, 3])) {
        atoms = points(raw, "positions")?;
    }
    let pbc = p.pbc("pbc", [false; 3])?;
    let periodic = pbc.iter().any(|&b| b);
    let mut cell = identity();
    if p.has("cell") {
        cell = cell_matrix(p.array("cell")?)?;
    } else if periodic || p.has("current_cell") {
        return Err("cell must be a nonempty finite array of the required dimension".into());
    }
    if p.has("current_cell") {
        let mapping = multiply(&inverse(&cell_matrix(p.array("current_cell")?)?)?, &cell);
        for a in atoms.iter_mut() {
            *a = row_times(*a, &mapping);
        }
    }
    let mic = MinimumImage::new(&cell, pbc)?;
    let (mut assignments, mut distances, mut ties) = (Vec::new(), Vec::new(), Vec::new());
    let mut occupancy = vec![0i32; sites.len()];
    for a in &atoms {
        let norms: Vec<f64> = sites.iter().map(|s| norm(mic.apply(sub(*s, *a)))).collect();
        let mut selected = 0;
        for s in 1..norms.len() {
            if norms[s] < norms[selected] {
                selected = s;
            }
        }
        let close = norms.iter().filter(|&&v| (v - norms[selected]).abs() <= 1e-10 + 1e-10 * norms[selected].abs()).count();
        assignments.push(selected as i32);
        distances.push(norms[selected]);
        ties.push(close > 1);
        occupancy[selected] += 1;
    }
    let (mut vacancy_sites, mut interstitial_sites) = (Vec::new(), Vec::new());
    let (mut vacancies, mut excess) = (0, 0);
    for (s, &o) in occupancy.iter().enumerate() {
        if o == 0 {
            vacancy_sites.push(s as i32);
            vacancies += 1;
        }
        if o > 1 {
            interstitial_sites.push(s as i32);
            excess += o - 1;
        }
    }
    let mut result = Json::object();
    result.set("occupancy", ji(&occupancy));
    result.set("site_index", ji(&assignments));
    result.set("distance_A", jv(&distances));
    result.set("ambiguous", jb(&ties));
    result.set("vacancy_sites", ji(&vacancy_sites));
    result.set("interstitial_sites", ji(&interstitial_sites));
    result.set("vacancies", vacancies as f64);
    result.set("interstitial_excess", excess as f64);
    Ok(result)
}

pub fn static_structure_factor(p: &Parameters) -> R {
    let mut c = finite_array(p.array("positions")?, "positions", -1)?;
    if c.ndim() == 2 {
        c.shape.insert(0, 1);
    }
    if c.ndim() != 3 || c.shape[2] != 3 {
        return Err("positions must have shape (atoms,3) or (frames,atoms,3)".into());
    }
    let wavevectors = points(p.array("q_vectors")?, "q_vectors")?;
    let (frames, atoms) = (c.shape[0], c.shape[1]);
    let factors = if p.has("weights") { p.array("weights")?.values } else { vec![1.0; atoms] };
    let mut total = 0.0;
    let mut finite = factors.len() == atoms;
    for &f in &factors {
        finite = finite && f.is_finite();
        total += f * f;
    }
    if !finite || total == 0.0 {
        return Err("weights must contain one finite scattering amplitude per atom, not all zero".into());
    }
    let mut intensity = vec![0.0; wavevectors.len()];
    for f in 0..frames {
        for (q, wv) in wavevectors.iter().enumerate() {
            let mut amplitude = Complex64::new(0.0, 0.0);
            for a in 0..atoms {
                let base = (f * atoms + a) * 3;
                let phase = wv[0] * c.values[base] + wv[1] * c.values[base + 1] + wv[2] * c.values[base + 2];
                amplitude += Complex64::new(factors[a] * phase.cos(), factors[a] * phase.sin());
            }
            intensity[q] += amplitude.norm_sqr();
        }
    }
    for value in intensity.iter_mut() {
        *value /= frames as f64 * total;
    }
    let mut result = Json::object();
    result.set("q_vectors_rad_per_A", Json::Array(wavevectors.iter().map(|q| vec3_json(*q)).collect()));
    result.set("S_q", jv(&intensity));
    Ok(result)
}

// ---------------------------------------------------------------- electronic
pub fn band_gap(p: &Parameters) -> R {
    let mut energy = finite_array(p.array("energies_eV")?, "energies_eV", -1)?;
    let fermi = p.number("fermi_eV")?;
    let tolerance = positive(p.number_or("tolerance_eV", 1e-6)?, "tolerance_eV")?;
    if energy.ndim() == 2 {
        energy.shape.insert(0, 1);
    }
    if energy.ndim() != 3 {
        return Err("energies_eV must have shape (k,band) or (spin,k,band)".into());
    }
    let (spins, kpoints, bands) = (energy.shape[0], energy.shape[1], energy.shape[2]);
    let e = |s: usize, k: usize, b: usize| energy.values[(s * kpoints + k) * bands + b];
    let (mut any_below, mut any_above, mut touching, mut crossing) = (false, false, false, false);
    let (mut valence, mut conduction) = (f64::NEG_INFINITY, f64::INFINITY);
    for &value in &energy.values {
        let below = value < fermi - tolerance;
        let above = value > fermi + tolerance;
        any_below = any_below || below;
        any_above = any_above || above;
        touching = touching || (value - fermi).abs() <= tolerance;
        if below {
            valence = fmax(valence, value);
        }
        if above {
            conduction = fmin(conduction, value);
        }
    }
    if !any_below || !any_above {
        return Err("Include both occupied and unoccupied bands".into());
    }
    for s in 0..spins {
        for b in 0..bands {
            let (mut below, mut above) = (false, false);
            for k in 0..kpoints {
                below = below || e(s, k, b) < fermi - tolerance;
                above = above || e(s, k, b) > fermi + tolerance;
            }
            crossing = crossing || (below && above);
        }
    }
    let metal = crossing || touching;
    let mut direct = f64::INFINITY;
    for s in 0..spins {
        for k in 0..kpoints {
            let (mut occupied, mut empty) = (f64::NEG_INFINITY, f64::INFINITY);
            for b in 0..bands {
                let value = e(s, k, b);
                if value < fermi - tolerance {
                    occupied = fmax(occupied, value);
                }
                if value > fermi + tolerance {
                    empty = fmin(empty, value);
                }
            }
            if occupied.is_finite() && empty.is_finite() {
                direct = fmin(direct, empty - occupied);
            }
        }
    }
    let (mut vbm, mut cbm) = (Vec::new(), Vec::new());
    for s in 0..spins {
        for k in 0..kpoints {
            for b in 0..bands {
                let value = e(s, k, b);
                let index = Json::Array(vec![Json::Number(s as f64), Json::Number(k as f64), Json::Number(b as f64)]);
                if value < fermi - tolerance && (value - valence).abs() <= tolerance {
                    vbm.push(index.clone());
                }
                if value > fermi + tolerance && (value - conduction).abs() <= tolerance {
                    cbm.push(index);
                }
            }
        }
    }
    let mut result = Json::object();
    result.set("metal_on_sampled_mesh", metal);
    result.set("gap_eV", if metal { 0.0 } else { conduction - valence });
    result.set(
        "direct_gap_eV",
        if metal {
            Json::Number(0.0)
        } else if direct.is_finite() {
            Json::Number(direct)
        } else {
            Json::Null
        },
    );
    result.set("vbm_eV", valence);
    result.set("cbm_eV", conduction);
    result.set("vbm_indices_spin_k_band", Json::Array(vbm));
    result.set("cbm_indices_spin_k_band", Json::Array(cbm));
    Ok(result)
}

pub fn effective_mass(p: &Parameters) -> R {
    let k = points(p.array("kpoints_inv_A")?, "kpoints_inv_A")?;
    let energy = vector1d(p.array("energies_eV")?, "energies_eV")?;
    let center = vector1d(p.array("center_inv_A")?, "center_inv_A")?;
    if center.len() != 3 || energy.len() != k.len() || k.len() < 10 {
        return Err("At least ten 3D k points, matching energies and one center are required".into());
    }
    let rows = k.len();
    let mut design = DMatrix::zeros(rows, 10);
    let mut values = DMatrix::zeros(rows, 1);
    for i in 0..rows {
        let (x, y, z) = (k[i][0] - center[0], k[i][1] - center[1], k[i][2] - center[2]);
        let row = [1.0, x, y, z, x * x / 2.0, y * y / 2.0, z * z / 2.0, x * y, x * z, y * z];
        for (j, v) in row.iter().enumerate() {
            design[(i, j)] = *v;
        }
        values[(i, 0)] = energy[i];
    }
    let scale: Vec<f64> = (0..10).map(|j| design.column(j).norm()).collect();
    let insufficient = "k points must span a full 3D quadratic fit; a line path is insufficient";
    if scale.iter().any(|&s| s == 0.0) {
        return Err(insufficient.into());
    }
    let scaled = DMatrix::from_fn(rows, 10, |i, j| design[(i, j)] * (1.0 / scale[j]));
    if matrix_rank(&scaled) != 10 {
        return Err(insufficient.into());
    }
    let solution = least_squares(&scaled, &values);
    let coefficients: Vec<f64> = (0..10).map(|j| (1.0 / scale[j]) * solution[(j, 0)]).collect();
    let c = &coefficients;
    let hessian = DMatrix::from_row_slice(3, 3, &[c[4], c[7], c[8], c[7], c[5], c[9], c[8], c[9], c[6]]);
    let (curvature, axes) = symmetric_eigen(&hessian);
    if curvature.iter().any(|v| v.abs() < 1e-10) {
        return Err("Singular band curvature: finite effective mass cannot be determined".into());
    }
    let constant = 7.619964231073853;
    let fitted = &design * DMatrix::from_column_slice(10, 1, c);
    let h3 = Matrix3::from_fn(|i, j| hessian[(i, j)]);
    let inverse = h3.try_inverse().unwrap_or_else(|| Matrix3::from_element(f64::NAN));
    let mass = DMatrix::from_fn(3, 3, |i, j| constant * inverse[(i, j)]);
    let mut result = Json::object();
    result.set("principal_masses_m_e", vector_json(&curvature.iter().map(|v| constant / v).collect::<Vec<_>>()));
    result.set("principal_axes_columns", matrix_json(&axes));
    result.set("mass_tensor_m_e", matrix_json(&mass));
    result.set("hessian_eV_A2", matrix_json(&hessian));
    result.set("gradient_eV_A", vector_json(&c[1..4]));
    result.set("energy_at_center_eV", c[0]);
    result.set("rms_fit_error_eV", ((&fitted - &values).norm_squared() / rows as f64).sqrt());
    Ok(result)
}

pub fn work_function(p: &Parameters) -> R {
    let position = vector1d(p.array("distance_A")?, "distance_A")?;
    let potential = vector1d(p.array("potential_eV")?, "potential_eV")?;
    let region = vector1d(p.array("vacuum_range_A")?, "vacuum_range_A")?;
    let fermi = p.number("fermi_eV")?;
    let limit = positive(p.number_or("max_slope_eV_per_A", 0.01)?, "max_slope_eV_per_A")?;
    let increasing = position.len() == potential.len() && position.windows(2).all(|w| w[1] > w[0]);
    if !increasing {
        return Err("Distances must increase and match the potential".into());
    }
    if region.len() != 2 || region[0] >= region[1] || region[0] < position[0] || region[1] > *position.last().unwrap() {
        return Err("Vacuum interval must lie within the supplied profile".into());
    }
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for i in 0..position.len() {
        if position[i] >= region[0] && position[i] <= region[1] {
            x.push(position[i]);
            y.push(potential[i]);
        }
    }
    if x.len() < 3 {
        return Err("At least three vacuum samples are required".into());
    }
    let slope = linear_regression(&x, &y).slope;
    let vacuum = y.iter().fold(0.0, |s, v| s + v) / y.len() as f64;
    let mut variance = 0.0;
    for v in &y {
        variance += (v - vacuum) * (v - vacuum);
    }
    let mut result = Json::object();
    result.set("work_function_eV", vacuum - fermi);
    result.set("vacuum_level_eV", vacuum);
    result.set("vacuum_std_eV", (variance / y.len() as f64).sqrt());
    result.set("slope_eV_per_A", slope);
    result.set("flat_vacuum", slope.abs() <= limit);
    result.set("samples", x.len() as f64);
    result.set("fermi_eV", fermi);
    result.set("distance_A", jv(&position));
    result.set("potential_eV", jv(&potential));
    Ok(result)
}

// ---------------------------------------------------------------- fits
fn birch_murnaghan(volume: f64, q: &Vector4<f64>, gradient: Option<&mut Vector4<f64>>) -> f64 {
    let x = (q[3] / volume).powf(2.0 / 3.0) - 1.0;
    let shape = 2.0 * x * x + (q[2] - 4.0) * x * x * x;
    let prefactor = 9.0 / 16.0;
    if let Some(g) = gradient {
        g[0] = 1.0;
        g[1] = prefactor * q[3] * shape;
        g[2] = prefactor * q[1] * q[3] * x * x * x;
        g[3] = prefactor * q[1] * (shape + 2.0 / 3.0 * (x + 1.0) * (4.0 * x + 3.0 * (q[2] - 4.0) * x * x));
    }
    q[0] + prefactor * q[1] * q[3] * shape
}

pub fn equation_of_state(p: &Parameters) -> R {
    let volumes = vector1d(p.array("volumes_A3")?, "volumes_A3")?;
    let energies = vector1d(p.array("energies_eV")?, "energies_eV")?;
    let mut unique = volumes.clone();
    sort_by(&mut unique, |a, b| a < b);
    let distinct = !unique.windows(2).any(|w| w[0] == w[1]);
    if volumes.len() != energies.len() || volumes.len() < 5 || unique[0] <= 0.0 || !distinct {
        return Err("Provide at least five distinct positive volumes and matching energies".into());
    }
    let n = volumes.len();
    let design = DMatrix::from_fn(n, 3, |i, j| match j {
        0 => volumes[i] * volumes[i],
        1 => volumes[i],
        _ => 1.0,
    });
    let e = DMatrix::from_column_slice(n, 1, &energies);
    let quadratic = least_squares(&design, &e);
    let (a, b, c) = (quadratic[(0, 0)], quadratic[(1, 0)], quadratic[(2, 0)]);
    if a <= 0.0 {
        return Err("Energy-volume samples must contain a convex minimum".into());
    }
    let (vmin, vmax) = (unique[0], *unique.last().unwrap());
    let v0 = -b / (2.0 * a);
    if !(vmin < v0 && v0 < vmax) {
        return Err("Volume samples must bracket the equilibrium minimum".into());
    }
    let mut q = Vector4::new(a * v0 * v0 + b * v0 + c, 2.0 * a * v0, 4.0, v0);
    let lower = Vector4::new(f64::NEG_INFINITY, 1e-12, -20.0, vmin);
    let upper = Vector4::new(f64::INFINITY, f64::INFINITY, 30.0, vmax);
    let evaluate = |parameters: &Vector4<f64>| {
        let mut residual = DVector::zeros(n);
        let mut jacobian = DMatrix::zeros(n, 4);
        for i in 0..n {
            let mut gradient = Vector4::zeros();
            residual[i] = energies[i] - birch_murnaghan(volumes[i], parameters, Some(&mut gradient));
            for k in 0..4 {
                jacobian[(i, k)] = gradient[k];
            }
        }
        let cost = residual.norm_squared();
        (residual, jacobian, cost)
    };
    let (mut residual, mut jacobian, mut cost) = evaluate(&q);
    let mut lambda = 1e-3;
    let mut iteration = 0i64;
    while iteration < 20000 {
        let normal: Matrix4<f64> = Matrix4::from_fn(|i, j| (jacobian.transpose() * &jacobian)[(i, j)]);
        let jt_r = jacobian.transpose() * &residual;
        let gradient = Vector4::new(jt_r[0], jt_r[1], jt_r[2], jt_r[3]);
        let mut improved = false;
        while lambda < 1e20 {
            let mut damped = normal;
            for k in 0..4 {
                damped[(k, k)] += lambda * fmax(normal[(k, k)], 1e-30);
            }
            let step = damped.cholesky().map(|ch| ch.solve(&gradient)).or_else(|| damped.lu().solve(&gradient)).unwrap_or(Vector4::from_element(f64::NAN));
            let mut candidate = q + step;
            for k in 0..4 {
                candidate[k] = candidate[k].clamp(lower[k], upper[k]);
            }
            let (candidate_residual, candidate_jacobian, candidate_cost) = evaluate(&candidate);
            if candidate_cost.is_finite() && candidate_cost < cost {
                let change = (candidate - q).abs().max();
                let relative = (cost - candidate_cost) / fmax(cost, 1e-300);
                q = candidate;
                residual = candidate_residual;
                jacobian = candidate_jacobian;
                cost = candidate_cost;
                lambda = fmax(lambda / 10.0, 1e-15);
                improved = true;
                if change < 1e-15 * (1.0 + q.abs().max()) || relative < 1e-15 {
                    iteration = 1 << 30;
                }
                break;
            }
            lambda *= 10.0;
        }
        if !improved || cost == 0.0 {
            break;
        }
        iteration += 1;
    }
    let (e0, modulus, derivative, volume) = (q[0], q[1], q[2], q[3]);
    if !(vmin + 1e-8 < volume && volume < vmax - 1e-8) {
        return Err("Fitted equilibrium lies on the volume boundary; sample a wider range".into());
    }
    let residuals: Vec<f64> = (0..n).map(|i| energies[i] - birch_murnaghan(volumes[i], &q, None)).collect();
    let mut squared = 0.0;
    for r in &residuals {
        squared += r * r;
    }
    let normal = jacobian.transpose() * &jacobian;
    let mut covariance = DMatrix::from_element(4, 4, f64::INFINITY);
    if n > 4 && matrix_rank(&normal) == 4 {
        if let Some(inv) = normal.clone().try_inverse() {
            covariance = inv * (squared / (n - 4) as f64);
        }
    }
    let fit_volumes: Vec<f64> = (0..=200).map(|i| vmin + (vmax - vmin) * i as f64 / 200.0).collect();
    let fit_energies: Vec<f64> = fit_volumes.iter().map(|&v| birch_murnaghan(v, &q, None)).collect();
    let mut result = Json::object();
    result.set("volume_A3", volume);
    result.set("energy_eV", e0);
    result.set("bulk_modulus_GPa", modulus * EV_PER_A3_TO_GPA);
    result.set("bulk_derivative", derivative);
    result.set("residuals_eV", jv(&residuals));
    result.set("rms_error_eV", (squared / residuals.len() as f64).sqrt());
    result.set("parameter_covariance", matrix_json(&covariance));
    result.set("volumes_A3", jv(&volumes));
    result.set("energies_eV", jv(&energies));
    result.set("fit_volumes_A3", jv(&fit_volumes));
    result.set("fit_energies_eV", jv(&fit_energies));
    Ok(result)
}

pub fn elastic_tensor(p: &Parameters) -> R {
    let strain = finite_array(p.array("strains")?, "strains", 2)?;
    let stress = finite_array(p.array("stresses_GPa")?, "stresses_GPa", 2)?;
    if strain.shape != stress.shape || strain.shape[1] != 6 || strain.shape[0] < 7 {
        return Err("Provide matching (samples>=7,6) engineering strains and stresses".into());
    }
    let n = strain.shape[0];
    let design = DMatrix::from_fn(n, 7, |i, j| if j == 0 { 1.0 } else { strain.values[i * 6 + j - 1] });
    let target = DMatrix::from_fn(n, 6, |i, j| stress.values[i * 6 + j]);
    if matrix_rank(&design) != 7 {
        return Err("Strains must independently sample all six strain components".into());
    }
    let coefficients = least_squares(&design, &target);
    let raw = DMatrix::from_fn(6, 6, |i, j| coefficients[(j + 1, i)]);
    let stiffness = (&raw + raw.transpose()) / 2.0;
    let (eigenvalues, _) = symmetric_eigen(&stiffness);
    let stable = eigenvalues.iter().all(|&v| v > 0.0);
    let mut moduli = Json::Null;
    if stable {
        let compliance = stiffness.clone().try_inverse().unwrap_or_else(|| DMatrix::from_element(6, 6, f64::NAN));
        let diagonal = |m: &DMatrix<f64>| {
            [m[(0, 0)] + m[(1, 1)] + m[(2, 2)], m[(0, 1)] + m[(0, 2)] + m[(1, 2)], m[(3, 3)] + m[(4, 4)] + m[(5, 5)]]
        };
        let [a, b, c] = diagonal(&stiffness);
        let bulk_v = (a + 2.0 * b) / 9.0;
        let shear_v = (a - b + 3.0 * c) / 15.0;
        let s = diagonal(&compliance);
        let bulk_r = 1.0 / (s[0] + 2.0 * s[1]);
        let shear_r = 15.0 / (4.0 * s[0] - 4.0 * s[1] + 3.0 * s[2]);
        let bulk = (bulk_v + bulk_r) / 2.0;
        let shear = (shear_v + shear_r) / 2.0;
        let mut m = Json::object();
        m.set("bulk_voigt_GPa", bulk_v);
        m.set("bulk_reuss_GPa", bulk_r);
        m.set("shear_voigt_GPa", shear_v);
        m.set("shear_reuss_GPa", shear_r);
        m.set("bulk_hill_GPa", bulk);
        m.set("shear_hill_GPa", shear);
        m.set("young_hill_GPa", 9.0 * bulk * shear / (3.0 * bulk + shear));
        m.set("poisson_hill", (3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear)));
        m.set("universal_anisotropy", 5.0 * shear_v / shear_r + bulk_v / bulk_r - 6.0);
        moduli = m;
    }
    let prestress: Vec<f64> = (0..6).map(|j| coefficients[(0, j)]).collect();
    let prediction = DMatrix::from_fn(n, 6, |i, j| {
        let mut v = design[(i, 0)] * coefficients[(0, j)];
        let mut lin = 0.0;
        for k in 0..6 {
            lin += design[(i, k + 1)] * stiffness[(j, k)];
        }
        v += lin;
        v
    });
    let residual = &target - prediction;
    let mut result = Json::object();
    result.set("stiffness_GPa", matrix_json(&stiffness));
    result.set("prestress_GPa", vector_json(&prestress));
    result.set("stable_zero_prestress", stable);
    result.set("stiffness_eigenvalues_GPa", vector_json(&eigenvalues));
    result.set("antisymmetric_norm_GPa", (&raw - raw.transpose()).norm());
    result.set("rms_stress_error_GPa", (residual.norm_squared() / residual.len() as f64).sqrt());
    result.set("moduli", moduli);
    Ok(result)
}

// ---------------------------------------------------------------- phonons
struct Modes {
    energies: NdArray,
    weights: Vec<f64>,
}

fn modes(p: &Parameters) -> Result<Modes, String> {
    let energies = finite_array(p.array("energies_eV")?, "energies_eV", 2)?;
    let rows = energies.shape[0];
    let mut weights = if p.has("weights") { vector1d(p.array("weights")?, "weights")? } else { vec![1.0; rows] };
    let mut total = 0.0;
    let mut valid = weights.len() == rows;
    for &w in &weights {
        valid = valid && w >= 0.0;
        total += w;
    }
    if !valid || total <= 0.0 {
        return Err("Supply one nonnegative q-point weight per row, with positive total".into());
    }
    for w in weights.iter_mut() {
        *w /= total;
    }
    Ok(Modes { energies, weights })
}

pub fn phonon_dos(p: &Parameters) -> R {
    let m = modes(p)?;
    let grid = vector1d(p.array("energy_grid_eV")?, "energy_grid_eV")?;
    if !(grid.len() >= 2 && grid.windows(2).all(|w| w[1] > w[0])) {
        return Err("Energy grid must contain at least two increasing samples".into());
    }
    let sigma = positive(p.number_or("sigma_eV", 0.001)?, "sigma_eV")?;
    let mut density = vec![0.0; grid.len()];
    let mut imaginary = 0.0;
    let (rows, columns) = (m.energies.shape[0], m.energies.shape[1]);
    for q in 0..rows {
        for b in 0..columns {
            let mode = m.energies.values[q * columns + b];
            if mode < 0.0 {
                imaginary += m.weights[q];
            }
            for g in 0..grid.len() {
                let z = (grid[g] - mode) / sigma;
                density[g] += m.weights[q] * (-0.5 * (z * z)).exp() / (sigma * (2.0 * PI).sqrt());
            }
        }
    }
    let mut result = Json::object();
    result.set("energy_eV", jv(&grid));
    result.set("dos_per_eV", jv(&density));
    result.set("total_modes", columns as f64);
    result.set("enclosed_modes", trapezoid(&density, &grid));
    result.set("imaginary_weight", imaginary);
    Ok(result)
}

pub fn harmonic_thermodynamics(p: &Parameters) -> R {
    let m = modes(p)?;
    let temperatures = vector1d(p.array("temperatures_K")?, "temperatures_K")?;
    let tolerance = positive(p.number_or("zero_tolerance_eV", 1e-8)?, "zero_tolerance_eV")?;
    let valid = temperatures.iter().all(|&t| t >= 0.0) && m.energies.values.iter().all(|&e| e >= -tolerance);
    if !valid {
        return Err("Nonnegative temperatures and dynamically stable phonons are required".into());
    }
    let (rows, columns) = (m.energies.shape[0], m.energies.shape[1]);
    let (mut mode_energies, mut mode_weights) = (Vec::new(), Vec::new());
    let mut omitted = 0.0;
    for q in 0..rows {
        for b in 0..columns {
            let e = m.energies.values[q * columns + b];
            if e > tolerance {
                mode_energies.push(e);
                mode_weights.push(m.weights[q]);
            } else {
                omitted += m.weights[q];
            }
        }
    }
    if mode_energies.is_empty() {
        return Err("No positive-frequency modes remain".into());
    }
    let mut zpe = 0.0;
    for i in 0..mode_energies.len() {
        zpe += mode_weights[i] * mode_energies[i] / 2.0;
    }
    let (mut free, mut internal, mut entropy, mut capacity) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for &temperature in &temperatures {
        let (mut f, mut u, mut s, mut cv) = (zpe, zpe, 0.0, 0.0);
        if temperature > 0.0 {
            let (mut log_sum, mut energy_sum, mut heat) = (0.0, 0.0, 0.0);
            for i in 0..mode_energies.len() {
                let x = mode_energies[i] / (BOLTZMANN * temperature);
                let decay = (-x).exp();
                let denominator = -(-x).exp_m1();
                log_sum += mode_weights[i] * denominator.ln();
                energy_sum += mode_weights[i] * mode_energies[i] * decay / denominator;
                if x < 700.0 {
                    heat += mode_weights[i] * x * x * decay / (denominator * denominator);
                }
            }
            f = zpe + BOLTZMANN * temperature * log_sum;
            u = zpe + energy_sum;
            s = (u - f) / temperature;
            cv = BOLTZMANN * heat;
        }
        free.push(f);
        internal.push(u);
        entropy.push(s);
        capacity.push(cv);
    }
    let mut result = Json::object();
    result.set("temperature_K", jv(&temperatures));
    result.set("free_energy_eV", jv(&free));
    result.set("internal_energy_eV", jv(&internal));
    result.set("entropy_eV_per_K", jv(&entropy));
    result.set("heat_capacity_eV_per_K", jv(&capacity));
    result.set("zero_point_energy_eV", zpe);
    result.set("omitted_zero_mode_weight", omitted);
    Ok(result)
}
