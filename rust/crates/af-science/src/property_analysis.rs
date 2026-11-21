//! Analyses of per-atom properties (the C++ science/PropertyAnalysis):
//! spatial binning, histograms, spatial correlation functions (by FFT),
//! scatter plots and ambient occlusion.

use af_core::cxx::sort_by;
use af_core::json::Json;
use num_complex::Complex64;

use crate::analysis::{jv, R};
use crate::core::*;
use crate::data::Parameters;

fn fmin(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}
fn fmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

/// Per-atom values (1 for every atom when absent).
fn per_atom(p: &Parameters, name: &str, count: usize) -> Result<Vec<f64>, String> {
    if !p.has(name) {
        return Ok(vec![1.0; count]);
    }
    let array = p.array(name)?;
    if array.ndim() != 1 {
        return Err(format!("{name} must be a list of numbers"));
    }
    if array.values.len() != count {
        return Err(format!("{name} needs one value per atom"));
    }
    Ok(array.values)
}

fn bin_counts(p: &Parameters, fallback: [i64; 3]) -> Result<[i64; 3], String> {
    if !p.has("bins") {
        return Ok(fallback);
    }
    let values = p.array("bins")?.values;
    if values.len() != 3 {
        return Err("bins needs three counts (1 = not binned along that cell vector)".into());
    }
    let mut out = [0i64; 3];
    for d in 0..3 {
        out[d] = integer(values[d], "bins", 1)?;
    }
    if out[0] * out[1] * out[2] > 50_000_000 {
        return Err("Too many bins".into());
    }
    Ok(out)
}

fn text_or(p: &Parameters, name: &str, fallback: &str) -> Result<String, String> {
    Ok(if p.has(name) { p.json(name)?.string()?.to_string() } else { fallback.to_string() })
}

/// 3D DFT of a row-major grid, axis by axis, in place.
fn fft3(grid: &mut [Complex64], n: [usize; 3], inverse: bool) {
    let stride = [n[1] * n[2], n[2], 1];
    for axis in 0..3 {
        if n[axis] == 1 {
            continue;
        }
        let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
        for i in 0..n[a] {
            for j in 0..n[b] {
                let base = i * stride[a] + j * stride[b];
                let line: Vec<Complex64> = (0..n[axis]).map(|k| grid[base + k * stride[axis]]).collect();
                let line = fft(line, inverse);
                for k in 0..n[axis] {
                    grid[base + k * stride[axis]] = line[k];
                }
            }
        }
    }
}

pub fn spatial_binning(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    if !p.has("cell") {
        return Err("Spatial binning needs the cell".into());
    }
    let cell = cell_matrix(p.array("cell")?)?;
    let pbc = p.pbc("pbc", [true; 3])?;
    let values = per_atom(p, "values", positions.len())?;
    let bins = bin_counts(p, [20, 1, 1])?;
    let operation = text_or(p, "operation", "mean")?;
    if !["mean", "sum", "density", "min", "max", "count"].contains(&operation.as_str()) {
        return Err("operation must be mean, sum, density, min, max or count".into());
    }
    let total = (bins[0] * bins[1] * bins[2]) as usize;
    let mut sum = vec![0.0; total];
    let mut low = vec![f64::INFINITY; total];
    let mut high = vec![f64::NEG_INFINITY; total];
    let mut count = vec![0i64; total];
    let mut outside = 0usize;
    for (i, x) in positions.iter().enumerate() {
        let mut f = fractional(*x, &cell)?;
        let mut index = [0i64; 3];
        let mut inside = true;
        for d in 0..3 {
            if pbc[d] {
                f[d] -= f[d].floor();
                if f[d] >= 1.0 {
                    f[d] = 0.0;
                }
            }
            if f[d] < 0.0 || f[d] >= 1.0 {
                inside = false;
                break;
            }
            index[d] = (bins[d] - 1).min((f[d] * bins[d] as f64) as i64);
        }
        if !inside {
            outside += 1;
            continue;
        }
        if values[i].is_nan() {
            continue;
        }
        let k = ((index[0] * bins[1] + index[1]) * bins[2] + index[2]) as usize;
        sum[k] += values[i];
        low[k] = fmin(low[k], values[i]);
        high[k] = fmax(high[k], values[i]);
        count[k] += 1;
    }
    let bin_volume = cell_volume(&cell) / total as f64;
    let out: Vec<f64> = (0..total)
        .map(|k| match operation.as_str() {
            "mean" => {
                if count[k] > 0 {
                    sum[k] / count[k] as f64
                } else {
                    f64::NAN
                }
            }
            "sum" => sum[k],
            "density" => sum[k] / bin_volume,
            "min" => {
                if count[k] > 0 {
                    low[k]
                } else {
                    f64::NAN
                }
            }
            "max" => {
                if count[k] > 0 {
                    high[k]
                } else {
                    f64::NAN
                }
            }
            _ => count[k] as f64,
        })
        .collect();
    let mut result = Json::object();
    result.set("bins", jv(&[bins[0] as f64, bins[1] as f64, bins[2] as f64]));
    result.set("values", jv(&out));
    result.set("atom_counts", jv(&count.iter().map(|&c| c as f64).collect::<Vec<_>>()));
    result.set("operation", operation.as_str());
    result.set("bin_volume_A3", bin_volume);
    result.set("atoms_outside", outside as f64);
    let mut centers = Json::object();
    for (d, axis_name) in ["a", "b", "c"].iter().enumerate() {
        if bins[d] == 1 {
            continue;
        }
        let length = norm(cell[d]);
        let fractional_centers: Vec<f64> = (0..bins[d]).map(|b| (b as f64 + 0.5) / bins[d] as f64).collect();
        let distance: Vec<f64> = fractional_centers.iter().map(|f| f * length).collect();
        let mut axis = Json::object();
        axis.set("fractional", jv(&fractional_centers));
        axis.set("distance_A", jv(&distance));
        centers.set(axis_name, axis);
    }
    result.set("centers", centers);
    Ok(result)
}

pub fn property_histogram(p: &Parameters) -> R {
    let input = p.array("values")?;
    if input.ndim() != 1 {
        return Err("values must be a list of numbers".into());
    }
    let mut values = input.values;
    if p.has("mask") {
        let mask = finite_array(p.array("mask")?, "mask", 1)?.values;
        if mask.len() != values.len() {
            return Err("mask needs one value per atom".into());
        }
        values = values.iter().zip(&mask).filter(|(_, &m)| m != 0.0).map(|(&v, _)| v).collect();
    }
    values.retain(|v| !v.is_nan());
    if values.is_empty() {
        return Err("No values to analyse".into());
    }
    let bins = integer(p.number_or("bins", 50.0)?, "bins", 1)?;
    let mut sorted = values.clone();
    sort_by(&mut sorted, |a, b| a < b);
    let (mut low, mut high) = (sorted[0], *sorted.last().unwrap());
    if p.has("range") {
        let range = p.array("range")?.values;
        if range.len() != 2 || !(range[1] > range[0]) {
            return Err("range must be [low, high] with high > low".into());
        }
        low = range[0];
        high = range[1];
    }
    if high - low <= 1e-9 * fmax(1.0, high.abs()) {
        low -= 0.5;
        high += 0.5;
    }
    let width = (high - low) / bins as f64;
    let mut counts = vec![0.0; bins as usize];
    let centers: Vec<f64> = (0..bins).map(|b| low + (b as f64 + 0.5) * width).collect();
    let mut outside = 0usize;
    for &v in &values {
        if v < low || v > high {
            outside += 1;
            continue;
        }
        let b = (bins - 1).min(((v - low) / width) as i64);
        counts[b as usize] += 1.0;
    }
    let n = values.len() as f64;
    let sum = values.iter().fold(0.0, |s, v| s + v);
    let mean = sum / n;
    let mut squares = 0.0;
    for &v in &values {
        squares += (v - mean) * (v - mean);
    }
    let half = sorted.len() / 2;
    let mut result = Json::object();
    result.set("centers", jv(&centers));
    result.set("counts", jv(&counts));
    result.set("bin_width", width);
    result.set("values_outside_range", outside as f64);
    result.set("count", n);
    result.set("sum", sum);
    result.set("mean", mean);
    result.set("std", (squares / n).sqrt());
    result.set("min", sorted[0]);
    result.set("max", *sorted.last().unwrap());
    result.set("median", if sorted.len() % 2 == 1 { sorted[half] } else { 0.5 * (sorted[half - 1] + sorted[half]) });
    Ok(result)
}

pub fn spatial_correlation(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    if !p.has("cell") {
        return Err("The spatial correlation function needs a periodic cell".into());
    }
    let cell = cell_matrix(p.array("cell")?)?;
    let atoms = positions.len();
    let a = per_atom(p, "values_a", atoms)?;
    let b = if p.has("values_b") { per_atom(p, "values_b", atoms)? } else { a.clone() };
    let spacing = positive(p.number_or("grid_spacing_A", 1.0)?, "grid_spacing_A")?;
    let radial_bins = integer(p.number_or("bins", 50.0)?, "bins", 1)?;
    let volume = cell_volume(&cell);
    let mut height = f64::INFINITY;
    for d in 0..3 {
        height = fmin(height, volume / norm(cross(cell[(d + 1) % 3], cell[(d + 2) % 3])));
    }
    let r_max = p.number_or("r_max_A", height / 2.0)?;
    if !(r_max > 0.0) {
        return Err("r_max_A must be positive".into());
    }
    let n: [usize; 3] = std::array::from_fn(|d| 1usize.max((norm(cell[d]) / spacing).round() as usize));
    let cells = n[0] * n[1] * n[2];
    if cells > 32_000_000 {
        return Err("FFT grid too fine; increase grid_spacing_A".into());
    }
    let zero = Complex64::new(0.0, 0.0);
    let (mut grid_a, mut grid_b, mut grid_n) = (vec![zero; cells], vec![zero; cells], vec![zero; cells]);
    let mut used = 0usize;
    for i in 0..atoms {
        if a[i].is_nan() || b[i].is_nan() {
            continue;
        }
        used += 1;
        let mut f = fractional(positions[i], &cell)?;
        let mut g = [0usize; 3];
        for d in 0..3 {
            f[d] -= f[d].floor();
            g[d] = ((f[d] * n[d] as f64).round() as i64 as usize) % n[d];
        }
        let k = (g[0] * n[1] + g[1]) * n[2] + g[2];
        grid_a[k].re += a[i];
        grid_b[k].re += b[i];
        grid_n[k].re += 1.0;
    }
    fft3(&mut grid_a, n, false);
    fft3(&mut grid_b, n, false);
    fft3(&mut grid_n, n, false);
    let mut ab: Vec<Complex64> = (0..cells).map(|k| grid_a[k].conj() * grid_b[k]).collect();
    let mut nn: Vec<Complex64> = (0..cells).map(|k| grid_n[k].conj() * grid_n[k]).collect();
    fft3(&mut ab, n, true);
    fft3(&mut nn, n, true);
    let bins = radial_bins as usize;
    let (mut numerator, mut pairs) = (vec![0.0; bins], vec![0.0; bins]);
    let width = r_max / radial_bins as f64;
    let wrapped = |index: usize, size: usize| {
        let v = index as f64 / size as f64;
        if v >= 0.5 {
            v - 1.0
        } else {
            v
        }
    };
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                let f = [wrapped(i, n[0]), wrapped(j, n[1]), wrapped(k, n[2])];
                let r = norm(row_times(f, &cell));
                let bin = (r / width * (1.0 + 1e-12) + 1e-12) as usize;
                if bin >= bins {
                    continue;
                }
                let index = (i * n[1] + j) * n[2] + k;
                numerator[bin] += ab[index].re;
                pairs[bin] += nn[index].re;
            }
        }
    }
    if used == 0 {
        return Err("No atom has values".into());
    }
    let (mut mean_a, mut mean_b, mut mean_ab) = (0.0, 0.0, 0.0);
    for i in 0..atoms {
        if !a[i].is_nan() && !b[i].is_nan() {
            mean_a += a[i];
            mean_b += b[i];
            mean_ab += a[i] * b[i];
        }
    }
    mean_a /= used as f64;
    mean_b /= used as f64;
    mean_ab /= used as f64;
    let covariance = mean_ab - mean_a * mean_b;
    let radius: Vec<f64> = (0..bins).map(|bin| (bin as f64 + 0.5) * width).collect();
    let (mut correlation, mut normalized) = (vec![f64::NAN; bins], vec![f64::NAN; bins]);
    for bin in 0..bins {
        if pairs[bin] > 0.5 {
            correlation[bin] = numerator[bin] / pairs[bin];
            if covariance.abs() > 0.0 {
                normalized[bin] = (correlation[bin] - mean_a * mean_b) / covariance;
            }
        }
    }
    let mut result = Json::object();
    result.set("r_A", jv(&radius));
    result.set("correlation", jv(&correlation));
    result.set("normalized_correlation", jv(&normalized));
    result.set("pair_counts", jv(&pairs));
    result.set("mean_a", mean_a);
    result.set("mean_b", mean_b);
    result.set("covariance", covariance);
    result.set("grid", jv(&[n[0] as f64, n[1] as f64, n[2] as f64]));
    result.set("r_max_A", r_max);
    Ok(result)
}

pub fn scatter_plot(p: &Parameters) -> R {
    let (xs, ys) = (p.array("values_x")?, p.array("values_y")?);
    if xs.ndim() != 1 || ys.ndim() != 1 || xs.values.len() != ys.values.len() {
        return Err("values_x and values_y must be lists of the same length".into());
    }
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for i in 0..xs.values.len() {
        if xs.values[i].is_finite() && ys.values[i].is_finite() {
            x.push(xs.values[i]);
            y.push(ys.values[i]);
        }
    }
    if x.len() < 2 {
        return Err("At least two atoms need both values".into());
    }
    let n = x.len() as f64;
    let mx = x.iter().fold(0.0, |s, v| s + v) / n;
    let my = y.iter().fold(0.0, |s, v| s + v) / n;
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    for i in 0..x.len() {
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
        sxy += (x[i] - mx) * (y[i] - my);
    }
    let mut result = Json::object();
    result.set("x", jv(&x));
    result.set("y", jv(&y));
    result.set("points", n);
    result.set("pearson_r", if sxx > 0.0 && syy > 0.0 { sxy / (sxx * syy).sqrt() } else { f64::NAN });
    result.set("slope", if sxx > 0.0 { sxy / sxx } else { f64::NAN });
    result.set("intercept", if sxx > 0.0 { my - sxy / sxx * mx } else { f64::NAN });
    Ok(result)
}

pub fn ambient_occlusion(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    let n = positions.len();
    if n == 0 {
        return Err("No atoms".into());
    }
    let mut radius = vec![0.0; n];
    if p.has("radii") {
        let radii = p.array("radii")?;
        if radii.values.len() != n {
            return Err("radii needs one value per atom".into());
        }
        radius = radii.values;
    } else {
        let (mut low, mut high) = (positions[0], positions[0]);
        for x in &positions {
            for d in 0..3 {
                low[d] = fmin(low[d], x[d]);
                high[d] = fmax(high[d], x[d]);
            }
        }
        let mut volume = 1.0;
        for d in 0..3 {
            volume *= fmax(high[d] - low[d], 1.0);
        }
        radius = vec![0.55 * (volume / n as f64).cbrt(); n];
    }
    if radius.iter().any(|r| !(*r > 0.0)) {
        return Err("radii must be positive".into());
    }
    let directions = integer(p.number_or("directions", 40.0)?, "directions", 1)?;
    let resolution = integer(p.number_or("resolution", 256.0)?, "resolution", 16)?;
    let golden = std::f64::consts::PI * (3.0 - 5.0f64.sqrt());
    let mut light = vec![0.0; n];
    let mut max_radius = radius[0];
    for &r in &radius[1..] {
        if max_radius < r {
            max_radius = r;
        }
    }
    for k in 0..directions {
        let z = 1.0 - 2.0 * (k as f64 + 0.5) / directions as f64;
        let ring = fmax(0.0, 1.0 - z * z).sqrt();
        let dir = [ring * (golden * k as f64).cos(), ring * (golden * k as f64).sin(), z];
        let helper = if dir[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let u = scale(cross(dir, helper), 1.0 / norm(cross(dir, helper)));
        let v = cross(dir, u);
        let (mut low_u, mut low_v, mut high_u, mut high_v) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for x in &positions {
            low_u = fmin(low_u, dot(*x, u));
            high_u = fmax(high_u, dot(*x, u));
            low_v = fmin(low_v, dot(*x, v));
            high_v = fmax(high_v, dot(*x, v));
        }
        let extent = fmax(high_u - low_u, high_v - low_v) + 2.0 * max_radius;
        let pixel = extent / resolution as f64;
        let side = resolution as usize;
        let mut depth = vec![f64::NEG_INFINITY; side * side];
        let mut owner = vec![-1i64; side * side];
        for i in 0..n {
            let cu = (dot(positions[i], u) - low_u + max_radius) / pixel;
            let cv = (dot(positions[i], v) - low_v + max_radius) / pixel;
            let r = radius[i] / pixel;
            let h = dot(positions[i], dir);
            let u0 = 0i64.max((cu - r).floor() as i64);
            let u1 = (side as i64 - 1).min((cu + r).ceil() as i64);
            let v0 = 0i64.max((cv - r).floor() as i64);
            let v1 = (side as i64 - 1).min((cv + r).ceil() as i64);
            for a in u0..=u1 {
                for b in v0..=v1 {
                    let du = a as f64 + 0.5 - cu;
                    let dv = b as f64 + 0.5 - cv;
                    let inside = r * r - du * du - dv * dv;
                    if inside < 0.0 {
                        continue;
                    }
                    let front = h + inside.sqrt() * pixel;
                    let k2 = a as usize * side + b as usize;
                    if front > depth[k2] {
                        depth[k2] = front;
                        owner[k2] = i as i64;
                    }
                }
            }
        }
        for &o in &owner {
            if o >= 0 {
                light[o as usize] += pixel * pixel;
            }
        }
    }
    let mut brightness = vec![0.0; n];
    let mut best = 0.0;
    for i in 0..n {
        brightness[i] = light[i] / (radius[i] * radius[i]);
        best = fmax(best, brightness[i]);
    }
    for b in brightness.iter_mut() {
        *b = if best > 0.0 { *b / best } else { 0.0 };
    }
    let mut result = Json::object();
    result.set("brightness", jv(&brightness));
    result.set("directions", directions as f64);
    Ok(result)
}
