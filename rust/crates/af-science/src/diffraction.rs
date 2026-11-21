//! Kinematic powder X-ray and electron diffraction (the C++
//! science/Diffraction).

use std::collections::BTreeMap;
use std::f64::consts::{LN_2, PI};

use af_core::cxx::sort_by;
use af_core::elements;
use af_core::json::Json;

use crate::analysis::jv;
use crate::core::*;
use crate::data::{Parameters, StructureInput};
use crate::scattering_factors::{Coefficients, TABLE};
use crate::symmetry::conventional_index_transform;
use crate::tools::ToolOutput;

type Hkl = [i32; 3];

fn coefficients(symbol: &str) -> Result<&'static Coefficients, String> {
    TABLE.iter().find(|c| c.symbol == symbol).ok_or_else(|| format!("No scattering factor coefficients for {symbol}"))
}

struct Crystal {
    cell: Mat3,
    /// Rows b_i with a_i . b_j = delta_ij (no 2 pi).
    reciprocal: Mat3,
    symbols: Vec<String>,
    fractional: Vec<Vec3>,
}

fn crystal_from(input: &StructureInput) -> Result<Crystal, String> {
    let s = &input.structure;
    if !s.has_unit_cell || s.atoms.is_empty() {
        return Err("Diffraction requires a nonempty periodic structure".into());
    }
    let cell = s.cell_vectors;
    let mut c = Crystal { cell, reciprocal: transpose(&inverse(&cell)?), symbols: Vec::new(), fractional: Vec::new() };
    for atom in &s.atoms {
        c.symbols.push(atom.symbol.clone());
        c.fractional.push(fractional([atom.x, atom.y, atom.z], &cell)?);
        coefficients(&atom.symbol)?;
    }
    Ok(c)
}

fn hkl_vec(hkl: &Hkl) -> Vec3 {
    [hkl[0] as f64, hkl[1] as f64, hkl[2] as f64]
}

/// The structure factor (real, imaginary).
fn structure_factor(c: &Crystal, hkl: &Hkl, s: f64, factor: fn(&str, f64) -> f64, debye_waller: f64) -> (f64, f64) {
    let (mut re, mut im) = (0.0, 0.0);
    let mut cache: BTreeMap<&str, f64> = BTreeMap::new();
    for (symbol, f) in c.symbols.iter().zip(&c.fractional) {
        let value = *cache.entry(symbol).or_insert_with(|| factor(symbol, s) * (-debye_waller * s * s).exp());
        let phase = 2.0 * PI * (hkl[0] as f64 * f[0] + hkl[1] as f64 * f[1] + hkl[2] as f64 * f[2]);
        re += value * phase.cos();
        im += value * phase.sin();
    }
    (re, im)
}

fn norm2((re, im): (f64, f64)) -> f64 {
    re * re + im * im
}

/// All (h, k, l) with |g| <= g_max (g = h b1 + k b2 + l b3, 1/Angstrom).
fn reflections(c: &Crystal, g_max: f64) -> Vec<Hkl> {
    let range: [i32; 3] = std::array::from_fn(|k| (g_max * norm(c.cell[k])).ceil() as i32 + 1);
    let mut result = Vec::new();
    for h in -range[0]..=range[0] {
        for k in -range[1]..=range[1] {
            for l in -range[2]..=range[2] {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                if norm(row_times([h as f64, k as f64, l as f64], &c.reciprocal)) <= g_max {
                    result.push([h, k, l]);
                }
            }
        }
    }
    result
}

fn hkl_json(hkl: &Hkl) -> Json {
    Json::Array(hkl.iter().map(|&v| Json::from(v)).collect())
}

fn scattering_sum(c: &Coefficients, s: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..4 {
        sum += c.a[i] * (-c.b[i] * s * s).exp();
    }
    sum
}

pub fn xray_scattering_factor(symbol: &str, s: f64) -> f64 {
    let c = coefficients(symbol).expect("checked by crystal_from");
    let sum = scattering_sum(c, s);
    elements::atomic_number(symbol) as f64 - 41.78214 * s * s * sum
}

/// Mott-Bethe: f_e = 0.023934 (Z - f_x) / s^2 = 0.023934 * 41.78214 * sum a_i exp(-b_i s^2).
pub fn electron_scattering_factor(symbol: &str, s: f64) -> f64 {
    let c = coefficients(symbol).expect("checked by crystal_from");
    0.023934 * 41.78214 * scattering_sum(c, s)
}

/// The representative of equivalent reflections reads most simply: fewest
/// negative indices, then the smallest |h|+|k|+|l|, then the largest indices
/// first. Returns the first maximum (std::max_element).
fn canonical(list: &[Hkl]) -> Hkl {
    let key = |h: &Hkl| {
        let negatives = (h[0] < 0) as i32 + (h[1] < 0) as i32 + (h[2] < 0) as i32;
        let sum = h[0].abs() + h[1].abs() + h[2].abs();
        (-negatives, -sum, h[0], h[1], h[2])
    };
    let mut best = list[0];
    for h in &list[1..] {
        if key(&best) < key(h) {
            best = *h;
        }
    }
    best
}

#[derive(Default)]
struct Peak {
    two_theta: f64,
    d: f64,
    intensity: f64,
    hkls: Vec<Hkl>,
}

pub fn powder_xrd(input: &StructureInput, p: &Parameters) -> Result<ToolOutput, String> {
    let c = crystal_from(input)?;
    let wavelength = positive(p.number_or("wavelength_A", 1.54184)?, "wavelength_A")?;
    let (mut low, mut high) = (10.0, 90.0);
    if p.has("two_theta_range_deg") {
        let values = finite_array(p.array("two_theta_range_deg")?, "two_theta_range_deg", 1)?.values;
        if values.len() != 2 || values[0] < 0.0 || values[1] > 180.0 || values[0] >= values[1] {
            return Err("two_theta_range_deg must be [low, high] within 0-180".into());
        }
        low = values[0];
        high = values[1];
    }
    let debye_waller = p.number_or("debye_waller_A2", 0.0)?;
    let fwhm = positive(p.number_or("fwhm_deg", 0.1)?, "fwhm_deg")?;
    let g_max = 2.0 * (high / 2.0 * PI / 180.0).sin() / wavelength;
    // Merge reflections sharing a Bragg angle (within 1e-5 degrees).
    let mut peaks: BTreeMap<i64, Peak> = BTreeMap::new();
    let all = reflections(&c, g_max);
    // Reflections whose structure factor vanishes are not counted in a peak's multiplicity.
    let mut electrons = 0.0;
    for atom in &input.structure.atoms {
        electrons += atom.atomic_number as f64;
    }
    let vanishing = 1e-10 * 1f64.max(electrons * electrons);
    for hkl in &all {
        let g = norm(row_times(hkl_vec(hkl), &c.reciprocal));
        let theta = (wavelength * g / 2.0).min(1.0).asin();
        let two_theta = 2.0 * theta * 180.0 / PI;
        if two_theta < low || two_theta > high {
            continue;
        }
        let f = structure_factor(&c, hkl, g / 2.0, xray_scattering_factor, debye_waller);
        let cos2 = (2.0 * theta).cos();
        let sin = theta.sin();
        let lorentz = (1.0 + cos2 * cos2) / (sin * sin * theta.cos());
        let key = (two_theta / 1e-5).round() as i64;
        let peak = peaks.entry(key).or_default();
        if peak.hkls.is_empty() {
            peak.two_theta = two_theta;
            peak.d = 1.0 / g;
        }
        peak.intensity += norm2(f) * lorentz;
        if norm2(f) > vanishing {
            peak.hkls.push(*hkl);
        }
    }
    let mut maximum = 0.0f64;
    for peak in peaks.values() {
        if maximum < peak.intensity {
            maximum = peak.intensity;
        }
    }
    // Indices of the conventional standard cell too, so a supercell's peaks read as the crystal's.
    let numbers: Vec<i32> = input.structure.atoms.iter().map(|a| a.atomic_number).collect();
    let to_conventional = conventional_index_transform(&c.cell, &c.fractional, &numbers, 1e-3);
    let mut list = Vec::new();
    let mut kept: Vec<(f64, f64)> = Vec::new();
    for peak in peaks.values() {
        let scaled = if maximum > 0.0 { 100.0 * peak.intensity / maximum } else { 0.0 };
        if scaled < 1e-3 || peak.hkls.is_empty() {
            continue;
        }
        let mut entry = Json::object();
        entry.set("two_theta_deg", peak.two_theta);
        entry.set("d_A", peak.d);
        entry.set("intensity", scaled);
        entry.set("hkl", hkl_json(&canonical(&peak.hkls)));
        entry.set("multiplicity", peak.hkls.len());
        if let Some(t) = &to_conventional {
            let mut indices = Vec::new();
            let mut whole = true;
            for hkl in &peak.hkls {
                let h = row_times(hkl_vec(hkl), t);
                let mut rounded = [0i32; 3];
                for k in 0..3 {
                    rounded[k] = h[k].round() as i32;
                    whole = whole && (h[k] - rounded[k] as f64).abs() < 1e-6;
                }
                indices.push(rounded);
            }
            if whole && !indices.is_empty() {
                entry.set("hkl_conventional", hkl_json(&canonical(&indices)));
            }
        }
        list.push(entry);
        kept.push((peak.two_theta, scaled));
    }
    // Gaussian-broadened profile on a 0.02 degree grid.
    let sigma = fwhm / (2.0 * (2.0 * LN_2).sqrt());
    let (mut x, mut y) = (Vec::new(), Vec::new());
    let mut t = low;
    while t <= high + 1e-9 {
        let mut value = 0.0;
        for &(position, intensity) in &kept {
            if (t - position).abs() < 6.0 * sigma {
                let u = (t - position) / sigma;
                value += intensity * (-0.5 * (u * u)).exp();
            }
        }
        x.push(t);
        y.push(value);
        t += 0.02;
    }
    let mut result = Json::object();
    result.set("wavelength_A", wavelength);
    let count = list.len();
    result.set("peaks", Json::Array(list));
    result.set("peak_count", count);
    result.set("two_theta_deg", jv(&x));
    result.set("profile", jv(&y));
    Ok(ToolOutput { result, ..Default::default() })
}

pub fn electron_diffraction(input: &StructureInput, p: &Parameters) -> Result<ToolOutput, String> {
    let c = crystal_from(input)?;
    let mut axis = vec![0.0, 0.0, 1.0];
    if p.has("zone_axis") {
        axis = finite_array(p.array("zone_axis")?, "zone_axis", 1)?.values;
    }
    if axis.len() != 3 {
        return Err("zone_axis needs three indices [u, v, w]".into());
    }
    let beam = row_times([axis[0], axis[1], axis[2]], &c.cell);
    if norm(beam) < 1e-9 {
        return Err("zone_axis cannot be zero".into());
    }
    let voltage = positive(p.number_or("voltage_kV", 200.0)?, "voltage_kV")? * 1000.0;
    // Relativistic electron wavelength (Angstrom).
    let wavelength = 12.2642598 / (voltage * (1.0 + 0.978475e-6 * voltage)).sqrt();
    let g_max = positive(p.number_or("g_max_inv_A", 1.0)?, "g_max_inv_A")?;
    let direction = scale(beam, 1.0 / norm(beam));
    let mut spots: Vec<(Hkl, Vec3, f64)> = Vec::new();
    let mut maximum = 0.0f64;
    for hkl in reflections(&c, g_max) {
        if (hkl[0] as f64 * axis[0] + hkl[1] as f64 * axis[1] + hkl[2] as f64 * axis[2]).abs() > 1e-9 {
            continue;
        }
        let g = row_times(hkl_vec(&hkl), &c.reciprocal);
        let intensity = norm2(structure_factor(&c, &hkl, norm(g) / 2.0, electron_scattering_factor, 0.0));
        if maximum < intensity {
            maximum = intensity;
        }
        spots.push((hkl, g, intensity));
    }
    sort_by(&mut spots, |a, b| norm(a.1) < norm(b.1));
    // Detector axes: e1 along the shortest allowed reflection, e2 = beam x e1.
    let mut e1 = [0.0; 3];
    for spot in &spots {
        if spot.2 > 1e-6 * maximum {
            e1 = scale(spot.1, 1.0 / norm(spot.1));
            break;
        }
    }
    if norm(e1) == 0.0 {
        return Err("No allowed reflections in this zone within g_max".into());
    }
    let e2 = cross(direction, e1);
    let mut list = Vec::new();
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for (hkl, g, intensity) in &spots {
        let scaled = 100.0 * intensity / maximum;
        if scaled < 1e-3 {
            continue;
        }
        let mut entry = Json::object();
        entry.set("hkl", hkl_json(hkl));
        entry.set("g_inv_A", norm(*g));
        entry.set("d_A", 1.0 / norm(*g));
        entry.set("x_inv_A", dot(*g, e1));
        entry.set("y_inv_A", dot(*g, e2));
        entry.set("intensity", scaled);
        entry.set("angle_mrad", 2.0 * (wavelength * norm(*g) / 2.0).asin() * 1000.0);
        list.push(entry);
        xs.push(dot(*g, e1));
        ys.push(dot(*g, e2));
    }
    let mut result = Json::object();
    result.set("wavelength_A", wavelength);
    result.set("zone_axis", jv(&axis));
    let count = list.len();
    result.set("spots", Json::Array(list));
    result.set("spot_count", count);
    result.set("spot_x_inv_A", jv(&xs));
    result.set("spot_y_inv_A", jv(&ys));
    Ok(ToolOutput { result, ..Default::default() })
}
