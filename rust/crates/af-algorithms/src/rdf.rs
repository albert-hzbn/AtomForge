//! Radial distribution function (the C++ algorithms/RadialDistributionAnalysis):
//! g(r) over all periodic images, first peak, minimum and width, and the
//! first-shell distortion metrics. Single precision as the C++.

use af_core::glm;
use af_core::neighbours;
use af_core::Structure;
use glam::Vec3;

const PI: f32 = std::f32::consts::PI;
const MIN_DISTANCE: f32 = 1e-6;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RdfBin {
    pub r_center: f32,
    pub g: f32,
    pub raw_count: f32,
    pub cumulative: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RdfResult {
    pub valid: bool,
    pub message: String,
    pub pbc_used: bool,
    pub normalized: bool,
    pub atom_count: i32,
    pub ref_count: i32,
    pub target_count: i32,
    pub volume: f32,
    pub density: f32,
    pub r_min: f32,
    pub r_max: f32,
    pub bin_count: i32,
    pub bin_width: f32,
    pub first_peak_r: f32,
    pub first_peak_value: f32,
    pub first_minimum_r: f32,
    pub first_minimum_value: f32,
    pub has_first_peak: bool,
    pub has_first_minimum: bool,
    pub has_distortion_metrics: bool,
    pub distortion_window_auto: bool,
    pub distortion_window_min: f32,
    pub distortion_window_max: f32,
    pub shell_mean_distance: f32,
    pub shell_std_dev: f32,
    pub shell_relative_distortion: f32,
    pub shell_coordination: f32,
    pub first_peak_fwhm: f32,
    pub has_first_peak_fwhm: bool,
    pub bins: Vec<RdfBin>,
}

impl Default for RdfResult {
    fn default() -> Self {
        Self {
            valid: false,
            message: String::new(),
            pbc_used: false,
            normalized: false,
            atom_count: 0,
            ref_count: 0,
            target_count: 0,
            volume: 0.0,
            density: 0.0,
            r_min: 0.0,
            r_max: 0.0,
            bin_count: 0,
            bin_width: 0.0,
            first_peak_r: 0.0,
            first_peak_value: 0.0,
            first_minimum_r: 0.0,
            first_minimum_value: 0.0,
            has_first_peak: false,
            has_first_minimum: false,
            has_distortion_metrics: false,
            distortion_window_auto: true,
            distortion_window_min: 0.0,
            distortion_window_max: 0.0,
            shell_mean_distance: 0.0,
            shell_std_dev: 0.0,
            shell_relative_distortion: 0.0,
            shell_coordination: 0.0,
            first_peak_fwhm: 0.0,
            has_first_peak_fwhm: false,
            bins: Vec::new(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RdfParams {
    pub ref_atomic_number_filter: i32,
    pub target_atomic_number_filter: i32,
    pub use_pbc: bool,
    pub normalize: bool,
    pub r_min: f32,
    pub r_max: f32,
    pub bin_count: i32,
    pub smoothing_passes: i32,
    pub enable_distortion_analysis: bool,
    pub auto_distortion_window: bool,
    pub manual_distortion_min: f32,
    pub manual_distortion_max: f32,
}

impl Default for RdfParams {
    fn default() -> Self {
        Self {
            ref_atomic_number_filter: 0,
            target_atomic_number_filter: 0,
            use_pbc: true,
            normalize: true,
            r_min: 0.0,
            r_max: 8.0,
            bin_count: 200,
            smoothing_passes: 0,
            enable_distortion_analysis: true,
            auto_distortion_window: true,
            manual_distortion_min: 1.8,
            manual_distortion_max: 3.2,
        }
    }
}

/// The volume of the atoms' bounding box (each side at least 1 A).
pub fn bounding_volume(structure: &Structure) -> f32 {
    let Some(first) = structure.atoms.first() else { return 0.0 };
    let p = |a: &af_core::AtomSite| Vec3::new(a.x as f32, a.y as f32, a.z as f32);
    let mut min_p = p(first);
    let mut max_p = min_p;
    for a in &structure.atoms[1..] {
        min_p = glm::min(min_p, p(a));
        max_p = glm::max(max_p, p(a));
    }
    let e = glm::max(max_p - min_p, Vec3::ONE);
    e.x * e.y * e.z
}

pub fn compute_rdf(structure: &Structure, p: &RdfParams) -> RdfResult {
    if !p.r_min.is_finite() || !p.r_max.is_finite() || p.r_min < 0.0 || !p.manual_distortion_min.is_finite() || !p.manual_distortion_max.is_finite() {
        return RdfResult { message: "Radius parameters must be finite and r_min non-negative.".into(), ..Default::default() };
    }
    let (r_min, r_max, bin_count) = (p.r_min, p.r_max, p.bin_count);
    let mut result = RdfResult {
        atom_count: structure.atoms.len() as i32,
        normalized: p.normalize,
        r_min,
        r_max,
        bin_count,
        ..Default::default()
    };
    if structure.atoms.is_empty() {
        result.message = "No atoms available.".into();
        return result;
    }
    if r_max <= r_min {
        result.message = "r_max must be greater than r_min.".into();
        return result;
    }
    if bin_count < 8 {
        result.message = "Bin count too small.".into();
        return result;
    }
    let mut use_pbc = false;
    let mut volume = 0.0f32;
    if p.use_pbc && structure.has_unit_cell {
        if let Some((cell, _)) = glm::cell_matrices(structure) {
            volume = glm::determinant(&cell).abs();
            use_pbc = true;
        }
    }
    if !use_pbc {
        volume = bounding_volume(structure);
    }
    result.pbc_used = use_pbc;
    result.volume = volume;
    result.bin_width = (r_max - r_min) / bin_count as f32;

    let is_ref: Vec<bool> =
        structure.atoms.iter().map(|a| p.ref_atomic_number_filter <= 0 || a.atomic_number == p.ref_atomic_number_filter).collect();
    let is_target: Vec<bool> =
        structure.atoms.iter().map(|a| p.target_atomic_number_filter <= 0 || a.atomic_number == p.target_atomic_number_filter).collect();
    let ref_count = is_ref.iter().filter(|&&b| b).count();
    let target_count = is_target.iter().filter(|&&b| b).count();
    result.ref_count = ref_count as i32;
    result.target_count = target_count as i32;
    if ref_count == 0 || target_count == 0 {
        result.message = "Reference or target species selection is empty.".into();
        return result;
    }
    if volume <= 0.0 {
        result.message = "Unable to determine analysis volume.".into();
        return result;
    }
    let n_bins = bin_count as usize;
    let mut histogram = vec![0.0f32; n_bins];
    let bin_width = result.bin_width as f64;
    // Every periodic image within r_max (up to 200 images of the cell).
    let visited = neighbours::visit_within(structure, r_max as f64, use_pbc, 200, |pair| {
        if !is_ref[pair.i] || !is_target[pair.j] {
            return;
        }
        let r = pair.r;
        if r < r_min as f64 || r >= r_max as f64 || r <= MIN_DISTANCE as f64 {
            return;
        }
        let bin = ((r - r_min as f64) / bin_width) as i32;
        if bin >= 0 && bin < bin_count {
            histogram[bin as usize] += 1.0;
        }
    });
    if let Err(e) = visited {
        result.message = e.to_string();
        return result;
    }
    for _ in 0..p.smoothing_passes.max(0) {
        let mut smoothed = histogram.clone();
        for i in 1..n_bins - 1 {
            smoothed[i] = 0.25 * histogram[i - 1] + 0.5 * histogram[i] + 0.25 * histogram[i + 1];
        }
        histogram = smoothed;
    }
    let rho_target = target_count as f32 / volume;
    result.density = rho_target;
    let mut cumulative = 0.0f32;
    for i in 0..n_bins {
        let r0 = r_min + i as f32 * result.bin_width;
        let r1 = r0 + result.bin_width;
        let rc = 0.5 * (r0 + r1);
        let shell_volume = (4.0f32 / 3.0) * PI * (r1 * r1 * r1 - r0 * r0 * r0);
        let expected = ref_count as f32 * rho_target * shell_volume;
        let g = if p.normalize && expected > 1e-12 { histogram[i] / expected } else { histogram[i] };
        cumulative += histogram[i] / ref_count as f32;
        result.bins.push(RdfBin { r_center: rc, raw_count: histogram[i], g, cumulative });
    }
    let bins = &result.bins;
    let mut best_peak = -f32::MAX;
    let mut peak_index: i32 = -1;
    for i in 1..n_bins - 1 {
        let y = bins[i].g;
        if y >= bins[i - 1].g && y >= bins[i + 1].g && y > best_peak {
            best_peak = y;
            peak_index = i as i32;
        }
    }
    let mut first_peak = None;
    if peak_index >= 0 {
        let pi = peak_index as usize;
        let peak_value = bins[pi].g;
        let mut minimum = None;
        for i in pi + 1..n_bins - 1 {
            let y = bins[i].g;
            if y <= bins[i - 1].g && y <= bins[i + 1].g {
                minimum = Some(i);
                break;
            }
        }
        let half = 0.5 * peak_value;
        let mut left = bins[pi].r_center;
        let mut right = bins[pi].r_center;
        let interpolate = |y0: f32, y1: f32, r0: f32, r1: f32| {
            let t = if (y1 - y0).abs() > 1e-12 { (half - y0) / (y1 - y0) } else { 0.0 };
            r0 + t * (r1 - r0)
        };
        for i in (1..=pi).rev() {
            let (y0, y1) = (bins[i - 1].g, bins[i].g);
            if y0 <= half && y1 >= half {
                left = interpolate(y0, y1, bins[i - 1].r_center, bins[i].r_center);
                break;
            }
        }
        for i in pi..n_bins - 1 {
            let (y0, y1) = (bins[i].g, bins[i + 1].g);
            if y0 >= half && y1 <= half {
                right = interpolate(y0, y1, bins[i].r_center, bins[i + 1].r_center);
                break;
            }
        }
        first_peak = Some((pi, peak_value, minimum, left, right));
    }
    if let Some((pi, peak_value, minimum, left, right)) = first_peak {
        result.has_first_peak = true;
        result.first_peak_r = result.bins[pi].r_center;
        result.first_peak_value = peak_value;
        if let Some(m) = minimum {
            result.has_first_minimum = true;
            result.first_minimum_r = result.bins[m].r_center;
            result.first_minimum_value = result.bins[m].g;
        }
        if right > left {
            result.has_first_peak_fwhm = true;
            result.first_peak_fwhm = right - left;
        }
    }
    if p.enable_distortion_analysis {
        let mut shell_min = p.manual_distortion_min;
        let mut shell_max = p.manual_distortion_max;
        let mut has_window = false;
        let bins = &result.bins;
        if p.auto_distortion_window {
            if result.has_first_peak {
                let pi = peak_index as usize;
                let is_minimum = |i: usize| bins[i].g <= bins[i - 1].g && bins[i].g <= bins[i + 1].g;
                let left_min = (1..pi).rev().find(|&i| is_minimum(i));
                let right_min = (pi + 1..n_bins - 1).find(|&i| is_minimum(i));
                if let Some(r) = right_min {
                    shell_min = left_min.map_or(r_min, |l| bins[l].r_center);
                    shell_max = bins[r].r_center;
                    has_window = shell_max > shell_min;
                }
            }
        } else {
            has_window = shell_max > shell_min;
        }
        if has_window {
            let in_window = |b: &&RdfBin| !(b.r_center < shell_min || b.r_center > shell_max);
            let (mut sum_w, mut sum_wr, mut raw_sum) = (0.0f32, 0.0f32, 0.0f32);
            for b in bins.iter().filter(in_window) {
                sum_w += b.raw_count;
                sum_wr += b.raw_count * b.r_center;
                raw_sum += b.raw_count;
            }
            if sum_w > 1e-8 && result.ref_count > 0 {
                let mean = sum_wr / sum_w;
                let mut var = 0.0f32;
                for b in bins.iter().filter(in_window) {
                    let d = b.r_center - mean;
                    var += b.raw_count * d * d;
                }
                var /= sum_w;
                let std_dev = (if var < 0.0 { 0.0 } else { var }).sqrt();
                result.has_distortion_metrics = true;
                result.distortion_window_auto = p.auto_distortion_window;
                result.distortion_window_min = shell_min;
                result.distortion_window_max = shell_max;
                result.shell_mean_distance = mean;
                result.shell_std_dev = std_dev;
                result.shell_relative_distortion = if mean > 1e-8 { (std_dev / mean) * 100.0 } else { 0.0 };
                result.shell_coordination = raw_sum / result.ref_count as f32;
            }
        }
    }
    result.valid = true;
    result.message = if p.normalize { "RDF computed.".into() } else { "Radial histogram computed.".into() };
    result
}
