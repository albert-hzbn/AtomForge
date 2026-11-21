//! Bond-angle distribution (the C++ algorithms/AngularDistributionAnalysis):
//! angles j-i-k over neighbour pairs within a cutoff, smoothed, with peaks
//! labelled by geometry and per-element coordination statistics.

use std::collections::BTreeMap;

use af_core::glm::{self, dot, length};
use af_core::neighbours;
use af_core::{AtomSite, Structure};
use glam::Vec3;

const PI: f32 = std::f32::consts::PI;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AdfCentreMode {
    /// Any atom is a centre.
    #[default]
    All,
    /// Centres of `centre_symbol` only.
    ByElement,
    /// Explicit i, j, k types.
    ByPair,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdfParams {
    pub r_cutoff: f32,
    pub bin_count: i32,
    pub use_pbc: bool,
    pub centre_mode: AdfCentreMode,
    pub centre_symbol: String,
    pub neigh_symbol1: String,
    pub neigh_symbol2: String,
    /// Normalise to peak = 1.
    pub normalize: bool,
    pub smooth_passes: i32,
}

impl Default for AdfParams {
    fn default() -> Self {
        Self {
            r_cutoff: 3.5,
            bin_count: 180,
            use_pbc: true,
            centre_mode: AdfCentreMode::All,
            centre_symbol: String::new(),
            neigh_symbol1: String::new(),
            neigh_symbol2: String::new(),
            normalize: true,
            smooth_passes: 2,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AdfBin {
    pub angle_deg: f32,
    pub count: f32,
    pub value: f32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AdfPeakInfo {
    pub angle_deg: f32,
    pub value: f32,
    pub label: String,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CoordStat {
    pub mean: f32,
    pub stddev: f32,
    pub count: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AdfResult {
    pub valid: bool,
    pub message: String,
    pub r_cutoff: f32,
    pub bin_count: i32,
    pub bin_width: f32,
    pub pbc_used: bool,
    pub normalized: bool,
    pub n_atoms: i32,
    pub n_centre_atoms: i32,
    pub n_triplets: i64,
    pub element_counts: BTreeMap<String, i32>,
    pub bins: Vec<AdfBin>,
    pub peaks: Vec<AdfPeakInfo>,
    pub coord_stats: BTreeMap<String, CoordStat>,
}

fn angle_deg(rij: Vec3, rik: Vec3) -> f32 {
    let d1 = length(rij);
    let d2 = length(rik);
    if d1 < 1e-6 || d2 < 1e-6 {
        return 0.0;
    }
    let mut cos_a = dot(rij, rik) / (d1 * d2);
    // std::max(-1, std::min(1, cosA))
    cos_a = if cos_a < 1.0 { cos_a } else { 1.0 };
    cos_a = if -1.0 < cos_a { cos_a } else { -1.0 };
    cos_a.acos() * (180.0 / PI)
}

/// Box smoothing (1/4, 1/2, 1/4), reflecting at the ends.
fn smooth_histogram(v: &mut Vec<f32>, passes: i32) {
    if v.len() < 3 {
        return;
    }
    let n = v.len();
    for _ in 0..passes {
        let mut tmp = vec![0.0f32; n];
        for i in 0..n {
            let mut s = v[i] * 0.5;
            s += v[if i > 0 { i - 1 } else { i }] * 0.25;
            s += v[if i + 1 < n { i + 1 } else { i }] * 0.25;
            tmp[i] = s;
        }
        *v = tmp;
    }
}

fn max_element(v: &[f32]) -> f32 {
    let mut best = v[0];
    for &x in &v[1..] {
        if best < x {
            best = x;
        }
    }
    best
}

/// Local maxima above `threshold` of the global maximum.
fn find_peaks(v: &[f32], threshold: f32) -> Vec<usize> {
    let gmax = max_element(v);
    if gmax < 1e-10 {
        return Vec::new();
    }
    (1..v.len().saturating_sub(1)).filter(|&i| v[i] > v[i - 1] && v[i] >= v[i + 1] && v[i] > threshold * gmax).collect()
}

fn geometry_label(deg: f32) -> &'static str {
    const KNOWN: [(f32, &str); 8] = [
        (60.0, "triangular 60\u{b0}"),
        (70.5, "FCC/HCP 70.5\u{b0}"),
        (90.0, "octahedral 90\u{b0}"),
        (109.47, "tetrahedral 109.5\u{b0}"),
        (120.0, "trigonal planar 120\u{b0}"),
        (135.0, "square antiprism 135\u{b0}"),
        (150.0, "BCC 150\u{b0}"),
        (180.0, "linear 180\u{b0}"),
    ];
    KNOWN.iter().find(|(angle, _)| (deg - angle).abs() < 3.5).map_or("", |(_, label)| label)
}

pub fn compute_adf(structure: &Structure, params: &AdfParams) -> AdfResult {
    let mut result = AdfResult {
        r_cutoff: params.r_cutoff,
        bin_count: params.bin_count,
        normalized: params.normalize,
        n_atoms: structure.atoms.len() as i32,
        ..Default::default()
    };
    if structure.atoms.is_empty() {
        result.message = "No atoms in structure.".into();
        return result;
    }
    if !params.r_cutoff.is_finite()
        || params.r_cutoff < 0.1
        || params.bin_count < 1
        || params.bin_count > 100000
        || params.smooth_passes < 0
        || params.smooth_passes > 1000
    {
        result.message = "Use a finite cutoff >= 0.1 A, 1-100000 bins and 0-1000 smoothing passes.".into();
        return result;
    }
    result.bin_width = 180.0 / params.bin_count as f32;
    if structure.atoms.iter().any(|a| !a.x.is_finite() || !a.y.is_finite() || !a.z.is_finite()) {
        result.message = "Atom coordinates must be finite.".into();
        return result;
    }
    for a in &structure.atoms {
        *result.element_counts.entry(a.symbol.clone()).or_default() += 1;
    }
    let use_pbc = params.use_pbc && structure.has_unit_cell && glm::cell_matrices(structure).is_some();
    result.pbc_used = use_pbc;
    let is_centre = |a: &AtomSite| params.centre_mode == AdfCentreMode::All || a.symbol == params.centre_symbol;
    let is_neigh1 = |a: &AtomSite| params.centre_mode != AdfCentreMode::ByPair || params.neigh_symbol1.is_empty() || a.symbol == params.neigh_symbol1;
    let is_neigh2 = |a: &AtomSite| params.centre_mode != AdfCentreMode::ByPair || params.neigh_symbol2.is_empty() || a.symbol == params.neigh_symbol2;

    let n_bins = params.bin_count as usize;
    let mut hist = vec![0.0f32; n_bins];
    let mut coord_counts: BTreeMap<String, Vec<i32>> = BTreeMap::new();
    let mut n_triplets = 0i64;
    let pairs = match neighbours::within(structure, params.r_cutoff as f64, use_pbc) {
        Ok(p) => p,
        Err(e) => {
            result.message = e.to_string();
            return result;
        }
    };
    let mut next = 0;
    for (i, ai) in structure.atoms.iter().enumerate() {
        let begin = next;
        while next < pairs.len() && pairs[next].i == i {
            next += 1;
        }
        if !is_centre(ai) {
            continue;
        }
        let neighs: Vec<(usize, Vec3)> =
            pairs[begin..next].iter().map(|p| (p.j, Vec3::new(p.d.x as f32, p.d.y as f32, p.d.z as f32))).collect();
        coord_counts.entry(ai.symbol.clone()).or_default().push(neighs.len() as i32);
        // Each unordered neighbour pair once, independent of atom order.
        for a in 0..neighs.len() {
            for b in a + 1..neighs.len() {
                let first = &structure.atoms[neighs[a].0];
                let second = &structure.atoms[neighs[b].0];
                if !((is_neigh1(first) && is_neigh2(second)) || (is_neigh1(second) && is_neigh2(first))) {
                    continue;
                }
                let angle = angle_deg(neighs[a].1, neighs[b].1);
                let mut bin = (angle / result.bin_width) as i32;
                if bin >= params.bin_count {
                    bin = params.bin_count - 1;
                }
                hist[bin as usize] += 1.0;
                n_triplets += 1;
            }
        }
    }
    result.n_triplets = n_triplets;
    result.n_centre_atoms = structure.atoms.iter().filter(|a| is_centre(a)).count() as i32;
    if n_triplets == 0 {
        result.message = "No triplets found \u{2014} try a larger cutoff radius.".into();
        return result;
    }
    let raw_counts = hist.clone();
    smooth_histogram(&mut hist, params.smooth_passes);
    let gmax = max_element(&hist);
    for b in 0..n_bins {
        result.bins.push(AdfBin {
            angle_deg: (b as f32 + 0.5) * result.bin_width,
            count: raw_counts[b],
            value: if params.normalize && gmax > 1e-10 { hist[b] / gmax } else { hist[b] },
        });
    }
    let plot: Vec<f32> = result.bins.iter().map(|b| b.value).collect();
    for idx in find_peaks(&plot, 0.08) {
        let angle = result.bins[idx].angle_deg;
        result.peaks.push(AdfPeakInfo { angle_deg: angle, value: result.bins[idx].value, label: geometry_label(angle).to_string() });
    }
    for (symbol, counts) in &coord_counts {
        let (mut sum, mut sum2) = (0.0f32, 0.0f32);
        for &c in counts {
            sum += c as f32;
            sum2 += (c * c) as f32;
        }
        let mean = sum / counts.len() as f32;
        let var = sum2 / counts.len() as f32 - mean * mean;
        let stddev = (if 0.0 < var { var } else { 0.0 }).sqrt();
        result.coord_stats.insert(symbol.clone(), CoordStat { mean, stddev, count: counts.len() as i32 });
    }
    result.valid = true;
    result.message = "ADF computed successfully.".into();
    result
}
