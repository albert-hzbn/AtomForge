//! Chemical short-range order (the C++ algorithms/ShortRangeOrderAnalysis):
//! neighbour shells clustered from all periodic pair distances, and the
//! Warren-Cowley and Rao-Curtin (2022) parameters per shell and pair.

use std::collections::BTreeMap;

use af_core::neighbours;
use af_core::Structure;
use glam::DVec3;

/// Warren-Cowley alpha_ij = 1 - P(j|i) / c_j (Cowley, Phys. Rev. 138 (1965) A1384).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SroEntry {
    pub shell: i32,
    pub shell_distance: f64,
    pub center_element: String,
    pub neighbor_element: String,
    pub nij: i32,
    pub zi: i32,
    pub pij: f64,
    pub cj: f64,
    pub alpha: f64,
}

/// Rao-Curtin alpha_ij = 1 - P_ij / (2 c_i c_j) (Acta Mater. 226 (2022) 117621).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SroEntryRC {
    pub shell: i32,
    pub shell_distance: f64,
    pub elem_i: String,
    pub elem_j: String,
    pub n_ij: i32,
    pub n_ji: i32,
    pub z_shell: f64,
    pub p_ij: f64,
    pub ci: f64,
    pub cj: f64,
    pub alpha: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SroReport {
    pub n_atoms: i32,
    pub n_shells: i32,
    pub shell_distances: Vec<f64>,
    pub composition: BTreeMap<String, f64>,
    pub entries: Vec<SroEntry>,
    pub message: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SroReportRC {
    pub n_atoms: i32,
    pub n_shells: i32,
    pub shell_distances: Vec<f64>,
    pub composition: BTreeMap<String, f64>,
    pub shell_coordination: Vec<f64>,
    pub entries: Vec<SroEntryRC>,
    pub message: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SroAnalysisResult {
    pub success: bool,
    pub message: String,
    pub warren_cowley: SroReport,
    pub rao_curtin: SroReportRC,
}

/// An ordered neighbour pair: atom j lies in shell `shell` around atom i.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SroPair {
    pub i: usize,
    pub j: usize,
    pub shell: usize,
}

/// Every periodic pair within the cutoff with its distance.
fn sro_neighbours(s: &Structure, cutoff: f64) -> Result<Vec<(usize, usize, f64)>, String> {
    Ok(neighbours::within(s, cutoff, true).map_err(|e| e.to_string())?.iter().map(|p| (p.i, p.j, p.r)).collect())
}

/// Sorted distances clustered into shells: each shell averages the distances
/// within `tolerance` of its first one.
fn cluster_shells(mut dists: Vec<f64>, tolerance: f64, max_shells: usize) -> Vec<f64> {
    dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut shells = Vec::new();
    let mut i = 0;
    while i < dists.len() && shells.len() < max_shells {
        let mut sum = 0.0;
        let mut j = i;
        while j < dists.len() && dists[j] - dists[i] < tolerance {
            sum += dists[j];
            j += 1;
        }
        if j == i {
            break; // a non-positive tolerance (the C++ loops forever here)
        }
        shells.push(sum / (j - i) as f64);
        i = j;
    }
    shells
}

pub fn compute_composition(structure: &Structure) -> BTreeMap<String, f64> {
    let mut comp = BTreeMap::new();
    let total = structure.atoms.len();
    if total == 0 {
        return comp;
    }
    let mut counts: BTreeMap<String, i32> = BTreeMap::new();
    for a in &structure.atoms {
        *counts.entry(a.symbol.clone()).or_default() += 1;
    }
    for (symbol, count) in counts {
        comp.insert(symbol, count as f64 / total as f64);
    }
    comp
}

/// The first `max_shells` neighbour shell distances.
pub fn build_shell_distances(structure: &Structure, max_shells: i32, tolerance: f64) -> Vec<f64> {
    if structure.atoms.is_empty() || !structure.has_unit_cell || max_shells <= 0 {
        return Vec::new();
    }
    let v = |r: usize| DVec3::from(structure.cell_vectors[r]);
    let volume = v(0).dot(v(1).cross(v(2))).abs();
    if volume < 1e-9 {
        return Vec::new();
    }
    // Start near the second-neighbour distance of a dense packing and widen
    // until the requested shells lie wholly inside the search.
    let mut cutoff = 1.6 * (volume / structure.atoms.len() as f64).cbrt();
    let mut complete = Vec::new();
    for _ in 0..8 {
        let Ok(pairs) = sro_neighbours(structure, cutoff) else { break };
        let dists = pairs.iter().map(|p| p.2).collect();
        complete = cluster_shells(dists, tolerance, max_shells as usize + 1).into_iter().filter(|&shell| shell + tolerance < cutoff).collect();
        if complete.len() as i32 >= max_shells {
            break;
        }
        cutoff *= 1.5;
    }
    complete.truncate(max_shells as usize);
    complete
}

pub fn build_pairs(structure: &Structure, shell_distances: &[f64], tolerance: f64) -> Vec<SroPair> {
    let mut pairs = Vec::new();
    if structure.atoms.is_empty() || !structure.has_unit_cell || shell_distances.is_empty() {
        return pairs;
    }
    let mut max = shell_distances[0];
    for &d in &shell_distances[1..] {
        if max < d {
            max = d;
        }
    }
    // The shells were found within a larger successful search, so this one succeeds too.
    for (i, j, d) in sro_neighbours(structure, max + tolerance).unwrap_or_default() {
        if let Some(shell) = shell_distances.iter().position(|&s| (d - s).abs() < tolerance) {
            pairs.push(SroPair { i, j, shell });
        }
    }
    pairs
}

pub fn analyze_warren_cowley(structure: &Structure, n_shells: i32, tolerance: f64) -> SroReport {
    let mut report = SroReport { n_atoms: structure.atoms.len() as i32, ..Default::default() };
    if structure.atoms.is_empty() || !structure.has_unit_cell {
        report.message = "No structure or unit cell defined.".into();
        return report;
    }
    report.composition = compute_composition(structure);
    report.shell_distances = build_shell_distances(structure, n_shells, tolerance);
    report.n_shells = report.shell_distances.len() as i32;
    if report.n_shells == 0 {
        report.message = "Could not identify neighbor shells.".into();
        return report;
    }
    let mut pairs: BTreeMap<(usize, &str, &str), i32> = BTreeMap::new();
    let mut coord: BTreeMap<(usize, &str), i32> = BTreeMap::new();
    for p in build_pairs(structure, &report.shell_distances, tolerance) {
        let (si, sj) = (structure.atoms[p.i].symbol.as_str(), structure.atoms[p.j].symbol.as_str());
        *pairs.entry((p.shell, si, sj)).or_default() += 1;
        *coord.entry((p.shell, si)).or_default() += 1;
    }
    for s in 0..report.n_shells as usize {
        for (ei, _) in &report.composition {
            let zi = coord.get(&(s, ei.as_str())).copied().unwrap_or(0);
            if zi == 0 {
                continue;
            }
            for (ej, &cj) in &report.composition {
                let nij = pairs.get(&(s, ei.as_str(), ej.as_str())).copied().unwrap_or(0);
                let pij = nij as f64 / zi as f64;
                report.entries.push(SroEntry {
                    shell: s as i32,
                    shell_distance: report.shell_distances[s],
                    center_element: ei.clone(),
                    neighbor_element: ej.clone(),
                    nij,
                    zi,
                    pij,
                    cj,
                    alpha: if cj > 1e-12 { 1.0 - pij / cj } else { 0.0 },
                });
            }
        }
    }
    report.message = "Warren-Cowley SRO completed.".into();
    report
}

pub fn analyze_rao_curtin(structure: &Structure, n_shells: i32, tolerance: f64) -> SroReportRC {
    let mut report = SroReportRC { n_atoms: structure.atoms.len() as i32, ..Default::default() };
    if structure.atoms.is_empty() || !structure.has_unit_cell {
        report.message = "No structure or unit cell defined.".into();
        return report;
    }
    report.composition = compute_composition(structure);
    report.shell_distances = build_shell_distances(structure, n_shells, tolerance);
    report.n_shells = report.shell_distances.len() as i32;
    if report.n_shells == 0 {
        report.message = "Could not identify neighbor shells.".into();
        return report;
    }
    let n = report.n_atoms;
    let mut pairs: BTreeMap<(usize, &str, &str), i32> = BTreeMap::new();
    let mut total_per_shell = vec![0i32; report.n_shells as usize];
    for p in build_pairs(structure, &report.shell_distances, tolerance) {
        *pairs.entry((p.shell, structure.atoms[p.i].symbol.as_str(), structure.atoms[p.j].symbol.as_str())).or_default() += 1;
        total_per_shell[p.shell] += 1;
    }
    report.shell_coordination = total_per_shell.iter().map(|&t| if n > 0 { t as f64 / n as f64 } else { 0.0 }).collect();
    // Unique pairs (i <= j), shell first.
    for s in 0..report.n_shells as usize {
        let zs = report.shell_coordination[s];
        for (ei, &ci) in &report.composition {
            for (ej, &cj) in &report.composition {
                if ei > ej || zs < 1e-12 {
                    continue;
                }
                let get = |a: &str, b: &str| pairs.get(&(s, a, b)).copied().unwrap_or(0);
                let nij = get(ei, ej);
                let nji = if ei != ej { get(ej, ei) } else { nij };
                let pij = (nij + nji) as f64 / (n as f64 * zs);
                let denom = 2.0 * ci * cj;
                report.entries.push(SroEntryRC {
                    shell: s as i32,
                    shell_distance: report.shell_distances[s],
                    elem_i: ei.clone(),
                    elem_j: ej.clone(),
                    n_ij: nij,
                    n_ji: nji,
                    z_shell: zs,
                    p_ij: pij,
                    ci,
                    cj,
                    alpha: if denom > 1e-12 { 1.0 - pij / denom } else { 0.0 },
                });
            }
        }
    }
    report.message = "Rao-Curtin SRO (2022) completed.".into();
    report
}

pub fn analyze_sro(structure: &Structure, n_shells: i32, tolerance: f64) -> SroAnalysisResult {
    if structure.atoms.is_empty() || !structure.has_unit_cell {
        return SroAnalysisResult { message: "Structure is empty or has no unit cell.".into(), ..Default::default() };
    }
    let warren_cowley = analyze_warren_cowley(structure, n_shells, tolerance);
    let rao_curtin = analyze_rao_curtin(structure, n_shells, tolerance);
    let success = warren_cowley.n_shells > 0 && rao_curtin.n_shells > 0;
    SroAnalysisResult {
        success,
        message: if success { "SRO analysis complete (Warren-Cowley + Rao-Curtin).".into() } else { "SRO analysis failed.".into() },
        warren_cowley,
        rao_curtin,
    }
}
