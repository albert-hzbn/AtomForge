//! Homogeneous strain of a cell (the C++ algorithms/StrainTool) and random
//! vacancies (algorithms/VacancyBuilder).

use af_core::cxx_random::{shuffle, Mt19937};
use af_core::Structure;
use af_science::core::{norm, MinimumImage};

/// Deformation gradient F: new lattice vector i = F * old vector i (row-major).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StrainParams {
    pub f: [[f64; 3]; 3],
}

impl Default for StrainParams {
    fn default() -> Self {
        Self { f: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StrainResult {
    pub success: bool,
    pub message: String,
    pub structure: Structure,
}

pub fn engineering_strain(exx: f64, eyy: f64, ezz: f64, exy: f64, exz: f64, eyz: f64) -> StrainParams {
    StrainParams { f: [[1.0 + exx, exy, exz], [exy, 1.0 + eyy, eyz], [exz, eyz, 1.0 + ezz]] }
}

pub fn apply_strain(source: &Structure, params: &StrainParams) -> StrainResult {
    let mut result = StrainResult::default();
    if !source.has_unit_cell {
        result.message = "Structure has no unit cell to strain.".into();
        return result;
    }
    let c = &source.cell_vectors;
    let det = c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0])
        + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0]);
    if det.abs() < 1e-12 {
        result.message = "Structure's cell is degenerate (zero volume).".into();
        return result;
    }
    // Inverse of the row-vector cell: frac = cart * inverse.
    let k = 1.0 / det;
    let inv = [
        [(c[1][1] * c[2][2] - c[1][2] * c[2][1]) * k, (c[0][2] * c[2][1] - c[0][1] * c[2][2]) * k, (c[0][1] * c[1][2] - c[0][2] * c[1][1]) * k],
        [(c[1][2] * c[2][0] - c[1][0] * c[2][2]) * k, (c[0][0] * c[2][2] - c[0][2] * c[2][0]) * k, (c[0][2] * c[1][0] - c[0][0] * c[1][2]) * k],
        [(c[1][0] * c[2][1] - c[1][1] * c[2][0]) * k, (c[0][1] * c[2][0] - c[0][0] * c[2][1]) * k, (c[0][0] * c[1][1] - c[0][1] * c[1][0]) * k],
    ];
    let mut out = source.clone();
    for i in 0..3 {
        for j in 0..3 {
            let mut v = 0.0;
            for kk in 0..3 {
                v += params.f[i][kk] * c[kk][j];
            }
            out.cell_vectors[i][j] = v;
        }
    }
    let new_cell = out.cell_vectors;
    for atom in out.atoms.iter_mut() {
        let cart = [atom.x, atom.y, atom.z];
        let mut frac = [0.0; 3];
        for j in 0..3 {
            for i in 0..3 {
                frac[j] += cart[i] * inv[i][j];
            }
        }
        let mut moved = [0.0; 3];
        for j in 0..3 {
            for i in 0..3 {
                moved[j] += frac[i] * new_cell[i][j];
            }
        }
        atom.x = moved[0];
        atom.y = moved[1];
        atom.z = moved[2];
    }
    result.message = format!("Applied deformation gradient to cell and {} atoms.", out.atoms.len());
    result.structure = out;
    result.success = true;
    result
}

#[derive(Clone, Debug, PartialEq)]
pub struct VacancyParams {
    /// Only atoms of this element (empty: all).
    pub element: String,
    /// Percent of the eligible atoms (when target_count is 0).
    pub target_percentage: f64,
    /// Exact count; wins when positive.
    pub target_count: i32,
    /// Minimum separation (A) between vacancies; 0 disables.
    pub min_separation: f64,
    pub seed: u32,
}

impl Default for VacancyParams {
    fn default() -> Self {
        Self { element: String::new(), target_percentage: 0.0, target_count: 0, min_separation: 0.0, seed: 1 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VacancyResult {
    pub success: bool,
    pub message: String,
    pub structure: Structure,
    pub eligible_count: i32,
    pub requested_count: i32,
    pub removed_count: i32,
}

pub fn build_vacancies(source: &Structure, params: &VacancyParams) -> VacancyResult {
    let mut result = VacancyResult { structure: source.clone(), ..Default::default() };
    let mut eligible: Vec<usize> = (0..source.atoms.len()).filter(|&i| params.element.is_empty() || source.atoms[i].symbol == params.element).collect();
    result.eligible_count = eligible.len() as i32;
    if eligible.is_empty() {
        result.message = if params.element.is_empty() { "Structure has no atoms.".into() } else { format!("No atoms of element '{}' found.", params.element) };
        return result;
    }
    let requested = if params.target_count > 0 {
        params.target_count
    } else {
        (params.target_percentage / 100.0 * eligible.len() as f64).round() as i64 as i32
    };
    result.requested_count = requested;
    if requested < 0 || requested as usize > eligible.len() {
        result.message = format!("Requested {requested} vacancies but only {} atoms are eligible.", eligible.len());
        return result;
    }
    if requested == 0 {
        result.success = true;
        result.message = "Requested zero vacancies; structure unchanged.".into();
        return result;
    }
    let mut rng = Mt19937::new(params.seed);
    shuffle(&mut eligible, &mut rng);
    let mic = if source.has_unit_cell { MinimumImage::new(&source.cell_vectors, [true; 3]).ok() } else { None };
    let mut removed: Vec<usize> = Vec::with_capacity(requested as usize);
    for &idx in &eligible {
        if removed.len() as i32 >= requested {
            break;
        }
        if params.min_separation > 0.0 {
            let p = &source.atoms[idx];
            let too_close = removed.iter().any(|&chosen| {
                let q = &source.atoms[chosen];
                let d = [p.x - q.x, p.y - q.y, p.z - q.z];
                // Periodic images count: the minimum image.
                let separation = match &mic {
                    Some(m) => norm(m.apply(d)),
                    None => (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt(),
                };
                separation < params.min_separation
            });
            if too_close {
                continue;
            }
        }
        removed.push(idx);
    }
    let mut remove = vec![false; source.atoms.len()];
    for &i in &removed {
        remove[i] = true;
    }
    let n = source.atoms.len();
    let mut out = source.clone();
    out.atoms.clear();
    out.grain_colors.clear();
    out.grain_region_ids.clear();
    for i in 0..n {
        if remove[i] {
            continue;
        }
        out.atoms.push(source.atoms[i].clone());
        if source.grain_colors.len() == n {
            out.grain_colors.push(source.grain_colors[i]);
        }
        if source.grain_region_ids.len() == n {
            out.grain_region_ids.push(source.grain_region_ids[i]);
        }
    }
    result.structure = out;
    result.removed_count = removed.len() as i32;
    result.success = true;
    result.message = format!(
        "Removed {} of {requested} requested vacancies ({} eligible sites).{}",
        result.removed_count,
        eligible.len(),
        if result.removed_count < requested { " Minimum separation prevented reaching the full request." } else { "" }
    );
    result
}
