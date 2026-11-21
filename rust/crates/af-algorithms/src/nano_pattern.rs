//! Nanowires and core-shell element assignment (the C++
//! algorithms/NanostructureTools) and neighbour-direction pattern matching
//! (algorithms/PatternMatch).

use af_core::cxx::sort_by;
use af_core::eigen;
use af_core::elements;
use af_core::{AtomSite, Structure};
use glam::DVec3;

const PI: f64 = 3.14159265358979323846;

fn normalize(v: DVec3) -> DVec3 {
    af_core::glm::normalize_d(v)
}

fn length(v: DVec3) -> f64 {
    v.dot(v).sqrt()
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NanowireParams {
    /// The periodic lattice vector (0 = a, 1 = b, 2 = c).
    pub axis: i32,
    /// Apothem (sides >= 3) or radius (A).
    pub radius: f64,
    /// Facets of the cross-section; < 3 is circular.
    pub sides: i32,
    pub vacuum: f64,
    pub axis_repeats: i32,
}

impl Default for NanowireParams {
    fn default() -> Self {
        Self { axis: 2, radius: 10.0, sides: 0, vacuum: 10.0, axis_repeats: 1 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NanowireResult {
    pub success: bool,
    pub message: String,
    pub structure: Structure,
    pub atom_count: i32,
}

pub fn build_nanowire(source: &Structure, params: &NanowireParams) -> NanowireResult {
    let mut result = NanowireResult::default();
    if !source.has_unit_cell || source.atoms.is_empty() {
        result.message = "Structure needs a unit cell and at least one atom.".into();
        return result;
    }
    if !(0..=2).contains(&params.axis) {
        result.message = "axis must be 0 (a), 1 (b), or 2 (c).".into();
        return result;
    }
    if params.radius <= 0.0 {
        result.message = "radius must be positive.".into();
        return result;
    }
    if params.axis_repeats < 1 {
        result.message = "axisRepeats must be at least 1.".into();
        return result;
    }
    let axis = params.axis as usize;
    let in_plane = [(axis + 1) % 3, (axis + 2) % 3];
    let cell: [DVec3; 3] = std::array::from_fn(|r| DVec3::from(source.cell_vectors[r]));
    let axis_len = length(cell[axis]);
    if axis_len < 1e-9 {
        result.message = "The wire-axis cell vector is degenerate.".into();
        return result;
    }
    let axis_dir = cell[axis] / axis_len;
    // An orthonormal frame of the cross-section plane.
    let e1 = normalize(if axis_dir.dot(DVec3::X).abs() < 0.9 { axis_dir.cross(DVec3::X) } else { axis_dir.cross(DVec3::Y) });
    let e2 = normalize(axis_dir.cross(e1));
    let coverage = params.radius + params.vacuum + axis_len;
    let reps = |v: DVec3| (coverage / length(v).max(1e-6)).ceil() as i32 + 1;
    let (reps1, reps2) = (reps(cell[in_plane[0]]), reps(cell[in_plane[1]]));
    let mut centroid = DVec3::ZERO;
    for a in &source.atoms {
        centroid += DVec3::new(a.x, a.y, a.z);
    }
    centroid /= source.atoms.len() as f64;
    let face_angle_step = if params.sides >= 3 { 2.0 * PI / params.sides as f64 } else { 0.0 };
    let mut kept: Vec<AtomSite> = Vec::new();
    for i1 in -reps1..=reps1 {
        for i2 in -reps2..=reps2 {
            for ia in 0..params.axis_repeats {
                let shift = i1 as f64 * cell[in_plane[0]] + i2 as f64 * cell[in_plane[1]] + ia as f64 * cell[axis];
                for atom in &source.atoms {
                    let rel = DVec3::new(atom.x, atom.y, atom.z) + shift - centroid;
                    let along = rel.dot(axis_dir);
                    if along < -1e-6 || along >= axis_len * params.axis_repeats as f64 - 1e-6 {
                        continue;
                    }
                    let (c1, c2) = (rel.dot(e1), rel.dot(e2));
                    let inside = if params.sides >= 3 {
                        (0..params.sides).all(|k| {
                            let angle = k as f64 * face_angle_step;
                            !(c1 * angle.cos() + c2 * angle.sin() > params.radius)
                        })
                    } else {
                        (c1 * c1 + c2 * c2).sqrt() <= params.radius
                    };
                    if !inside {
                        continue;
                    }
                    let mut placed = atom.clone();
                    let p = along * axis_dir + (c1 + params.radius + params.vacuum) * e1 + (c2 + params.radius + params.vacuum) * e2;
                    placed.x = p.x;
                    placed.y = p.y;
                    placed.z = p.z;
                    kept.push(placed);
                }
            }
        }
    }
    if kept.is_empty() {
        result.message = "No atoms fell inside the requested cross-section; increase radius.".into();
        return result;
    }
    let mut out = Structure { has_unit_cell: true, atoms: kept, ..Default::default() };
    let box_size = 2.0 * (params.radius + params.vacuum);
    let mut new_cell = [DVec3::ZERO; 3];
    new_cell[in_plane[0]] = e1 * box_size;
    new_cell[in_plane[1]] = e2 * box_size;
    new_cell[axis] = axis_dir * axis_len * params.axis_repeats as f64;
    out.cell_vectors = new_cell.map(|v| v.to_array());
    result.atom_count = out.atoms.len() as i32;
    result.structure = out;
    result.success = true;
    let kind = if params.sides >= 3 { format!("{}-sided", params.sides) } else { "circular".into() };
    result.message = format!("Built a {kind} nanowire with {} atoms, periodic along axis {axis}.", result.atom_count);
    result
}

#[derive(Clone, Debug, PartialEq)]
pub struct CoreShellParams {
    pub core_radius: f64,
    pub core_element: String,
    pub shell_element: String,
    /// Centre on the atoms' centroid, else on `center`.
    pub use_centroid: bool,
    pub center: [f64; 3],
}

impl Default for CoreShellParams {
    fn default() -> Self {
        Self { core_radius: 5.0, core_element: String::new(), shell_element: String::new(), use_centroid: true, center: [0.0; 3] }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CoreShellResult {
    pub success: bool,
    pub message: String,
    pub structure: Structure,
    pub core_count: i32,
    pub shell_count: i32,
}

pub fn apply_core_shell(source: &Structure, params: &CoreShellParams) -> CoreShellResult {
    let mut result = CoreShellResult::default();
    if source.atoms.is_empty() {
        result.message = "Structure has no atoms.".into();
        return result;
    }
    if params.core_element.is_empty() || params.shell_element.is_empty() {
        result.message = "Both coreElement and shellElement must be specified.".into();
        return result;
    }
    let (core_z, shell_z) = (elements::atomic_number(&params.core_element), elements::atomic_number(&params.shell_element));
    if core_z == 0 || shell_z == 0 {
        result.message = format!("Unknown element symbol ('{}' or '{}').", params.core_element, params.shell_element);
        return result;
    }
    let mut center = DVec3::from(params.center);
    if params.use_centroid {
        center = DVec3::ZERO;
        for a in &source.atoms {
            center += DVec3::new(a.x, a.y, a.z);
        }
        center /= source.atoms.len() as f64;
    }
    let mut out = source.clone();
    for atom in out.atoms.iter_mut() {
        let is_core = length(DVec3::new(atom.x, atom.y, atom.z) - center) < params.core_radius;
        let z = if is_core { core_z } else { shell_z };
        atom.atomic_number = z;
        atom.symbol = elements::symbol(z).to_string();
        [atom.r, atom.g, atom.b] = elements::default_colour(z);
        if is_core {
            result.core_count += 1;
        } else {
            result.shell_count += 1;
        }
    }
    result.structure = out;
    result.success = true;
    result.message = format!(
        "Assigned {} core ({}) and {} shell ({}) atoms.",
        result.core_count, params.core_element, result.shell_count, params.shell_element
    );
    result
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CrystalPattern {
    pub success: bool,
    pub message: String,
    /// Unit vectors from the representative atom to its neighbours.
    pub directions: Vec<DVec3>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PatternMatchRow {
    pub index: i32,
    pub matched: bool,
    pub neighbor_count: i32,
    pub max_angle_deviation_deg: f64,
}

impl Default for PatternMatchRow {
    fn default() -> Self {
        Self { index: -1, matched: false, neighbor_count: 0, max_angle_deviation_deg: 0.0 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PatternDetectResult {
    pub success: bool,
    pub message: String,
    pub rows: Vec<PatternMatchRow>,
}

type Cells = Option<(eigen::Matrix3, eigen::Matrix3)>;

/// The cell (columns a, b, c) and its inverse, as Eigen computes them.
fn cell_matrices(s: &Structure) -> Cells {
    if !s.has_unit_cell {
        return None;
    }
    let mut cell = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cell[j][i] = s.cell_vectors[i][j];
        }
    }
    if eigen::determinant3(&cell).abs() < 1e-12 {
        return None;
    }
    Some((cell, eigen::inverse3(&cell)))
}

fn neighbor_directions(pos: &[DVec3], i: usize, cutoff: f32, cells: &Cells) -> Vec<DVec3> {
    let cutoff2 = cutoff as f64 * cutoff as f64;
    let mut dirs = Vec::new();
    for (j, &pj) in pos.iter().enumerate() {
        if j == i {
            continue;
        }
        let mut delta = pj - pos[i];
        if let Some((cell, inverse)) = cells {
            let mut frac = eigen::mul3(inverse, delta.to_array());
            for f in frac.iter_mut() {
                *f -= f.round();
            }
            delta = DVec3::from(eigen::mul3(cell, frac));
        }
        let d2 = delta.dot(delta);
        if d2 > 1e-8 && d2 <= cutoff2 {
            dirs.push(delta / d2.sqrt());
        }
    }
    dirs
}

/// The neighbour directions of the atom nearest the centroid of `reference`.
pub fn build_pattern(reference: &Structure, cutoff: f32, use_pbc: bool) -> CrystalPattern {
    let mut result = CrystalPattern::default();
    let n = reference.atoms.len();
    if n < 2 {
        result.message = "Reference structure needs at least two atoms.".into();
        return result;
    }
    if cutoff <= 0.0 {
        result.message = "cutoff must be positive.".into();
        return result;
    }
    let pos: Vec<DVec3> = reference.atoms.iter().map(|a| DVec3::new(a.x, a.y, a.z)).collect();
    let mut center = DVec3::ZERO;
    for &p in &pos {
        center += p;
    }
    center /= n as f64;
    let cells = if use_pbc { cell_matrices(reference) } else { None };
    let mut representative = 0;
    let mut best = f64::MAX;
    for (i, &p) in pos.iter().enumerate() {
        let d2 = (p - center).dot(p - center);
        if d2 < best {
            best = d2;
            representative = i;
        }
    }
    result.directions = neighbor_directions(&pos, representative, cutoff, &cells);
    if result.directions.is_empty() {
        result.message = "No neighbors found within cutoff around the representative atom.".into();
        return result;
    }
    result.success = true;
    result.message = format!("Pattern built from atom {representative} with {} neighbor directions.", result.directions.len());
    result
}

/// Which atoms have the pattern's neighbour directions (greedy nearest-angle assignment).
pub fn detect_pattern(structure: &Structure, pattern: &CrystalPattern, cutoff: f32, angle_threshold_deg: f64, use_pbc: bool) -> PatternDetectResult {
    let mut result = PatternDetectResult::default();
    if !pattern.success || pattern.directions.is_empty() {
        result.message = "Pattern is not valid; build it with buildPattern first.".into();
        return result;
    }
    if cutoff <= 0.0 {
        result.message = "cutoff must be positive.".into();
        return result;
    }
    let n = structure.atoms.len();
    let pos: Vec<DVec3> = structure.atoms.iter().map(|a| DVec3::new(a.x, a.y, a.z)).collect();
    let cells = if use_pbc { cell_matrices(structure) } else { None };
    let threshold = angle_threshold_deg * PI / 180.0;
    for i in 0..n {
        let dirs = neighbor_directions(&pos, i, cutoff, &cells);
        let mut row = PatternMatchRow { index: i as i32, neighbor_count: dirs.len() as i32, ..Default::default() };
        if dirs.len() != pattern.directions.len() {
            result.rows.push(row);
            continue;
        }
        let mut candidates: Vec<(usize, usize, f64)> = Vec::with_capacity(dirs.len() * dirs.len());
        for (a, d) in dirs.iter().enumerate() {
            for (b, p) in pattern.directions.iter().enumerate() {
                let cos = d.dot(*p);
                let cos = if cos < -1.0 { -1.0 } else if 1.0 < cos { 1.0 } else { cos };
                candidates.push((a, b, cos.acos()));
            }
        }
        sort_by(&mut candidates, |x, y| x.2 < y.2);
        let mut used_a = vec![false; dirs.len()];
        let mut used_b = vec![false; pattern.directions.len()];
        let (mut assigned, mut max_angle, mut ok) = (0usize, 0.0f64, true);
        for &(a, b, angle) in &candidates {
            if used_a[a] || used_b[b] {
                continue;
            }
            if angle > threshold {
                ok = false;
                break;
            }
            used_a[a] = true;
            used_b[b] = true;
            if max_angle < angle {
                max_angle = angle;
            }
            assigned += 1;
            if assigned == dirs.len() {
                break;
            }
        }
        row.matched = ok && assigned == dirs.len();
        row.max_angle_deviation_deg = if row.matched { max_angle * 180.0 / PI } else { 0.0 };
        result.rows.push(row);
    }
    result.success = true;
    let matched = result.rows.iter().filter(|r| r.matched).count();
    result.message = format!("{matched} of {n} atoms matched the pattern.");
    result
}
