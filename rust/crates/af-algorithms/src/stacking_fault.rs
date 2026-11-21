//! Generalised stacking-fault sequences (the C++ algorithms/
//! StackingFaultBuilder): the crystal family from spglib (all atoms one
//! species) or common-neighbour signatures, the fault plane and slip
//! direction, and copies with the upper layers shifted step by step.

use std::collections::{HashMap, HashSet, VecDeque};

use af_core::cell::{BOND_TOLERANCE_FACTOR, MIN_BOND_DISTANCE};
use af_core::cxx::sort_by;
use af_core::elements;
use af_core::glm::{self, cross, dot, length, mul, normalize};
use af_core::supercell::build_supercell;
use af_core::Structure;
use glam::{IVec3, Mat3, Vec3};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StackingFaultFamily {
    #[default]
    Unknown = 0,
    Fcc,
    Hcp,
    Bcc,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StackingFaultPlane {
    #[default]
    Auto = 0,
    Fcc111,
    HcpBasal,
    HcpPrismatic,
    HcpPyramidal,
    Bcc110,
    Bcc112,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StackingFaultCellMode {
    #[default]
    SmallestUnitCell = 0,
    OrthogonalCell,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StackingFaultSequenceItem {
    pub structure: Structure,
    pub label: String,
    pub displacement_factor: f32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StackingFaultDetectionResult {
    pub success: bool,
    pub message: String,
    pub family: StackingFaultFamily,
    pub fcc_count: i32,
    pub hcp_count: i32,
    pub bcc_count: i32,
    pub recognized_count: i32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StackingFaultParams {
    pub plane: StackingFaultPlane,
    pub cell_mode: StackingFaultCellMode,
    pub layer_count: i32,
    pub interval: f32,
    pub max_displacement_factor: f32,
    pub use_pbc_for_detection: bool,
}

impl Default for StackingFaultParams {
    fn default() -> Self {
        Self {
            plane: StackingFaultPlane::Auto,
            cell_mode: StackingFaultCellMode::SmallestUnitCell,
            layer_count: 9,
            interval: 0.1,
            max_displacement_factor: 2.0,
            use_pbc_for_detection: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StackingFaultResult {
    pub success: bool,
    pub message: String,
    pub family: StackingFaultFamily,
    pub detected_phase: String,
    pub sequence: Vec<StackingFaultSequenceItem>,
    pub detected_layer_count: i32,
    pub shifted_atom_count: i32,
    pub partial_displacement: f32,
}

fn min_f(a: f32, b: f32) -> f32 {
    if b < a {
        b
    } else {
        a
    }
}
fn max_f(a: f32, b: f32) -> f32 {
    if a < b {
        b
    } else {
        a
    }
}

pub fn stacking_fault_family_name(family: StackingFaultFamily) -> &'static str {
    match family {
        StackingFaultFamily::Fcc => "FCC-like",
        StackingFaultFamily::Hcp => "HCP-like",
        StackingFaultFamily::Bcc => "BCC-like",
        StackingFaultFamily::Unknown => "Unknown",
    }
}

fn edge_key(a: usize, b: usize) -> u64 {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    ((a as u32 as u64) << 32) | b as u32 as u64
}

/// Bonds (not atoms) in the longest chain among the common neighbours.
fn longest_chain_length(common: &[usize], edges: &HashSet<u64>) -> i32 {
    if common.len() < 2 {
        return 0;
    }
    let n = common.len();
    let mut adjacency = vec![Vec::new(); n];
    for i in 0..n {
        for j in i + 1..n {
            if edges.contains(&edge_key(common[i], common[j])) {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
        }
    }
    let mut best = 0;
    for src in 0..n {
        let mut dist = vec![-1i32; n];
        let mut queue = VecDeque::new();
        dist[src] = 0;
        queue.push_back(src);
        while let Some(u) = queue.pop_front() {
            for &v in &adjacency[u] {
                if dist[v] >= 0 {
                    continue;
                }
                dist[v] = dist[u] + 1;
                best = best.max(dist[v]);
                queue.push_back(v);
            }
        }
    }
    best
}

#[derive(Default)]
struct DetectionSummary {
    family: StackingFaultFamily,
    fcc_count: i32,
    hcp_count: i32,
    bcc_count: i32,
    recognized_count: i32,
}

fn grid_cell(p: Vec3) -> (i32, i32, i32) {
    let c = |x: f32| (x / 4.0).floor() as i32;
    (c(p.x), c(p.y), c(p.z))
}

fn detect_family(structure: &Structure, use_pbc_request: bool) -> DetectionSummary {
    let mut summary = DetectionSummary::default();
    let n = structure.atoms.len();
    if n < 4 {
        return summary;
    }
    let radii = elements::literature_covalent_radii();
    let radius_of = |z: i32| if z >= 0 && (z as usize) < radii.len() { radii[z as usize] } else { 1.0 };
    let matrices = if use_pbc_request && structure.has_unit_cell { glm::cell_matrices(structure) } else { None };
    let use_pbc = matrices.is_some();
    let (cell, inv) = matrices.unwrap_or((Mat3::IDENTITY, Mat3::IDENTITY));
    let pairwise_pbc = use_pbc && n <= 5000;
    let positions: Vec<Vec3> = structure.atoms.iter().map(|a| Vec3::new(a.x as f32, a.y as f32, a.z as f32)).collect();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut edges: HashSet<u64> = HashSet::new();
    let mut try_add_bond = |i: usize, j: usize| {
        if j <= i {
            return;
        }
        let ri = radius_of(structure.atoms[i].atomic_number);
        let rj = radius_of(structure.atoms[j].atomic_number);
        let delta = glm::minimum_image_delta(positions[j] - positions[i], pairwise_pbc, &cell, &inv);
        let distance = length(delta);
        if distance <= MIN_BOND_DISTANCE {
            return;
        }
        if distance > (ri + rj) * BOND_TOLERANCE_FACTOR {
            return;
        }
        neighbors[i].push(j);
        neighbors[j].push(i);
        edges.insert(edge_key(i, j));
    };
    if !pairwise_pbc {
        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &p) in positions.iter().enumerate() {
            grid.entry(grid_cell(p)).or_default().push(i);
        }
        for (i, &p) in positions.iter().enumerate() {
            let c = grid_cell(p);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&(c.0 + dx, c.1 + dy, c.2 + dz)) {
                            for &j in list {
                                try_add_bond(i, j);
                            }
                        }
                    }
                }
            }
        }
    } else {
        // Candidate pairs over the periodic images, judged by the minimum image.
        let mut largest = 0.0f32;
        for a in &structure.atoms {
            if a.atomic_number >= 0 && (a.atomic_number as usize) < radii.len() {
                largest = max_f(largest, radii[a.atomic_number as usize]);
            }
        }
        largest = max_f(largest, 1.0);
        af_core::neighbours::for_each_candidate_pair(structure, 2.0 * largest as f64 * BOND_TOLERANCE_FACTOR as f64, true, |i, j| try_add_bond(i, j));
    }
    for list in neighbors.iter_mut() {
        list.sort();
    }
    // CNA signatures: FCC 4-2-1; HCP 4-2-1 and 4-2-2; BCC 6-6-3 and 4-4-2.
    let (mut n421, mut n422, mut n663, mut n442) = (vec![0; n], vec![0; n], vec![0; n], vec![0; n]);
    for i in 0..n {
        for &neighbor in &neighbors[i] {
            if neighbor <= i {
                continue;
            }
            // std::set_intersection of the sorted lists.
            let (a, b) = (&neighbors[i], &neighbors[neighbor]);
            let (mut x, mut y) = (0, 0);
            let mut common = Vec::new();
            while x < a.len() && y < b.len() {
                if a[x] < b[y] {
                    x += 1;
                } else if b[y] < a[x] {
                    y += 1;
                } else {
                    common.push(a[x]);
                    x += 1;
                    y += 1;
                }
            }
            let mut bonds = 0;
            for p in 0..common.len() {
                for q in p + 1..common.len() {
                    if edges.contains(&edge_key(common[p], common[q])) {
                        bonds += 1;
                    }
                }
            }
            let sig = (common.len() as i32, bonds, longest_chain_length(&common, &edges));
            let counter = match sig {
                (4, 2, 1) => Some(&mut n421),
                (4, 2, 2) => Some(&mut n422),
                (6, 6, 3) => Some(&mut n663),
                (4, 4, 2) => Some(&mut n442),
                _ => None,
            };
            if let Some(c) = counter {
                c[i] += 1;
                c[neighbor] += 1;
            }
        }
    }
    for i in 0..n {
        let close_packed = n421[i] + n422[i];
        let bcc = n663[i] + n442[i];
        if bcc > close_packed && bcc >= 4 {
            summary.bcc_count += 1;
        } else if close_packed >= 4 && n422[i] > 0 {
            summary.hcp_count += 1;
        } else if close_packed >= 4 {
            summary.fcc_count += 1;
        }
    }
    summary.recognized_count = summary.fcc_count + summary.hcp_count + summary.bcc_count;
    if summary.recognized_count <= 0 {
        return summary;
    }
    summary.family = if summary.fcc_count >= summary.hcp_count && summary.fcc_count >= summary.bcc_count {
        StackingFaultFamily::Fcc
    } else if summary.hcp_count >= summary.bcc_count {
        StackingFaultFamily::Hcp
    } else {
        StackingFaultFamily::Bcc
    };
    summary
}

fn direction_from_uvw(cell: &Mat3, uvw: IVec3) -> Vec3 {
    cell.x_axis * uvw.x as f32 + cell.y_axis * uvw.y as f32 + cell.z_axis * uvw.z as f32
}

fn plane_normal_from_hkl(cell: &Mat3, hkl: IVec3) -> Vec3 {
    let volume = dot(cell.x_axis, cross(cell.y_axis, cell.z_axis));
    if volume.abs() <= 1e-8 {
        return Vec3::new(0.0, 0.0, 1.0);
    }
    let a = cross(cell.y_axis, cell.z_axis) / volume;
    let b = cross(cell.z_axis, cell.x_axis) / volume;
    let c = cross(cell.x_axis, cell.y_axis) / volume;
    a * hkl.x as f32 + b * hkl.y as f32 + c * hkl.z as f32
}

fn resolve_plane_and_direction(family: StackingFaultFamily, plane: StackingFaultPlane) -> (IVec3, IVec3) {
    let v = IVec3::new;
    match family {
        StackingFaultFamily::Fcc => (v(1, 1, 2), v(1, 1, 1)),
        StackingFaultFamily::Hcp => match plane {
            StackingFaultPlane::HcpPrismatic => (v(1, 0, 0), v(0, 1, 0)),
            StackingFaultPlane::HcpPyramidal => (v(1, 0, 1), v(1, -1, 1)),
            _ => (v(0, 0, 1), v(1, -1, 0)),
        },
        StackingFaultFamily::Bcc => {
            if plane == StackingFaultPlane::Bcc112 {
                (v(1, 1, 2), v(1, 1, -1))
            } else {
                (v(1, 1, 0), v(1, -1, 1))
            }
        }
        StackingFaultFamily::Unknown => (v(0, 0, 1), v(1, 0, 0)),
    }
}

fn plane_label_for_selection(family: StackingFaultFamily, plane: StackingFaultPlane) -> &'static str {
    match family {
        StackingFaultFamily::Fcc => "(112) <111>",
        StackingFaultFamily::Bcc => {
            if plane == StackingFaultPlane::Bcc112 {
                "(112) <111>"
            } else {
                "(110) <111>"
            }
        }
        StackingFaultFamily::Hcp => match plane {
            StackingFaultPlane::HcpPrismatic => "prismatic",
            StackingFaultPlane::HcpPyramidal => "pyramidal",
            _ => "basal",
        },
        StackingFaultFamily::Unknown => "unknown",
    }
}

fn family_partial_scale(family: StackingFaultFamily) -> f32 {
    match family {
        StackingFaultFamily::Fcc => 1.0 / 6.0,
        StackingFaultFamily::Hcp => 1.0 / 3.0,
        StackingFaultFamily::Bcc => 0.5,
        StackingFaultFamily::Unknown => 0.25,
    }
}

/// Layer index of each projection (layers within 0.35 of the smallest gap merge).
fn assign_layer_ids(projections: &[f32]) -> (Vec<usize>, f32) {
    let mut sorted = projections.to_vec();
    sort_by(&mut sorted, |a, b| a < b);
    let mut min_diff = f32::MAX;
    for w in sorted.windows(2) {
        let diff = w[1] - w[0];
        if diff > 1e-4 {
            min_diff = min_f(min_diff, diff);
        }
    }
    if !min_diff.is_finite() {
        min_diff = 1.0;
    }
    let tolerance = max_f(1e-3, min_diff * 0.35);
    let mut centers = vec![sorted[0]];
    for &s in &sorted[1..] {
        if (s - centers[centers.len() - 1]).abs() > tolerance {
            centers.push(s);
        }
    }
    let spacing = if centers.len() >= 2 { max_f(1e-3, centers[1] - centers[0]) } else { min_diff };
    let ids = projections
        .iter()
        .map(|&p| {
            let mut best = 0;
            let mut best_dist = (p - centers[0]).abs();
            for (layer, &c) in centers.iter().enumerate().skip(1) {
                let dist = (p - c).abs();
                if dist < best_dist {
                    best_dist = dist;
                    best = layer;
                }
            }
            best
        })
        .collect();
    (ids, spacing)
}

fn plane_candidates(family: StackingFaultFamily, plane: StackingFaultPlane) -> Vec<IVec3> {
    let l = |v: &[[i32; 3]]| v.iter().map(|&[a, b, c]| IVec3::new(a, b, c)).collect::<Vec<_>>();
    const P112: [[i32; 3]; 9] = [[1, 1, 2], [1, 2, 1], [2, 1, 1], [-1, 1, 2], [1, -1, 2], [1, 1, -2], [-1, -1, 2], [-1, 1, -2], [1, -1, -2]];
    match family {
        StackingFaultFamily::Fcc => l(&P112),
        StackingFaultFamily::Bcc => {
            if plane == StackingFaultPlane::Bcc112 {
                l(&P112)
            } else {
                l(&[[1, 1, 0], [1, 0, 1], [0, 1, 1], [-1, 1, 0], [1, -1, 0], [-1, 0, 1], [1, 0, -1], [0, -1, 1], [0, 1, -1]])
            }
        }
        StackingFaultFamily::Hcp => match plane {
            StackingFaultPlane::HcpPrismatic => l(&[[1, 0, 0], [0, 1, 0], [-1, 1, 0], [-1, 0, 0], [0, -1, 0], [1, -1, 0]]),
            StackingFaultPlane::HcpPyramidal => l(&[[1, 0, 1], [0, 1, 1], [-1, 1, 1], [-1, 0, -1], [0, -1, -1], [1, -1, -1]]),
            _ => l(&[[0, 0, 1], [0, 0, -1]]),
        },
        StackingFaultFamily::Unknown => l(&[[0, 0, 1]]),
    }
}

fn direction_candidates(family: StackingFaultFamily, plane: StackingFaultPlane) -> Vec<IVec3> {
    let l = |v: &[[i32; 3]]| v.iter().map(|&[a, b, c]| IVec3::new(a, b, c)).collect::<Vec<_>>();
    const D111: [[i32; 3]; 8] = [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]];
    match family {
        StackingFaultFamily::Fcc | StackingFaultFamily::Bcc => l(&D111),
        StackingFaultFamily::Hcp => {
            if plane == StackingFaultPlane::HcpPyramidal {
                l(&[[1, 0, 1], [0, 1, 1], [-1, 1, 1], [-1, 0, -1], [0, -1, -1], [1, -1, -1]])
            } else {
                l(&[[1, 0, 0], [0, 1, 0], [1, -1, 0], [-1, 0, 0], [0, -1, 0], [-1, 1, 0]])
            }
        }
        StackingFaultFamily::Unknown => l(&[[1, 0, 0]]),
    }
}

fn projections(structure: &Structure, normal: Vec3) -> Vec<f32> {
    structure.atoms.iter().map(|a| dot(Vec3::new(a.x as f32, a.y as f32, a.z as f32), normal)).collect()
}

fn evaluate_plane(structure: &Structure, normal: Vec3) -> f32 {
    if structure.atoms.is_empty() || length(normal) <= 1e-6 {
        return -1e9;
    }
    let (ids, _) = assign_layer_ids(&projections(structure, normal));
    let layers = ids.iter().max().copied().unwrap_or(0) as i32 + 1;
    if layers < 2 {
        return -1e9;
    }
    let shift_start = (layers / 2) as usize;
    let shifted = ids.iter().filter(|&&id| id >= shift_start).count();
    let half = max_f(1.0, 0.5 * structure.atoms.len() as f32);
    let imbalance = (shifted as f32 - half).abs() / half;
    layers as f32 - 0.5 * imbalance
}

fn refine_plane_and_direction(
    structure: &Structure,
    cell: &Mat3,
    family: StackingFaultFamily,
    plane: StackingFaultPlane,
    hkl: &mut IVec3,
    uvw: &mut IVec3,
) {
    let mut best_plane = *hkl;
    let mut best_score = -1e9f32;
    let mut best_normal = Vec3::new(0.0, 0.0, 1.0);
    for candidate in plane_candidates(family, plane) {
        let normal = plane_normal_from_hkl(cell, candidate);
        if length(normal) <= 1e-6 {
            continue;
        }
        let normal = normalize(normal);
        let score = evaluate_plane(structure, normal);
        if score > best_score {
            best_score = score;
            best_plane = candidate;
            best_normal = normal;
        }
    }
    *hkl = best_plane;
    let mut best_direction = *uvw;
    let mut best_orthogonality = f32::MAX;
    let mut best_length = f32::MAX;
    for candidate in direction_candidates(family, plane) {
        let dir = direction_from_uvw(cell, candidate);
        let len = length(dir);
        if len <= 1e-6 {
            continue;
        }
        let orthogonality = dot(normalize(dir), best_normal).abs();
        if orthogonality < best_orthogonality - 1e-5 || ((orthogonality - best_orthogonality).abs() <= 1e-5 && len < best_length) {
            best_orthogonality = orthogonality;
            best_length = len;
            best_direction = candidate;
        }
    }
    *uvw = best_direction;
}

fn make_orthogonal_cell_structure(source: &Structure, slip: Vec3, plane_normal: Vec3) -> Structure {
    let mut output = source.clone();
    let e1 = normalize(slip);
    let e3 = normalize(plane_normal);
    let mut e2 = normalize(cross(e3, e1));
    if length(e2) <= 1e-6 {
        e2 = Vec3::new(0.0, 1.0, 0.0);
    }
    let mut min_c = [f32::MAX; 3];
    let mut max_c = [-f32::MAX; 3];
    let local: Vec<Vec3> = source
        .atoms
        .iter()
        .map(|a| {
            let p = Vec3::new(a.x as f32, a.y as f32, a.z as f32);
            let l = Vec3::new(dot(p, e1), dot(p, e2), dot(p, e3));
            for k in 0..3 {
                min_c[k] = min_f(min_c[k], l[k]);
                max_c[k] = max_f(max_c[k], l[k]);
            }
            l
        })
        .collect();
    let bx = e1 * max_f(1.0, max_c[0] - min_c[0]);
    let by = e2 * max_f(1.0, max_c[1] - min_c[1]);
    let bz = e3 * max_f(1.0, max_c[2] - min_c[2]);
    output.has_unit_cell = true;
    output.cell_vectors = [bx, by, bz].map(|v| [v.x as f64, v.y as f64, v.z as f64]);
    output.cell_offset = [0.0; 3];
    for (atom, l) in output.atoms.iter_mut().zip(&local) {
        let shifted = *l - Vec3::new(min_c[0], min_c[1], min_c[2]);
        let world = e1 * shifted.x + e2 * shifted.y + e3 * shifted.z;
        atom.x = world.x as f64;
        atom.y = world.y as f64;
        atom.z = world.z as f64;
    }
    output
}

fn build_sequence_label(family: StackingFaultFamily, factor: f32, interval: f32) -> String {
    let tol = max_f(0.02, interval * 0.5);
    let close_packed = family == StackingFaultFamily::Fcc || family == StackingFaultFamily::Hcp;
    if close_packed && (factor - 1.0).abs() <= tol {
        return "ISF".into();
    }
    if close_packed && (factor - 2.0).abs() <= tol {
        return "ESF".into();
    }
    format!("Shift {:.2}", factor)
}

/// FCC (F-centred cubic), BCC (I-centred cubic) or HCP (hexagonal) from
/// spglib with every atom one species.
fn classify_with_spglib(structure: &Structure) -> Result<(StackingFaultFamily, i32, String), String> {
    if !structure.has_unit_cell {
        return Err("No unit cell available for spglib detection.".into());
    }
    if structure.atoms.is_empty() {
        return Err("No atoms available for spglib detection.".into());
    }
    let Some((_, inv)) = glm::cell_matrices(structure) else {
        return Err("Invalid unit cell for spglib detection.".into());
    };
    let origin = Vec3::new(structure.cell_offset[0] as f32, structure.cell_offset[1] as f32, structure.cell_offset[2] as f32);
    let positions: Vec<[f64; 3]> = structure
        .atoms
        .iter()
        .map(|a| {
            let mut f = mul(&inv, Vec3::new(a.x as f32, a.y as f32, a.z as f32) - origin);
            f -= glm::floor(f);
            [f.x as f64, f.y as f64, f.z as f64]
        })
        .collect();
    let types = vec![1; positions.len()];
    let dataset = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
        .iter()
        .find_map(|&symprec| af_science::spg::get_dataset(&structure.cell_vectors, &positions, &types, symprec));
    let Some(dataset) = dataset else {
        return Err("spglib could not identify symmetry for this structure.".into());
    };
    let (number, symbol) = (dataset.spacegroup_number, dataset.international_symbol);
    let centering = symbol.chars().next().map(|c| c.to_ascii_uppercase());
    let family = if (195..=230).contains(&number) {
        match centering {
            Some('F') => StackingFaultFamily::Fcc,
            Some('I') => StackingFaultFamily::Bcc,
            _ => StackingFaultFamily::Unknown,
        }
    } else if (168..=194).contains(&number) {
        StackingFaultFamily::Hcp
    } else {
        StackingFaultFamily::Unknown
    };
    if family == StackingFaultFamily::Unknown {
        return Err(format!("spglib identified SG {number} ({symbol}), but it does not map to FCC/HCP/BCC."));
    }
    Ok((family, number, symbol))
}

pub fn detect_stacking_fault_family(structure: &Structure, use_pbc_for_detection: bool) -> StackingFaultDetectionResult {
    let mut result = StackingFaultDetectionResult::default();
    if structure.atoms.is_empty() {
        result.message = "Input structure has no atoms.".into();
        return result;
    }
    if structure.has_unit_cell {
        if let Ok((family, number, symbol)) = classify_with_spglib(structure) {
            result.success = true;
            result.family = family;
            result.recognized_count = structure.atoms.len() as i32;
            result.message = format!(
                "Detected {} using spglib (homogeneous species), SG {number} ({symbol}).",
                stacking_fault_family_name(family)
            );
            return result;
        }
    }
    let detection = detect_family(structure, use_pbc_for_detection);
    result.family = detection.family;
    result.fcc_count = detection.fcc_count;
    result.hcp_count = detection.hcp_count;
    result.bcc_count = detection.bcc_count;
    result.recognized_count = detection.recognized_count;
    if detection.family == StackingFaultFamily::Unknown {
        result.message = "Could not identify an FCC, HCP, or BCC-like environment.".into();
        return result;
    }
    result.message = format!(
        "Detected {} environment (FCC={}, HCP={}, BCC={}).",
        stacking_fault_family_name(detection.family),
        detection.fcc_count,
        detection.hcp_count,
        detection.bcc_count
    );
    result.success = true;
    result
}

/// The index of the cell vector most aligned with the normal.
fn normal_vector_index(cell: &Mat3, normal: Vec3) -> usize {
    let mut best = 0.0f32;
    let mut index = 0;
    for axis in 0..3 {
        let alignment = dot(normalize(cell.col(axis)), normal).abs();
        if alignment > best {
            best = alignment;
            index = axis;
        }
    }
    index
}

pub fn build_stacking_fault_sequence(base: &Structure, params: &StackingFaultParams) -> StackingFaultResult {
    let mut result = StackingFaultResult::default();
    if base.atoms.is_empty() {
        result.message = "Input structure has no atoms.".into();
        return result;
    }
    if !base.has_unit_cell {
        result.message = "Stacking faults builder requires a structure with a unit cell.".into();
        return result;
    }
    if params.layer_count < 2 || params.layer_count > 10000 {
        result.message = "Layer count must be at least 2.".into();
        return result;
    }
    if !params.interval.is_finite()
        || !params.max_displacement_factor.is_finite()
        || params.interval <= 0.0
        || params.max_displacement_factor < 0.0
        || params.max_displacement_factor / params.interval > 10000.0
    {
        result.message = "Interval and maximum displacement must be positive.".into();
        return result;
    }
    let detection = detect_stacking_fault_family(base, params.use_pbc_for_detection);
    result.family = detection.family;
    result.detected_phase = stacking_fault_family_name(detection.family).into();
    if detection.family == StackingFaultFamily::Unknown {
        result.message = "Could not identify an FCC, HCP, or BCC-like environment from the loaded structure.".into();
        return result;
    }
    let fcc = detection.family == StackingFaultFamily::Fcc;
    // FCC: a working cell with one axis along [111] and the shear axis along [11-2].
    let mut working = if fcc { build_supercell(base, &[[1, 1, -2], [1, -1, 0], [1, 1, 1]]) } else { base.clone() };
    let Some((mut cell, _)) = glm::cell_matrices(&working) else {
        result.message = "Input unit cell is singular.".into();
        return result;
    };
    let (mut hkl, mut uvw) = (IVec3::new(0, 0, 1), IVec3::new(1, 0, 0));
    let plane_label;
    let mut plane_normal;
    let mut raw_direction;
    let mut normal_index;
    if fcc {
        plane_normal = normalize(cell.z_axis);
        raw_direction = cell.x_axis;
        plane_label = "(111) <-1-12>";
        normal_index = 2;
    } else {
        (hkl, uvw) = resolve_plane_and_direction(detection.family, params.plane);
        refine_plane_and_direction(&working, &cell, detection.family, params.plane, &mut hkl, &mut uvw);
        plane_label = plane_label_for_selection(detection.family, params.plane);
        plane_normal = plane_normal_from_hkl(&cell, hkl);
        raw_direction = direction_from_uvw(&cell, uvw);
        if length(plane_normal) <= 1e-6 || length(raw_direction) <= 1e-6 {
            result.message = "Could not derive stacking-fault plane vectors from the current cell.".into();
            return result;
        }
        plane_normal = normalize(plane_normal);
        normal_index = normal_vector_index(&cell, plane_normal);
    }
    let mut slip = raw_direction - plane_normal * dot(raw_direction, plane_normal);
    if length(slip) <= 1e-6 {
        result.message = "Chosen slip direction is not compatible with the detected plane.".into();
        return result;
    }
    slip = normalize(slip);
    let (mut layer_ids, _) = assign_layer_ids(&projections(&working, plane_normal));
    let mut detected_layers = layer_ids.iter().max().copied().unwrap_or(0) as i32 + 1;
    if detected_layers >= 1 && detected_layers < params.layer_count {
        let repeats = 2.max((params.layer_count + detected_layers - 1) / detected_layers);
        let mut transform = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        transform[normal_index][normal_index] = repeats;
        working = build_supercell(&working, &transform);
        let Some((c, _)) = glm::cell_matrices(&working) else {
            result.message = "Could not build a valid supercell for stacking-fault generation.".into();
            return result;
        };
        cell = c;
        if fcc {
            plane_normal = normalize(cell.z_axis);
            raw_direction = cell.x_axis;
            normal_index = 2;
        } else {
            plane_normal = plane_normal_from_hkl(&cell, hkl);
            raw_direction = direction_from_uvw(&cell, uvw);
            if length(plane_normal) <= 1e-6 || length(raw_direction) <= 1e-6 {
                result.message = "Could not derive stacking-fault vectors after supercell expansion.".into();
                return result;
            }
            plane_normal = normalize(plane_normal);
            normal_index = normal_vector_index(&cell, plane_normal);
        }
        slip = raw_direction - plane_normal * dot(raw_direction, plane_normal);
        if length(slip) <= 1e-6 {
            result.message = "Slip direction became invalid after supercell expansion.".into();
            return result;
        }
        slip = normalize(slip);
        layer_ids = assign_layer_ids(&projections(&working, plane_normal)).0;
        detected_layers = layer_ids.iter().max().copied().unwrap_or(0) as i32 + 1;
    }
    result.detected_layer_count = detected_layers;
    if detected_layers < 2 {
        result.message = "The selected structure does not contain enough distinct layers along the fault normal.".into();
        return result;
    }
    let used = params.layer_count.min(detected_layers);
    let start = 0.max((detected_layers - used) / 2);
    let end = start + used;
    let candidates: Vec<usize> = (0..layer_ids.len()).filter(|&i| (layer_ids[i] as i32) >= start && (layer_ids[i] as i32) < end).collect();
    if candidates.is_empty() {
        result.message = "No atoms fell into the shifted layer window. Increase the number of layers or load a thicker slab.".into();
        return result;
    }
    let shift_start = start + 1.max(used / 2);
    let shifted_count = candidates.iter().filter(|&&i| layer_ids[i] as i32 >= shift_start).count() as i32;
    if shifted_count <= 0 {
        result.message = "No atoms available in the selected shift region.".into();
        return result;
    }
    result.partial_displacement = length(raw_direction) * family_partial_scale(detection.family);
    result.shifted_atom_count = shifted_count;
    let step_count = 1.max((params.max_displacement_factor / params.interval + 0.5).floor() as i32);
    for step in 0..=step_count {
        let factor = min_f(params.max_displacement_factor, step as f32 * params.interval);
        let displacement = slip * (result.partial_displacement * factor);
        let mut generated = working.clone();
        for &i in &candidates {
            if layer_ids[i] as i32 >= shift_start {
                let a = &mut generated.atoms[i];
                a.x += displacement.x as f64;
                a.y += displacement.y as f64;
                a.z += displacement.z as f64;
            }
        }
        for k in 0..3 {
            generated.cell_vectors[normal_index][k] += displacement[k] as f64;
        }
        if params.cell_mode == StackingFaultCellMode::OrthogonalCell {
            generated = make_orthogonal_cell_structure(&generated, slip, plane_normal);
        }
        result.sequence.push(StackingFaultSequenceItem {
            structure: generated,
            displacement_factor: factor,
            label: build_sequence_label(detection.family, factor, params.interval),
        });
    }
    result.message = format!(
        "Generated {} {} stacking-fault structures on {plane_label}. Shifted atoms: {}, detected layers: {detected_layers}.",
        result.sequence.len(),
        stacking_fault_family_name(detection.family),
        result.shifted_atom_count
    );
    result.success = true;
    result
}
