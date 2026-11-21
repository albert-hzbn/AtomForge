//! Dislocation analysis from displacement fields: Vitek differential
//! displacement maps (the C++ algorithms/VitekMap), the Nye tensor
//! (algorithms/NyeTensor, Hartley & Mishin 2005 as in BABEL) and the line
//! position and Burgers vector fitted from it (algorithms/DislocationFit).

use af_core::cxx::sort_by;
use af_core::eigen::{self, Matrix3};
use af_core::Structure;
use af_science::core::neighbor_list;
use glam::{DVec2, DVec3};

fn normalize(v: DVec3) -> DVec3 {
    af_core::glm::normalize_d(v)
}

/// std::to_string(double): "%f".
fn f6(x: f64) -> String {
    format!("{x:.6}")
}

/// The cell (columns a, b, c) and its inverse as Eigen computes them.
fn cell_matrices(s: &Structure) -> Option<(Matrix3, Matrix3)> {
    if !s.has_unit_cell {
        return None;
    }
    let cell = eigen::cell_columns(&s.cell_vectors);
    if eigen::determinant3(&cell).abs() < 1e-12 {
        return None;
    }
    Some((cell, eigen::inverse3(&cell)))
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VitekPair {
    pub index_i: i32,
    pub index_j: i32,
    pub position_i: DVec3,
    /// Possibly a periodic image.
    pub position_j: DVec3,
    /// (u_j - u_i) . line, wrapped to [-b/2, b/2).
    pub screw_component: f64,
    /// (u_j - u_i) in the plane perpendicular to the line.
    pub edge_component: DVec2,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VitekResult {
    pub success: bool,
    pub message: String,
    pub pairs: Vec<VitekPair>,
}

/// The nearest lattice translation along a and b only (BABEL's dsu(1)/dsu(2)).
fn wrap_in_plane(v: DVec3, cells: &Option<(Matrix3, Matrix3)>) -> DVec3 {
    let Some((cell, inverse)) = cells else { return v };
    let frac = eigen::mul3(inverse, v.to_array());
    let (n0, n1) = (frac[0].round(), frac[1].round());
    if n0 == 0.0 && n1 == 0.0 {
        return v;
    }
    let shift: [f64; 3] = std::array::from_fn(|r| n0 * cell[r][0] + n1 * cell[r][1]);
    v - DVec3::from(shift)
}

pub fn compute_vitek_map(
    reference: &Structure,
    deformed: &Structure,
    line_direction: DVec3,
    burgers_magnitude: f64,
    mut neighbor_cutoff: f32,
    use_pbc_xy: bool,
) -> VitekResult {
    let mut result = VitekResult::default();
    let n = reference.atoms.len();
    if n != deformed.atoms.len() {
        result.message = "Reference and deformed structures must have the same atom count and ordering.".into();
        return result;
    }
    if n < 2 {
        result.message = "Need at least two atoms.".into();
        return result;
    }
    let line = normalize(line_direction);
    let e1 = if line.dot(DVec3::X).abs() < 0.9 { normalize(line.cross(DVec3::X)) } else { normalize(line.cross(DVec3::Y)) };
    let e2 = normalize(line.cross(e1));
    let pos0: Vec<DVec3> = reference.atoms.iter().map(|a| DVec3::new(a.x, a.y, a.z)).collect();
    let u: Vec<DVec3> = deformed.atoms.iter().zip(&pos0).map(|(a, p)| DVec3::new(a.x, a.y, a.z) - *p).collect();
    let ref_cells = if use_pbc_xy { cell_matrices(reference) } else { None };
    let def_cells = cell_matrices(deformed);
    // A first/second-shell midpoint cutoff (perpendicular to the line) from up to 8 atoms.
    if neighbor_cutoff <= 0.0 {
        let mut shells: Vec<f64> = Vec::new();
        for s in 0..n.min(8) {
            let mut distances: Vec<f64> = Vec::new();
            for k in 0..n {
                if k == s {
                    continue;
                }
                let mut delta = pos0[k] - pos0[s];
                delta -= line * delta.dot(line);
                let d = delta.dot(delta).sqrt();
                if d > 1e-6 {
                    distances.push(d);
                }
            }
            sort_by(&mut distances, |a, b| a < b);
            for d in distances {
                if shells.last().is_none_or(|&last| d > last * 1.05) {
                    shells.push(d);
                }
                if shells.len() >= 2 {
                    break;
                }
            }
        }
        if shells.len() < 2 {
            result.message = "Could not auto-detect a neighbor cutoff; pass neighborCutoff explicitly.".into();
            return result;
        }
        sort_by(&mut shells, |a, b| a < b);
        neighbor_cutoff = (0.5 * (shells[0] + shells[1])) as f32;
    }
    let cutoff2 = neighbor_cutoff as f64 * neighbor_cutoff as f64;
    for i in 0..n {
        let ui = wrap_in_plane(u[i], &def_cells);
        for j in i + 1..n {
            let mut pos_j = pos0[j];
            let mut dr0 = pos_j - pos0[i];
            if let Some((cell, inverse)) = &ref_cells {
                let frac = eigen::mul3(inverse, dr0.to_array());
                let (n0, n1) = (frac[0].round(), frac[1].round());
                if n0 != 0.0 || n1 != 0.0 {
                    let shift = -DVec3::from(std::array::from_fn(|r| n0 * cell[r][0] + n1 * cell[r][1]));
                    dr0 += shift;
                    pos_j += shift;
                }
            }
            let perp = dr0 - line * dr0.dot(line);
            let r2 = perp.dot(perp);
            if r2 > cutoff2 || r2 < 1e-8 {
                continue;
            }
            let uj = wrap_in_plane(u[j], &def_cells);
            let du = uj - ui;
            let mut screw = du.dot(line);
            if burgers_magnitude > 0.0 {
                screw -= burgers_magnitude * (screw / burgers_magnitude).round();
            }
            result.pairs.push(VitekPair {
                index_i: i as i32,
                index_j: j as i32,
                position_i: pos0[i],
                position_j: pos_j,
                screw_component: screw,
                edge_component: DVec2::new(du.dot(e1), du.dot(e2)),
            });
        }
    }
    result.success = true;
    result.message = format!("Found {} neighbor pairs within cutoff {} A.", result.pairs.len(), f6(neighbor_cutoff as f64));
    result
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NyeTensorResult {
    pub success: bool,
    pub message: String,
    /// The 3x3 Nye tensor per atom (A^-1).
    pub alpha: Vec<Matrix3>,
    pub alpha_norm: Vec<f64>,
}

/// The minimum image of a displacement (columns of `cell` are the lattice vectors).
fn minimum_image(delta: DVec3, cells: &Option<(Matrix3, Matrix3)>) -> DVec3 {
    let Some((cell, inverse)) = cells else { return delta };
    let mut frac = eigen::mul3(inverse, delta.to_array());
    for f in frac.iter_mut() {
        *f -= f.round();
    }
    DVec3::from(eigen::mul3(cell, frac))
}

fn degenerate(b: &Matrix3) -> bool {
    eigen::determinant3(b).abs() < 1e-12 * eigen::norm3(b).powi(3) + 1e-300
}

pub fn compute_nye_tensor(reference: &Structure, deformed: &Structure, cutoff_radius: f32, use_pbc: bool) -> NyeTensorResult {
    let mut result = NyeTensorResult::default();
    let n = reference.atoms.len();
    if n != deformed.atoms.len() {
        result.message = "Reference and deformed structures must have the same atom count and ordering.".into();
        return result;
    }
    result.alpha = vec![[[0.0; 3]; 3]; n];
    result.alpha_norm = vec![0.0; n];
    if n == 0 {
        result.success = true;
        return result;
    }
    if cutoff_radius <= 0.0 {
        result.message = "cutoffRadius must be positive.".into();
        return result;
    }
    let use_pbc = use_pbc && reference.has_unit_cell;
    let mut cells = None;
    let mut cell_delta: Matrix3 = [[0.0; 3]; 3];
    if use_pbc {
        let cell0 = eigen::cell_columns(&reference.cell_vectors);
        if eigen::determinant3(&cell0).abs() < 1e-12 {
            result.message = "Reference unit cell is singular.".into();
            return result;
        }
        if deformed.has_unit_cell {
            cell_delta = eigen::sub3(&eigen::cell_columns(&deformed.cell_vectors), &cell0);
        }
        cells = Some((cell0, eigen::inverse3(&cell0)));
    }
    let pos0: Vec<DVec3> = reference.atoms.iter().map(|a| DVec3::new(a.x, a.y, a.z)).collect();
    let u: Vec<DVec3> = deformed.atoms.iter().zip(&pos0).map(|(a, p)| DVec3::new(a.x, a.y, a.z) - *p).collect();
    // Each neighbour once through its minimum image.
    let cutoff2 = cutoff_radius as f64 * cutoff_radius as f64;
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    {
        let points: Vec<[f64; 3]> = pos0.iter().map(|p| p.to_array()).collect();
        let cell = if use_pbc { reference.cell_vectors } else { af_science::core::identity() };
        let pairs = neighbor_list(&points, &cell, [use_pbc; 3], (cutoff_radius as f64).next_up()).unwrap_or_default();
        for pair in pairs {
            if pair.i == pair.j {
                continue;
            }
            let delta = minimum_image(pos0[pair.j] - pos0[pair.i], &cells);
            if delta.dot(delta) > cutoff2 {
                continue;
            }
            let list = &mut neighbors[pair.i];
            if !list.contains(&pair.j) {
                list.push(pair.j);
            }
        }
        for list in neighbors.iter_mut() {
            list.sort();
        }
    }
    // Step 1: displacement gradients du(a, b) = dU_a / dX_b by least squares.
    let mut du = vec![[[0.0f64; 3]; 3]; n];
    let mut defined = vec![false; n];
    let mut skipped = 0;
    for i in 0..n {
        let nbrs = &neighbors[i];
        if nbrs.len() < 4 {
            skipped += 1;
            continue;
        }
        let (mut a, mut b) = ([[0.0f64; 3]; 3], [[0.0f64; 3]; 3]);
        for &j in nbrs {
            let r0 = minimum_image(pos0[j] - pos0[i], &cells);
            let mut du_n = u[j] - u[i];
            if let Some((_, inverse)) = &cells {
                // The integer image count times the change of the cell (BABEL's ds0 = aNInt(ds0)).
                let ds0 = eigen::mul3(inverse, (pos0[j] - pos0[i]).to_array()).map(|x| x.round());
                du_n -= DVec3::from(eigen::mul3(&cell_delta, ds0));
            }
            a = eigen::add3(&a, &eigen::outer3(r0.to_array(), du_n.to_array()));
            b = eigen::add3(&b, &eigen::outer3(r0.to_array(), r0.to_array()));
        }
        if degenerate(&b) {
            skipped += 1;
            continue;
        }
        let f_total = eigen::matmul3(&a, &eigen::inverse3(&b));
        du[i] = eigen::transpose3(&f_total);
        defined[i] = true;
    }
    // Step 2: the lattice correspondence G = -[(I + du)^-1 du]^T.
    let mut g = vec![[[0.0f64; 3]; 3]; n];
    for i in 0..n {
        if defined[i] {
            let f_inv = eigen::inverse3(&eigen::add3(&eigen::IDENTITY3, &du[i]));
            g[i] = eigen::transpose3(&eigen::matmul3(&f_inv, &du[i])).map(|row| row.map(|x| -x));
        }
    }
    // Step 3: the gradient of G and its curl (with Olmsted's index order).
    for i in 0..n {
        let nbrs = &neighbors[i];
        if !defined[i] || nbrs.len() < 4 {
            continue;
        }
        let mut b = [[0.0f64; 3]; 3];
        let mut r0s = Vec::new();
        let mut dgs = Vec::new();
        let mut undefined_neighbor = false;
        for &j in nbrs {
            if !defined[j] {
                undefined_neighbor = true;
                continue;
            }
            let r0 = minimum_image(pos0[j] - pos0[i], &cells).to_array();
            r0s.push(r0);
            dgs.push(eigen::sub3(&g[j], &g[i]));
            b = eigen::add3(&b, &eigen::outer3(r0, r0));
        }
        if r0s.len() < 4 || degenerate(&b) {
            if undefined_neighbor || r0s.len() < 4 {
                skipped += 1;
            }
            continue;
        }
        let inv_b = eigen::inverse3(&b);
        let mut alpha = [[0.0f64; 3]; 3];
        for l in 0..3 {
            let mut v = [[0.0f64; 3]; 3];
            for (r0, dg) in r0s.iter().zip(&dgs) {
                v = eigen::add3(&v, &eigen::outer3(*r0, [dg[0][l], dg[1][l], dg[2][l]]));
            }
            let d = eigen::matmul3(&inv_b, &v); // d(a, b) = dG(b, l) / dx(a)
            alpha[l][0] = -d[1][2] + d[2][1];
            alpha[l][1] = -d[2][0] + d[0][2];
            alpha[l][2] = -d[0][1] + d[1][0];
        }
        result.alpha[i] = alpha;
        let mut norm_sq = 0.0;
        for row in &alpha {
            for x in row {
                norm_sq += x * x;
            }
        }
        result.alpha_norm[i] = norm_sq.sqrt();
    }
    result.success = true;
    result.message = if skipped > 0 {
        format!(
            "{skipped} of {n} atoms had too few/degenerate neighbors for a defined Nye tensor (left at zero); increase cutoffRadius if this is unexpected."
        )
    } else {
        "Nye tensor computed for all atoms.".into()
    };
    result
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DislocationFitResult {
    pub success: bool,
    pub message: String,
    /// The |alpha|-weighted centroid of the positions.
    pub line_position: DVec3,
    /// Nye's theorem: alpha contracted with the line, times the area per atom.
    pub burgers_vector: DVec3,
    pub total_weight: f64,
}

pub fn fit_dislocation_from_nye(positions: &[DVec3], nye: &NyeTensorResult, line_direction: DVec3, area_per_atom: f64) -> DislocationFitResult {
    let mut result = DislocationFitResult::default();
    if positions.len() != nye.alpha.len() || positions.is_empty() {
        result.message = "positions and the Nye tensor result must be the same non-empty size.".into();
        return result;
    }
    let line = normalize(line_direction);
    let mut total = 0.0;
    let mut centroid = DVec3::ZERO;
    let mut burgers = DVec3::ZERO;
    for (i, &p) in positions.iter().enumerate() {
        let w = nye.alpha_norm[i];
        total += w;
        centroid += w * p;
        let a = &nye.alpha[i];
        for row in 0..3 {
            let mut contracted = 0.0;
            for col in 0..3 {
                contracted += a[col][row] * line[col];
            }
            burgers[row] += contracted * area_per_atom;
        }
    }
    if total <= 1e-12 {
        result.message = "Nye tensor is zero everywhere; no dislocation signal to fit.".into();
        return result;
    }
    result.line_position = centroid / total;
    result.burgers_vector = burgers;
    result.total_weight = total;
    result.success = true;
    result.message = format!("Fitted from {} atoms, total |alpha| weight = {}.", positions.len(), f6(total));
    result
}
