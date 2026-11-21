//! The cell sculptor (the C++ algorithms/CellSculptorAlgo): supercells cut by
//! (hkl) slabs, periodic slabs snapped to the structural repeat of their
//! planes, cells from three bounding slabs, and coincident-atom removal.

use std::collections::HashMap;

use af_core::cxx::sort_by;
use af_core::glm::{self, cross, dot, length, mul};
use af_core::Structure;
use glam::{DMat3, DVec3, Mat3, Vec3};

/// One slab: two parallel planes at [d1, d2] (A) along the (hkl) normal from
/// the reference point; periodic slabs use multiples of the structural period.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellSlabPlane {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    pub d1: f32,
    pub d2: f32,
    pub use_periodic: bool,
    pub start_plane: i32,
    pub n_periods: i32,
    pub color: [f32; 3],
}

impl Default for CellSlabPlane {
    fn default() -> Self {
        Self { h: 1, k: 0, l: 0, d1: -5.0, d2: 5.0, use_periodic: true, start_plane: 0, n_periods: 1, color: [0.4, 0.7, 1.0] }
    }
}

/// Tolerance (A) added to slab bounds so atoms exactly on a plane are kept.
pub const SLAB_TOL: f32 = 0.02;

fn lattice_matrix(source: &Structure) -> Mat3 {
    glm::cell_matrix(source)
}

/// The reciprocal-lattice matrix B = A^-T (B * hkl is the plane normal).
fn recip_matrix(source: &Structure) -> Mat3 {
    glm::inverse(&lattice_matrix(source)).transpose()
}

fn hkl(slab: &CellSlabPlane) -> Vec3 {
    Vec3::new(slab.h as f32, slab.k as f32, slab.l as f32)
}

/// The smallest N >= 1 such that N d_hkl along the normal maps every source
/// atom onto one of the same species (mod the cell).
fn structural_period_multiplier(slab: &CellSlabPlane, source: &Structure) -> i32 {
    if !source.has_unit_cell || source.atoms.is_empty() {
        return 1;
    }
    let g = mul(&recip_matrix(source), hkl(slab));
    let gsq = dot(g, g);
    if gsq < 1e-18 {
        return 1;
    }
    let step = g / gsq;
    let a_inv = glm::inverse(&lattice_matrix(source));
    let frac_atoms: Vec<(Vec3, &str)> = source
        .atoms
        .iter()
        .map(|a| {
            let mut f = mul(&a_inv, Vec3::new(a.x as f32, a.y as f32, a.z as f32));
            for i in 0..3 {
                f[i] -= f[i].floor();
            }
            (f, a.symbol.as_str())
        })
        .collect();
    let frac_step = mul(&a_inv, step);
    let tol = 0.01f32;
    for n in 1..=24 {
        let ft_step = frac_step * n as f32;
        let symmetric = frac_atoms.iter().all(|&(f, symbol)| {
            let mut ft = Vec3::ZERO;
            for i in 0..3 {
                ft[i] = f[i] + ft_step[i];
                ft[i] -= ft[i].floor();
            }
            frac_atoms.iter().any(|&(gf, gsym)| {
                if gsym != symbol {
                    return false;
                }
                let mut d = ft - gf;
                for i in 0..3 {
                    d[i] -= d[i].round();
                }
                d.x.abs() < tol && d.y.abs() < tol && d.z.abs() < tol
            })
        });
        if symmetric {
            return n;
        }
    }
    1
}

/// Distinct atom-plane projections (mod the period) over the 27 nearest cells.
fn all_layer_projections(slab: &CellSlabPlane, source: &Structure, period: f32) -> Vec<f32> {
    if period < 1e-9 || source.atoms.is_empty() {
        return Vec::new();
    }
    let n = csc_normal(slab, source);
    let tol = period * 0.005;
    let cv = &source.cell_vectors;
    let mut projections = Vec::with_capacity(source.atoms.len() * 27);
    for ia in -1..=1 {
        for ib in -1..=1 {
            for ic in -1..=1 {
                let (fa, fb, fc) = (ia as f64, ib as f64, ic as f64);
                let t = Vec3::new(
                    (cv[0][0] * fa + cv[1][0] * fb + cv[2][0] * fc) as f32,
                    (cv[0][1] * fa + cv[1][1] * fb + cv[2][1] * fc) as f32,
                    (cv[0][2] * fa + cv[1][2] * fb + cv[2][2] * fc) as f32,
                );
                for a in &source.atoms {
                    let p = Vec3::new(a.x as f32 + t.x, a.y as f32 + t.y, a.z as f32 + t.z);
                    let mut proj = dot(p, n) % period;
                    if proj < 0.0 {
                        proj += period;
                    }
                    projections.push(proj);
                }
            }
        }
    }
    sort_by(&mut projections, |a, b| a < b);
    // std::unique with "b - a < tol" against the last kept value.
    let mut unique: Vec<f32> = Vec::with_capacity(projections.len());
    for p in projections {
        match unique.last() {
            Some(&kept) if p - kept < tol => {}
            _ => unique.push(p),
        }
    }
    while unique.last().is_some_and(|&p| p >= period - tol) {
        unique.pop();
    }
    unique
}

/// The reference point of a slab: periodic slabs are phase-aligned so that
/// multiples of the structural period land on equivalent atom planes.
pub fn csc_slab_ref(slab: &CellSlabPlane, source: &Structure, centroid: Vec3) -> Vec3 {
    if !slab.use_periodic || !source.has_unit_cell {
        return centroid;
    }
    let base = Vec3::new(source.cell_offset[0] as f32, source.cell_offset[1] as f32, source.cell_offset[2] as f32);
    if source.atoms.is_empty() {
        return base;
    }
    let n = csc_normal(slab, source);
    let period = structural_period_multiplier(slab, source) as f32 * csc_dhkl(slab, source);
    let phase = all_layer_projections(slab, source, period).first().copied().unwrap_or(0.0);
    base - n * phase
}

/// The unit (hkl) normal (reciprocal lattice with a cell, else (h, k, l)).
pub fn csc_normal(slab: &CellSlabPlane, source: &Structure) -> Vec3 {
    let n = if source.has_unit_cell { mul(&recip_matrix(source), hkl(slab)) } else { hkl(slab) };
    let len = length(n);
    if len > 1e-9 {
        n / len
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

/// The interplanar spacing d_hkl (1 without a cell).
pub fn csc_dhkl(slab: &CellSlabPlane, source: &Structure) -> f32 {
    if !source.has_unit_cell {
        return 1.0;
    }
    let len = length(mul(&recip_matrix(source), hkl(slab)));
    if len > 1e-9 {
        1.0 / len
    } else {
        1.0
    }
}

/// The smallest gap between distinct atom planes in one structural period.
pub fn csc_atom_layer_spacing(slab: &CellSlabPlane, source: &Structure) -> f32 {
    let period = structural_period_multiplier(slab, source) as f32 * csc_dhkl(slab, source);
    if period < 1e-9 {
        return period;
    }
    let projections = all_layer_projections(slab, source, period);
    if projections.len() < 2 {
        return period;
    }
    let tol = period * 0.005;
    let mut min_gap = period - projections[projections.len() - 1] + projections[0];
    for w in projections.windows(2) {
        let gap = w[1] - w[0];
        if gap < min_gap {
            min_gap = gap;
        }
    }
    if min_gap > tol {
        min_gap
    } else {
        period
    }
}

/// The structural period N d_hkl.
pub fn csc_structural_dhkl(slab: &CellSlabPlane, source: &Structure) -> f32 {
    structural_period_multiplier(slab, source) as f32 * csc_dhkl(slab, source)
}

pub fn csc_d1(slab: &CellSlabPlane, source: &Structure) -> f32 {
    if !slab.use_periodic {
        return slab.d1;
    }
    slab.start_plane as f32 * csc_structural_dhkl(slab, source)
}

pub fn csc_d2(slab: &CellSlabPlane, source: &Structure) -> f32 {
    if !slab.use_periodic {
        return slab.d2;
    }
    (slab.start_plane + slab.n_periods) as f32 * csc_structural_dhkl(slab, source)
}

/// The mean atom position.
pub fn csc_centroid(s: &Structure) -> Vec3 {
    if s.atoms.is_empty() {
        return Vec3::ZERO;
    }
    let mut sum = Vec3::ZERO;
    for a in &s.atoms {
        sum += Vec3::new(a.x as f32, a.y as f32, a.z as f32);
    }
    sum / s.atoms.len() as f32
}

/// An nx x ny x nz supercell tiled symmetrically: ia in [-(nx/2), nx - nx/2).
pub fn csc_build_supercell(source: &Structure, nx: i32, ny: i32, nz: i32) -> Structure {
    if source.atoms.is_empty() {
        return Structure::default();
    }
    if nx == 1 && ny == 1 && nz == 1 {
        return source.clone();
    }
    let cv = source.cell_vectors;
    let mut sc = Structure { has_unit_cell: source.has_unit_cell, ..Default::default() };
    for i in 0..3 {
        sc.cell_vectors[0][i] = cv[0][i] * nx as f64;
        sc.cell_vectors[1][i] = cv[1][i] * ny as f64;
        sc.cell_vectors[2][i] = cv[2][i] * nz as f64;
    }
    let (ia_min, ia_max) = (-(nx / 2), nx - nx / 2);
    let (ib_min, ib_max) = (-(ny / 2), ny - ny / 2);
    let (ic_min, ic_max) = (-(nz / 2), nz - nz / 2);
    for k in 0..3 {
        sc.cell_offset[k] = source.cell_offset[k] + cv[0][k] * ia_min as f64 + cv[1][k] * ib_min as f64 + cv[2][k] * ic_min as f64;
    }
    let (a_d, b_d, c_d) = (DVec3::from(cv[0]), DVec3::from(cv[1]), DVec3::from(cv[2]));
    let origin = DVec3::from(source.cell_offset);
    let lattice = DMat3::from_cols(a_d, b_d, c_d);
    let lattice_inv = if glm::determinant_d(&lattice).abs() > 1e-12 { glm::inverse_d(&lattice) } else { DMat3::IDENTITY };
    let n = source.atoms.len();
    for ia in ia_min..ia_max {
        for ib in ib_min..ib_max {
            for ic in ic_min..ic_max {
                let (fa, fb, fc) = (ia as f64, ib as f64, ic as f64);
                let t = [
                    cv[0][0] * fa + cv[1][0] * fb + cv[2][0] * fc,
                    cv[0][1] * fa + cv[1][1] * fb + cv[2][1] * fc,
                    cv[0][2] * fa + cv[1][2] * fb + cv[2][2] * fc,
                ];
                for (index, atom) in source.atoms.iter().enumerate() {
                    let mut a = atom.clone();
                    if source.has_unit_cell {
                        // Wrap by whole cell vectors in double precision.
                        let frac = glm::mul_d(&lattice_inv, DVec3::new(a.x, a.y, a.z) - origin);
                        let cells = frac.floor();
                        let shift = a_d * cells.x + b_d * cells.y + c_d * cells.z;
                        a.x -= shift.x;
                        a.y -= shift.y;
                        a.z -= shift.z;
                    }
                    a.x += t[0];
                    a.y += t[1];
                    a.z += t[2];
                    sc.atoms.push(a);
                    if source.grain_colors.len() == n {
                        sc.grain_colors.push(source.grain_colors[index]);
                    }
                    if source.grain_region_ids.len() == n {
                        sc.grain_region_ids.push(source.grain_region_ids[index]);
                    }
                }
            }
        }
    }
    sc
}

/// The atoms of `sc` inside every slab; three slabs also define the cell.
pub fn csc_apply_slabs(sc: &Structure, slabs: &[CellSlabPlane], source: &Structure) -> Structure {
    if sc.atoms.is_empty() || slabs.is_empty() {
        return sc.clone();
    }
    let centroid = csc_centroid(sc);
    struct SlabData {
        n: Vec3,
        reference: Vec3,
        d1: f32,
        d2: f32,
        periodic: bool,
    }
    let enabled: Vec<SlabData> = slabs
        .iter()
        .map(|s| SlabData {
            n: csc_normal(s, source),
            reference: csc_slab_ref(s, source, centroid),
            d1: csc_d1(s, source),
            d2: csc_d2(s, source),
            periodic: s.use_periodic && source.has_unit_cell,
        })
        .collect();
    let mut result = sc.clone();
    result.atoms.clear();
    result.grain_colors.clear();
    result.grain_region_ids.clear();
    result.dislocation_loop_points.clear();
    result.dislocation_detection_done = false;
    let n = sc.atoms.len();
    for (index, atom) in sc.atoms.iter().enumerate() {
        let p = Vec3::new(atom.x as f32, atom.y as f32, atom.z as f32);
        let keep = enabled.iter().all(|sd| {
            let proj = dot(p - sd.reference, sd.n);
            if proj < sd.d1 - SLAB_TOL {
                return false;
            }
            // Periodic: half-open (the next period's first layer is out).
            !(if sd.periodic { proj >= sd.d2 - SLAB_TOL } else { proj > sd.d2 + SLAB_TOL })
        });
        if keep {
            result.atoms.push(atom.clone());
            if sc.grain_colors.len() == n {
                result.grain_colors.push(sc.grain_colors[index]);
            }
            if sc.grain_region_ids.len() == n {
                result.grain_region_ids.push(sc.grain_region_ids[index]);
            }
        }
    }
    result.has_unit_cell = false;
    if enabled.len() == 3 {
        // The parallelepiped of the slab normals from the corner where the
        // three lower planes meet.
        let a = Mat3::from_cols(enabled[0].n, enabled[1].n, enabled[2].n).transpose();
        if glm::determinant(&a).abs() > 1e-6 {
            let rhs = Vec3::new(
                enabled[0].d1 + dot(enabled[0].reference, enabled[0].n),
                enabled[1].d1 + dot(enabled[1].reference, enabled[1].n),
                enabled[2].d1 + dot(enabled[2].reference, enabled[2].n),
            );
            let a_inv = glm::inverse(&a);
            let corner = mul(&a_inv, rhs);
            result.has_unit_cell = true;
            result.cell_offset = [corner.x as f64, corner.y as f64, corner.z as f64];
            for j in 0..3 {
                let v = a_inv.col(j) * (enabled[j].d2 - enabled[j].d1);
                result.cell_vectors[j] = [v.x as f64, v.y as f64, v.z as f64];
            }
        }
    }
    // A cell vector or corner within 1e-3 A of a lattice translation
    // (twelfths of the cell vectors) is that translation exactly.
    if result.has_unit_cell && source.has_unit_cell {
        let scv = source.cell_vectors;
        let lattice = DMat3::from_cols(DVec3::from(scv[0]), DVec3::from(scv[1]), DVec3::from(scv[2]));
        if glm::determinant_d(&lattice).abs() > 1e-12 {
            let to_lattice = glm::inverse_d(&lattice);
            let snap = |v: DVec3, tolerance: f64| {
                let f = glm::mul_d(&to_lattice, v) * 12.0;
                let exact = glm::mul_d(&lattice, DVec3::new(f.x.round(), f.y.round(), f.z.round()) / 12.0);
                let d = exact - v;
                if d.dot(d).sqrt() < tolerance {
                    exact
                } else {
                    v
                }
            };
            for j in 0..3 {
                result.cell_vectors[j] = snap(DVec3::from(result.cell_vectors[j]), 1e-3).to_array();
            }
            let source_origin = DVec3::from(source.cell_offset);
            let corner = source_origin + snap(DVec3::from(result.cell_offset) - source_origin, 1e-3);
            result.cell_offset = corner.to_array();
        }
    }
    csc_deduplicate_atoms(&result, 0.01)
}

/// Removes coincident atoms (nearer than `tol`), keeping the first; with a
/// cell, positions compare wrapped into it (fractional 0 equals 1).
pub fn csc_deduplicate_atoms(s: &Structure, tol: f32) -> Structure {
    if s.atoms.is_empty() {
        return s.clone();
    }
    let tol2 = tol * tol;
    let inv = 1.0 / tol;
    let mut periodic = None;
    if s.has_unit_cell {
        let cell = glm::cell_matrix(s);
        if glm::determinant(&cell).abs() > 1e-8 {
            let origin = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
            periodic = Some((cell, glm::inverse(&cell), origin));
        }
    }
    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    let mut result = s.clone();
    result.atoms.clear();
    result.grain_colors.clear();
    result.grain_region_ids.clear();
    result.dislocation_loop_points.clear();
    result.dislocation_detection_done = false;
    let mut canonical_positions: Vec<Vec3> = Vec::new();
    let n = s.atoms.len();
    for (index, a) in s.atoms.iter().enumerate() {
        let mut canonical = Vec3::new(a.x as f32, a.y as f32, a.z as f32);
        if let Some((cell, inverse, origin)) = periodic {
            let mut frac = mul(&inverse, canonical - origin);
            for i in 0..3 {
                frac[i] -= frac[i].floor();
            }
            canonical = origin + mul(&cell, frac);
        }
        let key = |v: f32| (v * inv).floor() as i32;
        let (cx, cy, cz) = (key(canonical.x), key(canonical.y), key(canonical.z));
        let mut duplicate = false;
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(list) = grid.get(&(cx + dx, cy + dy, cz + dz)) {
                        for &idx in list {
                            let b = canonical_positions[idx];
                            let (ex, ey, ez) = (canonical.x - b.x, canonical.y - b.y, canonical.z - b.z);
                            if ex * ex + ey * ey + ez * ez < tol2 {
                                duplicate = true;
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
        if !duplicate {
            grid.entry((cx, cy, cz)).or_default().push(result.atoms.len());
            result.atoms.push(a.clone());
            if s.grain_colors.len() == n {
                result.grain_colors.push(s.grain_colors[index]);
            }
            if s.grain_region_ids.len() == n {
                result.grain_region_ids.push(s.grain_region_ids[index]);
            }
            canonical_positions.push(canonical);
        }
    }
    result
}

/// Whether the slab normals span 3D (the intersection is bounded).
pub fn csc_is_bounded(slabs: &[CellSlabPlane], source: &Structure) -> bool {
    let normals: Vec<Vec3> = slabs.iter().map(|s| csc_normal(s, source)).collect();
    if normals.len() < 3 {
        return false;
    }
    let (mut n0, mut n1) = (Vec3::ZERO, Vec3::ZERO);
    let mut found = 0;
    for &n in &normals {
        if found == 0 {
            if length(n) < 1e-6 {
                continue;
            }
            n0 = n;
            found += 1;
        } else if found == 1 {
            if length(cross(n0, n)) > 1e-4 {
                n1 = n;
                found += 1;
            }
        } else if dot(cross(n0, n1), n).abs() > 1e-4 {
            return true;
        }
    }
    false
}

/// Sutherland-Hodgman: the part of a polygon where dot(p, n) >= d.
pub fn csc_clip_halfspace(poly: &[Vec3], n: Vec3, plane_d: f32) -> Vec<Vec3> {
    let mut result = Vec::with_capacity(poly.len() + 1);
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        let da = dot(a, n) - plane_d;
        let db = dot(b, n) - plane_d;
        if da >= 0.0 {
            result.push(a);
        }
        if (da > 0.0 && db < 0.0) || (da < 0.0 && db > 0.0) {
            result.push(a + (da / (da - db)) * (b - a));
        }
    }
    result
}
