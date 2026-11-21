//! Surface slabs (the C++ algorithms/SurfaceBuilder): the shortest in-plane
//! lattice vectors of (hkl), the smallest out-of-plane step, layers stacked
//! along it, the normal made perpendicular and vacuum centred.

use af_core::cxx::sort_by;
use af_core::Structure;

use crate::csl::{
    grain_to_structure, is_orthogonal, make_supercell, make_supercell_diag, reduce_to_primitive, set_orthogonal_grain, structure_to_grain, M3,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceParams {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    pub layers: i32,
    pub vacuum: f64,
    pub nmax: i32,
    pub primitive_input: bool,
    pub primitive_symprec: f64,
}

impl Default for SurfaceParams {
    fn default() -> Self {
        Self { h: 1, k: 1, l: 1, layers: 4, vacuum: 15.0, nmax: 8, primitive_input: true, primitive_symprec: 1e-3 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SurfaceResult {
    pub success: bool,
    pub message: String,
    pub structure: Structure,
    pub atoms_per_layer: i32,
}

fn det_m3(m: &[[i32; 3]; 3]) -> i64 {
    let e = |i: usize, j: usize| m[i][j] as i64;
    e(0, 0) * (e(1, 1) * e(2, 2) - e(1, 2) * e(2, 1)) - e(0, 1) * (e(1, 0) * e(2, 2) - e(1, 2) * e(2, 0)) + e(0, 2) * (e(1, 0) * e(2, 1) - e(1, 1) * e(2, 0))
}

fn to_cartesian(v: &[i32; 3], cell: &M3) -> [f64; 3] {
    std::array::from_fn(|c| v[0] as f64 * cell[0][c] + v[1] as f64 * cell[1][c] + v[2] as f64 * cell[2][c])
}

fn length_of(v: &[i32; 3], cell: &M3) -> f64 {
    let c = to_cartesian(v, cell);
    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}

fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Miller indices in a reduced cell sharing the original's orientation
/// (M = reduced * original^-1 rational in twelfths); false otherwise.
fn transform_miller_indices(original: &M3, reduced: &M3, miller: &mut [i32; 3]) -> bool {
    let c = original;
    let det = c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0])
        + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0]);
    if det.abs() < 1e-12 {
        return false;
    }
    let inv = [
        [(c[1][1] * c[2][2] - c[1][2] * c[2][1]) / det, (c[0][2] * c[2][1] - c[0][1] * c[2][2]) / det, (c[0][1] * c[1][2] - c[0][2] * c[1][1]) / det],
        [(c[1][2] * c[2][0] - c[1][0] * c[2][2]) / det, (c[0][0] * c[2][2] - c[0][2] * c[2][0]) / det, (c[0][2] * c[1][0] - c[0][0] * c[1][2]) / det],
        [(c[1][0] * c[2][1] - c[1][1] * c[2][0]) / det, (c[0][1] * c[2][0] - c[0][0] * c[2][1]) / det, (c[0][0] * c[1][1] - c[0][1] * c[1][0]) / det],
    ];
    let mut transformed = [0.0f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut m = 0.0;
            for k in 0..3 {
                m += reduced[i][k] * inv[k][j];
            }
            if (m * 12.0 - (m * 12.0).round()).abs() > 1e-5 {
                return false;
            }
            transformed[i] += m * miller[j] as f64;
        }
    }
    let scaled = transformed.map(|t| (t * 12.0).round() as i64);
    let mut g = 0;
    for v in scaled {
        g = gcd(g, v.abs());
    }
    if g == 0 {
        return false;
    }
    for i in 0..3 {
        miller[i] = (scaled[i] / g) as i32;
    }
    true
}

fn is_multiple(a: &[i32; 3], b: &[i32; 3]) -> bool {
    a[0] * b[1] - a[1] * b[0] == 0 && a[0] * b[2] - a[2] * b[0] == 0 && a[1] * b[2] - a[2] * b[1] == 0
}

/// `ostream << double` (6 significant digits).
fn g6(x: f64) -> String {
    af_core::format::format_g(x, 6)
}

pub fn build_surface(source: &Structure, params: &SurfaceParams) -> SurfaceResult {
    let mut result = SurfaceResult::default();
    if params.h == 0 && params.k == 0 && params.l == 0 {
        result.message = "Miller indices (h,k,l) cannot all be zero.".into();
        return result;
    }
    if params.layers < 1 {
        result.message = "layers must be at least 1.".into();
        return result;
    }
    if params.nmax < 1 {
        result.message = "nmax must be at least 1.".into();
        return result;
    }
    if !source.has_unit_cell || source.atoms.is_empty() {
        result.message = "Structure needs a unit cell and at least one atom.".into();
        return result;
    }
    // Miller indices transform like lattice vectors: (h' k' l') = M (h k l).
    let mut miller = [params.h, params.k, params.l];
    let mut prepped = source.clone();
    if params.primitive_input && reduce_to_primitive(&mut prepped, params.primitive_symprec) {
        if !transform_miller_indices(&source.cell_vectors, &prepped.cell_vectors, &mut miller) {
            // The reduction rotated the cell: build from the input cell.
            prepped = source.clone();
            miller = [params.h, params.k, params.l];
        }
    }
    let grain = structure_to_grain(&prepped);
    let nmax = params.nmax;
    // 1. The two shortest independent in-plane vectors, Gauss-reduced.
    let mut in_plane: Vec<[i32; 3]> = Vec::new();
    for n0 in -nmax..=nmax {
        for n1 in -nmax..=nmax {
            for n2 in -nmax..=nmax {
                if (n0 == 0 && n1 == 0 && n2 == 0) || miller[0] * n0 + miller[1] * n1 + miller[2] * n2 != 0 {
                    continue;
                }
                in_plane.push([n0, n1, n2]);
            }
        }
    }
    let independent_error = "Could not find two independent in-plane lattice vectors; raise nmax.";
    if in_plane.len() < 2 {
        result.message = independent_error.into();
        return result;
    }
    sort_by(&mut in_plane, |a, b| length_of(a, &grain.cell) < length_of(b, &grain.cell));
    let mut v1 = in_plane[0];
    let Some(mut v2) = in_plane.iter().find(|c| !is_multiple(c, &v1)).copied() else {
        result.message = independent_error.into();
        return result;
    };
    for _ in 0..32 {
        if length_of(&v2, &grain.cell) < length_of(&v1, &grain.cell) {
            std::mem::swap(&mut v1, &mut v2);
        }
        let (c1, c2) = (to_cartesian(&v1, &grain.cell), to_cartesian(&v2, &grain.cell));
        let dot12 = c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
        let dot11 = c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2];
        if dot11 < 1e-12 {
            break;
        }
        let m = (dot12 / dot11).round() as i32;
        if m == 0 {
            break;
        }
        for c in 0..3 {
            v2[c] -= m * v1[c];
        }
    }
    // 2. The out-of-plane step skipping fewest planes, then shortest.
    let mut v3 = None;
    let (mut best_dot, mut best_len) = (0i64, 0.0f64);
    for n0 in -nmax..=nmax {
        for n1 in -nmax..=nmax {
            for n2 in -nmax..=nmax {
                let dot = miller[0] as i64 * n0 as i64 + miller[1] as i64 * n1 as i64 + miller[2] as i64 * n2 as i64;
                if dot == 0 {
                    continue;
                }
                let cand = [n0, n1, n2];
                if det_m3(&[v1, v2, cand]) == 0 {
                    continue;
                }
                let len = length_of(&cand, &grain.cell);
                let abs_dot = dot.abs();
                if v3.is_none() || abs_dot < best_dot || (abs_dot == best_dot && len < best_len) {
                    v3 = Some(cand);
                    best_dot = abs_dot;
                    best_len = len;
                }
            }
        }
    }
    let Some(v3) = v3 else {
        result.message = "Could not find an out-of-plane lattice vector completing the basis; raise nmax.".into();
        return result;
    };
    let slab_unit = make_supercell(&grain, &[v1, v2, v3]);
    result.atoms_per_layer = slab_unit.atoms.len() as i32;
    // Stack along the true lattice translation, then make the normal perpendicular.
    let mut slab = make_supercell_diag(&slab_unit, 1, 1, params.layers);
    if !is_orthogonal(&slab.cell, 0.01) {
        slab = set_orthogonal_grain(&slab, 2);
    }
    // 3. Vacuum centred along the normal.
    if params.vacuum > 0.0 && !slab.atoms.is_empty() {
        let c = slab.cell[2];
        let len = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        if len > 1e-9 {
            let unit = c.map(|x| x / len);
            let (mut low, mut high) = (0.0f64, 0.0f64);
            for (i, a) in slab.atoms.iter().enumerate() {
                let h = a.x * unit[0] + a.y * unit[1] + a.z * unit[2];
                if i == 0 || h < low {
                    low = h;
                }
                if i == 0 || h > high {
                    high = h;
                }
            }
            let new_len = (high - low) + params.vacuum;
            slab.cell[2] = unit.map(|u| u * new_len);
            let offset = params.vacuum * 0.5 - low;
            let shift = unit.map(|u| u * offset);
            for a in slab.atoms.iter_mut() {
                a.x += shift[0];
                a.y += shift[1];
                a.z += shift[2];
            }
        }
    }
    result.structure = grain_to_structure(&slab);
    result.success = true;
    result.message = format!(
        "Built a ({} {} {}) surface: {} layers, {} atoms, {} A vacuum.",
        params.h,
        params.k,
        params.l,
        params.layers,
        result.structure.atoms.len(),
        g6(params.vacuum)
    );
    result
}
