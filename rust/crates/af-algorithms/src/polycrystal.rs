//! Polycrystals by Voronoi tessellation (the C++ algorithms/PolyCrystalBuilder):
//! seeds in a periodic box, each grain a rotated tiling of the reference cell
//! kept where its seed is the nearest, duplicates across grain boundaries
//! merged, IPF-Z colours per grain.

use std::collections::HashMap;

use af_core::cxx_random::{uniform_real_f32, Mt19937};
use af_core::glm;
use af_core::{AtomSite, Structure};
use glam::{Mat3, Vec3};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GrainOrientationMode {
    #[default]
    AllRandom,
    AllSpecified,
    PartialSpecified,
}

/// Bunge Euler angles (degrees) of one grain.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GrainOrientation {
    pub grain_index: i32,
    pub phi1: f32,
    pub phi: f32,
    pub phi2: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PolyParams {
    pub size_x: f32,
    pub size_y: f32,
    pub size_z: f32,
    pub num_grains: i32,
    pub seed: i32,
    pub orientation_mode: GrainOrientationMode,
    pub specified_orientations: Vec<GrainOrientation>,
}

impl Default for PolyParams {
    fn default() -> Self {
        Self {
            size_x: 50.0,
            size_y: 50.0,
            size_z: 50.0,
            num_grains: 8,
            seed: 42,
            orientation_mode: GrainOrientationMode::AllRandom,
            specified_orientations: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PolyBuildResult {
    pub success: bool,
    pub message: String,
    pub input_atoms: i32,
    pub output_atoms: i32,
    pub num_grains: i32,
}

fn radians(deg: f32) -> f32 {
    deg * 0.01745329251994329576923690768489f32
}

/// Bunge (Z1-X-Z2) rotation; glm columns.
fn euler_to_matrix(phi1_deg: f32, phi_deg: f32, phi2_deg: f32) -> Mat3 {
    let (p1, p, p2) = (radians(phi1_deg), radians(phi_deg), radians(phi2_deg));
    let (c1, s1) = (p1.cos(), p1.sin());
    let (cp, sp) = (p.cos(), p.sin());
    let (c2, s2) = (p2.cos(), p2.sin());
    Mat3::from_cols(
        Vec3::new(c1 * c2 - s1 * s2 * cp, -c1 * s2 - s1 * c2 * cp, s1 * sp),
        Vec3::new(s1 * c2 + c1 * s2 * cp, -s1 * s2 + c1 * c2 * cp, -c1 * sp),
        Vec3::new(s2 * sp, c2 * sp, cp),
    )
}

/// A uniformly random rotation: a random unit quaternion through glm's
/// normalize and mat3_cast.
fn random_rotation(rng: &mut Mt19937) -> Mat3 {
    let u1 = uniform_real_f32(rng, 0.0, 1.0);
    let u2 = uniform_real_f32(rng, 0.0, 1.0) * 2.0 * 3.14159265358979f32;
    let u3 = uniform_real_f32(rng, 0.0, 1.0) * 2.0 * 3.14159265358979f32;
    let q0 = (1.0 - u1).sqrt() * u2.sin();
    let q1 = (1.0 - u1).sqrt() * u2.cos();
    let q2 = u1.sqrt() * u3.sin();
    let q3 = u1.sqrt() * u3.cos();
    // glm::quat(w, x, y, z) = (q3, q0, q1, q2), normalised.
    let (mut w, mut x, mut y, mut z) = (q3, q0, q1, q2);
    let len = ((w * w + x * x) + (y * y + z * z)).sqrt();
    if len <= 0.0 {
        (w, x, y, z) = (1.0, 0.0, 0.0, 0.0);
    } else {
        let o = 1.0 / len;
        (w, x, y, z) = (w * o, x * o, y * o, z * o);
    }
    let (qxx, qyy, qzz) = (x * x, y * y, z * z);
    let (qxz, qxy, qyz) = (x * z, x * y, y * z);
    let (qwx, qwy, qwz) = (w * x, w * y, w * z);
    Mat3::from_cols(
        Vec3::new(1.0 - 2.0 * (qyy + qzz), 2.0 * (qxy + qwz), 2.0 * (qxz - qwy)),
        Vec3::new(2.0 * (qxy - qwz), 1.0 - 2.0 * (qxx + qzz), 2.0 * (qyz + qwx)),
        Vec3::new(2.0 * (qxz + qwy), 2.0 * (qyz - qwx), 1.0 - 2.0 * (qxx + qyy)),
    )
}

fn wrap_periodic_delta(delta: f32, period: f32) -> f32 {
    if period <= 1e-8 {
        return delta;
    }
    delta - (delta / period).round() * period
}

fn wrap_to_box(mut value: f32, period: f32) -> f32 {
    if period <= 1e-8 {
        return value;
    }
    value -= (value / period).floor() * period;
    if value >= period {
        value -= period;
    }
    if value < 0.0 {
        value += period;
    }
    value
}

fn nearest_seed_periodic(pos: Vec3, seeds: &[Vec3], box_size: Vec3) -> usize {
    let delta = |s: Vec3| {
        Vec3::new(wrap_periodic_delta(pos.x - s.x, box_size.x), wrap_periodic_delta(pos.y - s.y, box_size.y), wrap_periodic_delta(pos.z - s.z, box_size.z))
    };
    let mut best = 0;
    let mut best_dist2 = glm::dot(delta(seeds[0]), delta(seeds[0]));
    for (i, s) in seeds.iter().enumerate().skip(1) {
        let d = delta(*s);
        let d2 = glm::dot(d, d);
        if d2 < best_dist2 {
            best_dist2 = d2;
            best = i;
        }
    }
    best
}

/// IPF-Z colour for cubic symmetry: [001] red, [011] green, [111] blue.
fn compute_ipf_color(rotation: &Mat3) -> [f32; 3] {
    let d = glm::mul(&rotation.transpose(), Vec3::new(0.0, 0.0, 1.0));
    let (mut h, mut k, mut l) = (d.x.abs(), d.y.abs(), d.z.abs());
    if h < k {
        std::mem::swap(&mut h, &mut k);
    }
    if h < l {
        std::mem::swap(&mut h, &mut l);
    }
    if k < l {
        std::mem::swap(&mut k, &mut l);
    }
    let len = (h * h + k * k + l * l).sqrt();
    if len < 1e-10 {
        return [0.5, 0.5, 0.5];
    }
    h /= len;
    k /= len;
    l /= len;
    let (r, g, b) = (h - k, k - l, l * 1.7320508f32);
    // std::max({r, g, b, 1e-6f}): the first largest.
    let mut max_c = r;
    for v in [g, b, 1e-6f32] {
        if max_c < v {
            max_c = v;
        }
    }
    [r / max_c, g / max_c, b / max_c]
}

/// Builds the polycrystal into `structure` (its atoms, cell and grain data).
pub fn build_polycrystal(structure: &mut Structure, reference: &Structure, params: &PolyParams, element_colors: &[[f32; 3]]) -> PolyBuildResult {
    let mut result = PolyBuildResult { input_atoms: reference.atoms.len() as i32, num_grains: params.num_grains, ..Default::default() };
    if reference.atoms.is_empty() {
        result.message = "Reference structure has no atoms.".into();
        return result;
    }
    if !reference.has_unit_cell {
        result.message = "Reference structure needs a unit cell for polycrystal tiling.".into();
        return result;
    }
    if params.num_grains < 1 {
        result.message = "Number of grains must be at least 1.".into();
        return result;
    }
    if params.size_x <= 0.0 || params.size_y <= 0.0 || params.size_z <= 0.0 {
        result.message = "Box dimensions must be positive.".into();
        return result;
    }
    let cv = &reference.cell_vectors;
    let row = |i: usize| Vec3::new(cv[i][0] as f32, cv[i][1] as f32, cv[i][2] as f32);
    let (a, b, c) = (row(0), row(1), row(2));
    let (la, lb, lc) = (glm::length(a), glm::length(b), glm::length(c));
    if la < 1e-8 || lb < 1e-8 || lc < 1e-8 {
        result.message = "Reference structure has degenerate lattice vectors.".into();
        return result;
    }
    let n = params.num_grains as usize;
    let mut rng = Mt19937::new(params.seed as u32);
    // glm::vec3(distX(rng), distY(rng), distZ(rng)): GCC draws z first.
    let seeds: Vec<Vec3> = (0..n)
        .map(|_| {
            let z = uniform_real_f32(&mut rng, 0.0, params.size_z);
            let y = uniform_real_f32(&mut rng, 0.0, params.size_y);
            let x = uniform_real_f32(&mut rng, 0.0, params.size_x);
            Vec3::new(x, y, z)
        })
        .collect();
    let box_size = Vec3::new(params.size_x, params.size_y, params.size_z);
    let rotations: Vec<Mat3> = match params.orientation_mode {
        GrainOrientationMode::AllRandom => (0..n).map(|_| random_rotation(&mut rng)).collect(),
        GrainOrientationMode::AllSpecified => (0..n)
            .map(|i| match params.specified_orientations.get(i) {
                Some(o) => euler_to_matrix(o.phi1, o.phi, o.phi2),
                None => Mat3::IDENTITY,
            })
            .collect(),
        GrainOrientationMode::PartialSpecified => {
            let mut r: Vec<Mat3> = (0..n).map(|_| random_rotation(&mut rng)).collect();
            for o in &params.specified_orientations {
                if o.grain_index >= 0 && o.grain_index < params.num_grains {
                    r[o.grain_index as usize] = euler_to_matrix(o.phi1, o.phi, o.phi2);
                }
            }
            r
        }
    };
    let box_diag = (params.size_x * params.size_x + params.size_y * params.size_y + params.size_z * params.size_z).sqrt();
    // std::min({la, lb, lc}): the first smallest.
    let mut min_len = la;
    for v in [lb, lc] {
        if v < min_len {
            min_len = v;
        }
    }
    let max_rep = 100.min((box_diag / min_len).ceil() as i32 + 2);
    let total_estimate = (2 * max_rep + 1) as i64 * (2 * max_rep + 1) as i64 * (2 * max_rep + 1) as i64 * reference.atoms.len() as i64;
    if total_estimate > 50_000_000 {
        result.message = format!("Tiling would generate ~{total_estimate} candidate atoms (limit 50M). Reduce box size or use a larger unit cell.");
        return result;
    }
    let rot_a: Vec<Vec3> = rotations.iter().map(|r| glm::mul(r, a)).collect();
    let rot_b: Vec<Vec3> = rotations.iter().map(|r| glm::mul(r, b)).collect();
    let rot_c: Vec<Vec3> = rotations.iter().map(|r| glm::mul(r, c)).collect();
    let ipf: Vec<[f32; 3]> = rotations.iter().map(compute_ipf_color).collect();
    const IN_BOX_TOL: f32 = 1e-3;
    const DUP_TOL: f32 = 0.08;
    const DUP_TOL_SQ: f32 = DUP_TOL * DUP_TOL;
    let cell_size = DUP_TOL;
    let nx = 1.max((params.size_x / cell_size).floor() as i32);
    let ny = 1.max((params.size_y / cell_size).floor() as i32);
    let nz = 1.max((params.size_z / cell_size).floor() as i32);

    // Per-grain tiling and Voronoi filtering.
    let candidates: Vec<Vec<AtomSite>> = (0..n)
        .map(|g| {
            let seed = seeds[g];
            let rot = &rotations[g];
            let reps = |v: Vec3| 100.min((box_diag / glm::length(v)).ceil() as i32 + 2);
            let (na, nb, nc) = (reps(rot_a[g]), reps(rot_b[g]), reps(rot_c[g]));
            let mut out = Vec::new();
            for ia in -na..=na {
                for ib in -nb..=nb {
                    for ic in -nc..=nc {
                        let offset = ia as f32 * rot_a[g] + ib as f32 * rot_b[g] + ic as f32 * rot_c[g];
                        for atom in &reference.atoms {
                            let ref_pos = Vec3::new(atom.x as f32, atom.y as f32, atom.z as f32);
                            let pos = seed + glm::mul(rot, ref_pos) + offset;
                            if pos.x < -IN_BOX_TOL
                                || pos.x >= params.size_x + IN_BOX_TOL
                                || pos.y < -IN_BOX_TOL
                                || pos.y >= params.size_y + IN_BOX_TOL
                                || pos.z < -IN_BOX_TOL
                                || pos.z >= params.size_z + IN_BOX_TOL
                            {
                                continue;
                            }
                            let wrapped = Vec3::new(wrap_to_box(pos.x, params.size_x), wrap_to_box(pos.y, params.size_y), wrap_to_box(pos.z, params.size_z));
                            if nearest_seed_periodic(wrapped, &seeds, box_size) != g {
                                continue;
                            }
                            let mut o = atom.clone();
                            o.x = wrapped.x as f64;
                            o.y = wrapped.y as f64;
                            o.z = wrapped.z as f64;
                            let z = o.atomic_number;
                            let colour = if z >= 0 && (z as usize) < element_colors.len() {
                                element_colors[z as usize]
                            } else {
                                af_core::elements::default_colour(z)
                            };
                            [o.r, o.g, o.b] = colour;
                            out.push(o);
                        }
                    }
                }
            }
            out
        })
        .collect();

    // Serial merge, dropping sites that duplicate accepted ones.
    let wrap_index = |idx: i32, n: i32| {
        if n <= 0 {
            return 0;
        }
        let i = idx % n;
        if i < 0 {
            i + n
        } else {
            i
        }
    };
    let cell_of = |p: Vec3| {
        (
            wrap_index((p.x / cell_size).floor() as i32, nx),
            wrap_index((p.y / cell_size).floor() as i32, ny),
            wrap_index((p.z / cell_size).floor() as i32, nz),
        )
    };
    let mut generated: Vec<AtomSite> = Vec::new();
    let mut grain_of: Vec<i32> = Vec::new();
    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    for (g, list) in candidates.into_iter().enumerate() {
        for atom in list {
            let p = Vec3::new(atom.x as f32, atom.y as f32, atom.z as f32);
            let centre = cell_of(p);
            let mut duplicate = false;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let key = (wrap_index(centre.0 + dx, nx), wrap_index(centre.1 + dy, ny), wrap_index(centre.2 + dz, nz));
                        let Some(ids) = grid.get(&key) else { continue };
                        for &idx in ids {
                            let e = &generated[idx];
                            if e.atomic_number != atom.atomic_number {
                                continue;
                            }
                            let ddx = wrap_periodic_delta(p.x - e.x as f32, params.size_x);
                            let ddy = wrap_periodic_delta(p.y - e.y as f32, params.size_y);
                            let ddz = wrap_periodic_delta(p.z - e.z as f32, params.size_z);
                            if ddx * ddx + ddy * ddy + ddz * ddz <= DUP_TOL_SQ {
                                duplicate = true;
                                break 'search;
                            }
                        }
                    }
                }
            }
            if duplicate {
                continue;
            }
            generated.push(atom);
            grain_of.push(g as i32);
            grid.entry(centre).or_default().push(generated.len() - 1);
        }
    }
    if generated.is_empty() {
        result.message = "No atoms generated. Try increasing the box size or check the reference structure.".into();
        return result;
    }
    result.output_atoms = generated.len() as i32;
    structure.atoms = generated;
    structure.has_unit_cell = true;
    structure.cell_offset = [0.0; 3];
    // The box lengths as typed: the shortest decimal that reads back as the same float.
    let typed = |v: f32| af_core::format::strtod(&af_core::format::format_g(v as f64, 7)).0;
    structure.cell_vectors = [[typed(params.size_x), 0.0, 0.0], [0.0, typed(params.size_y), 0.0], [0.0, 0.0, typed(params.size_z)]];
    structure.grain_colors = grain_of.iter().map(|&g| ipf[g as usize]).collect();
    structure.grain_region_ids = grain_of;
    result.success = true;
    result.message = format!("Polycrystal built: {} atoms, {} grains.", result.output_atoms, params.num_grains);
    result
}
