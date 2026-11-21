//! Nanocrystals (the C++ algorithms/NanoCrystalBuilder): a tiled, optionally
//! rotated reference crystal clipped by a shape, a triangle mesh or a Wulff
//! polyhedron (Barmparis et al., Beilstein J. Nanotechnol. 6, 361 (2015)),
//! with duplicate sites merged. Single precision as glm computes it.

use std::collections::HashMap;

use af_core::cxx::sort_by;
use af_core::glm;
use af_core::{AtomSite, Structure};
use glam::{DMat3, DVec3, Mat3, Vec3};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NanoShape {
    #[default]
    Sphere,
    Ellipsoid,
    Box,
    Cylinder,
    Octahedron,
    TruncatedOctahedron,
    Cuboctahedron,
    MeshModel,
}

pub const NUM_SHAPES: usize = 8;
pub const SHAPES: [NanoShape; NUM_SHAPES] = [
    NanoShape::Sphere,
    NanoShape::Ellipsoid,
    NanoShape::Box,
    NanoShape::Cylinder,
    NanoShape::Octahedron,
    NanoShape::TruncatedOctahedron,
    NanoShape::Cuboctahedron,
    NanoShape::MeshModel,
];

pub fn shape_label(s: NanoShape) -> &'static str {
    match s {
        NanoShape::Sphere => "Sphere",
        NanoShape::Ellipsoid => "Ellipsoid",
        NanoShape::Box => "Box",
        NanoShape::Cylinder => "Cylinder",
        NanoShape::Octahedron => "Octahedron",
        NanoShape::TruncatedOctahedron => "Truncated Octahedron",
        NanoShape::Cuboctahedron => "Cuboctahedron",
        NanoShape::MeshModel => "3D Model Fill",
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NanoGenerationMode {
    #[default]
    Shape,
    WulffConstruction,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WulffPlaneInput {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    pub surface_energy: f32,
}

impl Default for WulffPlaneInput {
    fn default() -> Self {
        Self { h: 1, k: 0, l: 0, surface_energy: 1.0 }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WulffPlane {
    pub family_index: i32,
    pub h: i32,
    pub k: i32,
    pub l: i32,
    pub distance: f32,
    pub normal: Vec3,
}

#[derive(Clone, Debug, PartialEq)]
pub struct WulffFace {
    pub family_index: i32,
    pub h: i32,
    pub k: i32,
    pub l: i32,
    pub distance: f32,
    pub normal: Vec3,
    pub vertices: Vec<Vec3>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WulffPreview {
    pub success: bool,
    pub message: String,
    pub center: Vec3,
    pub min_point: Vec3,
    pub max_point: Vec3,
    pub max_plane_distance: f32,
    pub bounding_radius: f32,
    pub planes: Vec<WulffPlane>,
    pub faces: Vec<WulffFace>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct NanoParams {
    pub generation_mode: NanoGenerationMode,
    pub shape: NanoShape,
    pub sphere_radius: f32,
    pub ellip_rx: f32,
    pub ellip_ry: f32,
    pub ellip_rz: f32,
    pub box_hx: f32,
    pub box_hy: f32,
    pub box_hz: f32,
    pub cyl_radius: f32,
    pub cyl_height: f32,
    pub cyl_axis: i32,
    pub oct_radius: f32,
    pub trunc_oct_radius: f32,
    pub trunc_oct_trunc: f32,
    pub cubo_radius: f32,
    pub model_hx: f32,
    pub model_hy: f32,
    pub model_hz: f32,
    pub model_scale: f32,
    pub auto_center_from_atoms: bool,
    pub cx: f32,
    pub cy: f32,
    pub cz: f32,
    pub apply_crystal_orientation: bool,
    pub use_miller_orientation: bool,
    pub orient_x_deg: f32,
    pub orient_y_deg: f32,
    pub orient_z_deg: f32,
    pub miller_h: i32,
    pub miller_k: i32,
    pub miller_l: i32,
    pub auto_replicate: bool,
    pub rep_a: i32,
    pub rep_b: i32,
    pub rep_c: i32,
    pub set_output_cell: bool,
    pub vacuum_padding: f32,
    pub wulff_planes: Vec<WulffPlaneInput>,
    pub wulff_max_radius: f32,
}

impl Default for NanoParams {
    fn default() -> Self {
        Self {
            generation_mode: NanoGenerationMode::Shape,
            shape: NanoShape::Sphere,
            sphere_radius: 15.0,
            ellip_rx: 15.0,
            ellip_ry: 12.0,
            ellip_rz: 10.0,
            box_hx: 15.0,
            box_hy: 15.0,
            box_hz: 15.0,
            cyl_radius: 12.0,
            cyl_height: 30.0,
            cyl_axis: 2,
            oct_radius: 15.0,
            trunc_oct_radius: 18.0,
            trunc_oct_trunc: 12.0,
            cubo_radius: 15.0,
            model_hx: 15.0,
            model_hy: 15.0,
            model_hz: 15.0,
            model_scale: 1.0,
            auto_center_from_atoms: true,
            cx: 0.0,
            cy: 0.0,
            cz: 0.0,
            apply_crystal_orientation: false,
            use_miller_orientation: false,
            orient_x_deg: 0.0,
            orient_y_deg: 0.0,
            orient_z_deg: 0.0,
            miller_h: 1,
            miller_k: 0,
            miller_l: 0,
            auto_replicate: true,
            rep_a: 5,
            rep_b: 5,
            rep_c: 5,
            set_output_cell: true,
            vacuum_padding: 5.0,
            wulff_planes: Vec::new(),
            wulff_max_radius: 20.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NanoBuildResult {
    pub success: bool,
    pub message: String,
    pub mode: NanoGenerationMode,
    pub input_atoms: i32,
    pub output_atoms: i32,
    pub shape: NanoShape,
    pub estimated_diameter: f32,
    pub rep_a: i32,
    pub rep_b: i32,
    pub rep_c: i32,
    pub tiling_used: bool,
    pub rep_clamped: bool,
    pub wulff_family_count: i32,
    pub wulff_face_count: i32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfExtents {
    pub hx: f32,
    pub hy: f32,
    pub hz: f32,
}

// ---------------------------------------------------------------- glm helpers

fn fmax(a: f32, b: f32) -> f32 {
    if a < b {
        b
    } else {
        a
    }
}
fn fmin(a: f32, b: f32) -> f32 {
    if b < a {
        b
    } else {
        a
    }
}
fn radians(deg: f32) -> f32 {
    deg * 0.01745329251994329576923690768489f32
}

type M4 = [[f32; 4]; 4];
const I4: M4 = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]];

fn col_add(a: [f32; 4], b: [f32; 4]) -> [f32; 4] {
    std::array::from_fn(|i| a[i] + b[i])
}
fn col_scale(a: [f32; 4], s: f32) -> [f32; 4] {
    std::array::from_fn(|i| a[i] * s)
}

/// glm::rotate(m, angle, axis).
fn rotate(m: &M4, angle: f32, v: Vec3) -> M4 {
    let (c, s) = (angle.cos(), angle.sin());
    let axis = glm::normalize(v);
    let temp = (1.0 - c) * axis;
    let r = [
        [c + temp.x * axis.x, temp.x * axis.y + s * axis.z, temp.x * axis.z - s * axis.y],
        [temp.y * axis.x - s * axis.z, c + temp.y * axis.y, temp.y * axis.z + s * axis.x],
        [temp.z * axis.x + s * axis.y, temp.z * axis.y - s * axis.x, c + temp.z * axis.z],
    ];
    let col = |i: usize| col_add(col_add(col_scale(m[0], r[i][0]), col_scale(m[1], r[i][1])), col_scale(m[2], r[i][2]));
    [col(0), col(1), col(2), m[3]]
}

/// glm mat4 * mat4 (the non-SIMD path).
fn mul4(a: &M4, b: &M4) -> M4 {
    std::array::from_fn(|j| {
        let mut t = col_scale(a[0], b[j][0]);
        t = col_add(t, col_scale(a[1], b[j][1]));
        t = col_add(t, col_scale(a[2], b[j][2]));
        col_add(t, col_scale(a[3], b[j][3]))
    })
}

fn mat3_of(m: &M4) -> Mat3 {
    Mat3::from_cols(Vec3::new(m[0][0], m[0][1], m[0][2]), Vec3::new(m[1][0], m[1][1], m[1][2]), Vec3::new(m[2][0], m[2][1], m[2][2]))
}

fn row_f32(s: &Structure, i: usize) -> Vec3 {
    Vec3::new(s.cell_vectors[i][0] as f32, s.cell_vectors[i][1] as f32, s.cell_vectors[i][2] as f32)
}

// ---------------------------------------------------------------- geometry

fn ray_intersects_triangle(orig: Vec3, dir: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> bool {
    let eps = 1e-8f32;
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let p = glm::cross(dir, e2);
    let det = glm::dot(e1, p);
    if det.abs() < eps {
        return false;
    }
    let inv_det = 1.0 / det;
    let t = orig - v0;
    let u = glm::dot(t, p) * inv_det;
    if u < 0.0 || u > 1.0 {
        return false;
    }
    let q = glm::cross(t, e1);
    let v = glm::dot(dir, q) * inv_det;
    if v < 0.0 || (u + v) > 1.0 {
        return false;
    }
    glm::dot(e2, q) * inv_det > eps
}

fn is_inside_mesh_model(p: Vec3, params: &NanoParams, vertices: &[Vec3], indices: &[u32]) -> bool {
    if vertices.is_empty() || indices.len() < 3 {
        return false;
    }
    let inv_scale = if params.model_scale > 1e-8 { 1.0 / params.model_scale } else { 1.0 };
    let pm = p * inv_scale;
    let mut votes = 0;
    for dir in [Vec3::X, Vec3::Y, Vec3::Z] {
        let mut hits = 0;
        let mut i = 0;
        while i + 2 < indices.len() {
            let (i0, i1, i2) = (indices[i] as usize, indices[i + 1] as usize, indices[i + 2] as usize);
            i += 3;
            if i0 >= vertices.len() || i1 >= vertices.len() || i2 >= vertices.len() {
                continue;
            }
            if ray_intersects_triangle(pm, dir, vertices[i0], vertices[i1], vertices[i2]) {
                hits += 1;
            }
        }
        if hits % 2 == 1 {
            votes += 1;
        }
    }
    votes >= 2
}

fn build_reciprocal_basis(a: Vec3, b: Vec3, c: Vec3) -> (Vec3, Vec3, Vec3) {
    let volume = glm::dot(a, glm::cross(b, c));
    (glm::cross(b, c) / volume, glm::cross(c, a) / volume, glm::cross(a, b) / volume)
}

fn invert_integer_matrix(m: &[[i32; 3]; 3]) -> Option<[[i32; 3]; 3]> {
    let (a, b, c) = (m[0][0], m[0][1], m[0][2]);
    let (d, e, f) = (m[1][0], m[1][1], m[1][2]);
    let (g, h, i) = (m[2][0], m[2][1], m[2][2]);
    let c00 = e * i - f * h;
    let c01 = -(d * i - f * g);
    let c02 = d * h - e * g;
    let c10 = -(b * i - c * h);
    let c11 = a * i - c * g;
    let c12 = -(a * h - b * g);
    let c20 = b * f - c * e;
    let c21 = -(a * f - c * d);
    let c22 = a * e - b * d;
    let det = a * c00 + b * c01 + c * c02;
    if det != 1 && det != -1 {
        return None;
    }
    Some([[c00 / det, c10 / det, c20 / det], [c01 / det, c11 / det, c21 / det], [c02 / det, c12 / det, c22 / det]])
}

fn build_miller_plane_normal(reference: &Structure, h: i32, k: i32, l: i32) -> Result<Vec3, String> {
    if !reference.has_unit_cell {
        return Err("Miller orientation requires a reference structure with a unit cell.".into());
    }
    let (a, b, c) = (row_f32(reference, 0), row_f32(reference, 1), row_f32(reference, 2));
    let volume = glm::dot(a, glm::cross(b, c));
    if volume.abs() <= 1e-8 {
        return Err("Reference structure has a degenerate lattice.".into());
    }
    let (sa, sb, sc) = build_reciprocal_basis(a, b, c);
    let normal = h as f32 * sa + k as f32 * sb + l as f32 * sc;
    let len = glm::length(normal);
    if len <= 1e-8 {
        return Err("Miller indices do not define a valid plane normal for the current lattice.".into());
    }
    Ok(normal / len)
}

fn rotation_aligning_direction_to_z(from: Vec3) -> Mat3 {
    let source = glm::normalize(from);
    let target = Vec3::Z;
    let dot_st = glm::dot(source, target);
    if dot_st > 0.999999 {
        return Mat3::IDENTITY;
    }
    if dot_st < -0.999999 {
        let mut axis = glm::cross(source, Vec3::X);
        if glm::length(axis) < 1e-6 {
            axis = glm::cross(source, Vec3::Y);
        }
        axis = glm::normalize(axis);
        return mat3_of(&rotate(&I4, 3.14159265358979323846f32, axis));
    }
    let mut axis = glm::cross(source, target);
    let len = glm::length(axis);
    if len < 1e-6 {
        return Mat3::IDENTITY;
    }
    axis /= len;
    let clamped = fmax(-1.0, fmin(1.0, dot_st));
    mat3_of(&rotate(&I4, clamped.acos(), axis))
}

fn build_orientation_matrix(reference: &Structure, params: &NanoParams) -> Result<Mat3, String> {
    if !params.apply_crystal_orientation {
        return Ok(Mat3::IDENTITY);
    }
    if params.use_miller_orientation {
        if params.miller_h == 0 && params.miller_k == 0 && params.miller_l == 0 {
            return Err("Miller indices [h k l] cannot all be zero.".into());
        }
        let normal = build_miller_plane_normal(reference, params.miller_h, params.miller_k, params.miller_l)?;
        return Ok(rotation_aligning_direction_to_z(normal));
    }
    let (rx, ry, rz) = (radians(params.orient_x_deg), radians(params.orient_y_deg), radians(params.orient_z_deg));
    let r = mul4(&mul4(&rotate(&I4, rz, Vec3::Z), &rotate(&I4, ry, Vec3::Y)), &rotate(&I4, rx, Vec3::X));
    Ok(mat3_of(&r))
}

/// The mean atom position (single precision).
pub fn compute_atom_centroid(atoms: &[AtomSite]) -> Vec3 {
    if atoms.is_empty() {
        return Vec3::ZERO;
    }
    let mut sum = Vec3::ZERO;
    for a in atoms {
        sum += Vec3::new(a.x as f32, a.y as f32, a.z as f32);
    }
    sum / atoms.len() as f32
}

fn build_clip_center(reference: &Structure, params: &NanoParams) -> Vec3 {
    if params.auto_center_from_atoms {
        compute_atom_centroid(&reference.atoms)
    } else {
        Vec3::new(params.cx, params.cy, params.cz)
    }
}

fn gcd(a: i32, b: i32) -> i32 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn reduce_miller(h: [i32; 3]) -> [i32; 3] {
    let d = gcd(h[0].abs(), gcd(h[1].abs(), h[2].abs()));
    if d <= 1 {
        h
    } else {
        [h[0] / d, h[1] / d, h[2] / d]
    }
}

fn transform_miller(rotation: &[[i32; 3]; 3], h: [i32; 3]) -> [i32; 3] {
    let Some(inv) = invert_integer_matrix(rotation) else { return h };
    [
        inv[0][0] * h[0] + inv[1][0] * h[1] + inv[2][0] * h[2],
        inv[0][1] * h[0] + inv[1][1] * h[1] + inv[2][1] * h[2],
        inv[0][2] * h[0] + inv[1][2] * h[1] + inv[2][2] * h[2],
    ]
}

/// The distinct rotations of the lattice holohedry (spglib on one atom).
fn collect_symmetry_rotations(reference: &Structure) -> Result<Vec<[[i32; 3]; 3]>, String> {
    if !reference.has_unit_cell {
        return Err("Wulff construction requires a reference structure with a unit cell.".into());
    }
    if reference.atoms.is_empty() {
        return Err("Reference structure has no atoms.".into());
    }
    let lattice = Mat3::from_cols(row_f32(reference, 0), row_f32(reference, 1), row_f32(reference, 2));
    if glm::determinant(&lattice).abs() <= 1e-10 {
        return Err("Reference structure has a degenerate lattice.".into());
    }
    let Some(dataset) = af_science::spg::dataset(&reference.cell_vectors, &[[0.0; 3]], &[1], 1e-5) else {
        return Err("spglib could not determine symmetry for the reference structure.".into());
    };
    let mut rotations: Vec<[[i32; 3]; 3]> = Vec::new();
    for r in dataset.rotations {
        if !rotations.contains(&r) {
            rotations.push(r);
        }
    }
    if rotations.is_empty() {
        return Err("No symmetry rotations were returned by spglib.".into());
    }
    Ok(rotations)
}

fn append_unique_point(points: &mut Vec<Vec3>, candidate: Vec3, tolerance: f32) -> bool {
    let t2 = tolerance * tolerance;
    for p in points.iter() {
        let d = *p - candidate;
        if glm::dot(d, d) <= t2 {
            return false;
        }
    }
    points.push(candidate);
    true
}

fn sort_face_vertices(vertices: &mut [Vec3], normal: Vec3) {
    if vertices.len() < 3 {
        return;
    }
    let mut centroid = Vec3::ZERO;
    for v in vertices.iter() {
        centroid += *v;
    }
    centroid /= vertices.len() as f32;
    let reference_axis = if normal.z.abs() < 0.9 { Vec3::Z } else { Vec3::Y };
    let u = glm::normalize(glm::cross(reference_axis, normal));
    let v = glm::normalize(glm::cross(normal, u));
    let angle = |p: &Vec3| {
        let d = *p - centroid;
        glm::dot(d, v).atan2(glm::dot(d, u))
    };
    sort_by(vertices, |l, r| angle(l) < angle(r));
}

fn is_inside_wulff(local: Vec3, preview: &WulffPreview, tolerance: f32) -> bool {
    !preview.planes.iter().any(|p| glm::dot(p.normal, local) > p.distance + tolerance)
}

pub fn compute_bounding_radius(p: &NanoParams) -> f32 {
    if p.generation_mode == NanoGenerationMode::WulffConstruction {
        return fmax(p.wulff_max_radius, 0.1);
    }
    match p.shape {
        NanoShape::Sphere => p.sphere_radius,
        // std::max({rx, ry, rz}): the first largest.
        NanoShape::Ellipsoid => {
            let mut m = p.ellip_rx;
            for v in [p.ellip_ry, p.ellip_rz] {
                if m < v {
                    m = v;
                }
            }
            m
        }
        NanoShape::Box => (p.box_hx * p.box_hx + p.box_hy * p.box_hy + p.box_hz * p.box_hz).sqrt(),
        NanoShape::Cylinder => (p.cyl_radius * p.cyl_radius + (p.cyl_height * 0.5) * (p.cyl_height * 0.5)).sqrt(),
        NanoShape::Octahedron => p.oct_radius,
        NanoShape::TruncatedOctahedron => fmax(p.trunc_oct_radius, p.trunc_oct_trunc * 3.0f32.sqrt()),
        NanoShape::Cuboctahedron => p.cubo_radius,
        NanoShape::MeshModel => p.model_scale * (p.model_hx * p.model_hx + p.model_hy * p.model_hy + p.model_hz * p.model_hz).sqrt(),
    }
}

pub fn compute_shape_half_extents(p: &NanoParams) -> HalfExtents {
    let all = |r: f32| HalfExtents { hx: r, hy: r, hz: r };
    if p.generation_mode == NanoGenerationMode::WulffConstruction {
        return all(fmax(p.wulff_max_radius, 0.1));
    }
    match p.shape {
        NanoShape::Sphere => all(p.sphere_radius),
        NanoShape::Ellipsoid => HalfExtents { hx: p.ellip_rx, hy: p.ellip_ry, hz: p.ellip_rz },
        NanoShape::Box => HalfExtents { hx: p.box_hx, hy: p.box_hy, hz: p.box_hz },
        NanoShape::Cylinder => match p.cyl_axis {
            0 => HalfExtents { hx: p.cyl_height * 0.5, hy: p.cyl_radius, hz: p.cyl_radius },
            1 => HalfExtents { hx: p.cyl_radius, hy: p.cyl_height * 0.5, hz: p.cyl_radius },
            _ => HalfExtents { hx: p.cyl_radius, hy: p.cyl_radius, hz: p.cyl_height * 0.5 },
        },
        NanoShape::Octahedron => all(p.oct_radius),
        NanoShape::TruncatedOctahedron => all(fmin(p.trunc_oct_trunc, p.trunc_oct_radius)),
        NanoShape::Cuboctahedron => all(p.cubo_radius),
        NanoShape::MeshModel => HalfExtents { hx: p.model_hx * p.model_scale, hy: p.model_hy * p.model_scale, hz: p.model_hz * p.model_scale },
    }
}

pub fn is_inside_shape(p: Vec3, params: &NanoParams, vertices: &[Vec3], indices: &[u32]) -> bool {
    match params.shape {
        NanoShape::Sphere => glm::dot(p, p) <= params.sphere_radius * params.sphere_radius,
        NanoShape::Ellipsoid => {
            let (fx, fy, fz) = (p.x / params.ellip_rx, p.y / params.ellip_ry, p.z / params.ellip_rz);
            fx * fx + fy * fy + fz * fz <= 1.0
        }
        NanoShape::Box => p.x.abs() <= params.box_hx && p.y.abs() <= params.box_hy && p.z.abs() <= params.box_hz,
        NanoShape::Cylinder => {
            let (r2, ax) = match params.cyl_axis {
                0 => (p.y * p.y + p.z * p.z, p.x),
                1 => (p.x * p.x + p.z * p.z, p.y),
                _ => (p.x * p.x + p.y * p.y, p.z),
            };
            r2 <= params.cyl_radius * params.cyl_radius && ax.abs() <= params.cyl_height * 0.5
        }
        NanoShape::Octahedron => p.x.abs() + p.y.abs() + p.z.abs() <= params.oct_radius,
        NanoShape::TruncatedOctahedron => {
            p.x.abs() + p.y.abs() + p.z.abs() <= params.trunc_oct_radius
                && p.x.abs() <= params.trunc_oct_trunc
                && p.y.abs() <= params.trunc_oct_trunc
                && p.z.abs() <= params.trunc_oct_trunc
        }
        NanoShape::Cuboctahedron => {
            p.x.abs() + p.y.abs() <= params.cubo_radius && p.y.abs() + p.z.abs() <= params.cubo_radius && p.x.abs() + p.z.abs() <= params.cubo_radius
        }
        NanoShape::MeshModel => is_inside_mesh_model(p, params, vertices, indices),
    }
}

/// The Wulff polyhedron of the reference's symmetry-expanded planes.
pub fn compute_wulff_preview(reference: &Structure, params: &NanoParams) -> WulffPreview {
    let mut preview = WulffPreview { center: build_clip_center(reference, params), ..Default::default() };
    if params.wulff_planes.is_empty() {
        preview.message = "Add at least one Miller-index and surface-energy entry.".into();
        return preview;
    }
    if !reference.has_unit_cell {
        preview.message = "Wulff construction requires a periodic reference structure with a unit cell.".into();
        return preview;
    }
    let orientation = match build_orientation_matrix(reference, params) {
        Ok(m) => m,
        Err(e) => {
            preview.message = e;
            return preview;
        }
    };
    let a = glm::mul(&orientation, row_f32(reference, 0));
    let b = glm::mul(&orientation, row_f32(reference, 1));
    let c = glm::mul(&orientation, row_f32(reference, 2));
    if glm::dot(a, glm::cross(b, c)).abs() <= 1e-8 {
        preview.message = "Reference structure has a degenerate lattice.".into();
        return preview;
    }
    let mut max_energy = 0.0f32;
    for input in &params.wulff_planes {
        if input.surface_energy <= 0.0 {
            preview.message = "Surface energies must be positive.".into();
            return preview;
        }
        if input.h == 0 && input.k == 0 && input.l == 0 {
            preview.message = "Miller indices [h k l] cannot all be zero.".into();
            return preview;
        }
        max_energy = fmax(max_energy, input.surface_energy);
    }
    if params.wulff_max_radius <= 0.0 {
        preview.message = "Maximum Wulff radius must be positive.".into();
        return preview;
    }
    let rotations = match collect_symmetry_rotations(reference) {
        Ok(r) => r,
        Err(e) => {
            preview.message = e;
            return preview;
        }
    };
    let (sa, sb, sc) = build_reciprocal_basis(a, b, c);
    let mut index_of: HashMap<[i32; 3], usize> = HashMap::new();
    for (family, input) in params.wulff_planes.iter().enumerate() {
        let base = reduce_miller([input.h, input.k, input.l]);
        let distance = params.wulff_max_radius * (input.surface_energy / max_energy);
        for rotation in &rotations {
            let t = reduce_miller(transform_miller(rotation, base));
            if t == [0, 0, 0] {
                continue;
            }
            let reciprocal = t[0] as f32 * sa + t[1] as f32 * sb + t[2] as f32 * sc;
            let len = glm::length(reciprocal);
            if len <= 1e-8 {
                continue;
            }
            match index_of.get(&t) {
                None => {
                    index_of.insert(t, preview.planes.len());
                    preview.planes.push(WulffPlane { family_index: family as i32, h: t[0], k: t[1], l: t[2], distance, normal: reciprocal / len });
                }
                Some(&i) => {
                    if distance < preview.planes[i].distance {
                        let plane = &mut preview.planes[i];
                        plane.family_index = family as i32;
                        plane.distance = distance;
                        plane.normal = reciprocal / len;
                    }
                }
            }
        }
    }
    if preview.planes.len() < 4 {
        preview.message = "Wulff construction produced too few symmetry-equivalent planes.".into();
        return preview;
    }
    let intersection_tolerance = fmax(1e-3, params.wulff_max_radius * 1e-4);
    let face_tolerance = fmax(2e-3, params.wulff_max_radius * 2e-4);
    let mut vertices: Vec<Vec3> = Vec::new();
    let planes = &preview.planes;
    for i in 0..planes.len() {
        for j in i + 1..planes.len() {
            for k in j + 1..planes.len() {
                let (ni, nj, nk) = (planes[i].normal, planes[j].normal, planes[k].normal);
                // system[c][r]: row r is plane r's normal.
                let system = Mat3::from_cols(Vec3::new(ni.x, nj.x, nk.x), Vec3::new(ni.y, nj.y, nk.y), Vec3::new(ni.z, nj.z, nk.z));
                if glm::determinant(&system).abs() <= 1e-6 {
                    continue;
                }
                let distances = Vec3::new(planes[i].distance, planes[j].distance, planes[k].distance);
                let local = glm::mul(&glm::inverse(&system), distances);
                if !is_inside_wulff(local, &preview, intersection_tolerance) {
                    continue;
                }
                append_unique_point(&mut vertices, preview.center + local, intersection_tolerance);
            }
        }
    }
    if vertices.len() < 4 {
        preview.message = "No valid Wulff polyhedron could be formed from the supplied planes and energies.".into();
        return preview;
    }
    preview.min_point = Vec3::splat(f32::MAX);
    preview.max_point = Vec3::splat(-f32::MAX);
    preview.bounding_radius = 0.0;
    preview.max_plane_distance = 0.0;
    for v in &vertices {
        preview.min_point = glm::min(preview.min_point, *v);
        preview.max_point = glm::max(preview.max_point, *v);
        preview.bounding_radius = fmax(preview.bounding_radius, glm::length(*v - preview.center));
    }
    for p in &preview.planes {
        preview.max_plane_distance = fmax(preview.max_plane_distance, p.distance);
    }
    let mut faces = Vec::new();
    for plane in &preview.planes {
        let mut face_vertices: Vec<Vec3> =
            vertices.iter().copied().filter(|v| (glm::dot(plane.normal, *v - preview.center) - plane.distance).abs() <= face_tolerance).collect();
        if face_vertices.len() < 3 {
            continue;
        }
        sort_face_vertices(&mut face_vertices, plane.normal);
        faces.push(WulffFace {
            family_index: plane.family_index,
            h: plane.h,
            k: plane.k,
            l: plane.l,
            distance: plane.distance,
            normal: plane.normal,
            vertices: face_vertices,
        });
    }
    preview.faces = faces;
    if preview.faces.is_empty() {
        preview.message = "Wulff construction generated vertices but no valid faces.".into();
        return preview;
    }
    preview.success = true;
    preview.message = format!("Wulff preview: {} faces from {} symmetry-expanded planes.", preview.faces.len(), preview.planes.len());
    preview
}

/// The reference atom rotated about the centre and moved by whole cells,
/// in double precision.
fn place_exactly(atom: &AtomSite, center: Vec3, orientation: &Mat3, reference: &Structure, ia: i32, ib: i32, ic: i32) -> DVec3 {
    let origin = center.as_dvec3();
    let rotation = DMat3::from_cols(orientation.x_axis.as_dvec3(), orientation.y_axis.as_dvec3(), orientation.z_axis.as_dvec3());
    let mut offset = DVec3::ZERO;
    if ia != 0 || ib != 0 || ic != 0 {
        let cv = &reference.cell_vectors;
        offset = ia as f64 * DVec3::from(cv[0]) + ib as f64 * DVec3::from(cv[1]) + ic as f64 * DVec3::from(cv[2]);
    }
    origin + glm::mul_d(&rotation, DVec3::new(atom.x, atom.y, atom.z) - origin + offset)
}

fn coloured(atom: &AtomSite, exact: DVec3, colours: &[[f32; 3]]) -> AtomSite {
    let mut out = atom.clone();
    out.x = exact.x;
    out.y = exact.y;
    out.z = exact.z;
    let z = out.atomic_number;
    [out.r, out.g, out.b] = if z >= 0 && (z as usize) < colours.len() { colours[z as usize] } else { af_core::elements::default_colour(z) };
    out
}

/// Tiles the oriented reference, keeping atoms `inside` (relative to the centre).
#[allow(clippy::too_many_arguments)]
fn tile(
    reference: &Structure,
    orientation: &Mat3,
    center: Vec3,
    reps: (i32, i32, i32),
    cell: (Vec3, Vec3, Vec3),
    colours: &[[f32; 3]],
    inside: &dyn Fn(Vec3) -> bool,
) -> Vec<AtomSite> {
    let (na, nb, nc) = reps;
    let (a, b, c) = cell;
    let mut out = Vec::new();
    for ia in -na..=na {
        for ib in -nb..=nb {
            for ic in -nc..=nc {
                let offset = ia as f32 * a + ib as f32 * b + ic as f32 * c;
                for atom in &reference.atoms {
                    let local = Vec3::new(atom.x as f32 - center.x, atom.y as f32 - center.y, atom.z as f32 - center.z);
                    let pos = center + glm::mul(orientation, local) + offset;
                    if !inside(pos - center) {
                        continue;
                    }
                    out.push(coloured(atom, place_exactly(atom, center, orientation, reference, ia, ib, ic), colours));
                }
            }
        }
    }
    out
}

/// Without a cell: each reference atom once.
fn untiled(reference: &Structure, orientation: &Mat3, center: Vec3, colours: &[[f32; 3]], inside: &dyn Fn(Vec3) -> bool) -> Vec<AtomSite> {
    let mut out = Vec::new();
    for atom in &reference.atoms {
        let local = Vec3::new(atom.x as f32 - center.x, atom.y as f32 - center.y, atom.z as f32 - center.z);
        let pos = center + glm::mul(orientation, local);
        if !inside(pos - center) {
            continue;
        }
        out.push(coloured(atom, place_exactly(atom, center, orientation, reference, 0, 0, 0), colours));
    }
    out
}

/// Duplicates within 0.01 A merged (the first kept), through a 0.5 A hash grid.
fn deduplicate(atoms: Vec<AtomSite>) -> Vec<AtomSite> {
    let inv_cell = 1.0f32 / 0.5f32;
    let tol2 = 0.01f32 * 0.01f32;
    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    let mut unique: Vec<AtomSite> = Vec::with_capacity(atoms.len());
    for atom in atoms {
        let p = Vec3::new(atom.x as f32, atom.y as f32, atom.z as f32);
        let cell = ((p.x * inv_cell).floor() as i32, (p.y * inv_cell).floor() as i32, (p.z * inv_cell).floor() as i32);
        let mut duplicate = false;
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(ids) = grid.get(&(cell.0 + dx, cell.1 + dy, cell.2 + dz)) {
                        for &idx in ids {
                            let u = &unique[idx];
                            let d = p - Vec3::new(u.x as f32, u.y as f32, u.z as f32);
                            if glm::dot(d, d) < tol2 {
                                duplicate = true;
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
        if !duplicate {
            grid.entry(cell).or_default().push(unique.len());
            unique.push(atom);
        }
    }
    unique
}

fn clear_metadata(structure: &mut Structure) {
    structure.grain_colors.clear();
    structure.grain_region_ids.clear();
    structure.pbc_boundary_tol = 0.0;
    structure.ipf_load_status.clear();
    structure.dislocation_loop_points.clear();
    structure.dislocation_detection_done = false;
}

pub fn build_nanocrystal(
    structure: &mut Structure,
    reference: &Structure,
    params: &NanoParams,
    element_colors: &[[f32; 3]],
    model_vertices: &[Vec3],
    model_indices: &[u32],
) -> NanoBuildResult {
    let mut result = NanoBuildResult { mode: params.generation_mode, shape: params.shape, input_atoms: reference.atoms.len() as i32, ..Default::default() };
    const MAX_REPS: i32 = 40;

    if params.generation_mode == NanoGenerationMode::WulffConstruction {
        if reference.atoms.is_empty() {
            result.message = "Reference structure has no atoms.".into();
            return result;
        }
        let preview = compute_wulff_preview(reference, params);
        if !preview.success {
            result.message = preview.message;
            return result;
        }
        let orientation = match build_orientation_matrix(reference, params) {
            Ok(m) => m,
            Err(e) => {
                result.message = e;
                return result;
            }
        };
        let center = preview.center;
        let max_r = fmax(preview.bounding_radius, 0.1);
        let inside_tolerance = fmax(1e-3, preview.max_plane_distance * 1e-4);
        let inside = |p: Vec3| is_inside_wulff(p, &preview, inside_tolerance);
        let generated = if reference.has_unit_cell {
            let a = glm::mul(&orientation, row_f32(reference, 0));
            let b = glm::mul(&orientation, row_f32(reference, 1));
            let c = glm::mul(&orientation, row_f32(reference, 2));
            let (la, lb, lc) = (glm::length(a), glm::length(b), glm::length(c));
            if la < 1e-8 || lb < 1e-8 || lc < 1e-8 {
                result.message = "Reference structure has degenerate lattice vectors.".into();
                return result;
            }
            let mut clamped = false;
            let mut safe_rep = |len: f32| {
                let n = (max_r / len).ceil() as i32 + 2;
                if n > MAX_REPS {
                    clamped = true;
                    MAX_REPS
                } else {
                    n
                }
            };
            let (na, nb, nc) = if params.auto_replicate {
                (safe_rep(la), safe_rep(lb), safe_rep(lc))
            } else {
                (params.rep_a.max(1), params.rep_b.max(1), params.rep_c.max(1))
            };
            result.rep_a = na;
            result.rep_b = nb;
            result.rep_c = nc;
            result.tiling_used = true;
            result.rep_clamped = clamped;
            let total = (2 * na + 1) as i64 * (2 * nb + 1) as i64 * (2 * nc + 1) as i64 * reference.atoms.len() as i64;
            if total > 8_000_000 {
                result.message = format!("Tiling would test {total} atoms (limit 8M). Reduce the Wulff radius or use manual replication.");
                return result;
            }
            tile(reference, &orientation, center, (na, nb, nc), (a, b, c), element_colors, &inside)
        } else {
            result.tiling_used = false;
            untiled(reference, &orientation, center, element_colors, &inside)
        };
        if generated.is_empty() {
            result.message = "No atoms fell inside the generated Wulff polyhedron.".into();
            return result;
        }
        let generated = deduplicate(generated);
        result.estimated_diameter = 2.0 * preview.max_plane_distance;
        result.output_atoms = generated.len() as i32;
        result.wulff_family_count = params.wulff_planes.len() as i32;
        result.wulff_face_count = preview.faces.len() as i32;
        structure.atoms = generated;
        clear_metadata(structure);
        if params.set_output_cell {
            let pad = params.vacuum_padding;
            let (lo, hi) = (preview.min_point, preview.max_point);
            structure.has_unit_cell = true;
            structure.cell_offset = [(lo.x - pad) as f64, (lo.y - pad) as f64, (lo.z - pad) as f64];
            structure.cell_vectors = [
                [(hi.x - lo.x + 2.0 * pad) as f64, 0.0, 0.0],
                [0.0, (hi.y - lo.y + 2.0 * pad) as f64, 0.0],
                [0.0, 0.0, (hi.z - lo.z + 2.0 * pad) as f64],
            ];
        } else {
            structure.has_unit_cell = false;
        }
        result.success = true;
        result.message = format!("Wulff nanocrystal built: {} atoms, {} faces.", result.output_atoms, result.wulff_face_count);
        return result;
    }

    if params.shape == NanoShape::MeshModel {
        if model_vertices.is_empty() || model_indices.len() < 3 {
            result.message = "No 3D model loaded. Drop an OBJ/STL file in the builder first.".into();
            return result;
        }
        if params.model_scale <= 0.0 {
            result.message = "Model scale must be positive.".into();
            return result;
        }
    }
    if reference.atoms.is_empty() {
        result.message = "Reference structure has no atoms.".into();
        return result;
    }
    let center = build_clip_center(reference, params);
    let max_r = compute_bounding_radius(params);
    let orientation = match build_orientation_matrix(reference, params) {
        Ok(m) => m,
        Err(e) => {
            result.message = e;
            return result;
        }
    };
    let inside = |p: Vec3| is_inside_shape(p, params, model_vertices, model_indices);
    let generated = if reference.has_unit_cell {
        let a = glm::mul(&orientation, row_f32(reference, 0));
        let b = glm::mul(&orientation, row_f32(reference, 1));
        let c = glm::mul(&orientation, row_f32(reference, 2));
        let (la, lb, lc) = (glm::length(a), glm::length(b), glm::length(c));
        if la < 1e-8 || lb < 1e-8 || lc < 1e-8 {
            result.message = "Reference structure has degenerate lattice vectors.".into();
            return result;
        }
        let mut clamped = false;
        let (na, nb, nc) = if params.auto_replicate {
            let mut safe_rep = |len: f32| {
                let n = (max_r / len).ceil() as i32 + 2;
                if n > MAX_REPS {
                    clamped = true;
                    MAX_REPS
                } else {
                    n
                }
            };
            (safe_rep(la), safe_rep(lb), safe_rep(lc))
        } else {
            (params.rep_a.max(1), params.rep_b.max(1), params.rep_c.max(1))
        };
        result.rep_a = na;
        result.rep_b = nb;
        result.rep_c = nc;
        result.tiling_used = true;
        result.rep_clamped = clamped;
        let total = (2 * na + 1) as i64 * (2 * nb + 1) as i64 * (2 * nc + 1) as i64 * reference.atoms.len() as i64;
        if total > 8_000_000 {
            result.message = format!("Tiling would test {total} atoms (limit 8M). Reduce shape size or use manual replication.");
            return result;
        }
        tile(reference, &orientation, center, (na, nb, nc), (a, b, c), element_colors, &inside)
    } else {
        result.tiling_used = false;
        untiled(reference, &orientation, center, element_colors, &inside)
    };
    if generated.is_empty() {
        result.message = "No atoms within the specified shape. Try increasing the size parameter(s) or adjusting Model Scale.".into();
        return result;
    }
    let generated = deduplicate(generated);
    result.estimated_diameter = 2.0 * max_r;
    result.output_atoms = generated.len() as i32;
    structure.atoms = generated;
    clear_metadata(structure);
    if params.set_output_cell {
        let he = compute_shape_half_extents(params);
        let pad = params.vacuum_padding;
        structure.has_unit_cell = true;
        structure.cell_offset = [(center.x - he.hx - pad) as f64, (center.y - he.hy - pad) as f64, (center.z - he.hz - pad) as f64];
        structure.cell_vectors = [
            [2.0 * (he.hx + pad) as f64, 0.0, 0.0],
            [0.0, 2.0 * (he.hy + pad) as f64, 0.0],
            [0.0, 0.0, 2.0 * (he.hz + pad) as f64],
        ];
    } else {
        structure.has_unit_cell = false;
    }
    result.success = true;
    result.message = if params.shape == NanoShape::MeshModel {
        let mut m = format!("Custom structure built: {} atoms", result.output_atoms);
        if result.tiling_used {
            m += &format!(" (reps={}x{}x{})", result.rep_a, result.rep_b, result.rep_c);
        }
        m
    } else {
        "Nanocrystal built successfully.".into()
    };
    result
}
