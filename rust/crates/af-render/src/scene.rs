//! The scene the renderer draws (the C++ graphics/StructureInstanceBuilder,
//! the CPU half of graphics/SceneBuffers, app/SceneView and the camera
//! framing of app/EditorOps): atom instances with periodic boundary images,
//! the cell box, bonds, the default view and the frame matrices. All in
//! single precision with glm's formulas, so the numbers match the C++.

use std::collections::HashMap;

use af_core::display::{light_direction, DisplaySettings};
use af_core::glm;
use af_core::structure::Structure;
use glam::{IVec3, Mat3, Mat4, Vec3, Vec4};

const EPSILON: f32 = 1e-5;
const PBC_BOUNDARY_TOL: f32 = 1e-4;
const ATOM_VISUAL_SCALE: f32 = 0.60;
const MIN_BOND_DISTANCE: f32 = 0.10;

/// Per-instance data of the atoms to draw.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InstanceData {
    pub positions: Vec<Vec3>,
    pub colors: Vec<[f32; 3]>,
    pub scales: Vec<f32>,
    pub bond_radii: Vec<f32>,
    pub shininess: Vec<f32>,
    pub atomic_numbers: Vec<i32>,
    pub box_lines: Vec<Vec3>,
    /// Index of each instance in Structure::atoms.
    pub atom_indices: Vec<i32>,
    pub orbit_center: Vec3,
}

/// glm::radians for floats.
pub fn radians(degrees: f32) -> f32 {
    degrees * 0.017_453_292_519_943_295_f32
}

fn wrap_and_snap_fractional(mut value: f32) -> f32 {
    value -= value.floor();
    if value.abs() <= PBC_BOUNDARY_TOL || (1.0 - value).abs() <= PBC_BOUNDARY_TOL {
        return 0.0;
    }
    value
}

fn cell_vectors(s: &Structure) -> (Vec3, Vec3, Vec3, Vec3) {
    let v = |r: usize| Vec3::new(s.cell_vectors[r][0] as f32, s.cell_vectors[r][1] as f32, s.cell_vectors[r][2] as f32);
    let origin = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
    (origin, v(0), v(1), v(2))
}

/// origin + x a + y b + z c, in glm's order of operations.
fn combine(origin: Vec3, f: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    origin + f.x * a + f.y * b + f.z * c
}

fn append_pbc_boundary_images(s: &Structure, d: &mut InstanceData) {
    if !s.has_unit_cell || d.positions.is_empty() {
        return;
    }
    let (origin, a, b, c) = cell_vectors(s);
    let cell = Mat3::from_cols(a, b, c);
    if glm::determinant(&cell).abs() <= EPSILON {
        return;
    }
    let inverse = glm::inverse(&cell);
    let tol = if s.pbc_boundary_tol > 0.0 { s.pbc_boundary_tol } else { PBC_BOUNDARY_TOL };
    let base = d.positions.len();
    for i in 0..base {
        let mut frac = glm::mul(&inverse, d.positions[i] - origin);
        frac = Vec3::new(wrap_and_snap_fractional(frac.x), wrap_and_snap_fractional(frac.y), wrap_and_snap_fractional(frac.z));
        let canonical = combine(origin, frac, a, b, c);
        d.positions[i] = canonical;
        let shifts = |f: f32| if f.abs() <= tol { vec![0, 1] } else { vec![0] };
        for sx in shifts(frac.x) {
            for sy in shifts(frac.y) {
                for sz in shifts(frac.z) {
                    if sx == 0 && sy == 0 && sz == 0 {
                        continue;
                    }
                    let image = canonical + sx as f32 * a + sy as f32 * b + sz as f32 * c;
                    d.positions.push(image);
                    d.colors.push(d.colors[i]);
                    d.scales.push(d.scales[i]);
                    d.bond_radii.push(d.bond_radii[i]);
                    d.shininess.push(d.shininess[i]);
                    d.atomic_numbers.push(d.atomic_numbers[i]);
                    d.atom_indices.push(d.atom_indices[i]);
                }
            }
        }
    }
}

fn append_box_edges(corners: &[Vec3; 8], lines: &mut Vec<Vec3>) {
    const EDGES: [[usize; 2]; 12] = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]];
    for e in EDGES {
        lines.push(corners[e[0]]);
        lines.push(corners[e[1]]);
    }
}

/// The instances of a structure (optionally expanded by an integer
/// transform), with periodic boundary images and the box.
pub fn build_instance_data(s: &Structure, use_transform: bool, transform: &[[i32; 3]; 3], element_radii: &[f32], element_shininess: &[f32]) -> InstanceData {
    let mut d = InstanceData::default();
    for (i, atom) in s.atoms.iter().enumerate() {
        d.positions.push(Vec3::new(atom.x as f32, atom.y as f32, atom.z as f32));
        d.colors.push([atom.r, atom.g, atom.b]);
        let z = atom.atomic_number;
        let mut radius = 1.0;
        if z >= 0 && (z as usize) < element_radii.len() && element_radii[z as usize] > 0.0 {
            radius = element_radii[z as usize];
        }
        d.scales.push(radius * ATOM_VISUAL_SCALE);
        d.bond_radii.push(radius);
        let shine = if z >= 0 && (z as usize) < element_shininess.len() { element_shininess[z as usize] } else { 32.0 };
        d.shininess.push(shine);
        d.atomic_numbers.push(z);
        d.atom_indices.push(i as i32);
    }
    if use_transform && s.has_unit_cell {
        let base = std::mem::take(&mut d);
        let (origin, a, b, c) = cell_vectors(s);
        let inverse_cell = glm::inverse(&Mat3::from_cols(a, b, c));
        // transform[col][row] = matrix[row][col] in glm.
        let t = Mat3::from_cols(
            Vec3::new(transform[0][0] as f32, transform[1][0] as f32, transform[2][0] as f32),
            Vec3::new(transform[0][1] as f32, transform[1][1] as f32, transform[2][1] as f32),
            Vec3::new(transform[0][2] as f32, transform[1][2] as f32, transform[2][2] as f32),
        );
        if glm::determinant(&t).abs() <= EPSILON {
            d = base;
        } else {
            let inverse_t = glm::inverse(&t);
            let mut n_min = [0; 3];
            let mut n_max = [0; 3];
            for row in 0..3 {
                let (mut low, mut high) = (0, 0);
                for col in 0..3 {
                    let v = transform[row][col];
                    if v < 0 {
                        low += v;
                    }
                    if v > 0 {
                        high += v;
                    }
                }
                n_min[row] = low - 1;
                n_max[row] = high + 1;
            }
            for k in 0..base.positions.len() {
                let mut frac = glm::mul(&inverse_cell, base.positions[k] - origin);
                frac -= frac.floor();
                for ix in n_min[0]..=n_max[0] {
                    for iy in n_min[1]..=n_max[1] {
                        for iz in n_min[2]..=n_max[2] {
                            let shifted = frac + Vec3::new(ix as f32, iy as f32, iz as f32);
                            let mapped = glm::mul(&inverse_t, shifted);
                            if mapped.x < -EPSILON || mapped.y < -EPSILON || mapped.z < -EPSILON {
                                continue;
                            }
                            if mapped.x >= 1.0 - EPSILON || mapped.y >= 1.0 - EPSILON || mapped.z >= 1.0 - EPSILON {
                                continue;
                            }
                            d.positions.push(combine(origin, shifted, a, b, c));
                            d.colors.push(base.colors[k]);
                            d.scales.push(base.scales[k]);
                            d.bond_radii.push(base.bond_radii[k]);
                            d.shininess.push(base.shininess[k]);
                            d.atomic_numbers.push(base.atomic_numbers[k]);
                            d.atom_indices.push(base.atom_indices[k]);
                        }
                    }
                }
            }
        }
    }
    // For periodic display, boundary atoms are repeated on the opposite faces.
    if !use_transform && s.has_unit_cell {
        append_pbc_boundary_images(s, &mut d);
    }
    if d.positions.is_empty() {
        return d;
    }
    let (mut lo, mut hi) = (Vec3::splat(1e9), Vec3::splat(-1e9));
    for &p in &d.positions {
        lo = glm::min(lo, p);
        hi = glm::max(hi, p);
    }
    let corners: [Vec3; 8] = if s.has_unit_cell {
        let (origin, a, b, c) = cell_vectors(s);
        let (a, b, c) = if use_transform {
            let row = |r: usize| transform[r][0] as f32 * a + transform[r][1] as f32 * b + transform[r][2] as f32 * c;
            (row(0), row(1), row(2))
        } else {
            (a, b, c)
        };
        [origin, origin + a, origin + a + b, origin + b, origin + c, origin + a + c, origin + a + b + c, origin + b + c]
    } else {
        [
            Vec3::new(lo.x, lo.y, lo.z),
            Vec3::new(hi.x, lo.y, lo.z),
            Vec3::new(hi.x, hi.y, lo.z),
            Vec3::new(lo.x, hi.y, lo.z),
            Vec3::new(lo.x, lo.y, hi.z),
            Vec3::new(hi.x, lo.y, hi.z),
            Vec3::new(hi.x, hi.y, hi.z),
            Vec3::new(lo.x, hi.y, hi.z),
        ]
    };
    d.orbit_center = 0.5 * (corners[0] + corners[6]);
    append_box_edges(&corners, &mut d.box_lines);
    d
}

/// How bonds are drawn (Display Settings).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondStyle {
    /// Multiplies the automatic bond radius.
    pub thickness: f32,
    /// One colour instead of half of each atom's.
    pub uniform_colour: bool,
    pub colour: [f32; 3],
    /// Bonded below tolerance x (sum of element radii).
    pub tolerance: f32,
    /// Also bonds between two metals.
    pub metal_bonds: bool,
    /// No bond search from this many atoms.
    pub atom_limit: usize,
}

impl Default for BondStyle {
    fn default() -> Self {
        Self { thickness: 1.0, uniform_colour: false, colour: [0.62; 3], tolerance: 1.18, metal_bonds: false, atom_limit: 500000 }
    }
}

impl BondStyle {
    pub fn from_display(d: &DisplaySettings) -> Self {
        Self {
            thickness: d.bond_scale,
            uniform_colour: d.uniform_bond_colour,
            colour: d.bond_colour,
            tolerance: d.bond_tolerance,
            metal_bonds: d.metal_bonds,
            atom_limit: d.bond_atom_limit as usize,
        }
    }
}

/// Bonds as cylinders between two atoms, each half in its atom's colour.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Bonds {
    pub starts: Vec<Vec3>,
    pub ends: Vec<Vec3>,
    pub colors_a: Vec<[f32; 3]>,
    pub colors_b: Vec<[f32; 3]>,
    pub radii: Vec<f32>,
    pub shininess_a: Vec<f32>,
    pub shininess_b: Vec<f32>,
}

fn is_likely_covalent(z: i32) -> bool {
    matches!(z, 1 | 5 | 6 | 7 | 8 | 9 | 14 | 15 | 16 | 17 | 32 | 33 | 34 | 35 | 52 | 53)
}

fn grid_cell(p: Vec3, size: f32) -> IVec3 {
    IVec3::new((p.x / size).floor() as i32, (p.y / size).floor() as i32, (p.z / size).floor() as i32)
}

/// The bonds of a scene (SceneBuffers::upload): pairs within tolerance x the
/// sum of the element radii, from a spatial hash, in the C++ order.
pub fn find_bonds(d: &InstanceData, style: &BondStyle, element_filter: Option<&[bool; 119]>) -> Bonds {
    const MAX_BOND_COUNT: usize = 1_000_000;
    const CELL: f32 = 4.0;
    const MIN_RADIUS: f32 = 0.06;
    const MAX_RADIUS: f32 = 0.16;
    const INSET: f32 = 0.55;
    let mut out = Bonds::default();
    let n = d.positions.len();
    // CPU caches (and so bonds) only up to 100000 atoms.
    let positions: &[Vec3] = if n > 100_000 { &[] } else { &d.positions };
    if positions.len() < 2 || n >= style.atom_limit {
        return out;
    }
    let source_radii = if d.bond_radii.len() == n { &d.bond_radii } else { &d.scales };
    let mut grid: HashMap<IVec3, Vec<usize>> = HashMap::new();
    for (i, &p) in positions.iter().enumerate() {
        grid.entry(grid_cell(p, CELL)).or_default().push(i);
    }
    'atoms: for i in 0..positions.len() {
        let centre = grid_cell(positions[i], CELL);
        let radius_a = source_radii.get(i).copied().unwrap_or(1.0);
        let drawn_a = d.scales.get(i).map_or(radius_a, |&s| radius_a.min(s));
        for dx in -2..=2 {
            for dy in -2..=2 {
                for dz in -2..=2 {
                    let Some(list) = grid.get(&(centre + IVec3::new(dx, dy, dz))) else { continue };
                    for &j in list {
                        if j <= i {
                            continue;
                        }
                        let za = d.atomic_numbers.get(i).copied().unwrap_or(0);
                        let zb = d.atomic_numbers.get(j).copied().unwrap_or(0);
                        if let Some(mask) = element_filter {
                            let allowed = |z: i32| (1..=118).contains(&z) && mask[z as usize];
                            if !allowed(za) && !allowed(zb) {
                                continue;
                            }
                        } else if !style.metal_bonds && (1..=118).contains(&za) && (1..=118).contains(&zb) && !is_likely_covalent(za) && !is_likely_covalent(zb) {
                            continue;
                        }
                        let delta = positions[j] - positions[i];
                        let distance = glm::length(delta);
                        if distance <= MIN_BOND_DISTANCE {
                            continue;
                        }
                        let radius_b = source_radii.get(j).copied().unwrap_or(1.0);
                        if distance > (radius_a + radius_b) * style.tolerance {
                            continue;
                        }
                        let direction = delta / distance;
                        let drawn_b = d.scales.get(j).map_or(radius_b, |&s| radius_b.min(s));
                        let inset_a = (drawn_a * INSET).min(distance * 0.30);
                        let inset_b = (drawn_b * INSET).min(distance * 0.30);
                        let start = positions[i] + direction * inset_a;
                        let end = positions[j] - direction * inset_b;
                        if glm::length(end - start) <= MIN_BOND_DISTANCE {
                            continue;
                        }
                        let radius = style.thickness * (0.12 * (radius_a + radius_b) * 0.5).min(MAX_RADIUS).max(MIN_RADIUS);
                        out.starts.push(start);
                        out.ends.push(end);
                        if style.uniform_colour {
                            out.colors_a.push(style.colour);
                            out.colors_b.push(style.colour);
                        } else {
                            out.colors_a.push(d.colors.get(i).copied().unwrap_or([0.8; 3]));
                            out.colors_b.push(d.colors.get(j).copied().unwrap_or([0.8; 3]));
                        }
                        out.radii.push(radius);
                        out.shininess_a.push(d.shininess.get(i).copied().unwrap_or(32.0));
                        out.shininess_b.push(d.shininess.get(j).copied().unwrap_or(32.0));
                        if out.starts.len() >= MAX_BOND_COUNT {
                            break 'atoms;
                        }
                    }
                }
            }
        }
    }
    out
}

/// The orbit camera (the C++ Camera's view parameters).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrbitCamera {
    pub yaw: f32,
    pub pitch: f32,
    pub roll: f32,
    pub distance: f32,
    pub pan_offset: Vec3,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self { yaw: 45.0, pitch: 35.2643897, roll: 0.0, distance: 10.0, pan_offset: Vec3::ZERO }
    }
}

pub const MIN_DISTANCE: f32 = 2.0;
pub const MAX_DISTANCE: f32 = 1.0e6;

/// The atoms (or box) seen from the orbit centre: the largest distance plus radius, at least 1.
pub fn estimate_scene_radius(positions: &[Vec3], radii: &[f32], box_lines: &[Vec3], center: Vec3) -> f32 {
    let mut max_radius = 0.0f32;
    if !positions.is_empty() {
        for (i, &p) in positions.iter().enumerate() {
            let r = radii.get(i).copied().unwrap_or(0.0);
            max_radius = max_radius.max(glm::length(p - center) + r);
        }
    } else {
        for &p in box_lines {
            max_radius = max_radius.max(glm::length(p - center));
        }
    }
    max_radius.max(1.0)
}

/// The isometric default view, at a distance that frames the scene.
pub fn apply_default_view(camera: &mut OrbitCamera, scene: &Scene, width: i32, height: i32, fit: bool) {
    camera.yaw = 45.0;
    camera.pitch = 35.2643897;
    camera.roll = 0.0;
    camera.pan_offset = Vec3::ZERO;
    if !fit || scene.atom_count == 0 {
        camera.distance = 10.0;
        return;
    }
    let max_radius = scene.radius();
    let aspect = if height > 0 { width as f32 / height as f32 } else { 1.0 };
    let vertical = radians(45.0);
    let horizontal = 2.0 * ((vertical * 0.5).tan() * aspect).atan();
    let half = (0.5 * vertical.min(horizontal)).max(radians(10.0));
    let framed = (max_radius / half.sin()) * 1.15;
    camera.distance = MIN_DISTANCE.max(MAX_DISTANCE.min(framed));
}

/// The data the renderer draws (what SceneBuffers holds on the CPU).
#[derive(Clone, Debug, Default)]
pub struct Scene {
    pub instances: InstanceData,
    pub bonds: Bonds,
    pub atom_count: usize,
    pub orbit_center: Vec3,
    /// Selected instances, drawn with a halo (the atoms keep their colours).
    pub selected: Vec<usize>,
}

impl Scene {
    /// SceneBuffers::upload: the instances and their bonds.
    pub fn new(instances: InstanceData, style: &BondStyle) -> Self {
        let bonds = find_bonds(&instances, style, None);
        let atom_count = instances.positions.len();
        let orbit_center = if atom_count > 0 { instances.orbit_center } else { Vec3::ZERO };
        Self { instances, bonds, atom_count, orbit_center, selected: Vec::new() }
    }

    pub fn radius(&self) -> f32 {
        let (p, r): (&[Vec3], &[f32]) =
            if self.atom_count > 100_000 { (&[], &[]) } else { (&self.instances.positions, &self.instances.scales) };
        estimate_scene_radius(p, r, &self.instances.box_lines, self.orbit_center)
    }
}

/// glm::perspective (right-handed, clip depth -1..1).
pub fn perspective(fovy: f32, aspect: f32, near: f32, far: f32) -> Mat4 {
    let t = (fovy / 2.0).tan();
    let mut m = Mat4::ZERO;
    m.x_axis.x = 1.0 / (aspect * t);
    m.y_axis.y = 1.0 / t;
    m.z_axis.z = -(far + near) / (far - near);
    m.z_axis.w = -1.0;
    m.w_axis.z = -(2.0 * far * near) / (far - near);
    m
}

/// glm::ortho (right-handed, clip depth -1..1).
pub fn ortho(left: f32, right: f32, bottom: f32, top: f32, near: f32, far: f32) -> Mat4 {
    let mut m = Mat4::IDENTITY;
    m.x_axis.x = 2.0 / (right - left);
    m.y_axis.y = 2.0 / (top - bottom);
    m.z_axis.z = -2.0 / (far - near);
    m.w_axis.x = -(right + left) / (right - left);
    m.w_axis.y = -(top + bottom) / (top - bottom);
    m.w_axis.z = -(far + near) / (far - near);
    m
}

/// glm::lookAt (right-handed).
pub fn look_at(eye: Vec3, center: Vec3, up: Vec3) -> Mat4 {
    let f = glm::normalize(center - eye);
    let s = glm::normalize(glm::cross(f, up));
    let u = glm::cross(s, f);
    let mut m = Mat4::IDENTITY;
    m.x_axis.x = s.x;
    m.y_axis.x = s.y;
    m.z_axis.x = s.z;
    m.x_axis.y = u.x;
    m.y_axis.y = u.y;
    m.z_axis.y = u.z;
    m.x_axis.z = -f.x;
    m.y_axis.z = -f.y;
    m.z_axis.z = -f.z;
    m.w_axis.x = -glm::dot(s, eye);
    m.w_axis.y = -glm::dot(u, eye);
    m.w_axis.z = glm::dot(f, eye);
    m
}

/// glm's mat4 product: each column a left-to-right sum of scaled columns.
pub fn mat_mul(a: &Mat4, b: &Mat4) -> Mat4 {
    let column = |c: Vec4| a.x_axis * c.x + a.y_axis * c.y + a.z_axis * c.z + a.w_axis * c.w;
    Mat4::from_cols(column(b.x_axis), column(b.y_axis), column(b.z_axis), column(b.w_axis))
}

/// The matrices and lighting of one frame (FrameView).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FrameView {
    pub width: i32,
    pub height: i32,
    pub projection: Mat4,
    pub view: Mat4,
    pub light_mvp: Mat4,
    pub light_position: Vec3,
    pub camera_position: Vec3,
    /// Camera distances of the nearest and farthest scene points, for depth cueing.
    pub depth_near: f32,
    pub depth_far: f32,
}

const LIGHT_DISTANCE: f32 = 69.282_032_3;

/// buildFrameView: projection, view and the light for a camera.
pub fn build_frame_view(camera: &OrbitCamera, scene: &Scene, orthographic: bool, width: i32, height: i32, display: &DisplaySettings) -> FrameView {
    let aspect = width as f32 / height as f32;
    let fov = radians(45.0);
    let scene_radius = scene.radius();
    let padding = 10.0f32.max(scene_radius * 0.25);
    let projection = if orthographic {
        let half_height = 0.1f32.max(camera.distance * (fov * 0.5).tan());
        let half_width = half_height * aspect;
        let range = 1000.0f32.max(camera.distance + scene_radius + padding);
        ortho(-half_width, half_width, -half_height, half_height, -range, range)
    } else {
        let nearest = camera.distance - scene_radius;
        let near = if nearest > 0.0 { 0.01f32.max(nearest * 0.25) } else { 0.01 };
        let far = (near + 100.0).max(camera.distance + scene_radius + padding);
        perspective(fov, aspect, near, far)
    };
    let (yaw, pitch) = (radians(camera.yaw), radians(camera.pitch));
    let offset = Vec3::new(camera.distance * pitch.cos() * yaw.sin(), camera.distance * pitch.sin(), camera.distance * pitch.cos() * yaw.cos());
    let target = scene.orbit_center + camera.pan_offset;
    let camera_position = target + offset;
    let base_up = glm::normalize(Vec3::new(-pitch.sin() * yaw.sin(), pitch.cos(), -pitch.sin() * yaw.cos()));
    let forward = glm::normalize(target - camera_position);
    let right = glm::normalize(glm::cross(forward, base_up));
    let roll = radians(camera.roll);
    let rolled_up = roll.cos() * base_up + roll.sin() * right;
    let view = look_at(camera_position, target, rolled_up);
    let light_projection = perspective(radians(45.0), 1.0, 0.1, 1000.0);
    let mut to_light = Vec3::from(light_direction(display));
    if display.light_follows_camera {
        to_light = to_light.x * right + to_light.y * rolled_up - to_light.z * forward;
    }
    let to_light = glm::normalize(to_light);
    let light_position = scene.orbit_center + LIGHT_DISTANCE * to_light;
    let centre_distance = glm::length(camera_position - scene.orbit_center);
    let depth_near = 0.0f32.max(centre_distance - scene_radius);
    let depth_far = (depth_near + 1e-3).max(centre_distance + scene_radius);
    let light_up = if to_light.y.abs() > 0.98 { Vec3::Z } else { Vec3::Y };
    let light_mvp = mat_mul(&light_projection, &look_at(light_position, scene.orbit_center, light_up));
    FrameView { width, height, projection, view, light_mvp, light_position, camera_position, depth_near, depth_far }
}

/// The unit sphere (stacks x slices quads, CCW from outside): vertices and indices.
pub fn sphere_mesh(stacks: i32, slices: i32) -> (Vec<[f32; 3]>, Vec<u32>) {
    const PI: f32 = std::f32::consts::PI;
    let mut vertices = Vec::new();
    for i in 0..=stacks {
        let lat = PI * (-0.5 + i as f32 / stacks as f32);
        let (cos_lat, sin_lat) = (lat.cos(), lat.sin());
        for j in 0..=slices {
            let lng = 2.0 * PI * j as f32 / slices as f32;
            vertices.push([lng.cos() * cos_lat, lng.sin() * cos_lat, sin_lat]);
        }
    }
    let mut indices = Vec::new();
    for i in 0..stacks {
        for j in 0..slices {
            let row0 = (i * (slices + 1)) as u32;
            let row1 = ((i + 1) * (slices + 1)) as u32;
            let (v00, v01, v10, v11) = (row0 + j as u32, row0 + j as u32 + 1, row1 + j as u32, row1 + j as u32 + 1);
            indices.extend([v00, v10, v01, v01, v10, v11]);
        }
    }
    (vertices, indices)
}

/// The unit cylinder along z (-0.5..0.5) with caps, as triangles.
pub fn cylinder_mesh(slices: i32) -> Vec<[f32; 3]> {
    let two_pi = 6.283_185_307_179_586_476_92_f32;
    let mut v = Vec::new();
    for index in 0..slices {
        let a0 = two_pi * index as f32 / slices as f32;
        let a1 = two_pi * (index + 1) as f32 / slices as f32;
        let (x0, y0, x1, y1) = (a0.cos(), a0.sin(), a1.cos(), a1.sin());
        v.extend([[x0, y0, -0.5], [x0, y0, 0.5], [x1, y1, 0.5]]);
        v.extend([[x0, y0, -0.5], [x1, y1, 0.5], [x1, y1, -0.5]]);
        v.extend([[0.0, 0.0, -0.5], [x1, y1, -0.5], [x0, y0, -0.5]]);
        v.extend([[0.0, 0.0, 0.5], [x0, y0, 0.5], [x1, y1, 0.5]]);
    }
    v
}

/// glm::degrees for floats.
pub fn degrees(radians: f32) -> f32 {
    radians * 57.295_779_513_082_32_f32
}

/// glm::rotate(mat4(1), angle, axis) as a 3x3 matrix (columns).
fn glm_rotation(angle: f32, v: Vec3) -> Mat3 {
    let (c, s) = (angle.cos(), angle.sin());
    let axis = glm::normalize(v);
    let temp = axis * (1.0 - c);
    Mat3::from_cols(
        Vec3::new(c + temp.x * axis.x, temp.x * axis.y + s * axis.z, temp.x * axis.z - s * axis.y),
        Vec3::new(temp.y * axis.x - s * axis.z, c + temp.y * axis.y, temp.y * axis.z + s * axis.x),
        Vec3::new(temp.z * axis.x + s * axis.y, temp.z * axis.y - s * axis.x, c + temp.z * axis.z),
    )
}

/// The toolbar's Rotate: turning the crystal by `angle` degrees about x, y
/// or z equals orbiting the camera by -angle (rotateCrystalAroundAxis); the
/// camera keeps its roll relative to the rotated frame.
pub fn rotate_crystal_view(camera: &mut OrbitCamera, axis: usize, angle: f64) {
    let rotation = glm_rotation(radians(-angle as f32), [Vec3::X, Vec3::Y, Vec3::Z][axis.min(2)]);
    let (yaw, pitch) = (radians(camera.yaw), radians(camera.pitch));
    let forward = Vec3::new(pitch.cos() * yaw.sin(), pitch.sin(), pitch.cos() * yaw.cos());
    let canonical_up = Vec3::new(-pitch.sin() * yaw.sin(), pitch.cos(), -pitch.sin() * yaw.cos());
    let right = glm::normalize(glm::cross(forward, canonical_up));
    let roll = radians(camera.roll);
    let up = canonical_up * roll.cos() + right * roll.sin();
    let new_forward = glm::normalize(glm::mul(&rotation, forward));
    let new_up = glm::normalize(glm::mul(&rotation, up));
    let new_pitch = degrees(new_forward.y.clamp(-1.0, 1.0).asin());
    let new_yaw = degrees(new_forward.x.atan2(new_forward.z));
    let (ny, np) = (radians(new_yaw), radians(new_pitch));
    let new_canonical_up = Vec3::new(-np.sin() * ny.sin(), np.cos(), -np.sin() * ny.cos());
    let sin_roll = glm::dot(glm::cross(new_canonical_up, new_up), new_forward);
    let cos_roll = glm::dot(new_canonical_up, new_up);
    camera.yaw = new_yaw;
    camera.pitch = new_pitch;
    camera.roll = degrees(sin_roll.atan2(cos_roll));
}
