//! Dislocation insertion (the C++ algorithms/DislocationBuilder): the lattice
//! family, crystallographic directions, a local frame around the line, the
//! isotropic (or anisotropic Stroh) displacement field with shape and cutoff
//! masks, an optional dipole partner, centre-of-mass correction, the
//! displaced region's outline and a validation report.

use af_core::format::format_g;
use af_core::glm::{self, cross, dot, length, mul};
use af_core::Structure;
use glam::{DMat3, DVec3, IVec3, Mat3, Vec2, Vec3};

use crate::anisotropic::{anisotropic_displacement, solve_stroh, StrohSolution};
use crate::elastic::{make_cubic_elastic_tensor, make_hexagonal_elastic_tensor, make_voigt_elastic_tensor, rotate_elastic_tensor};
use crate::stacking_fault::{detect_stacking_fault_family, StackingFaultFamily};

const PI: f32 = std::f32::consts::PI;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DislocationLatticeFamily {
    #[default]
    Unknown = 0,
    Fcc,
    Hcp,
    Bcc,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DislocationCharacter {
    #[default]
    Edge = 0,
    Screw,
    Mixed,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DislocationShape {
    #[default]
    HalfPlane = 0,
    Cylinder,
    Sphere,
    Ellipsoid,
    Freeform2D,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ElasticSymmetry {
    #[default]
    Cubic,
    Hexagonal,
    General,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DislocationDetectionResult {
    pub success: bool,
    pub message: String,
    pub family: DislocationLatticeFamily,
    pub fcc_count: i32,
    pub hcp_count: i32,
    pub bcc_count: i32,
    pub recognized_count: i32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DislocationParams {
    pub character: DislocationCharacter,
    pub shape: DislocationShape,
    pub auto_directions: bool,
    pub plane_hkl: IVec3,
    pub burgers_uvw: IVec3,
    pub line_uvw: IVec3,
    pub use_fractional_line_point: bool,
    pub line_point_fractional: Vec3,
    pub line_point_cartesian_offset: Vec3,
    pub burgers_scale: f32,
    /// If > 0, |b| in A instead of the detected value.
    pub burgers_override_magnitude: f32,
    pub mixed_character_angle_deg: f32,
    pub poisson_ratio: f32,
    pub core_radius: f32,
    pub cutoff_radius: f32,
    pub line_half_length: f32,
    /// Anisotropic elasticity (Stroh), constants in the cell's crystal axes (GPa).
    pub anisotropic_elasticity: bool,
    pub elastic_symmetry: ElasticSymmetry,
    pub elastic_c11: f64,
    pub elastic_c12: f64,
    pub elastic_c13: f64,
    pub elastic_c33: f64,
    pub elastic_c44: f64,
    pub elastic_voigt: [[f64; 6]; 6],
    pub elastic_noise_amplitude: f64,
    pub elastic_noise_seed: u32,
    /// A partner with the opposite Burgers vector at dipole_offset (e1, e2).
    pub dipole: bool,
    pub dipole_offset: Vec2,
    pub cylinder_radius: f32,
    pub sphere_radius: f32,
    pub ellipsoid_radii: Vec3,
    pub freeform_points: Vec<Vec2>,
}

impl Default for DislocationParams {
    fn default() -> Self {
        Self {
            character: DislocationCharacter::Edge,
            shape: DislocationShape::HalfPlane,
            auto_directions: true,
            plane_hkl: IVec3::new(1, 1, 1),
            burgers_uvw: IVec3::new(1, -1, 0),
            line_uvw: IVec3::new(1, 1, -2),
            use_fractional_line_point: true,
            line_point_fractional: Vec3::splat(0.5),
            line_point_cartesian_offset: Vec3::ZERO,
            burgers_scale: 1.0,
            burgers_override_magnitude: -1.0,
            mixed_character_angle_deg: 45.0,
            poisson_ratio: 0.33,
            core_radius: 1.2,
            cutoff_radius: 20.0,
            line_half_length: 1_000_000.0,
            anisotropic_elasticity: false,
            elastic_symmetry: ElasticSymmetry::Cubic,
            elastic_c11: 0.0,
            elastic_c12: 0.0,
            elastic_c13: 0.0,
            elastic_c33: 0.0,
            elastic_c44: 0.0,
            elastic_voigt: [[0.0; 6]; 6],
            elastic_noise_amplitude: 1e-4,
            elastic_noise_seed: 1,
            dipole: false,
            dipole_offset: Vec2::new(15.0, 0.0),
            cylinder_radius: 12.0,
            sphere_radius: 12.0,
            ellipsoid_radii: Vec3::new(14.0, 10.0, 10.0),
            freeform_points: vec![Vec2::new(-8.0, 0.0), Vec2::new(8.0, 0.0), Vec2::new(8.0, 8.0), Vec2::new(-8.0, 8.0)],
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DislocationValidationReport {
    pub passed: bool,
    pub finite_coordinates: bool,
    pub lattice_family_preserved: bool,
    pub min_distance_safe: bool,
    pub family_before: DislocationLatticeFamily,
    pub family_after: DislocationLatticeFamily,
    pub min_interatomic_distance: f32,
    pub rms_displacement: f32,
    pub max_displacement: f32,
    pub message: String,
}

impl Default for DislocationValidationReport {
    fn default() -> Self {
        Self {
            passed: false,
            finite_coordinates: true,
            lattice_family_preserved: false,
            min_distance_safe: true,
            family_before: DislocationLatticeFamily::Unknown,
            family_after: DislocationLatticeFamily::Unknown,
            min_interatomic_distance: 0.0,
            rms_displacement: 0.0,
            max_displacement: 0.0,
            message: String::new(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DislocationResult {
    pub success: bool,
    pub message: String,
    pub output: Structure,
    pub family: DislocationLatticeFamily,
    pub shifted_atom_count: i32,
    pub burgers_magnitude: f32,
    pub line_point: Vec3,
    pub line_direction: Vec3,
    pub burgers_direction: Vec3,
    pub e1: Vec3,
    pub e2: Vec3,
    /// The displaced region's outline in the slip plane.
    pub loop_points: Vec<Vec3>,
    pub validation: DislocationValidationReport,
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

pub fn dislocation_lattice_family_name(family: DislocationLatticeFamily) -> &'static str {
    match family {
        DislocationLatticeFamily::Fcc => "FCC-like",
        DislocationLatticeFamily::Hcp => "HCP-like",
        DislocationLatticeFamily::Bcc => "BCC-like",
        DislocationLatticeFamily::Unknown => "Unknown",
    }
}

fn safe_normalize(v: Vec3, fallback: Vec3) -> Vec3 {
    let n = length(v);
    if n <= 1e-6 {
        fallback
    } else {
        v / n
    }
}

fn direction_from_uvw(cell: &Mat3, uvw: IVec3) -> Vec3 {
    safe_normalize(cell.x_axis * uvw.x as f32 + cell.y_axis * uvw.y as f32 + cell.z_axis * uvw.z as f32, Vec3::X)
}

fn plane_normal_from_hkl(cell: &Mat3, hkl: IVec3) -> Vec3 {
    let volume = dot(cell.x_axis, cross(cell.y_axis, cell.z_axis));
    if volume.abs() <= 1e-8 {
        return Vec3::Z;
    }
    let a = cross(cell.y_axis, cell.z_axis) / volume;
    let b = cross(cell.z_axis, cell.x_axis) / volume;
    let c = cross(cell.x_axis, cell.y_axis) / volume;
    safe_normalize(a * hkl.x as f32 + b * hkl.y as f32 + c * hkl.z as f32, Vec3::Z)
}

fn auto_lattice_vectors(family: DislocationLatticeFamily) -> (IVec3, IVec3, IVec3) {
    let v = IVec3::new;
    match family {
        DislocationLatticeFamily::Fcc => (v(1, 1, 1), v(1, -1, 0), v(1, 1, -2)),
        DislocationLatticeFamily::Hcp => (v(0, 0, 1), v(1, -1, 0), v(1, 1, 0)),
        DislocationLatticeFamily::Bcc => (v(1, 1, 0), v(1, 1, 1), v(1, -1, 0)),
        DislocationLatticeFamily::Unknown => (v(0, 0, 1), v(1, 0, 0), v(0, 1, 0)),
    }
}

/// Approximate nearest-neighbour distance from a 2.5 A grid.
fn nearest_distance_approx(structure: &Structure) -> f32 {
    if structure.atoms.len() < 2 {
        return 0.0;
    }
    let mut grid: std::collections::HashMap<(i32, i32, i32), Vec<usize>> = std::collections::HashMap::new();
    let mut best = f32::MAX;
    let key = |p: Vec3| ((p.x / 2.5).floor() as i32, (p.y / 2.5).floor() as i32, (p.z / 2.5).floor() as i32);
    let pos = |i: usize| Vec3::new(structure.atoms[i].x as f32, structure.atoms[i].y as f32, structure.atoms[i].z as f32);
    for i in 0..structure.atoms.len() {
        let p = pos(i);
        let k = key(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(list) = grid.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                        for &j in list {
                            best = min_f(best, length(p - pos(j)));
                        }
                    }
                }
            }
        }
        grid.entry(k).or_default().push(i);
    }
    if best.is_finite() {
        best
    } else {
        0.0
    }
}

fn estimate_lattice_parameter(structure: &Structure, family: DislocationLatticeFamily) -> f32 {
    let dnn = nearest_distance_approx(structure);
    if dnn > 1e-4 {
        return match family {
            DislocationLatticeFamily::Fcc => max_f(1e-3, dnn * 2.0f32.sqrt()),
            DislocationLatticeFamily::Bcc => max_f(1e-3, dnn * (2.0 / 3.0f32.sqrt())),
            _ => max_f(1e-3, dnn),
        };
    }
    if !structure.has_unit_cell {
        return 1.0;
    }
    let c = glm::cell_matrix(structure);
    max_f(1e-3, (length(c.x_axis) + length(c.y_axis) + length(c.z_axis)) / 3.0)
}

fn base_burgers_magnitude(family: DislocationLatticeFamily, a: f32) -> f32 {
    match family {
        DislocationLatticeFamily::Fcc => a / 2.0f32.sqrt(),
        DislocationLatticeFamily::Bcc => a * 3.0f32.sqrt() * 0.5,
        _ => a,
    }
}

/// Ray-casting point-in-polygon.
fn point_in_polygon(point: Vec2, polygon: &[Vec2]) -> bool {
    if polygon.len() < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = polygon.len() - 1;
    for i in 0..polygon.len() {
        let (pi, pj) = (polygon[i], polygon[j]);
        let crosses = (pi.y > point.y) != (pj.y > point.y) && point.x < (pj.x - pi.x) * (point.y - pi.y) / max_f(1e-8, pj.y - pi.y) + pi.x;
        if crosses {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn shape_mask(p: &DislocationParams, x: f32, y: f32, z: f32) -> f32 {
    let inside = match p.shape {
        DislocationShape::HalfPlane => y >= 0.0,
        DislocationShape::Cylinder => x * x + y * y <= p.cylinder_radius * p.cylinder_radius,
        DislocationShape::Sphere => x * x + y * y + z * z <= p.sphere_radius * p.sphere_radius,
        DislocationShape::Ellipsoid => {
            let rx = max_f(1e-4, p.ellipsoid_radii.x);
            let ry = max_f(1e-4, p.ellipsoid_radii.y);
            let rz = max_f(1e-4, p.ellipsoid_radii.z);
            (x * x) / (rx * rx) + (y * y) / (ry * ry) + (z * z) / (rz * rz) <= 1.0
        }
        DislocationShape::Freeform2D => point_in_polygon(Vec2::new(x, y), &p.freeform_points),
    };
    if inside {
        1.0
    } else {
        0.0
    }
}

fn attenuation_mask(p: &DislocationParams, radial: f32) -> f32 {
    if p.cutoff_radius <= 0.0 {
        return 1.0;
    }
    let r = radial / max_f(1e-4, p.cutoff_radius);
    let t = r * r;
    (-(t * t)).exp()
}

fn structure_center(structure: &Structure) -> Vec3 {
    if structure.atoms.is_empty() {
        return Vec3::ZERO;
    }
    let mut center = Vec3::ZERO;
    for a in &structure.atoms {
        center += Vec3::new(a.x as f32, a.y as f32, a.z as f32);
    }
    center / structure.atoms.len() as f32
}

struct Kernel {
    line_point: Vec3,
    line_dir: Vec3,
    e1: Vec3,
    e2: Vec3,
    pref_edge: f32,
    pref_screw: f32,
    poisson: f32,
    core: f32,
    line_half_length: f32,
    stroh: Option<StrohSolution>,
    dipole_line_point: Option<Vec3>,
}

/// One dislocation's displacement at `position` (`sign` -1 for a dipole partner).
fn one_displacement(p: &DislocationParams, k: &Kernel, position: Vec3, line_point: Vec3, sign: f32) -> Vec3 {
    let rel = position - line_point;
    let z = dot(rel, k.line_dir);
    if z.abs() > k.line_half_length {
        return Vec3::ZERO;
    }
    let x = dot(rel, k.e1);
    let y = dot(rel, k.e2);
    let r2 = x * x + y * y + k.core * k.core;
    let sh = shape_mask(p, x, y, z);
    if sh <= 0.0 {
        return Vec3::ZERO;
    }
    let weight = sh * attenuation_mask(p, (x * x + y * y).sqrt());
    if weight <= 1e-6 {
        return Vec3::ZERO;
    }
    let (mut ux, mut uy, mut uz) = (0.0f32, 0.0f32, 0.0f32);
    if let Some(stroh) = &k.stroh {
        // Clamp to the core radius along the same ray (a fixed ray on the line itself).
        let radial = (x * x + y * y).sqrt();
        let (mut ex, mut ey) = (x as f64, y as f64);
        if radial < k.core {
            if radial < 1e-6 {
                ex = k.core as f64;
                ey = 0.0;
            } else {
                ex = (x / radial * k.core) as f64;
                ey = (y / radial * k.core) as f64;
            }
        }
        let u = anisotropic_displacement(stroh, ex, ey);
        return sign * weight * (u.x as f32 * k.e1 + u.y as f32 * k.e2 + u.z as f32 * k.line_dir);
    }
    if k.pref_edge.abs() > 0.0 {
        let atan_term = y.atan2(x);
        let inv = 1.0 / max_f(1e-7, 2.0 * (1.0 - k.poisson) * r2);
        ux = k.pref_edge * (atan_term + (x * y) * inv);
        let log_term = max_f(1e-7, r2).ln();
        let q = max_f(1e-7, 4.0 * (1.0 - k.poisson));
        uy = -k.pref_edge * (((1.0 - 2.0 * k.poisson) / q) * log_term + (x * x - y * y) / (q * max_f(1e-7, r2)));
    }
    if k.pref_screw.abs() > 0.0 {
        uz = k.pref_screw * y.atan2(x);
    }
    sign * weight * (ux * k.e1 + uy * k.e2 + uz * k.line_dir)
}

fn displacement(p: &DislocationParams, k: &Kernel, position: Vec3) -> Vec3 {
    let mut total = one_displacement(p, k, position, k.line_point, 1.0);
    if let Some(partner) = k.dipole_line_point {
        total += one_displacement(p, k, position, partner, -1.0);
    }
    total
}

pub fn detect_dislocation_lattice(structure: &Structure, use_pbc_for_detection: bool) -> DislocationDetectionResult {
    let sf = detect_stacking_fault_family(structure, use_pbc_for_detection);
    DislocationDetectionResult {
        success: sf.success,
        message: sf.message,
        family: match sf.family {
            StackingFaultFamily::Fcc => DislocationLatticeFamily::Fcc,
            StackingFaultFamily::Hcp => DislocationLatticeFamily::Hcp,
            StackingFaultFamily::Bcc => DislocationLatticeFamily::Bcc,
            StackingFaultFamily::Unknown => DislocationLatticeFamily::Unknown,
        },
        fcc_count: sf.fcc_count,
        hcp_count: sf.hcp_count,
        bcc_count: sf.bcc_count,
        recognized_count: sf.recognized_count,
    }
}

/// `ostream << float` (6 significant digits).
fn g6(x: f32) -> String {
    format_g(x as f64, 6)
}

fn validation_report(before: &Structure, after: &Structure, family_before: DislocationLatticeFamily) -> DislocationValidationReport {
    let mut report = DislocationValidationReport { family_before, finite_coordinates: true, ..Default::default() };
    report.finite_coordinates = after.atoms.iter().all(|a| a.x.is_finite() && a.y.is_finite() && a.z.is_finite());
    report.min_interatomic_distance = nearest_distance_approx(after);
    report.min_distance_safe = report.min_interatomic_distance >= 0.20;
    // Without periodic images: a single dislocation's field is not periodic.
    report.family_after = detect_dislocation_lattice(after, false).family;
    report.lattice_family_preserved = if family_before == DislocationLatticeFamily::Unknown {
        report.family_after != DislocationLatticeFamily::Unknown
    } else {
        family_before == report.family_after
    };
    let (mut sum_sq, mut max_disp) = (0.0f32, 0.0f32);
    let n = before.atoms.len().min(after.atoms.len());
    for i in 0..n {
        let p0 = Vec3::new(before.atoms[i].x as f32, before.atoms[i].y as f32, before.atoms[i].z as f32);
        let p1 = Vec3::new(after.atoms[i].x as f32, after.atoms[i].y as f32, after.atoms[i].z as f32);
        let d = length(p1 - p0);
        sum_sq += d * d;
        max_disp = max_f(max_disp, d);
    }
    report.rms_displacement = if n > 0 { (sum_sq / n as f32).sqrt() } else { 0.0 };
    report.max_displacement = max_disp;
    report.passed = report.finite_coordinates && report.min_distance_safe && report.lattice_family_preserved;
    report.message = if report.passed {
        format!("Validation passed. Family preserved, coordinates finite, and min distance = {} A.", g6(report.min_interatomic_distance))
    } else {
        let mut parts = Vec::new();
        if !report.finite_coordinates {
            parts.push("non-finite coordinates detected".to_string());
        }
        if !report.min_distance_safe {
            parts.push(format!("very short interatomic distance ({} A)", g6(report.min_interatomic_distance)));
        }
        if !report.lattice_family_preserved {
            parts.push(format!(
                "lattice family changed from {} to {}",
                dislocation_lattice_family_name(report.family_before),
                dislocation_lattice_family_name(report.family_after)
            ));
        }
        format!("Validation warning: {}.", parts.join("; "))
    };
    report
}

pub fn build_dislocation(base: &Structure, params: &DislocationParams) -> DislocationResult {
    let mut result = DislocationResult { line_direction: Vec3::Z, burgers_direction: Vec3::X, e1: Vec3::X, e2: Vec3::Y, ..Default::default() };
    if base.atoms.is_empty() {
        result.message = "Input structure has no atoms.".into();
        return result;
    }
    if !base.has_unit_cell {
        result.message = "Dislocation insertion requires a structure with a unit cell.".into();
        return result;
    }
    let detection = detect_dislocation_lattice(base, true);
    result.family = detection.family;
    if params.auto_directions && detection.family == DislocationLatticeFamily::Unknown {
        result.message = "Could not identify lattice family (FCC/HCP/BCC). Disable auto vectors and set hkl/uvw manually.".into();
        return result;
    }
    let Some((cell, _)) = glm::cell_matrices(base) else {
        result.message = "Could not build cell matrices from the structure.".into();
        return result;
    };
    let (plane_hkl, burgers_uvw, line_uvw) =
        if params.auto_directions { auto_lattice_vectors(detection.family) } else { (params.plane_hkl, params.burgers_uvw, params.line_uvw) };
    let plane_normal = plane_normal_from_hkl(&cell, plane_hkl);
    let burgers_dir = direction_from_uvw(&cell, burgers_uvw);
    let screw = params.character == DislocationCharacter::Screw;
    let line_dir = if screw {
        burgers_dir
    } else if params.auto_directions {
        safe_normalize(cross(plane_normal, burgers_dir), direction_from_uvw(&cell, line_uvw))
    } else {
        direction_from_uvw(&cell, line_uvw)
    };
    // A screw's e2 is the plane normal made perpendicular to the line (Gram-Schmidt).
    let e2 = if screw {
        let mut reference = plane_normal;
        if length(cross(line_dir, reference)) <= 1e-6 {
            reference = if dot(line_dir, Vec3::Z).abs() < 0.9 { Vec3::Z } else { Vec3::Y };
        }
        safe_normalize(reference - line_dir * dot(line_dir, reference), Vec3::X)
    } else {
        safe_normalize(cross(line_dir, burgers_dir), plane_normal)
    };
    if !screw && length(cross(line_dir, burgers_dir)) <= 1e-6 {
        result.message = "For edge/mixed dislocations, line direction must not be parallel to Burgers direction.".into();
        return result;
    }
    let e1 = safe_normalize(cross(e2, line_dir), burgers_dir);
    let lattice_a = estimate_lattice_parameter(base, detection.family);
    let bmag = if params.burgers_override_magnitude > 0.0 {
        max_f(1e-5, params.burgers_scale * params.burgers_override_magnitude)
    } else {
        max_f(1e-5, params.burgers_scale * base_burgers_magnitude(detection.family, lattice_a))
    };
    let poisson = if params.poisson_ratio < 0.02 {
        0.02
    } else if 0.49 < params.poisson_ratio {
        0.49
    } else {
        params.poisson_ratio
    };
    let core = max_f(1e-4, params.core_radius);
    let (edge_fraction, screw_fraction) = match params.character {
        DislocationCharacter::Edge => (1.0f32, 0.0f32),
        DislocationCharacter::Screw => (0.0, 1.0),
        DislocationCharacter::Mixed => {
            let angle = params.mixed_character_angle_deg * (PI / 180.0);
            (angle.cos(), angle.sin())
        }
    };
    let line_point = if params.use_fractional_line_point && base.has_unit_cell {
        let origin = Vec3::new(base.cell_offset[0] as f32, base.cell_offset[1] as f32, base.cell_offset[2] as f32);
        origin + mul(&cell, params.line_point_fractional) + params.line_point_cartesian_offset
    } else {
        structure_center(base) + params.line_point_cartesian_offset
    };
    let mut kernel = Kernel {
        line_point,
        line_dir,
        e1,
        e2,
        pref_edge: (edge_fraction * bmag) / (2.0 * PI),
        pref_screw: (screw_fraction * bmag) / (2.0 * PI),
        poisson,
        core,
        line_half_length: params.line_half_length,
        stroh: None,
        dipole_line_point: params.dipole.then(|| line_point + params.dipole_offset.x * e1 + params.dipole_offset.y * e2),
    };
    if params.anisotropic_elasticity {
        let crystal = match params.elastic_symmetry {
            ElasticSymmetry::Cubic => make_cubic_elastic_tensor(params.elastic_c11, params.elastic_c12, params.elastic_c44),
            ElasticSymmetry::Hexagonal => {
                make_hexagonal_elastic_tensor(params.elastic_c11, params.elastic_c12, params.elastic_c13, params.elastic_c33, params.elastic_c44)
            }
            ElasticSymmetry::General => make_voigt_elastic_tensor(&params.elastic_voigt),
        };
        // Hexagonal: c is the 6-fold axis; otherwise the cell vectors are x, y, z.
        let col = |v: Vec3| v.as_dvec3();
        let (x_c, y_c, z_c) = if params.elastic_symmetry == ElasticSymmetry::Hexagonal {
            let z = glm::normalize_d(col(cell.z_axis));
            let a = col(cell.x_axis);
            let x = glm::normalize_d(a - a.dot(z) * z);
            (x, z.cross(x), z)
        } else {
            (glm::normalize_d(col(cell.x_axis)), glm::normalize_d(col(cell.y_axis)), glm::normalize_d(col(cell.z_axis)))
        };
        let frame = |v: DVec3| DVec3::new(v.dot(x_c), v.dot(y_c), v.dot(z_c));
        let local_frame = DMat3::from_cols(frame(col(e1)), frame(col(e2)), frame(col(line_dir)));
        let local = rotate_elastic_tensor(&crystal, &local_frame);
        let burgers_local = DVec3::new((edge_fraction * bmag) as f64, 0.0, (screw_fraction * bmag) as f64);
        let solution = solve_stroh(&local, burgers_local, params.elastic_noise_seed, params.elastic_noise_amplitude);
        if !solution.valid {
            result.message = format!("Anisotropic elasticity solve failed: {}", solution.error);
            return result;
        }
        kernel.stroh = Some(solution);
    }
    let mut output = base.clone();
    let atom_count = output.atoms.len();
    // The C++ splits 50000+ atoms into min(8, cores) contiguous chunks and
    // sums each chunk's displacement before adding the chunks.
    let hw = std::thread::available_parallelism().map_or(1, |n| n.get()).max(1);
    let workers = if atom_count >= 50000 { hw.min(8) } else { 1 };
    let chunk = atom_count.div_ceil(workers);
    let mut shifted = 0i32;
    let (mut sum_dx, mut sum_dy, mut sum_dz) = (0.0f64, 0.0f64, 0.0f64);
    for w in 0..workers {
        let (begin, end) = (w * chunk, ((w + 1) * chunk).min(atom_count));
        if begin >= end {
            break;
        }
        let (mut lx, mut ly, mut lz, mut ls) = (0.0f64, 0.0f64, 0.0f64, 0);
        for atom in &mut output.atoms[begin..end] {
            let p = Vec3::new(atom.x as f32, atom.y as f32, atom.z as f32);
            let d = displacement(params, &kernel, p);
            if dot(d, d) <= 1e-14 {
                continue;
            }
            atom.x += d.x as f64;
            atom.y += d.y as f64;
            atom.z += d.z as f64;
            lx += d.x as f64;
            ly += d.y as f64;
            lz += d.z as f64;
            ls += 1;
        }
        shifted += ls;
        sum_dx += lx;
        sum_dy += ly;
        sum_dz += lz;
    }
    // Centre-of-mass conservation for edge and mixed dislocations.
    if shifted > 0 && !output.atoms.is_empty() && !screw {
        let inv_n = 1.0 / output.atoms.len() as f64;
        let (mx, my, mz) = (sum_dx * inv_n, sum_dy * inv_n, sum_dz * inv_n);
        for atom in output.atoms.iter_mut() {
            atom.x -= mx;
            atom.y -= my;
            atom.z -= mz;
        }
    }
    result.success = shifted > 0;
    result.shifted_atom_count = shifted;
    result.burgers_magnitude = bmag;
    result.line_point = line_point;
    result.line_direction = line_dir;
    result.burgers_direction = burgers_dir;
    result.e1 = e1;
    result.e2 = e2;
    // The displaced region's outline: the convex hull (Jarvis march) of the
    // displaced atoms projected on the slip plane.
    if params.shape == DislocationShape::Freeform2D {
        result.loop_points = params.freeform_points.iter().map(|pt| line_point + pt.x * e1 + pt.y * e2).collect();
    } else {
        let threshold = max_f(0.05, 0.05 * bmag);
        let n = base.atoms.len().min(output.atoms.len());
        let mut pts: Vec<Vec2> = Vec::new();
        for i in 0..n {
            let b = Vec3::new(base.atoms[i].x as f32, base.atoms[i].y as f32, base.atoms[i].z as f32);
            let o = Vec3::new(output.atoms[i].x as f32, output.atoms[i].y as f32, output.atoms[i].z as f32);
            if length(o - b) < threshold {
                continue;
            }
            let rel = b - line_point;
            pts.push(Vec2::new(dot(rel, e1), dot(rel, line_dir)));
        }
        if pts.len() >= 3 {
            let n = pts.len();
            let mut start = 0;
            for i in 1..n {
                if pts[i].x < pts[start].x || (pts[i].x == pts[start].x && pts[i].y < pts[start].y) {
                    start = i;
                }
            }
            let mut hull = Vec::new();
            let mut cur = start;
            loop {
                hull.push(cur);
                let mut next = if cur == 0 { 1 } else { 0 };
                for i in 0..n {
                    if i == next {
                        continue;
                    }
                    let a = pts[next] - pts[cur];
                    let b = pts[i] - pts[cur];
                    let c = a.x * b.y - a.y * b.x;
                    if c > 0.0 || (c.abs() < 1e-9 && b.x * b.x + b.y * b.y > a.x * a.x + a.y * a.y) {
                        next = i;
                    }
                }
                cur = next;
                if cur == start || hull.len() > n {
                    break;
                }
            }
            result.loop_points = hull.iter().map(|&i| line_point + pts[i].x * e1 + pts[i].y * line_dir).collect();
        } else {
            result.loop_points = pts.iter().map(|p| line_point + p.x * e1 + p.y * line_dir).collect();
        }
    }
    result.output = output;
    if !result.success {
        result.message = "No atoms were shifted. Try increasing the shape radius/cutoff or moving the line position.".into();
        return result;
    }
    result.validation = validation_report(base, &result.output, detection.family);
    let character = match params.character {
        DislocationCharacter::Edge => "edge",
        DislocationCharacter::Screw => "screw",
        DislocationCharacter::Mixed => "mixed",
    };
    result.message = format!(
        "Inserted {character}{} on {}. Shifted atoms: {shifted}, |b| = {} A{} {}",
        if params.dipole { " dislocation dipole" } else { " dislocation" },
        dislocation_lattice_family_name(detection.family),
        g6(bmag),
        if params.dipole { " per partner (net Burgers vector zero)." } else { "." },
        result.validation.message
    );
    result
}
