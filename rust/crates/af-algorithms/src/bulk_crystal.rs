//! Bulk crystals from an asymmetric unit and a space group (the C++
//! algorithms/BulkCrystalBuilder): spglib's symmetry operations of the
//! group's Hall setting expand the atoms; duplicates within 0.01 A merge.

use af_core::elements;
use af_core::glm;
use af_core::{AtomSite, Structure};
use af_io::spacegroups;
use glam::{DMat3, DVec3};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CrystalSystem {
    Triclinic = 0,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BulkBuildResult {
    pub success: bool,
    pub message: String,
    pub generated_atoms: i32,
    pub space_group_number: i32,
    pub space_group_symbol: String,
}

pub struct SpaceGroupRange {
    pub system: CrystalSystem,
    pub label: &'static str,
    pub first: i32,
    pub last: i32,
}

pub const SPACE_GROUP_RANGES: [SpaceGroupRange; 7] = [
    SpaceGroupRange { system: CrystalSystem::Triclinic, label: "Triclinic", first: 1, last: 2 },
    SpaceGroupRange { system: CrystalSystem::Monoclinic, label: "Monoclinic", first: 3, last: 15 },
    SpaceGroupRange { system: CrystalSystem::Orthorhombic, label: "Orthorhombic", first: 16, last: 74 },
    SpaceGroupRange { system: CrystalSystem::Tetragonal, label: "Tetragonal", first: 75, last: 142 },
    SpaceGroupRange { system: CrystalSystem::Trigonal, label: "Trigonal", first: 143, last: 167 },
    SpaceGroupRange { system: CrystalSystem::Hexagonal, label: "Hexagonal", first: 168, last: 194 },
    SpaceGroupRange { system: CrystalSystem::Cubic, label: "Cubic", first: 195, last: 230 },
];

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LatticeParameters {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl Default for LatticeParameters {
    fn default() -> Self {
        Self { a: 4.0, b: 4.0, c: 4.0, alpha: 90.0, beta: 90.0, gamma: 90.0 }
    }
}

pub fn has_element_color(z: i32, colours: &[[f32; 3]]) -> bool {
    z >= 0 && (z as usize) < colours.len()
}

/// The colour of element z from a table, or its default colour.
pub fn element_colour(z: i32, colours: &[[f32; 3]]) -> [f32; 3] {
    if has_element_color(z, colours) {
        colours[z as usize]
    } else {
        elements::default_colour(z)
    }
}

pub fn apply_element_to_atom(atom: &mut AtomSite, z: i32, colours: &[[f32; 3]]) {
    atom.atomic_number = z;
    atom.symbol = elements::symbol(z).to_string();
    [atom.r, atom.g, atom.b] = element_colour(z, colours);
}

pub fn crystal_system_label(system: CrystalSystem) -> &'static str {
    SPACE_GROUP_RANGES.iter().find(|r| r.system == system).map_or("Unknown", |r| r.label)
}

pub fn current_range(system: CrystalSystem) -> &'static SpaceGroupRange {
    SPACE_GROUP_RANGES.iter().find(|r| r.system == system).unwrap_or(&SPACE_GROUP_RANGES[0])
}

/// Lattice vectors (columns a, b, c) from lengths and angles: a along x, b in the xy plane.
pub fn build_lattice_from_parameters(p: &LatticeParameters) -> DMat3 {
    let deg = 3.14159265358979323846 / 180.0;
    let (alpha, beta, gamma) = (p.alpha * deg, p.beta * deg, p.gamma * deg);
    let (cos_alpha, cos_beta, cos_gamma, sin_gamma) = (alpha.cos(), beta.cos(), gamma.cos(), gamma.sin());
    let a = DVec3::new(p.a, 0.0, 0.0);
    let b = DVec3::new(p.b * cos_gamma, p.b * sin_gamma, 0.0);
    let cx = p.c * cos_beta;
    let cy = if sin_gamma.abs() > 1e-10 { p.c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma } else { 0.0 };
    let mut cz2 = p.c * p.c - cx * cx - cy * cy;
    if cz2 < 0.0 && cz2 > -1e-8 {
        cz2 = 0.0;
    }
    let c = DVec3::new(cx, cy, if cz2 > 0.0 { cz2.sqrt() } else { 0.0 });
    DMat3::from_cols(a, b, c)
}

pub fn validate_parameters(system: CrystalSystem, p: &LatticeParameters) -> Result<(), String> {
    let near = |l: f64, r: f64| (l - r).abs() < 1e-6;
    if p.a <= 0.0 || p.b <= 0.0 || p.c <= 0.0 {
        return Err("Lattice lengths must be positive.".into());
    }
    let bad_angle = |x: f64| x <= 0.0 || x >= 180.0;
    if bad_angle(p.alpha) || bad_angle(p.beta) || bad_angle(p.gamma) {
        return Err("Angles must be in the range (0, 180).".into());
    }
    let right = |x: f64| near(x, 90.0);
    let error = match system {
        CrystalSystem::Monoclinic if !right(p.alpha) || !right(p.gamma) => Some("Monoclinic builder uses alpha=gamma=90 deg."),
        CrystalSystem::Orthorhombic if !right(p.alpha) || !right(p.beta) || !right(p.gamma) => {
            Some("Orthorhombic builder uses alpha=beta=gamma=90 deg.")
        }
        CrystalSystem::Tetragonal if !near(p.a, p.b) || !right(p.alpha) || !right(p.beta) || !right(p.gamma) => {
            Some("Tetragonal builder uses a=b and alpha=beta=gamma=90 deg.")
        }
        CrystalSystem::Trigonal if !near(p.a, p.b) || !right(p.alpha) || !right(p.beta) || !near(p.gamma, 120.0) => {
            Some("Trigonal builder uses the hexagonal setting: a=b, alpha=beta=90 deg, gamma=120 deg.")
        }
        CrystalSystem::Hexagonal if !near(p.a, p.b) || !right(p.alpha) || !right(p.beta) || !near(p.gamma, 120.0) => {
            Some("Hexagonal builder uses a=b, alpha=beta=90 deg, gamma=120 deg.")
        }
        CrystalSystem::Cubic if !near(p.a, p.b) || !near(p.b, p.c) || !right(p.alpha) || !right(p.beta) || !right(p.gamma) => {
            Some("Cubic builder uses a=b=c and alpha=beta=gamma=90 deg.")
        }
        _ => None,
    };
    if let Some(e) = error {
        return Err(e.into());
    }
    if glm::determinant_d(&build_lattice_from_parameters(p)).abs() < 1e-10 {
        return Err("Lattice parameters produce a zero-volume cell.".into());
    }
    Ok(())
}

pub fn apply_system_constraints(system: CrystalSystem, p: &mut LatticeParameters) {
    match system {
        CrystalSystem::Monoclinic => {
            p.alpha = 90.0;
            p.gamma = 90.0;
        }
        CrystalSystem::Orthorhombic => {
            p.alpha = 90.0;
            p.beta = 90.0;
            p.gamma = 90.0;
        }
        CrystalSystem::Tetragonal => {
            p.b = p.a;
            p.alpha = 90.0;
            p.beta = 90.0;
            p.gamma = 90.0;
        }
        CrystalSystem::Trigonal | CrystalSystem::Hexagonal => {
            p.b = p.a;
            p.alpha = 90.0;
            p.beta = 90.0;
            p.gamma = 120.0;
        }
        CrystalSystem::Cubic => {
            p.b = p.a;
            p.c = p.a;
            p.alpha = 90.0;
            p.beta = 90.0;
            p.gamma = 90.0;
        }
        CrystalSystem::Triclinic => {}
    }
}

/// The first Hall number of each space group (index 1-230; 0 = none).
pub fn hall_by_space_group() -> &'static [i32; 231] {
    static TABLE: std::sync::OnceLock<[i32; 231]> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut by_number = [0i32; 231];
        for hall in 1..=530 {
            let number = spacegroups::of_hall(hall).number;
            if (1..=230).contains(&number) && by_number[number as usize] == 0 {
                by_number[number as usize] = hall;
            }
        }
        by_number
    })
}

/// "225 Fm-3m".
pub fn space_group_label(number: i32) -> String {
    if (1..=230).contains(&number) && hall_by_space_group()[number as usize] > 0 {
        let t = spacegroups::of_hall(hall_by_space_group()[number as usize]);
        return format!("{} {}", t.number, t.international_short);
    }
    number.to_string()
}

/// The Hall number of a space group's origin choice 2, or 0.
fn second_origin_hall(number: i32) -> i32 {
    if !(1..=230).contains(&number) || hall_by_space_group()[number as usize] == 0 {
        return 0;
    }
    for hall in hall_by_space_group()[number as usize]..=530 {
        let t = spacegroups::of_hall(hall);
        if t.number != number {
            break;
        }
        if t.choice == "2" {
            return hall;
        }
    }
    0
}

/// Whether a space group has two origin choices (Fd-3m and 23 others).
pub fn has_two_origins(number: i32) -> bool {
    second_origin_hall(number) != 0
}

/// Fractions wrapped into [0, 1), snapped to simple fractions (1/3 as 0.3333333).
pub fn wrap_fractional(frac: DVec3) -> DVec3 {
    let mut w = frac;
    w.x -= w.x.floor();
    w.y -= w.y.floor();
    w.z -= w.z.floor();
    for k in 0..3 {
        if w[k] >= 1.0 - 1e-9 {
            w[k] = 0.0;
        }
        for d in [2.0, 3.0, 4.0, 6.0, 8.0, 12.0] {
            let n = (w[k] * d).round();
            if (w[k] * d - n).abs() < 1e-6 * d {
                w[k] = n / d;
                break;
            }
        }
        if w[k] >= 1.0 {
            w[k] = 0.0;
        }
    }
    w
}

pub fn same_fractional_position(lhs: DVec3, rhs: DVec3, tol: f64) -> bool {
    let mut d = lhs - rhs;
    d.x -= d.x.round();
    d.y -= d.y.round();
    d.z -= d.z.round();
    d.x.abs() <= tol && d.y.abs() <= tol && d.z.abs() <= tol
}

/// Expands the asymmetric atoms (fractional positions in x, y, z) by the
/// space group's operations into `structure`.
pub fn build_bulk_crystal(
    structure: &mut Structure,
    system: CrystalSystem,
    space_group_number: i32,
    params: &LatticeParameters,
    asymmetric_atoms: &[AtomSite],
    colours: &[[f32; 3]],
    origin_choice: i32,
) -> BulkBuildResult {
    let mut result =
        BulkBuildResult { space_group_number, space_group_symbol: space_group_label(space_group_number), ..Default::default() };
    if origin_choice != 1 && origin_choice != 2 {
        result.message = "The origin choice must be 1 or 2.".into();
        return result;
    }
    if has_two_origins(space_group_number) {
        result.space_group_symbol += &format!(", origin choice {origin_choice}");
    } else if origin_choice == 2 {
        result.message = format!("Space group {space_group_number} has only one origin choice.");
        return result;
    }
    if let Err(e) = validate_parameters(system, params) {
        result.message = e;
        return result;
    }
    if asymmetric_atoms.is_empty() {
        result.message = "Add at least one asymmetric-unit atom before building.".into();
        return result;
    }
    if !(1..=230).contains(&space_group_number) || hall_by_space_group()[space_group_number as usize] == 0 {
        result.message = "Selected space group is not available in the symmetry database.".into();
        return result;
    }
    let hall = if origin_choice == 2 { second_origin_hall(space_group_number) } else { hall_by_space_group()[space_group_number as usize] };
    let operations = &spacegroups::of_hall(hall).operations;
    if operations.is_empty() {
        result.message = "Could not load symmetry operations for the selected space group.".into();
        return result;
    }
    let lattice = build_lattice_from_parameters(params);
    let (a, b, c) = (lattice.x_axis, lattice.y_axis, lattice.z_axis);
    let mut generated: Vec<AtomSite> = Vec::new();
    let mut generated_frac: Vec<DVec3> = Vec::new();
    for basis in asymmetric_atoms {
        let frac = wrap_fractional(DVec3::new(basis.x, basis.y, basis.z));
        for (rotation, translation) in operations {
            let mut t = DVec3::ZERO;
            for row in 0..3 {
                t[row] = translation[row];
                for col in 0..3 {
                    t[row] += rotation.col(col)[row] * frac[col];
                }
            }
            let t = wrap_fractional(t);
            // Within 0.01 A: the same site, typed with rounded fractions.
            let duplicate = generated.iter().zip(&generated_frac).any(|(existing, &f)| {
                if existing.atomic_number != basis.atomic_number {
                    return false;
                }
                let mut d = t - f;
                d -= DVec3::new(d.x.round(), d.y.round(), d.z.round());
                let v = a * d.x + b * d.y + c * d.z;
                v.dot(v).sqrt() < 0.01
            });
            if duplicate {
                continue;
            }
            let cart = a * t.x + b * t.y + c * t.z;
            let z = basis.atomic_number;
            let [r, g, bl] = element_colour(z, colours);
            generated.push(AtomSite {
                atomic_number: z,
                symbol: elements::symbol(z).to_string(),
                x: cart.x,
                y: cart.y,
                z: cart.z,
                r,
                g,
                b: bl,
                ..Default::default()
            });
            generated_frac.push(t);
        }
    }
    if generated.is_empty() {
        result.message = "Symmetry expansion produced no atoms.".into();
        return result;
    }
    structure.atoms = generated;
    structure.has_unit_cell = true;
    structure.cell_offset = [0.0; 3];
    structure.cell_vectors = [a.to_array(), b.to_array(), c.to_array()];
    result.generated_atoms = structure.atoms.len() as i32;
    result.success = true;
    result.message = "Bulk crystal generated from asymmetric-unit atoms and selected space-group symmetry.".into();
    result
}

pub fn add_default_asymmetric_atom(atoms: &mut Vec<AtomSite>, colours: &[[f32; 3]]) {
    let mut atom = AtomSite::default();
    apply_element_to_atom(&mut atom, 1, colours);
    atoms.push(atom);
}
