//! Amorphous structures by random sequential addition (the C++
//! algorithms/AmorphousBuilder): atoms at uniform random positions in an
//! orthogonal periodic box, accepted when no placed atom is nearer than the
//! pair's minimum distance (std::mt19937 and float uniforms as the C++).

use std::collections::HashMap;

use af_core::cxx_random::{shuffle, uniform_real_f32, Mt19937};
use af_core::elements;
use af_core::{AtomSite, Structure};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AmorphousBoxMode {
    /// a, b, c given in Angstrom.
    Manual = 0,
    /// A cube whose volume gives the target density.
    #[default]
    AutoFromDensity,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AmorphousElementSpec {
    pub atomic_number: i32,
    pub count: i32,
}

impl Default for AmorphousElementSpec {
    fn default() -> Self {
        Self { atomic_number: 14, count: 100 }
    }
}

/// An explicit minimum separation for one element pair.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AmorphousPairDist {
    pub z1: i32,
    pub z2: i32,
    pub min_dist: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AmorphousParams {
    pub elements: Vec<AmorphousElementSpec>,
    pub box_mode: AmorphousBoxMode,
    pub box_a: f32,
    pub box_b: f32,
    pub box_c: f32,
    /// g/cm^3 (AutoFromDensity).
    pub target_density: f32,
    /// Applied to every box vector before packing.
    pub cell_scale_factor: f32,
    pub pair_distances: Vec<AmorphousPairDist>,
    /// Fraction of r_cov1 + r_cov2 for pairs without an explicit distance.
    pub covalent_tolerance: f32,
    /// 0 seeds from the clock.
    pub seed: u32,
    pub max_attempts: i32,
    pub periodic: bool,
}

impl Default for AmorphousParams {
    fn default() -> Self {
        Self {
            elements: Vec::new(),
            box_mode: AmorphousBoxMode::AutoFromDensity,
            box_a: 20.0,
            box_b: 20.0,
            box_c: 20.0,
            target_density: 2.0,
            cell_scale_factor: 1.0,
            pair_distances: Vec::new(),
            covalent_tolerance: 0.75,
            seed: 42,
            max_attempts: 1000,
            periodic: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AmorphousResult {
    pub success: bool,
    pub message: String,
    pub output: Structure,
    pub requested_atoms: i32,
    pub placed_atoms: i32,
    /// g/cm^3 of the placed atoms in the box.
    pub actual_density: f32,
    pub box_a: f32,
    pub box_b: f32,
    pub box_c: f32,
}

fn max_f(a: f32, b: f32) -> f32 {
    if a < b {
        b
    } else {
        a
    }
}

/// `covalent_radii` and `colours` are indexed by atomic number.
pub fn build_amorphous_structure(params: &AmorphousParams, covalent_radii: &[f32], colours: &[[f32; 3]]) -> AmorphousResult {
    let mut result = AmorphousResult::default();
    if params.elements.is_empty() {
        result.message = "No element species specified.".into();
        return result;
    }
    let mut total_atoms: i32 = 0;
    for e in &params.elements {
        if !(1..=118).contains(&e.atomic_number) {
            result.message = "Invalid atomic number in composition.".into();
            return result;
        }
        if e.count < 0 {
            result.message = format!("Negative atom count for element Z={}", e.atomic_number);
            return result;
        }
        if e.count > i32::MAX - total_atoms {
            result.message = "Total atom count exceeds the supported range.".into();
            return result;
        }
        total_atoms += e.count;
    }
    if total_atoms == 0 {
        result.message = "Total atom count is zero.".into();
        return result;
    }
    result.requested_atoms = total_atoms;

    let (mut box_a, mut box_b, mut box_c);
    if params.box_mode == AmorphousBoxMode::Manual {
        (box_a, box_b, box_c) = (params.box_a, params.box_b, params.box_c);
    } else {
        let mut total_mass = 0.0f64;
        for e in &params.elements {
            total_mass += e.count as f64 * elements::atomic_mass(e.atomic_number);
        }
        // volume [A^3] = mass [u] * 1.66054 / density [g/cm^3]
        let density = (params.target_density as f64).max(1e-3);
        let side = (total_mass * 1.66054 / density).cbrt() as f32;
        (box_a, box_b, box_c) = (side, side, side);
    }
    box_a *= params.cell_scale_factor;
    box_b *= params.cell_scale_factor;
    box_c *= params.cell_scale_factor;
    if !box_a.is_finite() || !box_b.is_finite() || !box_c.is_finite() || box_a < 0.1 || box_b < 0.1 || box_c < 0.1 {
        result.message = "Computed box dimensions must be finite and at least 0.1 \u{c5}.".into();
        return result;
    }
    result.box_a = box_a;
    result.box_b = box_b;
    result.box_c = box_c;

    let mut pair_dist: HashMap<i32, f32> = HashMap::new();
    for pd in &params.pair_distances {
        pair_dist.insert(pd.z1.min(pd.z2) * 200 + pd.z1.max(pd.z2), pd.min_dist);
    }
    let min_dist = |z1: i32, z2: i32| -> f32 {
        if let Some(&d) = pair_dist.get(&(z1.min(z2) * 200 + z1.max(z2))) {
            return d;
        }
        let radius = |z: i32| if z >= 1 && (z as usize) < covalent_radii.len() { covalent_radii[z as usize] } else { 1.5 };
        (radius(z1) + radius(z2)) * params.covalent_tolerance
    };
    let mut max_min_dist = 0.01f32;
    for e1 in &params.elements {
        for e2 in &params.elements {
            max_min_dist = max_f(max_min_dist, min_dist(e1.atomic_number, e2.atomic_number));
        }
    }
    let mut atom_list: Vec<i32> = Vec::with_capacity(total_atoms as usize);
    for e in &params.elements {
        atom_list.extend(std::iter::repeat(e.atomic_number).take(e.count as usize));
    }
    let seed = if params.seed != 0 {
        params.seed
    } else {
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(1, |d| d.as_secs() as u32)
    };
    let mut rng = Mt19937::new(seed);
    shuffle(&mut atom_list, &mut rng);

    // Grid cells of size maxMinDist: the 3x3x3 neighbourhood holds every conflict.
    let cell_size = max_min_dist;
    let num = |b: f32| ((b / cell_size).floor() as i32).max(1);
    let (num_a, num_b, num_c) = (num(box_a), num(box_b), num(box_c));
    let wrap = |v: i32, n: i32| ((v % n) + n) % n;
    let cell_of = |x: f32, n: i32| wrap((x / cell_size).floor() as i32, n);
    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    let mut positions: Vec<[f32; 3]> = Vec::with_capacity(total_atoms as usize);
    let mut placed_z: Vec<i32> = Vec::with_capacity(total_atoms as usize);
    let (mut placed, mut skipped) = (0, 0);
    for &z in &atom_list {
        let mut success = false;
        for _ in 0..params.max_attempts {
            let cx = uniform_real_f32(&mut rng, 0.0, box_a);
            let cy = uniform_real_f32(&mut rng, 0.0, box_b);
            let cz = uniform_real_f32(&mut rng, 0.0, box_c);
            let (gx, gy, gz) = (cell_of(cx, num_a), cell_of(cy, num_b), cell_of(cz, num_c));
            let mut conflict = false;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let key = (wrap(gx + dx, num_a), wrap(gy + dy, num_b), wrap(gz + dz, num_c));
                        let Some(list) = grid.get(&key) else { continue };
                        for &idx in list {
                            // Minimum-image distance.
                            let mut ddx = cx - positions[idx][0];
                            let mut ddy = cy - positions[idx][1];
                            let mut ddz = cz - positions[idx][2];
                            ddx -= (ddx / box_a).round() * box_a;
                            ddy -= (ddy / box_b).round() * box_b;
                            ddz -= (ddz / box_c).round() * box_c;
                            let dist2 = ddx * ddx + ddy * ddy + ddz * ddz;
                            let d = min_dist(z, placed_z[idx]);
                            if dist2 < d * d {
                                conflict = true;
                                break 'search;
                            }
                        }
                    }
                }
            }
            if !conflict {
                grid.entry((gx, gy, gz)).or_default().push(positions.len());
                positions.push([cx, cy, cz]);
                placed_z.push(z);
                placed += 1;
                success = true;
                break;
            }
        }
        if !success {
            skipped += 1;
        }
    }
    let mut output = Structure::default();
    for (p, &z) in positions.iter().zip(&placed_z) {
        let mut site = AtomSite {
            atomic_number: z,
            symbol: elements::symbol(z).to_string(),
            x: p[0] as f64,
            y: p[1] as f64,
            z: p[2] as f64,
            ..Default::default()
        };
        [site.r, site.g, site.b] = elements::default_colour(z);
        if z >= 0 && (z as usize) < colours.len() {
            [site.r, site.g, site.b] = colours[z as usize];
        }
        output.atoms.push(site);
    }
    if params.periodic {
        output.has_unit_cell = true;
        output.cell_vectors = [[box_a as f64, 0.0, 0.0], [0.0, box_b as f64, 0.0], [0.0, 0.0, box_c as f64]];
        output.cell_offset = [0.0; 3];
    }
    let mut placed_mass = 0.0f64;
    for &z in &placed_z {
        placed_mass += elements::atomic_mass(z);
    }
    let volume = box_a as f64 * box_b as f64 * box_c as f64;
    result.actual_density = (placed_mass * 1.66054 / volume) as f32;
    result.output = output;
    result.placed_atoms = placed;
    if skipped == 0 {
        result.success = true;
        result.message = format!("Successfully placed {placed} atoms.");
    } else {
        result.success = placed > 0;
        result.message = format!(
            "Placed {placed} / {total_atoms} atoms. {skipped} atom(s) could not be placed \u{2014} try increasing box size, scale factor, or max attempts."
        );
    }
    result
}
