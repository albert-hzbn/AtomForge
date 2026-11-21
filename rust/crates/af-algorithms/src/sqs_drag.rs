//! Special quasirandom structures by simulated annealing of species swaps
//! (the C++ algorithms/SQSBuilder) and drag-method interpolation between two
//! structures (algorithms/DragPrep).

use std::collections::BTreeMap;

use af_core::cxx::sort_by;
use af_core::cxx_random::{shuffle, uniform_i64, uniform_real, Mt19937};
use af_core::eigen;
use af_core::elements;
use af_core::{AtomSite, Structure};
use glam::DVec3;

use crate::sro::{build_pairs, build_shell_distances, SroPair};

#[derive(Clone, Debug, PartialEq)]
pub struct SqsParams {
    /// Element symbol -> fraction (summing to about 1).
    pub composition: BTreeMap<String, f64>,
    pub shells: i32,
    pub shell_tolerance: f64,
    pub steps: i32,
    pub start_temperature: f64,
    pub end_temperature: f64,
    pub seed: u32,
}

impl Default for SqsParams {
    fn default() -> Self {
        Self { composition: BTreeMap::new(), shells: 2, shell_tolerance: 0.2, steps: 3000, start_temperature: 1.0, end_temperature: 0.02, seed: 1 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SqsResult {
    pub success: bool,
    pub message: String,
    pub structure: Structure,
    pub initial_objective: f64,
    pub final_objective: f64,
    pub accepted_swaps: i32,
}

fn set_species(atom: &mut AtomSite, z: i32) {
    atom.atomic_number = z;
    atom.symbol = elements::symbol(z).to_string();
    [atom.r, atom.g, atom.b] = elements::default_colour(z);
}

/// Sum of squared Warren-Cowley parameters over the neighbour pairs.
fn objective(pairs: &[SroPair], shells: usize, kind: &[usize], concentration: &[f64]) -> f64 {
    let m = concentration.len();
    let mut n = vec![0.0f64; shells * m * m];
    for p in pairs {
        n[(p.shell * m + kind[p.i]) * m + kind[p.j]] += 1.0;
    }
    let mut sum = 0.0;
    for s in 0..shells {
        for a in 0..m {
            let mut z = 0.0;
            for b in 0..m {
                z += n[(s * m + a) * m + b];
            }
            if z == 0.0 {
                continue;
            }
            for b in 0..m {
                let alpha = if concentration[b] > 1e-12 { 1.0 - n[(s * m + a) * m + b] / z / concentration[b] } else { 0.0 };
                sum += alpha * alpha;
            }
        }
    }
    sum
}

/// std::to_string(double): "%f".
fn f6(x: f64) -> String {
    format!("{x:.6}")
}

pub fn build_sqs(source: &Structure, params: &SqsParams) -> SqsResult {
    let mut result = SqsResult::default();
    if source.atoms.is_empty() {
        result.message = "Structure has no atoms.".into();
        return result;
    }
    if params.composition.len() < 2 {
        result.message = "composition needs at least two elements.".into();
        return result;
    }
    let mut total = 0.0;
    let mut fractions: Vec<(i32, f64)> = Vec::new();
    for (symbol, &fraction) in &params.composition {
        let z = elements::atomic_number(symbol);
        if z == 0 {
            result.message = format!("Unknown element symbol '{symbol}'.");
            return result;
        }
        if fraction < 0.0 {
            result.message = "Composition fractions must be non-negative.".into();
            return result;
        }
        fractions.push((z, fraction));
        total += fraction;
    }
    if (total - 1.0).abs() > 1e-2 {
        result.message = format!("Composition fractions must sum to 1.0 (got {}).", f6(total));
        return result;
    }
    if params.steps < 0 {
        result.message = "steps must be non-negative.".into();
        return result;
    }
    let n = source.atoms.len();
    let mut rng = Mt19937::new(params.seed);
    // Largest-remainder counts, shuffled onto the fixed sites.
    let mut counts = vec![0i32; fractions.len()];
    let mut remainders = vec![0.0f64; fractions.len()];
    let mut assigned = 0;
    for (i, &(_, fraction)) in fractions.iter().enumerate() {
        let exact = fraction * n as f64;
        counts[i] = exact.floor() as i32;
        remainders[i] = exact - counts[i] as f64;
        assigned += counts[i];
    }
    let mut order: Vec<usize> = (0..fractions.len()).collect();
    sort_by(&mut order, |&a, &b| remainders[a] > remainders[b]);
    for &k in &order {
        if assigned >= n as i32 {
            break;
        }
        counts[k] += 1;
        assigned += 1;
    }
    let mut species: Vec<i32> = Vec::with_capacity(n);
    for (i, &(z, _)) in fractions.iter().enumerate() {
        species.extend(std::iter::repeat(z).take(counts[i].max(0) as usize));
    }
    shuffle(&mut species, &mut rng);
    let mut structure = source.clone();
    for (atom, &z) in structure.atoms.iter_mut().zip(&species) {
        set_species(atom, z);
    }
    // The geometry never changes: the neighbour pairs once, species recounted per swap.
    let shell_distances = build_shell_distances(&structure, params.shells, params.shell_tolerance);
    let pairs = build_pairs(&structure, &shell_distances, params.shell_tolerance);
    let shells = shell_distances.len();
    let concentration: Vec<f64> = counts.iter().map(|&c| c as f64 / n as f64).collect();
    let mut kind = vec![0usize; n];
    for i in 0..n {
        for (k, &(z, _)) in fractions.iter().enumerate() {
            if z == species[i] {
                kind[i] = k;
            }
        }
    }
    result.initial_objective = objective(&pairs, shells, &kind, &concentration);
    let mut current = result.initial_objective;
    let mut best = structure.clone();
    let mut best_objective = current;
    for step in 0..params.steps {
        let i = uniform_i64(&mut rng, 0, n as i64 - 1) as usize;
        let mut j = uniform_i64(&mut rng, 0, n as i64 - 1) as usize;
        let mut guard = 0;
        while j == i || structure.atoms[j].atomic_number == structure.atoms[i].atomic_number {
            let keep_going = guard < 64;
            guard += 1;
            if !keep_going {
                break;
            }
            j = uniform_i64(&mut rng, 0, n as i64 - 1) as usize;
        }
        if structure.atoms[j].atomic_number == structure.atoms[i].atomic_number {
            continue;
        }
        let (zi, zj) = (structure.atoms[i].atomic_number, structure.atoms[j].atomic_number);
        set_species(&mut structure.atoms[i], zj);
        set_species(&mut structure.atoms[j], zi);
        kind.swap(i, j);
        let candidate = objective(&pairs, shells, &kind, &concentration);
        let delta = candidate - current;
        let fraction = if params.steps > 1 { step as f64 / (params.steps - 1) as f64 } else { 1.0 };
        let t = params.start_temperature + (params.end_temperature - params.start_temperature) * fraction;
        let accept = delta <= 0.0 || (t > 1e-9 && (-delta / t).exp() > uniform_real(&mut rng, 0.0, 1.0));
        if accept {
            current = candidate;
            result.accepted_swaps += 1;
            if current < best_objective {
                best_objective = current;
                best = structure.clone();
            }
        } else {
            set_species(&mut structure.atoms[i], zi);
            set_species(&mut structure.atoms[j], zj);
            kind.swap(i, j);
        }
    }
    result.structure = best;
    result.final_objective = best_objective;
    result.success = true;
    result.message = format!(
        "SQS optimization: sum(alpha^2) over {} shells went from {} to {} over {} steps ({} accepted swaps).",
        params.shells,
        f6(result.initial_objective),
        f6(result.final_objective),
        params.steps,
        result.accepted_swaps
    );
    result
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DragPrepResult {
    pub success: bool,
    pub message: String,
    /// initial + zeta (final - initial).
    pub interpolated: Structure,
    /// final - initial per atom (minimum image with a cell).
    pub constraint_direction: Vec<DVec3>,
}

pub fn prepare_drag(initial: &Structure, final_: &Structure, zeta: f64, clip_displacement: bool, use_pbc: bool) -> DragPrepResult {
    let mut result = DragPrepResult::default();
    let n = initial.atoms.len();
    if n != final_.atoms.len() {
        result.message = "Initial and final structures must have the same atom count and ordering.".into();
        return result;
    }
    if n == 0 {
        result.message = "Structures have no atoms.".into();
        return result;
    }
    if !(0.0..=1.0).contains(&zeta) {
        result.message = "zeta must be in [0, 1].".into();
        return result;
    }
    // Columns of the cell are a, b and c.
    let mut cell = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cell[j][i] = initial.cell_vectors[i][j];
        }
    }
    let have_pbc = clip_displacement && use_pbc && initial.has_unit_cell && eigen::determinant3(&cell).abs() >= 1e-12;
    let inverse = if have_pbc { eigen::inverse3(&cell) } else { [[0.0; 3]; 3] };
    result.interpolated = initial.clone();
    for i in 0..n {
        let (a, b) = (&initial.atoms[i], &final_.atoms[i]);
        let mut delta = [b.x - a.x, b.y - a.y, b.z - a.z];
        if have_pbc {
            let mut frac = eigen::mul3(&inverse, delta);
            for f in frac.iter_mut() {
                *f -= f.round();
            }
            delta = eigen::mul3(&cell, frac);
        }
        result.constraint_direction.push(DVec3::from(delta));
        let atom = &mut result.interpolated.atoms[i];
        atom.x = a.x + zeta * delta[0];
        atom.y = a.y + zeta * delta[1];
        atom.z = a.z + zeta * delta[2];
    }
    result.success = true;
    result.message = format!("Interpolated {n} atoms at zeta={}.", f6(zeta));
    result
}
