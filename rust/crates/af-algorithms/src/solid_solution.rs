//! Substitutional solid solutions (the C++ algorithms/
//! SubstitutionalSolidSolutionBuilder): the sites of a host structure get
//! random species in the requested proportions (largest-remainder counts,
//! a seeded std::mt19937 shuffle).

use af_core::cxx::sort_by;
use af_core::cxx_random::{shuffle, Mt19937};
use af_core::elements;
use af_core::Structure;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ElementFraction {
    pub atomic_number: i32,
    /// Target fraction in [0, 1]; fractions are normalised.
    pub fraction: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SssParams {
    pub composition: Vec<ElementFraction>,
    /// 0 seeds from the clock.
    pub seed: u32,
}

impl Default for SssParams {
    fn default() -> Self {
        Self { composition: Vec::new(), seed: 12345 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SssResult {
    pub success: bool,
    pub message: String,
    pub output: Structure,
    /// (atomic number, count) in the composition's order.
    pub element_counts: Vec<(i32, i32)>,
}

pub fn build_substitutional_solid_solution(base: &Structure, params: &SssParams) -> SssResult {
    let mut result = SssResult::default();
    if base.atoms.is_empty() {
        result.message = "Base structure has no atoms.".into();
        return result;
    }
    if params.composition.is_empty() {
        result.message = "No composition entries specified.".into();
        return result;
    }
    for ef in &params.composition {
        if !(1..=118).contains(&ef.atomic_number) {
            result.message = "Invalid atomic number in composition.".into();
            return result;
        }
        if !ef.fraction.is_finite() || ef.fraction < 0.0 {
            result.message = "Composition fractions must be finite and non-negative.".into();
            return result;
        }
    }
    let mut total = 0.0f64;
    for ef in &params.composition {
        total += ef.fraction as f64;
    }
    if total < 1e-6f32 as f64 {
        result.message = "All composition fractions are zero.".into();
        return result;
    }
    let n = base.atoms.len() as i32;
    let num_el = params.composition.len();
    let mut counts = vec![0i32; num_el];
    let mut remainders = vec![0.0f64; num_el];
    let mut assigned = 0;
    for i in 0..num_el {
        let exact = (params.composition[i].fraction as f64 / total) * n as f64;
        counts[i] = exact.floor() as i32;
        remainders[i] = exact - counts[i] as f64;
        assigned += counts[i];
    }
    // Leftover sites by largest remainder.
    let leftover = n - assigned;
    if leftover > 0 {
        let mut order: Vec<usize> = (0..num_el).collect();
        sort_by(&mut order, |&a, &b| remainders[a] > remainders[b]);
        for k in 0..leftover as usize {
            counts[order[k % num_el]] += 1;
        }
    }
    let seed = if params.seed != 0 {
        params.seed
    } else {
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(1, |d| d.as_secs() as u32)
    };
    let mut rng = Mt19937::new(seed);
    let mut indices: Vec<usize> = (0..n as usize).collect();
    shuffle(&mut indices, &mut rng);
    let mut atomic_numbers = vec![0i32; n as usize];
    let mut pos = 0;
    for (i, ef) in params.composition.iter().enumerate() {
        for _ in 0..counts[i] {
            atomic_numbers[indices[pos]] = ef.atomic_number;
            pos += 1;
        }
    }
    result.output = base.clone();
    result.output.grain_colors.clear();
    result.output.grain_region_ids.clear();
    let colours = elements::default_colours();
    for (atom, &z) in result.output.atoms.iter_mut().zip(&atomic_numbers) {
        atom.atomic_number = z;
        atom.symbol = elements::symbol(z).to_string();
        if z >= 1 && (z as usize) < colours.len() {
            [atom.r, atom.g, atom.b] = colours[z as usize];
        }
    }
    result.element_counts = params.composition.iter().zip(&counts).map(|(ef, &c)| (ef.atomic_number, c)).collect();
    let mut message = format!("Built substitutional solid solution: {n} atoms");
    for (ef, c) in params.composition.iter().zip(&counts) {
        message.push_str(&format!(",  {} {c}", elements::symbol(ef.atomic_number)));
    }
    result.message = message;
    result.success = true;
    result
}
