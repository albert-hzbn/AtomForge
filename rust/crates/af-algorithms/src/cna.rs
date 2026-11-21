//! Common neighbour analysis (the C++ algorithms/CommonNeighbourAnalysis):
//! for every bonded pair the number of common neighbours, the bonds among
//! them and their longest chain; per atom the dominant signature and an
//! FCC/HCP/BCC/ICO environment.

use std::collections::BTreeMap;

use af_core::cell::MIN_BOND_DISTANCE;
use af_core::elements;
use af_core::glm;
use af_core::neighbours;
use af_core::Structure;
use glam::{DMat3, DVec3};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Signature {
    pub common: i32,
    pub bonds: i32,
    pub chain: i32,
}

impl Signature {
    pub const fn new(common: i32, bonds: i32, chain: i32) -> Self {
        Self { common, bonds, chain }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AtomRow {
    pub index: i32,
    pub atomic_number: i32,
    pub symbol: String,
    pub coordination: i32,
    pub dominant_signature: Signature,
    pub dominant_signature_count: i32,
    pub environment: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CnaResult {
    pub valid: bool,
    pub message: String,
    pub atom_count: i32,
    pub pair_count: i32,
    pub pbc_used: bool,
    pub signature_counts: BTreeMap<Signature, i32>,
    pub environment_counts: BTreeMap<String, i32>,
    pub atom_rows: Vec<AtomRow>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CnaParams {
    pub cutoff_scale: f32,
    pub use_pbc: bool,
    /// Angstrom; zero selects scaled covalent radii.
    pub cutoff_radius: f32,
}

impl Default for CnaParams {
    fn default() -> Self {
        Self { cutoff_scale: 1.18, use_pbc: true, cutoff_radius: 0.0 }
    }
}

/// Atom j in periodic image n, relative to the central atom in image 0.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Neighbour {
    j: i32,
    n: [i32; 3],
}

fn classify_environment(signatures: &BTreeMap<Signature, i32>, coordination: i32) -> &'static str {
    let count = |s: Signature| signatures.get(&s).copied().unwrap_or(0);
    if coordination == 12 && count(Signature::new(4, 2, 1)) == 12 {
        return "FCC-like";
    }
    if coordination == 12 && count(Signature::new(4, 2, 1)) == 6 && count(Signature::new(4, 2, 2)) == 6 {
        return "HCP-like";
    }
    if coordination == 14 && count(Signature::new(4, 4, 4)) == 6 && count(Signature::new(6, 6, 6)) == 8 {
        return "BCC-like";
    }
    if coordination == 12 && count(Signature::new(5, 5, 5)) == 12 {
        return "ICO-like";
    }
    "Unknown"
}

/// The longest chain of bonds among the common neighbours (a ring counts its
/// closing bond); -1 when the bounded search gives up.
fn longest_chain_length(node_count: usize, bonded: impl Fn(usize, usize) -> bool) -> i32 {
    if node_count > 16 {
        return -1;
    }
    let mut adjacency = vec![Vec::new(); node_count];
    for i in 0..node_count {
        for j in i + 1..node_count {
            if bonded(i, j) {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
        }
    }
    struct Search<'a> {
        adjacency: &'a [Vec<usize>],
        start: usize,
        best: i32,
        visits: i32,
        exhausted: bool,
    }
    fn visit(s: &mut Search, node: usize, visited: u32, length: i32) {
        s.visits += 1;
        if s.visits > 200000 {
            s.exhausted = true;
            return;
        }
        s.best = s.best.max(length);
        for &next in &s.adjacency[node] {
            if s.exhausted {
                return;
            }
            if next == s.start && length >= 2 {
                s.best = s.best.max(length + 1);
            } else if visited & (1 << next) == 0 {
                visit(s, next, visited | (1 << next), length + 1);
            }
        }
    }
    let mut search = Search { adjacency: &adjacency, start: 0, best: 0, visits: 0, exhausted: false };
    for start in 0..node_count {
        search.start = start;
        visit(&mut search, start, 1 << start, 0);
    }
    if search.exhausted {
        -1
    } else {
        search.best
    }
}

fn run_cna(structure: &Structure, cutoff_scale: f32, use_pbc_request: bool, cutoff_radius: f32) -> CnaResult {
    let mut result = CnaResult { atom_count: structure.atoms.len() as i32, ..Default::default() };
    if structure.atoms.is_empty() {
        result.message = "No atoms available.".into();
        return result;
    }
    let radii = elements::literature_covalent_radii();
    let radius_of = |z: i32| if z >= 0 && (z as usize) < radii.len() { radii[z as usize] } else { 1.0 };
    let use_pbc = use_pbc_request && structure.has_unit_cell && glm::cell_matrices(structure).is_some();
    result.pbc_used = use_pbc;

    let mut maximum_cutoff = cutoff_radius;
    if maximum_cutoff == 0.0 {
        for atom in &structure.atoms {
            let candidate = 2.0 * radius_of(atom.atomic_number) * cutoff_scale;
            if maximum_cutoff < candidate {
                maximum_cutoff = candidate;
            }
        }
    }
    // Every image of every neighbour within the cutoff.
    let pairs = match neighbours::within(structure, maximum_cutoff as f64, use_pbc) {
        Ok(p) => p,
        Err(e) => {
            result.message = e.to_string();
            return result;
        }
    };
    let inverse = if use_pbc {
        let v = |r: usize| DVec3::from(structure.cell_vectors[r]);
        glm::inverse_d(&DMat3::from_cols(v(0), v(1), v(2)))
    } else {
        DMat3::IDENTITY
    };
    let n_atoms = structure.atoms.len();
    let mut neighbors: Vec<Vec<Neighbour>> = vec![Vec::new(); n_atoms];
    for pair in &pairs {
        if pair.r <= MIN_BOND_DISTANCE as f64 {
            continue;
        }
        let ri = radius_of(structure.atoms[pair.i].atomic_number);
        let rj = radius_of(structure.atoms[pair.j].atomic_number);
        let cutoff = if cutoff_radius > 0.0 { cutoff_radius as f64 } else { ((ri + rj) * cutoff_scale) as f64 };
        if pair.r > cutoff {
            continue;
        }
        let mut neighbour = Neighbour { j: pair.j as i32, n: [0; 3] };
        if use_pbc {
            let (from, to) = (&structure.atoms[pair.i], &structure.atoms[pair.j]);
            let image = glm::mul_d(&inverse, pair.d - DVec3::new(to.x - from.x, to.y - from.y, to.z - from.z));
            neighbour.n = [image.x.round() as i32, image.y.round() as i32, image.z.round() as i32];
        }
        neighbors[pair.i].push(neighbour);
    }
    for list in neighbors.iter_mut() {
        list.sort();
        list.dedup();
    }
    let sub = |a: [i32; 3], b: [i32; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    // Is atom b in image n a neighbour of atom a in image 0?
    let bonded = |a: i32, b: i32, n: [i32; 3]| neighbors[a as usize].binary_search(&Neighbour { j: b, n }).is_ok();

    let mut atom_signatures: Vec<BTreeMap<Signature, i32>> = vec![BTreeMap::new(); n_atoms];
    let mut common: Vec<Neighbour> = Vec::new();
    for i in 0..n_atoms {
        for partner in &neighbors[i] {
            // Each bonded pair once: from the lower index (or, for an atom and
            // its own image, from the image with the positive offset).
            let j = partner.j as usize;
            if j < i {
                continue;
            }
            if j == i && !(Neighbour { j: i as i32, n: [0; 3] } < Neighbour { j: i as i32, n: partner.n }) {
                continue;
            }
            common.clear();
            for k in &neighbors[i] {
                if bonded(partner.j, k.j, sub(k.n, partner.n)) {
                    common.push(*k);
                }
            }
            let mut bond_count = 0;
            for x in 0..common.len() {
                for y in x + 1..common.len() {
                    if bonded(common[x].j, common[y].j, sub(common[y].n, common[x].n)) {
                        bond_count += 1;
                    }
                }
            }
            let chain = longest_chain_length(common.len(), |x, y| bonded(common[x].j, common[y].j, sub(common[y].n, common[x].n)));
            let sig = Signature { common: common.len() as i32, bonds: bond_count, chain };
            *result.signature_counts.entry(sig).or_default() += 1;
            *atom_signatures[i].entry(sig).or_default() += 1;
            *atom_signatures[j].entry(sig).or_default() += 1;
            result.pair_count += 1;
        }
    }
    for (i, atom) in structure.atoms.iter().enumerate() {
        let mut dominant = Signature::default();
        let mut dominant_count = 0;
        for (&sig, &count) in &atom_signatures[i] {
            if count > dominant_count {
                dominant = sig;
                dominant_count = count;
            }
        }
        let coordination = neighbors[i].len() as i32;
        let environment = classify_environment(&atom_signatures[i], coordination).to_string();
        *result.environment_counts.entry(environment.clone()).or_default() += 1;
        result.atom_rows.push(AtomRow {
            index: i as i32,
            atomic_number: atom.atomic_number,
            symbol: atom.symbol.clone(),
            coordination,
            dominant_signature: dominant,
            dominant_signature_count: dominant_count,
            environment,
        });
    }
    result.valid = true;
    result.message = "CNA completed.".into();
    result
}

pub fn compute_cna(structure: &Structure, params: &CnaParams) -> CnaResult {
    if !params.cutoff_scale.is_finite()
        || params.cutoff_scale <= 0.0
        || !params.cutoff_radius.is_finite()
        || params.cutoff_radius < 0.0
        || params.cutoff_radius > 100000.0
        || params.cutoff_scale > 100.0
    {
        return CnaResult { message: "Cutoff scale must be finite and positive.".into(), ..Default::default() };
    }
    if structure.atoms.iter().any(|a| {
        !a.x.is_finite() || !a.y.is_finite() || !a.z.is_finite() || a.x.abs() > 1e8 || a.y.abs() > 1e8 || a.z.abs() > 1e8
    }) {
        return CnaResult { message: "Atom coordinates must be finite and within 1e8 A.".into(), ..Default::default() };
    }
    run_cna(structure, params.cutoff_scale, params.use_pbc, params.cutoff_radius)
}
