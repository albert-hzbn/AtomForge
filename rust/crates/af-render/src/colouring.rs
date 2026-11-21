//! Atom colouring modes of View > Colour Atoms By (the C++ app/EditorOps
//! applyDisplayColors and applyGrainBoundaryColors): element colours, crystal
//! orientation (per-atom grain colours), grain boundaries and a per-atom
//! property on the viridis scale.

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use af_core::cell::{BOND_TOLERANCE_FACTOR, MIN_BOND_DISTANCE};
use af_core::elements;
use af_core::glm::{self, length};
use af_core::Structure;
use glam::{Mat3, Vec3};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ColourMode {
    #[default]
    Element,
    CrystalOrientation,
    GrainBoundary,
    AtomProperty,
}

impl ColourMode {
    pub fn label(self) -> &'static str {
        match self {
            ColourMode::Element => "Element",
            ColourMode::CrystalOrientation => "Crystal orientation",
            ColourMode::GrainBoundary => "Grain boundaries",
            ColourMode::AtomProperty => "Atom property",
        }
    }
}

pub const CRYSTAL_COLOUR: [f32; 3] = [0.20, 0.72, 0.98];
pub const BOUNDARY_COLOUR: [f32; 3] = [1.00, 0.45, 0.12];

/// Colours the atoms for a mode (applyDisplayColors): element colours first,
/// then the mode's colours over them.
pub fn apply_colour_mode(structure: &mut Structure, element_colours: &[[f32; 3]], mode: ColourMode, property_colours: Option<&[[f32; 3]]>) {
    for atom in structure.atoms.iter_mut() {
        let z = atom.atomic_number;
        if z >= 0 && (z as usize) < element_colours.len() {
            [atom.r, atom.g, atom.b] = element_colours[z as usize];
        }
    }
    match mode {
        ColourMode::Element => {}
        ColourMode::CrystalOrientation => {
            if structure.grain_colors.len() == structure.atoms.len() {
                for (atom, c) in structure.atoms.iter_mut().zip(structure.grain_colors.iter()) {
                    [atom.r, atom.g, atom.b] = [c[0], c[1], c[2]];
                }
            }
        }
        ColourMode::GrainBoundary => apply_grain_boundary_colours(structure),
        ColourMode::AtomProperty => {
            if let Some(colours) = property_colours {
                if colours.len() == structure.atoms.len() {
                    for (atom, c) in structure.atoms.iter_mut().zip(colours.iter()) {
                        [atom.r, atom.g, atom.b] = *c;
                    }
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct Signature {
    common: i32,
    bonds: i32,
    chain: i32,
}

fn edge_key(a: usize, b: usize) -> u64 {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    ((a as u32 as u64) << 32) | b as u32 as u64
}

fn grid_cell(p: Vec3) -> (i32, i32, i32) {
    let c = |x: f32| (x / 4.0).floor() as i32;
    (c(p.x), c(p.y), c(p.z))
}

/// Atoms (not bonds) in the longest shortest path among the common neighbours.
fn longest_chain_length(common: &[usize], edges: &HashSet<u64>) -> i32 {
    match common.len() {
        0 => return 0,
        1 => return 1,
        _ => {}
    }
    let n = common.len();
    let mut adjacency = vec![Vec::new(); n];
    for i in 0..n {
        for j in i + 1..n {
            if edges.contains(&edge_key(common[i], common[j])) {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
        }
    }
    let mut best = 1;
    for source in 0..n {
        let mut dist = vec![-1i32; n];
        let mut queue = VecDeque::new();
        dist[source] = 0;
        queue.push_back(source);
        while let Some(u) = queue.pop_front() {
            for &v in &adjacency[u] {
                if dist[v] >= 0 {
                    continue;
                }
                dist[v] = dist[u] + 1;
                best = best.max(dist[v] + 1);
                queue.push_back(v);
            }
        }
    }
    best
}

fn crystal_like(s: Signature) -> bool {
    let sig = |common, bonds, chain| Signature { common, bonds, chain };
    [sig(4, 2, 1), sig(4, 2, 2), sig(4, 4, 1), sig(6, 6, 1), sig(5, 5, 1)].contains(&s)
}

fn paint(structure: &mut Structure, boundary: &[bool]) {
    for (atom, &b) in structure.atoms.iter_mut().zip(boundary) {
        [atom.r, atom.g, atom.b] = if b { BOUNDARY_COLOUR } else { CRYSTAL_COLOUR };
    }
}

/// Which atoms are grain-boundary atoms (applyGrainBoundaryColors): across
/// grain regions when the structure has them, else atoms whose coordination
/// differs from their element's usual one (away from the box edges), else
/// atoms without a crystal-like common-neighbour signature.
pub fn grain_boundary_atoms(structure: &Structure) -> Vec<bool> {
    let n = structure.atoms.len();
    if n < 4 {
        return vec![false; n];
    }
    let radii = elements::literature_covalent_radii();
    let radius_of = |z: i32| if z >= 0 && (z as usize) < radii.len() { radii[z as usize] } else { 1.0 };
    let matrices = if structure.has_unit_cell { glm::cell_matrices(structure) } else { None };
    let use_pbc = matrices.is_some();
    let (cell, inv) = matrices.unwrap_or((Mat3::IDENTITY, Mat3::IDENTITY));
    let positions: Vec<Vec3> = structure.atoms.iter().map(|a| Vec3::new(a.x as f32, a.y as f32, a.z as f32)).collect();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut edges: HashSet<u64> = HashSet::new();
    let mut try_add_bond = |i: usize, j: usize| {
        if j <= i {
            return;
        }
        let delta = glm::minimum_image_delta(positions[j] - positions[i], use_pbc, &cell, &inv);
        let d = length(delta);
        if d <= MIN_BOND_DISTANCE {
            return;
        }
        if d > (radius_of(structure.atoms[i].atomic_number) + radius_of(structure.atoms[j].atomic_number)) * BOND_TOLERANCE_FACTOR {
            return;
        }
        neighbors[i].push(j);
        neighbors[j].push(i);
        edges.insert(edge_key(i, j));
    };
    if !use_pbc {
        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &p) in positions.iter().enumerate() {
            grid.entry(grid_cell(p)).or_default().push(i);
        }
        for (i, &p) in positions.iter().enumerate() {
            let c = grid_cell(p);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&(c.0 + dx, c.1 + dy, c.2 + dz)) {
                            for &j in list {
                                try_add_bond(i, j);
                            }
                        }
                    }
                }
            }
        }
    } else {
        let mut largest = 0.0f32;
        for a in &structure.atoms {
            if a.atomic_number >= 0 && (a.atomic_number as usize) < radii.len() {
                largest = largest.max(radii[a.atomic_number as usize]);
            }
        }
        largest = largest.max(1.0);
        af_core::neighbours::for_each_candidate_pair(structure, 2.0 * largest as f64 * BOND_TOLERANCE_FACTOR as f64, true, |i, j| try_add_bond(i, j));
    }
    for list in neighbors.iter_mut() {
        list.sort_unstable();
    }

    let mut boundary = vec![false; n];
    if structure.grain_region_ids.len() == n {
        for i in 0..n {
            let region = structure.grain_region_ids[i];
            for &j in &neighbors[i] {
                if structure.grain_region_ids[j] == region {
                    continue;
                }
                // Bonds across a periodic boundary join opposite faces of the
                // box (artefacts of the periodic tessellation).
                if use_pbc {
                    let raw = positions[j] - positions[i];
                    let min = glm::minimum_image_delta(raw, true, &cell, &inv);
                    if length(raw - min) > 0.1 {
                        continue;
                    }
                }
                boundary[i] = true;
                boundary[j] = true;
            }
        }
        return boundary;
    }

    let coordination: Vec<i32> = neighbors.iter().map(|l| l.len() as i32).collect();
    let mut histogram: BTreeMap<i32, BTreeMap<i32, i32>> = BTreeMap::new();
    for i in 0..n {
        *histogram.entry(structure.atoms[i].atomic_number).or_default().entry(coordination[i]).or_default() += 1;
    }
    let modal: HashMap<i32, i32> = histogram
        .iter()
        .map(|(&z, counts)| {
            let (mut best, mut best_count) = (0, -1);
            for (&c, &k) in counts {
                if k > best_count {
                    best = c;
                    best_count = k;
                }
            }
            (z, best)
        })
        .collect();

    // Atoms near the box (or bounding box) edges are surface atoms.
    let (lo, hi) = if use_pbc {
        let lo = Vec3::new(structure.cell_offset[0] as f32, structure.cell_offset[1] as f32, structure.cell_offset[2] as f32);
        (lo, lo + cell.x_axis + cell.y_axis + cell.z_axis)
    } else {
        let mut lo = Vec3::splat(f32::MAX);
        let mut hi = Vec3::splat(f32::MIN);
        for &p in &positions {
            lo = glm::min(lo, p);
            hi = glm::max(hi, p);
        }
        (lo, hi)
    };
    let mut max_radius = 0.0f32;
    for a in &structure.atoms {
        let r = radius_of(a.atomic_number);
        if r > max_radius {
            max_radius = r;
        }
    }
    let edge = max_radius * 2.0 * BOND_TOLERANCE_FACTOR;
    let near_edge = |p: Vec3| p.x - lo.x < edge || hi.x - p.x < edge || p.y - lo.y < edge || hi.y - p.y < edge || p.z - lo.z < edge || hi.z - p.z < edge;

    let mut seeds = 0;
    for i in 0..n {
        let usual = modal.get(&structure.atoms[i].atomic_number).copied().unwrap_or(coordination[i]);
        if (coordination[i] - usual).abs() >= 1 && !near_edge(positions[i]) {
            boundary[i] = true;
            seeds += 1;
        }
    }
    if seeds > 0 {
        return boundary;
    }

    let mut signatures: Vec<BTreeMap<Signature, i32>> = vec![BTreeMap::new(); n];
    for i in 0..n {
        for &j in &neighbors[i] {
            if j <= i {
                continue;
            }
            let common: Vec<usize> = intersection(&neighbors[i], &neighbors[j]);
            let mut bonds = 0;
            for a in 0..common.len() {
                for b in a + 1..common.len() {
                    if edges.contains(&edge_key(common[a], common[b])) {
                        bonds += 1;
                    }
                }
            }
            let sig = Signature { common: common.len() as i32, bonds, chain: longest_chain_length(&common, &edges) };
            *signatures[i].entry(sig).or_default() += 1;
            *signatures[j].entry(sig).or_default() += 1;
        }
    }
    for i in 0..n {
        let (mut dominant, mut count) = (Signature::default(), 0);
        for (&sig, &k) in &signatures[i] {
            if k > count {
                dominant = sig;
                count = k;
            }
        }
        boundary[i] = !(count > 0 && crystal_like(dominant));
    }
    boundary
}

fn intersection(a: &[usize], b: &[usize]) -> Vec<usize> {
    let (mut i, mut j) = (0, 0);
    let mut out = Vec::new();
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

/// Grain-boundary atoms orange, the rest blue.
pub fn apply_grain_boundary_colours(structure: &mut Structure) {
    if structure.atoms.is_empty() {
        return;
    }
    let boundary = grain_boundary_atoms(structure);
    paint(structure, &boundary);
}
