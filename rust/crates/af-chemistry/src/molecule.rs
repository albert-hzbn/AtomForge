//! Molecular graphs (the C++ chemistry/Molecule): atoms with formal charges
//! and hydrogen counts, bonds with orders, stereo descriptors, and the
//! perception that recovers them from coordinates (bonds, bond orders, rings,
//! aromaticity, hybridisation, CIP R/S and E/Z labels).

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::sync::OnceLock;

use af_core::cxx::sort_by;
use af_core::elements;
use af_core::structure::{AtomSite, Structure};

pub type Vec3 = [f64; 3];

const PI: f64 = 3.14159265358979323846;

#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    pub z: i32,
    /// Formal charge.
    pub charge: i32,
    /// Implicit hydrogens (not yet explicit atoms).
    pub hydrogens: i32,
    /// Aromatic as written (SMILES) or perceived.
    pub aromatic: bool,
    pub isotope: i32,
    /// SMILES atom class ([C:1]); 0 = none.
    pub map_number: i32,
    /// '*': an attachment point.
    pub dummy: bool,
    pub position: Vec3,
}

impl Default for Atom {
    fn default() -> Self {
        Self { z: 6, charge: 0, hydrogens: 0, aromatic: false, isotope: 0, map_number: 0, dummy: false, position: [0.0; 3] }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    pub a: i32,
    pub b: i32,
    /// 1, 2, 3; aromatic bonds have their Kekule order.
    pub order: i32,
    pub aromatic: bool,
    /// To a metal: not counted in valences.
    pub dative: bool,
    /// The vector from a to b is position[b] - position[a] + image (a periodic
    /// bond across the cell); usually zero.
    pub image: Vec3,
}

impl Default for Bond {
    fn default() -> Self {
        Self { a: 0, b: 0, order: 1, aromatic: false, dative: false, image: [0.0; 3] }
    }
}

/// A tetrahedral stereocentre: sign of the triple product of the vectors from
/// the centre to neighbours n[1], n[2], n[3] (+1 is SMILES '@@' with n[0] first).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TetrahedralStereo {
    pub centre: i32,
    pub n: [i32; 4],
    pub sign: i32,
}

/// A double bond x=y with a neighbour a of x and b of y on the same side (cis) or not.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DoubleBondStereo {
    pub a: i32,
    pub x: i32,
    pub y: i32,
    pub b: i32,
    pub cis: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub tetrahedral: Vec<TetrahedralStereo>,
    pub double_bonds: Vec<DoubleBondStereo>,
    atom_bonds: Vec<Vec<i32>>,
}

pub(crate) fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(crate) fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(crate) fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(crate) fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(crate) fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub(crate) fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}
pub(crate) fn unit(a: Vec3) -> Vec3 {
    let n = norm(a);
    if n > 1e-12 {
        scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 1.0]
    }
}
/// Any unit vector perpendicular to a.
pub(crate) fn perpendicular(a: Vec3) -> Vec3 {
    let trial = if a[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    unit(cross(a, trial))
}
pub(crate) fn fmin(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}
pub(crate) fn fmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

/// Group (1, 2, 13-18; 3-12 for the d-block, 0 for the f-block).
pub(crate) fn group_of(z: i32) -> i32 {
    if z == 1 {
        return 1;
    }
    if z == 2 {
        return 18;
    }
    const STARTS: [i32; 6] = [3, 11, 19, 37, 55, 87];
    for p in (0..6).rev() {
        let s = STARTS[p];
        if z < s {
            continue;
        }
        let k = z - s; // position within the period
        if p <= 1 {
            return if k < 2 { k + 1 } else { k + 11 }; // periods 2 and 3: groups 1, 2, 13-18
        }
        if p <= 3 {
            return k + 1; // periods 4 and 5: groups 1-18
        }
        // Periods 6 and 7: groups 1, 2, then 14 f-block elements, then 3-18.
        if k < 2 {
            return k + 1;
        }
        if k < 17 {
            return if k == 2 { 3 } else { 0 };
        }
        return k - 13;
    }
    0
}

pub(crate) fn period_of(z: i32) -> i32 {
    const ENDS: [i32; 7] = [2, 10, 18, 36, 54, 86, 118];
    for (p, &end) in ENDS.iter().enumerate() {
        if z <= end {
            return p as i32 + 1;
        }
    }
    7
}

fn main_group(z: i32) -> bool {
    let g = group_of(z);
    g == 1 || g == 2 || g >= 13
}

fn valences_of_group(group: i32, period: i32) -> &'static [i32] {
    match group {
        1 => &[1],
        2 => &[2],
        13 => &[3],
        14 => &[4],
        15 => {
            if period <= 2 {
                &[3]
            } else {
                &[3, 5]
            }
        }
        16 => {
            if period <= 2 {
                &[2]
            } else {
                &[2, 4, 6]
            }
        }
        17 => {
            if period <= 2 {
                &[1]
            } else {
                &[1, 3, 5, 7]
            }
        }
        18 => {
            if period <= 2 {
                &[0]
            } else {
                &[0, 2, 4, 6]
            }
        }
        _ => &[],
    }
}

/// Maximum-cardinality matching in a general graph (Edmonds' blossom
/// algorithm), extending an initial matching; matched vertices stay matched.
/// Only vertices with start[v] begin augmenting paths.
pub(crate) fn augment_matching(adjacency: &[Vec<i32>], matching: &mut [i32], start: &[u8]) {
    let n = adjacency.len();
    let mut parent = vec![-1i32; n];
    let mut base: Vec<i32> = (0..n as i32).collect();
    let mut used = vec![0u8; n];
    let mut blossom = vec![0u8; n];
    let mut queue: Vec<i32> = Vec::new();
    let lca = |a0: i32, b0: i32, base: &[i32], parent: &[i32], matching: &[i32]| -> i32 {
        let mut seen = vec![0u8; n];
        let mut a = a0;
        loop {
            a = base[a as usize];
            seen[a as usize] = 1;
            if matching[a as usize] == -1 {
                break;
            }
            a = parent[matching[a as usize] as usize];
        }
        let mut b = b0;
        loop {
            b = base[b as usize];
            if seen[b as usize] != 0 {
                return b;
            }
            b = parent[matching[b as usize] as usize];
        }
    };
    let mark_path = |mut v: i32, b: i32, mut child: i32, base: &[i32], parent: &mut [i32], blossom: &mut [u8], matching: &[i32]| {
        while base[v as usize] != b {
            let mv = matching[v as usize];
            blossom[base[v as usize] as usize] = 1;
            blossom[base[mv as usize] as usize] = 1;
            parent[v as usize] = child;
            child = mv;
            v = parent[mv as usize];
        }
    };
    for v in 0..n {
        if matching[v] != -1 || start[v] == 0 || adjacency[v].is_empty() {
            continue;
        }
        // findPath(v)
        let root = v as i32;
        used.iter_mut().for_each(|u| *u = 0);
        parent.iter_mut().for_each(|p| *p = -1);
        for (i, b) in base.iter_mut().enumerate() {
            *b = i as i32;
        }
        used[v] = 1;
        queue.clear();
        queue.push(root);
        let mut found = -1;
        let mut head = 0;
        'search: while head < queue.len() {
            let vv = queue[head];
            head += 1;
            for &to in &adjacency[vv as usize] {
                if base[vv as usize] == base[to as usize] || matching[vv as usize] == to {
                    continue;
                }
                if to == root || (matching[to as usize] != -1 && parent[matching[to as usize] as usize] != -1) {
                    let current = lca(vv, to, &base, &parent, matching);
                    blossom.iter_mut().for_each(|x| *x = 0);
                    mark_path(vv, current, to, &base, &mut parent, &mut blossom, matching);
                    mark_path(to, current, vv, &base, &mut parent, &mut blossom, matching);
                    for i in 0..n {
                        if blossom[base[i] as usize] != 0 {
                            base[i] = current;
                            if used[i] == 0 {
                                used[i] = 1;
                                queue.push(i as i32);
                            }
                        }
                    }
                } else if parent[to as usize] == -1 {
                    parent[to as usize] = vv;
                    if matching[to as usize] == -1 {
                        found = to;
                        break 'search;
                    }
                    used[matching[to as usize] as usize] = 1;
                    queue.push(matching[to as usize]);
                }
            }
        }
        let mut u = found;
        while u != -1 {
            let pv = parent[u as usize];
            let next = matching[pv as usize];
            matching[u as usize] = pv;
            matching[pv as usize] = u;
            u = next;
        }
    }
}

fn radii() -> &'static [f32] {
    static TABLE: OnceLock<Vec<f32>> = OnceLock::new();
    TABLE.get_or_init(elements::literature_covalent_radii)
}

impl Molecule {
    pub fn add_atom(&mut self, atom: Atom) -> i32 {
        self.atoms.push(atom);
        self.atom_bonds.push(Vec::new());
        self.atoms.len() as i32 - 1
    }

    pub fn add_bond(&mut self, a: i32, b: i32, order: i32) -> Result<i32, String> {
        if a == b {
            return Err("a bond needs two different atoms".into());
        }
        if self.atom_bonds.len() < self.atoms.len() {
            self.atom_bonds.resize(self.atoms.len(), Vec::new());
        }
        self.bonds.push(Bond { a, b, order, ..Default::default() });
        let index = self.bonds.len() as i32 - 1;
        self.atom_bonds[a as usize].push(index);
        self.atom_bonds[b as usize].push(index);
        Ok(index)
    }

    /// Rebuilds the per-atom bond lists after bonds were edited directly.
    pub fn reindex(&mut self) {
        self.atom_bonds = vec![Vec::new(); self.atoms.len()];
        for (k, bond) in self.bonds.iter().enumerate() {
            self.atom_bonds[bond.a as usize].push(k as i32);
            self.atom_bonds[bond.b as usize].push(k as i32);
        }
    }

    /// Bond index between a and b, or -1.
    pub fn bond_between(&self, a: i32, b: i32) -> i32 {
        self.bonds_of(a).iter().copied().find(|&k| self.other(k, a) == b).unwrap_or(-1)
    }

    pub fn other(&self, bond: i32, atom: i32) -> i32 {
        let b = &self.bonds[bond as usize];
        if b.a == atom {
            b.b
        } else {
            b.a
        }
    }

    /// Bond indices of an atom, in the order they were added.
    pub fn bonds_of(&self, atom: i32) -> &[i32] {
        if atom >= 0 && (atom as usize) < self.atom_bonds.len() {
            &self.atom_bonds[atom as usize]
        } else {
            &[]
        }
    }

    pub fn neighbours(&self, atom: i32) -> Vec<i32> {
        self.bonds_of(atom).iter().map(|&k| self.other(k, atom)).collect()
    }

    pub fn degree(&self, atom: i32) -> i32 {
        self.bonds_of(atom).len() as i32
    }

    /// Sum of bond orders (dative bonds excluded) plus implicit hydrogens.
    pub fn valence(&self, atom: i32) -> i32 {
        let mut sum = self.atoms[atom as usize].hydrogens;
        for &k in self.bonds_of(atom) {
            if !self.bonds[k as usize].dative {
                sum += self.bonds[k as usize].order;
            }
        }
        sum
    }

    /// Explicit H neighbours plus implicit hydrogens.
    pub fn hydrogen_count(&self, atom: i32) -> i32 {
        self.atoms[atom as usize].hydrogens + self.bonds_of(atom).iter().filter(|&&k| self.atoms[self.other(k, atom) as usize].z == 1).count() as i32
    }

    pub fn heavy_degree(&self, atom: i32) -> i32 {
        self.bonds_of(atom).iter().filter(|&&k| self.atoms[self.other(k, atom) as usize].z != 1).count() as i32
    }

    /// Removes atoms (keep[i] == 0) with their bonds and stereo; returns the old -> new map (-1 removed).
    pub fn remove_atoms(&mut self, keep: &[u8]) -> Vec<i32> {
        let mut map = vec![-1i32; self.atoms.len()];
        let mut kept = Vec::new();
        for (i, atom) in self.atoms.iter().enumerate() {
            if keep[i] != 0 {
                map[i] = kept.len() as i32;
                kept.push(atom.clone());
            }
        }
        let mut kept_bonds = Vec::new();
        for b in &self.bonds {
            if map[b.a as usize] >= 0 && map[b.b as usize] >= 0 {
                let mut c = b.clone();
                c.a = map[b.a as usize];
                c.b = map[b.b as usize];
                kept_bonds.push(c);
            }
        }
        let mut kept_tetra = Vec::new();
        for t in &self.tetrahedral {
            let mut c = *t;
            let mut ok = map[t.centre as usize] >= 0;
            c.centre = if ok { map[t.centre as usize] } else { 0 };
            for n in c.n.iter_mut() {
                if *n < 0 {
                    continue; // an implicit hydrogen
                }
                if map[*n as usize] < 0 {
                    ok = false;
                } else {
                    *n = map[*n as usize];
                }
            }
            if ok {
                kept_tetra.push(c);
            }
        }
        let mut kept_double = Vec::new();
        for d in &self.double_bonds {
            if map[d.a as usize] >= 0 && map[d.x as usize] >= 0 && map[d.y as usize] >= 0 && map[d.b as usize] >= 0 {
                kept_double.push(DoubleBondStereo { a: map[d.a as usize], x: map[d.x as usize], y: map[d.y as usize], b: map[d.b as usize], cis: d.cis });
            }
        }
        self.atoms = kept;
        self.bonds = kept_bonds;
        self.tetrahedral = kept_tetra;
        self.double_bonds = kept_double;
        self.reindex();
        map
    }

    /// The vector of bond k from atom `from` to its other atom (image included).
    pub(crate) fn bond_vector(&self, k: i32, from: i32) -> Vec3 {
        let b = &self.bonds[k as usize];
        let d = add(sub(self.atoms[b.b as usize].position, self.atoms[b.a as usize].position), b.image);
        if b.a == from {
            d
        } else {
            scale(d, -1.0)
        }
    }
}

// ------------------------------------------------------------------ elements

/// Atomic number of a symbol (any case; D and T as hydrogen), 0 when unknown.
pub fn atomic_number(symbol: &str) -> i32 {
    static TABLE: OnceLock<HashMap<String, i32>> = OnceLock::new();
    let table = TABLE.get_or_init(|| {
        let mut t: HashMap<String, i32> = (1..=118).map(|z| (elements::symbol(z).to_string(), z)).collect();
        t.insert("D".into(), 1);
        t.insert("T".into(), 1);
        t
    });
    if symbol.is_empty() {
        return 0;
    }
    let mut s: Vec<u8> = symbol.bytes().collect();
    s[0] = s[0].to_ascii_uppercase();
    for c in s.iter_mut().skip(1) {
        *c = c.to_ascii_lowercase();
    }
    table.get(&String::from_utf8_lossy(&s).into_owned()).copied().unwrap_or(0)
}

pub fn symbol_of(z: i32) -> String {
    if z == 0 {
        "*".into()
    } else {
        elements::symbol(z).to_string()
    }
}

pub fn is_metal(z: i32) -> bool {
    if z <= 0 || z == 1 {
        return false;
    }
    let g = group_of(z);
    if g == 1 || g == 2 || (3..=12).contains(&g) || g == 0 {
        return true;
    }
    z == 13 || z == 31 || z == 49 || z == 81 // Al, Ga, In, Tl
}

/// The usual valence of an atom with the given bond-order sum (the smallest
/// allowed valence >= it, e.g. S 2, 4, 6); -1 for no fixed valence.
pub fn default_valence(z: i32, charge: i32, bond_order_sum: i32) -> i32 {
    if z <= 0 {
        return -1;
    }
    if z == 1 {
        return if charge == 0 { 1 } else { 0 };
    }
    if !main_group(z) || is_metal(z) {
        return -1;
    }
    let mut group = group_of(z);
    if group == 1 || group == 2 {
        return if charge == 0 { group } else { -1 };
    }
    // Isoelectronic shift: N+ behaves as C, O- as F, C- as N, B- as C, C+ as B.
    if group == 14 && charge != 0 {
        group = 13 + if charge < 0 { 2 } else { 0 };
    } else {
        group -= charge;
    }
    if !(13..=18).contains(&group) {
        return -1;
    }
    for &v in valences_of_group(group, period_of(z)) {
        if v >= bond_order_sum {
            return v;
        }
    }
    bond_order_sum // hypervalent beyond the table
}

pub fn covalent_radius(z: i32) -> f64 {
    let r = radii();
    if z > 0 && (z as usize) < r.len() && r[z as usize] > 0.0 {
        f64::from(r[z as usize])
    } else {
        1.5
    }
}

/// Bond length to hydrogen (Angstrom).
pub fn hydrogen_bond_length(z: i32) -> f64 {
    match z {
        1 => 0.74,
        5 => 1.19,
        6 => 1.09,
        7 => 1.01,
        8 => 0.96,
        9 => 0.92,
        14 => 1.48,
        15 => 1.42,
        16 => 1.34,
        17 => 1.27,
        35 => 1.41,
        53 => 1.61,
        _ => covalent_radius(z) + 0.31,
    }
}

// ------------------------------------------------------------------ perception

/// Minimum-image vector from atom i to atom j of a (periodic) structure.
pub fn bond_vector(s: &Structure, i: usize, j: usize) -> Vec3 {
    let (ai, aj) = (&s.atoms[i], &s.atoms[j]);
    let d = [aj.x - ai.x, aj.y - ai.y, aj.z - ai.z];
    let periodic = s.has_unit_cell && s.pbc.iter().any(|&p| p);
    if !periodic {
        return d;
    }
    // Fractional coordinates through the inverse cell (rows are a, b, c).
    let c = &s.cell_vectors;
    let m = [[c[0][0], c[1][0], c[2][0]], [c[0][1], c[1][1], c[2][1]], [c[0][2], c[1][2], c[2][2]]];
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-12 {
        return d;
    }
    let inv = [
        [(m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det, (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det, (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det],
        [(m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det, (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det, (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det],
        [(m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det, (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det, (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det],
    ];
    let mut f = [0.0; 3];
    for r in 0..3 {
        f[r] = inv[r][0] * d[0] + inv[r][1] * d[1] + inv[r][2] * d[2];
    }
    for r in 0..3 {
        if s.pbc[r] {
            f[r] -= f[r].round();
        }
    }
    // The nearest image may also be a neighbouring combination in skewed cells.
    let mut best = [0.0; 3];
    let mut best_norm = 1e300;
    for a in -1..=1 {
        for b in -1..=1 {
            for e in -1..=1 {
                if (a != 0 && !s.pbc[0]) || (b != 0 && !s.pbc[1]) || (e != 0 && !s.pbc[2]) {
                    continue;
                }
                let g = [f[0] + a as f64, f[1] + b as f64, f[2] + e as f64];
                let mut v = [0.0; 3];
                for r in 0..3 {
                    v[r] = m[r][0] * g[0] + m[r][1] * g[1] + m[r][2] * g[2];
                }
                let n = dot(v, v);
                if n < best_norm {
                    best_norm = n;
                    best = v;
                }
            }
        }
    }
    best
}

fn angle_sum(m: &Molecule, atom: i32) -> f64 {
    let nb = m.bonds_of(atom);
    let mut sum = 0.0;
    for x in 0..nb.len() {
        for y in x + 1..nb.len() {
            let u = unit(m.bond_vector(nb[x], atom));
            let v = unit(m.bond_vector(nb[y], atom));
            sum += dot(u, v).clamp(-1.0, 1.0).acos() * 180.0 / PI;
        }
    }
    sum
}

fn bond_length(m: &Molecule, k: usize) -> f64 {
    let b = &m.bonds[k];
    norm(add(sub(m.atoms[b.b as usize].position, m.atoms[b.a as usize].position), b.image))
}

fn assign_bond_orders(m: &mut Molecule) {
    let n = m.atoms.len();
    let mut deficit = vec![0i32; n];
    // Formal charges of over-coordinated first-row atoms (ammonium N, oxonium O, borate B).
    for i in 0..n {
        let z = m.atoms[i].z;
        if is_metal(z) || z == 1 {
            continue;
        }
        let ii = i as i32;
        let sigma = m.valence(ii);
        let terminal = |m: &Molecule, zs: &[i32]| m.bonds_of(ii).iter().any(|&k| {
            let o = m.other(k, ii);
            zs.contains(&m.atoms[o as usize].z) && m.degree(o) == 1
        });
        if z == 7 && sigma == 4 {
            m.atoms[i].charge = 1;
        } else if z == 15 && sigma == 4 {
            if !terminal(m, &[8, 16]) {
                m.atoms[i].charge = 1; // phosphonium; phosphates keep P(V)
            }
        } else if z == 8 && sigma == 3 {
            m.atoms[i].charge = 1;
        } else if z == 5 && sigma == 4 {
            m.atoms[i].charge = -1;
        } else if z == 7 && sigma == 3 && m.degree(ii) == 3 {
            // Planar N with a terminal O: nitro groups and N-oxides (N+).
            if terminal(m, &[8]) && angle_sum(m, ii) > 350.0 {
                m.atoms[i].charge = 1;
            }
        }
    }
    for i in 0..n {
        let (z, charge) = (m.atoms[i].z, m.atoms[i].charge);
        if is_metal(z) || z == 1 {
            continue;
        }
        let ii = i as i32;
        let sigma = m.valence(ii);
        let target = default_valence(z, charge, sigma);
        if target > sigma {
            deficit[i] = target - sigma;
        }
        // A tetrahedral (pyramidal) atom takes no pi bond.
        if deficit[i] > 0 && m.degree(ii) == 4 && z == 6 {
            deficit[i] = 0;
        }
        if deficit[i] > 0 && m.degree(ii) == 3 && (z == 6 || z == 5) && angle_sum(m, ii) < 340.0 {
            deficit[i] = 0;
        }
    }
    // Rounds of maximum matching over bonds between atoms that still lack valence.
    for _ in 0..3 {
        let mut adjacency = vec![Vec::new(); n];
        let mut candidates: Vec<(f64, usize)> = Vec::new();
        for (k, b) in m.bonds.iter().enumerate() {
            if b.dative || deficit[b.a as usize] <= 0 || deficit[b.b as usize] <= 0 {
                continue;
            }
            let ratio = bond_length(m, k) / (covalent_radius(m.atoms[b.a as usize].z) + covalent_radius(m.atoms[b.b as usize].z));
            // Multiple bonds are shorter than single ones: a long bond stays single.
            if ratio > if b.order == 1 { 0.97 } else { 0.86 } {
                continue;
            }
            candidates.push((ratio, k));
        }
        if candidates.is_empty() {
            break;
        }
        candidates.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        let mut matching = vec![-1i32; n];
        // Shortest (relative) bonds first, then augmenting paths.
        for &(_, k) in &candidates {
            let (a, b) = (m.bonds[k].a, m.bonds[k].b);
            adjacency[a as usize].push(b);
            adjacency[b as usize].push(a);
            if matching[a as usize] == -1 && matching[b as usize] == -1 {
                matching[a as usize] = b;
                matching[b as usize] = a;
            }
        }
        augment_matching(&adjacency, &mut matching, &vec![1u8; n]);
        let mut changed = false;
        for i in 0..n {
            if matching[i] > i as i32 {
                let k = m.bond_between(i as i32, matching[i]) as usize;
                m.bonds[k].order += 1;
                deficit[i] -= 1;
                deficit[matching[i] as usize] -= 1;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    // A carbon left short next to a planar three-connected N: an iminium N+=C.
    for i in 0..n {
        if deficit[i] <= 0 || m.atoms[i].z != 6 {
            continue;
        }
        let ii = i as i32;
        let mut best = -1i32;
        let mut best_length = 1e300;
        for &k in m.bonds_of(ii) {
            let x = m.other(k, ii);
            let bond = &m.bonds[k as usize];
            let xa = &m.atoms[x as usize];
            if bond.dative || bond.order != 1 || xa.z != 7 || xa.charge != 0 {
                continue;
            }
            if m.valence(x) != 3 || m.degree(x) != 3 || angle_sum(m, x) < 350.0 {
                continue;
            }
            let length = bond_length(m, k as usize);
            if length < best_length {
                best_length = length;
                best = k;
            }
        }
        if best < 0 {
            continue;
        }
        m.bonds[best as usize].order += 1;
        deficit[i] -= 1;
        let o = m.other(best, ii);
        m.atoms[o as usize].charge = 1;
    }
    // Atoms left short: anions (O-, S-, N-) or radicals.
    for i in 0..n {
        if deficit[i] <= 0 {
            continue;
        }
        let z = m.atoms[i].z;
        let g = group_of(z);
        let ii = i as i32;
        if (15..=17).contains(&g) && m.valence(ii) < (if g == 15 { 3 } else if g == 16 { 2 } else { 1 }) + 1 {
            m.atoms[i].charge -= deficit[i];
        } else if z == 6 {
            let next_to_metal = m.bonds_of(ii).iter().any(|&k| m.bonds[k as usize].dative);
            if next_to_metal {
                m.atoms[i].charge -= deficit[i];
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct PerceiveOptions {
    /// Bonded when closer than this times the sum of covalent radii.
    pub bond_scale: f64,
    /// Assign double/triple bonds and formal charges.
    pub bond_orders: bool,
}

impl Default for PerceiveOptions {
    fn default() -> Self {
        Self { bond_scale: 1.2, bond_orders: true }
    }
}

/// Bonds (minimum image in periodic cells), bond orders and formal charges,
/// then aromaticity, from a structure's coordinates.
pub fn perceive(s: &Structure, options: &PerceiveOptions) -> Result<Molecule, String> {
    let mut m = Molecule::default();
    let n = s.atoms.len();
    for site in &s.atoms {
        m.add_atom(Atom { z: if site.atomic_number > 0 { site.atomic_number } else { atomic_number(&site.symbol) }, position: [site.x, site.y, site.z], ..Default::default() });
    }
    let mut max_radius = 0.5f64;
    for a in &m.atoms {
        max_radius = fmax(max_radius, covalent_radius(a.z));
    }
    let periodic = s.has_unit_cell && s.pbc.iter().any(|&p| p);
    let mut pairs = Vec::new();
    af_core::neighbours::for_each_candidate_pair(s, 2.0 * max_radius * options.bond_scale, periodic, |i, j| pairs.push((i, j)));
    pairs.sort();
    let mut vectors = Vec::new();
    for (i, j) in pairs {
        let (zi, zj) = (m.atoms[i].z, m.atoms[j].z);
        if zi <= 0 || zj <= 0 {
            continue;
        }
        let d = bond_vector(s, i, j);
        let r = norm(d);
        if r < 0.1 || r > options.bond_scale * (covalent_radius(zi) + covalent_radius(zj)) {
            continue;
        }
        if is_metal(zi) && is_metal(zj) {
            continue; // metallic contacts are not bonds of a molecule
        }
        let k = m.add_bond(i as i32, j as i32, 1)? as usize;
        m.bonds[k].dative = is_metal(zi) || is_metal(zj);
        vectors.push(d);
    }
    // H bonded to more than one atom keeps its shortest bond only.
    {
        let mut keep_bond = vec![true; m.bonds.len()];
        for i in 0..n as i32 {
            if m.atoms[i as usize].z != 1 || m.degree(i) <= 1 {
                continue;
            }
            let mut best = -1i32;
            for &k in m.bonds_of(i) {
                if best < 0 || norm(vectors[k as usize]) < norm(vectors[best as usize]) {
                    best = k;
                }
            }
            for &k in m.bonds_of(i) {
                if k != best && m.atoms[m.other(k, i) as usize].z != 5 {
                    keep_bond[k as usize] = false;
                }
            }
        }
        let (mut kept, mut kept_vectors) = (Vec::new(), Vec::new());
        for (k, bond) in m.bonds.iter().enumerate() {
            if keep_bond[k] {
                kept.push(bond.clone());
                kept_vectors.push(vectors[k]);
            }
        }
        m.bonds = kept;
        vectors = kept_vectors;
        m.reindex();
    }
    // Unwrap each molecule along the bond vectors and record the bonds that
    // still cross the cell (infinite chains and networks).
    if periodic {
        let mut placed = vec![false; n];
        for root in 0..n {
            if placed[root] {
                continue;
            }
            placed[root] = true;
            let mut queue = VecDeque::from([root as i32]);
            while let Some(a) = queue.pop_front() {
                for &k in m.bonds_of(a).to_vec().iter() {
                    let b = m.other(k, a);
                    if placed[b as usize] {
                        continue;
                    }
                    let d = if m.bonds[k as usize].a == a { vectors[k as usize] } else { scale(vectors[k as usize], -1.0) };
                    m.atoms[b as usize].position = add(m.atoms[a as usize].position, d);
                    placed[b as usize] = true;
                    queue.push_back(b);
                }
            }
        }
        for k in 0..m.bonds.len() {
            let (a, b) = (m.bonds[k].a as usize, m.bonds[k].b as usize);
            let image = sub(vectors[k], sub(m.atoms[b].position, m.atoms[a].position));
            m.bonds[k].image = if norm(image) > 1e-6 { image } else { [0.0; 3] };
        }
    }
    if options.bond_orders {
        assign_bond_orders(&mut m);
        perceive_aromaticity(&mut m);
    }
    Ok(m)
}

/// The atoms of a molecule as a structure (symbols, positions, colours), with
/// formal charges in the per-atom charge field when `charges_to_field` is set.
pub fn to_structure(m: &Molecule, charges_to_field: bool) -> Structure {
    let mut s = Structure::default();
    for a in &m.atoms {
        let [r, g, b] = elements::default_colour(a.z);
        s.atoms.push(AtomSite {
            atomic_number: a.z,
            symbol: if a.dummy || a.z == 0 { "X".into() } else { elements::symbol(a.z).to_string() },
            x: a.position[0],
            y: a.position[1],
            z: a.position[2],
            r,
            g,
            b,
            charge: if charges_to_field { f64::from(a.charge) } else { 0.0 },
            ..Default::default()
        });
    }
    s.has_unit_cell = false;
    s.pbc = [false; 3];
    s
}

// ------------------------------------------------------------------ rings

/// Smallest set of smallest rings (atom cycles in ring order).
pub fn smallest_rings(m: &Molecule) -> Vec<Vec<i32>> {
    let n = m.atoms.len();
    let nb = m.bonds.len();
    // Number of independent rings: E - V + components.
    let mut component = vec![-1i32; n];
    let mut components = 0;
    for i in 0..n {
        if component[i] >= 0 {
            continue;
        }
        let mut stack = vec![i as i32];
        component[i] = components;
        while let Some(a) = stack.pop() {
            for b in m.neighbours(a) {
                if component[b as usize] < 0 {
                    component[b as usize] = components;
                    stack.push(b);
                }
            }
        }
        components += 1;
    }
    let wanted = nb as i64 - n as i64 + i64::from(components);
    if wanted <= 0 {
        return Vec::new();
    }
    // Candidate rings: for every bond, the shortest path between its atoms that avoids it.
    struct Candidate {
        atoms: Vec<i32>,
        edges: Vec<u64>,
    }
    let mut candidates: Vec<Candidate> = Vec::new();
    let mut seen: BTreeSet<Vec<i32>> = BTreeSet::new();
    let words = nb.div_ceil(64);
    for k in 0..nb {
        let (from, to) = (m.bonds[k].a, m.bonds[k].b);
        let mut previous = vec![-2i32; n];
        let mut previous_bond = vec![-1i32; n];
        let mut queue = VecDeque::from([from]);
        previous[from as usize] = -1;
        while previous[to as usize] == -2 {
            let Some(a) = queue.pop_front() else { break };
            for &e in m.bonds_of(a) {
                if e as usize == k {
                    continue;
                }
                let b = m.other(e, a);
                if previous[b as usize] != -2 {
                    continue;
                }
                previous[b as usize] = a;
                previous_bond[b as usize] = e;
                queue.push_back(b);
            }
        }
        if previous[to as usize] == -2 {
            continue; // a bridge
        }
        let mut c = Candidate { atoms: Vec::new(), edges: vec![0; words] };
        c.edges[k / 64] |= 1u64 << (k % 64);
        let mut a = to;
        while a != from {
            c.atoms.push(a);
            let e = previous_bond[a as usize] as usize;
            c.edges[e / 64] |= 1u64 << (e % 64);
            a = previous[a as usize];
        }
        c.atoms.push(from);
        let mut key = c.atoms.clone();
        key.sort();
        if seen.insert(key) {
            candidates.push(c);
        }
    }
    candidates.sort_by_key(|c| c.atoms.len()); // stable, as std::stable_sort
    // Keep linearly independent rings (over GF(2) on the bond sets), smallest first.
    let mut basis: Vec<Vec<u64>> = Vec::new();
    let mut pivots: Vec<usize> = Vec::new();
    let mut rings = Vec::new();
    for c in &candidates {
        let mut v = c.edges.clone();
        for (r, row) in basis.iter().enumerate() {
            if (v[pivots[r] / 64] >> (pivots[r] % 64)) & 1 != 0 {
                for w in 0..words {
                    v[w] ^= row[w];
                }
            }
        }
        let mut pivot = None;
        'find: for w in 0..words {
            if v[w] != 0 {
                for bit in 0..64 {
                    if (v[w] >> bit) & 1 != 0 {
                        pivot = Some(w * 64 + bit);
                        break 'find;
                    }
                }
            }
        }
        let Some(pivot) = pivot else { continue };
        // Keep the basis reduced so the pivot test above stays valid.
        for row in basis.iter_mut() {
            if (row[pivot / 64] >> (pivot % 64)) & 1 != 0 {
                for w in 0..words {
                    row[w] ^= v[w];
                }
            }
        }
        basis.push(v);
        pivots.push(pivot);
        rings.push(c.atoms.clone());
        if rings.len() as i64 == wanted {
            break;
        }
    }
    rings
}

/// Per atom: size of the smallest ring through it (0 = not in a ring).
pub fn smallest_ring_size(m: &Molecule, rings: &[Vec<i32>]) -> Vec<i32> {
    let mut size = vec![0i32; m.atoms.len()];
    for ring in rings {
        for &a in ring {
            let r = ring.len() as i32;
            if size[a as usize] == 0 || r < size[a as usize] {
                size[a as usize] = r;
            }
        }
    }
    size
}

pub fn bond_in_ring(m: &Molecule, k: i32) -> bool {
    let (from, to) = (m.bonds[k as usize].a, m.bonds[k as usize].b);
    let mut seen = vec![false; m.atoms.len()];
    let mut stack = vec![from];
    seen[from as usize] = true;
    while let Some(a) = stack.pop() {
        for &e in m.bonds_of(a) {
            if e == k {
                continue;
            }
            let b = m.other(e, a);
            if b == to {
                return true;
            }
            if !seen[b as usize] {
                seen[b as usize] = true;
                stack.push(b);
            }
        }
    }
    false
}

/// Marks aromatic atoms and bonds (Hueckel 4n+2 on rings and fused ring pairs).
pub fn perceive_aromaticity(m: &mut Molecule) {
    for a in m.atoms.iter_mut() {
        a.aromatic = false;
    }
    for b in m.bonds.iter_mut() {
        b.aromatic = false;
    }
    let rings = smallest_rings(m);
    if rings.is_empty() {
        return;
    }
    // Ring systems to test: every ring and every pair of rings sharing a bond.
    let mut systems: Vec<Vec<i32>> = rings.iter().filter(|r| r.len() <= 8).cloned().collect();
    for x in 0..rings.len() {
        for y in x + 1..rings.len() {
            if rings[x].len() > 8 || rings[y].len() > 8 {
                continue;
            }
            let mut a: BTreeSet<i32> = rings[x].iter().copied().collect();
            let shared = rings[y].iter().filter(|v| a.contains(v)).count();
            if shared < 2 {
                continue;
            }
            a.extend(rings[y].iter().copied());
            systems.push(a.into_iter().collect());
        }
    }
    for _ in 0..4 {
        let mut changed = false;
        for system in &systems {
            let inside_set: BTreeSet<i32> = system.iter().copied().collect();
            if system.iter().all(|&a| m.atoms[a as usize].aromatic) {
                continue;
            }
            let mut electrons = 0;
            let mut ok = true;
            for &a in system {
                let atom = &m.atoms[a as usize];
                let z = atom.z;
                if ![5, 6, 7, 8, 15, 16, 33, 34].contains(&z) {
                    ok = false;
                    break;
                }
                let (mut inside, mut exo_double, mut exo_aromatic) = (0, 0, 0);
                for &k in m.bonds_of(a) {
                    let b = m.other(k, a);
                    let order = m.bonds[k as usize].order;
                    if order == 2 {
                        if inside_set.contains(&b) {
                            inside += 1;
                        } else if m.atoms[b as usize].aromatic {
                            exo_aromatic += 1;
                        } else {
                            exo_double += 1;
                        }
                    } else if order == 3 {
                        ok = false;
                    }
                }
                if !ok {
                    break;
                }
                if inside == 1 || exo_aromatic > 0 {
                    electrons += 1;
                } else if inside > 1 {
                    ok = false;
                    break;
                } else if exo_double > 0 {
                    // e.g. the C of a ring C=O
                } else if (z == 7 || z == 15 || z == 33) && m.valence(a) - atom.charge == 3 {
                    electrons += 2;
                } else if (z == 8 || z == 16 || z == 34) && m.valence(a) == 2 && atom.charge == 0 {
                    electrons += 2;
                } else if z == 6 && atom.charge == -1 {
                    electrons += 2;
                } else if (z == 6 && atom.charge == 1) || z == 5 {
                } else {
                    ok = false;
                    break;
                }
            }
            if !ok || electrons % 4 != 2 {
                continue;
            }
            for &a in system {
                m.atoms[a as usize].aromatic = true;
            }
            for &a in system {
                for &k in m.bonds_of(a).to_vec().iter() {
                    if inside_set.contains(&m.other(k, a)) {
                        m.bonds[k as usize].aromatic = true;
                    }
                }
            }
            changed = true;
        }
        if !changed {
            break;
        }
    }
    // A bond between aromatic atoms of different rings (biphenyl) is not aromatic.
    for k in 0..m.bonds.len() {
        if m.bonds[k].aromatic && !bond_in_ring(m, k as i32) {
            m.bonds[k].aromatic = false;
        }
    }
}

/// Aromatic bonds as alternating single and double bonds; an aromatic N left
/// without a double bond gets a hydrogen (pyrrole written as c1ccnc1).
pub fn kekulize(m: &mut Molecule) -> Result<(), String> {
    let n = m.atoms.len();
    let (mut required, mut optional) = (vec![0u8; n], vec![0u8; n]);
    for i in 0..n {
        let a = &m.atoms[i];
        if !a.aromatic {
            continue;
        }
        let ii = i as i32;
        if m.bonds_of(ii).iter().any(|&k| !m.bonds[k as usize].aromatic && m.bonds[k as usize].order >= 2) {
            continue;
        }
        let connections = a.hydrogens + m.bonds_of(ii).iter().filter(|&&k| !m.bonds[k as usize].dative).count() as i32;
        match a.z {
            6 => {
                if (a.charge == 0 && connections <= 3) || (a.charge == 1 && connections <= 2) {
                    required[i] = 1;
                }
            }
            5 => {
                if connections <= 2 {
                    optional[i] = 1;
                }
            }
            7 | 15 | 33 => {
                if a.charge == 1 && connections <= 3 {
                    required[i] = 1;
                } else if a.charge == 0 && connections == 2 {
                    optional[i] = 1;
                }
            }
            8 | 16 | 34 => {
                if a.charge == 1 && connections == 2 {
                    required[i] = 1;
                }
            }
            _ => {}
        }
    }
    let mut adjacency = vec![Vec::new(); n];
    let mut required_adjacency = vec![Vec::new(); n];
    for b in &m.bonds {
        if !b.aromatic {
            continue;
        }
        let (a, c) = (b.a as usize, b.b as usize);
        if !((required[a] != 0 || optional[a] != 0) && (required[c] != 0 || optional[c] != 0)) {
            continue;
        }
        adjacency[a].push(b.b);
        adjacency[c].push(b.a);
        if required[a] != 0 && required[c] != 0 {
            required_adjacency[a].push(b.b);
            required_adjacency[c].push(b.a);
        }
    }
    let mut matching = vec![-1i32; n];
    augment_matching(&required_adjacency, &mut matching, &required);
    augment_matching(&adjacency, &mut matching, &vec![1u8; n]);
    for i in 0..n {
        if required[i] != 0 && matching[i] == -1 {
            return Err(format!("cannot kekulize the aromatic system at atom {} (check the ring atoms, e.g. [nH] in pyrrole)", i + 1));
        }
    }
    for b in m.bonds.iter_mut() {
        if b.aromatic {
            b.order = if matching[b.a as usize] == b.b { 2 } else { 1 };
        }
    }
    // An aromatic n written without a hydrogen but needing one (pyrrole as c1ccnc1).
    for i in 0..n {
        if optional[i] != 0 && matching[i] == -1 && m.atoms[i].z == 7 && m.atoms[i].charge == 0 {
            m.atoms[i].hydrogens = m.atoms[i].hydrogens.max(1);
        }
    }
    Ok(())
}

// ------------------------------------------------------------------ hydrogens

/// Implicit hydrogens of atoms written without brackets (organic subset).
pub fn assign_implicit_hydrogens(m: &mut Molecule, bracketed: &[u8]) {
    for i in 0..m.atoms.len() {
        if (i < bracketed.len() && bracketed[i] != 0) || m.atoms[i].dummy {
            continue;
        }
        let mut sum = 0;
        for &k in m.bonds_of(i as i32) {
            if !m.bonds[k as usize].dative {
                sum += m.bonds[k as usize].order;
            }
        }
        let a = &m.atoms[i];
        let v = default_valence(a.z, a.charge, sum + a.hydrogens);
        if v > sum + a.hydrogens {
            m.atoms[i].hydrogens = v - sum;
        }
    }
}

/// 1 sp, 2 sp2, 3 sp3; 0 for metals and hydrogen.
pub fn hybridisation(m: &Molecule, atom: i32) -> i32 {
    let a = &m.atoms[atom as usize];
    if a.z <= 1 || is_metal(a.z) {
        return 0;
    }
    let (mut doubles, mut triples) = (0, 0);
    for &k in m.bonds_of(atom) {
        let b = &m.bonds[k as usize];
        if b.dative {
            continue;
        }
        if b.order == 2 {
            doubles += 1;
        }
        if b.order >= 3 {
            triples += 1;
        }
    }
    if triples > 0 || doubles >= 2 {
        return 1;
    }
    if doubles > 0 || a.aromatic {
        return 2;
    }
    let connections = m.degree(atom) + a.hydrogens;
    if a.z == 6 && a.charge == 1 && connections == 3 {
        return 2;
    }
    if a.z == 5 && connections <= 3 {
        return 2;
    }
    if (a.z == 7 || (a.z == 6 && a.charge == -1)) && connections == 3 {
        // Next to a pi system (amide, aniline, enamine, enolate): planar.
        for &k in m.bonds_of(atom) {
            let b = m.other(k, atom);
            for &e in m.bonds_of(b) {
                let be = &m.bonds[e as usize];
                if e != k && !be.dative && (be.order >= 2 || be.aromatic) {
                    return 2;
                }
            }
        }
    }
    3
}

/// Directions for `count` new bonds on `atom` from its existing bonds and hybridisation.
pub fn new_bond_directions(m: &Molecule, atom: i32, count: i32, hybrid: i32) -> Vec<Vec3> {
    let existing: Vec<Vec3> = m.bonds_of(atom).iter().map(|&k| unit(m.bond_vector(k, atom))).collect();
    let mut out: Vec<Vec3> = Vec::new();
    if count <= 0 {
        return out;
    }
    let total = existing.len() as i32 + count;
    let ideal = match hybrid {
        1 => 2,
        2 => 3,
        _ => 4,
    };
    // A reference direction for staggered and in-plane placements: a neighbour's other neighbour.
    let reference = |u: Vec3| -> Vec3 {
        for &k in m.bonds_of(atom) {
            let nb = m.other(k, atom);
            for &e in m.bonds_of(nb) {
                if m.other(e, nb) == atom {
                    continue;
                }
                let d = m.bond_vector(e, nb);
                let p = sub(d, scale(u, dot(d, u)));
                if norm(p) > 1e-3 {
                    return unit(p);
                }
            }
        }
        perpendicular(u)
    };
    let s3 = 3f64.sqrt();
    if total <= ideal {
        if existing.is_empty() {
            out = match ideal {
                2 => vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
                3 => vec![[1.0, 0.0, 0.0], [-0.5, s3 / 2.0, 0.0], [-0.5, -s3 / 2.0, 0.0]],
                _ => {
                    let s = 1.0 / s3;
                    vec![[s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]]
                }
            };
        } else if existing.len() == 1 {
            let u = existing[0];
            if ideal == 2 {
                out = vec![scale(u, -1.0)];
            } else if ideal == 3 {
                let p = reference(u);
                // In the plane of the reference, the first trans to it.
                out = vec![add(scale(u, -0.5), scale(p, -s3 / 2.0)), add(scale(u, -0.5), scale(p, s3 / 2.0))];
            } else {
                let p = reference(u);
                let q = cross(u, p);
                let s = (8.0f64 / 9.0).sqrt();
                for phi in [180.0f64, 60.0, -60.0] {
                    let r = phi * PI / 180.0;
                    out.push(add(scale(u, -1.0 / 3.0), add(scale(p, s * r.cos()), scale(q, s * r.sin()))));
                }
            }
        } else if existing.len() == 2 {
            let bis = scale(unit(add(existing[0], existing[1])), -1.0);
            if ideal == 3 {
                out = vec![bis];
            } else if ideal == 4 {
                let mut nrm = cross(existing[0], existing[1]);
                nrm = if norm(nrm) < 1e-6 { perpendicular(existing[0]) } else { unit(nrm) };
                let (ca, sa) = (1.0 / s3, (2.0f64 / 3.0).sqrt());
                out = vec![add(scale(bis, ca), scale(nrm, sa)), add(scale(bis, ca), scale(nrm, -sa))];
            }
        } else if existing.len() == 3 {
            let mut sum = add(add(existing[0], existing[1]), existing[2]);
            if norm(sum) < 0.1 {
                sum = cross(sub(existing[1], existing[0]), sub(existing[2], existing[0]));
            }
            out = vec![scale(unit(sum), -1.0)];
        }
    }
    if (out.len() as i32) < count {
        // Fallback: directions farthest from the existing bonds (Fibonacci sphere).
        out.clear();
        let mut taken = existing.clone();
        for _ in 0..count {
            let mut best = [0.0, 0.0, 1.0];
            let mut best_score = -2.0;
            let samples = 600;
            for s in 0..samples {
                let y = 1.0 - 2.0 * (s as f64 + 0.5) / samples as f64;
                let r = (1.0 - y * y).sqrt();
                let phi = s as f64 * 2.39996322972865332;
                let d = [r * phi.cos(), y, r * phi.sin()];
                let mut closest = 2.0;
                for t in &taken {
                    closest = fmin(closest, 1.0 - dot(d, *t));
                }
                if taken.is_empty() {
                    closest = 1.0;
                }
                if closest > best_score {
                    best_score = closest;
                    best = d;
                }
            }
            out.push(best);
            taken.push(best);
        }
    }
    out.resize(count as usize, [0.0; 3]);
    out.into_iter().map(unit).collect()
}

/// Makes implicit hydrogens explicit atoms (placed by hybridisation when
/// `place`); stereo neighbour lists refer to the new atoms. Returns the number added.
pub fn make_hydrogens_explicit(m: &mut Molecule, place: bool) -> i32 {
    let mut added = 0;
    let n = m.atoms.len();
    for i in 0..n {
        let count = m.atoms[i].hydrogens;
        if count <= 0 {
            continue;
        }
        let ii = i as i32;
        let hybrid = hybridisation(m, ii);
        m.atoms[i].hydrogens = 0;
        let directions = if place { new_bond_directions(m, ii, count, if hybrid == 0 { 3 } else { hybrid }) } else { Vec::new() };
        let length = hydrogen_bond_length(m.atoms[i].z);
        for h in 0..count as usize {
            let position = if place { add(m.atoms[i].position, scale(directions[h], length)) } else { m.atoms[i].position };
            let index = m.add_atom(Atom { z: 1, position, ..Default::default() });
            m.add_bond(ii, index, 1).expect("distinct atoms");
            added += 1;
            for t in m.tetrahedral.iter_mut() {
                if t.centre == ii {
                    if let Some(x) = t.n.iter_mut().find(|x| **x == -1) {
                        *x = index;
                    }
                }
            }
        }
    }
    added
}

// ------------------------------------------------------------------ properties

/// Hill-order formula (implicit hydrogens included).
pub fn formula(m: &Molecule) -> String {
    let mut counts: BTreeMap<String, i32> = BTreeMap::new();
    let mut carbon = false;
    for a in &m.atoms {
        if a.dummy || a.z <= 0 {
            continue;
        }
        *counts.entry(symbol_of(a.z)).or_insert(0) += 1;
        if a.hydrogens != 0 {
            *counts.entry("H".into()).or_insert(0) += a.hydrogens;
        }
        carbon = carbon || a.z == 6;
    }
    let mut out = String::new();
    let mut put = |counts: &mut BTreeMap<String, i32>, symbol: &str| {
        if let Some(count) = counts.remove(symbol) {
            out += symbol;
            if count > 1 {
                out += &count.to_string();
            }
        }
    };
    if carbon {
        put(&mut counts, "C");
        put(&mut counts, "H");
    }
    while let Some(first) = counts.keys().next().cloned() {
        put(&mut counts, &first);
    }
    out
}

pub fn molar_mass(m: &Molecule) -> f64 {
    let mut sum = 0.0;
    for a in &m.atoms {
        if a.dummy || a.z <= 0 {
            continue;
        }
        sum += elements::atomic_mass(a.z) + f64::from(a.hydrogens) * elements::atomic_mass(1);
    }
    sum
}

pub fn monoisotopic_mass(m: &Molecule) -> f64 {
    const TABLE: [(i32, f64); 15] = [
        (1, 1.00782503207),
        (5, 11.0093054),
        (6, 12.0),
        (7, 14.0030740048),
        (8, 15.99491461956),
        (9, 18.99840322),
        (11, 22.9897692809),
        (14, 27.9769265325),
        (15, 30.97376163),
        (16, 31.972071),
        (17, 34.96885268),
        (19, 38.96370668),
        (34, 79.9165213),
        (35, 78.9183371),
        (53, 126.904473),
    ];
    let mass = |z: i32| TABLE.iter().find(|e| e.0 == z).map(|e| e.1);
    let mut sum = 0.0;
    for a in &m.atoms {
        if a.dummy || a.z <= 0 {
            continue;
        }
        sum += mass(a.z).unwrap_or_else(|| elements::atomic_mass(a.z)) + f64::from(a.hydrogens) * TABLE[0].1;
    }
    sum
}

pub fn total_charge(m: &Molecule) -> i32 {
    m.atoms.iter().map(|a| a.charge).sum()
}

/// Sum of Z (with implicit hydrogens) minus the charge.
pub fn electron_count(m: &Molecule) -> i32 {
    m.atoms.iter().map(|a| a.z + a.hydrogens).sum::<i32>() - total_charge(m)
}

/// Spin multiplicity 2S+1 from radicals (valence deficits) or the electron count's parity.
pub fn multiplicity(m: &Molecule) -> i32 {
    // O2 is the familiar exception: a triplet ground state.
    if m.atoms.len() == 2 && m.atoms[0].z == 8 && m.atoms[1].z == 8 && total_charge(m) == 0 {
        return 3;
    }
    let mut unpaired = 0;
    let mut metal = false;
    for (i, a) in m.atoms.iter().enumerate() {
        if is_metal(a.z) {
            metal = true;
            continue;
        }
        if a.z <= 0 {
            continue;
        }
        let sum = m.valence(i as i32);
        let v = default_valence(a.z, a.charge, sum);
        if v > sum {
            unpaired += v - sum;
        }
    }
    let electrons = electron_count(m);
    if metal || unpaired == 0 {
        return if electrons % 2 != 0 { 2 } else { 1 };
    }
    if unpaired % 2 != electrons % 2 {
        unpaired = electrons % 2;
    }
    unpaired + 1
}

/// Connected components (sorted atom lists).
pub fn fragments(m: &Molecule) -> Vec<Vec<i32>> {
    let n = m.atoms.len();
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for i in 0..n {
        if seen[i] {
            continue;
        }
        let mut part = Vec::new();
        let mut stack = vec![i as i32];
        seen[i] = true;
        while let Some(a) = stack.pop() {
            part.push(a);
            for b in m.neighbours(a) {
                if !seen[b as usize] {
                    seen[b as usize] = true;
                    stack.push(b);
                }
            }
        }
        part.sort();
        out.push(part);
    }
    out
}

/// Single, non-ring bonds between two heavy atoms each with another heavy
/// neighbour, excluding amide C-N and triple-bond neighbours.
pub fn rotatable_bonds(m: &Molecule) -> Vec<i32> {
    let mut out = Vec::new();
    for (k, b) in m.bonds.iter().enumerate() {
        if b.order != 1 || b.aromatic || b.dative {
            continue;
        }
        if m.atoms[b.a as usize].z == 1 || m.atoms[b.b as usize].z == 1 {
            continue;
        }
        if m.heavy_degree(b.a) < 2 || m.heavy_degree(b.b) < 2 {
            continue;
        }
        if bond_in_ring(m, k as i32) {
            continue;
        }
        if [b.a, b.b].iter().any(|&x| m.bonds_of(x).iter().any(|&e| m.bonds[e as usize].order == 3)) {
            continue;
        }
        // Amide C-N: C bonded to O= (or S=).
        let amide = |c: i32, n_atom: i32| {
            m.atoms[c as usize].z == 6
                && m.atoms[n_atom as usize].z == 7
                && m.bonds_of(c).iter().any(|&e| m.bonds[e as usize].order == 2 && [8, 16].contains(&m.atoms[m.other(e, c) as usize].z))
        };
        if amide(b.a, b.b) || amide(b.b, b.a) {
            continue;
        }
        out.push(k as i32);
    }
    out
}

// ------------------------------------------------------------------ ranks and CIP

/// Graph-invariant ranks (equal for symmetry-equivalent atoms; Morgan-style refinement).
pub fn canonical_ranks(m: &Molecule) -> Vec<i32> {
    let n = m.atoms.len();
    let rings = smallest_rings(m);
    let ring_size = smallest_ring_size(m, &rings);
    let mut invariant = vec![0i64; n];
    for i in 0..n {
        let a = &m.atoms[i];
        let ii = i as i32;
        invariant[i] = ((((i64::from(a.z) * 16 + i64::from(m.degree(ii))) * 16 + i64::from(m.hydrogen_count(ii))) * 16 + i64::from(a.charge + 8)) * 64 + i64::from(ring_size[i])) * 2
            + i64::from(a.aromatic);
        invariant[i] = invariant[i] * 512 + i64::from(a.isotope);
    }
    let to_ranks = |keys: &[Vec<i64>]| -> Vec<i32> {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&x, &y| keys[x].cmp(&keys[y]));
        let mut rank = vec![0i32; n];
        for k in 1..n {
            rank[order[k]] = rank[order[k - 1]] + i32::from(keys[order[k]] != keys[order[k - 1]]);
        }
        rank
    };
    let mut keys: Vec<Vec<i64>> = invariant.iter().map(|&v| vec![v]).collect();
    let mut rank = to_ranks(&keys);
    for _ in 0..n + 1 {
        for i in 0..n {
            let mut nb: Vec<i64> = m.bonds_of(i as i32).iter().map(|&k| i64::from(rank[m.other(k, i as i32) as usize]) * 8 + i64::from(m.bonds[k as usize].order)).collect();
            nb.sort();
            keys[i] = std::iter::once(i64::from(rank[i])).chain(nb).collect();
        }
        let next = to_ranks(&keys);
        let before = rank.iter().copied().max().unwrap_or(0);
        let after = next.iter().copied().max().unwrap_or(0);
        rank = next;
        if after == before {
            break;
        }
    }
    rank
}

/// Signed volume (n1-c).((n2-c)x(n3-c)) of a centre's neighbours.
pub fn signed_volume(c: Vec3, n1: Vec3, n2: Vec3, n3: Vec3) -> f64 {
    dot(sub(n1, c), cross(sub(n2, c), sub(n3, c)))
}

/// A node of the CIP hierarchical digraph.
#[derive(Clone, Debug)]
struct CipNode {
    z: i32,
    /// -1 for duplicated and phantom atoms.
    atom: i32,
    /// Atom it was reached from.
    from: i32,
    /// Atoms on the way from the root (to stop at ring closures).
    path: Vec<i32>,
}

/// Compares two branches (root atoms a and b seen from centre) sphere by
/// sphere; > 0 when a ranks higher.
fn compare_branches(m: &Molecule, centre: i32, a: i32, b: i32, depth: i32) -> i32 {
    let children = |node: &CipNode| -> Vec<CipNode> {
        let mut out = Vec::new();
        if node.atom < 0 {
            return out; // duplicated atoms carry phantom substituents only
        }
        for &k in m.bonds_of(node.atom) {
            let x = m.other(k, node.atom);
            let bond = &m.bonds[k as usize];
            let order = if bond.dative { 1 } else { bond.order };
            let zx = m.atoms[x as usize].z;
            if x == node.from {
                // The way back: duplicates for the extra bond order only.
                for _ in 1..order {
                    out.push(CipNode { z: zx, atom: -1, from: node.atom, path: Vec::new() });
                }
                continue;
            }
            let closure = node.path.contains(&x);
            let mut path = node.path.clone();
            path.push(node.atom);
            out.push(CipNode { z: zx, atom: if closure { -1 } else { x }, from: node.atom, path });
            for _ in 1..order {
                out.push(CipNode { z: zx, atom: -1, from: node.atom, path: Vec::new() });
            }
        }
        for _ in 0..m.atoms[node.atom as usize].hydrogens {
            out.push(CipNode { z: 1, atom: -1, from: node.atom, path: Vec::new() });
        }
        sort_by(&mut out, |x, y| x.z > y.z);
        out
    };
    let mut la = vec![CipNode { z: m.atoms[a as usize].z, atom: a, from: centre, path: vec![centre] }];
    let mut lb = vec![CipNode { z: m.atoms[b as usize].z, atom: b, from: centre, path: vec![centre] }];
    if la[0].z != lb[0].z {
        return if la[0].z > lb[0].z { 1 } else { -1 };
    }
    for _ in 0..depth {
        let (mut na, mut nb) = (Vec::new(), Vec::new());
        // Compare the substituent sets of each node in turn (hierarchically).
        for i in 0..la.len().min(lb.len()) {
            let (ca, cb) = (children(&la[i]), children(&lb[i]));
            let mut za: Vec<i32> = ca.iter().map(|c| c.z).collect();
            let mut zb: Vec<i32> = cb.iter().map(|c| c.z).collect();
            let len = za.len().max(zb.len());
            za.resize(len, 0);
            zb.resize(len, 0);
            if za != zb {
                return if za > zb { 1 } else { -1 };
            }
            na.extend(ca);
            nb.extend(cb);
        }
        if la.len() != lb.len() {
            return if la.len() > lb.len() { 1 } else { -1 };
        }
        if na.is_empty() && nb.is_empty() {
            return 0;
        }
        la = na;
        lb = nb;
    }
    0
}

/// Branches of an atom ordered by CIP priority (highest first); None on a tie.
fn cip_order(m: &Molecule, centre: i32, mut branches: Vec<i32>) -> Option<Vec<i32>> {
    for x in 0..branches.len() {
        for y in x + 1..branches.len() {
            if compare_branches(m, centre, branches[x], branches[y], 12) == 0 {
                return None;
            }
        }
    }
    sort_by(&mut branches, |&x, &y| compare_branches(m, centre, x, y, 12) > 0);
    Some(branches)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoLabel {
    /// Stereocentre, or first atom of the double bond.
    pub atom: i32,
    /// Second atom of the double bond, else -1.
    pub other: i32,
    /// "R", "S", "E", "Z".
    pub label: String,
}

/// CIP labels of tetrahedral centres with four different branches and of
/// double bonds with two different branches at each end, from coordinates.
pub fn cip_labels(m: &Molecule) -> Vec<StereoLabel> {
    let mut out = Vec::new();
    let n = m.atoms.len();
    // Neighbour position as bonded (periodic images).
    let pos = |centre: i32, x: i32| add(m.atoms[centre as usize].position, m.bond_vector(m.bond_between(centre, x), centre));
    for i in 0..n {
        let a = &m.atoms[i];
        let ii = i as i32;
        if a.z <= 1 || is_metal(a.z) || a.hydrogens > 0 {
            continue;
        }
        let nb = m.neighbours(ii);
        if nb.len() != 4 || hybridisation(m, ii) != 3 {
            continue;
        }
        let Some(sorted) = cip_order(m, ii, nb) else { continue };
        let v = signed_volume(a.position, pos(ii, sorted[0]), pos(ii, sorted[1]), pos(ii, sorted[2]));
        if v.abs() < 1e-3 {
            continue;
        }
        out.push(StereoLabel { atom: ii, other: -1, label: if v < 0.0 { "R" } else { "S" }.into() });
    }
    for (k, b) in m.bonds.iter().enumerate() {
        if b.order != 2 || b.aromatic || bond_in_ring(m, k as i32) {
            // Ring double bonds of eight or more atoms may be E; smaller rings are always Z.
            if !(b.order == 2 && !b.aromatic) {
                continue;
            }
            let rings = smallest_rings(m);
            if rings.iter().any(|r| r.len() < 8 && r.contains(&b.a) && r.contains(&b.b)) {
                continue;
            }
        }
        let ea: Vec<i32> = m.neighbours(b.a).into_iter().filter(|&x| x != b.b).collect();
        let eb: Vec<i32> = m.neighbours(b.b).into_iter().filter(|&x| x != b.a).collect();
        if ea.is_empty() || eb.is_empty() || ea.len() > 2 || eb.len() > 2 {
            continue;
        }
        if m.atoms[b.a as usize].hydrogens + ea.len() as i32 != 2 || m.atoms[b.b as usize].hydrogens + eb.len() as i32 != 2 {
            continue;
        }
        let sa = if ea.len() == 2 { cip_order(m, b.a, ea) } else { Some(ea) };
        let Some(sa) = sa else { continue };
        let sb = if eb.len() == 2 { cip_order(m, b.b, eb) } else { Some(eb) };
        let Some(sb) = sb else { continue };
        // Same side when the dihedral sa-a-b-sb is below 90 degrees.
        let pa = m.atoms[b.a as usize].position;
        let pb = pos(b.a, b.b);
        let x1 = pos(b.a, sa[0]);
        let x2 = add(pb, sub(pos(b.b, sb[0]), m.atoms[b.b as usize].position));
        let axis = unit(sub(pb, pa));
        let mut u = sub(x1, pa);
        let mut w = sub(x2, pb);
        u = sub(u, scale(axis, dot(u, axis)));
        w = sub(w, scale(axis, dot(w, axis)));
        if norm(u) < 1e-6 || norm(w) < 1e-6 {
            continue;
        }
        out.push(StereoLabel { atom: b.a, other: b.b, label: if dot(u, w) > 0.0 { "Z" } else { "E" }.into() });
    }
    out
}
