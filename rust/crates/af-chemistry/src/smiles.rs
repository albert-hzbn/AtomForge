//! SMILES reading and writing (the C++ chemistry/Smiles): organic subset and
//! bracket atoms, branches, ring closures, aromaticity, tetrahedral '@'/'@@'
//! and double-bond '/'/'\' stereo; writing from canonical ranks with stereo
//! from coordinates.

use std::collections::{BTreeMap, BTreeSet};

use crate::molecule::*;

/// A SMILES error, at a character position when the text caused it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmilesError {
    pub message: String,
    /// None for errors of the molecule rather than of the text (kekulization).
    pub position: Option<usize>,
}

impl std::fmt::Display for SmilesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

fn error(message: impl Into<String>, at: usize) -> SmilesError {
    SmilesError { message: message.into(), position: Some(at) }
}

fn organic_subset(z: i32) -> bool {
    matches!(z, 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53)
}

fn aromatic_symbol(z: i32) -> bool {
    matches!(z, 5 | 6 | 7 | 8 | 15 | 16 | 33 | 34)
}

#[derive(Clone, Copy, Debug, Default)]
struct PendingBond {
    /// 0: default (single, or aromatic between aromatic atoms).
    order: i32,
    /// ':'.
    aromatic: bool,
    /// b'/' or b'\\', or 0.
    direction: u8,
}

#[derive(Clone, Copy, Debug)]
struct RingOpen {
    atom: i32,
    bond: PendingBond,
    /// Position in the atom's stereo neighbour order.
    slot: usize,
}

struct DirectedBond {
    bond: i32,
    first: i32,
    direction: u8,
}

struct Parser {
    m: Molecule,
    bracketed: Vec<u8>,
    /// 0, 1 '@', 2 '@@'.
    chirality: Vec<i32>,
    /// Stereo neighbour order (-1 implicit H, < -1 an open ring slot).
    order: Vec<Vec<i32>>,
    directed: Vec<DirectedBond>,
    previous: i32,
    pending: PendingBond,
    dot: bool,
}

impl Parser {
    fn connect(&mut self, a: i32, b: i32, bond: PendingBond, at: usize) -> Result<i32, SmilesError> {
        if self.m.bond_between(a, b) >= 0 {
            return Err(error(format!("atoms {} and {} are bonded twice", a + 1, b + 1), at));
        }
        let mut order = bond.order;
        let mut aromatic = bond.aromatic;
        if order == 0 {
            if self.m.atoms[a as usize].aromatic && self.m.atoms[b as usize].aromatic && bond.direction == 0 {
                aromatic = true;
            }
            order = 1;
        }
        let k = self.m.add_bond(a, b, order).map_err(|e| error(e, at))?;
        self.m.bonds[k as usize].aromatic = aromatic && self.m.atoms[a as usize].aromatic && self.m.atoms[b as usize].aromatic;
        if bond.direction != 0 {
            self.directed.push(DirectedBond { bond: k, first: a, direction: bond.direction });
        }
        Ok(k)
    }

    fn add_atom(&mut self, mut atom: Atom, bracket: bool, chiral: i32, hydrogens: i32, at: usize) -> Result<(), SmilesError> {
        atom.hydrogens = if bracket { hydrogens } else { 0 };
        let index = self.m.add_atom(atom);
        self.bracketed.push(u8::from(bracket));
        self.chirality.push(chiral);
        self.order.push(Vec::new());
        if self.previous >= 0 && !self.dot {
            let previous = self.previous;
            self.connect(previous, index, self.pending, at)?;
            self.order[previous as usize].push(index);
            self.order[index as usize].push(previous);
        }
        if bracket && hydrogens > 0 && chiral != 0 {
            self.order[index as usize].push(-1);
        }
        self.pending = PendingBond::default();
        self.dot = false;
        self.previous = index;
        Ok(())
    }
}

fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\x0B' | b'\x0C' | b'\r')
}

/// Reads a SMILES string (it ends at whitespace; the rest is a title).
pub fn parse_smiles(text: &str) -> Result<Molecule, SmilesError> {
    let t = text.as_bytes();
    let n = t.len();
    let mut p = Parser {
        m: Molecule::default(),
        bracketed: Vec::new(),
        chirality: Vec::new(),
        order: Vec::new(),
        directed: Vec::new(),
        previous: -1,
        pending: PendingBond::default(),
        dot: false,
    };
    let mut rings: BTreeMap<i32, RingOpen> = BTreeMap::new();
    let mut branches: Vec<i32> = Vec::new();
    let mut i = 0usize;
    while i < n {
        let c = t[i];
        if is_space(c) {
            break;
        }
        match c {
            b'(' => {
                if p.previous < 0 {
                    return Err(error("a branch '(' needs an atom before it", i));
                }
                branches.push(p.previous);
                i += 1;
                continue;
            }
            b')' => {
                let Some(b) = branches.pop() else {
                    return Err(error("')' without a matching '('", i));
                };
                p.previous = b;
                i += 1;
                continue;
            }
            b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\' => {
                match c {
                    b'-' => p.pending.order = 1,
                    b'=' => p.pending.order = 2,
                    b'#' => p.pending.order = 3,
                    b'$' => p.pending.order = 4,
                    b':' => p.pending.aromatic = true,
                    _ => p.pending.direction = c,
                }
                if c == b'/' || c == b'\\' {
                    p.pending.order = if p.pending.order != 0 { p.pending.order } else { 1 };
                }
                i += 1;
                continue;
            }
            b'.' => {
                p.dot = true;
                p.pending = PendingBond::default();
                i += 1;
                continue;
            }
            _ => {}
        }
        if c.is_ascii_digit() || c == b'%' {
            if p.previous < 0 {
                return Err(error("a ring-closure digit needs an atom before it", i));
            }
            let at = i;
            let digit;
            if c == b'%' {
                if i + 2 >= n || !t[i + 1].is_ascii_digit() || !t[i + 2].is_ascii_digit() {
                    return Err(error("'%' must be followed by two digits", i));
                }
                digit = i32::from(t[i + 1] - b'0') * 10 + i32::from(t[i + 2] - b'0');
                i += 3;
            } else {
                digit = i32::from(c - b'0');
                i += 1;
            }
            match rings.remove(&digit) {
                None => {
                    let previous = p.previous;
                    let slot = p.order[previous as usize].len();
                    p.order[previous as usize].push(-100 - digit);
                    rings.insert(digit, RingOpen { atom: previous, bond: p.pending, slot });
                }
                Some(open) => {
                    if open.atom == p.previous {
                        return Err(error(format!("ring bond {digit} closes on its own atom"), at));
                    }
                    let bond = if p.pending.order != 0 || p.pending.aromatic || p.pending.direction != 0 { p.pending } else { open.bond };
                    if open.bond.order != 0 && p.pending.order != 0 && open.bond.order != p.pending.order {
                        return Err(error(format!("ring bond {digit} has two different bond orders"), at));
                    }
                    let previous = p.previous;
                    p.connect(open.atom, previous, bond, at)?;
                    p.order[open.atom as usize][open.slot] = previous;
                    p.order[previous as usize].push(open.atom);
                }
            }
            p.pending = PendingBond::default();
            continue;
        }
        // Atoms.
        let mut atom = Atom::default();
        if c == b'[' {
            let Some(close) = text[i..].find(']').map(|k| k + i) else {
                return Err(error("'[' is not closed with ']'", i));
            };
            let mut j = i + 1;
            let mut isotope = 0;
            while j < close && t[j].is_ascii_digit() {
                isotope = isotope * 10 + i32::from(t[j] - b'0');
                j += 1;
            }
            if j >= close {
                return Err(error("missing element symbol in brackets", i));
            }
            let symbol: String;
            let mut aromatic = false;
            if t[j] == b'*' {
                symbol = "*".into();
                j += 1;
            } else if t[j].is_ascii_lowercase() {
                // Aromatic: c, n, o, p, s, b, se, as.
                if j + 1 < close && (&t[j..j + 2] == b"se" || &t[j..j + 2] == b"as") {
                    symbol = text[j..j + 2].to_string();
                    j += 2;
                } else {
                    symbol = (t[j] as char).to_string();
                    j += 1;
                }
                aromatic = true;
            } else {
                let mut s = (t[j] as char).to_string();
                j += 1;
                if j < close && t[j].is_ascii_lowercase() {
                    let two = format!("{s}{}", t[j] as char);
                    if atomic_number(&two) > 0 {
                        s = two;
                        j += 1;
                    }
                }
                symbol = s;
            }
            let z = if symbol == "*" { 0 } else { atomic_number(&symbol) };
            if symbol != "*" && z <= 0 {
                return Err(error(format!("unknown element '{symbol}'"), i + 1));
            }
            if aromatic && !aromatic_symbol(z) {
                return Err(error(format!("'{symbol}' cannot be aromatic"), i + 1));
            }
            let mut chiral = 0;
            if j < close && t[j] == b'@' {
                j += 1;
                chiral = 1;
                if j < close && t[j] == b'@' {
                    chiral = 2;
                    j += 1;
                } else if j + 2 < close && &t[j..j + 2] == b"TH" {
                    chiral = if t[j + 2] == b'2' { 2 } else { 1 };
                    j += 3;
                } else if j + 1 < close && t[j].is_ascii_uppercase() {
                    // @SP, @TB, @OH: not supported; skip the descriptor.
                    chiral = 0;
                    while j < close && t[j].is_ascii_alphanumeric() {
                        j += 1;
                    }
                }
            }
            let mut hydrogens = 0;
            if j < close && t[j] == b'H' {
                j += 1;
                hydrogens = 1;
                if j < close && t[j].is_ascii_digit() {
                    hydrogens = i32::from(t[j] - b'0');
                    j += 1;
                }
            }
            let mut charge = 0;
            while j < close && (t[j] == b'+' || t[j] == b'-') {
                let sign = if t[j] == b'+' { 1 } else { -1 };
                j += 1;
                if j < close && t[j].is_ascii_digit() {
                    let mut value = 0;
                    while j < close && t[j].is_ascii_digit() {
                        value = value * 10 + i32::from(t[j] - b'0');
                        j += 1;
                    }
                    charge += sign * value;
                } else {
                    charge += sign;
                }
            }
            let mut map_number = 0;
            if j < close && t[j] == b':' {
                j += 1;
                while j < close && t[j].is_ascii_digit() {
                    map_number = map_number * 10 + i32::from(t[j] - b'0');
                    j += 1;
                }
            }
            if j != close {
                return Err(error(format!("unexpected '{}' in brackets", t[j] as char), j));
            }
            atom.z = z;
            atom.dummy = z == 0;
            atom.aromatic = aromatic;
            atom.isotope = isotope;
            atom.charge = charge;
            atom.map_number = map_number;
            p.add_atom(atom, true, chiral, hydrogens, i)?;
            i = close + 1;
            continue;
        }
        if c == b'*' {
            atom.z = 0;
            atom.dummy = true;
            p.add_atom(atom, true, 0, 0, i)?;
            i += 1;
            continue;
        }
        let symbol = if c == b'C' && i + 1 < n && t[i + 1] == b'l' {
            "Cl".to_string()
        } else if c == b'B' && i + 1 < n && t[i + 1] == b'r' {
            "Br".to_string()
        } else {
            (c as char).to_string()
        };
        let aromatic = c.is_ascii_lowercase();
        let z;
        if aromatic {
            z = match c {
                b'b' => 5,
                b'c' => 6,
                b'n' => 7,
                b'o' => 8,
                b'p' => 15,
                b's' => 16,
                _ => return Err(error(format!("'{symbol}' is not an atom (use brackets, e.g. [Na])"), i)),
            };
        } else {
            z = atomic_number(&symbol);
            if !organic_subset(z) {
                return Err(error(format!("'{symbol}' must be written in brackets, e.g. [{symbol}]"), i));
            }
        }
        atom.z = z;
        atom.aromatic = aromatic;
        p.add_atom(atom, false, 0, 0, i)?;
        i += symbol.len();
    }
    if let Some((digit, _)) = rings.iter().next() {
        return Err(error(format!("ring bond {digit} is not closed"), n));
    }
    if !branches.is_empty() {
        return Err(error("'(' is not closed with ')'", n));
    }
    if p.m.atoms.is_empty() {
        return Err(error("the SMILES has no atoms", 0));
    }
    let Parser { mut m, bracketed, chirality, order, directed, .. } = p;
    // Stereocentres with four neighbours (an implicit H counts as one).
    for a in 0..m.atoms.len() {
        if chirality[a] == 0 || order[a].len() != 4 {
            continue; // lone-pair centres (e.g. sulfoxides) are not handled
        }
        let o = &order[a];
        m.tetrahedral.push(TetrahedralStereo { centre: a as i32, n: [o[0], o[1], o[2], o[3]], sign: if chirality[a] == 2 { 1 } else { -1 } });
    }
    kekulize(&mut m).map_err(|message| SmilesError { message, position: None })?;
    assign_implicit_hydrogens(&mut m, &bracketed);
    // Double-bond stereo from bond directions: equal "up" on both ends means cis.
    let mut stereo = Vec::new();
    for b in &m.bonds {
        if b.order != 2 || b.aromatic {
            continue;
        }
        let side = |x: i32, y: i32| -> Option<(i32, bool)> {
            for d in &directed {
                let e = &m.bonds[d.bond as usize];
                if e.a != x && e.b != x {
                    continue;
                }
                let other = if e.a == x { e.b } else { e.a };
                if other == y {
                    continue;
                }
                return Some((other, (d.direction == b'/') == (d.first == x)));
            }
            None
        };
        if let (Some((a, up_a)), Some((c2, up_b))) = (side(b.a, b.b), side(b.b, b.a)) {
            stereo.push(DoubleBondStereo { a, x: b.a, y: b.b, b: c2, cis: up_a == up_b });
        }
    }
    m.double_bonds.extend(stereo);
    Ok(m)
}

/// Writes canonical-ranked SMILES with stereo (from coordinates) when `stereo`.
pub fn write_smiles(m: &Molecule, stereo: bool) -> String {
    let n = m.atoms.len();
    // Hydrogens folded into their heavy atom.
    let mut kept = vec![true; n];
    let mut hcount = vec![0i32; n];
    let mut folded_h = vec![-1i32; n]; // one folded H per atom, for stereo
    for i in 0..n {
        let a = &m.atoms[i];
        hcount[i] += a.hydrogens;
        let ii = i as i32;
        if a.z != 1 || a.charge != 0 || a.isotope != 0 || m.degree(ii) != 1 {
            continue;
        }
        let heavy = m.neighbours(ii)[0];
        if m.atoms[heavy as usize].z == 1 {
            continue;
        }
        kept[i] = false;
        hcount[heavy as usize] += 1;
        folded_h[heavy as usize] = ii;
    }
    let rank = canonical_ranks(m);
    let position = |centre: i32, x: i32| add(m.atoms[centre as usize].position, m.bond_vector(m.bond_between(centre, x), centre));
    let kept_neighbours = |a: i32| -> Vec<i32> {
        let mut out: Vec<i32> = m.neighbours(a).into_iter().filter(|&x| kept[x as usize]).collect();
        out.sort_by(|&x, &y| (rank[x as usize], x).cmp(&(rank[y as usize], y)));
        out
    };

    // Pass 1: depth-first tree, ring closures and the order atoms are written.
    let mut parent = vec![-1i32; n];
    let mut preorder = vec![-1i32; n];
    let mut children: Vec<Vec<i32>> = vec![Vec::new(); n];
    let mut closures: Vec<Vec<i32>> = vec![Vec::new(); n];
    let mut roots = Vec::new();
    let mut counter = 0;
    let mut state = vec![0u8; n]; // 0 new, 1 open, 2 done
    let mut closure_bonds: BTreeSet<i32> = BTreeSet::new();
    struct Dfs<'a> {
        m: &'a Molecule,
        parent: &'a mut Vec<i32>,
        preorder: &'a mut Vec<i32>,
        children: &'a mut Vec<Vec<i32>>,
        closures: &'a mut Vec<Vec<i32>>,
        state: &'a mut Vec<u8>,
        closure_bonds: &'a mut BTreeSet<i32>,
        counter: &'a mut i32,
    }
    fn visit(d: &mut Dfs, a: i32, kept_neighbours: &dyn Fn(i32) -> Vec<i32>) {
        d.state[a as usize] = 1;
        d.preorder[a as usize] = *d.counter;
        *d.counter += 1;
        for x in kept_neighbours(a) {
            if x == d.parent[a as usize] {
                continue;
            }
            let k = d.m.bond_between(a, x);
            if d.state[x as usize] == 1 {
                if !d.closure_bonds.insert(k) {
                    continue;
                }
                d.closures[x as usize].push(k); // opened at the ancestor x
                d.closures[a as usize].push(k); // closed here
            } else if d.state[x as usize] == 0 {
                d.parent[x as usize] = a;
                d.children[a as usize].push(x);
                visit(d, x, kept_neighbours);
            }
        }
        d.state[a as usize] = 2;
    }
    {
        let mut order: Vec<i32> = (0..n as i32).collect();
        order.sort_by(|&x, &y| {
            // Start fragments at a low-rank end atom.
            let (dx, dy) = (kept_neighbours(x).len(), kept_neighbours(y).len());
            if (dx <= 1) != (dy <= 1) {
                return if dx <= 1 { std::cmp::Ordering::Less } else { std::cmp::Ordering::Greater };
            }
            (rank[x as usize], x).cmp(&(rank[y as usize], y))
        });
        let mut d = Dfs {
            m,
            parent: &mut parent,
            preorder: &mut preorder,
            children: &mut children,
            closures: &mut closures,
            state: &mut state,
            closure_bonds: &mut closure_bonds,
            counter: &mut counter,
        };
        for a in order {
            if kept[a as usize] && d.state[a as usize] == 0 {
                roots.push(a);
                visit(&mut d, a, &kept_neighbours);
            }
        }
    }
    // Closures at an atom: closing ones first, then openings, by the partner's preorder.
    for a in 0..n {
        let pa = preorder[a];
        closures[a].sort_by_key(|&k| {
            let partner = preorder[m.other(k, a as i32) as usize];
            (u8::from(partner >= pa), partner)
        });
    }

    // Double-bond stereo marks: bond index -> direction character.
    let mut directions: BTreeMap<i32, u8> = BTreeMap::new();
    if stereo {
        for label in cip_labels(m) {
            if label.other < 0 {
                continue;
            }
            let (mut x, mut y) = (label.atom, label.other);
            if preorder[x as usize] > preorder[y as usize] {
                std::mem::swap(&mut x, &mut y);
            }
            if preorder[x as usize] < 0 || preorder[y as usize] < 0 {
                continue;
            }
            // A written, non-closure neighbour bond.
            let pick = |end: i32, partner: i32| kept_neighbours(end).into_iter().find(|&z| z != partner && !closure_bonds.contains(&m.bond_between(end, z))).unwrap_or(-1);
            let (a, b) = (pick(x, y), pick(y, x));
            if a < 0 || b < 0 {
                continue;
            }
            // Geometric cis/trans of a and b.
            let px = m.atoms[x as usize].position;
            let py = position(x, y);
            let pa = position(x, a);
            let py_atom = m.atoms[y as usize].position;
            let pb_raw = position(y, b);
            let pb = [py[0] + pb_raw[0] - py_atom[0], py[1] + pb_raw[1] - py_atom[1], py[2] + pb_raw[2] - py_atom[2]];
            let mut axis = [py[0] - px[0], py[1] - px[1], py[2] - px[2]];
            let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
            for v in axis.iter_mut() {
                *v /= len;
            }
            let mut u = [pa[0] - px[0], pa[1] - px[1], pa[2] - px[2]];
            let mut w = [pb[0] - py[0], pb[1] - py[1], pb[2] - py[2]];
            let du = u[0] * axis[0] + u[1] * axis[1] + u[2] * axis[2];
            let dw = w[0] * axis[0] + w[1] * axis[1] + w[2] * axis[2];
            for r in 0..3 {
                u[r] -= du * axis[r];
                w[r] -= dw * axis[r];
            }
            let cis = u[0] * w[0] + u[1] * w[1] + u[2] * w[2] > 0.0;
            let (ka, kb) = (m.bond_between(x, a), m.bond_between(y, b));
            let first_of = |k: i32| {
                let e = &m.bonds[k as usize];
                if preorder[e.a as usize] < preorder[e.b as usize] {
                    e.a
                } else {
                    e.b
                }
            };
            let mut up_a = true;
            if let Some(&d) = directions.get(&ka) {
                up_a = (d == b'/') == (first_of(ka) == x);
            } else {
                directions.insert(ka, if first_of(ka) == x { b'/' } else { b'\\' });
            }
            let want = if cis { up_a } else { !up_a };
            let cb = if want == (first_of(kb) == y) { b'/' } else { b'\\' };
            if directions.get(&kb).is_some_and(|&d| d != cb) {
                continue;
            }
            directions.insert(kb, cb);
        }
    }

    // Stereocentres: four different neighbours (canonical ranks), sp3.
    let mut chiral_sign = vec![0i32; n];
    if stereo {
        for a in 0..n {
            let ai = a as i32;
            if !kept[a] || m.atoms[a].z <= 1 || hcount[a] > 1 {
                continue;
            }
            let nb = m.neighbours(ai);
            if nb.len() as i32 + m.atoms[a].hydrogens != 4 || m.atoms[a].hydrogens > 0 || hybridisation(m, ai) != 3 {
                continue;
            }
            let distinct: BTreeSet<i32> = nb.iter().map(|&x| rank[x as usize]).collect();
            if distinct.len() != 4 {
                continue;
            }
            // Order as written.
            let mut written = Vec::new();
            if parent[a] >= 0 {
                written.push(parent[a]);
            }
            if folded_h[a] >= 0 {
                written.push(folded_h[a]);
            }
            for &k in &closures[a] {
                written.push(m.other(k, ai));
            }
            written.extend(children[a].iter().copied());
            if written.len() != 4 {
                continue;
            }
            let v = signed_volume(m.atoms[a].position, position(ai, written[1]), position(ai, written[2]), position(ai, written[3]));
            if v.abs() < 1e-4 {
                continue;
            }
            chiral_sign[a] = if v > 0.0 { 2 } else { 1 };
        }
    }

    let atom_token = |a: usize| -> String {
        let atom = &m.atoms[a];
        let z = atom.z;
        let mut sum = 0;
        for &k in m.bonds_of(a as i32) {
            if kept[m.other(k, a as i32) as usize] && !m.bonds[k as usize].dative {
                sum += m.bonds[k as usize].order;
            }
        }
        let aromatic = atom.aromatic && aromatic_symbol(z);
        let mut symbol = if atom.dummy || z == 0 { "*".to_string() } else { symbol_of(z) };
        if aromatic {
            let first = symbol[..1].to_ascii_lowercase();
            symbol.replace_range(..1, &first);
        }
        let mut plain = organic_subset(z) && atom.charge == 0 && atom.isotope == 0 && chiral_sign[a] == 0 && !atom.dummy;
        if plain {
            if aromatic {
                let degree = kept_neighbours(a as i32).len();
                if z == 6 {
                    plain = hcount[a] == default_valence(6, 0, sum) - sum;
                } else if z == 7 || z == 15 {
                    plain = hcount[a] == 0 && (degree == 2 || degree == 3);
                } else if z == 8 || z == 16 {
                    plain = hcount[a] == 0 && degree == 2;
                } else if z == 5 {
                    plain = hcount[a] == 0;
                }
                if z == 7 && degree == 3 {
                    plain = plain && sum <= 3;
                }
            } else {
                let v = default_valence(z, 0, sum);
                plain = v >= sum && hcount[a] == v - sum;
            }
        }
        if atom.dummy && atom.charge == 0 && atom.isotope == 0 && hcount[a] == 0 {
            return "*".into();
        }
        if plain {
            return symbol;
        }
        let mut out = String::from("[");
        if atom.isotope != 0 {
            out += &atom.isotope.to_string();
        }
        out += &symbol;
        if chiral_sign[a] != 0 {
            out += if chiral_sign[a] == 2 { "@@" } else { "@" };
        }
        if hcount[a] == 1 {
            out += "H";
        } else if hcount[a] > 1 {
            out += &format!("H{}", hcount[a]);
        }
        if atom.charge > 0 {
            out += &if atom.charge == 1 { "+".to_string() } else { format!("+{}", atom.charge) };
        }
        if atom.charge < 0 {
            out += &if atom.charge == -1 { "-".to_string() } else { format!("-{}", -atom.charge) };
        }
        if atom.map_number != 0 {
            out += &format!(":{}", atom.map_number);
        }
        out + "]"
    };
    let bond_token = |k: i32, from: i32| -> String {
        let b = &m.bonds[k as usize];
        if let Some(&d) = directions.get(&k) {
            // Directions were defined for the bond written from its first atom.
            let first = if preorder[b.a as usize] < preorder[b.b as usize] { b.a } else { b.b };
            let d = if from == first { d } else if d == b'/' { b'\\' } else { b'/' };
            return (d as char).to_string();
        }
        let both_aromatic = m.atoms[b.a as usize].aromatic && m.atoms[b.b as usize].aromatic;
        if b.aromatic && both_aromatic {
            return String::new();
        }
        match b.order {
            2 => "=".into(),
            3 => "#".into(),
            o if o >= 4 => "$".into(),
            _ if both_aromatic => "-".into(),
            _ => String::new(),
        }
    };

    struct Writer<'a> {
        out: String,
        free_digits: BTreeSet<i32>,
        digit_of: BTreeMap<i32, i32>,
        closures: &'a [Vec<i32>],
        children: &'a [Vec<i32>],
        m: &'a Molecule,
    }
    fn digit_text(d: i32) -> String {
        if d < 10 {
            d.to_string()
        } else {
            format!("%{d}")
        }
    }
    fn write(w: &mut Writer, a: i32, atom_token: &dyn Fn(usize) -> String, bond_token: &dyn Fn(i32, i32) -> String) {
        w.out += &atom_token(a as usize);
        for &k in &w.closures[a as usize] {
            match w.digit_of.remove(&k) {
                None => {
                    let d = *w.free_digits.iter().next().expect("ring-closure digits");
                    w.free_digits.remove(&d);
                    w.digit_of.insert(k, d);
                    w.out += &bond_token(k, a);
                    w.out += &digit_text(d);
                }
                Some(d) => {
                    w.out += &digit_text(d);
                    w.free_digits.insert(d);
                }
            }
        }
        let ch = &w.children[a as usize];
        for c in 0..ch.len() {
            let branch = c + 1 < ch.len();
            if branch {
                w.out += "(";
            }
            w.out += &bond_token(w.m.bond_between(a, ch[c]), a);
            write(w, ch[c], atom_token, bond_token);
            if branch {
                w.out += ")";
            }
        }
    }
    let mut w = Writer { out: String::new(), free_digits: (1..100).collect(), digit_of: BTreeMap::new(), closures: &closures, children: &children, m };
    for (r, &root) in roots.iter().enumerate() {
        if r > 0 {
            w.out += ".";
        }
        write(&mut w, root, &atom_token, &bond_token);
    }
    w.out
}
