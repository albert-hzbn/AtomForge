//! Editing operations on molecules with coordinates (the C++
//! chemistry/Editing): substituting a hydrogen with a functional group,
//! changing bond orders with the hydrogens that follow, adding and removing
//! hydrogens, (de)protonation, inverting stereocentres and double bonds,
//! conformer searches, geometry scans, Gasteiger charges, superposition and RMSD.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use af_core::crmath;
use af_core::cxx::sort_by;
use af_core::cxx_random::{uniform_real, Mt19937};
use af_core::eigen::{determinant3, eigen_matmul3, eigen_mul3, jacobi_svd3, transpose3, IDENTITY3};

use crate::embedding::{molecule_from_smiles, EmbedOptions};
use crate::forcefield::{minimize, uff_bond_length, uff_type, EnergyTerms, ForceFieldOptions, Restraint, UffForceField};
use crate::molecule::*;

const PI: f64 = 3.14159265358979323846;

/// Rotation of v about the unit axis through the origin (Rodrigues).
fn rotate_about(v: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let (c, s) = if angle == PI { (crmath::cos(PI), crmath::sin(PI)) } else { (angle.cos(), angle.sin()) };
    add(add(scale(v, c), scale(cross(axis, v), s)), scale(axis, dot(axis, v) * (1.0 - c)))
}

/// Rotation taking unit vector a onto unit vector b, applied to v.
fn align_vector(v: Vec3, a: Vec3, b: Vec3) -> Vec3 {
    let axis = cross(a, b);
    let (s, c) = (norm(axis), dot(a, b));
    if s < 1e-9 {
        if c > 0.0 {
            return v;
        }
        // Opposite: half a turn about any perpendicular axis.
        let p = unit(cross(a, if a[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] }));
        return rotate_about(v, p, PI);
    }
    rotate_about(v, scale(axis, 1.0 / s), s.atan2(c))
}

/// Position of atom x as bonded to atom from (periodic images).
fn bonded_position(m: &Molecule, from: i32, x: i32) -> Vec3 {
    let k = m.bond_between(from, x);
    if k < 0 {
        return m.atoms[x as usize].position;
    }
    let b = &m.bonds[k as usize];
    let d = add(sub(m.atoms[b.b as usize].position, m.atoms[b.a as usize].position), b.image);
    add(m.atoms[from as usize].position, if b.a == from { d } else { scale(d, -1.0) })
}

fn contact_radius(z: i32) -> f64 {
    match z {
        1 => 1.1,
        6 => 1.7,
        7 => 1.55,
        8 => 1.52,
        _ => 1.8,
    }
}

fn heavy_bond_sum(m: &Molecule, atom: i32) -> i32 {
    m.bonds_of(atom)
        .iter()
        .filter(|&&k| !m.bonds[k as usize].dative && m.atoms[m.other(k, atom) as usize].z != 1)
        .map(|&k| m.bonds[k as usize].order)
        .sum()
}

fn hydrogen_neighbours(m: &Molecule, atom: i32) -> i32 {
    m.neighbours(atom).iter().filter(|&&x| m.atoms[x as usize].z == 1).count() as i32
}

fn place_hydrogens(m: &mut Molecule, atom: i32, count: i32) -> Result<(), String> {
    if count <= 0 {
        return Ok(());
    }
    let mut hyb = hybridisation(m, atom);
    if hyb == 0 {
        hyb = 3;
    }
    let dirs = new_bond_directions(m, atom, count, hyb);
    let length = hydrogen_bond_length(m.atoms[atom as usize].z);
    for d in dirs {
        let h = Atom { z: 1, position: add(m.atoms[atom as usize].position, scale(d, length)), ..Default::default() };
        let index = m.add_atom(h);
        m.add_bond(atom, index, 1)?;
    }
    Ok(())
}

/// Functional groups for substitution: name -> SMILES with '*' at the
/// attachment (methyl "*C", phenyl "*c1ccccc1", ...).
pub fn functional_groups() -> &'static [(&'static str, &'static str)] {
    &[
        ("methyl", "*C"), ("ethyl", "*CC"), ("propyl", "*CCC"), ("isopropyl", "*C(C)C"), ("butyl", "*CCCC"), ("isobutyl", "*CC(C)C"),
        ("sec-butyl", "*C(C)CC"), ("tert-butyl", "*C(C)(C)C"), ("cyclopropyl", "*C1CC1"), ("cyclopentyl", "*C1CCCC1"), ("cyclohexyl", "*C1CCCCC1"),
        ("vinyl", "*C=C"), ("allyl", "*CC=C"), ("ethynyl", "*C#C"), ("phenyl", "*c1ccccc1"), ("benzyl", "*Cc1ccccc1"),
        ("naphthyl", "*c1cccc2ccccc12"), ("2-pyridyl", "*c1ccccn1"), ("3-pyridyl", "*c1cccnc1"), ("4-pyridyl", "*c1ccncc1"),
        ("2-thienyl", "*c1cccs1"), ("2-furyl", "*c1ccco1"), ("imidazolyl", "*n1ccnc1"), ("tolyl", "*c1ccc(C)cc1"),
        ("mesityl", "*c1c(C)cc(C)cc1C"), ("hydroxyl", "*O"), ("methoxy", "*OC"), ("ethoxy", "*OCC"), ("phenoxy", "*Oc1ccccc1"),
        ("trifluoromethoxy", "*OC(F)(F)F"), ("amino", "*N"), ("methylamino", "*NC"), ("dimethylamino", "*N(C)C"), ("ammonium", "*[NH3+]"),
        ("anilino", "*Nc1ccccc1"), ("hydrazino", "*NN"), ("guanidino", "*NC(N)=N"), ("nitro", "*[N+](=O)[O-]"), ("nitroso", "*N=O"),
        ("cyano", "*C#N"), ("isocyano", "*[N+]#[C-]"), ("azido", "*N=[N+]=[N-]"), ("isocyanato", "*N=C=O"), ("isothiocyanato", "*N=C=S"),
        ("fluoro", "*F"), ("chloro", "*Cl"), ("bromo", "*Br"), ("iodo", "*I"), ("trifluoromethyl", "*C(F)(F)F"),
        ("trichloromethyl", "*C(Cl)(Cl)Cl"), ("pentafluorophenyl", "*c1c(F)c(F)c(F)c(F)c1F"), ("formyl", "*C=O"), ("acetyl", "*C(C)=O"),
        ("benzoyl", "*C(=O)c1ccccc1"), ("carboxyl", "*C(=O)O"), ("carboxylate", "*C(=O)[O-]"), ("methyl-ester", "*C(=O)OC"),
        ("ethyl-ester", "*C(=O)OCC"), ("acetoxy", "*OC(C)=O"), ("amide", "*C(N)=O"), ("acetamido", "*NC(C)=O"), ("boc", "*C(=O)OC(C)(C)C"),
        ("hydroxymethyl", "*CO"), ("aminomethyl", "*CN"), ("thiol", "*S"), ("methylthio", "*SC"), ("sulfonic-acid", "*S(=O)(=O)O"),
        ("sulfonate", "*S(=O)(=O)[O-]"), ("sulfonamide", "*S(N)(=O)=O"), ("methylsulfonyl", "*S(C)(=O)=O"),
        ("tosyl", "*S(=O)(=O)c1ccc(C)cc1"), ("sulfonyl-chloride", "*S(Cl)(=O)=O"), ("phosphate", "*OP(=O)(O)O"), ("phosphonate", "*P(=O)(O)O"),
        ("diphenylphosphino", "*P(c1ccccc1)c1ccccc1"), ("boronic-acid", "*B(O)O"), ("trimethylsilyl", "*[Si](C)(C)C"),
        ("morpholino", "*N1CCOCC1"), ("piperidino", "*N1CCCCC1"), ("pyrrolidino", "*N1CCCC1"), ("piperazinyl", "*N1CCNCC1"),
        ("hydrogen", "*[H]"),
    ]
}

/// The SMILES of a group name, or the text itself when it contains '*'.
pub fn group_smiles(name: &str) -> Result<String, String> {
    if name.contains('*') {
        return Ok(name.to_string());
    }
    let mut key: String = name.bytes().map(|c| (if c == b'_' || c == b' ' { b'-' } else { c }).to_ascii_lowercase() as char).collect();
    static ALIASES: OnceLock<BTreeMap<&'static str, &'static str>> = OnceLock::new();
    let aliases = ALIASES.get_or_init(|| {
        [
            ("me", "methyl"), ("et", "ethyl"), ("pr", "propyl"), ("ipr", "isopropyl"), ("bu", "butyl"), ("tbu", "tert-butyl"),
            ("t-butyl", "tert-butyl"), ("ph", "phenyl"), ("bn", "benzyl"), ("oh", "hydroxyl"), ("ome", "methoxy"), ("oet", "ethoxy"),
            ("nh2", "amino"), ("nme2", "dimethylamino"), ("no2", "nitro"), ("cn", "cyano"), ("f", "fluoro"), ("cl", "chloro"),
            ("br", "bromo"), ("i", "iodo"), ("cf3", "trifluoromethyl"), ("cho", "formyl"), ("ac", "acetyl"), ("cooh", "carboxyl"),
            ("coo-", "carboxylate"), ("cooh-", "carboxylate"), ("co2me", "methyl-ester"), ("oac", "acetoxy"), ("conh2", "amide"),
            ("nhac", "acetamido"), ("sh", "thiol"), ("sme", "methylthio"), ("so3h", "sulfonic-acid"), ("ts", "tosyl"),
            ("tms", "trimethylsilyl"), ("h", "hydrogen"), ("hydroxy", "hydroxyl"), ("ch3", "methyl"), ("c2h5", "ethyl"),
            ("ocf3", "trifluoromethoxy"), ("nco", "isocyanato"), ("ncs", "isothiocyanato"), ("n3", "azido"), ("pph2", "diphenylphosphino"),
            ("b(oh)2", "boronic-acid"), ("pyridyl", "4-pyridyl"), ("thienyl", "2-thienyl"), ("furyl", "2-furyl"), ("tolyl", "tolyl"),
        ]
        .into_iter()
        .collect()
    });
    if let Some(alias) = aliases.get(key.as_str()) {
        key = alias.to_string();
    }
    if let Some((_, smiles)) = functional_groups().iter().find(|(n, _)| *n == key) {
        return Ok(smiles.to_string());
    }
    let names: Vec<&str> = functional_groups().iter().map(|g| g.0).collect();
    Err(format!("unknown group '{name}'; give a SMILES with * at the attachment or one of: {}", names.join(", ")))
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SubstituteResult {
    /// Indices of the group's atoms in the result.
    pub new_atoms: Vec<i32>,
    /// The group atom bonded to the host.
    pub attachment: i32,
}

/// Replaces terminal atom `atom` (usually a hydrogen) by the group, bonded to
/// the atom it was attached to with the UFF bond length; the group is turned
/// about the new bond to the position with the fewest clashes.
pub fn substitute(m: &mut Molecule, atom: i32, group: &str) -> Result<SubstituteResult, String> {
    if atom < 0 || atom as usize >= m.atoms.len() {
        return Err(format!("no atom {atom}"));
    }
    if m.degree(atom) != 1 {
        return Err(format!("atom {atom} is not terminal (it has {} bonds); substitute a hydrogen or another end atom", m.degree(atom)));
    }
    let host = m.neighbours(atom)[0];
    let host_bond = m.bond_between(host, atom);
    let host_order = m.bonds[host_bond as usize].order;
    let smiles = group_smiles(group)?;
    let options = EmbedOptions { attempts: 2, ..Default::default() };
    let (g, _) = molecule_from_smiles(&smiles, &options).map_err(|e| e.message)?;
    let mut dummy = -1i32;
    for (i, a) in g.atoms.iter().enumerate() {
        if a.dummy {
            if dummy >= 0 {
                return Err("the group has more than one '*'".into());
            }
            dummy = i as i32;
        }
    }
    if dummy < 0 {
        return Err("the group needs a '*' at its attachment".into());
    }
    if g.degree(dummy) != 1 {
        return Err("'*' must have one bond".into());
    }
    let attach = g.neighbours(dummy)[0];
    let group_order = g.bonds[g.bond_between(dummy, attach) as usize].order;
    if group_order != host_order {
        return Err(format!("the group's bond order ({group_order}) differs from the bond it replaces ({host_order})"));
    }
    // Bond length from UFF types (the host's type before the change).
    let mut host_type = String::new();
    let mut attach_type = String::new();
    if let Ok(t) = uff_type(m, host) {
        host_type = t;
        if let Ok(t) = uff_type(&g, attach) {
            attach_type = t;
        }
    }
    let mut length = covalent_radius(m.atoms[host as usize].z) + covalent_radius(g.atoms[attach as usize].z);
    if !host_type.is_empty() && !attach_type.is_empty() {
        length = uff_bond_length(&host_type, &attach_type, f64::from(group_order))?;
    }
    let p = m.atoms[host as usize].position;
    let u = unit(sub(bonded_position(m, host, atom), p));
    let v = unit(sub(g.atoms[attach as usize].position, g.atoms[dummy as usize].position));
    let a0 = g.atoms[attach as usize].position;
    let placed: Vec<Vec3> = g.atoms.iter().map(|a| add(add(p, scale(u, length)), align_vector(sub(a.position, a0), v, u))).collect();
    // Turn about the new bond to avoid clashes with the host.
    let pivot = add(p, scale(u, length));
    let (mut best_score, mut best_angle) = (1e300, 0.0);
    for step in 0..36 {
        let angle = f64::from(step) * PI / 18.0;
        let mut score = 0.0;
        for i in 0..g.atoms.len() {
            if i as i32 == dummy || i as i32 == attach {
                continue;
            }
            let x = add(pivot, rotate_about(sub(placed[i], pivot), u, angle));
            for h in 0..m.atoms.len() {
                if h as i32 == atom || h as i32 == host {
                    continue;
                }
                let d = norm(sub(x, m.atoms[h].position));
                let r = 0.8 * (contact_radius(g.atoms[i].z) + contact_radius(m.atoms[h].z));
                if d < r {
                    score += (r - d) * (r - d);
                }
            }
        }
        if score < best_score - 1e-9 {
            best_score = score;
            best_angle = angle;
        }
    }
    let mut map = vec![-1i32; g.atoms.len()];
    for i in 0..g.atoms.len() {
        if i as i32 == dummy {
            continue;
        }
        let mut a = g.atoms[i].clone();
        a.dummy = false;
        a.position = add(pivot, rotate_about(sub(placed[i], pivot), u, best_angle));
        map[i] = m.add_atom(a);
    }
    for b in &g.bonds {
        if b.a == dummy || b.b == dummy {
            continue;
        }
        let k = m.add_bond(map[b.a as usize], map[b.b as usize], b.order)?;
        m.bonds[k as usize].aromatic = b.aromatic;
    }
    m.add_bond(host, map[attach as usize], group_order)?;
    // Remove the replaced atom.
    let mut keep = vec![1u8; m.atoms.len()];
    keep[atom as usize] = 0;
    let renumber = m.remove_atoms(&keep);
    let mut result = SubstituteResult::default();
    for &x in &map {
        if x >= 0 {
            result.new_atoms.push(renumber[x as usize]);
        }
    }
    result.attachment = renumber[map[attach as usize] as usize];
    Ok(result)
}

/// Sets the order of the bond a-b (creating it when absent; 0 removes it) and
/// adds or removes hydrogens on a and b to keep their valences. Returns the
/// change in the number of hydrogens.
pub fn set_bond_order(m: &mut Molecule, mut a: i32, mut b: i32, order: i32) -> Result<i32, String> {
    let n = m.atoms.len() as i32;
    if a < 0 || b < 0 || a >= n || b >= n || a == b {
        return Err("give two different atoms".into());
    }
    if (m.atoms[a as usize].z == 1 || m.atoms[b as usize].z == 1) && order > 1 {
        return Err("hydrogen takes single bonds only".into());
    }
    let k = m.bond_between(a, b);
    if order <= 0 {
        if k < 0 {
            return Ok(0);
        }
        m.bonds.remove(k as usize);
        m.reindex();
    } else if k < 0 {
        m.add_bond(a, b, order)?;
    } else {
        m.bonds[k as usize].order = order;
        m.bonds[k as usize].aromatic = false;
    }
    // Hydrogens to match the valences of a and b: removals in one pass, then additions.
    let mut change = 0;
    let mut keep = vec![1u8; m.atoms.len()];
    for x in [a, b] {
        let ax = &m.atoms[x as usize];
        if ax.z == 1 || is_metal(ax.z) {
            continue;
        }
        let heavy = heavy_bond_sum(m, x);
        let target = default_valence(ax.z, ax.charge, heavy);
        if target < 0 {
            continue;
        }
        let mut surplus = hydrogen_neighbours(m, x) - 0.max(target - heavy);
        for y in m.neighbours(x) {
            if surplus > 0 && m.atoms[y as usize].z == 1 && m.degree(y) == 1 && keep[y as usize] != 0 {
                keep[y as usize] = 0;
                surplus -= 1;
                change -= 1;
            }
        }
    }
    if change < 0 {
        let map = m.remove_atoms(&keep);
        a = map[a as usize];
        b = map[b as usize];
    }
    for x in [a, b] {
        let (z, charge) = (m.atoms[x as usize].z, m.atoms[x as usize].charge);
        if z == 1 || is_metal(z) {
            continue;
        }
        let heavy = heavy_bond_sum(m, x);
        let target = default_valence(z, charge, heavy);
        if target < 0 {
            continue;
        }
        let have = hydrogen_neighbours(m, x);
        if target - heavy > have {
            place_hydrogens(m, x, target - heavy - have)?;
            change += target - heavy - have;
        }
    }
    Ok(change)
}

/// Adds hydrogens to atoms (mask, or all) that have fewer bonds than their
/// valence; returns the number added.
pub fn add_hydrogens(m: &mut Molecule, mask: &[u8]) -> Result<i32, String> {
    let mut added = 0;
    let n = m.atoms.len();
    for i in 0..n {
        if !mask.is_empty() && mask[i] == 0 {
            continue;
        }
        let (z, charge) = (m.atoms[i].z, m.atoms[i].charge);
        if z <= 1 || is_metal(z) {
            continue;
        }
        let sum = m.valence(i as i32);
        let target = default_valence(z, charge, sum);
        if target > sum {
            place_hydrogens(m, i as i32, target - sum)?;
            added += target - sum;
        }
    }
    Ok(added)
}

/// Removes hydrogens bonded to the atoms in the mask (or all hydrogens).
pub fn remove_hydrogens(m: &mut Molecule, mask: &[u8]) -> i32 {
    let mut keep = vec![1u8; m.atoms.len()];
    let mut removed = 0;
    for i in 0..m.atoms.len() {
        if m.atoms[i].z != 1 {
            continue;
        }
        let mut hit = mask.is_empty();
        if !mask.is_empty() && m.neighbours(i as i32).iter().any(|&x| mask[x as usize] != 0) {
            hit = true;
        }
        if !mask.is_empty() && mask[i] != 0 {
            hit = true;
        }
        if hit {
            keep[i] = 0;
            removed += 1;
        }
    }
    if removed > 0 {
        m.remove_atoms(&keep);
    }
    removed
}

/// Adds a proton to an atom with a lone pair (N, O, S, P): charge +1.
pub fn protonate(m: &mut Molecule, atom: i32) -> Result<(), String> {
    let z = m.atoms[atom as usize].z;
    if !matches!(z, 7 | 8 | 15 | 16 | 6 | 34 | 33) {
        return Err(format!("atom {atom} ({}) has no lone pair to protonate", symbol_of(z)));
    }
    let mut hyb = hybridisation(m, atom);
    // An sp3 amine becomes ammonium (tetrahedral); an sp2 N (pyridine, imine) keeps its plane.
    if hyb == 0 {
        hyb = 3;
    }
    if m.atoms[atom as usize].aromatic || hyb == 2 {
        hyb = 2;
    }
    let degree = m.degree(atom);
    if (hyb == 2 && degree >= 3) || degree >= 4 {
        hyb = 3;
    }
    let dirs = new_bond_directions(m, atom, 1, hyb);
    let h = Atom { z: 1, position: add(m.atoms[atom as usize].position, scale(dirs[0], hydrogen_bond_length(z))), ..Default::default() };
    let index = m.add_atom(h);
    m.add_bond(atom, index, 1)?;
    m.atoms[atom as usize].charge += 1;
    Ok(())
}

/// Removes a proton (a bonded H) from an atom: charge -1.
pub fn deprotonate(m: &mut Molecule, atom: i32) -> Result<(), String> {
    let mut hydrogen = -1;
    for x in m.neighbours(atom) {
        if m.atoms[x as usize].z == 1 && m.degree(x) == 1 {
            hydrogen = x;
        }
    }
    if hydrogen < 0 {
        return Err(format!("atom {atom} has no hydrogen to remove"));
    }
    let mut keep = vec![1u8; m.atoms.len()];
    keep[hydrogen as usize] = 0;
    m.atoms[atom as usize].charge -= 1;
    m.remove_atoms(&keep);
    Ok(())
}

/// Ionisation at about pH 7: carboxylic, sulfonic and phosphoric acids lose a
/// proton, aliphatic amines, amidines and guanidines gain one. Returns a
/// description of the changes.
pub fn ionize_at_neutral_ph(m: &mut Molecule) -> Result<String, String> {
    let (mut acids, mut sulfonic, mut phosphoric, mut amines, mut amidines) = (0, 0, 0, 0, 0);
    // Collect the sites first, then edit (indices shift on removal).
    let mut deprotonate_at: Vec<i32> = Vec::new();
    let mut protonate_at: Vec<i32> = Vec::new();
    let mut done_p: BTreeSet<i32> = BTreeSet::new();
    let n = m.atoms.len() as i32;
    for o in 0..n {
        let ao = &m.atoms[o as usize];
        if ao.z != 8 || ao.charge != 0 || m.degree(o) != 2 || hydrogen_neighbours(m, o) != 1 {
            continue;
        }
        let mut centre = -1;
        for x in m.neighbours(o) {
            if m.atoms[x as usize].z != 1 {
                centre = x;
            }
        }
        if centre < 0 {
            continue;
        }
        let oxo = m.bonds_of(centre).iter().any(|&k| {
            let y = m.other(k, centre);
            y != o && m.bonds[k as usize].order == 2 && matches!(m.atoms[y as usize].z, 8 | 16)
        });
        if !oxo {
            continue;
        }
        match m.atoms[centre as usize].z {
            6 => {
                deprotonate_at.push(o);
                acids += 1;
            }
            16 => {
                deprotonate_at.push(o);
                sulfonic += 1;
            }
            15 if !done_p.contains(&centre) => {
                deprotonate_at.push(o);
                done_p.insert(centre);
                phosphoric += 1;
            }
            _ => {}
        }
    }
    for n_atom in 0..n {
        let a = &m.atoms[n_atom as usize];
        if a.z != 7 || a.charge != 0 || a.aromatic {
            continue;
        }
        let (mut skip, mut imine) = (false, false);
        let mut heavy = 0;
        for &k in m.bonds_of(n_atom) {
            let y = m.other(k, n_atom);
            if m.atoms[y as usize].z == 1 {
                continue;
            }
            heavy += 1;
            let bk = &m.bonds[k as usize];
            if bk.order == 2 && m.atoms[y as usize].z == 6 {
                imine = true;
            }
            if bk.order >= 2 && m.atoms[y as usize].z != 6 {
                skip = true;
            }
            if m.atoms[y as usize].z != 6 {
                skip = true; // N-N, N-O, N-S
            }
            if bk.order == 1 {
                for &e in m.bonds_of(y) {
                    let w = m.other(e, y);
                    if w == n_atom {
                        continue;
                    }
                    // Amides, anilines, sulfonamides, enamines: not basic.
                    if (m.bonds[e as usize].order >= 2 && !imine) || m.bonds[e as usize].aromatic || m.atoms[y as usize].aromatic {
                        skip = true;
                    }
                }
            }
        }
        if heavy == 0 {
            continue;
        }
        if imine {
            // Amidines and guanidines: C(=N)N.
            let mut carbon = -1;
            for &k in m.bonds_of(n_atom) {
                if m.bonds[k as usize].order == 2 {
                    carbon = m.other(k, n_atom);
                }
            }
            let mut amidine = false;
            if carbon >= 0 {
                for &k in m.bonds_of(carbon) {
                    let y = m.other(k, carbon);
                    if y != n_atom && m.atoms[y as usize].z == 7 && m.bonds[k as usize].order == 1 {
                        amidine = true;
                    }
                }
            }
            if amidine && !skip {
                protonate_at.push(n_atom);
                amidines += 1;
            }
            continue;
        }
        if skip || m.degree(n_atom) != 3 || hybridisation(m, n_atom) != 3 {
            continue;
        }
        protonate_at.push(n_atom);
        amines += 1;
    }
    for &x in &protonate_at {
        protonate(m, x)?;
    }
    deprotonate_at.sort_unstable_by(|a, b| b.cmp(a));
    for &x in &deprotonate_at {
        deprotonate(m, x)?;
    }
    let mut out = String::new();
    let mut put = |count: i32, what: &str| {
        if count == 0 {
            return;
        }
        if !out.is_empty() {
            out += ", ";
        }
        out += &format!("{count} {what}");
    };
    put(acids, "carboxylic acid(s) deprotonated");
    put(sulfonic, "sulfonic or sulfuric acid(s) deprotonated");
    put(phosphoric, "phosphate or phosphonate group(s) deprotonated once");
    put(amines, "amine(s) protonated");
    put(amidines, "amidine or guanidine group(s) protonated");
    Ok(if out.is_empty() { "no ionisable groups changed".into() } else { out })
}

/// Atoms on b's side of the bond a-b (excluding a); empty when a-b is in a ring.
pub fn side_of(m: &Molecule, a: i32, b: i32) -> Vec<i32> {
    let mut seen = vec![false; m.atoms.len()];
    let mut out = vec![b];
    let mut stack = vec![b];
    seen[a as usize] = true;
    seen[b as usize] = true;
    while let Some(x) = stack.pop() {
        for &k in m.bonds_of(x) {
            let y = m.other(k, x);
            if x == b && y == a {
                continue;
            }
            if y == a {
                return Vec::new(); // a ring through a-b
            }
            if seen[y as usize] {
                continue;
            }
            seen[y as usize] = true;
            out.push(y);
            stack.push(y);
        }
    }
    out
}

/// Dihedral a-b-c-d in degrees.
pub fn dihedral_angle(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f64 {
    let (b1, b2, b3) = (sub(b, a), sub(c, b), sub(d, c));
    let (n1, n2) = (cross(b1, b2), cross(b2, b3));
    (norm(b2) * dot(b1, n2)).atan2(dot(n1, n2)) * 180.0 / PI
}

/// Sets a dihedral i-j-k-l by turning l's side about j-k (rigid).
pub fn set_dihedral(m: &mut Molecule, i: i32, j: i32, k: i32, l: i32, degrees: f64) -> Result<(), String> {
    let side = side_of(m, j, k);
    if side.is_empty() {
        return Err(format!("the bond {j}-{k} is in a ring"));
    }
    let p = |x: i32| m.atoms[x as usize].position;
    let current = dihedral_angle(p(i), p(j), p(k), p(l));
    let axis = unit(sub(p(k), p(j)));
    let turn = (degrees - current) * PI / 180.0;
    let origin = p(k);
    for x in side {
        let pos = &mut m.atoms[x as usize].position;
        *pos = add(origin, rotate_about(sub(*pos, origin), axis, turn));
    }
    Ok(())
}

/// Inverts a tetrahedral centre: two branches that are not in a ring swap by a
/// 180-degree turn about their bisector, or (in rings) the centre and a
/// terminal substituent are reflected through the plane of the other three.
pub fn invert_stereocentre(m: &mut Molecule, atom: i32) -> Result<(), String> {
    let nb = m.neighbours(atom);
    if nb.len() != 4 {
        return Err(format!("atom {atom} is not a tetrahedral centre (it has {} bonds)", nb.len()));
    }
    let c = m.atoms[atom as usize].position;
    let mut acyclic: Vec<(usize, i32)> = Vec::new();
    for &x in &nb {
        let side = side_of(m, atom, x);
        if !side.is_empty() {
            acyclic.push((side.len(), x));
        }
    }
    if acyclic.len() >= 2 {
        acyclic.sort();
        let (p, q) = (acyclic[0].1, acyclic[1].1);
        let axis = unit(add(unit(sub(m.atoms[p as usize].position, c)), unit(sub(m.atoms[q as usize].position, c))));
        for branch in [p, q] {
            for x in side_of(m, atom, branch) {
                let pos = &mut m.atoms[x as usize].position;
                *pos = add(c, rotate_about(sub(*pos, c), axis, PI));
            }
        }
        return Ok(());
    }
    // In rings: reflect the centre and a terminal neighbour through the plane of the other three.
    let mut terminal = -1;
    for &x in &nb {
        if m.degree(x) == 1 {
            terminal = x;
        }
    }
    if terminal < 0 {
        return Err(format!("atom {atom} is a bridgehead (no free branch); it cannot be inverted alone"));
    }
    let others: Vec<Vec3> = nb.iter().filter(|&&x| x != terminal).map(|&x| m.atoms[x as usize].position).collect();
    let normal = unit(cross(sub(others[1], others[0]), sub(others[2], others[0])));
    let reflect = |x: Vec3| sub(x, scale(normal, 2.0 * dot(sub(x, others[0]), normal)));
    m.atoms[atom as usize].position = reflect(m.atoms[atom as usize].position);
    m.atoms[terminal as usize].position = reflect(m.atoms[terminal as usize].position);
    Ok(())
}

/// Turns one side of the double bond a=b by 180 degrees (E <-> Z).
pub fn flip_double_bond(m: &mut Molecule, a: i32, b: i32) -> Result<(), String> {
    if m.bond_between(a, b) < 0 {
        return Err(format!("atoms {a} and {b} are not bonded"));
    }
    let side = side_of(m, a, b);
    if side.is_empty() {
        return Err("the bond is in a ring; it cannot be flipped".into());
    }
    let origin = m.atoms[b as usize].position;
    let axis = unit(sub(origin, m.atoms[a as usize].position));
    for x in side {
        if x != b {
            let pos = &mut m.atoms[x as usize].position;
            *pos = add(origin, rotate_about(sub(*pos, origin), axis, PI));
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Default)]
struct Abc {
    a: f64,
    b: f64,
    c: f64,
}

/// Gasteiger-Marsili partial charges (PEOE; Tetrahedron 36, 3219 (1980)),
/// starting from the formal charges; atoms without parameters keep theirs.
pub fn gasteiger_charges(m: &Molecule, iterations: i32) -> Vec<f64> {
    // chi(q) = a + b q + c q^2 for each element and hybridisation.
    let conjugated_lone_pair = |i: i32| {
        // O or S with single bonds next to a pi system (esters, acids, phenols): sp2-like.
        m.bonds_of(i).iter().any(|&k| {
            let x = m.other(k, i);
            m.bonds_of(x).iter().any(|&e| e != k && (m.bonds[e as usize].order >= 2 || m.bonds[e as usize].aromatic))
        })
    };
    let abc = |a, b, c| Abc { a, b, c };
    let parameters = |i: i32| -> Option<Abc> {
        let z = m.atoms[i as usize].z;
        let mut hyb = hybridisation(m, i);
        if (z == 8 || z == 16) && hyb == 3 && conjugated_lone_pair(i) {
            hyb = 2;
        }
        Some(match z {
            1 => abc(7.17, 6.24, -0.56),
            5 => abc(5.98, 6.82, 1.605),
            6 => match hyb {
                1 => abc(10.39, 9.45, 0.73),
                2 => abc(8.79, 9.32, 1.51),
                _ => abc(7.98, 9.18, 1.88),
            },
            7 => match hyb {
                1 => abc(15.68, 11.70, -0.27),
                2 => abc(12.87, 11.15, 0.85),
                _ => abc(11.54, 10.82, 1.36),
            },
            8 => {
                if hyb == 2 || hyb == 1 {
                    abc(17.07, 13.79, 0.47)
                } else {
                    abc(14.18, 12.92, 1.39)
                }
            }
            9 => abc(14.66, 13.85, 2.31),
            14 => abc(7.30, 6.567, 0.657),
            15 => abc(8.90, 8.24, 0.96),
            16 => {
                if hyb == 2 {
                    abc(10.88, 9.485, 1.325)
                } else {
                    abc(10.14, 9.13, 1.38)
                }
            }
            17 => abc(11.00, 9.69, 1.35),
            35 => abc(10.08, 8.47, 1.16),
            53 => abc(9.90, 7.96, 0.96),
            _ => return None,
        })
    };
    let n = m.atoms.len();
    let mut p = vec![Abc::default(); n];
    let mut known = vec![false; n];
    let mut q = vec![0.0; n];
    let mut cation = vec![0.0; n];
    for i in 0..n {
        if let Some(v) = parameters(i as i32) {
            p[i] = v;
            known[i] = true;
        }
        q[i] = f64::from(m.atoms[i].charge);
        cation[i] = if m.atoms[i].z == 1 { 20.02 } else { p[i].a + p[i].b + p[i].c };
    }
    // Resonance: terminal atoms of one element on the same centre share their
    // formal charges (both O of a carboxylate or nitro group, the O of a sulfonate).
    for c in 0..n as i32 {
        let mut terminal: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
        for x in m.neighbours(c) {
            if m.degree(x) == 1 && m.atoms[x as usize].z != 1 {
                terminal.entry(m.atoms[x as usize].z).or_default().push(x);
            }
        }
        for atoms in terminal.values() {
            if atoms.len() < 2 {
                continue;
            }
            let mut sum = 0.0;
            let mut charged = false;
            for &x in atoms {
                sum += q[x as usize];
                charged = charged || q[x as usize] != 0.0;
            }
            if !charged {
                continue;
            }
            for &x in atoms {
                q[x as usize] = sum / atoms.len() as f64;
            }
        }
    }
    let mut damping = 1.0;
    for _ in 0..iterations {
        damping *= 0.5;
        let chi: Vec<f64> = (0..n).map(|i| p[i].a + p[i].b * q[i] + p[i].c * q[i] * q[i]).collect();
        let mut dq = vec![0.0; n];
        for b in &m.bonds {
            let (ba, bb) = (b.a as usize, b.b as usize);
            if b.dative || !known[ba] || !known[bb] {
                continue;
            }
            let lo = if chi[ba] < chi[bb] { ba } else { bb };
            let hi = if lo == ba { bb } else { ba };
            // Electrons move to the more electronegative atom.
            let transfer = (chi[hi] - chi[lo]) / cation[lo] * damping;
            dq[lo] += transfer;
            dq[hi] -= transfer;
        }
        for i in 0..n {
            q[i] += dq[i];
        }
    }
    q
}

/// Root-mean-square deviation after optimal superposition (Kabsch) of the
/// given atoms (all when empty); `aligned` receives b moved onto a.
pub fn rmsd(a: &[Vec3], b: &[Vec3], atoms: &[i32], aligned: Option<&mut Vec<Vec3>>) -> f64 {
    let mut use_: Vec<i32> = atoms.to_vec();
    if use_.is_empty() {
        use_ = (0..a.len() as i32).collect();
    }
    if use_.is_empty() {
        return 0.0;
    }
    let (mut ca, mut cb) = ([0.0; 3], [0.0; 3]);
    for &i in &use_ {
        ca = add(ca, a[i as usize]);
        cb = add(cb, b[i as usize]);
    }
    let count = use_.len() as f64;
    ca = scale(ca, 1.0 / count);
    cb = scale(cb, 1.0 / count);
    let mut h = [[0.0; 3]; 3];
    for &i in &use_ {
        let (p, q) = (sub(b[i as usize], cb), sub(a[i as usize], ca));
        for r in 0..3 {
            for s in 0..3 {
                h[r][s] += p[r] * q[s];
            }
        }
    }
    let (u, _, v) = jacobi_svd3(&h);
    let ut = transpose3(&u);
    let d = if determinant3(&eigen_matmul3(&v, &ut)) > 0.0 { 1.0 } else { -1.0 };
    let mut diag = IDENTITY3;
    diag[2][2] = d;
    let r = eigen_matmul3(&eigen_matmul3(&v, &diag), &ut);
    let mv = |x: Vec3| {
        let q = eigen_mul3(&r, [x[0] - cb[0], x[1] - cb[1], x[2] - cb[2]]);
        [q[0] + ca[0], q[1] + ca[1], q[2] + ca[2]]
    };
    let mut sum = 0.0;
    for &i in &use_ {
        let diff = sub(mv(b[i as usize]), a[i as usize]);
        sum += dot(diff, diff);
    }
    if let Some(out) = aligned {
        *out = b.iter().map(|&x| mv(x)).collect();
    }
    (sum / count).sqrt()
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Conformer {
    pub positions: Vec<Vec3>,
    /// kcal/mol.
    pub energy: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ConformerOptions {
    /// random or systematic.
    pub method: String,
    /// random: starting points.
    pub trials: i32,
    /// systematic: angles per rotatable bond (3 = 60, 180, 300).
    pub steps: i32,
    pub max_conformers: i32,
    /// A, heavy atoms after superposition.
    pub rmsd_threshold: f64,
    /// kcal/mol above the lowest.
    pub energy_window: f64,
    pub seed: u32,
    pub optimize_steps: i32,
    pub fixed: Vec<u8>,
}

impl Default for ConformerOptions {
    fn default() -> Self {
        Self {
            method: "random".into(),
            trials: 100,
            steps: 3,
            max_conformers: 10,
            rmsd_threshold: 0.5,
            energy_window: 10.0,
            seed: 1,
            optimize_steps: 1000,
            fixed: Vec::new(),
        }
    }
}

/// Distinct low-energy conformers (lowest first) from torsion changes about
/// the rotatable bonds and UFF minimisation; also the number of rotors.
pub fn conformer_search(molecule: &Molecule, options: &ConformerOptions) -> Result<(Vec<Conformer>, i32), String> {
    let m = molecule.clone();
    let field = UffForceField::new(&m, None, ForceFieldOptions::default())?;
    // Rotatable bonds with reference atoms and the smaller side to turn.
    let mut rotors: Vec<[i32; 4]> = Vec::new();
    for kb in rotatable_bonds(&m) {
        let (mut j, mut k) = (m.bonds[kb as usize].a, m.bonds[kb as usize].b);
        if side_of(&m, j, k).len() > side_of(&m, k, j).len() {
            std::mem::swap(&mut j, &mut k);
        }
        // Turn the side of j (smaller): dihedral l-k-j-i with i on j's side.
        let (mut i, mut l) = (-1i32, -1i32);
        for x in m.neighbours(j) {
            if x != k && (i < 0 || m.atoms[x as usize].z > m.atoms[i as usize].z) {
                i = x;
            }
        }
        for x in m.neighbours(k) {
            if x != j && (l < 0 || m.atoms[x as usize].z > m.atoms[l as usize].z) {
                l = x;
            }
        }
        if i < 0 || l < 0 {
            continue;
        }
        rotors.push([l, k, j, i]);
    }
    let rotatable = rotors.len() as i32;
    let heavy: Vec<i32> = (0..m.atoms.len() as i32).filter(|&a| m.atoms[a as usize].z != 1).collect();
    let mut found: Vec<Conformer> = Vec::new();
    let mut relax = |trial: &Molecule| -> Result<(), String> {
        let mut x: Vec<Vec3> = trial.atoms.iter().map(|a| a.position).collect();
        let r = minimize(&field, &mut x, &options.fixed, options.optimize_steps, 0.02)?;
        found.push(Conformer { positions: x, energy: r.energy });
        Ok(())
    };
    relax(&m)?;
    let mut random = Mt19937::new(options.seed);
    let systematic = options.method == "systematic";
    let mut combinations: usize = 1;
    if systematic {
        let mut r = 0;
        while r < rotors.len() && combinations <= 100000 {
            combinations *= options.steps.max(1) as usize;
            r += 1;
        }
    }
    if systematic && combinations <= 100000 {
        for c in 0..combinations {
            let mut trial = m.clone();
            let mut code = c;
            for r in &rotors {
                let choice = (code % options.steps as usize) as i32;
                code /= options.steps as usize;
                let angle = 60.0 + 360.0 * f64::from(choice) / f64::from(options.steps);
                set_dihedral(&mut trial, r[0], r[1], r[2], r[3], angle)?;
            }
            relax(&trial)?;
        }
    } else {
        let trials = (if systematic { options.trials.max(1000) } else { options.trials }).max(1);
        let mut t = 0;
        while t < trials && !rotors.is_empty() {
            let mut trial = m.clone();
            for r in &rotors {
                let angle = uniform_real(&mut random, -180.0, 180.0);
                set_dihedral(&mut trial, r[0], r[1], r[2], r[3], angle)?;
            }
            relax(&trial)?;
            t += 1;
        }
    }
    sort_by(&mut found, |a, b| a.energy < b.energy);
    let mut kept: Vec<Conformer> = Vec::new();
    for c in &found {
        if c.energy - found[0].energy > options.energy_window {
            break;
        }
        if !kept.iter().any(|k| rmsd(&k.positions, &c.positions, &heavy, None) < options.rmsd_threshold) {
            kept.push(c.clone());
        }
        if kept.len() as i32 >= options.max_conformers {
            break;
        }
    }
    Ok((kept, rotatable))
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScanPoint {
    /// Distance (A) or angle (degrees).
    pub value: f64,
    /// kcal/mol.
    pub energy: f64,
    pub positions: Vec<Vec3>,
}

/// Scans a distance (2 atoms), angle (3) or dihedral (4) from `start` to
/// `stop` in `step`; relaxed (UFF with the coordinate restrained) or rigid.
pub fn geometry_scan(molecule: &Molecule, atoms: &[i32], start: f64, stop: f64, step: f64, relaxed: bool, optimize_steps: i32) -> Result<Vec<ScanPoint>, String> {
    if atoms.len() < 2 || atoms.len() > 4 {
        return Err("a scan needs 2 (distance), 3 (angle) or 4 (dihedral) atoms".into());
    }
    if step == 0.0 {
        return Err("the scan step must not be zero".into());
    }
    let s = (if stop >= start { 1.0 } else { -1.0 }) * step.abs();
    let count = ((stop - start) / s + 1e-9).floor() as i32 + 1;
    if count > 2000 {
        return Err("the scan has more than 2000 points".into());
    }
    let mut m = molecule.clone();
    let mut out = Vec::new();
    let field = UffForceField::new(&m, None, ForceFieldOptions::default())?;
    let original = m.clone();
    let a = atoms;
    for p in 0..count {
        let value = start + f64::from(p) * s;
        if !relaxed {
            m = original.clone();
        }
        // Rigid move of the last atom's side.
        let moved = (|| -> Result<(), String> {
            let pos = |m: &Molecule, x: i32| m.atoms[x as usize].position;
            if a.len() == 2 {
                let side = side_of(&m, a[0], a[1]);
                let dir = unit(sub(pos(&m, a[1]), pos(&m, a[0])));
                let current = norm(sub(pos(&m, a[1]), pos(&m, a[0])));
                let moving = if side.is_empty() { vec![a[1]] } else { side };
                for x in moving {
                    m.atoms[x as usize].position = add(pos(&m, x), scale(dir, value - current));
                }
            } else if a.len() == 3 {
                let side = side_of(&m, a[1], a[2]);
                let (u, v) = (sub(pos(&m, a[0]), pos(&m, a[1])), sub(pos(&m, a[2]), pos(&m, a[1])));
                let mut axis = cross(u, v);
                if norm(axis) < 1e-9 {
                    axis = cross(u, [0.3, 0.5, 0.8]);
                }
                axis = unit(axis);
                let current = dot(unit(u), unit(v)).clamp(-1.0, 1.0).acos() * 180.0 / PI;
                let origin = pos(&m, a[1]);
                let moving = if side.is_empty() { vec![a[2]] } else { side };
                for x in moving {
                    m.atoms[x as usize].position = add(origin, rotate_about(sub(pos(&m, x), origin), axis, (value - current) * PI / 180.0));
                }
            } else {
                set_dihedral(&mut m, a[0], a[1], a[2], a[3], value)?;
            }
            Ok(())
        })();
        if let Err(e) = moved {
            if !relaxed {
                return Err(e);
            }
        }
        let mut x: Vec<Vec3> = m.atoms.iter().map(|at| at.position).collect();
        let energy;
        if relaxed {
            let mut restrained = field.clone();
            let mut r = Restraint { kind: a.len() as i32, target: value, k: 2000.0, ..Default::default() };
            for (q, &atom) in a.iter().enumerate() {
                r.atoms[q] = atom;
            }
            restrained.restraints.push(r);
            minimize(&restrained, &mut x, &[], optimize_steps, 0.02)?;
            let mut terms = EnergyTerms::default();
            restrained.energy(&x, None, Some(&mut terms))?;
            energy = terms.total() - terms.restraint;
            for (i, &xi) in x.iter().enumerate() {
                m.atoms[i].position = xi;
            }
        } else {
            energy = field.energy(&x, None, None)?;
        }
        out.push(ScanPoint { value, energy, positions: x });
    }
    Ok(out)
}
