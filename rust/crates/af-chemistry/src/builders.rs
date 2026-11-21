//! Builders for organic, biological and coordination chemistry (the C++
//! chemistry/Builders): a library of named compounds (as SMILES), peptides
//! from a sequence with backbone dihedrals, polymers from a repeat unit,
//! coordination complexes from a metal, a geometry and ligands, Z-matrices
//! (internal coordinates) and quantum-chemistry input files.

use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fmt::Write as _;

use af_core::crmath;
use af_core::cxx_random::{uniform_real, Mt19937};
use af_core::format::stod;
use af_core::structure::Structure;

use crate::compound_library::COMPOUNDS;
use crate::editing::{dihedral_angle, rmsd, set_dihedral, side_of, substitute};
use crate::embedding::{molecule_from_smiles, EmbedOptions};
use crate::forcefield::{minimize, uff_bond_length, uff_type, ForceFieldOptions, UffForceField};
use crate::molecule::*;

const PI: f64 = 3.14159265358979323846;

fn rotate_about(v: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let (c, s) = if angle == PI { (crmath::cos(PI), crmath::sin(PI)) } else { (angle.cos(), angle.sin()) };
    add(add(scale(v, c), scale(cross(axis, v), s)), scale(axis, dot(axis, v) * (1.0 - c)))
}

fn align_vector(v: Vec3, a: Vec3, b: Vec3) -> Vec3 {
    let axis = cross(a, b);
    let (s, c) = (norm(axis), dot(a, b));
    if s < 1e-9 {
        if c > 0.0 {
            return v;
        }
        let p = unit(cross(a, if a[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] }));
        return rotate_about(v, p, PI);
    }
    rotate_about(v, scale(axis, 1.0 / s), s.atan2(c))
}

/// Lower case without spaces, '-', '_', ',' and '\''.
fn squeeze(s: &str) -> String {
    s.bytes().filter(|c| !matches!(c, b' ' | b'-' | b'_' | b',' | b'\'')).map(|c| c.to_ascii_lowercase() as char).collect()
}

/// C isspace.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r')
}

/// istream >> std::string words.
fn words(s: &str) -> Vec<String> {
    s.as_bytes().split(|&c| is_space(c)).filter(|w| !w.is_empty()).map(|w| String::from_utf8_lossy(w).into_owned()).collect()
}

/// std::stoi with the characters used; None where it throws.
fn stoi_used(s: &str) -> Option<(i32, usize)> {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && is_space(b[i]) {
        i += 1;
    }
    let start = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let digits = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == digits {
        return None;
    }
    let v: i64 = s[start..i].parse().unwrap_or(i64::MAX);
    if v < i64::from(i32::MIN) || v > i64::from(i32::MAX) {
        return None;
    }
    Some((v as i32, i))
}

/// Hydrogens for every atom short of its valence, placed with the coordination
/// they will have (an amide N with its H is planar).
fn complete_hydrogens(m: &mut Molecule) {
    for i in 0..m.atoms.len() {
        let (z, charge, dummy) = (m.atoms[i].z, m.atoms[i].charge, m.atoms[i].dummy);
        if z <= 1 || is_metal(z) || dummy {
            continue;
        }
        let sum = m.valence(i as i32);
        let target = default_valence(z, charge, sum);
        if target > sum {
            m.atoms[i].hydrogens = target - sum;
        }
    }
    make_hydrogens_explicit(m, true);
}

/// NeRF: the position d with |cd| = bond, angle bcd and dihedral abcd (degrees).
pub fn place_atom(a: Vec3, b: Vec3, c: Vec3, bond: f64, angle: f64, torsion: f64) -> Vec3 {
    let (t, p) = (angle * PI / 180.0, torsion * PI / 180.0);
    let bc = unit(sub(c, b));
    let mut n = cross(sub(b, a), bc);
    if norm(n) < 1e-9 {
        n = cross(bc, if bc[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] });
    }
    n = unit(n);
    let m = cross(n, bc);
    let (dx, dy, dz) = (-bond * t.cos(), bond * t.sin() * p.cos(), bond * t.sin() * p.sin());
    add(c, add(add(scale(bc, dx), scale(m, dy)), scale(n, dz)))
}

// ------------------------------------------------------------------ compound library

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NamedCompound {
    pub name: &'static str,
    pub category: &'static str,
    pub smiles: &'static str,
}

pub fn compound_library() -> &'static [NamedCompound] {
    COMPOUNDS
}

/// SMILES for a name (case-insensitive; spaces, '-' and '_' ignored), or "".
pub fn compound_smiles(name: &str) -> String {
    let key = squeeze(name);
    if let Some(c) = COMPOUNDS.iter().find(|c| squeeze(c.name) == key) {
        return c.smiles.to_string();
    }
    const ALIASES: &[(&str, &str)] = &[
        ("isopropanol", "2-propanol"), ("ipa", "2-propanol"), ("meoh", "methanol"), ("etoh", "ethanol"), ("dcm", "dichloromethane"),
        ("methylenechloride", "dichloromethane"), ("ether", "diethyl ether"), ("etoac", "ethyl acetate"), ("acoh", "acetic acid"),
        ("mecn", "acetonitrile"), ("acn", "acetonitrile"), ("dioxane", "1,4-dioxane"), ("ec", "ethylene carbonate"),
        ("pc", "propylene carbonate"), ("dmc", "dimethyl carbonate"), ("dec", "diethyl carbonate"), ("acetaminophen", "paracetamol"),
        ("tylenol", "paracetamol"), ("gly", "glycine"), ("ala", "alanine"), ("val", "valine"), ("leu", "leucine"), ("ile", "isoleucine"),
        ("ser", "serine"), ("thr", "threonine"), ("cys", "cysteine"), ("met", "methionine"), ("pro", "proline"), ("phe", "phenylalanine"),
        ("tyr", "tyrosine"), ("trp", "tryptophan"), ("his", "histidine"), ("asp", "aspartic acid"), ("glu", "glutamic acid"),
        ("asn", "asparagine"), ("gln", "glutamine"), ("lys", "lysine"), ("arg", "arginine"), ("bipy", "2,2'-bipyridine"),
        ("bpy", "2,2'-bipyridine"), ("phen", "1,10-phenanthroline"), ("pph3", "triphenylphosphine"), ("bf4", "tetrafluoroborate"),
        ("pf6", "hexafluorophosphate"), ("otf", "triflate"), ("ntf2", "TFSI"), ("co2", "carbon dioxide"), ("co", "carbon monoxide"),
        ("h2o", "water"), ("nh3", "ammonia"), ("ch4", "methane"), ("h2", "hydrogen"), ("n2", "nitrogen"), ("h2s", "hydrogen sulfide"),
        ("so2", "sulfur dioxide"), ("hcn", "hydrogen cyanide"), ("c6h6", "benzene"), ("tea", "trimethylamine"),
    ];
    match ALIASES.iter().find(|(a, _)| *a == key) {
        Some((_, target)) => compound_smiles(target),
        None => String::new(),
    }
}

// ------------------------------------------------------------------ peptides

struct ZAtom {
    name: &'static str,
    z: i32,
    /// Bonded to, angle with, dihedral with.
    a: &'static str,
    b: &'static str,
    c: &'static str,
    bond: f64,
    angle: f64,
    torsion: f64,
}

struct Residue {
    one: u8,
    three: &'static str,
    side: Vec<ZAtom>,
    /// Ring closures and double bonds (orders).
    extra_bonds: Vec<(&'static str, &'static str, i32)>,
    /// Atom -> charge when side chains are ionised.
    ionized: Vec<(&'static str, i32)>,
}

#[allow(clippy::too_many_arguments)]
fn za(name: &'static str, z: i32, a: &'static str, b: &'static str, c: &'static str, bond: f64, angle: f64, torsion: f64) -> ZAtom {
    ZAtom { name, z, a, b, c, bond, angle, torsion }
}

fn residues() -> &'static [Residue] {
    static TABLE: std::sync::OnceLock<Vec<Residue>> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let r = |one: u8, three, side, extra_bonds, ionized| Residue { one, three, side, extra_bonds, ionized };
        vec![
            r(b'G', "GLY", vec![], vec![], vec![]),
            r(b'A', "ALA", vec![], vec![], vec![]),
            r(b'S', "SER", vec![za("OG", 8, "CB", "CA", "N", 1.417, 110.8, -63.3)], vec![], vec![]),
            r(b'C', "CYS", vec![za("SG", 16, "CB", "CA", "N", 1.808, 113.8, -62.2)], vec![], vec![]),
            r(b'V', "VAL", vec![za("CG1", 6, "CB", "CA", "N", 1.527, 110.7, 177.2), za("CG2", 6, "CB", "CA", "N", 1.527, 110.4, -63.3)], vec![], vec![]),
            r(b'T', "THR", vec![za("OG1", 8, "CB", "CA", "N", 1.43, 109.2, 60.0), za("CG2", 6, "CB", "CA", "N", 1.53, 111.1, -60.0)], vec![], vec![]),
            r(
                b'L',
                "LEU",
                vec![za("CG", 6, "CB", "CA", "N", 1.53, 116.1, -60.1), za("CD1", 6, "CG", "CB", "CA", 1.524, 110.3, 174.9), za("CD2", 6, "CG", "CB", "CA", 1.525, 110.6, 66.7)],
                vec![],
                vec![],
            ),
            r(
                b'I',
                "ILE",
                vec![za("CG1", 6, "CB", "CA", "N", 1.527, 110.7, 59.7), za("CG2", 6, "CB", "CA", "N", 1.527, 110.4, -61.6), za("CD1", 6, "CG1", "CB", "CA", 1.52, 113.97, 169.8)],
                vec![],
                vec![],
            ),
            r(
                b'M',
                "MET",
                vec![za("CG", 6, "CB", "CA", "N", 1.52, 113.68, -64.4), za("SD", 16, "CG", "CB", "CA", 1.81, 112.69, -179.6), za("CE", 6, "SD", "CG", "CB", 1.79, 100.61, 70.0)],
                vec![],
                vec![],
            ),
            r(b'P', "PRO", vec![za("CG", 6, "CB", "CA", "N", 1.495, 104.21, 29.6), za("CD", 6, "CG", "CB", "CA", 1.507, 105.03, -34.8)], vec![("CD", "N", 1)], vec![]),
            r(
                b'F',
                "PHE",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.50, 113.85, -64.7),
                    za("CD1", 6, "CG", "CB", "CA", 1.39, 120.0, 93.3),
                    za("CD2", 6, "CG", "CB", "CA", 1.39, 120.0, -86.7),
                    za("CE1", 6, "CD1", "CG", "CB", 1.39, 120.0, 180.0),
                    za("CE2", 6, "CD2", "CG", "CB", 1.39, 120.0, 180.0),
                    za("CZ", 6, "CE1", "CD1", "CG", 1.39, 120.0, 0.0),
                ],
                vec![("CZ", "CE2", 1), ("CG", "CD1", 2), ("CE1", "CZ", 2), ("CD2", "CE2", 2)],
                vec![],
            ),
            r(
                b'Y',
                "TYR",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.51, 113.8, -64.3),
                    za("CD1", 6, "CG", "CB", "CA", 1.39, 120.98, 93.1),
                    za("CD2", 6, "CG", "CB", "CA", 1.39, 120.82, -86.9),
                    za("CE1", 6, "CD1", "CG", "CB", 1.39, 120.0, 180.0),
                    za("CE2", 6, "CD2", "CG", "CB", 1.39, 120.0, 180.0),
                    za("CZ", 6, "CE1", "CD1", "CG", 1.39, 120.0, 0.0),
                    za("OH", 8, "CZ", "CE1", "CD1", 1.39, 119.78, 180.0),
                ],
                vec![("CZ", "CE2", 1), ("CG", "CD1", 2), ("CE1", "CZ", 2), ("CD2", "CE2", 2)],
                vec![],
            ),
            r(
                b'W',
                "TRP",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.50, 114.1, -66.4),
                    za("CD1", 6, "CG", "CB", "CA", 1.37, 127.07, 96.3),
                    za("CD2", 6, "CG", "CB", "CA", 1.43, 126.66, -83.7),
                    za("NE1", 7, "CD1", "CG", "CB", 1.38, 108.5, 180.0),
                    za("CE2", 6, "CD2", "CG", "CB", 1.40, 108.5, 180.0),
                    za("CE3", 6, "CD2", "CG", "CB", 1.40, 133.83, 0.0),
                    za("CZ2", 6, "CE2", "CD2", "CG", 1.40, 120.0, 180.0),
                    za("CZ3", 6, "CE3", "CD2", "CG", 1.39, 120.0, 180.0),
                    za("CH2", 6, "CZ2", "CE2", "CD2", 1.39, 120.0, 0.0),
                ],
                vec![("NE1", "CE2", 1), ("CH2", "CZ3", 1), ("CG", "CD1", 2), ("CD2", "CE3", 2), ("CE2", "CZ2", 2), ("CZ3", "CH2", 2)],
                vec![],
            ),
            r(
                b'H',
                "HIS",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.50, 113.74, -63.2),
                    za("ND1", 7, "CG", "CB", "CA", 1.38, 122.85, -75.7),
                    za("CD2", 6, "CG", "CB", "CA", 1.36, 130.61, 104.3),
                    za("CE1", 6, "ND1", "CG", "CB", 1.32, 108.5, 180.0),
                    za("NE2", 7, "CD2", "CG", "CB", 1.35, 108.5, 180.0),
                ],
                vec![("CE1", "NE2", 1), ("CG", "CD2", 2), ("ND1", "CE1", 2)],
                vec![],
            ),
            r(
                b'D',
                "ASP",
                vec![za("CG", 6, "CB", "CA", "N", 1.52, 113.06, -66.4), za("OD1", 8, "CG", "CB", "CA", 1.25, 119.22, -46.7), za("OD2", 8, "CG", "CB", "CA", 1.25, 118.22, 133.3)],
                vec![("CG", "OD1", 2)],
                vec![("OD2", -1)],
            ),
            r(
                b'N',
                "ASN",
                vec![za("CG", 6, "CB", "CA", "N", 1.52, 112.62, -65.5), za("OD1", 8, "CG", "CB", "CA", 1.23, 120.85, -58.3), za("ND2", 7, "CG", "CB", "CA", 1.33, 116.48, 121.7)],
                vec![("CG", "OD1", 2)],
                vec![],
            ),
            r(
                b'E',
                "GLU",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.52, 113.82, -63.4),
                    za("CD", 6, "CG", "CB", "CA", 1.52, 113.31, -179.8),
                    za("OE1", 8, "CD", "CG", "CB", 1.25, 119.02, -6.2),
                    za("OE2", 8, "CD", "CG", "CB", 1.25, 118.08, 173.8),
                ],
                vec![("CD", "OE1", 2)],
                vec![("OE2", -1)],
            ),
            r(
                b'Q',
                "GLN",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.52, 113.75, -60.2),
                    za("CD", 6, "CG", "CB", "CA", 1.52, 112.78, -69.6),
                    za("OE1", 8, "CD", "CG", "CB", 1.24, 120.86, -50.5),
                    za("NE2", 7, "CD", "CG", "CB", 1.33, 116.5, 129.5),
                ],
                vec![("CD", "OE1", 2)],
                vec![],
            ),
            r(
                b'K',
                "LYS",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.52, 113.83, -64.5),
                    za("CD", 6, "CG", "CB", "CA", 1.52, 111.79, -178.1),
                    za("CE", 6, "CD", "CG", "CB", 1.52, 111.68, -179.6),
                    za("NZ", 7, "CE", "CD", "CG", 1.49, 111.36, 179.6),
                ],
                vec![],
                vec![("NZ", 1)],
            ),
            r(
                b'R',
                "ARG",
                vec![
                    za("CG", 6, "CB", "CA", "N", 1.52, 113.83, -65.2),
                    za("CD", 6, "CG", "CB", "CA", 1.52, 111.79, -179.2),
                    za("NE", 7, "CD", "CG", "CB", 1.46, 111.68, -179.3),
                    za("CZ", 6, "NE", "CD", "CG", 1.33, 124.79, -178.7),
                    za("NH1", 7, "CZ", "NE", "CD", 1.33, 120.64, 0.0),
                    za("NH2", 7, "CZ", "NE", "CD", 1.33, 119.63, 180.0),
                ],
                vec![("CZ", "NH2", 2)],
                vec![("NH2", 1)],
            ),
        ]
    })
}

fn residue_for(code: &str) -> Result<&'static Residue, String> {
    let mut c = code.to_ascii_uppercase();
    if matches!(c.as_str(), "HID" | "HIE" | "HIP" | "HSD" | "HSE") {
        c = "HIS".into();
    }
    if matches!(c.as_str(), "CYX" | "CYM") {
        c = "CYS".into();
    }
    if c == "ASH" {
        c = "ASP".into();
    }
    if c == "GLH" {
        c = "GLU".into();
    }
    if c == "LYN" {
        c = "LYS".into();
    }
    residues()
        .iter()
        .find(|r| (c.len() == 1 && r.one == c.as_bytes()[0]) || c == r.three)
        .ok_or_else(|| format!("unknown amino acid '{code}' (use one-letter codes such as ACDEFGHIKLMNPQRSTVWY or three-letter codes)"))
}

#[derive(Clone, Debug, PartialEq)]
pub struct PeptideOptions {
    /// One-letter (ACDE...) or three-letter codes (Ala-Gly-...).
    pub sequence: String,
    /// alpha, 3-10, pi, beta, parallel-beta, ppii, extended, custom.
    pub conformation: String,
    /// Custom backbone dihedrals.
    pub phi: f64,
    pub psi: f64,
    pub omega: f64,
    /// charged (NH3+, COO-), neutral (NH2, COOH), capped (ACE, NME).
    pub termini: String,
    /// Asp, Glu deprotonated; Lys, Arg protonated.
    pub ionize_side_chains: bool,
    /// D-amino acids (mirror image).
    pub d_amino: bool,
}

impl Default for PeptideOptions {
    fn default() -> Self {
        Self {
            sequence: String::new(),
            conformation: "alpha".into(),
            phi: -57.0,
            psi: -47.0,
            omega: 180.0,
            termini: "charged".into(),
            ionize_side_chains: true,
            d_amino: false,
        }
    }
}

/// Heavy atoms from standard geometry (Engh and Huber) by the NeRF method,
/// hydrogens from the valences; also the residue label of every atom.
pub fn build_peptide(o: &PeptideOptions) -> Result<(Molecule, Vec<String>), String> {
    // Sequence: three-letter codes with separators, else one letter per residue.
    let mut sequence: Vec<&Residue> = Vec::new();
    if o.sequence.bytes().any(|c| matches!(c, b'-' | b' ' | b',')) {
        let s: String = o.sequence.chars().map(|c| if c == '-' || c == ',' { ' ' } else { c }).collect();
        for token in words(&s) {
            sequence.push(residue_for(&token)?);
        }
    } else {
        for &c in o.sequence.as_bytes() {
            if !is_space(c) {
                sequence.push(residue_for(&String::from_utf8_lossy(&[c]))?);
            }
        }
    }
    if sequence.is_empty() {
        return Err("give an amino-acid sequence, e.g. ACDKLG or Ala-Gly-Ser".into());
    }
    if sequence.len() > 500 {
        return Err("at most 500 residues".into());
    }
    let (mut phi, mut psi, mut omega) = (o.phi, o.psi, o.omega);
    let conf = squeeze(&o.conformation);
    match conf.as_str() {
        "alpha" | "alphahelix" => (phi, psi, omega) = (-57.0, -47.0, 180.0),
        "310" | "310helix" => (phi, psi, omega) = (-49.0, -26.0, 180.0),
        "pi" | "pihelix" => (phi, psi, omega) = (-57.0, -70.0, 180.0),
        "beta" | "antiparallelbeta" | "betasheet" => (phi, psi, omega) = (-139.0, 135.0, 180.0),
        "parallelbeta" => (phi, psi, omega) = (-119.0, 113.0, 180.0),
        "ppii" | "polyproline" => (phi, psi, omega) = (-75.0, 145.0, 180.0),
        "extended" => (phi, psi, omega) = (180.0, 180.0, 180.0),
        "custom" => {}
        _ => return Err("conformation is alpha, 3-10, pi, beta, parallel-beta, ppii, extended or custom".into()),
    }
    let termini = squeeze(&o.termini);
    if termini != "charged" && termini != "neutral" && termini != "capped" {
        return Err("termini are charged, neutral or capped".into());
    }

    let mut m = Molecule::default();
    let mut owner: Vec<String> = Vec::new();
    let add_atom = |m: &mut Molecule, owner: &mut Vec<String>, z: i32, p: Vec3, label: &str| {
        owner.push(label.to_string());
        m.add_atom(Atom { z, position: p, ..Default::default() })
    };
    let (mut n_idx, mut ca_idx, mut c_idx): (Vec<i32>, Vec<i32>, Vec<i32>) = (Vec::new(), Vec::new(), Vec::new());
    let count = sequence.len();
    let pos = |m: &Molecule, i: i32| m.atoms[i as usize].position;
    for (r, res) in sequence.iter().enumerate() {
        let label = format!("{}{}", res.three, r + 1);
        let phi_r = if res.one == b'P' && conf != "custom" { -63.0 } else { phi };
        let (pn, pca, pc);
        if r == 0 {
            pn = [0.0, 0.0, 0.0];
            pca = [1.458, 0.0, 0.0];
            let t = (180.0 - 111.2) * PI / 180.0;
            pc = [1.458 + 1.525 * crmath::cos(t), 1.525 * crmath::sin(t), 0.0];
        } else {
            let (n0, ca0, c0) = (pos(&m, *n_idx.last().unwrap()), pos(&m, *ca_idx.last().unwrap()), pos(&m, *c_idx.last().unwrap()));
            pn = place_atom(n0, ca0, c0, 1.329, 116.2, psi);
            pca = place_atom(ca0, c0, pn, 1.458, 121.7, omega);
            pc = place_atom(c0, pn, pca, 1.525, 111.2, phi_r);
        }
        let in_ = add_atom(&mut m, &mut owner, 7, pn, &label);
        let ica = add_atom(&mut m, &mut owner, 6, pca, &label);
        let ic = add_atom(&mut m, &mut owner, 6, pc, &label);
        m.add_bond(in_, ica, 1)?;
        m.add_bond(ica, ic, 1)?;
        if r > 0 {
            m.add_bond(*c_idx.last().unwrap(), in_, 1)?;
        }
        n_idx.push(in_);
        ca_idx.push(ica);
        c_idx.push(ic);
        // Carbonyl O anti to the next N (psi + 180).
        let io = add_atom(&mut m, &mut owner, 8, place_atom(pn, pca, pc, 1.231, 120.5, psi + 180.0), &label);
        m.add_bond(ic, io, 2)?;
        let mut named: HashMap<&str, i32> = HashMap::from([("N", in_), ("CA", ica), ("C", ic), ("O", io)]);
        if res.one != b'G' {
            let icb = add_atom(&mut m, &mut owner, 6, place_atom(pc, pn, pca, 1.53, 110.5, -122.6), &label);
            m.add_bond(ica, icb, 1)?;
            named.insert("CB", icb);
        }
        for z in &res.side {
            let p = place_atom(pos(&m, named[z.c]), pos(&m, named[z.b]), pos(&m, named[z.a]), z.bond, z.angle, z.torsion);
            let idx = add_atom(&mut m, &mut owner, z.z, p, &label);
            named.insert(z.name, idx);
            m.add_bond(named[z.a], idx, 1)?;
        }
        for &(na, nb, order) in &res.extra_bonds {
            let (a, b) = (named[na], named[nb]);
            let k = m.bond_between(a, b);
            if k >= 0 {
                m.bonds[k as usize].order = order;
            } else {
                m.add_bond(a, b, order)?;
            }
        }
        if o.ionize_side_chains {
            for &(atom, charge) in &res.ionized {
                m.atoms[named[atom] as usize].charge = charge;
            }
        }
    }
    // Termini.
    if termini == "charged" {
        m.atoms[n_idx[0] as usize].charge = 1;
    } else if termini == "capped" {
        // Acetyl before the first N (C0 at phi, O0 cis to CA1, CH3 trans), N-methyl after the last C.
        let (n1, ca1, c1) = (pos(&m, n_idx[0]), pos(&m, ca_idx[0]), pos(&m, c_idx[0]));
        let pc0 = place_atom(c1, ca1, n1, 1.329, 121.7, phi);
        let c0 = add_atom(&mut m, &mut owner, 6, pc0, "ACE0");
        let o0 = add_atom(&mut m, &mut owner, 8, place_atom(ca1, n1, pc0, 1.231, 122.7, 0.0), "ACE0");
        let ch3 = add_atom(&mut m, &mut owner, 6, place_atom(ca1, n1, pc0, 1.52, 116.2, 180.0), "ACE0");
        m.add_bond(c0, n_idx[0], 1)?;
        m.add_bond(c0, o0, 2)?;
        m.add_bond(c0, ch3, 1)?;
    }
    let (last_n, last_ca, last_c) = (*n_idx.last().unwrap(), *ca_idx.last().unwrap(), *c_idx.last().unwrap());
    let (nl, cal, cl) = (pos(&m, last_n), pos(&m, last_ca), pos(&m, last_c));
    let last_label = owner[last_c as usize].clone();
    if termini == "capped" {
        let pn = place_atom(nl, cal, cl, 1.329, 116.2, psi);
        let label = format!("NME{}", count + 1);
        let nn = add_atom(&mut m, &mut owner, 7, pn, &label);
        let cm = add_atom(&mut m, &mut owner, 6, place_atom(cal, cl, pn, 1.458, 121.7, omega), &label);
        m.add_bond(last_c, nn, 1)?;
        m.add_bond(nn, cm, 1)?;
    } else {
        let oxt = add_atom(&mut m, &mut owner, 8, place_atom(nl, cal, cl, 1.25, 117.0, psi), &last_label);
        m.add_bond(last_c, oxt, 1)?;
        if termini == "charged" {
            m.atoms[oxt as usize].charge = -1;
        }
    }
    let heavy = m.atoms.len();
    complete_hydrogens(&mut m);
    for i in heavy..m.atoms.len() {
        let host = m.neighbours(i as i32)[0];
        owner.push(owner[host as usize].clone());
    }
    if o.d_amino {
        for a in m.atoms.iter_mut() {
            a.position[0] = -a.position[0];
        }
    }
    perceive_aromaticity(&mut m);
    Ok((m, owner))
}

// ------------------------------------------------------------------ polymers

pub fn polymer_library() -> &'static [(&'static str, &'static str)] {
    &[
        ("PE", "*CC*"), ("PP", "*CC(C)*"), ("PS", "*CC(c1ccccc1)*"), ("PVC", "*CC(Cl)*"), ("PTFE", "*C(F)(F)C(F)(F)*"),
        ("PVDF", "*CC(F)(F)*"), ("PMMA", "*CC(C)(C(=O)OC)*"), ("PEO", "*CCO*"), ("PEG", "*CCO*"), ("PPO", "*CC(C)O*"),
        ("POM", "*CO*"), ("PVA", "*CC(O)*"), ("PVAc", "*CC(OC(C)=O)*"), ("PAN", "*CC(C#N)*"), ("PAA", "*CC(C(=O)O)*"),
        ("PIB", "*CC(C)(C)*"), ("PB", "*CC=CC*"), ("PI", "*CC(C)=CC*"), ("PDMS", "*[Si](C)(C)O*"), ("PLA", "*OC(C)C(=O)*"),
        ("PGA", "*OCC(=O)*"), ("PCL", "*OCCCCCC(=O)*"), ("PET", "*OCCOC(=O)c1ccc(cc1)C(=O)*"), ("PBT", "*OCCCCOC(=O)c1ccc(cc1)C(=O)*"),
        ("nylon-6", "*NCCCCCC(=O)*"), ("nylon-6,6", "*NCCCCCCNC(=O)CCCCC(=O)*"), ("PC", "*Oc1ccc(cc1)C(C)(C)c1ccc(cc1)OC(=O)*"),
        ("PA", "*C=C*"), ("PPV", "*c1ccc(cc1)C=C*"), ("PPP", "*c1ccc(cc1)*"), ("PT", "*c1ccc(s1)*"), ("P3HT", "*c1cc(CCCCCC)c(s1)*"),
        ("PEDOT", "*c1c2OCCOc2c(s1)*"), ("PPy", "*c1ccc([nH]1)*"), ("PANI", "*c1ccc(cc1)N*"), ("PVP", "*CC(N1CCCC1=O)*"),
        ("PEI", "*CCN*"), ("PVK", "*CC(n1c2ccccc2c2ccccc12)*"), ("PEEK", "*Oc1ccc(cc1)Oc1ccc(cc1)C(=O)c1ccc(cc1)*"),
        ("Kevlar", "*NC(=O)c1ccc(cc1)C(=O)Nc1ccc(cc1)*"), ("cellulose", "*O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]1CO)*"),
    ]
}

#[derive(Clone, Debug, PartialEq)]
pub struct PolymerOptions {
    /// SMILES with two '*' (head, tail) or a library name (PE, PS, ...).
    pub monomer: String,
    pub units: i32,
    /// isotactic, syndiotactic, atactic.
    pub tacticity: String,
    /// trans (zig-zag) or random (gauche and trans).
    pub conformation: String,
    /// End groups: H, CH3, or a group name / SMILES with '*'.
    pub head_group: String,
    pub tail_group: String,
    /// Infinite chain along c (no end groups).
    pub periodic: bool,
    pub seed: u32,
}

impl Default for PolymerOptions {
    fn default() -> Self {
        Self {
            monomer: "*CC*".into(),
            units: 10,
            tacticity: "isotactic".into(),
            conformation: "trans".into(),
            head_group: "H".into(),
            tail_group: "H".into(),
            periodic: false,
            seed: 1,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PolymerResult {
    pub molecule: Molecule,
    /// Periodic: c along the chain.
    pub cell: [[f64; 3]; 3],
    /// Periodic: RMS mismatch of the closing unit (A).
    pub misfit: f64,
    pub monomer_smiles: String,
}

/// The growing chain of build_polymer.
struct Chain<'a> {
    p: &'a Molecule,
    last_tail: i32,
    last_tail_dummy: i32,
    first_head_dummy: i32,
    backbone: &'a [i32],
}

struct RepeatUnit<'a> {
    head: usize,
    d_head: usize,
    d_tail: usize,
    path: &'a [i32],
    link: f64,
}

/// Places a unit (copy u) so its head bonds to the current last tail.
fn place_unit(chain: &Chain, ru: &RepeatUnit, u: &Molecule, torsion: f64) -> Vec<Vec3> {
    if chain.last_tail < 0 {
        return u.atoms.iter().map(|a| a.position).collect();
    }
    let p = chain.p;
    let t = p.atoms[chain.last_tail as usize].position;
    let dir = unit(sub(p.atoms[chain.last_tail_dummy as usize].position, t));
    let v = unit(sub(u.atoms[ru.head].position, u.atoms[ru.d_head].position));
    let h_new = add(t, scale(dir, ru.link));
    let mut positions: Vec<Vec3> = u.atoms.iter().map(|a| add(h_new, align_vector(sub(a.position, u.atoms[ru.head].position), v, dir))).collect();
    // Backbone torsion X - T - H' - Y.
    let x = if chain.backbone.len() >= 2 { chain.backbone[chain.backbone.len() - 2] } else { chain.first_head_dummy };
    let px = p.atoms[x as usize].position;
    let py = if ru.path.len() > 1 { positions[ru.path[1] as usize] } else { positions[ru.d_tail] };
    let current = dihedral_angle(px, t, h_new, py);
    let turn = (torsion - current) * PI / 180.0;
    for q in positions.iter_mut() {
        *q = add(h_new, rotate_about(sub(*q, h_new), dir, turn));
    }
    positions
}

pub fn build_polymer(o: &PolymerOptions) -> Result<PolymerResult, String> {
    if o.units < 1 {
        return Err("a polymer needs at least one repeat unit".into());
    }
    if o.units > 2000 {
        return Err("at most 2000 repeat units".into());
    }
    let mut smiles = o.monomer.clone();
    for (name, s) in polymer_library() {
        if squeeze(name) == squeeze(&o.monomer) {
            smiles = s.to_string();
        }
    }
    if smiles.bytes().filter(|&c| c == b'*').count() != 2 {
        return Err("the repeat unit needs two '*' (head and tail), e.g. *CC(C)* for polypropylene; or a name: PE, PP, PS, ...".into());
    }
    let (tact, conf) = (squeeze(&o.tacticity), squeeze(&o.conformation));
    if tact != "isotactic" && tact != "syndiotactic" && tact != "atactic" {
        return Err("tacticity is isotactic, syndiotactic or atactic".into());
    }
    if conf != "trans" && conf != "random" {
        return Err("conformation is trans or random".into());
    }
    let eo = EmbedOptions { attempts: 3, seed: o.seed, ..Default::default() };
    let (mut repeat, _) = molecule_from_smiles(&smiles, &eo).map_err(|e| e.message)?;
    let (mut d_head, mut d_tail) = (-1i32, -1i32);
    for (i, a) in repeat.atoms.iter().enumerate() {
        if a.dummy {
            if d_head < 0 {
                d_head = i as i32;
            } else {
                d_tail = i as i32;
            }
        }
    }
    let head = repeat.neighbours(d_head)[0];
    let tail = repeat.neighbours(d_tail)[0];
    // Backbone path head -> tail.
    let mut path: Vec<i32> = Vec::new();
    {
        let mut previous = vec![-2i32; repeat.atoms.len()];
        let mut q = VecDeque::from([head]);
        previous[head as usize] = -1;
        while let Some(a) = q.pop_front() {
            if a == tail {
                break;
            }
            for b in repeat.neighbours(a) {
                if previous[b as usize] == -2 && !repeat.atoms[b as usize].dummy {
                    previous[b as usize] = a;
                    q.push_back(b);
                }
            }
        }
        let mut a = tail;
        while a != -1 {
            if a < 0 {
                return Err("the repeat unit's head and tail are not connected".into());
            }
            path.push(a);
            a = previous[a as usize];
        }
        path.reverse();
    }
    // Backbone torsions inside the repeat unit: trans (zig-zag), or random
    // gauche/trans for a random coil (bonds in rings and double bonds stay).
    let mut torsion_random = Mt19937::new(o.seed.wrapping_mul(7919).wrapping_add(1));
    {
        let mut ext = vec![d_head];
        ext.extend_from_slice(&path);
        ext.push(d_tail);
        let mut q = 0;
        while q + 3 < ext.len() {
            let k = repeat.bond_between(ext[q + 1], ext[q + 2]);
            if k >= 0 && repeat.bonds[k as usize].order == 1 && !repeat.bonds[k as usize].aromatic {
                let mut target = 180.0;
                if conf == "random" {
                    let r = uniform_real(&mut torsion_random, 0.0, 1.0);
                    target = if r < 0.5 {
                        180.0
                    } else if r < 0.75 {
                        60.0
                    } else {
                        -60.0
                    };
                }
                let _ = set_dihedral(&mut repeat, ext[q], ext[q + 1], ext[q + 2], ext[q + 3], target); // in a ring: kept
            }
            q += 1;
        }
    }
    // GCC evaluates the arguments right to left: the head's type first.
    let head_type = uff_type(&repeat, head)?;
    let tail_type = uff_type(&repeat, tail)?;
    let link = uff_bond_length(&tail_type, &head_type, 1.0)?;
    let mut mirror = repeat.clone();
    for a in mirror.atoms.iter_mut() {
        a.position[0] = -a.position[0];
    }
    let mut random = Mt19937::new(o.seed);
    let mut result = PolymerResult { monomer_smiles: smiles.clone(), ..Default::default() };
    let mut p = Molecule::default();
    let mut remove: Vec<u8> = Vec::new();
    let mut backbone: Vec<i32> = Vec::new();
    let (mut last_tail, mut last_tail_dummy, mut first_head, mut first_head_dummy) = (-1i32, -1i32, -1i32, -1i32);
    let ru = RepeatUnit { head: head as usize, d_head: d_head as usize, d_tail: d_tail as usize, path: &path, link };
    for k in 0..o.units {
        let use_mirror = match tact.as_str() {
            "syndiotactic" => k % 2 == 1,
            "atactic" => random.next_u32() % 2 == 1,
            _ => false,
        };
        let u = if use_mirror { &mirror } else { &repeat };
        let mut options = vec![180.0];
        if conf == "random" && k > 0 {
            let r = uniform_real(&mut random, 0.0, 1.0);
            options = if r < 0.5 {
                vec![180.0, 60.0, -60.0]
            } else if r < 0.75 {
                vec![60.0, 180.0, -60.0]
            } else {
                vec![-60.0, 180.0, 60.0]
            };
        }
        let mut best = 1i32 << 30;
        let mut best_positions: Vec<Vec3> = Vec::new();
        {
            let chain = Chain { p: &p, last_tail, last_tail_dummy, first_head_dummy, backbone: &backbone };
            for &torsion in &options {
                let positions = place_unit(&chain, &ru, u, torsion);
                let c = if k == 0 {
                    0
                } else {
                    // Clashes of the unit's heavy atoms with the chain.
                    let mut bad = 0;
                    for (i, ua) in u.atoms.iter().enumerate() {
                        if ua.z == 1 || ua.dummy || i == ru.head {
                            continue;
                        }
                        for (j, pa) in p.atoms.iter().enumerate() {
                            if pa.z == 1 || remove[j] != 0 || j as i32 == last_tail {
                                continue;
                            }
                            if norm(sub(positions[i], pa.position)) < 2.6 {
                                bad += 1;
                            }
                        }
                    }
                    bad
                };
                if c < best {
                    best = c;
                    best_positions = positions;
                }
                if c == 0 {
                    break;
                }
            }
        }
        let mut map = vec![0i32; u.atoms.len()];
        for (i, ua) in u.atoms.iter().enumerate() {
            let mut a = ua.clone();
            a.position = best_positions[i];
            map[i] = p.add_atom(a);
            remove.push(0);
        }
        for b in &u.bonds {
            let e = p.add_bond(map[b.a as usize], map[b.b as usize], b.order)?;
            p.bonds[e as usize].aromatic = b.aromatic;
        }
        if last_tail >= 0 {
            p.add_bond(last_tail, map[head as usize], 1)?;
            remove[last_tail_dummy as usize] = 1;
            remove[map[d_head as usize] as usize] = 1;
        } else {
            first_head = map[head as usize];
            first_head_dummy = map[d_head as usize];
        }
        for &a in &path {
            backbone.push(map[a as usize]);
        }
        last_tail = map[tail as usize];
        last_tail_dummy = map[d_tail as usize];
    }
    if o.periodic {
        // Where the next unit's head would go: the chain repeats with that period.
        if tact == "syndiotactic" && o.units % 2 == 1 {
            return Err("a periodic syndiotactic chain needs an even number of units".into());
        }
        let ghost = place_unit(&Chain { p: &p, last_tail, last_tail_dummy, first_head_dummy, backbone: &backbone }, &ru, &repeat, 180.0);
        let period = sub(ghost[head as usize], p.atoms[first_head as usize].position);
        // Misfit: the ghost unit against the first unit shifted by the period.
        let mut sum = 0.0;
        let mut n = 0;
        for (i, a) in repeat.atoms.iter().enumerate() {
            if a.dummy {
                continue;
            }
            let d = sub(ghost[i], add(p.atoms[i].position, period)); // unit 0 holds atoms 0..n-1
            sum += dot(d, d);
            n += 1;
        }
        result.misfit = if n > 0 { (sum / f64::from(n)).sqrt() } else { 0.0 };
        remove[first_head_dummy as usize] = 1;
        remove[last_tail_dummy as usize] = 1;
        let e = p.add_bond(last_tail, first_head, 1)?;
        p.bonds[e as usize].image = period; // from lastTail: position[firstHead] + period - position[lastTail]
        // Orient the chain along z and build the cell.
        let cz = unit(period);
        let length = norm(period);
        for a in p.atoms.iter_mut() {
            a.position = align_vector(a.position, cz, [0.0, 0.0, 1.0]);
        }
        for b in p.bonds.iter_mut() {
            b.image = align_vector(b.image, cz, [0.0, 0.0, 1.0]);
        }
        let mut rmax = 0.0;
        let mut centre = [0.0; 3];
        let mut heavy = 0;
        for (i, a) in p.atoms.iter().enumerate() {
            if remove[i] == 0 {
                centre = add(centre, a.position);
                heavy += 1;
            }
        }
        centre = scale(centre, 1.0 / f64::from(heavy.max(1)));
        for (i, a) in p.atoms.iter().enumerate() {
            if remove[i] == 0 {
                let d = sub(a.position, centre);
                rmax = fmax(rmax, (d[0] * d[0] + d[1] * d[1]).sqrt());
            }
        }
        let lateral = fmax(10.0, 2.0 * rmax + 10.0);
        for a in p.atoms.iter_mut() {
            a.position = add(sub(a.position, [centre[0], centre[1], 0.0]), [lateral / 2.0, lateral / 2.0, 0.0]);
        }
        result.cell = [[lateral, 0.0, 0.0], [0.0, lateral, 0.0], [0.0, 0.0, length]];
    } else {
        // End groups: hydrogens (the dummies become H at the C-H length) or groups.
        for dummy in [first_head_dummy, last_tail_dummy] {
            let host = p.neighbours(dummy)[0];
            let hp = p.atoms[host as usize].position;
            let hz = p.atoms[host as usize].z;
            let a = &mut p.atoms[dummy as usize];
            a.dummy = false;
            a.z = 1;
            let d = unit(sub(a.position, hp));
            a.position = add(hp, scale(d, hydrogen_bond_length(hz)));
        }
    }
    let keep: Vec<u8> = remove.iter().map(|&r| u8::from(r == 0)).collect();
    let map = p.remove_atoms(&keep);
    if !o.periodic {
        // End groups replace the terminal hydrogens (head first; the tail index
        // moves down by one when the removed head H came before it).
        let (head_h, mut tail_h) = (map[first_head_dummy as usize], map[last_tail_dummy as usize]);
        let group_name = |group: &str| {
            let g = squeeze(group);
            if g.is_empty() || g == "h" || g == "hydrogen" {
                String::new()
            } else if g == "ch3" {
                "methyl".to_string()
            } else {
                group.to_string()
            }
        };
        let (head_name, tail_name) = (group_name(&o.head_group), group_name(&o.tail_group));
        if !head_name.is_empty() {
            substitute(&mut p, head_h, &head_name)?;
            if tail_h > head_h {
                tail_h -= 1;
            }
        }
        if !tail_name.is_empty() {
            substitute(&mut p, tail_h, &tail_name)?;
        }
    }
    perceive_aromaticity(&mut p);
    result.molecule = p;
    Ok(result)
}

// ------------------------------------------------------------------ complexes

#[derive(Clone, Debug, PartialEq)]
pub struct LigandInfo {
    pub name: &'static str,
    pub smiles: &'static str,
    /// Heavy-atom indices in the SMILES.
    pub donors: &'static [i32],
    pub description: &'static str,
}

pub fn ligand_library() -> &'static [LigandInfo] {
    const fn l(name: &'static str, smiles: &'static str, donors: &'static [i32], description: &'static str) -> LigandInfo {
        LigandInfo { name, smiles, donors, description }
    }
    const LIBRARY: &[LigandInfo] = &[
        l("H2O", "O", &[0], "aqua"), l("NH3", "N", &[0], "ammine"), l("CO", "[C-]#[O+]", &[0], "carbonyl (C-bound)"),
        l("CN", "[C-]#N", &[0], "cyanido (C-bound)"), l("NCS", "[N-]=C=S", &[0], "isothiocyanato (N-bound)"),
        l("SCN", "[S-]C#N", &[0], "thiocyanato (S-bound)"), l("NO2", "O=N[O-]", &[1], "nitrito-N"), l("F", "[F-]", &[0], "fluorido"),
        l("Cl", "[Cl-]", &[0], "chlorido"), l("Br", "[Br-]", &[0], "bromido"), l("I", "[I-]", &[0], "iodido"), l("OH", "[OH-]", &[0], "hydroxido"),
        l("O", "[O-2]", &[0], "oxido"), l("H", "[H-]", &[0], "hydrido"), l("MeCN", "N#CC", &[0], "acetonitrile"), l("py", "c1ccncc1", &[3], "pyridine"),
        l("PMe3", "CP(C)C", &[1], "trimethylphosphine"), l("PPh3", "P(c1ccccc1)(c1ccccc1)c1ccccc1", &[0], "triphenylphosphine"),
        l("THF", "C1CCOC1", &[3], "tetrahydrofuran"), l("DMSO", "CS(C)=O", &[3], "dimethyl sulfoxide (O-bound)"),
        l("en", "NCCN", &[0, 3], "ethylenediamine (bidentate)"), l("bipy", "c1ccc(nc1)-c1ccccn1", &[4, 11], "2,2'-bipyridine (bidentate)"),
        l("phen", "c1cnc2c(c1)ccc1cccnc12", &[2, 11], "1,10-phenanthroline (bidentate)"),
        l("acac", "CC(=O)/C=C(\\[O-])C", &[2, 5], "acetylacetonate (bidentate)"), l("ox", "[O-]C(=O)C(=O)[O-]", &[0, 5], "oxalate (bidentate)"),
        l("gly", "NCC(=O)[O-]", &[0, 4], "glycinate (bidentate)"), l("tmeda", "CN(C)CCN(C)C", &[1, 5], "TMEDA (bidentate)"),
        l("dppe", "P(CCP(c1ccccc1)c1ccccc1)(c1ccccc1)c1ccccc1", &[0, 3], "dppe (bidentate)"),
        l("CO3", "O=C([O-])[O-]", &[2, 3], "carbonato (bidentate)"),
    ];
    LIBRARY
}

fn geometry_sites(g: &str) -> Result<Vec<Vec3>, String> {
    let s3 = 3.0f64.sqrt() / 2.0;
    let t = 1.0 / 3.0f64.sqrt();
    Ok(match g {
        "linear" => vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
        "trigonalplanar" => vec![[1.0, 0.0, 0.0], [-0.5, s3, 0.0], [-0.5, -s3, 0.0]],
        "tetrahedral" => vec![[t, t, t], [t, -t, -t], [-t, t, -t], [-t, -t, t]],
        "squareplanar" => vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        "trigonalbipyramidal" => vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [-0.5, s3, 0.0], [-0.5, -s3, 0.0]],
        "squarepyramidal" => vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        "octahedral" => vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
        "pentagonalbipyramidal" => {
            let mut v = vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
            for k in 0..5 {
                let a = 2.0 * PI * f64::from(k) / 5.0;
                v.push([crmath::cos(a), crmath::sin(a), 0.0]);
            }
            v
        }
        _ => {
            return Err("geometry is linear, trigonal-planar, tetrahedral, square-planar, trigonal-bipyramidal, square-pyramidal, octahedral or pentagonal-bipyramidal".into())
        }
    })
}

#[derive(Clone, Debug, PartialEq)]
pub struct ComplexOptions {
    pub metal: String,
    pub oxidation_state: i32,
    /// linear, trigonal-planar, tetrahedral, square-planar, trigonal-bipyramidal,
    /// square-pyramidal, octahedral, pentagonal-bipyramidal.
    pub geometry: String,
    /// Names (or smiles:SMILES:donor[,donor2]) with counts.
    pub ligands: String,
    /// Metal-donor distance (0 = covalent radii).
    pub distance: f64,
    /// high or low.
    pub spin: String,
}

impl Default for ComplexOptions {
    fn default() -> Self {
        Self { metal: "Fe".into(), oxidation_state: 2, geometry: "octahedral".into(), ligands: "NH3 6".into(), distance: 0.0, spin: "high".into() }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ComplexResult {
    pub molecule: Molecule,
    pub charge: i32,
    pub d_electrons: i32,
    pub multiplicity: i32,
}

impl Default for ComplexResult {
    fn default() -> Self {
        Self { molecule: Molecule::default(), charge: 0, d_electrons: -1, multiplicity: 1 }
    }
}

#[derive(Clone)]
struct Wanted {
    smiles: String,
    donors: Vec<i32>,
}

pub fn build_complex(o: &ComplexOptions) -> Result<ComplexResult, String> {
    let zm = atomic_number(&o.metal);
    if zm <= 0 {
        return Err(format!("unknown metal '{}'", o.metal));
    }
    let geometry = squeeze(&o.geometry);
    let sites = geometry_sites(&geometry)?;
    let mut used = vec![false; sites.len()];
    // Ligands with counts.
    let mut wanted: Vec<Wanted> = Vec::new();
    {
        let list = words(&o.ligands);
        let mut i = 0;
        while i < list.len() {
            let word = &list[i];
            let item;
            if let Some(rest) = word.strip_prefix("smiles:") {
                let colon = rest.rfind(':').ok_or("custom ligands are smiles:SMILES:donor (heavy-atom index, e.g. smiles:NCCN:0,3)")?;
                let mut donors = Vec::new();
                for token in words(&rest[colon + 1..].replace(',', " ")) {
                    // istream >> int: stops at the first word that is not a number.
                    match stoi_used(&token) {
                        Some((x, used)) if !token.starts_with(char::is_whitespace) => {
                            donors.push(x);
                            if used != token.len() {
                                break;
                            }
                        }
                        _ => break,
                    }
                }
                item = Wanted { smiles: rest[..colon].to_string(), donors };
            } else {
                let mut found = None;
                for l in ligand_library() {
                    if squeeze(l.name) == squeeze(word) {
                        found = Some(Wanted { smiles: l.smiles.to_string(), donors: l.donors.to_vec() });
                    }
                }
                item = found.ok_or_else(|| {
                    let names: Vec<&str> = ligand_library().iter().map(|l| l.name).collect();
                    format!("unknown ligand '{word}'; use one of {} or smiles:SMILES:donor", names.join(", "))
                })?;
            }
            let mut count = 1;
            if i + 1 < list.len() && list[i + 1].bytes().all(|c| c.is_ascii_digit()) {
                i += 1;
                count = stoi_used(&list[i]).map(|(v, _)| v).ok_or("stoi")?;
            }
            for _ in 0..count {
                wanted.push(item.clone());
            }
            i += 1;
        }
    }
    let dentate: usize = wanted.iter().map(|w| w.donors.len()).sum();
    if dentate != sites.len() {
        return Err(format!("the ligands give {dentate} donor atoms; a {} complex has {} sites", o.geometry, sites.len()));
    }
    // Bidentate ligands first (they need cis site pairs).
    wanted.sort_by(|a, b| b.donors.len().cmp(&a.donors.len()));
    let mut result = ComplexResult::default();
    let mut m = Molecule::default();
    m.add_atom(Atom { z: zm, charge: o.oxidation_state, ..Default::default() });
    for w in &wanted {
        let eo = EmbedOptions { attempts: 2, ..Default::default() };
        let (mut lig, _) = molecule_from_smiles(&w.smiles, &eo).map_err(|e| e.message)?;
        // Donor indices count heavy atoms in SMILES order, which are the first atoms.
        for &d in &w.donors {
            if d < 0 || d as usize >= lig.atoms.len() || lig.atoms[d as usize].z == 1 {
                return Err(format!("donor index {d} is not a heavy atom of {}", w.smiles));
            }
        }
        let distance_to = |lig: &Molecule, donor: i32| if o.distance > 0.0 { o.distance } else { covalent_radius(zm) + covalent_radius(lig.atoms[donor as usize].z) };
        let mut placed = vec![[0.0; 3]; lig.atoms.len()];
        if w.donors.len() == 1 {
            let site = used.iter().position(|&u| !u).ok_or("no free site")?;
            used[site] = true;
            let d = w.donors[0];
            let pd = lig.atoms[d as usize].position;
            let mut lone = [0.0; 3];
            for x in lig.neighbours(d) {
                lone = sub(lone, unit(sub(lig.atoms[x as usize].position, pd)));
            }
            if norm(lone) < 1e-6 {
                lone = if lig.degree(d) > 0 { unit(sub(pd, lig.atoms[lig.neighbours(d)[0] as usize].position)) } else { [0.0, 0.0, 1.0] };
            }
            lone = unit(lone);
            let v = unit(sites[site]);
            let target = scale(v, distance_to(&lig, d));
            for (i, a) in lig.atoms.iter().enumerate() {
                placed[i] = add(target, align_vector(sub(a.position, pd), lone, scale(v, -1.0)));
            }
            // Turn about the M-L axis away from the atoms already placed.
            let (mut best_score, mut best_angle) = (1e300, 0.0);
            for step in 0..24 {
                let angle = f64::from(step) * PI / 12.0;
                let mut score = 0.0;
                for &pi in &placed {
                    let x = add(target, rotate_about(sub(pi, target), v, angle));
                    for a in &m.atoms[1..] {
                        let r = norm(sub(x, a.position));
                        if r < 3.0 {
                            score += (3.0 - r) * (3.0 - r);
                        }
                    }
                }
                if score < best_score - 1e-9 {
                    best_score = score;
                    best_angle = angle;
                }
            }
            for x in placed.iter_mut() {
                *x = add(target, rotate_about(sub(*x, target), v, best_angle));
            }
        } else if w.donors.len() == 2 {
            let (d1, d2) = (w.donors[0], w.donors[1]);
            // Turn the single bonds between the donors until they face each
            // other (a chelating conformation: syn bipyridine, gauche en).
            {
                let mut previous = vec![-2i32; lig.atoms.len()];
                let mut q = VecDeque::from([d1]);
                previous[d1 as usize] = -1;
                while let Some(a) = q.pop_front() {
                    for b in lig.neighbours(a) {
                        if previous[b as usize] == -2 {
                            previous[b as usize] = a;
                            q.push_back(b);
                        }
                    }
                }
                let mut path = Vec::new();
                let mut a = d2;
                while a >= 0 {
                    path.push(a);
                    a = previous[a as usize];
                }
                path.reverse();
                let mut rotors: Vec<[i32; 4]> = Vec::new();
                let mut s2 = 1;
                while s2 + 2 < path.len() + 1 && s2 + 1 < path.len() {
                    let (j, k) = (path[s2], path[s2 + 1]);
                    let e = lig.bond_between(j, k);
                    if e >= 0 && lig.bonds[e as usize].order == 1 && !side_of(&lig, j, k).is_empty() && s2 + 2 < path.len() {
                        rotors.push([path[s2 - 1], j, k, path[s2 + 2]]);
                    }
                    s2 += 1;
                }
                rotors.truncate(3);
                let start = lig.clone();
                let mut best = lig.clone();
                let mut best_score = 1e300;
                let combos = 12usize.pow(rotors.len() as u32);
                for c in 0..combos {
                    let mut trial = start.clone();
                    let mut code = c;
                    for r in &rotors {
                        set_dihedral(&mut trial, r[0], r[1], r[2], r[3], -180.0 + 30.0 * (code % 12) as f64)?;
                        code /= 12;
                    }
                    let tp = |i: i32| trial.atoms[i as usize].position;
                    let mut score = (norm(sub(tp(d1), tp(d2))) - 2.75).abs();
                    let n = trial.atoms.len() as i32;
                    for x in 0..n {
                        for y in x + 1..n {
                            if trial.atoms[x as usize].z == 1 && trial.atoms[y as usize].z == 1 {
                                continue;
                            }
                            let r = norm(sub(tp(x), tp(y)));
                            let mut close = trial.bond_between(x, y) >= 0;
                            for z in trial.neighbours(x) {
                                close = close || trial.bond_between(z, y) >= 0;
                            }
                            if !close && r < 2.4 {
                                score += 2.4 - r;
                            }
                        }
                    }
                    if score < best_score {
                        best_score = score;
                        best = trial;
                    }
                }
                lig = best;
            }
            let (r1, r2) = (distance_to(&lig, d1), distance_to(&lig, d2));
            let (p1, p2) = (lig.atoms[d1 as usize].position, lig.atoms[d2 as usize].position);
            let bite = norm(sub(p1, p2));
            let bite_angle = 2.0 * fmin(1.0, bite / (r1 + r2)).asin() * 180.0 / PI;
            let (mut s1, mut s2) = (-1i32, -1i32);
            let mut best_diff = 1e300;
            for a in 0..sites.len() {
                for b in 0..sites.len() {
                    if a == b || used[a] || used[b] {
                        continue;
                    }
                    let angle = dot(unit(sites[a]), unit(sites[b])).clamp(-1.0, 1.0).acos() * 180.0 / PI;
                    let diff = (angle - bite_angle).abs();
                    if diff < best_diff - 1e-9 {
                        best_diff = diff;
                        s1 = a as i32;
                        s2 = b as i32;
                    }
                }
            }
            if s1 < 0 {
                return Err("no free pair of sites for a bidentate ligand".into());
            }
            let (s1, s2) = (s1 as usize, s2 as usize);
            used[s1] = true;
            used[s2] = true;
            // Superpose (donor 1, donor 2, centroid) onto (site 1, site 2, outward point).
            let mut centroid = [0.0; 3];
            let mut heavy = 0;
            for a in &lig.atoms {
                if a.z != 1 {
                    centroid = add(centroid, a.position);
                    heavy += 1;
                }
            }
            centroid = scale(centroid, 1.0 / f64::from(heavy));
            let mid = scale(add(p1, p2), 0.5);
            let out = norm(sub(centroid, mid));
            let (t1, t2) = (scale(unit(sites[s1]), r1), scale(unit(sites[s2]), r2));
            let t3 = add(scale(add(t1, t2), 0.5), scale(unit(add(sites[s1], sites[s2])), fmax(out, 0.5)));
            let n = lig.atoms.len();
            let mut a = vec![[0.0; 3]; n + 1];
            let mut b: Vec<Vec3> = lig.atoms.iter().map(|x| x.position).collect();
            b.push(if norm(sub(centroid, mid)) > 0.3 { centroid } else { add(mid, [0.0, 0.0, 1.0]) });
            a[d1 as usize] = t1;
            a[d2 as usize] = t2;
            a[n] = t3;
            let mut best_clash = 1e300;
            for swap in 0..2 {
                let mut ta = a.clone();
                if swap == 1 {
                    ta[d1 as usize] = t2;
                    ta[d2 as usize] = t1;
                }
                let mut aligned = Vec::new();
                rmsd(&ta, &b, &[d1, d2, n as i32], Some(&mut aligned));
                let mut clash = 0.0;
                for al in aligned.iter().take(n) {
                    for ma in &m.atoms {
                        let r = norm(sub(*al, ma.position));
                        if r < 2.5 {
                            clash += 2.5 - r;
                        }
                    }
                }
                if clash < best_clash {
                    best_clash = clash;
                    placed[..n].copy_from_slice(&aligned[..n]);
                }
            }
        } else {
            return Err("only mono- and bidentate ligands are supported".into());
        }
        let mut map = vec![0i32; lig.atoms.len()];
        for (i, la) in lig.atoms.iter().enumerate() {
            let mut a = la.clone();
            a.position = placed[i];
            map[i] = m.add_atom(a);
        }
        for b in &lig.bonds {
            let e = m.add_bond(map[b.a as usize], map[b.b as usize], b.order)?;
            m.bonds[e as usize].aromatic = b.aromatic;
        }
        for &d in &w.donors {
            let e = m.add_bond(0, map[d as usize], 1)?;
            m.bonds[e as usize].dative = true;
        }
    }
    // Relax the ligands (UFF) with the metal and the donor atoms held, which
    // removes contacts between neighbouring ligands.
    let relaxed = (|| -> Result<Vec<Vec3>, String> {
        let field = UffForceField::new(&m, None, ForceFieldOptions::default())?;
        let mut x: Vec<Vec3> = m.atoms.iter().map(|a| a.position).collect();
        let mut fixed = vec![0u8; m.atoms.len()];
        fixed[0] = 7;
        for &k in m.bonds_of(0) {
            fixed[m.other(k, 0) as usize] = 7;
        }
        minimize(&field, &mut x, &fixed, 600, 0.05)?;
        Ok(x)
    })();
    if let Ok(x) = relaxed {
        for (a, p) in m.atoms.iter_mut().zip(x) {
            a.position = p;
        }
    }
    result.charge = total_charge(&m);
    // Group of the metal (d block 3-12): Sc, Y, Lu start group 3.
    let group = [21, 39, 71].iter().find(|&&s| zm >= s && zm < s + 10).map(|&s| 3 + (zm - s)).unwrap_or(if zm == 57 { 3 } else { 0 });
    if group >= 3 {
        result.d_electrons = group - o.oxidation_state;
        let d = result.d_electrons;
        let mut unpaired = 0;
        if (0..=10).contains(&d) {
            let low = squeeze(&o.spin) == "low";
            if geometry == "squareplanar" && d == 8 {
                unpaired = 0;
            } else if !low || geometry == "tetrahedral" {
                unpaired = if d <= 5 { d } else { 10 - d };
            } else {
                const LOW_SPIN: [i32; 11] = [0, 1, 2, 3, 2, 1, 0, 1, 2, 1, 0];
                unpaired = LOW_SPIN[d as usize];
            }
        }
        result.multiplicity = unpaired + 1;
    } else {
        result.multiplicity = if electron_count(&m) % 2 != 0 { 2 } else { 1 };
    }
    perceive_aromaticity(&mut m);
    result.molecule = m;
    Ok(result)
}

// ------------------------------------------------------------------ Z-matrices

fn trim(s: &str) -> &str {
    s.trim_matches(|c| c == ' ' || c == '\t')
}

/// Reads a Z-matrix (Gaussian style: element or label, then bonded atom and
/// distance, angle atom and angle, dihedral atom and dihedral; references by
/// number or label; values or variables defined after a blank line or a
/// "Variables:" line; X dummy atoms dropped).
pub fn parse_z_matrix(text: &str) -> Result<Molecule, String> {
    // std::getline lines (no empty last line after a final newline).
    let mut lines: Vec<&str> = text.split('\n').collect();
    if lines.last() == Some(&"") {
        lines.pop();
    }
    let lines: Vec<&str> = lines.into_iter().map(|l| l.strip_suffix('\r').unwrap_or(l)).collect();
    let mut start = 0;
    // A Gaussian input: skip the Link 0 and route sections, the title and the charge line.
    let gaussian = lines.iter().map(|l| trim(l)).find(|l| !l.is_empty()).is_some_and(|l| l.starts_with('#') || l.starts_with('%'));
    if gaussian {
        let mut blanks = 0;
        let mut i = 0;
        while i < lines.len() && trim(lines[i]).is_empty() {
            i += 1;
        }
        while i < lines.len() && blanks < 2 {
            if trim(lines[i]).is_empty() {
                blanks += 1;
                while i + 1 < lines.len() && trim(lines[i + 1]).is_empty() {
                    i += 1;
                }
            }
            i += 1;
        }
        start = i;
    }
    // Atom lines until a blank line or "Variables:"; then variable definitions.
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut variables: BTreeMap<String, f64> = BTreeMap::new();
    let mut in_variables = false;
    for line in &lines[start.min(lines.len())..] {
        let mut l = trim(line).to_string();
        if let Some(comment) = l.find('!') {
            l = trim(&l[..comment]).to_string();
        }
        let lower = l.to_ascii_lowercase();
        if lower.starts_with("variables") || lower.starts_with("constants") {
            in_variables = true;
            continue;
        }
        if l.is_empty() {
            if !rows.is_empty() {
                in_variables = true;
            }
            continue;
        }
        let l: String = l.chars().map(|c| if matches!(c, ',' | '=' | '\t') { ' ' } else { c }).collect();
        let w = words(&l);
        if in_variables {
            if w.len() >= 2 {
                let v = stod(&w[1]).map_err(|_| "stod".to_string())?.0;
                variables.insert(w[0].clone(), v);
            }
            continue;
        }
        // A "charge multiplicity" line before the first atom.
        if rows.is_empty() && w.len() == 2 && w[0].bytes().all(|c| c.is_ascii_digit() || c == b'-') {
            continue;
        }
        rows.push(w);
    }
    if rows.is_empty() {
        return Err("the Z-matrix has no atoms".into());
    }
    let mut labels: Vec<String> = Vec::new();
    let mut elements: Vec<i32> = Vec::new();
    let mut positions: Vec<Vec3> = Vec::new();
    let value = |word: &str, row: usize| -> Result<f64, String> {
        if let Ok((v, used)) = stod(word) {
            if used == word.len() {
                return Ok(v);
            }
        }
        let negative = word.starts_with('-');
        let name = if negative { &word[1..] } else { word };
        match variables.get(name) {
            Some(&v) => Ok(if negative { -v } else { v }),
            None => Err(format!("line {}: undefined variable '{name}'", row + 1)),
        }
    };
    for (r, w) in rows.iter().enumerate() {
        let reference = |word: &str, labels: &[String]| -> Result<usize, String> {
            let mut index = -1i64;
            if let Some((v, used)) = stoi_used(word) {
                if used == word.len() {
                    index = i64::from(v) - 1;
                }
            }
            if index < 0 {
                for (k, l) in labels.iter().enumerate() {
                    if l == word {
                        index = k as i64;
                    }
                }
            }
            if index < 0 || index >= r as i64 {
                return Err(format!("line {}: '{word}' is not an earlier atom", r + 1));
            }
            Ok(index as usize)
        };
        let symbol: String = w[0].chars().take_while(|c| c.is_ascii_alphabetic()).collect();
        let z = if symbol == "X" || symbol == "x" || symbol == "Xx" {
            -1
        } else {
            let mut z = atomic_number(&symbol);
            if z <= 0 && symbol.len() > 1 {
                z = atomic_number(&symbol[..1]);
            }
            if z <= 0 {
                return Err(format!("line {}: unknown element '{}'", r + 1, w[0]));
            }
            z
        };
        let mut p = [0.0; 3];
        let needed = match r {
            0 => 1,
            1 => 3,
            2 => 5,
            _ => 7,
        };
        if w.len() == 4 && r >= 1 && w.len() < needed {
            // Cartesian line: element x y z.
            p = [value(&w[1], r)?, value(&w[2], r)?, value(&w[3], r)?];
        } else {
            if w.len() < needed {
                return Err(format!("line {}: expected {needed} fields", r + 1));
            }
            if r == 1 {
                p = [0.0, 0.0, value(&w[2], r)?];
            }
            if r >= 1 {
                let a = reference(&w[1], &labels)?;
                let bond = value(&w[2], r)?;
                if r == 1 {
                    p = add(positions[a], [0.0, 0.0, bond]);
                } else {
                    let b = reference(&w[3], &labels)?;
                    let angle = value(&w[4], r)?;
                    if r == 2 {
                        // In the xz plane.
                        let u = unit(sub(positions[b], positions[a]));
                        let mut perp = cross(u, [0.0, 1.0, 0.0]);
                        if norm(perp) < 1e-6 {
                            perp = cross(u, [1.0, 0.0, 0.0]);
                        }
                        perp = unit(cross(perp, u));
                        let t = angle * PI / 180.0;
                        p = add(positions[a], scale(add(scale(u, t.cos()), scale(perp, t.sin())), bond));
                    } else {
                        let c = reference(&w[5], &labels)?;
                        let torsion = value(&w[6], r)?;
                        p = place_atom(positions[c], positions[b], positions[a], bond, angle, torsion);
                    }
                }
            }
        }
        labels.push(w[0].clone());
        elements.push(z);
        positions.push(p);
    }
    let mut m = Molecule::default();
    for (i, &z) in elements.iter().enumerate() {
        if z < 0 {
            continue;
        }
        m.add_atom(Atom { z, position: positions[i], ..Default::default() });
    }
    Ok(m)
}

/// Writes a Z-matrix from coordinates (references along bonds where possible);
/// with variables, the values follow as R1, A1, D1 definitions.
pub fn write_z_matrix(m: &Molecule, variables: bool) -> String {
    let n = m.atoms.len() as i32;
    let mut out = String::new();
    let mut vars = String::new();
    let p = |i: i32| m.atoms[i as usize].position;
    let distance = |a: i32, b: i32| norm(sub(p(a), p(b)));
    let angle_of = |a: i32, b: i32, c: i32| {
        let (u, v) = (unit(sub(p(a), p(b))), unit(sub(p(c), p(b))));
        dot(u, v).clamp(-1.0, 1.0).acos() * 180.0 / PI
    };
    // References: a bonded earlier atom when there is one, else the nearest earlier atom.
    let pick = |i: i32, exclude: &[i32]| {
        let mut best = -1;
        let mut best_score = 1e300;
        for j in 0..i {
            if exclude.contains(&j) {
                continue;
            }
            let mut score = distance(i, j);
            if m.bond_between(i, j) >= 0 {
                score -= 100.0;
            }
            if score < best_score {
                best_score = score;
                best = j;
            }
        }
        best
    };
    for i in 0..n {
        out += &symbol_of(m.atoms[i as usize].z);
        if i >= 1 {
            let a = pick(i, &[]);
            let (mut b, mut c) = (-1, -1);
            if i >= 2 {
                // b bonded to a where possible, not collinear with i-a.
                let mut best = 1e300;
                for j in 0..i {
                    if j == a {
                        continue;
                    }
                    let ang = angle_of(i, a, j);
                    if ang < 5.0 || ang > 175.0 {
                        continue;
                    }
                    let score = distance(a, j) - if m.bond_between(a, j) >= 0 { 100.0 } else { 0.0 };
                    if score < best {
                        best = score;
                        b = j;
                    }
                }
                if b < 0 {
                    b = pick(i, &[a]);
                }
            }
            if i >= 3 {
                let mut best = 1e300;
                for j in 0..i {
                    if j == a || j == b {
                        continue;
                    }
                    let ang = angle_of(a, b, j);
                    if ang < 5.0 || ang > 175.0 {
                        continue;
                    }
                    let score = distance(b, j) - if m.bond_between(b, j) >= 0 { 100.0 } else { 0.0 };
                    if score < best {
                        best = score;
                        c = j;
                    }
                }
                if c < 0 {
                    c = pick(i, &[a, b]);
                }
            }
            let r = distance(i, a);
            if variables {
                let _ = write!(out, "  {}  R{}", a + 1, i + 1);
                let _ = writeln!(vars, "R{} = {r:.6}", i + 1);
            } else {
                let _ = write!(out, "  {}  {r:.6}", a + 1);
            }
            if b >= 0 {
                let t = angle_of(i, a, b);
                if variables {
                    let _ = write!(out, "  {}  A{}", b + 1, i + 1);
                    let _ = writeln!(vars, "A{} = {t:.4}", i + 1);
                } else {
                    let _ = write!(out, "  {}  {t:.4}", b + 1);
                }
            }
            if c >= 0 {
                let d = dihedral_angle(p(i), p(a), p(b), p(c));
                if variables {
                    let _ = write!(out, "  {}  D{}", c + 1, i + 1);
                    let _ = writeln!(vars, "D{} = {d:.4}", i + 1);
                } else {
                    let _ = write!(out, "  {}  {d:.4}", c + 1);
                }
            }
        }
        out += "\n";
    }
    if variables {
        out += "\n";
        out += &vars;
    }
    out
}

// ------------------------------------------------------------------ quantum-chemistry inputs

/// Formats: orca, nwchem, psi4, qchem, mopac, turbomole, gzmat, gamess, zmat.
pub fn is_quantum_format(format: &str) -> bool {
    matches!(format.to_ascii_lowercase().as_str(), "orca" | "nwchem" | "psi4" | "qchem" | "mopac" | "turbomole" | "gzmat" | "gamess" | "zmat")
}

/// Charge and multiplicity from the per-atom charge field when it holds integer
/// formal charges, else from the perceived formal charges and radicals.
pub fn charge_and_multiplicity(s: &Structure) -> (i32, i32) {
    let mut sum = 0.0;
    let (mut integral, mut any) = (true, false);
    for a in &s.atoms {
        sum += a.charge;
        if a.charge != 0.0 {
            any = true;
        }
        if (a.charge - a.charge.round()).abs() > 1e-6 {
            integral = false;
        }
    }
    let m = perceive(s, &PerceiveOptions::default()).unwrap_or_default();
    if any && integral {
        let charge = sum.round() as i32;
        let electrons: i32 = s.atoms.iter().map(|a| a.atomic_number).sum::<i32>() - charge;
        return (charge, if electrons % 2 != 0 { 2 } else { 1 });
    }
    (total_charge(&m), if m.atoms.is_empty() { 1 } else { multiplicity(&m) })
}

/// The input file for a quantum-chemistry code; Err with a message for an
/// empty structure or an unknown format.
pub fn write_quantum_input(s: &Structure, format: &str) -> Result<String, String> {
    let f = format.to_ascii_lowercase();
    if s.atoms.is_empty() {
        return Err("The structure has no atoms.".into());
    }
    let (q, mult) = charge_and_multiplicity(s);
    let periodic = s.has_unit_cell && s.pbc.iter().any(|&p| p);
    let symbol = |a: &af_core::structure::AtomSite| if a.atomic_number > 0 { af_core::elements::symbol(a.atomic_number).to_string() } else { a.symbol.clone() };
    let atom_lines = |indent: &str| {
        let mut o = String::new();
        for a in &s.atoms {
            let _ = writeln!(o, "{indent}{:<3}{:>16.8}{:>16.8}{:>16.8}", symbol(a), a.x, a.y, a.z);
        }
        o
    };
    let note = format!("charge {q} and multiplicity {mult} from the formal charges (check them); method and basis are placeholders");
    let mut out = String::new();
    match f.as_str() {
        "orca" => {
            let _ = writeln!(out, "# ORCA input written by AtomForge: {note}");
            if periodic {
                out += "# The structure is periodic; ORCA treats it as a molecule.\n";
            }
            let _ = write!(out, "! B3LYP D3BJ def2-SVP Opt{}\n\n%pal nprocs 4 end\n%maxcore 2000\n\n", if mult > 1 { " UKS" } else { "" });
            let _ = write!(out, "* xyz {q} {mult}\n{}*\n", atom_lines("  "));
        }
        "nwchem" => {
            let _ = write!(out, "# NWChem input written by AtomForge: {note}\nstart molecule\ntitle \"AtomForge structure\"\ncharge {q}\n\n");
            let _ = write!(out, "geometry units angstrom noautosym\n{}", atom_lines("  "));
            if periodic {
                out += "  system crystal\n";
                let c = &s.cell_vectors;
                let len = |v: &[f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                let ang = |u: &[f64; 3], v: &[f64; 3]| ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (len(u) * len(v))).acos() * 180.0 / PI;
                let _ = write!(
                    out,
                    "    lat_a {:.6} lat_b {:.6} lat_c {:.6}\n    alpha {:.6} beta {:.6} gamma {:.6}\n  end\n",
                    len(&c[0]),
                    len(&c[1]),
                    len(&c[2]),
                    ang(&c[1], &c[2]),
                    ang(&c[0], &c[2]),
                    ang(&c[0], &c[1])
                );
            }
            let _ = write!(out, "end\n\nbasis\n  * library def2-svp\nend\n\ndft\n  xc b3lyp\n  mult {mult}\nend\n\ntask dft optimize\n");
        }
        "psi4" => {
            let _ = write!(
                out,
                "# Psi4 input written by AtomForge: {note}\nmolecule {{\n{q} {mult}\n{}  units angstrom\n  no_reorient\n}}\n\nset {{\n  basis def2-svp\n  reference {}\n}}\n\noptimize('b3lyp-d3bj')\n",
                atom_lines("  "),
                if mult > 1 { "uks" } else { "rks" }
            );
        }
        "qchem" => {
            let _ = write!(
                out,
                "$comment\nQ-Chem input written by AtomForge: {note}\n$end\n\n$molecule\n{q} {mult}\n{}$end\n\n$rem\n  JOBTYPE   opt\n  METHOD    b3lyp\n  BASIS     def2-svp\n  DFT_D     d3_bj\n{}$end\n",
                atom_lines("  "),
                if mult > 1 { "  UNRESTRICTED true\n" } else { "" }
            );
        }
        "mopac" => {
            const SPINS: [&str; 7] = ["", "SINGLET", "DOUBLET", "TRIPLET", "QUARTET", "QUINTET", "SEXTET"];
            let spin = if mult > 1 && mult <= 6 { format!(" UHF {}", SPINS[mult as usize]) } else { String::new() };
            let _ = write!(out, "PM7 CHARGE={q}{spin} PRECISE\nAtomForge structure\n{note}\n");
            for a in &s.atoms {
                let _ = writeln!(out, "{:<3}{:>16.8} 1{:>16.8} 1{:>16.8} 1", symbol(a), a.x, a.y, a.z);
            }
            if periodic {
                for k in 0..3 {
                    if s.pbc[k] {
                        let v = s.cell_vectors[k];
                        let _ = writeln!(out, "Tv {:>16.8} 1{:>16.8} 1{:>16.8} 1", v[0], v[1], v[2]);
                    }
                }
            }
        }
        "turbomole" => {
            let bohr = 1.0 / 0.52917721092;
            out += "$coord\n";
            for a in &s.atoms {
                let _ = writeln!(out, "{:>20.10}{:>20.10}{:>20.10}  {}", a.x * bohr, a.y * bohr, a.z * bohr, symbol(a).to_ascii_lowercase());
            }
            if periodic {
                let dims = s.pbc.iter().filter(|&&p| p).count();
                let _ = write!(out, "$periodic {dims}\n$lattice angs\n");
                for k in 0..3 {
                    if s.pbc[k] {
                        let v = s.cell_vectors[k];
                        let _ = writeln!(out, "  {:.10}  {:.10}  {:.10}", v[0], v[1], v[2]);
                    }
                }
            }
            out += "$end\n";
        }
        "gzmat" | "zmat" => {
            let mut m = Molecule::default();
            for a in &s.atoms {
                let z = if a.atomic_number > 0 { a.atomic_number } else { atomic_number(&a.symbol) };
                m.add_atom(Atom { z, position: [a.x, a.y, a.z], ..Default::default() });
            }
            if let Ok(p) = perceive(s, &PerceiveOptions { bond_orders: false, bond_scale: 1.2 }) {
                m.bonds = p.bonds;
                m.reindex();
            }
            if f == "gzmat" {
                let _ = write!(out, "%nprocshared=4\n%mem=4GB\n# B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt\n\nAtomForge structure: {note}\n\n{q} {mult}\n");
            }
            out += &write_z_matrix(&m, true);
            if f == "gzmat" {
                out += "\n";
            }
        }
        "gamess" => {
            let _ = writeln!(out, "! GAMESS input written by AtomForge: {note}");
            let _ = writeln!(out, " $CONTRL SCFTYP={} DFTTYP=B3LYP RUNTYP=OPTIMIZE ICHARG={q} MULT={mult} $END", if mult > 1 { "UHF" } else { "RHF" });
            out += " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n $DATA\nAtomForge structure\nC1\n";
            for a in &s.atoms {
                let _ = writeln!(out, "{:<3}{:>6.1}{:>16.8}{:>16.8}{:>16.8}", symbol(a), f64::from(a.atomic_number), a.x, a.y, a.z);
            }
            out += " $END\n";
        }
        _ => return Err(format!("Unknown quantum-chemistry format {format}")),
    }
    Ok(out)
}
