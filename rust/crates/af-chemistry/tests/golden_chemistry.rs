//! af-chemistry against the C++ chemistry library: the cases of
//! tools/golden/golden_chemistry.cpp (tests/golden/chemistry.json), run in the
//! same order and compared exactly.

use std::path::PathBuf;

use af_chemistry::*;
use af_core::json::Json;
use af_core::structure::{AtomSite, Structure};

fn molecule_json(m: &Molecule) -> Json {
    let mut o = Json::object();
    o.set(
        "atoms",
        Json::Array(
            m.atoms
                .iter()
                .map(|a| {
                    Json::Array(vec![
                        a.z.into(),
                        a.charge.into(),
                        a.hydrogens.into(),
                        a.aromatic.into(),
                        a.isotope.into(),
                        a.map_number.into(),
                        a.dummy.into(),
                        a.position[0].into(),
                        a.position[1].into(),
                        a.position[2].into(),
                    ])
                })
                .collect(),
        ),
    );
    o.set(
        "bonds",
        Json::Array(
            m.bonds
                .iter()
                .map(|b| Json::Array(vec![b.a.into(), b.b.into(), b.order.into(), b.aromatic.into(), b.dative.into(), b.image[0].into(), b.image[1].into(), b.image[2].into()]))
                .collect(),
        ),
    );
    o.set("tetrahedral", Json::Array(m.tetrahedral.iter().map(|t| Json::Array(vec![t.centre.into(), t.n[0].into(), t.n[1].into(), t.n[2].into(), t.n[3].into(), t.sign.into()])).collect()));
    o.set("double_bonds", Json::Array(m.double_bonds.iter().map(|d| Json::Array(vec![d.a.into(), d.x.into(), d.y.into(), d.b.into(), d.cis.into()])).collect()));
    o
}

fn ints(v: &[i32]) -> Json {
    Json::Array(v.iter().map(|&x| Json::from(x)).collect())
}

fn rings(list: &[Vec<i32>]) -> Json {
    Json::Array(list.iter().map(|r| ints(r)).collect())
}

fn properties(m: &Molecule) -> Json {
    let mut o = Json::object();
    o.set("formula", formula(m));
    o.set("molar_mass", molar_mass(m));
    o.set("monoisotopic_mass", monoisotopic_mass(m));
    o.set("charge", total_charge(m));
    o.set("electrons", electron_count(m));
    o.set("multiplicity", multiplicity(m));
    o.set("fragments", rings(&fragments(m)));
    o.set("rotatable", ints(&rotatable_bonds(m)));
    o.set("ranks", ints(&canonical_ranks(m)));
    let r = smallest_rings(m);
    o.set("rings", rings(&r));
    o.set("ring_sizes", ints(&smallest_ring_size(m, &r)));
    let n = m.atoms.len() as i32;
    o.set("hybridisation", ints(&(0..n).map(|i| hybridisation(m, i)).collect::<Vec<_>>()));
    o.set("valences", ints(&(0..n).map(|i| m.valence(i)).collect::<Vec<_>>()));
    o.set("hydrogen_counts", ints(&(0..n).map(|i| m.hydrogen_count(i)).collect::<Vec<_>>()));
    o.set("cip", Json::Array(cip_labels(m).iter().map(|l| Json::Array(vec![l.atom.into(), l.other.into(), l.label.as_str().into()])).collect()));
    o.set("smiles", write_smiles(m, false));
    o.set("smiles_stereo", write_smiles(m, true));
    o
}

fn read_xyz(path: &str) -> Structure {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let n: usize = lines.next().unwrap().trim().parse().unwrap();
    let comment = lines.next().unwrap();
    let mut s = Structure::default();
    if let Some(values) = comment.strip_prefix("cell ") {
        let v: Vec<f64> = values.split_whitespace().map(|x| x.parse().unwrap()).collect();
        s.cell_vectors = [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]];
        s.has_unit_cell = true;
        s.pbc = [true; 3];
    }
    for _ in 0..n {
        let w: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        s.atoms.push(AtomSite {
            symbol: w[0].to_string(),
            x: w[1].parse().unwrap(),
            y: w[2].parse().unwrap(),
            z: w[3].parse().unwrap(),
            atomic_number: atomic_number(w[0]),
            ..Default::default()
        });
    }
    s
}

struct Cases(Vec<Json>);

impl Cases {
    fn attempt(&mut self, name: &str, result: Result<Json, SmilesError>) {
        let mut entry = Json::object();
        entry.set("name", name);
        match result {
            Ok(r) => entry.set("result", r),
            Err(e) => {
                entry.set("error", e.message.as_str());
                if let Some(p) = e.position {
                    entry.set("position", p as f64);
                }
            }
        }
        self.0.push(entry);
    }
}

fn plain(e: String) -> SmilesError {
    SmilesError { message: e, position: None }
}

#[test]
fn chemistry_matches_cpp() {
    std::env::set_current_dir(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")).unwrap();
    let input = "rust/tests/fixtures/chemistry/";
    let mut c = Cases(Vec::new());
    let smiles = [
        "CCO", "c1ccccc1", "C1=CC=CC=C1", "c1ccncc1", "c1cc[nH]c1", "c1ccnc1", "o1cccc1", "s1cccc1", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O",
        "N[C@@H](C)C(=O)O", "F/C=C/F", "F/C=C\\F", "C/C=C(/C)C", "[NH4+]", "CC(=O)[O-]", "[Na+].[Cl-]", "[13CH4]", "[CH3:1]C", "C%10CC%10",
        "c1ccc2ccccc2c1", "O=C1NC(=O)C=C1", "[Fe+2]", "*CC*", "CC(C)(C)C", "C#N", "O=S(=O)(O)O", "P(=O)(O)(O)O", "C[N+](=O)[O-]",
        "[O-][N+](=O)c1ccccc1", "OC[C@@H](O)[C@H](O)CO", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C1CCCCCCC1", "C/1=C/CCCCCC1", "CC.CC",
        "C-1CC1", "[2H]O", "C$C", "c1ccccc1-c2ccccc2", "[C@H](F)(Cl)Br", "CCO ethanol", "[CH]", "[C]", "[O]", "O=O", "[Cu+2].[O-2]",
        "[se]1cccc1", "c1ccc[se]1", "[as]1cccc1", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "[nH]1cccc1", "C[S@](=O)c1ccccc1", "[C@@](F)(Cl)(Br)I",
        "N#[N+][O-]", "[B-](F)(F)(F)F", "[Si](C)(C)(C)C", "c1cccc1", "C1CC", "C(", "C)", "(C)", "[Xx]", "[c]1ccccc1", "Q", "%1", "C%1",
        "C%1a", "1C", "C11", "[CH3", "[C@SP1H4]", "CC==C", "", " ", "C.", "C1.C1", "C=1CC-1", "C=1CCC#1", "C[NH3+]", "[O-]C(=O)C",
    ];
    for s in smiles {
        let result = parse_smiles(s).map(|m| {
            let mut o = Json::object();
            o.set("molecule", molecule_json(&m));
            o.set("properties", properties(&m));
            let mut h = m.clone();
            o.set("explicit_added", make_hydrogens_explicit(&mut h, true));
            o.set("explicit", molecule_json(&h));
            o.set("explicit_smiles", write_smiles(&h, true));
            o
        });
        c.attempt(&format!("smiles {s}"), result);
    }
    for name in ["aspirin", "caffeine", "benzene", "ethanol", "chfclbr", "chfclbr_mirror", "z_butene", "e_butene", "ammonium_acetate", "polyethylene"] {
        let s = read_xyz(&format!("{input}{name}.xyz"));
        let result = (|| -> Result<Json, String> {
            let mut o = Json::object();
            let m = perceive(&s, &PerceiveOptions::default())?;
            o.set("molecule", molecule_json(&m));
            o.set("properties", properties(&m));
            o.set("plain", molecule_json(&perceive(&s, &PerceiveOptions { bond_orders: false, bond_scale: 1.3 })?));
            // Hydrogens removed, counted back and placed by hybridisation.
            let mut heavy = m.clone();
            let keep: Vec<u8> = heavy.atoms.iter().map(|a| u8::from(a.z != 1)).collect();
            o.set("map", ints(&heavy.remove_atoms(&keep)));
            assign_implicit_hydrogens(&mut heavy, &[]);
            o.set("placed_added", make_hydrogens_explicit(&mut heavy, true));
            o.set("placed", molecule_json(&heavy));
            o.set(
                "structure",
                Json::Array(
                    to_structure(&m, true)
                        .atoms
                        .iter()
                        .map(|a| Json::Array(vec![a.symbol.as_str().into(), a.atomic_number.into(), a.x.into(), a.y.into(), a.z.into(), a.charge.into(), a.r.into(), a.g.into(), a.b.into()]))
                        .collect(),
                ),
            );
            o.set(
                "bond_vectors",
                Json::Array((1..s.atoms.len().min(4)).map(|j| {
                    let v = bond_vector(&s, 0, j);
                    Json::Array(vec![v[0].into(), v[1].into(), v[2].into()])
                }).collect()),
            );
            Ok(o)
        })();
        c.attempt(&format!("perceive {name}"), result.map_err(plain));
    }
    let mut elements = Vec::new();
    for z in 0..=118 {
        let mut row = vec![symbol_of(z).into(), is_metal(z).into(), covalent_radius(z).into(), hydrogen_bond_length(z).into(), atomic_number(&symbol_of(z)).into()];
        for s in 0..8 {
            row.push(default_valence(z, s % 3 - 1, s).into());
        }
        elements.push(Json::Array(row));
    }
    for s in ["cl", "CL", "D", "T", "Xx", ""] {
        elements.push(Json::Array(vec![s.into(), atomic_number(s).into()]));
    }
    c.attempt("elements", Ok(Json::Array(elements)));

    let golden = Json::parse(&std::fs::read_to_string("rust/tests/golden/chemistry.json").unwrap()).unwrap();
    let expected = golden.items().unwrap();
    assert_eq!(c.0.len(), expected.len(), "case count");
    let mut differences = Vec::new();
    for (got, want) in c.0.iter().zip(expected) {
        if got != want {
            let name = want.at("name").unwrap().string().unwrap();
            let (g, w) = (got.dump(-1), want.dump(-1));
            let first = g.bytes().zip(w.bytes()).position(|(a, b)| a != b).unwrap_or(g.len().min(w.len()));
            let from = first.saturating_sub(80);
            differences.push(format!("{name}:\n  {}\n  vs {}", &g[from..(first + 120).min(g.len())], &w[from..(first + 120).min(w.len())]));
        }
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.iter().take(15).cloned().collect::<Vec<_>>().join("\n"));
}
