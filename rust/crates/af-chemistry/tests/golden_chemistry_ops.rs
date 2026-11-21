//! af-chemistry operations against the C++ chemistry library: the cases of
//! tools/golden/golden_chemistry_ops.cpp (tests/golden/chemistry_ops.json),
//! run in the same order. Everything is compared exactly except the
//! distance-geometry embeddings, whose coordinates come from a different
//! eigensolver: there the graph, the unspecified stereocentres and the
//! stereochemistry must match, and the UFF energy of the kept embedding must
//! be close (docs/rust-port-plan.md section 6a).

use std::path::PathBuf;

use af_chemistry::*;
use af_core::json::Json;
use af_core::structure::{AtomSite, Structure};

fn vec3(v: Vec3) -> Json {
    Json::Array(vec![v[0].into(), v[1].into(), v[2].into()])
}

fn positions(x: &[Vec3]) -> Json {
    Json::Array(x.iter().map(|&v| vec3(v)).collect())
}

fn graph_json(m: &Molecule) -> Json {
    let mut o = Json::object();
    o.set("atoms", Json::Array(m.atoms.iter().map(|a| Json::Array(vec![a.z.into(), a.charge.into(), a.hydrogens.into(), a.aromatic.into(), a.dummy.into()])).collect()));
    o.set("bonds", Json::Array(m.bonds.iter().map(|b| Json::Array(vec![b.a.into(), b.b.into(), b.order.into(), b.aromatic.into(), b.dative.into()])).collect()));
    o.set("tetrahedral", Json::Array(m.tetrahedral.iter().map(|t| Json::Array(vec![t.centre.into(), t.n[0].into(), t.n[1].into(), t.n[2].into(), t.n[3].into(), t.sign.into()])).collect()));
    o.set("double_bonds", Json::Array(m.double_bonds.iter().map(|d| Json::Array(vec![d.a.into(), d.x.into(), d.y.into(), d.b.into(), d.cis.into()])).collect()));
    o.set("formula", formula(m));
    o
}

fn molecule_json(m: &Molecule) -> Json {
    let mut o = graph_json(m);
    o.set("positions", positions(&coordinates(m)));
    o.set("images", Json::Array(m.bonds.iter().map(|b| vec3(b.image)).collect()));
    o
}

fn ints(v: &[i32]) -> Json {
    Json::Array(v.iter().map(|&x| x.into()).collect())
}

fn doubles(v: &[f64]) -> Json {
    Json::Array(v.iter().map(|&x| x.into()).collect())
}

fn strings<S: AsRef<str>>(v: &[S]) -> Json {
    Json::Array(v.iter().map(|s| s.as_ref().into()).collect())
}

fn terms_json(t: &EnergyTerms) -> Json {
    Json::Array(vec![t.bond.into(), t.angle.into(), t.torsion.into(), t.inversion.into(), t.vdw.into(), t.electrostatic.into(), t.restraint.into(), t.total().into()])
}

fn minimize_json(r: &MinimizeResult) -> Json {
    Json::Array(vec![r.initial_energy.into(), r.energy.into(), r.steps.into(), r.converged.into(), r.max_force.into()])
}

fn group_json(g: &PointGroupResult) -> Json {
    let mut o = Json::object();
    o.set("name", g.name.as_str());
    o.set("order", g.order);
    o.set("rotational", g.rotational_symmetry_number);
    o.set("operations", Json::Array(g.operations.iter().map(|m| Json::Array(m.iter().flat_map(|r| r.iter().map(|&v| v.into())).collect())).collect()));
    o.set("centre", vec3(g.centre));
    o.set("principal", vec3(g.principal_axis));
    o.set("secondary", vec3(g.secondary_axis));
    o.set("deviation", g.max_deviation);
    o
}

fn coordinates(m: &Molecule) -> Vec<Vec3> {
    m.atoms.iter().map(|a| a.position).collect()
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

fn load(name: &str) -> Structure {
    read_xyz(&format!("rust/tests/fixtures/chemistry/{name}.xyz"))
}

fn perceived(name: &str) -> Molecule {
    perceive(&load(name), &PerceiveOptions::default()).unwrap()
}

type Outcome = Result<Json, String>;

struct Cases(Vec<Json>);

impl Cases {
    fn attempt(&mut self, name: &str, result: Outcome) {
        let mut entry = Json::object();
        entry.set("name", name);
        match result {
            Ok(r) => entry.set("result", r),
            Err(e) => entry.set("error", e.as_str()),
        }
        self.0.push(entry);
    }
}

fn smiles_error(e: SmilesError) -> String {
    e.message
}

fn force_field_cases(c: &mut Cases) {
    for name in ["aspirin", "caffeine", "benzene", "ethanol", "chfclbr", "z_butene", "ammonium_acetate", "polyethylene", "butane", "water", "sf6"] {
        let s = load(name);
        let result = (|| -> Outcome {
            let m = perceive(&s, &PerceiveOptions::default())?;
            let cell = if s.has_unit_cell { Some(&s) } else { None };
            let mut o = Json::object();
            let field = UffForceField::new(&m, cell, ForceFieldOptions::default())?;
            o.set("types", strings(field.types()));
            let x = coordinates(&m);
            let mut g = Vec::new();
            let mut terms = EnergyTerms::default();
            o.set("energy", field.energy(&x, Some(&mut g), Some(&mut terms))?);
            o.set("terms", terms_json(&terms));
            o.set("gradient", positions(&g));
            let eo = ForceFieldOptions { electrostatics: true, dielectric: 4.0, cutoff: 9.0 };
            let charged = UffForceField::new(&m, cell, eo)?;
            let mut ct = EnergyTerms::default();
            let mut cg = Vec::new();
            o.set("charged_energy", charged.energy(&x, Some(&mut cg), Some(&mut ct))?);
            o.set("charged_terms", terms_json(&ct));
            o.set("charged_gradient", positions(&cg));
            let mut restrained = UffForceField::new(&m, cell, ForceFieldOptions::default())?;
            let n = m.atoms.len();
            restrained.restraints.push(Restraint { kind: 2, atoms: [0, 1, 0, 0], target: 1.7, k: 300.0 });
            if n >= 3 {
                restrained.restraints.push(Restraint { kind: 3, atoms: [0, 1, 2, 0], target: 100.0, k: 50.0 });
            }
            if n >= 4 {
                restrained.restraints.push(Restraint { kind: 4, atoms: [0, 1, 2, 3], target: 30.0, k: 20.0 });
            }
            let mut rt = EnergyTerms::default();
            let mut rg = Vec::new();
            o.set("restrained_energy", restrained.energy(&x, Some(&mut rg), Some(&mut rt))?);
            o.set("restrained_terms", terms_json(&rt));
            o.set("restrained_gradient", positions(&rg));
            let mut y = x.clone();
            o.set("minimize", minimize_json(&minimize(&field, &mut y, &[], 300, 0.01)?));
            o.set("minimized", positions(&y));
            let mut z = x.clone();
            let mut fixed = vec![0u8; n];
            fixed[0] = 7;
            fixed[1] = 1;
            o.set("minimize_fixed", minimize_json(&minimize(&restrained, &mut z, &fixed, 100, 0.05)?));
            o.set("minimized_fixed", positions(&z));
            Ok(o)
        })();
        c.attempt(&format!("uff {name}"), result);
    }
    let tables = (|| -> Outcome {
        let mut o = Json::object();
        o.set("has", Json::Array((0..=120).map(|z| uff_has_element(z).into()).collect()));
        let mut rows = Vec::new();
        for t in ["C_3", "C_R", "C_2", "C_1", "N_3", "N_R", "O_3", "O_2", "H_", "S_3+2", "Fe3+2", "Cl", "P_3+5", "Si3"] {
            rows.push(Json::Array(vec![
                t.into(),
                uff_angle(t).into(),
                uff_vdw_radius(t).into(),
                uff_bond_length(t, "C_3", 1.0)?.into(),
                uff_bond_length(t, "O_2", 2.0)?.into(),
                uff_bond_length(t, "C_R", 1.5)?.into(),
            ]));
        }
        o.set("rows", Json::Array(rows));
        Ok(o)
    })();
    c.attempt("uff tables", tables);
    c.attempt("uff unknown type", uff_bond_length("Qq", "C_3", 1.0).map(Json::from));
    c.attempt(
        "uff implicit hydrogens",
        (|| -> Outcome {
            UffForceField::new(&parse_smiles("CCO").map_err(smiles_error)?, None, ForceFieldOptions::default())?;
            Ok(0.into())
        })(),
    );
}

fn point_group_cases(c: &mut Cases) {
    for name in [
        "methane", "sf6", "water", "ammonia", "pcl5", "ethane", "trans_dichloroethylene", "allene", "co2", "hcn", "benzene", "benzene_distorted", "ethanol", "aspirin",
        "caffeine", "chfclbr",
    ] {
        let s = load(name);
        let z: Vec<i32> = s.atoms.iter().map(|a| a.atomic_number).collect();
        let x: Vec<Vec3> = s.atoms.iter().map(|a| [a.x, a.y, a.z]).collect();
        let mut o = Json::object();
        let g = point_group(&z, &x, 0.1);
        o.set("group", group_json(&g));
        o.set("tight", group_json(&point_group(&z, &x, 0.01)));
        let mut sym = x.clone();
        o.set("symmetrized_group", group_json(&symmetrize(&z, &mut sym, 0.1)));
        o.set("symmetrized", positions(&sym));
        let mut oriented = x.clone();
        standard_orientation(&mut oriented, &g);
        o.set("oriented", positions(&oriented));
        c.attempt(&format!("point group {name}"), Ok(o));
    }
    c.attempt("point group single", Ok(group_json(&point_group(&[6], &[[1.0, 2.0, 3.0]], 0.1))));
    c.attempt("point group empty", Ok(group_json(&point_group(&[], &[], 0.1))));
}

fn edit_case(c: &mut Cases, name: &str, file: &str, work: impl FnOnce(&mut Molecule) -> Outcome) {
    let mut m = perceived(file);
    let result = work(&mut m).map(|value| {
        let mut o = Json::object();
        o.set("value", value);
        o.set("molecule", molecule_json(&m));
        o
    });
    c.attempt(name, result);
}

fn editing_cases(c: &mut Cases) {
    for name in ["aspirin", "caffeine", "ethanol", "ammonium_acetate", "z_butene", "benzene"] {
        let m = perceived(name);
        let mut o = Json::object();
        o.set("six", doubles(&gasteiger_charges(&m, 6)));
        o.set("two", doubles(&gasteiger_charges(&m, 2)));
        c.attempt(&format!("gasteiger {name}"), Ok(o));
    }
    {
        let (a, b) = (perceived("chfclbr"), perceived("chfclbr_mirror"));
        let mut o = Json::object();
        let mut aligned = Vec::new();
        o.set("mirror", rmsd(&coordinates(&a), &coordinates(&b), &[], Some(&mut aligned)));
        o.set("aligned", positions(&aligned));
        o.set("subset", rmsd(&coordinates(&a), &coordinates(&b), &[0, 1, 2], None));
        let cm = perceived("butane");
        let turned: Vec<Vec3> = coordinates(&cm).iter().map(|p| [0.3 * p[0] - 0.954 * p[1] + 1.0, 0.954 * p[0] + 0.3 * p[1] - 2.0, p[2] + 0.5]).collect();
        o.set("turned", rmsd(&coordinates(&cm), &turned, &[], Some(&mut aligned)));
        o.set("turned_aligned", positions(&aligned));
        c.attempt("rmsd", Ok(o));
    }
    c.attempt(
        "butane torsions",
        (|| -> Outcome {
            let mut m = perceived("butane");
            let mut o = Json::object();
            o.set("side12", ints(&side_of(&m, 1, 2)));
            o.set("side21", ints(&side_of(&m, 2, 1)));
            let p = coordinates(&m);
            o.set("dihedral", dihedral_angle(p[0], p[1], p[2], p[3]));
            set_dihedral(&mut m, 0, 1, 2, 3, 60.0)?;
            o.set("set60", positions(&coordinates(&m)));
            set_dihedral(&mut m, 0, 1, 2, 3, -75.5)?;
            o.set("set-75", positions(&coordinates(&m)));
            Ok(o)
        })(),
    );
    c.attempt(
        "dihedral in ring",
        (|| -> Outcome {
            let mut m = perceived("benzene");
            set_dihedral(&mut m, 0, 1, 2, 3, 10.0)?;
            Ok(0.into())
        })(),
    );
    edit_case(c, "bond order 2", "ethanol", |m| Ok(set_bond_order(m, 0, 1, 2)?.into()));
    edit_case(c, "bond order 3", "ethanol", |m| Ok(set_bond_order(m, 0, 1, 3)?.into()));
    edit_case(c, "bond order 0", "ethanol", |m| Ok(set_bond_order(m, 0, 1, 0)?.into()));
    edit_case(c, "bond order new", "ethanol", |m| Ok(set_bond_order(m, 0, 2, 1)?.into()));
    edit_case(c, "bond order hydrogen", "ethanol", |m| {
        let mut h = -1;
        for x in m.neighbours(0) {
            if m.atoms[x as usize].z == 1 {
                h = x;
            }
        }
        Ok(set_bond_order(m, 0, h, 2)?.into())
    });
    edit_case(c, "bond order same", "ethanol", |m| Ok(set_bond_order(m, 1, 1, 1)?.into()));
    edit_case(c, "bond order butene", "z_butene", |m| Ok(set_bond_order(m, 0, 1, 1)?.into()));
    edit_case(c, "add hydrogens", "aspirin", |m| {
        let keep: Vec<u8> = m.atoms.iter().map(|a| u8::from(a.z != 1)).collect();
        m.remove_atoms(&keep);
        Ok(add_hydrogens(m, &[])?.into())
    });
    edit_case(c, "add hydrogens masked", "caffeine", |m| {
        let keep: Vec<u8> = m.atoms.iter().map(|a| u8::from(a.z != 1)).collect();
        m.remove_atoms(&keep);
        let mut mask = vec![0u8; m.atoms.len()];
        mask[0] = 1;
        mask[2] = 1;
        Ok(add_hydrogens(m, &mask)?.into())
    });
    edit_case(c, "remove hydrogens", "aspirin", |m| Ok(remove_hydrogens(m, &[]).into()));
    edit_case(c, "remove hydrogens masked", "ethanol", |m| {
        let mut mask = vec![0u8; m.atoms.len()];
        mask[0] = 1;
        Ok(remove_hydrogens(m, &mask).into())
    });
    edit_case(c, "protonate ammonia", "ammonia", |m| protonate(m, 0).map(|_| 0.into()));
    edit_case(c, "protonate water", "water", |m| protonate(m, 0).map(|_| 0.into()));
    edit_case(c, "protonate hydrogen", "water", |m| protonate(m, 1).map(|_| 0.into()));
    edit_case(c, "deprotonate ethanol", "ethanol", |m| {
        let o = m.atoms.iter().rposition(|a| a.z == 8).map_or(-1, |i| i as i32);
        deprotonate(m, o)?;
        Ok(o.into())
    });
    edit_case(c, "deprotonate none", "sf6", |m| deprotonate(m, 0).map(|_| 0.into()));
    edit_case(c, "ionize aspirin", "aspirin", |m| Ok(ionize_at_neutral_ph(m)?.into()));
    edit_case(c, "ionize caffeine", "caffeine", |m| Ok(ionize_at_neutral_ph(m)?.into()));
    edit_case(c, "invert chfclbr", "chfclbr", |m| {
        invert_stereocentre(m, 0)?;
        Ok(stereo_matches(m).into())
    });
    edit_case(c, "invert three bonds", "ammonia", |m| invert_stereocentre(m, 0).map(|_| 0.into()));
    edit_case(c, "flip butene", "z_butene", |m| flip_double_bond(m, 0, 1).map(|_| 0.into()));
    edit_case(c, "flip unbonded", "z_butene", |m| flip_double_bond(m, 2, 3).map(|_| 0.into()));
    edit_case(c, "flip ring", "benzene", |m| flip_double_bond(m, 0, 1).map(|_| 0.into()));
    for s in ["NCC(=O)O", "CCN", "NC(N)=N", "OS(=O)(=O)c1ccccc1", "OP(=O)(O)O", "CN(C)C", "Nc1ccccc1", "CC(N)=O", "C[NH3+]", "OCCN"] {
        let result = (|| -> Outcome {
            let mut m = parse_smiles(s).map_err(smiles_error)?;
            make_hydrogens_explicit(&mut m, true);
            let mut o = Json::object();
            o.set("text", ionize_at_neutral_ph(&mut m)?);
            o.set("molecule", molecule_json(&m));
            Ok(o)
        })();
        c.attempt(&format!("ionize {s}"), result);
    }
    c.attempt(
        "group smiles",
        (|| -> Outcome {
            let mut o = Vec::new();
            for g in ["Me", "methyl", "tBu", "TBU", "ph", "COOH", "coo-", "*CC", "nitro", "Tert_Butyl", "hydrogen", "b(oh)2", "OAc"] {
                o.push(group_smiles(g)?.into());
            }
            o.push(Json::Array(functional_groups().iter().map(|(n, s)| Json::Array(vec![(*n).into(), (*s).into()])).collect()));
            Ok(Json::Array(o))
        })(),
    );
    c.attempt("group unknown", group_smiles("unobtainium").map(Json::from));
    for (file, atom, group) in [
        ("ethanol", 3, "methyl"),
        ("ethanol", 3, "phenyl"),
        ("ethanol", 4, "*C(F)(F)F"),
        ("benzene", 6, "nitro"),
        ("benzene", 7, "COOH"),
        ("ethanol", 0, "methyl"),
        ("ethanol", 3, "unknown-group"),
        ("ethanol", 3, "*C*"),
        ("ethanol", 3, "*=C"),
        ("ethanol", 99, "methyl"),
    ] {
        let mut m = perceived(file);
        let result = substitute(&mut m, atom, group).map(|r| {
            let mut o = Json::object();
            o.set("new_atoms", ints(&r.new_atoms));
            o.set("attachment", r.attachment);
            o.set("graph", graph_json(&m));
            o
        });
        c.attempt(&format!("substitute {file} {atom} {group}"), result);
    }
    for (file, method, trials) in [("butane", "random", 12), ("butane", "systematic", 0), ("ethanol", "random", 8), ("aspirin", "random", 4), ("benzene", "random", 5)] {
        let m = perceived(file);
        let co = ConformerOptions { method: method.into(), trials, seed: 3, optimize_steps: 150, max_conformers: 4, rmsd_threshold: 0.3, ..Default::default() };
        let result = conformer_search(&m, &co).map(|(found, rotors)| {
            let mut o = Json::object();
            o.set("rotors", rotors);
            o.set("conformers", Json::Array(found.iter().map(|c| Json::Array(vec![c.energy.into(), positions(&c.positions)])).collect()));
            o
        });
        c.attempt(&format!("conformers {file} {method}"), result);
    }
    let scans: [(&[i32], f64, f64, f64, bool); 8] = [
        (&[0, 1, 2, 3], 0.0, 180.0, 60.0, true),
        (&[0, 1, 2, 3], 180.0, -180.0, 90.0, false),
        (&[0, 1], 1.4, 1.7, 0.1, false),
        (&[0, 1], 1.4, 1.7, 0.15, true),
        (&[0, 1, 2], 100.0, 120.0, 10.0, false),
        (&[0, 1, 2], 100.0, 120.0, 10.0, true),
        (&[0], 0.0, 1.0, 1.0, false),
        (&[0, 1], 0.0, 1.0, 0.0, false),
    ];
    for (atoms, start, stop, step, relaxed) in scans {
        let m = perceived("butane");
        let mut name = "scan".to_string();
        for a in atoms {
            name += &format!(" {a}");
        }
        name += &format!("{}{start:.6} {step:.6}", if relaxed { " relaxed " } else { " rigid " });
        let result = geometry_scan(&m, atoms, start, stop, step, relaxed, 80)
            .map(|points| Json::Array(points.iter().map(|p| Json::Array(vec![p.value.into(), p.energy.into(), positions(&p.positions)])).collect()));
        c.attempt(&name, result);
    }
}

fn builder_cases(c: &mut Cases) {
    c.attempt(
        "compound library",
        Ok(Json::Array(compound_library().iter().map(|x| Json::Array(vec![x.name.into(), x.category.into(), x.smiles.into()])).collect())),
    );
    c.attempt(
        "compound smiles",
        Ok(Json::Array(
            ["water", "Ethanol", "diethyl-ether", "DIETHYL_ETHER", "meoh", "Tylenol", "bpy", "lys", "nothing", "", "1,4-dioxane", "14dioxane"]
                .iter()
                .map(|n| compound_smiles(n).into())
                .collect(),
        )),
    );
    c.attempt("polymer library", Ok(Json::Array(polymer_library().iter().map(|(n, s)| Json::Array(vec![(*n).into(), (*s).into()])).collect())));
    c.attempt(
        "ligand library",
        Ok(Json::Array(ligand_library().iter().map(|l| Json::Array(vec![l.name.into(), l.smiles.into(), ints(l.donors), l.description.into()])).collect())),
    );
    let peptides: [(&str, &str, &str, bool, bool, f64, f64, f64); 12] = [
        ("ACDEFGHIKLMNPQRSTVWY", "alpha", "charged", true, false, 0.0, 0.0, 0.0),
        ("Ala-Gly-Ser", "beta", "neutral", true, false, 0.0, 0.0, 0.0),
        ("GPG", "3-10", "capped", true, false, 0.0, 0.0, 0.0),
        ("KRDE", "extended", "charged", false, false, 0.0, 0.0, 0.0),
        ("AW", "custom", "neutral", true, false, -60.0, -40.0, 175.0),
        ("AA", "ppii", "charged", true, true, 0.0, 0.0, 0.0),
        ("HIE-CYX,ash glh LYN", "pi", "capped", true, false, 0.0, 0.0, 0.0),
        ("y q", "parallel-beta", "Charged", true, false, 0.0, 0.0, 0.0),
        ("", "alpha", "charged", true, false, 0.0, 0.0, 0.0),
        ("AXZ", "alpha", "charged", true, false, 0.0, 0.0, 0.0),
        ("AG", "spiral", "charged", true, false, 0.0, 0.0, 0.0),
        ("AG", "alpha", "zwitter", true, false, 0.0, 0.0, 0.0),
    ];
    for (sequence, conformation, termini, ionize, d, phi, psi, omega) in peptides {
        let mut o = PeptideOptions {
            sequence: sequence.into(),
            conformation: conformation.into(),
            termini: termini.into(),
            ionize_side_chains: ionize,
            d_amino: d,
            ..Default::default()
        };
        if conformation == "custom" {
            (o.phi, o.psi, o.omega) = (phi, psi, omega);
        }
        let result = build_peptide(&o).map(|(m, labels)| {
            let mut r = Json::object();
            r.set("molecule", molecule_json(&m));
            r.set("labels", strings(&labels));
            r
        });
        c.attempt(&format!("peptide {sequence} {conformation} {termini}{}", if d { " D" } else { "" }), result);
    }
    let polymers: [(&str, i32, &str, &str, &str, &str, bool, u32); 12] = [
        ("PE", 5, "isotactic", "trans", "H", "H", false, 1),
        ("PS", 4, "syndiotactic", "trans", "CH3", "H", false, 1),
        ("pp", 6, "atactic", "trans", "H", "phenyl", false, 7),
        ("PE", 6, "isotactic", "random", "H", "H", false, 5),
        ("PE", 4, "isotactic", "trans", "H", "H", true, 1),
        ("*CC(Cl)*", 3, "isotactic", "trans", "methyl", "methyl", false, 2),
        ("PS", 3, "syndiotactic", "trans", "H", "H", true, 1),
        ("*CC", 3, "isotactic", "trans", "H", "H", false, 1),
        ("PE", 0, "isotactic", "trans", "H", "H", false, 1),
        ("PE", 3, "random", "trans", "H", "H", false, 1),
        ("PE", 3, "isotactic", "helix", "H", "H", false, 1),
        ("PET", 2, "isotactic", "trans", "H", "H", false, 1),
    ];
    for (monomer, units, tacticity, conformation, head, tail, periodic, seed) in polymers {
        let o = PolymerOptions {
            monomer: monomer.into(),
            units,
            tacticity: tacticity.into(),
            conformation: conformation.into(),
            head_group: head.into(),
            tail_group: tail.into(),
            periodic,
            seed,
        };
        let result = build_polymer(&o).map(|r| {
            let mut j = Json::object();
            j.set("graph", graph_json(&r.molecule));
            j.set("monomer", r.monomer_smiles.as_str());
            j
        });
        c.attempt(&format!("polymer {monomer} {units} {tacticity} {conformation} {head} {tail}{}", if periodic { " periodic" } else { "" }), result);
    }
    let complexes: [(&str, i32, &str, &str, &str, f64); 19] = [
        ("Fe", 2, "octahedral", "NH3 6", "high", 0.0),
        ("Co", 3, "octahedral", "en 3", "low", 0.0),
        ("Ni", 2, "square-planar", "Cl 4", "high", 0.0),
        ("Cu", 2, "tetrahedral", "Cl 4", "low", 0.0),
        ("Pt", 2, "square planar", "NH3 2 Cl 2", "high", 2.05),
        ("Ru", 2, "octahedral", "bipy 3", "low", 0.0),
        ("Ni", 2, "octahedral", "smiles:NCCN:0,3 2 H2O 2", "high", 0.0),
        ("Zn", 2, "tetrahedral", "py 2 Cl 2", "high", 0.0),
        ("Fe", 3, "trigonal-bipyramidal", "CO 5", "high", 0.0),
        ("V", 4, "square-pyramidal", "O 1 F 4", "high", 0.0),
        ("Na", 1, "linear", "H2O 2", "high", 0.0),
        ("Cr", 3, "octahedral", "acac 3", "high", 0.0),
        ("Xx", 2, "octahedral", "NH3 6", "high", 0.0),
        ("Fe", 2, "octahedral", "NH3 5", "high", 0.0),
        ("Fe", 2, "hexagonal", "NH3 6", "high", 0.0),
        ("Fe", 2, "octahedral", "foo 6", "high", 0.0),
        ("Fe", 2, "octahedral", "smiles:NCCN", "high", 0.0),
        ("Fe", 2, "octahedral", "smiles:NCC(N)CN:0,4,5 2", "high", 0.0),
        ("Fe", 2, "octahedral", "smiles:NCCN:0,9 3", "high", 0.0),
    ];
    for (metal, ox, geometry, ligands, spin, distance) in complexes {
        let o = ComplexOptions { metal: metal.into(), oxidation_state: ox, geometry: geometry.into(), ligands: ligands.into(), spin: spin.into(), distance };
        let result = build_complex(&o).map(|r| {
            let mut j = Json::object();
            j.set("graph", graph_json(&r.molecule));
            j.set("charge", r.charge);
            j.set("d", r.d_electrons);
            j.set("multiplicity", r.multiplicity);
            j
        });
        c.attempt(&format!("complex {metal} {geometry} {ligands} {spin}"), result);
    }
    let zmatrices = [
        "O\nH 1 0.96\nH 1 0.96 2 104.5\n",
        "O\nH 1 r\nH 1 r 2 a\n\nr 0.97\na=105.0\n",
        "C\nC 1 1.54\nH 1 1.09 2 109.5\nH 1 1.09 2 109.5 3 120\nH 1 1.09 2 109.5 3 -120\nH 2 1.09 1 109.5 3 60\nH 2 1.09 1 109.5 6 120\nH 2 1.09 1 109.5 6 -120\n",
        "%chk=x.chk\n# HF/STO-3G\n\ntitle line\n\n0 1\nC1\nO2 C1 rco\nH3 C1 rch O2 ang\nH4 C1 rch O2 ang H3 180.0\nVariables:\nrco 1.2\nrch 1.1\nang 121.0\n",
        "X\nN 1 1.0\nH 2 1.01 1 90.0\nH 2 1.01 1 90.0 3 120.0\nH 2 1.01 1 90.0 3 -120.0 ! comment\n",
        "C\nO 1 1.2\nH 0.5 0.9 0.1\nCl 1 1.8 2 -r 3 90\n\nr 120\n",
        "O\nH 1 r\n",
        "O\nH 2 0.96\n",
        "Q\nH 1 1\n",
        "O\nH 1\n",
        "",
        "  \n\n",
        "Fe\nCL 1 2.2\nCa 1 2.4 2 90\n",
        "O\nH 1 r\n\nr abc\n",
        "C\r\nH 1 1.09\r\nH 1 1.09 2 109.47\r\nH 1 1.09 2 109.47 3 120\r\n",
        "C\nC 1 1.5\nC 2 1.5 1 180\nC 3 1.5 2 120 1 0\n",
    ];
    for (k, text) in zmatrices.iter().enumerate() {
        let result = parse_z_matrix(text).map(|m| {
            let mut o = Json::object();
            o.set("molecule", molecule_json(&m));
            o.set("written", write_z_matrix(&m, false));
            o
        });
        c.attempt(&format!("zmatrix {k}"), result);
    }
    for name in ["aspirin", "ethanol", "water", "co2", "hcn", "sf6"] {
        let m = perceived(name);
        c.attempt(&format!("write zmatrix {name}"), Ok(Json::Array(vec![write_z_matrix(&m, false).into(), write_z_matrix(&m, true).into()])));
    }
    c.attempt(
        "quantum formats",
        Ok(Json::Array(["orca", "ORCA", "nwchem", "psi4", "qchem", "mopac", "turbomole", "gzmat", "gamess", "zmat", "xyz", ""].iter().map(|f| is_quantum_format(f).into()).collect())),
    );
    let mut ions = load("ammonium_acetate");
    ions.atoms[0].charge = 1.0;
    ions.atoms[8].charge = -1.0;
    let mut partial = load("ethanol");
    partial.atoms[0].charge = 0.25;
    let mut radical = load("water");
    radical.atoms.pop();
    let mut unknown = load("water");
    unknown.atoms[0].atomic_number = 0;
    unknown.atoms[0].symbol = "O".into();
    let structures = [
        ("ethanol", load("ethanol")),
        ("polyethylene", load("polyethylene")),
        ("ammonium_acetate", load("ammonium_acetate")),
        ("ions", ions),
        ("partial", partial),
        ("radical", radical),
        ("empty", Structure::default()),
        ("unknown", unknown),
    ];
    for (name, s) in &structures {
        let (q, mult) = charge_and_multiplicity(s);
        c.attempt(&format!("charge {name}"), Ok(Json::Array(vec![q.into(), mult.into()])));
        for f in ["orca", "nwchem", "psi4", "qchem", "mopac", "turbomole", "gzmat", "gamess", "zmat", "cp2k"] {
            let row = match write_quantum_input(s, f) {
                Ok(text) => vec![true.into(), text.into(), "".into()],
                Err(e) => vec![false.into(), "".into(), e.into()],
            };
            c.attempt(&format!("quantum {name} {f}"), Ok(Json::Array(row)));
        }
    }
}

const EMBED_SMILES: [&str; 24] = [
    "CCO",
    "c1ccccc1",
    "C[C@H](N)C(=O)O",
    "C[C@@H](N)C(=O)O",
    "F/C=C/F",
    "F/C=C\\F",
    "CC(C)(C)C",
    "OC[C@@H](O)[C@H](O)CO",
    "C1CCCCC1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "C",
    "[Na+].[Cl-]",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "C1CC1",
    "N#N",
    "C=C=C",
    "CC(N)C(C)O",
    "C/C=C/C=C\\C",
    "[C@H](F)(Cl)Br",
    "*CC*",
    "[Fe+2]",
    "C1CC2CCC1C2",
    "[Xx]",
    "C(",
];

fn embedding_cases(c: &mut Cases) {
    for s in EMBED_SMILES {
        for optimize in [true, false] {
            let eo = EmbedOptions { attempts: 3, seed: 11, optimize, ..Default::default() };
            let result = molecule_from_smiles(s, &eo).map_err(smiles_error).map(|(m, r)| {
                let mut o = Json::object();
                o.set("graph", graph_json(&m));
                o.set("energy", r.energy);
                o.set("succeeded", r.succeeded);
                o.set("unspecified", r.unspecified_stereo);
                o.set("stereo", stereo_matches(&m));
                o
            });
            c.attempt(&format!("embed {s}{}", if optimize { "" } else { " plain" }), result);
        }
    }
    c.attempt(
        "embed implicit hydrogens",
        (|| -> Outcome {
            let mut m = parse_smiles("CCO").map_err(smiles_error)?;
            embed(&mut m, &EmbedOptions::default())?;
            Ok(0.into())
        })(),
    );
    c.attempt("embed empty", embed(&mut Molecule::default(), &EmbedOptions::default()).map(|_| 0.into()));
}

/// Embedding results: exact graph, stereo flags and counts; energy (after the
/// UFF minimisation) within a few kcal/mol, since each run may settle in a
/// different conformer from its own random metric-matrix start.
fn embedding_matches(got: &Json, want: &Json) -> Result<(), String> {
    let (Ok(g), Ok(w)) = (got.at("result"), want.at("result")) else {
        return if got == want { Ok(()) } else { Err("error differs".into()) };
    };
    for key in ["graph", "unspecified", "stereo"] {
        if g.at(key) != w.at(key) {
            return Err(format!("{key} differs"));
        }
    }
    let (gs, ws) = (g.at("succeeded").unwrap().number().unwrap(), w.at("succeeded").unwrap().number().unwrap());
    if (ws > 0.0) != (gs > 0.0) {
        return Err(format!("succeeded {gs} vs {ws}"));
    }
    let name = got.at("name").unwrap().string().unwrap();
    if !name.ends_with(" plain") {
        let (ge, we) = (g.at("energy").unwrap().number().unwrap(), w.at("energy").unwrap().number().unwrap());
        if std::env::var_os("AF_VERBOSE").is_some() {
            eprintln!("{name}: energy {ge} (C++ {we})");
        }
        // A higher energy than the C++ embedding: within 2 kcal/mol or 10 %; a
        // lower one (another conformer of a flexible molecule) within 5 or 25 %.
        if ge > we + f64::max(2.0, 0.1 * we.abs()) || ge < we - f64::max(5.0, 0.25 * we.abs()) {
            return Err(format!("energy {ge} vs {we}"));
        }
    }
    Ok(())
}

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn chemistry_operations_match_cpp() {
    std::env::set_current_dir(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")).unwrap();
    let mut c = Cases(Vec::new());
    force_field_cases(&mut c);
    point_group_cases(&mut c);
    editing_cases(&mut c);
    builder_cases(&mut c);
    embedding_cases(&mut c);

    let golden = Json::parse(&std::fs::read_to_string("rust/tests/golden/chemistry_ops.json").unwrap()).unwrap();
    let expected = golden.items().unwrap();
    assert_eq!(c.0.len(), expected.len(), "case count");
    let mut differences = Vec::new();
    for (got, want) in c.0.iter().zip(expected) {
        let name = want.at("name").unwrap().string().unwrap();
        assert_eq!(got.at("name").unwrap().string().unwrap(), name, "case order");
        if name.starts_with("embed ") {
            if let Err(e) = embedding_matches(got, want) {
                differences.push(format!("{name}: {e}"));
            }
            continue;
        }
        if got != want {
            let (g, w) = (got.dump(-1), want.dump(-1));
            let first = g.bytes().zip(w.bytes()).position(|(a, b)| a != b).unwrap_or(g.len().min(w.len()));
            let from = first.saturating_sub(100);
            differences.push(format!("{name}:\n  {}\n  vs {}", &g[from..(first + 100).min(g.len())], &w[from..(first + 100).min(w.len())]));
        }
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.iter().take(25).cloned().collect::<Vec<_>>().join("\n"));
}
