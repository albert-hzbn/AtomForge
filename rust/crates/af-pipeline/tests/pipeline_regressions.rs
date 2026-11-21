//! Port of the C++ tests/pipeline_regressions.cpp (the parts covering the
//! steps ported so far; insert-interstitials, simulation steps and the
//! options reader follow with their crates).

use af_core::json::Json;
use af_core::{AtomSite, Structure};
use af_pipeline::*;

fn atom(symbol: &str, z: i32, x: f64, y: f64, w: f64) -> AtomSite {
    AtomSite { symbol: symbol.into(), atomic_number: z, x, y, z: w, ..Default::default() }
}

/// Conventional fcc cell(s) of copper, a = 3.615.
fn copper(repeat: i32) -> Structure {
    let a = 3.615;
    let r = repeat as f64;
    let mut s = Structure { has_unit_cell: true, cell_vectors: [[a * r, 0.0, 0.0], [0.0, a * r, 0.0], [0.0, 0.0, a * r]], ..Default::default() };
    for i in 0..repeat {
        for j in 0..repeat {
            for k in 0..repeat {
                for f in [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]] {
                    s.atoms.push(atom("Cu", 29, (i as f64 + f[0]) * a, (j as f64 + f[1]) * a, (k as f64 + f[2]) * a));
                }
            }
        }
    }
    s
}

fn run(text: &str, input: &Structure) -> Result<PipelineData, String> {
    let (out, stages) = Pipeline::parse(text)?.evaluate(input)?;
    if let Some(stage) = stages.iter().find(|s| !s.error.is_empty()) {
        return Err(stage.error.clone());
    }
    Ok(out)
}

fn ok(text: &str, input: &Structure) -> PipelineData {
    run(text, input).unwrap_or_else(|e| panic!("{text}: {e}"))
}

fn close(a: f64, b: f64, tolerance: f64, what: &str) {
    assert!((a - b).abs() <= tolerance, "{what}: {a} vs {b}");
}

fn evaluate(expression: &str, v: &AtomVariables) -> f64 {
    Expression::new(expression).unwrap().evaluate(v).unwrap()
}

#[test]
fn selection_expressions() {
    let v = AtomVariables { x: 2.0, y: -1.0, z: 4.0, fz: 0.75, element: "O".into(), atomic_number: 8.0, index: 3.0, ..Default::default() };
    let none = AtomVariables::default();
    close(evaluate("1 + 2 * 3 ^ 2", &v), 19.0, 0.0, "precedence");
    close(evaluate("-2 ^ 2", &none), -4.0, 0.0, "unary minus binds weaker than power");
    close(evaluate("(1 + 2) * 3", &none), 9.0, 0.0, "parentheses");
    close(evaluate("7 % 3", &none), 1.0, 0.0, "remainder");
    close(evaluate("-7 % 3", &none), 2.0, 0.0, "floored remainder of a negative value");
    close(evaluate("2 * 7 % 4", &none), 2.0, 0.0, "remainder at the precedence of * and /");
    let data = AtomVariables { magmom: -2.2, charge: 1.5, vx: 0.01, vy: 0.02, vz: -0.03, ..Default::default() };
    close(evaluate("magmom + charge + vx + vy + vz", &data), -0.7, 1e-12, "per-atom simulation data variables");
    close(evaluate("x * y + z", &v), 2.0, 0.0, "variables");
    assert_eq!(evaluate("fz > 0.5 && element == O", &v), 1.0, "element without quotes");
    assert_eq!(evaluate("element == \"O\" and Z == 8", &v), 1.0, "quoted element and keywords");
    assert_eq!(evaluate("element != 'O' or x > 100", &v), 0.0, "single quotes and or");
    assert!(evaluate("!(x > 1)", &v) == 0.0 && evaluate("not x > 5", &v) == 1.0, "negation");
    assert_eq!(evaluate("x != 2", &v), 0.0, "!= is not negation");
    assert_eq!(evaluate("index >= 3 && index <= 3", &v), 1.0, "comparisons");
    assert!(!evaluate("property", &v).is_nan(), "property defaults to zero in an empty variable set");
    assert!(evaluate("element == o", &v) == 1.0 && evaluate("element == O", &v) == 1.0, "element symbols in any case");
    assert!(Expression::new("elemnt == O").is_err(), "a misspelt variable is an error, not an element name");
    assert!(Expression::new("unknown(1)").is_err(), "unknown function");
    assert!(Expression::new("(x > 1").is_err(), "missing parenthesis");
    assert!(Expression::new("x > ").is_err(), "incomplete");
    assert!(Expression::new("").is_err(), "empty");
    assert!(Expression::new("element").unwrap().evaluate(&v).is_err(), "text result is not a condition");
    assert!(Expression::new("element < 3").unwrap().evaluate(&v).is_err(), "ordering element names");
}

#[test]
fn text_syntax_json_and_round_trips() {
    let p = Pipeline::parse("replicate 2 2 1 | select-expression fz > 0.5 || element == O | delete-selected").unwrap();
    assert_eq!(p.modifiers.len(), 3, "three steps (|| is not a pipe)");
    assert_eq!(p.modifiers[0].parameters["counts"].items().unwrap()[2].number().unwrap(), 1.0, "positional vector");
    assert_eq!(p.modifiers[1].parameters["expression"].string().unwrap(), "fz > 0.5 || element == O", "expression takes the rest of its step");
    let keyed = Pipeline::parse("select-slab axis=c minimum=0.25 maximum=0.75 mode=add | replicate counts=3,1,1").unwrap();
    assert!(keyed.modifiers[0].parameters["axis"].string().unwrap() == "c" && keyed.modifiers[0].parameters["mode"].string().unwrap() == "add", "key=value");
    assert_eq!(keyed.modifiers[1].parameters["counts"].items().unwrap()[0].number().unwrap(), 3.0, "comma vector");
    for misspelt in ["replicate count=2,2,2", "replicat 2 2 2"] {
        let message = Pipeline::parse(misspelt).unwrap_err();
        assert!(message.contains("did you mean"), "suggestion for {misspelt}: {message}");
    }
    let quoted = Pipeline::parse("select-expression 'element == \"O\" && x > 1' mode=intersect").unwrap();
    assert_eq!(quoted.modifiers[0].parameters["expression"].string().unwrap(), "element == \"O\" && x > 1", "quoted expression");
    let named_rest = Pipeline::parse("select-expression expression=fz > 0.5 mode=intersect").unwrap();
    assert!(
        named_rest.modifiers[0].parameters["expression"].string().unwrap() == "fz > 0.5"
            && named_rest.modifiers[0].parameters["mode"].string().unwrap() == "intersect",
        "named expression takes the following words"
    );
    let named_list = Pipeline::parse("select-element mode=add elements=Cu Ni").unwrap();
    let elements = named_list.modifiers[0].parameters["elements"].string().unwrap();
    assert!(elements.contains("Cu") && elements.contains("Ni") && named_list.modifiers[0].parameters["mode"].string().unwrap() == "add");
    for original in [&p, &keyed, &quoted] {
        let again = Pipeline::parse(&original.to_text().unwrap()).unwrap();
        assert_eq!(again.to_json().dump(-1), original.to_json().dump(-1), "text round trip of {}", original.to_text().unwrap());
        assert_eq!(Pipeline::from_json(&original.to_json()).unwrap().to_json().dump(-1), original.to_json().dump(-1), "JSON round trip");
    }
    assert!(Pipeline::parse("").unwrap().modifiers.is_empty(), "empty pipeline");
    for bad in ["no-such-modifier", "replicate 2 2", "replicate 2.5 1 1", "select-slab axis=w", "wrap | | wrap", "clear-selection extra"] {
        assert!(Pipeline::parse(bad).is_err(), "{bad} should fail");
    }
    assert!(Pipeline::from_json(&Json::parse(r#"{"modifiers": [{"type": "replicate", "parameters": {"bogus": 1}}]}"#).unwrap()).is_err());
    assert!(find_modifier_type("replicate").is_some() && find_modifier_type("nope").is_none(), "type lookup");
    for t in modifier_types() {
        let m = make_modifier(&t.id).unwrap(); // every default parses
        assert_eq!(m.parameters.size(), t.parameters.len(), "defaults of {}", t.id);
    }
}

#[test]
fn the_input_is_never_modified_disabled_steps_and_errors() {
    let input = copper(2);
    let mut p = Pipeline::parse("select-element Cu | delete-selected").unwrap();
    let (out, stages) = p.evaluate(&input).unwrap();
    assert!(out.structure.atoms.is_empty() && input.atoms.len() == 32, "input untouched");
    assert!(stages[0].selected == 32 && stages[1].atoms == 0, "stage results");
    p.modifiers[1].enabled = false;
    let (kept, stages) = p.evaluate(&input).unwrap();
    assert!(kept.structure.atoms.len() == 32 && kept.selected_count() == 32 && stages[1].skipped, "disabled step skipped");
    let broken = Pipeline::parse("replicate 2 1 1 | select-property 0 1 | delete-selected").unwrap();
    let (partial, stages) = broken.evaluate(&input).unwrap();
    assert!(!stages[1].error.is_empty() && stages[2].skipped, "error stops later steps");
    assert_eq!(partial.structure.atoms.len(), 64, "output before the failing step");
}

#[test]
fn per_stage_cache_matches_full_evaluation_after_edits_and_moves() {
    let input = copper(2);
    let mut p = Pipeline::parse("replicate 2 1 1 | select-random 0.25 seed=3 | random-displacement 0.1 | delete-selected").unwrap();
    let mut cache = PipelineCache::default();
    p.evaluate_cached(&input, &mut cache).unwrap();
    p.modifiers[2].parameters.set("amplitude", 0.2);
    cache.invalidate_from(2);
    let cached = p.evaluate_cached(&input, &mut cache).unwrap();
    let (fresh, _) = p.evaluate(&input).unwrap();
    assert_eq!(cached.structure.atoms.len(), fresh.structure.atoms.len(), "same atoms");
    for (a, b) in cached.structure.atoms.iter().zip(&fresh.structure.atoms) {
        assert_eq!(a.x, b.x, "same positions");
    }
    p.modifiers.swap(0, 1);
    cache.invalidate_from(0);
    assert_eq!(p.evaluate_cached(&input, &mut cache).unwrap().structure.atoms.len(), p.evaluate(&input).unwrap().0.structure.atoms.len(), "after a move");
    assert_eq!(cache.stages.len(), 4, "one stage result per step");
}

#[test]
fn selection_modifiers() {
    let fcc = copper(2);
    let count = |text: &str, s: &Structure| ok(text, s).selected_count();
    assert_eq!(count("select-expression fz < 0.5", &fcc), 16, "half the cell by expression");
    assert_eq!(count("select-slab z 0 3", &fcc), 16, "Cartesian slab");
    assert_eq!(count("select-slab c 0.49 1", &fcc), 16, "fractional slab");
    assert_eq!(count("select-sphere 0 0 0 2.6", &fcc), 13, "corner atom and its 12 neighbours");
    assert_eq!(count("select-sphere center=0,0,0 radius=0.1 | expand-selection 2.6", &fcc), 13, "expand to neighbours");
    assert_eq!(count("select-sphere 0 0 0 0.1 | expand-selection 2.6 steps=2", &fcc), 1 + 12 + 6 + 24 - 12, "second shell growth");
    let random = ok("select-random 0.25 seed=7", &fcc);
    assert!(random.selected_count() == 8 && ok("select-random 0.25 seed=7", &fcc).selected == random.selected, "exact, reproducible random fraction");
    assert_ne!(ok("select-random 0.25 seed=7", &fcc).selected, ok("select-random 0.25 seed=8", &fcc).selected, "seed changes the selection");
    assert_eq!(count("select-expression fz < 0.5 | select-expression fx < 0.5 mode=add", &fcc), 24, "add");
    assert_eq!(count("select-expression fz < 0.5 | select-expression fx < 0.5 mode=subtract", &fcc), 8, "subtract");
    assert_eq!(count("select-expression fz < 0.5 | select-expression fx < 0.5 mode=intersect", &fcc), 8, "intersect");
    assert_eq!(count("select-expression fz < 0.5 | invert-selection", &fcc), 16, "invert");
    assert_eq!(count("select-element Cu | clear-selection", &fcc), 0, "clear");
    let mut oxide = fcc.clone();
    oxide.atoms[0].symbol = "O".into();
    oxide.atoms[0].atomic_number = 8;
    oxide.atoms[5].symbol = "Ni".into();
    oxide.atoms[5].atomic_number = 28;
    assert!(count("select-element O Ni", &oxide) == 2 && count("select-element elements=O,Ni", &oxide) == 2, "element lists");
    assert!(run("select-element Xx", &oxide).is_err(), "unknown element");
    assert!(run("select-slab c 0.6 0.2", &fcc).is_err(), "inverted slab");
}

#[test]
fn modification_modifiers() {
    let fcc = copper(2);
    let deleted = ok("select-expression index < 4 | delete-selected", &fcc);
    assert!(deleted.structure.atoms.len() == 28 && deleted.selected_count() == 0, "delete");
    assert_eq!(deleted.structure.atoms[0].x, fcc.atoms[4].x, "remaining atoms keep their order");
    let alloy = ok("select-random 0.5 seed=1 | assign-element Ni", &fcc);
    assert_eq!(alloy.structure.atoms.iter().filter(|a| a.symbol == "Ni" && a.atomic_number == 28).count(), 16, "assign element to the selection");
    assert!(run("assign-element Ni", &fcc).is_err(), "assign without a selection");
    assert_eq!(ok("assign-element Ni target=all", &fcc).structure.atoms[31].symbol, "Ni", "assign to all");
    let moved = ok("select-expression index == 0 | displace 1 2 3 target=selected", &fcc);
    close(moved.structure.atoms[0].z, 3.0, 1e-12, "displace selected");
    assert_eq!(moved.structure.atoms[1].z, fcc.atoms[1].z, "unselected atoms stay");
    let noisy = ok("replicate 4 4 4 | random-displacement 0.1 seed=5", &fcc);
    let big = ok("replicate 4 4 4", &fcc);
    let sum2: f64 = noisy.structure.atoms.iter().zip(&big.structure.atoms).map(|(a, b)| (a.x - b.x).powi(2)).sum();
    close((sum2 / big.structure.atoms.len() as f64).sqrt(), 0.1, 0.005, "displacement standard deviation");
    assert_eq!(ok("replicate 4 4 4 | random-displacement 0.1 seed=5", &fcc).structure.atoms[7].y, noisy.structure.atoms[7].y, "reproducible");
    assert_eq!(ok("slice 0 0 1 3.0", &fcc).structure.atoms.len(), 16, "half-space above z = 3 removed");
    assert_eq!(ok("slice normal=0,0,1 distance=3 invert=true", &fcc).structure.atoms.len(), 16, "inverted half-space");
    assert_eq!(ok("slice 0 0 1 3.615 width=1", &fcc).structure.atoms.len(), 24, "slab of one (001) plane removed");
    let sliced = ok("slice 0 0 1 3.0 action=select", &fcc);
    assert!(sliced.structure.atoms.len() == 32 && sliced.selected_count() == 16, "slice can select instead");
}

#[test]
fn cell_modifiers() {
    let mut fcc = copper(1);
    fcc.atom_property = vec![1.0, 2.0, 3.0, 4.0];
    fcc.atom_property_name = "test".into();
    let big = ok("replicate 2 3 1", &fcc);
    assert_eq!(big.structure.atoms.len(), 24, "replicated atoms");
    close(big.structure.cell_vectors[1][1], 3.0 * 3.615, 1e-12, "replicated cell");
    assert!(big.structure.atom_property.len() == 24 && big.structure.atom_property[5] == 2.0, "per-atom property follows the copies");
    let mut open = fcc.clone();
    open.has_unit_cell = false;
    assert!(run("replicate 2 2 2", &open).is_err(), "replicate needs a cell");
    let doubled = ok("transform 2 0 0 0 2 0 0 0 2", &fcc);
    close(af_science::core::cell_volume(&doubled.structure.cell_vectors), 8.0 * 3.615f64.powi(3), 1e-9, "transform scales the cell");
    close(doubled.structure.atoms[1].y, 3.615, 1e-12, "and the atoms");
    assert!(run("transform 1 0 0 1 0 0 0 0 1", &fcc).is_err(), "singular matrix");
    close(ok("strain 0.01 0 0", &fcc).structure.cell_vectors[0][0], 3.615 * 1.01, 1e-12, "normal strain");
    close(ok("strain shear=0,0,0.02", &fcc).structure.cell_vectors[0][1], 3.615 * 0.01, 1e-12, "engineering shear is twice the tensor component");
    let mut outside = fcc.clone();
    outside.atoms[0].x = -0.5;
    outside.atoms[1].z = 3.615 + 0.25;
    let wrapped = ok("wrap", &outside);
    close(wrapped.structure.atoms[0].x, 3.115, 1e-12, "wrapped below");
    close(wrapped.structure.atoms[1].z, 0.25, 1e-12, "wrapped above");
    let slab = ok("add-vacuum c 10", &fcc);
    close(slab.structure.cell_vectors[2][2], 13.615, 1e-12, "vacuum lengthens c");
    close(slab.structure.atoms[0].z, 5.0, 1e-12, "atoms moved to the middle");
    let centred = ok("add-vacuum c 10 center=false | center", &fcc);
    let (top, middle) = (0.5 * 3.615, 0.5 * 13.615);
    close(centred.structure.atoms[0].z, middle - top / 2.0, 1e-9, "center puts the atoms' extent in the middle");
    close(centred.structure.atoms[1].z, middle + top / 2.0, 1e-9, "both ends of the extent");
}

#[test]
fn edit_modifiers_add_atom_set_cell_merge() {
    let fcc = copper(1);
    let added = ok("add-atom O 0.5 0.5 0.5 fractional=true", &fcc);
    assert!(added.structure.atoms.len() == 5 && added.structure.atoms[4].symbol == "O" && added.selected_count() == 1, "atom added and selected");
    close(added.structure.atoms[4].x, 1.8075, 1e-12, "fractional position");
    close(ok("set-cell 4 0 0 0 4 0 0 0 4", &fcc).structure.atoms[1].y, 2.0, 1e-12, "atoms scale with the cell");
    close(ok("set-cell a=4,0,0 b=0,4,0 c=0,0,4 scale_atoms=false", &fcc).structure.atoms[1].y, 1.8075, 1e-12, "atoms stay without scaling");
    assert!(run("set-cell 1 0 0 2 0 0 0 0 1", &fcc).is_err(), "coplanar cell");
    let file = std::env::temp_dir().join(format!("atomforge_merge_{}.extxyz", std::process::id()));
    std::fs::write(&file, "2\nmerge test\nNi 0 0 0\nNi 1 0 0\n").unwrap();
    let merged = ok(&format!("merge file={} offset=10,0,0", file.display()), &fcc);
    assert!(merged.structure.atoms.len() == 6 && merged.structure.atoms[5].symbol == "Ni" && merged.selected_count() == 2, "merged and selected");
    close(merged.structure.atoms[5].x, 11.0, 1e-12, "offset");
    std::fs::remove_file(&file).unwrap();
    assert!(run("merge", &fcc).is_err(), "merge needs a file");
}

#[test]
fn supercell_transform_structure() {
    let fcc = copper(1);
    let doubled = ok("supercell \"2 0 0 0 2 0 0 0 2\"", &fcc);
    assert_eq!(doubled.structure.atoms.len(), 32, "2x2x2 supercell");
    close(doubled.structure.cell_vectors[0][0], 2.0 * 3.615, 1e-9, "doubled a");
    let rotated = ok("supercell matrix=1,1,0,-1,1,0,0,0,1", &fcc);
    assert_eq!(rotated.structure.atoms.len(), 8, "rotated supercell");
    close(rotated.structure.cell_vectors[1][0], -3.615, 1e-9, "second vector is -a+b");
    let shrunk = ok("supercell \"1 0 0 0 1 0 0 0 1\" | select-random 0.5 seed=2", &fcc);
    assert!(shrunk.structure.atoms.len() == 4 && shrunk.selected_count() == 2, "identity keeps the cell");
    assert!(run("supercell \"1 0 0 1 0 0 0 0 1\"", &fcc).is_err(), "singular matrix");
    assert!(run("supercell \"1.5 0 0 0 1 0 0 0 1\"", &fcc).is_err(), "non-integer matrix");
}

#[test]
fn insert_interstitials_fills_fcc_voids() {
    let fcc = copper(2);
    let octahedral = ok("insert-interstitials H kind=octahedral count=0", &fcc);
    let tetrahedral = ok("insert-interstitials H kind=tetrahedral count=0", &fcc);
    // An fcc cell has as many octahedral sites as atoms and twice as many tetrahedral sites.
    assert_eq!(octahedral.structure.atoms.len(), 32 + 32, "octahedral sites");
    assert_eq!(tetrahedral.structure.atoms.len(), 32 + 64, "tetrahedral sites");
    let three = ok("insert-interstitials C count=3 order=random seed=4", &fcc);
    assert!(three.structure.atoms.len() == 35 && three.selected_count() == 3 && three.structure.atoms.last().unwrap().symbol == "C", "three new, selected atoms");
    assert_eq!(ok("insert-interstitials C count=3 order=random seed=4", &fcc).structure.atoms.last().unwrap().x, three.structure.atoms.last().unwrap().x, "reproducible");
}
