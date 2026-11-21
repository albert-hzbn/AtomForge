//! The pipeline against the C++ reference outputs (rust/tests/golden/pipeline.json,
//! written by tools/golden/golden_pipeline.cpp): text and JSON forms,
//! expression values, and outputs and stage results of pipelines.

use std::path::PathBuf;

use af_core::json::Json;
use af_core::{AtomSite, Structure};
use af_pipeline::support::evaluate_for_atoms;
use af_pipeline::*;

fn attempt(f: impl FnOnce() -> Result<Json, String>) -> Json {
    f().unwrap_or_else(|e| {
        let mut error = Json::object();
        error.set("error", e);
        error
    })
}

fn vec3(v: [f64; 3]) -> Json {
    Json::Array(v.iter().map(|&x| Json::Number(x)).collect())
}

fn maybe_nan(v: f64) -> Json {
    if v.is_nan() {
        Json::String("nan".into())
    } else {
        Json::Number(v)
    }
}

fn structure_dump(s: &Structure) -> Json {
    let mut o = Json::object();
    let mut atoms = Vec::new();
    for a in &s.atoms {
        let mut atom = Json::object();
        atom.set("symbol", a.symbol.as_str());
        atom.set("z", a.atomic_number);
        atom.set("position", vec3([a.x, a.y, a.z]));
        atom.set("colour", vec3([a.r as f64, a.g as f64, a.b as f64]));
        atom.set("fixed", a.fixed as i32);
        atom.set("magmom", a.magmom);
        atom.set("charge", a.charge);
        atom.set("velocity", vec3(a.velocity));
        atoms.push(atom);
    }
    o.set("atoms", Json::Array(atoms));
    o.set("has_cell", s.has_unit_cell);
    o.set("cell", Json::Array(s.cell_vectors.iter().map(|r| vec3(*r)).collect()));
    o.set("pbc", Json::Array(s.pbc.iter().map(|&b| Json::Bool(b)).collect()));
    o.set("offset", vec3(s.cell_offset));
    o.set("grain_colors", Json::Array(s.grain_colors.iter().map(|c| vec3([c[0] as f64, c[1] as f64, c[2] as f64])).collect()));
    o.set("grain_region_ids", Json::Array(s.grain_region_ids.iter().map(|&r| Json::Number(r as f64)).collect()));
    o.set("atom_property", Json::Array(s.atom_property.iter().map(|&v| maybe_nan(v)).collect()));
    o.set("frozen_property", Json::Array(s.frozen_property.iter().map(|&v| maybe_nan(v)).collect()));
    o.set("atom_property_name", s.atom_property_name.as_str());
    o.set("frozen_property_name", s.frozen_property_name.as_str());
    o.set("loop_points", s.dislocation_loop_points.len());
    o.set("detection_done", s.dislocation_detection_done);
    o
}

fn numbers(j: &Json) -> Vec<f64> {
    j.items()
        .unwrap()
        .iter()
        .map(|v| match v {
            Json::String(_) => f64::NAN,
            other => other.number().unwrap(),
        })
        .collect()
}

fn v3(j: &Json) -> [f64; 3] {
    let n = numbers(j);
    [n[0], n[1], n[2]]
}

fn structure_from_dump(o: &Json) -> Structure {
    let mut s = Structure::default();
    for a in o.at("atoms").unwrap().items().unwrap() {
        let colour = v3(a.at("colour").unwrap());
        let position = v3(a.at("position").unwrap());
        s.atoms.push(AtomSite {
            symbol: a.at("symbol").unwrap().string().unwrap().into(),
            atomic_number: a.at("z").unwrap().number().unwrap() as i32,
            x: position[0],
            y: position[1],
            z: position[2],
            r: colour[0] as f32,
            g: colour[1] as f32,
            b: colour[2] as f32,
            fixed: a.at("fixed").unwrap().number().unwrap() as u8,
            magmom: a.at("magmom").unwrap().number().unwrap(),
            charge: a.at("charge").unwrap().number().unwrap(),
            velocity: v3(a.at("velocity").unwrap()),
        });
    }
    s.has_unit_cell = o.at("has_cell").unwrap().boolean().unwrap();
    let cell = o.at("cell").unwrap().items().unwrap();
    s.cell_vectors = [v3(&cell[0]), v3(&cell[1]), v3(&cell[2])];
    let pbc = o.at("pbc").unwrap().items().unwrap();
    s.pbc = [pbc[0].boolean().unwrap(), pbc[1].boolean().unwrap(), pbc[2].boolean().unwrap()];
    s.cell_offset = v3(o.at("offset").unwrap());
    s.grain_colors = o.at("grain_colors").unwrap().items().unwrap().iter().map(|c| v3(c).map(|x| x as f32)).collect();
    s.grain_region_ids = numbers(o.at("grain_region_ids").unwrap()).iter().map(|&r| r as i32).collect();
    s.atom_property = numbers(o.at("atom_property").unwrap());
    s.frozen_property = numbers(o.at("frozen_property").unwrap());
    s.atom_property_name = o.at("atom_property_name").unwrap().string().unwrap().to_string();
    s.frozen_property_name = o.at("frozen_property_name").unwrap().string().unwrap().to_string();
    s
}

fn diff(path: &str, a: &Json, b: &Json, out: &mut Vec<String>) {
    if out.len() > 40 {
        return;
    }
    match (a, b) {
        (Json::Object(x), Json::Object(y)) => {
            let kx: Vec<_> = x.iter().map(|(k, _)| k).collect();
            let ky: Vec<_> = y.iter().map(|(k, _)| k).collect();
            if kx != ky {
                out.push(format!("{path}: keys {kx:?} vs {ky:?} ({} vs {})", a.dump(-1), b.dump(-1)));
                return;
            }
            for ((k, va), (_, vb)) in x.iter().zip(y) {
                diff(&format!("{path}.{k}"), va, vb, out);
            }
        }
        (Json::Array(x), Json::Array(y)) => {
            if x.len() != y.len() {
                out.push(format!("{path}: length {} vs {}", x.len(), y.len()));
                return;
            }
            for (i, (va, vb)) in x.iter().zip(y).enumerate() {
                diff(&format!("{path}[{i}]"), va, vb, out);
            }
        }
        (Json::Number(x), Json::Number(y)) if x == y => {}
        _ if a == b => {}
        _ => out.push(format!("{path}: {} vs {}", a.dump(-1), b.dump(-1))),
    }
}

const CPP_MESH: &str = "rust/tests/golden/shape_";
const RUST_MESH: &str = "rust/target/golden/rust_shape_";

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn pipeline_matches_cpp() {
    // Paths (merge files) are relative to the repository root, as for the dumper.
    std::env::set_current_dir(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")).unwrap();
    let golden = Json::parse(&std::fs::read_to_string("rust/tests/golden/pipeline.json").unwrap()).unwrap();
    let mut differences = Vec::new();

    // The registry (the types ported so far, in the C++ menu order).
    let cpp_types: Vec<String> = golden.at("types").unwrap().items().unwrap().iter().map(|t| t.string().unwrap().to_string()).collect();
    let rust_types: Vec<String> = modifier_types().iter().map(|t| t.id.clone()).collect();
    let in_cpp_order: Vec<&String> = cpp_types.iter().filter(|t| rust_types.contains(t)).collect();
    assert_eq!(rust_types.iter().collect::<Vec<_>>(), in_cpp_order, "registry order differs from the C++ menu order");
    assert!(rust_types.iter().all(|t| cpp_types.contains(t)), "a Rust step is not a C++ step");

    for (i, expected) in golden.at("texts").unwrap().items().unwrap().iter().enumerate() {
        let text = expected.at("text").unwrap().string().unwrap();
        let mut item = Json::object();
        item.set("text", text);
        item.set(
            "parse",
            attempt(|| {
                let p = Pipeline::parse(text)?;
                let mut r = Json::object();
                r.set("json", p.to_json());
                r.set("text", p.to_text()?);
                r.set("round_trip", Pipeline::from_json(&Json::parse(&p.to_json().dump(-1))?)?.to_text()?);
                Ok(r)
            }),
        );
        diff(&format!("texts[{i}] {text:?}"), &item, expected, &mut differences);
    }

    for (i, expected) in golden.at("jsons").unwrap().items().unwrap().iter().enumerate() {
        let text = expected.at("text").unwrap().string().unwrap();
        let mut item = Json::object();
        item.set("text", text);
        item.set(
            "result",
            attempt(|| {
                let p = Pipeline::from_json(&Json::parse(text)?)?;
                let mut r = Json::object();
                r.set("json", p.to_json());
                r.set("text", p.to_text()?);
                Ok(r)
            }),
        );
        diff(&format!("jsons[{i}] {text}"), &item, expected, &mut differences);
    }

    let inputs: Vec<Structure> = golden.at("inputs").unwrap().items().unwrap().iter().map(structure_from_dump).collect();
    for (i, (s, expected)) in inputs.iter().zip(golden.at("inputs").unwrap().items().unwrap()).enumerate() {
        diff(&format!("inputs[{i}]"), &structure_dump(s), expected, &mut differences);
    }

    for (i, expected) in golden.at("expressions").unwrap().items().unwrap().iter().enumerate() {
        let text = expected.at("text").unwrap().string().unwrap();
        let mut item = Json::object();
        item.set("text", text);
        let values = inputs
            .iter()
            .map(|s| {
                attempt(|| {
                    let e = Expression::new(text)?;
                    let mut data = PipelineData::new(s.clone());
                    if !data.selected.is_empty() {
                        data.selected[0] = 1;
                    }
                    Ok(Json::Array(evaluate_for_atoms(&e, &data)?.into_iter().map(maybe_nan).collect()))
                })
            })
            .collect();
        item.set("values", Json::Array(values));
        diff(&format!("expressions[{i}] {text:?}"), &item, expected, &mut differences);
    }

    for (i, expected) in golden.at("runs").unwrap().items().unwrap().iter().enumerate() {
        let text = expected.at("text").unwrap().string().unwrap();
        let disabled: Vec<usize> = numbers(expected.at("disabled").unwrap()).iter().map(|&d| d as usize).collect();
        let mut item = Json::object();
        item.set("text", text);
        item.set("disabled", expected.at("disabled").unwrap().clone());
        let outputs = inputs
            .iter()
            .map(|s| {
                attempt(|| {
                    // Shape surfaces go to the target directory, to compare with the C++ files.
                    let mut p = Pipeline::parse(&text.replace(CPP_MESH, RUST_MESH))?;
                    for &d in &disabled {
                        p.modifiers[d].enabled = false;
                    }
                    let (mut out, mut stages) = p.evaluate(s)?;
                    for note in out.notes.iter_mut().chain(stages.iter_mut().flat_map(|st| st.notes.iter_mut())) {
                        *note = note.replace(RUST_MESH, CPP_MESH);
                    }
                    let mut r = Json::object();
                    r.set("structure", structure_dump(&out.structure));
                    r.set("selected", Json::Array(out.selected.iter().map(|&c| Json::Number(c as f64)).collect()));
                    r.set("notes", Json::Array(out.notes.iter().map(|n| Json::String(n.clone())).collect()));
                    let stage_list = stages
                        .iter()
                        .map(|st| {
                            let mut o = Json::object();
                            o.set("atoms", st.atoms);
                            o.set("selected", st.selected);
                            o.set("error", st.error.as_str());
                            o.set("notes", Json::Array(st.notes.iter().map(|n| Json::String(n.clone())).collect()));
                            o.set("skipped", st.skipped);
                            o
                        })
                        .collect();
                    r.set("stages", Json::Array(stage_list));
                    Ok(r)
                })
            })
            .collect();
        item.set("outputs", Json::Array(outputs));
        diff(&format!("runs[{i}] {text:?}"), &item, expected, &mut differences);
        // The OBJ surfaces written by build-shape, byte for byte.
        for word in text.split_whitespace().filter_map(|w| w.strip_prefix("mesh_file=")) {
            let cpp = std::fs::read(word).ok();
            let rust = std::fs::read(word.replace(CPP_MESH, RUST_MESH)).ok();
            if cpp != rust {
                differences.push(format!("runs[{i}]: {word} differs from the Rust surface"));
            }
        }
    }

    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn cached_evaluation_equals_a_full_one() {
    let mut s = Structure { has_unit_cell: true, cell_vectors: [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]], ..Default::default() };
    s.atoms.push(AtomSite::new("Cu", [0.0, 0.0, 0.0]));
    s.atoms.push(AtomSite::new("O", [1.5, 1.5, 1.5]));
    let mut p = Pipeline::parse("replicate 2 2 2 | select-element O | delete-selected | displace 0.1 0 0").unwrap();
    let mut cache = PipelineCache::default();
    let first = p.evaluate_cached(&s, &mut cache).unwrap();
    assert_eq!(first.structure.atoms.len(), 8);
    assert_eq!(cache.valid, 4);
    // Editing step 1 re-runs steps 1 and later only.
    p.modifiers[1].parameters.set("elements", "Cu");
    cache.invalidate_from(1);
    let second = p.evaluate_cached(&s, &mut cache).unwrap();
    let (full, stages) = p.evaluate(&s).unwrap();
    assert_eq!(second, full);
    assert_eq!(second.structure.atoms.len(), 8);
    assert!(second.structure.atoms.iter().all(|a| a.symbol == "O"));
    assert_eq!(cache.stages.iter().map(|s| (s.atoms, s.error.clone())).collect::<Vec<_>>(), stages.iter().map(|s| (s.atoms, s.error.clone())).collect::<Vec<_>>());
    // A cancelled evaluation fails as a whole.
    let control = std::sync::Arc::new(af_core::task::TaskControl::default());
    control.cancel();
    let _scope = af_core::task::Scope::new(control.clone());
    cache.clear();
    assert_eq!(p.evaluate_cached(&s, &mut cache).unwrap_err(), "Calculation cancelled");
}
