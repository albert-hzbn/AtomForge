//! The science core and data readers against the C++ reference outputs
//! (rust/tests/golden/science.json, written by tools/golden/golden_science.cpp).

use std::path::{Path, PathBuf};

use af_core::json::Json;
use af_science::core::*;
use af_science::data::*;
use num_complex::Complex64;

fn num(v: f64) -> Json {
    Json::Number(v)
}

fn attempt(f: impl FnOnce() -> Result<Json, String>) -> Json {
    f().unwrap_or_else(|e| {
        let mut error = Json::object();
        error.set("error", Json::String(e));
        error
    })
}

fn bools(v: &[bool]) -> Json {
    Json::Array(v.iter().map(|&b| Json::Bool(b)).collect())
}

fn numbers(v: &[f64]) -> Json {
    Json::Array(v.iter().map(|&x| num(x)).collect())
}

fn vec3_of(j: &Json) -> Vec3 {
    let items = j.items().unwrap();
    [items[0].number().unwrap(), items[1].number().unwrap(), items[2].number().unwrap()]
}

fn mat3_of(j: &Json) -> Mat3 {
    let rows = j.items().unwrap();
    [vec3_of(&rows[0]), vec3_of(&rows[1]), vec3_of(&rows[2])]
}

fn frame_json(frame: &FrameData) -> Json {
    let mut result = Json::object();
    result.set("structure", structure_json(&frame.structure));
    result.set("has_cell", Json::Bool(frame.structure.has_unit_cell));
    result.set("pbc", bools(&frame.pbc));
    result.set("velocities", Json::Array(frame.velocities.iter().map(|v| vec3_json(*v)).collect()));
    result.set("masses", numbers(&frame.masses));
    result.set("time_fs", frame.time_fs.map_or(Json::Null, num));
    result.set("atom_property", numbers(&frame.structure.atom_property));
    result.set(
        "colours",
        Json::Array(frame.structure.atoms.iter().map(|a| numbers(&[a.r as f64, a.g as f64, a.b as f64])).collect()),
    );
    result
}

fn frames_json(frames: &[FrameData]) -> Json {
    Json::Array(frames.iter().map(frame_json).collect())
}

fn describe_parameters(p: &Parameters, names: &[&str]) -> Json {
    let mut r = Json::object();
    for &n in names {
        let mut item = Json::object();
        item.set("has", Json::Bool(p.has(n)));
        item.set("array", attempt(|| Ok(array_json(&p.array(n)?))));
        item.set("json", attempt(|| Ok(p.json(n)?.clone())));
        item.set(
            "structure",
            attempt(|| {
                let s = p.structure(n)?;
                let mut o = Json::object();
                o.set("structure", structure_json(&s.structure));
                o.set("pbc", bools(&s.pbc));
                Ok(o)
            }),
        );
        item.set("frames", attempt(|| Ok(frames_json(p.frames(n)?))));
        r.set(n, item);
    }
    r.set("pbc_fallback", bools(&p.pbc("pbc", [false, true, false]).unwrap()));
    r
}

fn parameter_cases(base: &Path, name: &str) -> Json {
    let requests = [
        ("x", r#"{"data": {"file": "@"}}"#),
        ("x", r#"{"data": {"file": "@", "column": 1}}"#),
        ("x", r#"{"data": {"file": "@", "column": 0}}"#),
        ("x", r#"{"data": {"file": "@", "field": "positions"}}"#),
        ("x", r#"{"data": {"file": "@", "field": "cell.1"}}"#),
        ("x", r#"{"energies_eV": {"file": "@"}}"#),
        ("x", r#"{"velocities": {"file": "@"}, "masses": {"file": "@"}, "cell": {"file": "@"}}"#),
        ("msd", r#"{"positions": {"file": "@"}, "timestep_fs": 2.5}"#),
        ("x", r#"{"structure": {"file": "@"}, "initial": {"file": "./@"}}"#),
        ("x", r#"{"restart_images": {"file": "@"}, "out_file": "sub/../@", "other": {"file": "@", "bogus": 1}}"#),
    ];
    let names = ["data", "energies_eV", "velocities", "masses", "cell", "positions", "structure", "initial", "restart_images", "out_file"];
    let mut cases = Json::Array(Vec::new());
    for (tool, pattern) in requests {
        let text = pattern.replace('@', name);
        let mut item = Json::object();
        item.set("tool", Json::String(tool.into()));
        item.set("request", Json::String(text.clone()));
        item.set(
            "result",
            attempt(|| {
                let p = Parameters::new(tool, &Json::parse(&text)?, base, None)?;
                Ok(describe_parameters(&p, &names))
            }),
        );
        cases.push(item);
    }
    cases
}

/// First difference between two JSON values (exact numbers).
fn diff(path: &str, a: &Json, b: &Json, out: &mut Vec<String>) {
    if out.len() > 30 {
        return;
    }
    match (a, b) {
        (Json::Object(x), Json::Object(y)) => {
            let kx: Vec<_> = x.iter().map(|(k, _)| k).collect();
            let ky: Vec<_> = y.iter().map(|(k, _)| k).collect();
            if kx != ky {
                out.push(format!("{path}: keys {kx:?} vs {ky:?}"));
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
        (Json::Number(x), Json::Number(y)) if x.to_bits() == y.to_bits() || (x == y) => {}
        _ if a == b => {}
        _ => out.push(format!("{path}: {} vs {}", a.dump(-1), b.dump(-1))),
    }
}

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn science_core_and_data_match_cpp() {
    // Paths in messages are relative to the repository root, as for the dumper.
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    std::env::set_current_dir(&root).unwrap();
    let golden = Json::parse(&std::fs::read_to_string("rust/tests/golden/science.json").unwrap()).unwrap();
    let mut differences = Vec::new();

    // Algebra, minimum image and neighbour lists.
    for (c, expected) in golden.at("algebra").unwrap().items().unwrap().iter().enumerate() {
        let cell = mat3_of(expected.at("cell").unwrap());
        let pbc_items = expected.at("pbc").unwrap().items().unwrap();
        let pbc = [pbc_items[0].boolean().unwrap(), pbc_items[1].boolean().unwrap(), pbc_items[2].boolean().unwrap()];
        let mut item = Json::object();
        item.set("cell", mat3_json(&cell));
        item.set("pbc", bools(&pbc));
        item.set("det", num(determinant(&cell)));
        item.set("inverse", mat3_json(&inverse(&cell).unwrap()));
        item.set("product", mat3_json(&multiply(&cell, &transpose(&cell))));
        let image = MinimumImage::new(&cell, pbc).unwrap();
        let vectors: Vec<Vec3> = expected.at("vectors").unwrap().items().unwrap().iter().map(vec3_of).collect();
        item.set("vectors", Json::Array(vectors.iter().map(|v| vec3_json(*v)).collect()));
        item.set("images", Json::Array(vectors.iter().map(|v| vec3_json(image.apply(*v))).collect()));
        item.set("fractional", Json::Array(vectors.iter().map(|v| vec3_json(fractional(*v, &cell).unwrap())).collect()));
        let positions: Vec<Vec3> = expected.at("positions").unwrap().items().unwrap().iter().map(vec3_of).collect();
        item.set("positions", Json::Array(positions.iter().map(|v| vec3_json(*v)).collect()));
        let cutoff = expected.at("cutoff").unwrap().number().unwrap();
        item.set("cutoff", num(cutoff));
        item.set(
            "neighbors",
            attempt(|| {
                Ok(Json::Array(
                    neighbor_list(&positions, &cell, pbc, cutoff)?
                        .iter()
                        .map(|n| {
                            numbers(&[
                                n.i as f64,
                                n.j as f64,
                                n.shift[0] as f64,
                                n.shift[1] as f64,
                                n.shift[2] as f64,
                                n.vector[0],
                                n.vector[1],
                                n.vector[2],
                            ])
                        })
                        .collect(),
                ))
            }),
        );
        diff(&format!("algebra[{c}]"), &item, expected, &mut differences);
    }

    // FFTs.
    for (c, expected) in golden.at("fft").unwrap().items().unwrap().iter().enumerate() {
        let input: Vec<Complex64> = expected
            .at("input")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|p| Complex64::new(p.items().unwrap()[0].number().unwrap(), p.items().unwrap()[1].number().unwrap()))
            .collect();
        let mut item = Json::object();
        item.set("input", Json::Array(input.iter().map(|d| numbers(&[d.re, d.im])).collect()));
        for inverse in [false, true] {
            let out = fft(input.clone(), inverse);
            item.set(if inverse { "inverse" } else { "forward" }, Json::Array(out.iter().map(|d| numbers(&[d.re, d.im])).collect()));
        }
        diff(&format!("fft[{c}]"), &item, expected, &mut differences);
    }

    // Arrays and validation helpers.
    for (c, expected) in golden.at("arrays").unwrap().items().unwrap().iter().enumerate() {
        let text = expected.at("text").unwrap().string().unwrap().to_string();
        let mut item = Json::object();
        item.set("text", Json::String(text.clone()));
        if expected.contains("array") {
            item.set(
                "array",
                attempt(|| {
                    let a = to_array(&Json::parse(&text)?, "x")?;
                    let mut r = Json::object();
                    r.set("shape", numbers(&a.shape.iter().map(|&s| s as f64).collect::<Vec<_>>()));
                    r.set("values", numbers(&a.values));
                    r.set("json", array_json(&a));
                    Ok(r)
                }),
            );
        }
        item.set("cell", attempt(|| Ok(mat3_json(&cell_matrix(to_array(&Json::parse(&text)?, "cell")?)?))));
        if expected.contains("pbc") {
            item.set("pbc", attempt(|| Ok(bools(&parse_pbc(&Json::parse(&text)?)?))));
        }
        diff(&format!("arrays[{c}]"), &item, expected, &mut differences);
    }
    for (c, expected) in golden.at("checks").unwrap().items().unwrap().iter().enumerate() {
        let v = expected.at("value").unwrap().number().unwrap();
        let mut item = Json::object();
        item.set("value", num(v));
        item.set("positive", attempt(|| Ok(num(positive(v, "x")?))));
        item.set("integer", attempt(|| Ok(num(integer(v, "n", 1)? as f64))));
        item.set("integer0", attempt(|| Ok(num(integer(v, "n", 0)? as f64))));
        diff(&format!("checks[{c}]"), &item, expected, &mut differences);
    }

    // Frame readers, streams and parameters.
    let fixtures = Path::new("rust/tests/fixtures/science");
    let files = golden.at("files").unwrap();
    let mut names: Vec<String> = std::fs::read_dir(fixtures).unwrap().map(|e| e.unwrap().file_name().to_string_lossy().into_owned()).collect();
    names.sort();
    let expected_names: Vec<&String> = files.members().unwrap().iter().map(|(k, _)| k).collect();
    assert_eq!(names.iter().collect::<Vec<_>>(), expected_names, "fixture files differ from the golden run");
    for name in &names {
        let path = fixtures.join(name);
        let mut item = Json::object();
        item.set("frames", attempt(|| Ok(frames_json(&read_frames(&path, None)?))));
        item.set(
            "stream",
            attempt(|| {
                let stream = TrajectoryStream::new(&path, None)?;
                let mut r = Json::object();
                r.set("streamed", Json::Bool(stream.streamed()));
                let frames = (0..stream.len()).map(|i| stream.frame(i)).collect::<Result<Vec<_>, _>>()?;
                r.set("frames", frames_json(&frames));
                Ok(r)
            }),
        );
        item.set("parameters", parameter_cases(fixtures, name));
        diff(&format!("files.{name}"), &item, files.at(name).unwrap(), &mut differences);
    }

    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}
