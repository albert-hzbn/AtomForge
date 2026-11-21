//! The science tools against the C++ ones: tests/fixtures/tools/requests.json
//! run through run_tool and compared with tests/golden/tools.json (from
//! tools/golden/golden_tools.cpp). Tools built on Eigen decompositions are
//! compared to 1e-9 relative plus 1e-10 absolute (docs/rust-port-plan.md 6a), eigenvector
//! columns up to sign; phonons, whose Gamma acoustic modes are square roots of
//! rounding noise that also set the DOS grid, to 1e-5 relative plus 1e-6
//! absolute; the rest exactly, reports included.

use af_core::json::Json;

// Tools whose results pass through an Eigen least-squares fit or a dense
// (complex) eigen-solver, compared to 1e-9 relative (docs/rust-port-plan.md 6a).
const EIGEN_TOOLS: [&str; 11] = [
    "local-strain",
    "effective-mass",
    "equation-of-state",
    "elastic-tensor",
    "ptm",
    "grain-segmentation",
    "elastic-strain",
    "fcc-planar-faults",
    "eos-scan",
    "elastic-constants",
    "phonons",
];

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn read_json(path: &std::path::Path) -> Json {
    Json::parse(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// (relative, absolute) tolerance; None: exact.
type Tolerance = Option<(f64, f64)>;

fn tolerance_of(tool: &str) -> Tolerance {
    if tool == "phonons" {
        Some((1e-5, 1e-6))
    } else if EIGEN_TOOLS.contains(&tool) {
        Some((1e-9, 1e-10))
    } else {
        None
    }
}

fn same(a: f64, b: f64, tolerance: Tolerance) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    match tolerance {
        Some((relative, absolute)) => a == b || (a - b).abs() <= relative * a.abs().max(b.abs()) + absolute,
        None => a.to_bits() == b.to_bits() || a == b,
    }
}

fn diff(path: &str, a: &Json, b: &Json, relative: Tolerance, out: &mut Vec<String>) {
    match (a, b) {
        // C++ Json keeps NaN as a number and writes it as null.
        (Json::Number(x), Json::Null) | (Json::Null, Json::Number(x)) if x.is_nan() => {}
        (Json::Number(x), Json::Number(y)) => {
            if !same(*x, *y, relative) {
                out.push(format!("{path}: {x:e} vs {y:e}"));
            }
        }
        (Json::Array(x), Json::Array(y)) => {
            if x.len() != y.len() {
                out.push(format!("{path}: length {} vs {}", x.len(), y.len()));
                return;
            }
            for (i, (p, q)) in x.iter().zip(y).enumerate() {
                diff(&format!("{path}[{i}]"), p, q, relative, out);
            }
        }
        (Json::Object(x), Json::Object(y)) => {
            let kx: Vec<&String> = x.iter().map(|(k, _)| k).collect();
            let ky: Vec<&String> = y.iter().map(|(k, _)| k).collect();
            if kx != ky {
                out.push(format!("{path}: keys {kx:?} vs {ky:?}"));
                return;
            }
            for ((k, p), (_, q)) in x.iter().zip(y) {
                diff(&format!("{path}.{k}"), p, q, relative, out);
            }
        }
        _ => {
            if a != b {
                out.push(format!("{path}: {} vs {}", a.dump(-1), b.dump(-1)));
            }
        }
    }
}

/// Eigenvector columns with their largest component positive.
fn canonical_columns(matrix: &mut Json) {
    let Json::Array(rows) = matrix else { return };
    let n = rows.len();
    for c in 0..n {
        let column: Vec<f64> = rows.iter().map(|r| r.items().unwrap()[c].number().unwrap()).collect();
        let largest = column.iter().cloned().fold(0.0f64, |m, v| if v.abs() > m.abs() { v } else { m });
        if largest < 0.0 {
            for r in rows.iter_mut() {
                if let Json::Array(items) = r {
                    if let Json::Number(v) = &mut items[c] {
                        *v = -*v;
                    }
                }
            }
        }
    }
}

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn science_tools_match_cpp() {
    let root = repo_root();
    std::env::set_current_dir(&root).unwrap();
    let requests = read_json(&root.join("rust/tests/fixtures/tools/requests.json"));
    let golden = read_json(&root.join("rust/tests/golden/tools.json"));
    let mut differences = Vec::new();
    let mut compared = 0;
    for (k, (case, expected)) in requests.items().unwrap().iter().zip(golden.items().unwrap()).enumerate() {
        let tool = case.at("tool").unwrap().string().unwrap();
        let path = format!("[{k}] {tool}");
        let relative = tolerance_of(tool);
        match af_science::tools::run_tool(tool, case.at("request").unwrap(), std::path::Path::new("."), None) {
            Ok(output) => {
                let Some(want) = expected.find("result") else {
                    differences.push(format!("{path}: Rust ran, C++ failed with {}", expected.at("error").unwrap().dump(-1)));
                    continue;
                };
                let (mut got, mut want) = (output.result.clone(), want.clone());
                if tool == "effective-mass" {
                    for j in [&mut got, &mut want] {
                        if let Json::Object(members) = j {
                            for (key, value) in members.iter_mut() {
                                if key == "principal_axes_columns" {
                                    canonical_columns(value);
                                }
                            }
                        }
                    }
                }
                diff(&path, &got, &want, relative, &mut differences);
                if relative.is_none() {
                    let report = af_science::tools::result_report(tool, &output.result).unwrap();
                    if report != expected.at("report").unwrap().string().unwrap() {
                        differences.push(format!("{path}: report differs:\n{report}\nvs\n{}", expected.at("report").unwrap().string().unwrap()));
                    }
                }
                let frames = Json::Array(output.frames.iter().map(af_science::data::structure_json).collect());
                diff(&format!("{path}.frames"), &frames, expected.at("frames").unwrap(), relative, &mut differences);
                let times = Json::Array(output.times.iter().map(|&t| Json::Number(t)).collect());
                diff(&format!("{path}.times"), &times, expected.at("times").unwrap(), relative, &mut differences);
                compared += 1;
            }
            Err(error) => {
                let want = expected.find("error").map(|e| e.string().unwrap().to_string()).unwrap_or_else(|| "(C++ succeeded)".into());
                if error != want {
                    differences.push(format!("{path}: error {error:?} vs {want:?}"));
                }
            }
        }
    }
    println!("{compared} results compared");
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}
