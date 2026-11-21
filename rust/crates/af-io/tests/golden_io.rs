//! af-io against the C++ loader and writers (tests/golden/io.json, written by
//! tools/golden/golden_io.cpp with Open Babel and spglib).

use std::path::{Path, PathBuf};

use af_core::Structure;
use serde_json::Value;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn golden() -> Value {
    let text = std::fs::read_to_string(root().join("tests/golden/io.json")).expect("io.json");
    serde_json::from_str(&text).expect("JSON")
}

fn fixture_path(relative: &str) -> String {
    let (dir, name) = relative.split_once('/').unwrap();
    let base = if dir == "fixtures" { root().join("tests/fixtures/io") } else { root().join("../docs/manual/examples") };
    base.join(name).to_string_lossy().into_owned()
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

/// Loads where the port is intentionally better than the C++ build.
const IMPROVED: [&str; 1] = [
    // Open Babel only writes Gaussian input; the port reads it.
    "fixtures/ethanol.gjf",
];

fn compare(name: &str, s: &Structure, e: &Value) {
    assert_eq!(s.has_unit_cell, e["has_cell"].as_bool().unwrap(), "{name}: has_cell");
    for r in 0..3 {
        for c in 0..3 {
            assert!((s.cell_vectors[r][c] - f(&e["cell"][r][c])).abs() < 1e-9, "{name}: cell {r} {c}: {} vs {}", s.cell_vectors[r][c], e["cell"][r][c]);
        }
        assert_eq!(s.pbc[r], e["pbc"][r].as_i64().unwrap() != 0, "{name}: pbc");
        assert!((s.cell_offset[r] - f(&e["offset"][r])).abs() < 1e-12, "{name}: offset");
    }
    assert_eq!(s.load_note, e["load_note"].as_str().unwrap(), "{name}: load note");
    assert_eq!(s.ipf_load_status, e["ipf_load_status"].as_str().unwrap(), "{name}: IPF status");
    let atoms = e["atoms"].as_array().unwrap();
    assert_eq!(s.atoms.len(), atoms.len(), "{name}: atom count");
    for (i, (a, x)) in s.atoms.iter().zip(atoms).enumerate() {
        assert_eq!(a.symbol, x["symbol"].as_str().unwrap(), "{name}: atom {i} symbol");
        assert_eq!(a.atomic_number as i64, x["z"].as_i64().unwrap(), "{name}: atom {i} Z");
        let p = &x["position"];
        for (k, v) in [a.x, a.y, a.z].iter().enumerate() {
            assert!((v - f(&p[k])).abs() < 1e-9, "{name}: atom {i} position {k}: {v} vs {}", p[k]);
        }
        assert_eq!(a.fixed as i64, x["fixed"].as_i64().unwrap(), "{name}: atom {i} fixed");
        assert!((a.magmom - f(&x["magmom"])).abs() < 1e-12, "{name}: atom {i} magmom");
        assert!((a.charge - f(&x["charge"])).abs() < 1e-12, "{name}: atom {i} charge");
        for k in 0..3 {
            assert!((a.velocity[k] - f(&x["velocity"][k])).abs() < 1e-15, "{name}: atom {i} velocity");
        }
        let colour = [a.r, a.g, a.b];
        for k in 0..3 {
            assert_eq!(colour[k], f(&x["colour"][k]) as f32, "{name}: atom {i} colour");
        }
    }
    let colours = e["grain_colors"].as_array().unwrap();
    assert_eq!(s.grain_colors.len(), colours.len(), "{name}: grain colour count");
    let mut differing = 0;
    for (i, (c, x)) in s.grain_colors.iter().zip(colours).enumerate() {
        let same = (0..3).all(|k| (c[k] - f(&x[k]) as f32).abs() <= 2e-6);
        if !same {
            differing += 1;
            if differing <= 3 {
                eprintln!("{name}: grain colour {i}: {c:?} vs {x}");
            }
        }
    }
    assert_eq!(differing, 0, "{name}: {differing} grain colours differ");
    let regions: Vec<i64> = e["grain_region_ids"].as_array().unwrap().iter().map(|v| v.as_i64().unwrap()).collect();
    assert_eq!(s.grain_region_ids.iter().map(|&r| r as i64).collect::<Vec<_>>(), regions, "{name}: regions");
}

#[test]
fn loads_equal_the_cpp_loader() {
    let g = golden();
    let mut checked = 0;
    for case in g["loads"].as_array().unwrap() {
        let name = case["file"].as_str().unwrap();
        let path = fixture_path(name);
        let result = af_io::load_structure_from_file(&path);
        if IMPROVED.contains(&name) {
            assert!(result.is_ok(), "{name}: {result:?}");
            continue;
        }
        if case["ok"].as_bool().unwrap() {
            let s = result.unwrap_or_else(|e| panic!("{name}: {e}"));
            compare(name, &s, &case["structure"]);
        } else {
            let expected = case["error"].as_str().unwrap();
            let error = result.err().unwrap_or_else(|| panic!("{name}: loaded, but the C++ failed: {expected}"));
            // Messages naming the file path differ only by the path.
            if expected.starts_with("File not found") {
                assert!(error.starts_with("File not found"), "{name}: {error}");
            } else {
                assert_eq!(error, expected, "{name}");
            }
        }
        checked += 1;
    }
    assert!(checked > 80, "only {checked} cases");
}

#[test]
fn native_writers_equal_the_cpp_text() {
    let g = golden();
    let native = ["cif", "xyz", "extxyz", "vasp", "data", "pwi", "cell", "gen", "xsf", "cp2k"];
    for case in g["saves"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        let format = case["format"].as_str().unwrap();
        let s = af_io::load_structure_from_file(&fixture_path(&format!("fixtures/{source}"))).unwrap();
        let result = af_io::structure_text(&s, format);
        if !case["ok"].as_bool().unwrap() {
            assert_eq!(result.unwrap_err(), case["error"].as_str().unwrap(), "{source} as {format}");
            continue;
        }
        let text = result.unwrap_or_else(|e| panic!("{source} as {format}: {e}"));
        if native.contains(&format) {
            assert_eq!(text, case["text"].as_str().unwrap(), "{source} as {format}");
        }
    }
}

#[test]
fn every_format_round_trips() {
    // Structures saved in each format read back with the same atoms (to the
    // precision the format keeps) and cell.
    let dir = std::env::temp_dir().join("af_io_round_trip");
    std::fs::create_dir_all(&dir).unwrap();
    for source in ["cell.extxyz", "nacl.vasp", "selective.vasp", "abs.cell", "charge_triclinic.data", "labels.xyz", "aspirin.mol"] {
        let s = af_io::load_structure_from_file(&fixture_path(&format!("fixtures/{source}"))).unwrap();
        for (format, ext, precision, keeps_cell) in [
            ("cif", "cif", 1e-9, true),
            ("xyz", "xyz", 1e-9, true),
            ("vasp", "vasp", 1e-9, true),
            ("data", "data", 1e-8, true),
            ("pwi", "pwi", 1e-9, true),
            ("cell", "cell", 1e-9, true),
            ("gen", "gen", 1e-9, true),
            ("xsf", "xsf", 1e-9, true),
            ("pdb", "pdb", 1e-3, true),
            ("sdf", "sdf", 1e-4, false),
            ("mol", "mol", 1e-4, false),
            ("mol2", "mol2", 1e-4, false),
            ("gjf", "gjf", 1e-5, true),
        ] {
            let path = dir.join(format!("out.{ext}")).to_string_lossy().into_owned();
            if af_io::save_structure(&s, &path, format).is_err() {
                continue;  // formats that need a cell
            }
            let back = af_io::load_structure_from_file(&path).unwrap_or_else(|e| panic!("{source} as {format}: {e}"));
            assert_eq!(back.atoms.len(), s.atoms.len(), "{source} as {format}");
            // A structure without a cell is written into a LAMMPS box with room to spare.
            if keeps_cell && format != "pdb" && !(format == "data" && !s.has_unit_cell) {
                assert_eq!(back.has_unit_cell, s.has_unit_cell, "{source} as {format}: cell");
            }
            // Same multiset of species and positions (writers may group by element).
            let key = |st: &Structure| {
                let mut v: Vec<(i32, [i64; 3])> =
                    st.atoms.iter().map(|a| (a.atomic_number, [a.x, a.y, a.z].map(|x| (x / (precision * 10.0)).round() as i64))).collect();
                v.sort();
                v
            };
            // LAMMPS data moves a cell-less structure into its new box.
            let shifted = format == "data" && !s.has_unit_cell;
            if !shifted && (!s.has_unit_cell || keeps_cell && format != "pdb") {
                let (a, b) = (key(&s), key(&back));
                let mismatches = a.iter().zip(&b).filter(|(x, y)| x.0 != y.0 || (0..3).any(|k| (x.1[k] - y.1[k]).abs() > 1)).count();
                assert_eq!(mismatches, 0, "{source} as {format}");
            }
        }
    }
}

#[test]
fn orientation_sidecar_round_trips_as_in_cpp() {
    let g = golden();
    let mut s = af_io::load_structure_from_file(&fixture_path("fixtures/nacl.vasp")).unwrap();
    s.grain_colors.clear();
    s.grain_region_ids.clear();
    for i in 0..s.atoms.len() {
        s.grain_colors.push([0.1 * i as f32, 0.5, 1.0 - 0.1 * i as f32]);
        s.grain_region_ids.push((i % 3) as i32);
    }
    s.ipf_load_status.clear();
    let dir = std::env::temp_dir().join("af_io_sidecar");
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("with_ipf.xyz").to_string_lossy().into_owned();
    af_io::save_structure(&s, &file, "xyz").unwrap();
    let sidecar = std::fs::read_to_string(dir.join("with_ipf.atomforge-ipf")).unwrap();
    // Byte for byte: both write in text mode (CRLF on Windows, where the
    // reference was made; LF elsewhere).
    let expected = g["sidecar"]["sidecar_text"].as_str().unwrap();
    let expected = if cfg!(windows) { expected.to_string() } else { expected.replace("\r\n", "\n") };
    assert_eq!(sidecar, expected);
    let back = af_io::load_structure_from_file(&file).unwrap();
    compare("sidecar reload", &back, &g["sidecar"]["reloaded"]);
}

/// Open Babel (the reader the C++ build used) reads the files the port
/// writes in its place: same atoms. Skipped when obabel is not installed.
#[test]
fn open_babel_reads_the_port_molecular_writers() {
    let candidates = ["C:/msys64/ucrt64/bin/obabel.exe", "/usr/bin/obabel", "/usr/local/bin/obabel"];
    let Some(obabel) = candidates.iter().find(|p| Path::new(p).exists()) else {
        eprintln!("obabel not found; skipped");
        return;
    };
    let plugins = |base: &str| std::fs::read_dir(base).ok().and_then(|mut d| d.next()).and_then(|e| e.ok()).map(|e| e.path());
    let dir = std::env::temp_dir().join("af_io_obabel");
    std::fs::create_dir_all(&dir).unwrap();
    for source in ["aspirin.mol", "cell.extxyz", "labels.xyz"] {
        let s = af_io::load_structure_from_file(&fixture_path(&format!("fixtures/{source}"))).unwrap();
        for format in ["pdb", "sdf", "mol2", "mol"] {
            let path = dir.join(format!("port.{format}"));
            af_io::save_structure(&s, &path.to_string_lossy(), format).unwrap();
            let mut command = std::process::Command::new(obabel);
            command.arg(&path).arg("-oxyz");
            if let Some(lib) = plugins("C:/msys64/ucrt64/lib/openbabel") {
                command.env("BABEL_LIBDIR", lib);
            }
            if let Some(data) = plugins("C:/msys64/ucrt64/share/openbabel") {
                command.env("BABEL_DATADIR", data);
            }
            let output = command.output().unwrap();
            let text = String::from_utf8_lossy(&output.stdout);
            let lines: Vec<&str> = text.lines().collect();
            assert!(!lines.is_empty(), "{source} as {format}: obabel read nothing: {}", String::from_utf8_lossy(&output.stderr));
            assert_eq!(lines[0].trim().parse::<usize>().unwrap(), s.atoms.len(), "{source} as {format}");
            for (a, line) in s.atoms.iter().zip(&lines[2..]) {
                let w: Vec<&str> = line.split_whitespace().collect();
                assert_eq!(w[0], a.symbol, "{source} as {format}");
                let tolerance = if format == "pdb" { 1e-3 } else { 1e-4 };
                for (k, v) in [a.x, a.y, a.z].iter().enumerate() {
                    assert!((w[1 + k].parse::<f64>().unwrap() - v).abs() <= tolerance, "{source} as {format}: {line}");
                }
            }
        }
    }
}

#[test]
fn trajectories_equal_the_cpp_reader() {
    let g = golden();
    for case in g["trajectories"].as_array().unwrap() {
        let file = case["file"].as_str().unwrap();
        let types = case["types"].as_str().unwrap();
        af_io::trajectory::set_lammps_type_elements(types).unwrap();
        let path = root().join("tests/fixtures/traj").join(file).to_string_lossy().into_owned();
        let result = af_io::trajectory::load_trajectory(&path);
        af_io::trajectory::set_lammps_type_elements("").unwrap();
        let expected_error = case["error"].as_str().unwrap();
        if !expected_error.is_empty() {
            assert_eq!(result.unwrap_err(), expected_error, "{file} ({types})");
            continue;
        }
        let frames = result.unwrap_or_else(|e| panic!("{file}: {e}"));
        let expected = case["frames"].as_array().unwrap();
        assert_eq!(frames.len(), expected.len(), "{file}: frames");
        for (k, (frame, e)) in frames.iter().zip(expected).enumerate() {
            compare(&format!("{file} frame {k}"), frame, e);
        }
    }
}
