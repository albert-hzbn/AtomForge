//! af-cloud against the C++ cloud library: the cases of
//! tools/golden/golden_cloud.cpp (tests/golden/cloud.json), with every built
//! .afcloud compared byte for byte with the C++ file (one build thread).

use std::path::{Path, PathBuf};

use af_cloud::*;
use af_core::json::Json;

fn vec(v: &[f64; 3]) -> Json {
    Json::Array(vec![v[0].into(), v[1].into(), v[2].into()])
}

fn info_json(info: &CloudInfo) -> Json {
    let mut o = Json::object();
    o.set("atoms", info.atoms as f64);
    o.set("lower", vec(&info.lower));
    o.set("upper", vec(&info.upper));
    o.set("has_cell", info.has_cell);
    o.set("cell", Json::Array(info.cell.iter().map(vec).collect()));
    o.set("cell_origin", vec(&info.cell_origin));
    o.set(
        "species",
        Json::Array(
            info.species
                .iter()
                .map(|s| Json::Array(vec![s.name.as_str().into(), s.color[0].into(), s.color[1].into(), s.color[2].into(), s.radius.into()]))
                .collect(),
        ),
    );
    o.set(
        "chunks",
        Json::Array(
            info.chunks
                .iter()
                .map(|c| {
                    Json::Array(vec![
                        (c.offset as f64).into(),
                        (c.count as f64).into(),
                        c.lower[0].into(),
                        c.lower[1].into(),
                        c.lower[2].into(),
                        c.upper[0].into(),
                        c.upper[1].into(),
                        c.upper[2].into(),
                    ])
                })
                .collect(),
        ),
    );
    o.set("source", info.source.as_str());
    o
}

fn options(chunk: u64, seed: u64) -> BuildOptions {
    BuildOptions { chunk_atoms: chunk, threads: 1, seed, temporary: None }
}

#[test]
fn clouds_match_cpp() {
    std::env::set_current_dir(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")).unwrap();
    let input = "rust/tests/fixtures/cloud/";
    let out = "rust/target/golden/cloud_rust/";
    std::fs::create_dir_all(out).unwrap();
    let fixture = |name: &str| PathBuf::from(format!("{input}{name}"));
    let mut cases: Vec<Json> = Vec::new();
    let mut attempt = |name: &str, result: Result<Json, String>| {
        let mut entry = Json::object();
        entry.set("name", name);
        match result {
            Ok(r) => entry.set("result", r),
            Err(e) => entry.set("error", e),
        }
        cases.push(entry);
    };
    let path = |name: &str| PathBuf::from(format!("{out}{name}.afcloud"));
    let names = [String::from("Cu"), String::from("Ni")];
    let builds: Vec<(&str, Box<dyn Fn(&Path) -> Result<CloudInfo, String>>)> = vec![
        ("fcc", Box::new(|p| generate_crystal(p, "fcc", 3.615, 0.0, [10, 10, 10], "Cu", "", None, options(1024, 12345)))),
        ("hcp", Box::new(|p| generate_crystal(p, "hcp", 3.21, 0.0, [6, 4, 5], "Mg", "Zn", None, options(2000, 7)))),
        ("diamond", Box::new(|p| generate_crystal(p, "diamond", 5.43, 0.0, [5, 5, 5], "si", "", None, options(1024, 12345)))),
        ("bcc", Box::new(|p| generate_crystal(p, "bcc", 2.87, 0.0, [3, 3, 3], "Fe", "type 9", None, options(1024, 12345)))),
        ("sc", Box::new(|p| generate_crystal(p, "sc", 2.5, 0.0, [1, 1, 1], "Po", "", None, options(1024, 12345)))),
        ("hcp_c", Box::new(|p| generate_crystal(p, "hcp", 2.95, 4.68, [3, 2, 2], "Ti", "", None, options(1024, 12345)))),
        ("bad_a", Box::new(|p| generate_crystal(p, "fcc", 0.0, 0.0, [1, 1, 1], "Cu", "", None, options(1024, 12345)))),
        ("bad_lattice", Box::new(|p| generate_crystal(p, "a15", 3.0, 0.0, [1, 1, 1], "Cu", "", None, options(1024, 12345)))),
        ("bad_cells", Box::new(|p| generate_crystal(p, "fcc", 3.0, 0.0, [1, 0, 1], "Cu", "", None, options(1024, 12345)))),
        ("xyz", Box::new(|p| build_from_xyz(&fixture("big.extxyz"), p, 0, None, options(1024, 12345)))),
        ("xyz_frame1", Box::new(|p| build_from_xyz(&fixture("two_frames.xyz"), p, 1, None, options(1500, 99)))),
        ("xyz_bad", Box::new(|p| build_from_xyz(&fixture("bad_coordinate.xyz"), p, 0, None, options(1024, 12345)))),
        ("xyz_short", Box::new(|p| build_from_xyz(&fixture("short.xyz"), p, 0, None, options(1024, 12345)))),
        ("dump_elements", Box::new(|p| build_from_lammps_dump(&fixture("elements.dump"), p, &[], 0, None, options(1024, 12345)))),
        ("dump_triclinic", Box::new(|p| build_from_lammps_dump(&fixture("triclinic.dump"), p, &names, 0, None, options(1024, 12345)))),
        ("dump_scaled_frame1", Box::new(|p| build_from_lammps_dump(&fixture("scaled_two_frames.dump"), p, &[], 1, None, options(1024, 12345)))),
        ("dump_no_positions", Box::new(|p| build_from_lammps_dump(&fixture("no_positions.dump"), p, &[], 0, None, options(1024, 12345)))),
        ("dump_short_line", Box::new(|p| build_from_lammps_dump(&fixture("short_line.dump"), p, &[], 0, None, options(1024, 12345)))),
        ("dump_not", Box::new(|p| build_from_lammps_dump(&fixture("not_a_dump.dump"), p, &[], 0, None, options(1024, 12345)))),
    ];
    let mut built = Vec::new();
    for (name, build) in &builds {
        let result = build(&path(name)).map(|info| info_json(&info));
        if result.is_ok() {
            built.push(*name);
        }
        attempt(&format!("build {name}"), result);
    }
    for name in ["fcc", "hcp", "xyz", "dump_triclinic"] {
        let result = (|| -> Result<Json, String> {
            let file = CloudFile::open(&path(name))?;
            let mut o = Json::object();
            o.set("info", info_json(file.info()));
            let mut atoms = Vec::new();
            for (c, chunk) in file.info().chunks.iter().enumerate() {
                let count = chunk.count.min(6);
                for a in file.read(c, chunk.count - count, count)? {
                    let p = file.position(c, &a)?;
                    atoms.push(Json::Array(vec![i32::from(a.x).into(), i32::from(a.y).into(), i32::from(a.z).into(), i32::from(a.species).into(), p[0].into(), p[1].into(), p[2].into()]));
                }
            }
            o.set("atoms", Json::Array(atoms));
            Ok(o)
        })();
        attempt(&format!("read {name}"), result);
    }
    attempt("read past end", CloudFile::open(&path("sc")).and_then(|f| f.read(0, 0, 5)).map(|_| Json::Bool(true)));
    let incomplete = format!("{out}incomplete.afcloud");
    std::fs::write(&incomplete, b"AFCLOUD1\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0").unwrap();
    // The C++ names its own output folder in this message.
    let cpp_incomplete = "rust/tests/golden/cloud/incomplete.afcloud".to_string();
    for p in [cpp_incomplete.clone(), format!("{input}text.txt"), format!("{input}missing.afcloud")] {
        let name = Path::new(&p).file_name().unwrap().to_string_lossy().into_owned();
        attempt(&format!("open {name}"), CloudFile::open(Path::new(&p)).map(|f| Json::from(f.info().atoms as f64)));
    }
    for p in [format!("{out}fcc.afcloud"), format!("{input}text.txt"), format!("{input}missing.afcloud"), incomplete.clone()] {
        let name = Path::new(&p).file_name().unwrap().to_string_lossy().into_owned();
        attempt(&format!("is cloud {name}"), Ok(Json::Bool(is_cloud_file(Path::new(&p)))));
    }
    for name in ["big.extxyz", "two_frames.xyz", "elements.dump", "triclinic.dump", "text.txt", "missing.xyz"] {
        attempt(&format!("peek {name}"), Ok(Json::from(peek_atom_count(&fixture(name)) as f64)));
    }

    let golden = Json::parse(&std::fs::read_to_string("rust/tests/golden/cloud.json").unwrap()).unwrap();
    let expected = golden.items().unwrap();
    assert_eq!(cases.len(), expected.len(), "case count");
    let mut differences = Vec::new();
    for (got, want) in cases.iter().zip(expected) {
        if got != want {
            differences.push(format!("{}:\n  {}\n  vs {}", want.at("name").unwrap().string().unwrap(), &got.dump(-1)[..got.dump(-1).len().min(300)], &want.dump(-1)[..want.dump(-1).len().min(300)]));
        }
    }
    for name in built {
        let cpp = std::fs::read(format!("rust/tests/golden/cloud/{name}.afcloud")).ok();
        let rust = std::fs::read(path(name)).ok();
        if cpp != rust {
            differences.push(format!("{name}.afcloud differs from the C++ file"));
        }
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}
