//! Projects against the C++ model/Workspace: the files of
//! tools/golden/golden_project.cpp (formats 1-7) are read, checked field by
//! field and written again byte for byte.

use std::path::PathBuf;

use af_electronic::project::{parse_project, project_bytes, save_project, load_project, Workspace};

fn golden(format: i32) -> Vec<u8> {
    std::fs::read(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/project_v{format}.afproject"))).unwrap()
}

#[test]
fn projects_round_trip_byte_for_byte() {
    for format in 1..=7 {
        let bytes = golden(format);
        let tabs = parse_project(&bytes).unwrap_or_else(|e| panic!("format {format}: {e}"));
        assert_eq!(tabs.len(), 2);
        assert_eq!(project_bytes(&tabs, format).unwrap(), bytes, "format {format}");
    }
}

#[test]
fn project_fields() {
    let tabs = parse_project(&golden(7)).unwrap();
    let w = &tabs[0];
    let s = &w.structure;
    assert_eq!(s.atoms.len(), 5);
    assert_eq!(s.atoms[1].symbol, "Al");
    assert_eq!(s.atoms[3].y, 1.0);
    assert_eq!(s.atoms[3].fixed, 3);
    assert_eq!(s.atoms[2].charge, -0.3 * 2.0);
    assert_eq!(s.atoms[4].velocity, [4e-3, -2.5, 1.0 / 7.0]);
    assert_eq!(s.pbc, [true, false, true]);
    assert_eq!(s.cell_offset, [0.5, -0.25, 0.0]);
    assert_eq!(s.grain_region_ids, vec![0, 0, 1, 1, 2]);
    assert!(s.atom_property[1].is_nan());
    assert_eq!(s.atom_property_name, "coordination");
    assert_eq!(s.frozen_property[4], -4.0);
    assert_eq!(s.dislocation_loop_points[1], [-1.5, 0.25, 9.0]);
    assert!(s.dislocation_detection_done);
    assert_eq!(w.volume.fields[0].values.len(), 12);
    assert_eq!(w.volume.fields[0].cell.z_axis.x, 0.5);
    assert_eq!(w.volume.sites[1].number, 13);
    assert_eq!(w.surface.vertices.len(), 3);
    assert_eq!(w.columns, 3);
    assert_eq!(w.settings["isovalue"], 0.05);
    assert_eq!(w.camera, vec![30.0, -20.0, 5.0, 12.5, 0.1, 0.2, 0.3]);
    assert_eq!(w.history, vec!["Opened", "Replicate 2 1 1"]);
    assert!(w.has_pipeline_input);
    assert_eq!(w.pipeline_input.atoms.len(), 3);
    assert_eq!(w.pipeline_input.atoms[2].fixed, 2);
    assert_eq!(w.display, "ambient 0.3\n");
    assert_eq!(tabs[1].title, "second");
    // Older formats leave the later fields at their defaults.
    let old = parse_project(&golden(1)).unwrap();
    assert!(old[0].structure.atom_property.is_empty());
    assert_eq!(old[0].structure.atoms[3].fixed, 0);
    assert_eq!(old[0].structure.atoms[2].charge, 0.0);
}

#[test]
fn project_errors() {
    let mut bytes = golden(7);
    bytes.push(0);
    assert_eq!(parse_project(&bytes).unwrap_err(), "Unexpected trailing project data");
    let bytes = golden(7);
    assert_eq!(parse_project(&bytes[..bytes.len() - 3]).unwrap_err(), "Truncated project");
    let mut newer = golden(7);
    newer[8 + 18] = b'8';
    assert_eq!(parse_project(&newer).unwrap_err(), "Unsupported project version (saved by a newer AtomForge?)");
    assert_eq!(project_bytes(&[], 7).unwrap_err(), "A project requires 1-1000 tabs");
    let bad = Workspace { camera: vec![1.0; 3], ..Default::default() };
    assert_eq!(project_bytes(&[bad], 7).unwrap_err(), "Invalid project camera");
    // Saving replaces the file atomically and reads back.
    let dir = std::env::temp_dir().join(format!("af_project_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("p.afproject").to_string_lossy().into_owned();
    let tabs = parse_project(&golden(7)).unwrap();
    save_project(&tabs, &path, 7).unwrap();
    save_project(&tabs, &path, 7).unwrap();
    assert_eq!(std::fs::read(&path).unwrap(), golden(7));
    assert_eq!(load_project(&path).unwrap().len(), 2);
    let _ = std::fs::remove_dir_all(&dir);
}
