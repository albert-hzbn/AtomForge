//! Measurements and Atom Info against the C++ desktop: the selections of
//! tools/golden/golden_measure.cpp (tests/golden/measure.json) over the
//! render fixtures, compared as text.

use std::path::PathBuf;

use af_core::display::default_look;
use af_core::json::Json;
use af_core::structure::{AtomSite, Structure};
use af_render::measure::{angle, atom_info, dihedral, distance};
use af_render::scene::{build_instance_data, BondStyle, Scene};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn read_xyz(path: &PathBuf) -> Structure {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let n: usize = lines.next().unwrap().trim().parse().unwrap();
    let comment = lines.next().unwrap();
    let mut s = Structure::default();
    if let Some(values) = comment.strip_prefix("cell ") {
        let words: Vec<&str> = values.split_whitespace().collect();
        let v: Vec<f64> = words[..9].iter().map(|x| x.parse().unwrap()).collect();
        s.cell_vectors = [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]];
        if words.len() > 9 && words[9] == "origin" {
            s.cell_offset = [words[10].parse().unwrap(), words[11].parse().unwrap(), words[12].parse().unwrap()];
        }
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
            atomic_number: af_core::elements::atomic_number(w[0]),
            ..Default::default()
        });
    }
    s
}

fn strings(j: &Json) -> (String, String) {
    let items = j.items().unwrap();
    (items[0].string().unwrap().to_string(), items[1].string().unwrap().to_string())
}

fn indices(j: &Json) -> Vec<usize> {
    j.items().unwrap().iter().map(|v| v.number().unwrap() as usize).collect()
}

#[test]
fn measurements_match_cpp() {
    let golden = Json::parse(&std::fs::read_to_string(root().join("rust/tests/golden/measure.json")).unwrap()).unwrap();
    let look = default_look();
    let shininess = vec![32.0f32; 119];
    let mut compared = 0;
    for case in golden.items().unwrap() {
        let name = case.at("name").unwrap().string().unwrap();
        let mut s = read_xyz(&root().join(format!("rust/tests/fixtures/render/{name}.xyz")));
        if name == "molecule" {
            s.atoms[0].fixed = 5;
            s.atoms[0].magmom = 1.23456;
            s.atoms[0].charge = -0.5;
            s.atoms[0].velocity = [0.001, -2.5, 1e-7];
        }
        let data = build_instance_data(&s, false, &[[1, 0, 0], [0, 1, 0], [0, 0, 1]], &look.radii, &shininess);
        let scene = Scene::new(data, &BondStyle::default());
        for r in case.at("results").unwrap().items().unwrap() {
            let selection = indices(r.at("selection").unwrap());
            let d = distance(&selection, &scene, &s);
            assert_eq!((d.value.clone(), d.message.clone()), strings(r.at("distance").unwrap()), "{name} {selection:?} distance");
            assert_eq!(d.lines, indices(r.at("distance_lines").unwrap()), "{name} {selection:?} distance lines");
            let a = angle(&selection, &scene, &s);
            assert_eq!((a.value.clone(), a.message.clone()), strings(r.at("angle").unwrap()), "{name} {selection:?} angle");
            assert_eq!(a.lines, indices(r.at("angle_lines").unwrap()), "{name} {selection:?} angle lines");
            let t = dihedral(&selection, &scene, &s);
            assert_eq!((t.value, t.message), strings(r.at("dihedral").unwrap()), "{name} {selection:?} dihedral");
            assert_eq!(atom_info(&selection, &scene, &s, &look.radii), r.at("atom_info").unwrap().string().unwrap(), "{name} {selection:?} atom info");
            compared += 1;
        }
    }
    assert_eq!(compared, 6 * 17);
}

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn crystal_rotation_matches_cpp() {
    use af_render::scene::{rotate_crystal_view, OrbitCamera};
    let golden = Json::parse(&std::fs::read_to_string(root().join("rust/tests/golden/camera_rotate.json")).unwrap()).unwrap();
    let mut compared = 0;
    for case in golden.items().unwrap() {
        let start: Vec<f32> = case.at("start").unwrap().items().unwrap().iter().map(|v| v.number().unwrap() as f32).collect();
        let axis = case.at("axis").unwrap().number().unwrap() as usize;
        let angle = case.at("angle").unwrap().number().unwrap();
        let mut camera = OrbitCamera { yaw: start[0], pitch: start[1], roll: start[2], ..Default::default() };
        for want in case.at("cameras").unwrap().items().unwrap() {
            rotate_crystal_view(&mut camera, axis, angle);
            let want: Vec<f32> = want.items().unwrap().iter().map(|v| v.number().unwrap() as f32).collect();
            assert_eq!([camera.yaw, camera.pitch, camera.roll], [want[0], want[1], want[2]], "start {start:?} axis {axis} angle {angle}");
            compared += 1;
        }
    }
    assert_eq!(compared, 6 * 3 * 6 * 4);
}
