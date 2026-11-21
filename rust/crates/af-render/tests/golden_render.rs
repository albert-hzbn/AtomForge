//! af-render's scene data against the C++ renderer: the cases of
//! tools/golden/golden_render.cpp (tests/golden/render.json), compared
//! exactly (instances, ambient occlusion, bonds, default view, frame matrices).

use std::path::PathBuf;

use af_core::display::{default_look, DisplaySettings};
use af_core::json::Json;
use af_core::occlusion::darken_occluded;
use af_core::structure::{AtomSite, Structure};
use af_render::scene::*;
use glam::{Mat4, Vec3};

fn f3(v: Vec3) -> Json {
    Json::Array(vec![v.x.into(), v.y.into(), v.z.into()])
}

fn c3(v: [f32; 3]) -> Json {
    Json::Array(vec![v[0].into(), v[1].into(), v[2].into()])
}

fn vecs(list: &[Vec3]) -> Json {
    Json::Array(list.iter().map(|&v| f3(v)).collect())
}

fn colours(list: &[[f32; 3]]) -> Json {
    Json::Array(list.iter().map(|&v| c3(v)).collect())
}

fn floats(list: &[f32]) -> Json {
    Json::Array(list.iter().map(|&v| v.into()).collect())
}

fn ints(list: &[i32]) -> Json {
    Json::Array(list.iter().map(|&v| v.into()).collect())
}

fn mat(m: &Mat4) -> Json {
    Json::Array(m.to_cols_array().iter().map(|&v| v.into()).collect())
}

fn instances_json(d: &InstanceData) -> Json {
    let mut o = Json::object();
    o.set("positions", vecs(&d.positions));
    o.set("colors", colours(&d.colors));
    o.set("scales", floats(&d.scales));
    o.set("bond_radii", floats(&d.bond_radii));
    o.set("shininess", floats(&d.shininess));
    o.set("atomic_numbers", ints(&d.atomic_numbers));
    o.set("box_lines", vecs(&d.box_lines));
    o.set("atom_indices", ints(&d.atom_indices));
    o.set("orbit_center", f3(d.orbit_center));
    o
}

fn bonds_json(b: &Bonds) -> Json {
    let mut o = Json::object();
    o.set("starts", vecs(&b.starts));
    o.set("ends", vecs(&b.ends));
    o.set("colors_a", colours(&b.colors_a));
    o.set("colors_b", colours(&b.colors_b));
    o.set("radii", floats(&b.radii));
    o
}

fn frame_json(f: &FrameView) -> Json {
    let mut o = Json::object();
    o.set("projection", mat(&f.projection));
    o.set("view", mat(&f.view));
    o.set("light_mvp", mat(&f.light_mvp));
    o.set("light_position", f3(f.light_position));
    o.set("camera_position", f3(f.camera_position));
    o.set("depth", Json::Array(vec![f.depth_near.into(), f.depth_far.into()]));
    o
}

fn read_xyz(path: &str) -> Structure {
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

#[test]
fn scene_data_matches_cpp() {
    std::env::set_current_dir(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")).unwrap();
    let look = default_look();
    let shininess = vec![32.0f32; 119];
    let mut cases = Vec::new();
    for name in ["cu_fcc", "si_diamond", "nacl", "mg_hcp", "molecule", "si_super"] {
        let mut s = read_xyz(&format!("rust/tests/fixtures/render/{name}.xyz"));
        for a in s.atoms.iter_mut() {
            [a.r, a.g, a.b] = look.colours[a.atomic_number as usize];
        }
        let mut entry = Json::object();
        entry.set("name", name);
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let mut data = build_instance_data(&s, false, &identity, &look.radii, &shininess);
        entry.set("instances", instances_json(&data));
        entry.set("doubled", instances_json(&build_instance_data(&s, true, &[[2, 0, 0], [0, 1, 0], [0, 0, 1]], &look.radii, &shininess)));
        entry.set("skew", instances_json(&build_instance_data(&s, true, &[[1, 1, 0], [-1, 1, 0], [0, 0, 1]], &look.radii, &shininess)));
        for scale in data.scales.iter_mut() {
            *scale *= 1.3;
        }
        darken_occluded(&data.positions, &data.scales, &mut data.colors, 0.6);
        entry.set("occluded", colours(&data.colors));
        let mut bonds = Vec::new();
        for style in 0..3 {
            let bond_style = BondStyle {
                metal_bonds: style == 1,
                uniform_colour: style == 2,
                thickness: if style == 2 { 1.7 } else { 1.0 },
                tolerance: if style == 2 { 1.3 } else { 1.18 },
                ..Default::default()
            };
            bonds.push(bonds_json(&find_bonds(&data, &bond_style, None)));
        }
        entry.set("bonds", Json::Array(bonds));
        let scene = Scene::new(data, &BondStyle::default());
        let mut views = Vec::new();
        for (w, h) in [(1600, 1200), (64, 48), (300, 900)] {
            let mut camera = OrbitCamera::default();
            apply_default_view(&mut camera, &scene, w, h, true);
            let mut view = Json::object();
            view.set("distance", camera.distance);
            for variant in 0..3 {
                let mut c = camera;
                let mut display = DisplaySettings::default();
                if variant == 1 {
                    (c.yaw, c.pitch, c.roll) = (20.0, -15.0, 30.0);
                    display.light_follows_camera = !display.light_follows_camera;
                }
                if variant == 2 {
                    c.distance = 7.5;
                    display.light_azimuth = 200.0;
                    display.light_elevation = 80.0;
                }
                let frame = build_frame_view(&c, &scene, variant == 2, w, h, &display);
                view.set(&format!("frame{variant}"), frame_json(&frame));
            }
            views.push(view);
        }
        entry.set("views", Json::Array(views));
        cases.push(entry);
    }
    let golden = Json::parse(&std::fs::read_to_string("rust/tests/golden/render.json").unwrap()).unwrap();
    let expected = golden.items().unwrap();
    assert_eq!(cases.len(), expected.len());
    let mut differences = Vec::new();
    for (got, want) in cases.iter().zip(expected) {
        for (key, value) in want.members().unwrap() {
            let mine = got.find(key);
            if mine != Some(value) {
                let (g, w) = (mine.map(|m| m.dump(-1)).unwrap_or_default(), value.dump(-1));
                let first = g.bytes().zip(w.bytes()).position(|(a, b)| a != b).unwrap_or(g.len().min(w.len()));
                let from = first.saturating_sub(60);
                differences.push(format!(
                    "{} {key}:\n  {}\n  vs {}",
                    want.at("name").unwrap().string().unwrap(),
                    &g[from..(first + 80).min(g.len())],
                    &w[from..(first + 80).min(w.len())]
                ));
            }
        }
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}
