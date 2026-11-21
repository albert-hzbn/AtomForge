//! The grain-boundary colouring against the C++ desktop: the cases of
//! tools/golden/golden_colouring.cpp (tests/golden/colouring.json), compared
//! atom by atom.

use std::path::PathBuf;

use af_core::json::Json;
use af_core::structure::{AtomSite, Structure};
use af_render::colouring::{apply_colour_mode, grain_boundary_atoms, ColourMode, BOUNDARY_COLOUR, CRYSTAL_COLOUR};

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
        let v: Vec<f64> = values.split_whitespace().take(9).map(|x| x.parse().unwrap()).collect();
        s.cell_vectors = [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]];
        s.has_unit_cell = true;
        s.pbc = [true; 3];
    }
    for _ in 0..n {
        let w: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        if w.len() > 4 {
            s.grain_region_ids.push(w[4].parse().unwrap());
        }
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
fn grain_boundary_colouring_matches_cpp() {
    let golden = Json::parse(&std::fs::read_to_string(root().join("rust/tests/golden/colouring.json")).unwrap()).unwrap();
    let Json::Object(cases) = golden else { panic!("colouring.json is not an object") };
    assert_eq!(cases.len(), 8);
    for (name, flags) in cases.iter() {
        let s = read_xyz(&root().join(format!("rust/tests/fixtures/colouring/{name}.xyz")));
        let want: Vec<bool> = flags.items().unwrap().iter().map(|f| f.boolean().unwrap()).collect();
        let got = grain_boundary_atoms(&s);
        assert_eq!(got.len(), want.len(), "{name}");
        let differing: Vec<usize> = (0..want.len()).filter(|&i| got[i] != want[i]).collect();
        assert!(differing.is_empty(), "{name}: atoms {differing:?} differ ({} boundary atoms in C++)", want.iter().filter(|&&b| b).count());
        // The colours as painted.
        let mut painted = s.clone();
        apply_colour_mode(&mut painted, &af_core::elements::default_colours(), ColourMode::GrainBoundary, None);
        for (atom, &b) in painted.atoms.iter().zip(&want) {
            assert_eq!([atom.r, atom.g, atom.b], if b { BOUNDARY_COLOUR } else { CRYSTAL_COLOUR }, "{name}");
        }
    }
}

#[test]
fn colour_modes() {
    let mut s = read_xyz(&root().join("rust/tests/fixtures/colouring/three.xyz"));
    let colours = af_core::elements::default_colours();
    apply_colour_mode(&mut s, &colours, ColourMode::Element, None);
    assert_eq!([s.atoms[0].r, s.atoms[0].g, s.atoms[0].b], colours[29]);
    s.grain_colors = vec![[0.1, 0.2, 0.3]; 3];
    apply_colour_mode(&mut s, &colours, ColourMode::CrystalOrientation, None);
    assert_eq!([s.atoms[2].r, s.atoms[2].g, s.atoms[2].b], [0.1, 0.2, 0.3]);
    let property = vec![[0.5, 0.5, 0.0]; 3];
    apply_colour_mode(&mut s, &colours, ColourMode::AtomProperty, Some(&property));
    assert_eq!([s.atoms[1].r, s.atoms[1].g, s.atoms[1].b], [0.5, 0.5, 0.0]);
    // Grain colours of the wrong length leave the element colours.
    s.grain_colors.pop();
    apply_colour_mode(&mut s, &colours, ColourMode::CrystalOrientation, None);
    assert_eq!([s.atoms[0].r, s.atoms[0].g, s.atoms[0].b], colours[29]);
}
