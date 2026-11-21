//! af-core against the C++ reference outputs (rust/tests/golden, written by
//! tools/golden/golden_dump.cpp from the C++ libraries).

use af_core::display::{self, DisplaySettings};
use af_core::{elements, format, neighbours, AtomSite, Structure};
use serde_json::Value;

fn golden(name: &str) -> Value {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden").join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    serde_json::from_str(&text).expect("golden JSON")
}

fn f32s(v: &Value) -> Vec<f32> {
    v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap() as f32).collect()
}

pub fn structure_from(v: &Value) -> Structure {
    let mut s = Structure::default();
    s.has_unit_cell = v["has_cell"].as_bool().unwrap();
    for r in 0..3 {
        for c in 0..3 {
            s.cell_vectors[r][c] = v["cell"][r][c].as_f64().unwrap();
        }
        s.pbc[r] = v["pbc"][r].as_i64().unwrap() != 0;
        s.cell_offset[r] = v["offset"][r].as_f64().unwrap();
    }
    for a in v["atoms"].as_array().unwrap() {
        let p = &a["position"];
        let mut site = AtomSite::new(a["symbol"].as_str().unwrap(), [p[0].as_f64().unwrap(), p[1].as_f64().unwrap(), p[2].as_f64().unwrap()]);
        site.atomic_number = a["z"].as_i64().unwrap() as i32;
        s.atoms.push(site);
    }
    s
}

#[test]
fn element_tables_equal_the_cpp_ones() {
    let g = golden("elements.json");
    let radii = elements::literature_covalent_radii();
    for e in g["elements"].as_array().unwrap() {
        let z = e["z"].as_i64().unwrap() as i32;
        assert_eq!(elements::symbol(z), e["symbol"].as_str().unwrap(), "symbol {z}");
        assert_eq!(elements::name(z), e["name"].as_str().unwrap(), "name {z}");
        assert_eq!(radii[z as usize], e["covalent_radius"].as_f64().unwrap() as f32, "radius {z}");
        assert_eq!(elements::atomic_mass(z), e["mass"].as_f64().unwrap(), "mass {z}");
        assert_eq!(elements::default_colour(z).to_vec(), f32s(&e["colour"]), "colour {z}");
        if z > 0 {
            assert_eq!(elements::atomic_number(elements::symbol(z)), z);
        }
    }
    for (k, scheme) in g["schemes"].as_array().unwrap().iter().enumerate() {
        let schemes = elements::colour_schemes();
        assert_eq!(schemes[k].name, scheme["name"].as_str().unwrap());
        let colours = elements::scheme_colours(k as i32);
        for (z, c) in scheme["colours"].as_array().unwrap().iter().enumerate() {
            assert_eq!(colours[z].to_vec(), f32s(c), "scheme {k} element {z}");
        }
        assert_eq!(elements::matching_scheme(&colours), k as i32);
    }
}

#[test]
fn neighbour_pairs_equal_the_cpp_ones_in_order() {
    let g = golden("neighbours.json");
    for case in g["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let s = structure_from(&case["structure"]);
        let pairs = neighbours::within(&s, case["cutoff"].as_f64().unwrap(), case["periodic"].as_bool().unwrap()).unwrap();
        let expected = case["pairs"].as_array().unwrap();
        assert_eq!(pairs.len(), expected.len(), "{name}: pair count");
        for (p, e) in pairs.iter().zip(expected) {
            assert_eq!(p.i as i64, e[0].as_i64().unwrap(), "{name}");
            assert_eq!(p.j as i64, e[1].as_i64().unwrap(), "{name}");
            for c in 0..3 {
                assert_eq!(p.d[c], e[2 + c].as_f64().unwrap(), "{name}: d");
            }
            assert_eq!(p.r, e[5].as_f64().unwrap(), "{name}: r");
        }
    }
}

#[test]
fn number_formatting_equals_cpp_streams() {
    let g = golden("formatting.json");
    for case in g["cases"].as_array().unwrap() {
        let value = case["value"].as_f64().unwrap();
        let precision = case["precision"].as_u64().unwrap() as usize;
        assert_eq!(format::format_g(value, precision), case["text"].as_str().unwrap(), "{value} at {precision}");
    }
}

#[test]
fn display_settings_text_equals_cpp() {
    let g = golden("display.json");
    assert_eq!(display::to_text(&DisplaySettings::default()), g["default_text"].as_str().unwrap());
    for (i, p) in g["presets"].as_array().unwrap().iter().enumerate() {
        let mut s = DisplaySettings::default();
        display::apply_preset(&mut s, &display::PRESETS[i]);
        assert_eq!(display::PRESETS[i].name, p["name"].as_str().unwrap());
        assert_eq!(display::to_text(&s), p["text"].as_str().unwrap(), "preset {i}");
        assert_eq!(display::matching_preset(&s) as i64, p["matching"].as_i64().unwrap());
    }
    for case in g["parsed"].as_array().unwrap() {
        let s = display::from_text(case["input"].as_str().unwrap(), &DisplaySettings::default());
        assert_eq!(display::to_text(&s), case["text"].as_str().unwrap());
        assert_eq!(display::matching_preset(&s) as i64, case["matching"].as_i64().unwrap());
        let light = display::light_direction(&s);
        for (c, e) in light.iter().zip(f32s(&case["light"])) {
            assert!((c - e).abs() < 1e-6, "light {light:?}");
        }
    }
    let mut look = display::default_look();
    look.colours[29] = [0.1, 0.2, 0.3];
    look.radii[8] = 0.7;
    look.display.atom_scale = 1.5;
    let text = display::look_text(&look);
    assert_eq!(text, g["look"]["text"].as_str().unwrap());
    let back = display::read_look(&(text + "radius.Xx=2\ncolour.Fe=0.5 0.5\nradius.Au=9\n"));
    assert_eq!(display::look_text(&back), g["look"]["reread_text"].as_str().unwrap());
}

#[test]
fn sort_order_equals_libstdcxx_std_sort() {
    let g = golden("sort.json");
    for (k, case) in g["cases"].as_array().unwrap().iter().enumerate() {
        let keys: Vec<f32> = case["keys"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap() as f32).collect();
        let mut v: Vec<(f32, usize)> = keys.iter().enumerate().map(|(i, &k)| (k, i)).collect();
        af_core::cxx::sort_by(&mut v, |a, b| a.0 < b.0);
        let order: Vec<usize> = v.iter().map(|p| p.1).collect();
        let expected: Vec<usize> = case["order"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as usize).collect();
        assert_eq!(order, expected, "case {k} ({} elements)", keys.len());
    }
}

#[test]
fn selection_order_equals_libstdcxx_nth_element_and_partial_sort() {
    let g = golden("select.json");
    let indices = |v: &Value| -> Vec<usize> { v.as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as usize).collect() };
    for (k, case) in g["cases"].as_array().unwrap().iter().enumerate() {
        let keys: Vec<f32> = case["keys"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap() as f32).collect();
        let nth = case["nth"].as_u64().unwrap() as usize;
        let input: Vec<(f32, usize)> = keys.iter().enumerate().map(|(i, &k)| (k, i)).collect();
        let mut v = input.clone();
        af_core::cxx::nth_element_by(&mut v, nth, |a, b| a.0 < b.0);
        assert_eq!(v.iter().map(|p| p.1).collect::<Vec<_>>(), indices(&case["nth_element"]), "nth_element case {k} (n {}, nth {nth})", keys.len());
        let mut v = input.clone();
        af_core::cxx::partial_sort_by(&mut v, nth, |a, b| a.0 < b.0);
        assert_eq!(v.iter().map(|p| p.1).collect::<Vec<_>>(), indices(&case["partial_sort"]), "partial_sort case {k} (n {}, nth {nth})", keys.len());
        let mut v = input.clone();
        af_core::cxx::stable_sort_by(&mut v, |a, b| {
            if (a.0 - b.0).abs() > 1.5 {
                return a.0 < b.0;
            }
            (a.1 * 31) % 7 < (b.1 * 31) % 7
        });
        assert_eq!(v.iter().map(|p| p.1).collect::<Vec<_>>(), indices(&case["stable_sort"]), "stable_sort case {k} (n {})", keys.len());
    }
}

fn check_random<E: af_core::cxx_random::Engine>(mut e: E, case: &Value) {
    use af_core::cxx_random::*;
    let name = format!("{} seed {}", case["engine"].as_str().unwrap(), case["seed"]);
    let u64s = |v: &Value| -> Vec<u64> { v.as_array().unwrap().iter().map(|x| x.as_u64().unwrap()).collect() };
    let raw: Vec<u64> = (0..20).map(|_| e.next_u64()).collect();
    assert_eq!(raw, u64s(&case["raw"]), "{name}: raw");
    for (i, x) in case["real"].as_array().unwrap().iter().enumerate() {
        assert_eq!(uniform_real(&mut e, -2.5, 7.0), x.as_f64().unwrap(), "{name}: real {i}");
    }
    let mut k = 0;
    let ints = case["int"].as_array().unwrap();
    for b in [1i64, 2, 5, 6, 100, 1000003, 4294967295, 4294967296, 9007199254740993] {
        for _ in 0..6 {
            assert_eq!(uniform_i64(&mut e, -3, b), ints[k].as_i64().unwrap(), "{name}: int {k}");
            k += 1;
        }
    }
    let uints = u64s(&case["uint"]);
    k = 0;
    for b in [1u64, 7, u64::MAX] {
        for _ in 0..4 {
            assert_eq!(uniform_u64(&mut e, 0, b), uints[k], "{name}: uint {k}");
            k += 1;
        }
    }
    let mut normal = Normal::new(0.3, 1.7);
    for (i, x) in case["normal"].as_array().unwrap().iter().enumerate() {
        assert_eq!(normal.sample(&mut e), x.as_f64().unwrap(), "{name}: normal {i}");
    }
    for (s, n) in [0usize, 1, 2, 3, 4, 5, 10, 11, 100, 1001, 70000].iter().enumerate() {
        let mut v: Vec<usize> = (0..*n).collect();
        shuffle(&mut v, &mut e);
        let expected: Vec<usize> = u64s(&case["shuffles"][s]).into_iter().map(|x| x as usize).collect();
        assert_eq!(v, expected, "{name}: shuffle of {n}");
    }
    assert_eq!(e.next_u64(), case["after"].as_u64().unwrap(), "{name}: engine state after");
}

#[test]
fn random_numbers_equal_libstdcxx() {
    let g = golden("random.json");
    for case in g["cases"].as_array().unwrap() {
        let seed = case["seed"].as_u64().unwrap();
        if case["engine"] == "mt19937" {
            check_random(af_core::cxx_random::Mt19937::new(seed as u32), case);
        } else {
            check_random(af_core::cxx_random::Mt19937_64::new(seed), case);
        }
    }
}
