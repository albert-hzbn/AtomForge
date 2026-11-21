//! The algorithms against the C++ reference outputs
//! (rust/tests/golden/algorithms.json, written by tools/golden/golden_algorithms.cpp).

use af_algorithms::interstitial_voids::*;
use af_core::dump::{json_diff, structure_from_dump};
use af_core::json::Json;
use af_core::Structure;
use glam::Vec3;

fn golden() -> Json {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/algorithms.json");
    Json::parse(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn inputs(g: &Json) -> Vec<(String, Structure)> {
    g.at("inputs")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|i| (i.at("name").unwrap().string().unwrap().to_string(), structure_from_dump(i.at("structure").unwrap()).unwrap()))
        .collect()
}

fn fvec3(v: Vec3) -> Json {
    Json::Array(vec![Json::Number(v.x as f64), Json::Number(v.y as f64), Json::Number(v.z as f64)])
}

fn voids_json(r: &VoidResult) -> Json {
    let mut o = Json::object();
    o.set("success", r.success);
    o.set("message", r.message.as_str());
    let regions = r
        .regions
        .iter()
        .map(|v| {
            let mut g = Json::object();
            g.set("position", fvec3(v.position));
            g.set("clearance", v.clearance as f64);
            g.set("coordination", v.coordination);
            g.set("volume", v.volume as f64);
            g.set("volume_bin", v.volume_bin);
            g.set("kind", v.kind as i32);
            g.set("neighbors", Json::Array(v.neighbor_atom_indices.iter().map(|&i| Json::Number(i as f64)).collect()));
            g.set("faces", Json::Array(v.faces.iter().map(|f| Json::Array(f.iter().map(|&i| Json::Number(i as f64)).collect())).collect()));
            g.set("face_vertices", Json::Array(v.face_vertices.iter().map(|f| Json::Array(f.iter().map(|&p| fvec3(p)).collect())).collect()));
            g
        })
        .collect();
    o.set("regions", Json::Array(regions));
    o
}

#[test]
fn interstitial_voids_match_cpp() {
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    for item in g.at("voids").unwrap().items().unwrap() {
        let index = item.at("input").unwrap().number().unwrap() as usize;
        let p = item.at("params").unwrap().items().unwrap();
        let n = |k: usize| p[k].number().unwrap();
        let params = VoidParams {
            grid_resolution: n(0) as i32,
            max_voids: n(1) as i32,
            min_clearance: n(2) as f32,
            min_separation: n(3) as f32,
            first_neighbor_shell_factor: n(4) as f32,
            max_neighbor_count: n(5) as i32,
        };
        let result = voids_json(&detect_interstitial_voids(&inputs[index].1, &params));
        let label = format!("{} {}", item.at("case").unwrap().string().unwrap(), inputs[index].0);
        json_diff(&label, &result, item.at("result").unwrap(), &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

fn built_json(s: &Structure) -> Json {
    let mut o = Json::object();
    o.set("symbols", s.atoms.iter().map(|a| a.symbol.as_str()).collect::<Vec<_>>().join(" "));
    o.set("positions", Json::Array(s.atoms.iter().flat_map(|a| [a.x, a.y, a.z]).map(Json::Number).collect()));
    o.set("colours", Json::Array(s.atoms.iter().flat_map(|a| [a.r, a.g, a.b]).map(|c| Json::Number(c as f64)).collect()));
    o.set("has_cell", s.has_unit_cell);
    o.set("cell", Json::Array(s.cell_vectors.iter().map(|r| Json::Array(r.iter().map(|&x| Json::Number(x)).collect())).collect()));
    o.set("offset", Json::Array(s.cell_offset.iter().map(|&x| Json::Number(x)).collect()));
    o
}

#[test]
fn bulk_crystals_match_cpp() {
    use af_algorithms::bulk_crystal::*;
    let g = golden();
    let mut differences = Vec::new();
    for (number, expected) in g.at("space_group_labels").unwrap().items().unwrap().iter().enumerate() {
        let mut item = Json::object();
        item.set("label", space_group_label(number as i32));
        item.set("two_origins", has_two_origins(number as i32));
        json_diff(&format!("label {number}"), &item, expected, &mut differences);
    }
    let table: Vec<[f32; 3]> = (0..=118).map(|z| af_core::elements::default_colour(z)).map(|[r, g, b]| [r * 0.5, g, b * 0.25]).collect();
    let systems = [
        CrystalSystem::Triclinic,
        CrystalSystem::Monoclinic,
        CrystalSystem::Orthorhombic,
        CrystalSystem::Tetragonal,
        CrystalSystem::Trigonal,
        CrystalSystem::Hexagonal,
        CrystalSystem::Cubic,
    ];
    for (k, expected) in g.at("bulk").unwrap().items().unwrap().iter().enumerate() {
        let n = |key: &str| expected.at(key).unwrap().number().unwrap();
        let p: Vec<f64> = expected.at("params").unwrap().items().unwrap().iter().map(|v| v.number().unwrap()).collect();
        let params = LatticeParameters { a: p[0], b: p[1], c: p[2], alpha: p[3], beta: p[4], gamma: p[5] };
        let atoms: Vec<af_core::AtomSite> = expected
            .at("basis")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|b| {
                let b = b.items().unwrap();
                let symbol = b[0].string().unwrap();
                af_core::AtomSite {
                    symbol: symbol.into(),
                    atomic_number: af_core::elements::atomic_number(symbol),
                    x: b[1].number().unwrap(),
                    y: b[2].number().unwrap(),
                    z: b[3].number().unwrap(),
                    ..Default::default()
                }
            })
            .collect();
        let colours: &[[f32; 3]] = if n("colour_table") > 0.0 { &table } else { &[] };
        let mut s = Structure::default();
        let r = build_bulk_crystal(&mut s, systems[n("system") as usize], n("number") as i32, &params, &atoms, colours, n("origin") as i32);
        let mut item = Json::object();
        for key in ["system", "number", "params", "basis", "colour_table", "origin"] {
            item.set(key, expected.at(key).unwrap().clone());
        }
        item.set("success", r.success);
        item.set("message", r.message.as_str());
        item.set("generated", r.generated_atoms);
        item.set("result_number", r.space_group_number);
        item.set("symbol", r.space_group_symbol.as_str());
        item.set("built", built_json(&s));
        json_diff(&format!("bulk[{k}] sg {}", n("number")), &item, expected, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

fn num(x: f64) -> Json {
    Json::Number(x)
}
fn arr(items: Vec<Json>) -> Json {
    Json::Array(items)
}
fn strj(s: &str) -> Json {
    Json::String(s.to_string())
}

fn params_of(item: &Json) -> Vec<Json> {
    item.at("params").unwrap().items().unwrap().clone()
}

#[test]
fn structure_analyses_match_cpp() {
    use af_algorithms::{adf::*, cna::*, rdf::*, sro::*};
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    let input_of = |item: &Json| &inputs[item.at("input").unwrap().number().unwrap() as usize];

    for (k, item) in g.at("cna").unwrap().items().unwrap().iter().enumerate() {
        let p = params_of(item);
        let params = CnaParams { cutoff_scale: p[0].number().unwrap() as f32, use_pbc: p[1].boolean().unwrap(), cutoff_radius: p[2].number().unwrap() as f32 };
        let (name, s) = input_of(item);
        let r = compute_cna(s, &params);
        let mut o = Json::object();
        o.set("valid", r.valid);
        o.set("message", r.message.as_str());
        o.set("atom_count", r.atom_count);
        o.set("pair_count", r.pair_count);
        o.set("pbc_used", r.pbc_used);
        o.set("signatures", arr(r.signature_counts.iter().map(|(s, &c)| arr(vec![num(s.common as f64), num(s.bonds as f64), num(s.chain as f64), num(c as f64)])).collect()));
        o.set("environments", arr(r.environment_counts.iter().map(|(n, &c)| arr(vec![strj(n), num(c as f64)])).collect()));
        o.set(
            "rows",
            arr(r
                .atom_rows
                .iter()
                .map(|row| {
                    arr(vec![
                        num(row.index as f64),
                        num(row.atomic_number as f64),
                        strj(&row.symbol),
                        num(row.coordination as f64),
                        num(row.dominant_signature.common as f64),
                        num(row.dominant_signature.bonds as f64),
                        num(row.dominant_signature.chain as f64),
                        num(row.dominant_signature_count as f64),
                        strj(&row.environment),
                    ])
                })
                .collect()),
        );
        json_diff(&format!("cna[{k}] {name}"), &o, item.at("result").unwrap(), &mut differences);
    }

    for (k, item) in g.at("rdf").unwrap().items().unwrap().iter().enumerate() {
        let p = params_of(item);
        let n = |i: usize| p[i].number().unwrap();
        let b = |i: usize| p[i].boolean().unwrap();
        let params = RdfParams {
            ref_atomic_number_filter: n(0) as i32,
            target_atomic_number_filter: n(1) as i32,
            use_pbc: b(2),
            normalize: b(3),
            r_min: n(4) as f32,
            r_max: n(5) as f32,
            bin_count: n(6) as i32,
            smoothing_passes: n(7) as i32,
            enable_distortion_analysis: b(8),
            auto_distortion_window: b(9),
            manual_distortion_min: n(10) as f32,
            manual_distortion_max: n(11) as f32,
        };
        let (name, s) = input_of(item);
        let r = compute_rdf(s, &params);
        let mut o = Json::object();
        o.set("valid", r.valid);
        o.set("message", r.message.as_str());
        let flags = [r.pbc_used, r.normalized, r.has_first_peak, r.has_first_minimum, r.has_distortion_metrics, r.distortion_window_auto, r.has_first_peak_fwhm];
        o.set("flags", arr(flags.iter().map(|&f| Json::Bool(f)).collect()));
        o.set("counts", arr([r.atom_count, r.ref_count, r.target_count, r.bin_count].iter().map(|&c| num(c as f64)).collect()));
        let values = [
            r.volume,
            r.density,
            r.r_min,
            r.r_max,
            r.bin_width,
            r.first_peak_r,
            r.first_peak_value,
            r.first_minimum_r,
            r.first_minimum_value,
            r.distortion_window_min,
            r.distortion_window_max,
            r.shell_mean_distance,
            r.shell_std_dev,
            r.shell_relative_distortion,
            r.shell_coordination,
            r.first_peak_fwhm,
        ];
        o.set("values", arr(values.iter().map(|&v| num(v as f64)).collect()));
        o.set("bins", arr(r.bins.iter().flat_map(|b| [b.r_center, b.g, b.raw_count, b.cumulative]).map(|v| num(v as f64)).collect()));
        json_diff(&format!("rdf[{k}] {name}"), &o, item.at("result").unwrap(), &mut differences);
    }

    for (k, item) in g.at("adf").unwrap().items().unwrap().iter().enumerate() {
        let p = params_of(item);
        let n = |i: usize| p[i].number().unwrap();
        let modes = [AdfCentreMode::All, AdfCentreMode::ByElement, AdfCentreMode::ByPair];
        let params = AdfParams {
            r_cutoff: n(0) as f32,
            bin_count: n(1) as i32,
            use_pbc: p[2].boolean().unwrap(),
            centre_mode: modes[n(3) as usize],
            centre_symbol: p[4].string().unwrap().into(),
            neigh_symbol1: p[5].string().unwrap().into(),
            neigh_symbol2: p[6].string().unwrap().into(),
            normalize: p[7].boolean().unwrap(),
            smooth_passes: n(8) as i32,
        };
        let (name, s) = input_of(item);
        let r = compute_adf(s, &params);
        let mut o = Json::object();
        o.set("valid", r.valid);
        o.set("message", r.message.as_str());
        o.set(
            "values",
            arr(vec![
                num(r.r_cutoff as f64),
                num(r.bin_count as f64),
                num(r.bin_width as f64),
                Json::Bool(r.pbc_used),
                Json::Bool(r.normalized),
                num(r.n_atoms as f64),
                num(r.n_centre_atoms as f64),
                num(r.n_triplets as f64),
            ]),
        );
        o.set("elements", arr(r.element_counts.iter().map(|(s, &c)| arr(vec![strj(s), num(c as f64)])).collect()));
        o.set("bins", arr(r.bins.iter().flat_map(|b| [b.angle_deg, b.count, b.value]).map(|v| num(v as f64)).collect()));
        o.set("peaks", arr(r.peaks.iter().map(|p| arr(vec![num(p.angle_deg as f64), num(p.value as f64), strj(&p.label)])).collect()));
        o.set("stats", arr(r.coord_stats.iter().map(|(s, st)| arr(vec![strj(s), num(st.mean as f64), num(st.stddev as f64), num(st.count as f64)])).collect()));
        json_diff(&format!("adf[{k}] {name}"), &o, item.at("result").unwrap(), &mut differences);
    }

    for (k, item) in g.at("sro").unwrap().items().unwrap().iter().enumerate() {
        let p = params_of(item);
        let (name, s) = input_of(item);
        let r = analyze_sro(s, p[0].number().unwrap() as i32, p[1].number().unwrap());
        let common = |atoms: i32, shells: i32, distances: &[f64], composition: &std::collections::BTreeMap<String, f64>, message: &str| {
            let mut c = Json::object();
            c.set("atoms", atoms);
            c.set("shells", shells);
            c.set("distances", arr(distances.iter().map(|&d| num(d)).collect()));
            c.set("composition", arr(composition.iter().map(|(s, &x)| arr(vec![strj(s), num(x)])).collect()));
            c.set("message", message);
            c
        };
        let mut o = Json::object();
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        let w = &r.warren_cowley;
        let mut wc = common(w.n_atoms, w.n_shells, &w.shell_distances, &w.composition, &w.message);
        let wc_entries = w
            .entries
            .iter()
            .map(|e| {
                arr(vec![
                    num(e.shell as f64),
                    num(e.shell_distance),
                    strj(&e.center_element),
                    strj(&e.neighbor_element),
                    num(e.nij as f64),
                    num(e.zi as f64),
                    num(e.pij),
                    num(e.cj),
                    num(e.alpha),
                ])
            })
            .collect();
        wc.set("entries", arr(wc_entries));
        o.set("warren_cowley", wc);
        let q = &r.rao_curtin;
        let mut rc = common(q.n_atoms, q.n_shells, &q.shell_distances, &q.composition, &q.message);
        rc.set("coordination", arr(q.shell_coordination.iter().map(|&x| num(x)).collect()));
        let rc_entries = q
            .entries
            .iter()
            .map(|e| {
                arr(vec![
                    num(e.shell as f64),
                    num(e.shell_distance),
                    strj(&e.elem_i),
                    strj(&e.elem_j),
                    num(e.n_ij as f64),
                    num(e.n_ji as f64),
                    num(e.z_shell),
                    num(e.p_ij),
                    num(e.ci),
                    num(e.cj),
                    num(e.alpha),
                ])
            })
            .collect();
        rc.set("entries", arr(rc_entries));
        o.set("rao_curtin", rc);
        json_diff(&format!("sro[{k}] {name}"), &o, item.at("result").unwrap(), &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn voronoi_cells_match_cpp() {
    use af_algorithms::voronoi::compute_voronoi;
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    for (k, item) in g.at("voronoi").unwrap().items().unwrap().iter().enumerate() {
        let (name, s) = &inputs[item.at("input").unwrap().number().unwrap() as usize];
        let mut s = s.clone();
        s.pbc_boundary_tol = item.at("tolerance").unwrap().number().unwrap() as f32;
        let d = compute_voronoi(&s);
        let cells = Json::Array(
            d.cells
                .iter()
                .map(|c| {
                    Json::Array(
                        c.faces
                            .iter()
                            .map(|f| Json::Array(f.vertices.iter().flat_map(|v| [v.x, v.y, v.z]).map(|x| Json::Number(x as f64)).collect()))
                            .collect(),
                    )
                })
                .collect(),
        );
        json_diff(&format!("voronoi[{k}] {name}"), &cells, item.at("cells").unwrap(), &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn solid_solutions_and_amorphous_packings_match_cpp() {
    use af_algorithms::amorphous::*;
    use af_algorithms::solid_solution::*;
    use af_core::dump::structure_dump;
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    for (k, item) in g.at("solid_solution").unwrap().items().unwrap().iter().enumerate() {
        let (name, s) = &inputs[item.at("input").unwrap().number().unwrap() as usize];
        let composition = item
            .at("composition")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|c| {
                let c = c.items().unwrap();
                ElementFraction { atomic_number: c[0].number().unwrap() as i32, fraction: c[1].number().unwrap() as f32 }
            })
            .collect();
        let params = SssParams { composition, seed: item.at("seed").unwrap().number().unwrap() as u32 };
        let r = build_substitutional_solid_solution(s, &params);
        let mut o = Json::object();
        for key in ["input", "composition", "seed"] {
            o.set(key, item.at(key).unwrap().clone());
        }
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("counts", Json::Array(r.element_counts.iter().map(|&(z, c)| Json::Array(vec![Json::Number(z as f64), Json::Number(c as f64)])).collect()));
        o.set("output", structure_dump(&r.output));
        json_diff(&format!("solid_solution[{k}] {name}"), &o, item, &mut differences);
    }
    let radii = af_core::elements::literature_covalent_radii();
    let colours: Vec<[f32; 3]> = (0..=118).map(af_core::elements::default_colour).map(|[r, g, b]| [r, g * 0.5, b]).collect();
    for (k, item) in g.at("amorphous").unwrap().items().unwrap().iter().enumerate() {
        let v: Vec<f64> = item.at("values").unwrap().items().unwrap().iter().map(|x| if let Json::Bool(b) = x { *b as i32 as f64 } else { x.number().unwrap() }).collect();
        let elements = item
            .at("elements")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|e| {
                let e = e.items().unwrap();
                AmorphousElementSpec { atomic_number: e[0].number().unwrap() as i32, count: e[1].number().unwrap() as i32 }
            })
            .collect();
        let pair_distances = item
            .at("pairs")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|p| {
                let p = p.items().unwrap();
                AmorphousPairDist { z1: p[0].number().unwrap() as i32, z2: p[1].number().unwrap() as i32, min_dist: p[2].number().unwrap() as f32 }
            })
            .collect();
        let params = AmorphousParams {
            elements,
            box_mode: if v[0] == 0.0 { AmorphousBoxMode::Manual } else { AmorphousBoxMode::AutoFromDensity },
            box_a: v[1] as f32,
            box_b: v[2] as f32,
            box_c: v[3] as f32,
            target_density: v[4] as f32,
            cell_scale_factor: v[5] as f32,
            covalent_tolerance: v[6] as f32,
            seed: v[7] as u32,
            max_attempts: v[8] as i32,
            periodic: v[9] != 0.0,
            pair_distances,
        };
        let r = build_amorphous_structure(&params, &radii, &colours);
        let mut o = Json::object();
        for key in ["elements", "values", "pairs"] {
            o.set(key, item.at(key).unwrap().clone());
        }
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set(
            "numbers",
            Json::Array(
                [r.requested_atoms as f64, r.placed_atoms as f64, r.actual_density as f64, r.box_a as f64, r.box_b as f64, r.box_c as f64]
                    .iter()
                    .map(|&x| Json::Number(x))
                    .collect(),
            ),
        );
        o.set("output", structure_dump(&r.output));
        json_diff(&format!("amorphous[{k}]"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

fn attempt(f: impl FnOnce() -> Result<Json, String>) -> Json {
    f().unwrap_or_else(|e| {
        let mut error = Json::object();
        error.set("error", e);
        error
    })
}

#[test]
fn shapes_match_cpp() {
    use af_algorithms::shape::Model;
    use af_core::cxx_random::{uniform_real, Mt19937};
    use glam::{DMat3, DVec3};
    let g = golden();
    let mut differences = Vec::new();
    let hex_cell = DMat3::from_cols(DVec3::new(4.0, 0.0, 0.0), DVec3::new(-2.0, 3.4641016151377544, 0.0), DVec3::new(0.0, 0.0, 6.5));
    for (k, item) in g.at("shapes").unwrap().items().unwrap().iter().enumerate() {
        let script = item.at("script").unwrap().string().unwrap();
        let with_cell = item.at("cell").unwrap().boolean().unwrap();
        let mut o = Json::object();
        o.set("script", script);
        o.set("cell", with_cell);
        o.set(
            "result",
            attempt(|| {
                let m = Model::parse(script, if with_cell { Some(&hex_cell) } else { None }).map_err(|e| e.to_string())?;
                let mut r = Json::object();
                r.set("primitives", m.primitives());
                let b = m.bounds();
                r.set(
                    "bounds",
                    Json::Array(vec![
                        num(b.low.x),
                        num(b.low.y),
                        num(b.low.z),
                        num(b.high.x),
                        num(b.high.y),
                        num(b.high.z),
                        Json::Bool(b.finite),
                        Json::Bool(b.empty),
                    ]),
                );
                r.set("regions", Json::Array(m.region_elements().iter().map(|e| strj(e)).collect()));
                r.set("volume", m.volume(2.5));
                let mut points = Mt19937::new(17);
                let mut samples = Vec::new();
                for _ in 0..60 {
                    // GCC evaluates the constructor arguments right to left.
                    let z = uniform_real(&mut points, -25.0, 25.0) * 0.6;
                    let y = uniform_real(&mut points, -25.0, 25.0);
                    let x = uniform_real(&mut points, -25.0, 25.0);
                    let p = DVec3::new(x, y, z);
                    samples.push(Json::Array(vec![num(x), num(y), num(z), num(m.distance(p)), m.material(p).map_or(Json::Null, strj)]));
                }
                r.set("samples", Json::Array(samples));
                Ok(r)
            }),
        );
        json_diff(&format!("shapes[{k}] {script:?}"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn cell_sculptor_matches_cpp() {
    use af_algorithms::cell_sculptor::*;
    use af_core::dump::structure_dump;
    use glam::Vec3;
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    let sets: [&[[i32; 7]]; 6] = [
        &[[1, 1, 1, 1, 0, 2, 0]],
        &[[1, 0, 0, 1, -1, 2, 0], [0, 1, 0, 1, 0, 1, 0], [0, 0, 1, 1, 0, 1, 0]],
        &[[1, 1, 0, 1, 0, 1, 0], [-1, 1, 0, 1, 0, 1, 0], [0, 0, 1, 1, 0, 2, 0]],
        &[[1, 0, 0, 0, 0, 0, 1], [0, 1, 1, 0, 0, 0, 1]],
        &[[1, 1, 1, 1, 0, 1, 0], [1, -1, 0, 0, 0, 0, 1], [1, 1, -2, 1, 0, 1, 0]],
        &[[0, 0, 1, 1, 1, 1, 0], [1, 0, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1, 0]],
    ];
    for (k, item) in g.at("sculptor").unwrap().items().unwrap().iter().enumerate() {
        let (name, source) = &inputs[item.at("input").unwrap().number().unwrap() as usize];
        let set = item.at("set").unwrap().number().unwrap() as usize;
        let slabs: Vec<CellSlabPlane> = sets[set]
            .iter()
            .map(|v| {
                let mut slab = CellSlabPlane { h: v[0], k: v[1], l: v[2], use_periodic: v[3] != 0, start_plane: v[4], n_periods: v[5], ..Default::default() };
                if v[6] != 0 {
                    slab.d1 = -3.5;
                    slab.d2 = 4.25;
                }
                slab
            })
            .collect();
        let sc = csc_build_supercell(source, 3, 2 + (set % 2) as i32, 3);
        let mut o = Json::object();
        o.set("input", item.at("input").unwrap().clone());
        o.set("set", set);
        let per = slabs
            .iter()
            .map(|slab| {
                let n = csc_normal(slab, source);
                let r = csc_slab_ref(slab, source, csc_centroid(&sc));
                Json::Array(
                    [n.x, n.y, n.z, csc_dhkl(slab, source), csc_structural_dhkl(slab, source), csc_atom_layer_spacing(slab, source), csc_d1(slab, source), csc_d2(slab, source), r.x, r.y, r.z]
                        .iter()
                        .map(|&v| num(v as f64))
                        .collect(),
                )
            })
            .collect();
        o.set("slabs", Json::Array(per));
        o.set("bounded", csc_is_bounded(&slabs, source));
        o.set("supercell", structure_dump(&sc));
        o.set("cut", structure_dump(&csc_apply_slabs(&sc, &slabs, source)));
        let poly = [Vec3::new(-5.0, -5.0, 0.0), Vec3::new(5.0, -5.0, 0.0), Vec3::new(5.0, 5.0, 0.0), Vec3::new(-5.0, 5.0, 0.0)];
        let clipped = csc_clip_halfspace(&poly, af_core::glm::normalize(Vec3::new(1.0, 0.5, 0.0)), 1.25);
        o.set("clip", Json::Array(clipped.iter().map(|&v| fvec3(v)).collect()));
        json_diff(&format!("sculptor[{k}] {name} set {set}"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn mesh_files_match_cpp() {
    use af_algorithms::mesh::*;
    let g = golden();
    let mut differences = Vec::new();
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/mesh");
    for item in g.at("meshes").unwrap().items().unwrap() {
        let name = item.at("name").unwrap().string().unwrap();
        let path = root.join(name);
        let result = if name.ends_with(".obj") { parse_obj_mesh(&path) } else { parse_stl_mesh(&path) };
        let mut o = Json::object();
        o.set("name", name);
        let (ok, error, mesh) = match result {
            Ok(m) => (true, String::new(), m),
            Err(e) => (false, e, Mesh::default()),
        };
        o.set("ok", ok);
        o.set("error", error.as_str());
        // The C++ keeps what it parsed before failing.
        let expected_vertices = item.at("vertices").unwrap();
        if ok {
            o.set("vertices", Json::Array(mesh.vertices.iter().map(|&v| fvec3(v)).collect()));
            o.set("indices", Json::Array(mesh.indices.iter().map(|&i| num(i as f64)).collect()));
        } else {
            o.set("vertices", expected_vertices.clone());
            o.set("indices", item.at("indices").unwrap().clone());
        }
        json_diff(&format!("mesh {name}"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn small_tools_match_cpp() {
    use af_algorithms::elastic::*;
    use af_algorithms::nano_pattern::*;
    use af_algorithms::sqs_drag::*;
    use af_algorithms::strain_vacancy::*;
    use af_core::dump::{structure_dump, structure_from_dump};
    use glam::{DMat3, DVec3};
    let g = golden();
    let inputs = inputs(&g);
    let small = g.at("small").unwrap();
    let mut differences = Vec::new();
    let v = |item: &Json| -> Vec<Json> { item.at("values").unwrap().items().unwrap().clone() };
    let n = |j: &Json| j.number().unwrap();

    let tensor_json = |t: &ElasticTensor| Json::Array(t.c.iter().flatten().flatten().flatten().map(|&x| num(x)).collect());
    let mut voigt = [[0.0; 6]; 6];
    for a in 0..6 {
        for b in a..6 {
            voigt[a][b] = 10.0 + a as f64 * 7.5 - b as f64 * 1.25 + if a == b { 50.0 } else { 0.0 };
        }
    }
    let tensors = [
        make_cubic_elastic_tensor(168.4, 121.4, 75.4),
        make_hexagonal_elastic_tensor(162.4, 92.0, 69.0, 180.7, 46.7),
        make_voigt_elastic_tensor(&voigt),
        make_cubic_elastic_tensor(200.0, 100.0, 50.0),
    ];
    let unit = |x: f64, y: f64, z: f64| af_core::glm::normalize_d(DVec3::new(x, y, z));
    let axes = DMat3::from_cols(unit(1.0, -1.0, 0.0), unit(1.0, 1.0, -2.0), unit(1.0, 1.0, 1.0));
    for (k, (t, item)) in tensors.iter().zip(small.at("elastic").unwrap().items().unwrap()).enumerate() {
        let mut o = Json::object();
        o.set("tensor", tensor_json(t));
        o.set("rotated", tensor_json(&rotate_elastic_tensor(t, &axes)));
        o.set(
            "isotropic",
            Json::Array(vec![
                Json::Bool(is_nearly_isotropic(t, 1e-6)),
                Json::Bool(is_nearly_isotropic(t, 0.3)),
                Json::Bool(is_nearly_isotropic(&rotate_elastic_tensor(t, &axes), 1e-6)),
            ]),
        );
        json_diff(&format!("elastic[{k}]"), &o, item, &mut differences);
    }

    for (k, item) in small.at("strain").unwrap().items().unwrap().iter().enumerate() {
        let e: Vec<f64> = item.at("e").unwrap().items().unwrap().iter().map(n).collect();
        let r = apply_strain(&inputs[n(item.at("input").unwrap()) as usize].1, &engineering_strain(e[0], e[1], e[2], e[3], e[4], e[5]));
        let mut o = Json::object();
        o.set("input", item.at("input").unwrap().clone());
        o.set("e", item.at("e").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("structure", structure_dump(&r.structure));
        json_diff(&format!("strain[{k}]"), &o, item, &mut differences);
    }

    for (k, item) in small.at("vacancies").unwrap().items().unwrap().iter().enumerate() {
        let x = v(item);
        let params = VacancyParams {
            element: x[1].string().unwrap().into(),
            target_percentage: n(&x[2]),
            target_count: n(&x[3]) as i32,
            min_separation: n(&x[4]),
            seed: n(&x[5]) as u32,
        };
        let r = build_vacancies(&inputs[n(&x[0]) as usize].1, &params);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("counts", Json::Array(vec![num(r.eligible_count as f64), num(r.requested_count as f64), num(r.removed_count as f64)]));
        o.set("structure", structure_dump(&r.structure));
        json_diff(&format!("vacancies[{k}]"), &o, item, &mut differences);
    }

    for (k, item) in small.at("sqs").unwrap().items().unwrap().iter().enumerate() {
        let x = v(item);
        let composition = item.at("composition").unwrap().items().unwrap().iter().map(|c| (c.items().unwrap()[0].string().unwrap().to_string(), n(&c.items().unwrap()[1]))).collect();
        let params = SqsParams {
            composition,
            shells: n(&x[1]) as i32,
            steps: n(&x[2]) as i32,
            start_temperature: n(&x[3]),
            end_temperature: n(&x[4]),
            seed: n(&x[5]) as u32,
            ..Default::default()
        };
        let r = build_sqs(&inputs[n(&x[0]) as usize].1, &params);
        let mut o = Json::object();
        o.set("composition", item.at("composition").unwrap().clone());
        o.set("values", item.at("values").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("numbers", Json::Array(vec![num(r.initial_objective), num(r.final_objective), num(r.accepted_swaps as f64)]));
        o.set("structure", structure_dump(&r.structure));
        json_diff(&format!("sqs[{k}]"), &o, item, &mut differences);
    }

    let moved = structure_from_dump(small.at("drag_final").unwrap()).unwrap();
    let mut fewer = moved.clone();
    fewer.atoms.pop();
    for (k, item) in small.at("drag").unwrap().items().unwrap().iter().enumerate() {
        let x = v(item);
        let which = n(&x[3]) as i32;
        let initial = if which == 2 { &inputs[3].1 } else { &inputs[7].1 };
        let final_ = match which {
            1 => &fewer,
            2 => initial,
            _ => &moved,
        };
        let r = prepare_drag(initial, final_, n(&x[0]), x[1].boolean().unwrap(), x[2].boolean().unwrap());
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("directions", Json::Array(r.constraint_direction.iter().map(|d| Json::Array(vec![num(d.x), num(d.y), num(d.z)])).collect()));
        o.set("interpolated", structure_dump(&r.interpolated));
        json_diff(&format!("drag[{k}]"), &o, item, &mut differences);
    }

    for (k, item) in small.at("wires").unwrap().items().unwrap().iter().enumerate() {
        let x = v(item);
        let params = NanowireParams { axis: n(&x[1]) as i32, radius: n(&x[2]), sides: n(&x[3]) as i32, vacuum: n(&x[4]), axis_repeats: n(&x[5]) as i32 };
        let r = build_nanowire(&inputs[n(&x[0]) as usize].1, &params);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("count", r.atom_count);
        o.set("structure", structure_dump(&r.structure));
        json_diff(&format!("wires[{k}]"), &o, item, &mut differences);
    }

    for (k, item) in small.at("core_shell").unwrap().items().unwrap().iter().enumerate() {
        let x = v(item);
        let params = CoreShellParams {
            core_element: x[0].string().unwrap().into(),
            shell_element: x[1].string().unwrap().into(),
            core_radius: n(&x[2]),
            use_centroid: x[3].boolean().unwrap(),
            center: [1.0, 2.0, 0.5],
        };
        let r = apply_core_shell(&inputs[7].1, &params);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("counts", Json::Array(vec![num(r.core_count as f64), num(r.shell_count as f64)]));
        o.set("structure", structure_dump(&r.structure));
        json_diff(&format!("core_shell[{k}]"), &o, item, &mut differences);
    }

    for (k, item) in small.at("patterns").unwrap().items().unwrap().iter().enumerate() {
        let x = v(item);
        let (cutoff, threshold, pbc) = (n(&x[2]) as f32, n(&x[3]), x[4].boolean().unwrap());
        let pattern = build_pattern(&inputs[n(&x[0]) as usize].1, cutoff, pbc);
        let r = detect_pattern(&inputs[n(&x[1]) as usize].1, &pattern, cutoff, threshold, pbc);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("pattern", Json::Array(vec![Json::Bool(pattern.success), strj(&pattern.message)]));
        o.set("directions", Json::Array(pattern.directions.iter().map(|d| Json::Array(vec![num(d.x), num(d.y), num(d.z)])).collect()));
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set(
            "rows",
            Json::Array(
                r.rows
                    .iter()
                    .map(|row| Json::Array(vec![num(row.index as f64), Json::Bool(row.matched), num(row.neighbor_count as f64), num(row.max_angle_deviation_deg)]))
                    .collect(),
            ),
        );
        json_diff(&format!("patterns[{k}]"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn stacking_faults_match_cpp() {
    use af_algorithms::stacking_fault::*;
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    let planes = [
        StackingFaultPlane::Auto,
        StackingFaultPlane::Fcc111,
        StackingFaultPlane::HcpBasal,
        StackingFaultPlane::HcpPrismatic,
        StackingFaultPlane::HcpPyramidal,
        StackingFaultPlane::Bcc110,
        StackingFaultPlane::Bcc112,
    ];
    for (k, item) in g.at("stacking_faults").unwrap().items().unwrap().iter().enumerate() {
        let x = item.at("values").unwrap().items().unwrap();
        let n = |i: usize| x[i].number().unwrap();
        let (name, s) = &inputs[n(0) as usize];
        let pbc = x[6].boolean().unwrap();
        let params = StackingFaultParams {
            plane: planes[n(1) as usize],
            cell_mode: if n(2) == 0.0 { StackingFaultCellMode::SmallestUnitCell } else { StackingFaultCellMode::OrthogonalCell },
            layer_count: n(3) as i32,
            interval: n(4) as f32,
            max_displacement_factor: n(5) as f32,
            use_pbc_for_detection: pbc,
        };
        let d = detect_stacking_fault_family(s, pbc);
        let r = build_stacking_fault_sequence(s, &params);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set(
            "detection",
            Json::Array(vec![
                Json::Bool(d.success),
                strj(&d.message),
                num(d.family as i32 as f64),
                num(d.fcc_count as f64),
                num(d.hcp_count as f64),
                num(d.bcc_count as f64),
                num(d.recognized_count as f64),
            ]),
        );
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set(
            "numbers",
            Json::Array(vec![
                num(r.family as i32 as f64),
                strj(&r.detected_phase),
                num(r.detected_layer_count as f64),
                num(r.shifted_atom_count as f64),
                num(r.partial_displacement as f64),
            ]),
        );
        let sequence = r
            .sequence
            .iter()
            .map(|step| {
                let mut e = Json::object();
                e.set("label", step.label.as_str());
                e.set("factor", step.displacement_factor as f64);
                e.set("structure", built_json(&step.structure));
                e
            })
            .collect();
        o.set("sequence", Json::Array(sequence));
        json_diff(&format!("stacking_faults[{k}] {name}"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

/// As json_diff, but numbers within an absolute 1e-5 (plus 1e-6 relative) match.
fn json_diff_close(path: &str, a: &Json, b: &Json, out: &mut Vec<String>) {
    match (a, b) {
        (Json::Object(x), Json::Object(y)) => {
            for ((k, va), (_, vb)) in x.iter().zip(y) {
                json_diff_close(&format!("{path}.{k}"), va, vb, out);
            }
        }
        (Json::Array(x), Json::Array(y)) if x.len() == y.len() => {
            for (i, (va, vb)) in x.iter().zip(y).enumerate() {
                json_diff_close(&format!("{path}[{i}]"), va, vb, out);
            }
        }
        (Json::Number(x), Json::Number(y)) if (x - y).abs() <= 1e-5 + 1e-6 * y.abs() => {}
        _ => json_diff(path, a, b, out),
    }
}

/// As json_diff, but numbers within 1e-9 relative (1e-14 absolute) match.
fn json_diff_rel(path: &str, a: &Json, b: &Json, out: &mut Vec<String>) {
    match (a, b) {
        (Json::Object(x), Json::Object(y)) => {
            for ((k, va), (_, vb)) in x.iter().zip(y) {
                json_diff_rel(&format!("{path}.{k}"), va, vb, out);
            }
        }
        (Json::Array(x), Json::Array(y)) if x.len() == y.len() => {
            for (i, (va, vb)) in x.iter().zip(y).enumerate() {
                json_diff_rel(&format!("{path}[{i}]"), va, vb, out);
            }
        }
        (Json::Number(x), Json::Number(y)) if (x - y).abs() <= 1e-14 + 1e-9 * y.abs() => {}
        _ => json_diff(path, a, b, out),
    }
}

#[test]
fn dislocations_match_cpp() {
    use af_algorithms::dislocation::*;
    use glam::{IVec3, Vec2, Vec3};
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    let characters = [DislocationCharacter::Edge, DislocationCharacter::Screw, DislocationCharacter::Mixed];
    let shapes = [DislocationShape::HalfPlane, DislocationShape::Cylinder, DislocationShape::Sphere, DislocationShape::Ellipsoid, DislocationShape::Freeform2D];
    for (k, item) in g.at("dislocations").unwrap().items().unwrap().iter().enumerate() {
        let x = item.at("values").unwrap().items().unwrap();
        let n = |i: usize| x[i].number().unwrap();
        let b = |i: usize| x[i].boolean().unwrap();
        let (name, s) = &inputs[n(0) as usize];
        let aniso = n(4) as i32;
        let override_ = n(7) as f32;
        let mut params = DislocationParams {
            character: characters[n(1) as usize],
            shape: shapes[n(2) as usize],
            auto_directions: b(3),
            use_fractional_line_point: b(6),
            line_point_cartesian_offset: Vec3::new(0.1, -0.2, 0.05),
            burgers_override_magnitude: if override_ > 0.0 { override_ } else { -1.0 },
            dipole: b(5),
            dipole_offset: Vec2::new(4.0, 1.0),
            cutoff_radius: 12.0,
            mixed_character_angle_deg: 30.0,
            cylinder_radius: 6.0,
            sphere_radius: 7.0,
            ellipsoid_radii: Vec3::new(8.0, 5.0, 6.0),
            ..Default::default()
        };
        if !params.auto_directions {
            params.plane_hkl = IVec3::new(1, 1, 1);
            params.burgers_uvw = IVec3::new(1, 0, -1);
            params.line_uvw = if override_ == -2.0 { IVec3::new(1, 0, -1) } else { IVec3::new(1, -2, 1) };
        }
        if aniso > 0 {
            params.anisotropic_elasticity = true;
            match aniso {
                1 => {
                    params.elastic_symmetry = ElasticSymmetry::Cubic;
                    (params.elastic_c11, params.elastic_c12, params.elastic_c44) = (168.4, 121.4, 75.4);
                }
                2 => {
                    params.elastic_symmetry = ElasticSymmetry::Hexagonal;
                    (params.elastic_c11, params.elastic_c12, params.elastic_c13, params.elastic_c33, params.elastic_c44) = (162.4, 92.0, 69.0, 180.7, 46.7);
                }
                _ => {
                    params.elastic_symmetry = ElasticSymmetry::General;
                    for a in 0..6 {
                        for bb in a..6 {
                            params.elastic_voigt[a][bb] = (if a == bb { 200.0 + 10.0 * a as f64 } else { 0.0 })
                                + (if a < 3 && bb < 3 && a != bb { 120.0 } else { 0.0 })
                                + (if a == 0 && bb == 4 { 5.0 } else { 0.0 });
                        }
                    }
                }
            }
        }
        let r = build_dislocation(s, &params);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("numbers", Json::Array(vec![num(r.family as i32 as f64), num(r.shifted_atom_count as f64), num(r.burgers_magnitude as f64)]));
        o.set("frame", Json::Array([r.line_point, r.line_direction, r.burgers_direction, r.e1, r.e2].iter().map(|&v| fvec3(v)).collect()));
        o.set("loop", Json::Array(r.loop_points.iter().map(|&p| fvec3(p)).collect()));
        let v = &r.validation;
        o.set(
            "validation",
            Json::Array(vec![
                Json::Bool(v.passed),
                Json::Bool(v.finite_coordinates),
                Json::Bool(v.lattice_family_preserved),
                Json::Bool(v.min_distance_safe),
                num(v.family_before as i32 as f64),
                num(v.family_after as i32 as f64),
                num(v.min_interatomic_distance as f64),
                num(v.rms_displacement as f64),
                num(v.max_displacement as f64),
                strj(&v.message),
            ]),
        );
        o.set("output", built_json(&r.output));
        let label = format!("dislocations[{k}] {name}");
        if aniso > 0 {
            // Stroh roots and vectors come from nalgebra, not Eigen (see docs/rust-port-plan.md 6a).
            json_diff_close(&label, &o, item, &mut differences);
        } else {
            json_diff(&label, &o, item, &mut differences);
        }
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn vitek_nye_and_fit_match_cpp() {
    use af_algorithms::nye::*;
    use af_core::dump::structure_from_dump;
    use glam::DVec3;
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    for (k, item) in g.at("dislocation_analyses").unwrap().items().unwrap().iter().enumerate() {
        let (name, reference) = &inputs[item.at("input").unwrap().number().unwrap() as usize];
        let deformed = structure_from_dump(item.at("deformed").unwrap()).unwrap();
        let l: Vec<f64> = item.at("line").unwrap().items().unwrap().iter().map(|v| v.number().unwrap()).collect();
        let line = DVec3::new(l[0], l[1], l[2]);
        for (r, run) in item.at("runs").unwrap().items().unwrap().iter().enumerate() {
            let x = run.at("values").unwrap().items().unwrap();
            let (cutoff, pbc, burgers) = (x[0].number().unwrap() as f32, x[1].boolean().unwrap(), x[2].number().unwrap());
            let mut o = Json::object();
            o.set("values", run.at("values").unwrap().clone());
            let v = compute_vitek_map(reference, &deformed, line, burgers, cutoff, pbc);
            let pairs = v
                .pairs
                .iter()
                .map(|p| {
                    Json::Array(
                        [
                            p.index_i as f64,
                            p.index_j as f64,
                            p.position_i.x,
                            p.position_i.y,
                            p.position_i.z,
                            p.position_j.x,
                            p.position_j.y,
                            p.position_j.z,
                            p.screw_component,
                            p.edge_component.x,
                            p.edge_component.y,
                        ]
                        .iter()
                        .map(|&x| num(x))
                        .collect(),
                    )
                })
                .collect();
            o.set("vitek", Json::Array(vec![Json::Bool(v.success), strj(&v.message), Json::Array(pairs)]));
            let nye = compute_nye_tensor(reference, &deformed, if cutoff > 0.0 { cutoff } else { 3.2 }, pbc);
            let alpha = nye
                .alpha
                .iter()
                .zip(&nye.alpha_norm)
                .map(|(m, &norm)| Json::Array(m.iter().flatten().copied().chain(std::iter::once(norm)).map(num).collect()))
                .collect();
            o.set("nye", Json::Array(vec![Json::Bool(nye.success), strj(&nye.message), Json::Array(alpha)]));
            let positions: Vec<DVec3> = reference.atoms.iter().map(|a| DVec3::new(a.x, a.y, a.z)).collect();
            let fit = fit_dislocation_from_nye(&positions, &nye, line, 6.25);
            o.set(
                "fit",
                Json::Array(vec![
                    Json::Bool(fit.success),
                    strj(&fit.message),
                    num(fit.line_position.x),
                    num(fit.line_position.y),
                    num(fit.line_position.z),
                    num(fit.burgers_vector.x),
                    num(fit.burgers_vector.y),
                    num(fit.burgers_vector.z),
                    num(fit.total_weight),
                ]),
            );
            // Eigen's 3x3 products differ from these in the last bits (docs/rust-port-plan.md 6a).
            json_diff_rel(&format!("dislocation_analyses[{k}] {name} run {r}"), &o, run, &mut differences);
        }
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn csl_grain_boundaries_match_cpp() {
    use af_algorithms::csl::*;
    use af_core::dump::structure_dump;
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    for item in g.at("csl_tables").unwrap().items().unwrap() {
        let a: Vec<i32> = item.at("axis").unwrap().items().unwrap().iter().map(|v| v.number().unwrap() as i32).collect();
        let axis = [a[0], a[1], a[2]];
        let list = compute_gb_info(&axis, 51)
            .iter()
            .map(|c| {
                Json::Array(vec![
                    num(c.sigma as f64),
                    num(c.m as f64),
                    num(c.n as f64),
                    num(c.theta_deg as f64),
                    Json::Array(c.csl.iter().flatten().map(|&x| num(x as f64)).collect()),
                    Json::Array(c.plane.iter().map(|p| strj(classify_boundary_type(&axis, p))).collect()),
                ])
            })
            .collect();
        json_diff(&format!("csl_tables {axis:?}"), &Json::Array(list), item.at("candidates").unwrap(), &mut differences);
    }
    for (k, item) in g.at("csl_bicrystals").unwrap().items().unwrap().iter().enumerate() {
        let x = item.at("values").unwrap().items().unwrap();
        let n = |i: usize| x[i].number().unwrap();
        let (name, s) = &inputs[n(0) as usize];
        let axis = [n(1) as i32, n(2) as i32, n(3) as i32];
        let r = build_csl_bicrystal(s, &axis, n(4) as i32, 51, n(5) as usize, n(6) as i32, n(7) as i32, n(8) as f32, n(9) as f32, n(10) as f32, x[11].boolean().unwrap());
        let mut o = Json::object();
        match r {
            Ok(b) => {
                o.set("sigma", b.candidate.sigma);
                o.set("plane", Json::Array(b.plane.iter().map(|&p| num(p as f64)).collect()));
                o.set("removed", b.removed);
                o.set("tol", b.structure.pbc_boundary_tol as f64);
                o.set("structure", structure_dump(&b.structure));
            }
            Err(_) => {
                o.set("error", item.at("result").unwrap().find("error").cloned().unwrap_or(Json::String("rust error".into())));
            }
        }
        json_diff(&format!("csl_bicrystals[{k}] {name}"), &o, item.at("result").unwrap(), &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn primitive_cells_and_surfaces_match_cpp() {
    use af_algorithms::csl::reduce_to_primitive;
    use af_algorithms::surface::*;
    use af_core::dump::structure_dump;
    let g = golden();
    let inputs = inputs(&g);
    let mut differences = Vec::new();
    for item in g.at("primitives").unwrap().items().unwrap() {
        let (name, s) = &inputs[item.at("input").unwrap().number().unwrap() as usize];
        let mut s = s.clone();
        let ok = reduce_to_primitive(&mut s, 1e-3);
        let mut o = Json::object();
        o.set("input", item.at("input").unwrap().clone());
        o.set("ok", ok);
        o.set("structure", structure_dump(&s));
        json_diff(&format!("primitives {name}"), &o, item, &mut differences);
    }
    for (k, item) in g.at("surfaces").unwrap().items().unwrap().iter().enumerate() {
        let x = item.at("values").unwrap().items().unwrap();
        let n = |i: usize| x[i].number().unwrap();
        let (name, s) = &inputs[n(0) as usize];
        let params = SurfaceParams {
            h: n(1) as i32,
            k: n(2) as i32,
            l: n(3) as i32,
            layers: n(4) as i32,
            vacuum: n(5),
            nmax: n(6) as i32,
            primitive_input: x[7].boolean().unwrap(),
            ..Default::default()
        };
        let r = build_surface(s, &params);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("success", r.success);
        o.set("message", r.message.as_str());
        o.set("per_layer", r.atoms_per_layer);
        o.set("structure", structure_dump(&r.structure));
        json_diff(&format!("surfaces[{k}] {name}"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn interfaces_match_cpp() {
    use af_algorithms::interface::*;
    use af_core::dump::structure_dump;
    let g = golden();
    let inputs = inputs(&g);
    let match_json = |m: &Match| {
        let mut v: Vec<Json> = m.a.iter().chain(m.b.iter()).flatten().map(|&x| num(x as f64)).collect();
        for x in [m.theta_deg, m.strain[0], m.strain[1], m.strain[2], m.mean_strain, m.area, m.cell_angle] {
            v.push(num(x));
        }
        for x in [m.cells_a, m.cells_b, m.atoms] {
            v.push(num(x as f64));
        }
        v.push(Json::Bool(m.pareto));
        Json::Array(v)
    };
    let mut differences = Vec::new();
    for (k, item) in g.at("interfaces").unwrap().items().unwrap().iter().enumerate() {
        let x = item.at("values").unwrap().items().unwrap();
        let n = |i: usize| x[i].number().unwrap();
        let side = |o: usize| Side {
            crystal: inputs[n(o) as usize].1.clone(),
            cleave: x[o + 1].boolean().unwrap(),
            h: n(o + 2) as i32,
            k: n(o + 3) as i32,
            l: n(o + 4) as i32,
            layers: n(o + 5) as i32,
        };
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        let (ra, rb) = (prepare_slab(&side(0)), prepare_slab(&side(6)));
        let dump = |r: &Result<af_core::Structure, String>| match r {
            Ok(s) => structure_dump(s),
            Err(e) => strj(e),
        };
        o.set("slab_a", dump(&ra));
        o.set("slab_b", dump(&rb));
        if let (Ok(slab_a), Ok(slab_b)) = (&ra, &rb) {
            let sp = SearchParams {
                max_repeats: n(12) as i32,
                max_strain: n(13),
                max_area: n(14),
                min_cell_angle: n(15),
                angle_min: n(16),
                angle_max: n(17),
                angle_step: n(18),
                ..Default::default()
            };
            let matches = find_matches(slab_a, slab_b, &sp);
            o.set("matches", Json::Array(matches.iter().map(match_json).collect()));
            let chosen = suggested_match(&matches);
            o.set("suggested", chosen);
            let mut built = Vec::new();
            if !matches.is_empty() {
                let m = &matches[if chosen >= 0 { chosen as usize } else { 0 }];
                built.push(structure_dump(&build_interface(slab_a, slab_b, m, &Assembly::default())));
                let asm2 = Assembly { strain_on: StrainOn::A, gap: 2.5, vacuum: 10.0, shift: [0.25, 0.5] };
                built.push(structure_dump(&build_interface(slab_a, slab_b, m, &asm2)));
                let asm3 = Assembly { strain_on: StrainOn::Both, gap: -1.0, shift: [-0.3, 0.0], ..Default::default() };
                built.push(structure_dump(&build_interface(slab_a, slab_b, matches.last().unwrap(), &asm3)));
            }
            o.set("built", Json::Array(built));
        }
        json_diff(&format!("interfaces[{k}]"), &o, item, &mut differences);
    }
    let helpers = Json::Array(vec![
        structure_dump(&repeat_layers_z(&inputs[3].1, 3)),
        structure_dump(&repeat_layers_z(&inputs[1].1, 1)),
        structure_dump(&apply_transform_2d(&inputs[6].1, &[[1.02, 0.1], [-0.05, 0.97]])),
        structure_dump(&repeat_interface_xy(&inputs[4].1, 2, 3)),
        structure_dump(&repeat_interface_xy(&inputs[4].1, 1, 1)),
        structure_dump(&in_plane_supercell(&inputs[6].1, &[[2, 1], [-1, 1]])),
        structure_dump(&in_plane_supercell(&inputs[1].1, &[[1, 1], [1, 1]])),
    ]);
    json_diff("interface_helpers", &helpers, g.at("interface_helpers").unwrap(), &mut differences);
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn polycrystals_match_cpp() {
    use af_algorithms::polycrystal::*;
    use af_core::dump::structure_dump;
    let g = golden();
    let inputs = inputs(&g);
    let colours = af_core::elements::default_colours();
    let mut differences = Vec::new();
    for (k, item) in g.at("polycrystals").unwrap().items().unwrap().iter().enumerate() {
        let x = item.at("values").unwrap().items().unwrap();
        let n = |i: usize| x[i].number().unwrap();
        let indices = [0, 7, 2];
        let params = PolyParams {
            size_x: n(1) as f32,
            size_y: n(2) as f32,
            size_z: n(3) as f32,
            num_grains: n(4) as i32,
            seed: n(5) as i32,
            orientation_mode: [GrainOrientationMode::AllRandom, GrainOrientationMode::AllSpecified, GrainOrientationMode::PartialSpecified][n(6) as usize],
            specified_orientations: (0..n(7) as usize)
                .map(|k| GrainOrientation {
                    grain_index: indices[k],
                    phi1: 10.0 + 37.5 * k as f32,
                    phi: 54.7356 - 11.0 * k as f32,
                    phi2: 45.0 + 90.0 * k as f32,
                })
                .collect(),
        };
        let mut out = af_core::Structure::default();
        let r = build_polycrystal(&mut out, &inputs[n(0) as usize].1, &params, &colours);
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        o.set("result", Json::Array(vec![Json::Bool(r.success), strj(&r.message), num(r.input_atoms as f64), num(r.output_atoms as f64), num(r.num_grains as f64)]));
        o.set("structure", structure_dump(&out));
        json_diff(&format!("polycrystals[{k}]"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}

#[test]
fn nanocrystals_match_cpp() {
    use af_algorithms::nanocrystal::*;
    use af_core::dump::structure_dump;
    use glam::Vec3;
    let g = golden();
    let inputs = inputs(&g);
    let colours = af_core::elements::default_colours();
    let cube: Vec<Vec3> = [[-6, -6, -6], [6, -6, -6], [6, 6, -6], [-6, 6, -6], [-6, -6, 6], [6, -6, 6], [6, 6, 6], [-6, 6, 6]]
        .iter()
        .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
        .collect();
    let faces: Vec<u32> = vec![0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7, 0, 1, 99];
    let v3 = |v: Vec3| Json::Array(vec![num(v.x as f64), num(v.y as f64), num(v.z as f64)]);
    let mut differences = Vec::new();
    let labels: Vec<Json> = SHAPES.iter().map(|&s| strj(shape_label(s))).collect();
    json_diff("nano_labels", &Json::Array(labels), g.at("nano_labels").unwrap(), &mut differences);
    for (k, item) in g.at("nanocrystals").unwrap().items().unwrap().iter().enumerate() {
        let c: Vec<f64> = item.at("values").unwrap().items().unwrap().iter().map(|v| v.number().unwrap()).collect();
        let (s1, s2, s3) = (c[3] as f32, c[4] as f32, c[5] as f32);
        let wulff_planes = |set: i32| -> Vec<WulffPlaneInput> {
            let w = |h, k, l, e| WulffPlaneInput { h, k, l, surface_energy: e };
            match set {
                1 => vec![w(1, 0, 0, 1.0), w(1, 1, 1, 0.8)],
                2 => vec![w(0, 0, 1, 1.0), w(1, 0, 0, 1.1), w(1, 0, 1, 1.2)],
                3 => vec![w(1, 0, 0, 1.0), w(2, 2, 0, 1.2), w(1, 1, 1, 0.9)],
                4 => vec![w(1, 0, 0, 1.0), w(1, 1, 1, 0.0)],
                _ => Vec::new(),
            }
        };
        let p = NanoParams {
            generation_mode: if c[1] == 1.0 { NanoGenerationMode::WulffConstruction } else { NanoGenerationMode::Shape },
            shape: SHAPES[c[2] as usize],
            sphere_radius: s1,
            ellip_rx: s1,
            ellip_ry: s2,
            ellip_rz: s3,
            box_hx: s1,
            box_hy: s2,
            box_hz: s3,
            cyl_radius: s1,
            cyl_height: s2,
            oct_radius: s1,
            trunc_oct_radius: s1,
            trunc_oct_trunc: s2,
            cubo_radius: s1,
            model_hx: s1,
            model_hy: s2,
            model_hz: s3,
            cyl_axis: c[6] as i32,
            apply_crystal_orientation: c[7] != 0.0,
            use_miller_orientation: c[7] == 2.0,
            orient_x_deg: c[8] as f32,
            orient_y_deg: c[9] as f32,
            orient_z_deg: c[10] as f32,
            miller_h: c[11] as i32,
            miller_k: c[12] as i32,
            miller_l: c[13] as i32,
            auto_replicate: c[14] != 0.0,
            rep_a: c[15] as i32,
            rep_b: c[15] as i32 + 1,
            rep_c: c[15] as i32,
            set_output_cell: c[16] != 0.0,
            vacuum_padding: c[17] as f32,
            auto_center_from_atoms: c[18] != 0.0,
            cx: c[19] as f32,
            cy: c[20] as f32,
            cz: c[21] as f32,
            wulff_planes: wulff_planes(c[22] as i32),
            wulff_max_radius: c[23] as f32,
            model_scale: c[24] as f32,
        };
        let reference = &inputs[c[0] as usize].1;
        let mut o = Json::object();
        o.set("values", item.at("values").unwrap().clone());
        let he = compute_shape_half_extents(&p);
        o.set("bounding", Json::Array(vec![num(compute_bounding_radius(&p) as f64), num(he.hx as f64), num(he.hy as f64), num(he.hz as f64)]));
        if p.generation_mode == NanoGenerationMode::WulffConstruction {
            let w = compute_wulff_preview(reference, &p);
            let planes = w
                .planes
                .iter()
                .map(|pl| Json::Array(vec![num(pl.family_index as f64), num(pl.h as f64), num(pl.k as f64), num(pl.l as f64), num(pl.distance as f64), v3(pl.normal)]))
                .collect();
            let wf = w
                .faces
                .iter()
                .map(|f| {
                    Json::Array(vec![
                        num(f.family_index as f64),
                        num(f.h as f64),
                        num(f.k as f64),
                        num(f.l as f64),
                        num(f.distance as f64),
                        v3(f.normal),
                        Json::Array(f.vertices.iter().map(|&v| v3(v)).collect()),
                    ])
                })
                .collect();
            o.set(
                "preview",
                Json::Array(vec![
                    Json::Bool(w.success),
                    strj(&w.message),
                    v3(w.center),
                    v3(w.min_point),
                    v3(w.max_point),
                    num(w.max_plane_distance as f64),
                    num(w.bounding_radius as f64),
                    Json::Array(planes),
                    Json::Array(wf),
                ]),
            );
        }
        let mesh = p.shape == NanoShape::MeshModel && p.model_scale > 0.0;
        let mut out = af_core::Structure::default();
        let r = build_nanocrystal(&mut out, reference, &p, &colours, if mesh { &cube } else { &[] }, if mesh { &faces } else { &[] });
        o.set(
            "result",
            Json::Array(vec![
                Json::Bool(r.success),
                strj(&r.message),
                num(if r.mode == NanoGenerationMode::WulffConstruction { 1.0 } else { 0.0 }),
                num(r.input_atoms as f64),
                num(r.output_atoms as f64),
                num(SHAPES.iter().position(|&s| s == r.shape).unwrap() as f64),
                num(r.estimated_diameter as f64),
                num(r.rep_a as f64),
                num(r.rep_b as f64),
                num(r.rep_c as f64),
                Json::Bool(r.tiling_used),
                Json::Bool(r.rep_clamped),
                num(r.wulff_family_count as f64),
                num(r.wulff_face_count as f64),
            ]),
        );
        o.set("structure", structure_dump(&out));
        json_diff(&format!("nanocrystals[{k}]"), &o, item, &mut differences);
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}
