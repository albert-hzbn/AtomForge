//! af-electronic against the C++ electronic library: the cases of
//! tools/golden/golden_electronic.cpp (tests/golden/electronic.json and the
//! saved volumes in tests/golden/electronic_saved), run in the same order.
//! Everything is compared exactly except the Wannier90 band, curvature and
//! Chern results, whose dense Hermitian diagonalisation is nalgebra's rather
//! than Eigen's (1e-9 relative plus 1e-10 absolute; docs/rust-port-plan.md 6a).

use std::path::PathBuf;

use af_core::json::Json;
use af_electronic::bader::{bader_on_grid, bader_populations};
use af_electronic::charge::*;
use af_electronic::density::*;
use af_electronic::field::*;
use af_electronic::fourier::*;
use af_electronic::lobster::{read_cohpcar, read_icohplist};
use af_electronic::topology::{betti_curve, betti_numbers};
use af_electronic::volume::*;
use af_electronic::volume_io::{load_volume, save_volume};
use af_electronic::wannier::*;
use glam::{DMat3, DVec2, DVec3, IVec3};

fn vec3(v: DVec3) -> Json {
    Json::Array(vec![v.x.into(), v.y.into(), v.z.into()])
}

fn numbers(values: &[f64]) -> Json {
    Json::Array(values.iter().map(|&v| Json::Number(v)).collect())
}

fn grid_json(g: &Grid) -> Json {
    let mut o = Json::object();
    o.set("shape", Json::Array(g.shape.iter().map(|&v| Json::from(v)).collect()));
    o.set("cell_columns", Json::Array((0..3).flat_map(|c| (0..3).map(move |r| (c, r))).map(|(c, r)| Json::from(g.cell.col(c)[r])).collect()));
    o.set("origin", vec3(g.origin));
    o.set("periodic", g.periodic);
    o.set("name", g.name.as_str());
    o.set("unit", g.unit.as_str());
    o.set("values", numbers(&g.values));
    o
}

fn volume_json(v: &Volume) -> Json {
    let mut o = Json::object();
    o.set(
        "sites",
        Json::Array(v.sites.iter().map(|s| Json::Array(vec![s.number.into(), s.position.x.into(), s.position.y.into(), s.position.z.into()])).collect()),
    );
    o.set("fields", Json::Array(v.fields.iter().map(grid_json).collect()));
    o
}

fn vec2s(list: &[DVec2]) -> Json {
    Json::Array(list.iter().map(|p| Json::Array(vec![p.x.into(), p.y.into()])).collect())
}

fn reflections(list: &[Reflection]) -> Json {
    Json::Array(list.iter().map(|r| Json::Array(vec![r.hkl.x.into(), r.hkl.y.into(), r.hkl.z.into(), r.value.re.into(), r.value.im.into()])).collect())
}

fn grids(list: &[Grid]) -> Json {
    Json::Array(list.iter().map(grid_json).collect())
}

struct Cases(Vec<Json>);

impl Cases {
    fn attempt(&mut self, name: &str, work: impl FnOnce() -> Result<Json, String>) {
        let mut entry = Json::object();
        entry.set("name", name);
        match work() {
            Ok(r) => entry.set("result", r),
            Err(e) => entry.set("error", e),
        }
        self.0.push(entry);
    }
}

fn same(a: f64, b: f64, relative: bool) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    if relative {
        a == b || (a - b).abs() <= 1e-9 * a.abs().max(b.abs()) + 1e-10
    } else {
        a.to_bits() == b.to_bits() || a == b
    }
}

fn diff(path: &str, a: &Json, b: &Json, relative: bool, out: &mut Vec<String>) {
    match (a, b) {
        (Json::Number(x), Json::Null) | (Json::Null, Json::Number(x)) if x.is_nan() || x.is_infinite() => {}
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

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn electronic_matches_cpp() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    std::env::set_current_dir(&root).unwrap();
    let input = "rust/tests/fixtures/electronic/";
    let saved_rust = "rust/target/golden/electronic_saved_rust/";
    std::fs::create_dir_all(saved_rust).unwrap();
    let path = |name: &str| format!("{input}{name}");
    let mut c = Cases(Vec::new());

    // Reading.
    let loads = [
        ("density.cube", "auto", "bohr"),
        ("density.cube", "density", "bohr"),
        ("density.cube", "raw", "angstrom"),
        ("negative_axis.cube", "density", "bohr"),
        ("two_channel.cube", "potential", "bohr"),
        ("orbitals.cube", "auto", "bohr"),
        ("truncated.cube", "auto", "bohr"),
        ("extra.cube", "auto", "bohr"),
        ("D_exponent.cube", "auto", "bohr"),
        ("density.cube", "elf", "bohr"),
        ("density.cube", "charge", "bohr"),
        ("density.cube", "auto", "nm"),
        ("CHGCAR", "auto", "bohr"),
        ("CHGCAR", "raw", "bohr"),
        ("LOCPOT", "auto", "bohr"),
        ("ELFCAR", "auto", "bohr"),
        ("CHGCAR_four", "auto", "bohr"),
        ("CHGCAR_three", "auto", "bohr"),
        ("CHGCAR_traj", "auto", "bohr"),
        ("grid.raw", "auto", "bohr"),
        ("grid.xsf", "auto", "bohr"),
        ("grid.xsf", "density", "bohr"),
        ("nogrid.xsf", "auto", "bohr"),
        ("anim.xsf", "auto", "bohr"),
        ("missing.cube", "auto", "bohr"),
    ];
    for (file, quantity, unit) in loads {
        c.attempt(&format!("load {file} {quantity} {unit}"), || Ok(volume_json(&load_volume(&path(file), quantity, unit)?)));
    }
    let cube_volume = load_volume(&path("density.cube"), "density", "bohr").unwrap();
    let chgcar = load_volume(&path("CHGCAR"), "auto", "bohr").unwrap();
    let locpot = load_volume(&path("LOCPOT"), "auto", "bohr").unwrap();
    let (d, p, s, l) = (&cube_volume.fields[0], &chgcar.fields[0], &chgcar.fields[1], &locpot.fields[0]);
    let cube_sites: Vec<DVec3> = cube_volume.sites.iter().map(|s| s.position).collect();
    let chg_sites: Vec<DVec3> = chgcar.sites.iter().map(|s| s.position).collect();

    // Writing.
    let saves = [
        ("cube", "cube", "density_cube"),
        ("cube", "xsf", "density_xsf"),
        ("chg", "vasp", "chgcar_vasp"),
        ("chg", "cube", "chgcar_cube"),
        ("chg", "xsf", "chgcar_xsf"),
        ("loc", "cube", "locpot_cube"),
        ("cube", "vasp", "error_vasp"),
        ("chg", "png", "error_png"),
    ];
    for (which, format, name) in saves {
        let volume = match which {
            "cube" => &cube_volume,
            "chg" => &chgcar,
            _ => &locpot,
        };
        c.attempt(&format!("save {name}"), || {
            save_volume(volume, &format!("{saved_rust}{name}"), format)?;
            Ok(Json::Bool(true))
        });
    }

    // Grid algebra.
    for op in ["add", "subtract", "multiply", "divide", "power"] {
        c.attempt(&format!("arithmetic {op}"), || Ok(grid_json(&arithmetic(p, s, op)?)));
    }
    c.attempt("arithmetic divide P", || Ok(grid_json(&arithmetic(p, &scale(p, 0.25)?, "divide")?)));
    c.attempt("arithmetic units", || Ok(grid_json(&arithmetic(p, l, "add")?)));
    c.attempt("arithmetic aligned", || Ok(grid_json(&arithmetic(d, p, "multiply")?)));
    c.attempt("scale", || Ok(grid_json(&scale(p, 2.5)?)));
    c.attempt("scale nan", || Ok(grid_json(&scale(p, f64::NAN)?)));
    c.attempt("resample finite", || {
        let target = Grid {
            shape: [5, 4, 3],
            cell: DMat3::from_cols(DVec3::new(1.5, 0.0, 0.0), DVec3::new(0.1, 1.2, 0.0), DVec3::new(0.0, 0.05, 0.9)),
            origin: d.origin + DVec3::new(0.4, 0.5, 0.3),
            values: vec![0.0; 60],
            ..Default::default()
        };
        Ok(grid_json(&resample(d, &target)?))
    });
    c.attempt("resample periodic", || {
        let target = Grid { shape: [5, 6, 4], values: vec![0.0; 120], ..p.clone() };
        Ok(grid_json(&resample(p, &target)?))
    });
    c.attempt("resample outside", || {
        let target = Grid { origin: d.origin + DVec3::new(50.0, 0.0, 0.0), ..d.clone() };
        Ok(grid_json(&resample(d, &target)?))
    });
    c.attempt("periodic endpoints", || {
        let mut finite = p.clone();
        finite.periodic = false;
        finite.shape = [p.shape[0] + 1, p.shape[1] + 1, p.shape[2] + 1];
        finite.values = vec![0.0; (finite.shape[0] * finite.shape[1] * finite.shape[2]) as usize];
        for z in 0..finite.shape[2] {
            for y in 0..finite.shape[1] {
                for x in 0..finite.shape[0] {
                    let i = finite.index(x, y, z);
                    finite.values[i] = p.values[p.index(x % p.shape[0], y % p.shape[1], z % p.shape[2])];
                }
            }
        }
        Ok(grid_json(&periodic_endpoints(&finite, 1e-8)?))
    });
    c.attempt("periodic endpoints mismatch", || Ok(grid_json(&periodic_endpoints(d, 1e-8)?)));
    c.attempt("minimum images", || {
        let skew = DMat3::from_cols(DVec3::new(4.0, 0.0, 0.0), DVec3::new(3.5, 1.2, 0.0), DVec3::new(0.7, -0.4, 2.0));
        let deltas = [DVec3::new(3.9, 0.2, -0.1), DVec3::new(-7.1, 2.5, 1.9), DVec3::new(0.1, 0.1, 0.1), DVec3::new(12.3, -8.8, 4.4)];
        Ok(Json::Array(deltas.iter().map(|&delta| minimum_image(delta, &skew).map(vec3)).collect::<Result<_, _>>()?))
    });

    // Integration.
    c.attempt("integrate D", || Ok(Json::from(integrate(d)?)));
    c.attempt("integrate P", || Ok(Json::from(integrate(p)?)));
    c.attempt("integrate sphere D", || Ok(Json::from(integrate_sphere(d, cube_sites[0], 0.9)?)));
    c.attempt("integrate sphere P", || Ok(Json::from(integrate_sphere(p, chg_sites[1], 1.3)?)));
    c.attempt("integrate sphere bad", || Ok(Json::from(integrate_sphere(p, chg_sites[1], -1.0)?)));
    let basins = |b: Vec<Basin>| Json::Array(b.iter().map(|b| Json::Array(vec![b.integral.into(), b.volume.into()])).collect());
    c.attempt("voronoi D", || Ok(basins(voronoi_integrate(d, &cube_sites)?)));
    c.attempt("voronoi P", || Ok(basins(voronoi_integrate(p, &chg_sites)?)));
    for limit in [0, 3] {
        for which in 0..2 {
            c.attempt(&format!("peaks {which} {limit}"), || {
                let list = peaks(if which == 1 { p } else { d }, limit)?;
                Ok(Json::Array(list.iter().map(|q| Json::Array(vec![q.position.x.into(), q.position.y.into(), q.position.z.into(), q.value.into()])).collect()))
            });
        }
    }

    // Field analysis.
    for which in 0..2 {
        let g = if which == 1 { p } else { d };
        let tag = if which == 1 { " P" } else { " D" };
        c.attempt(&format!("gradient{tag}"), || Ok(grids(&gradient(g)?)));
        c.attempt(&format!("laplacian{tag}"), || Ok(grid_json(&laplacian(g)?)));
        c.attempt(&format!("energy density{tag}"), || Ok(grids(&energy_density(g, 1e-6)?)));
        c.attempt(&format!("rdg{tag}"), || Ok(grid_json(&reduced_density_gradient(g, 1e-6)?)));
        c.attempt(&format!("signed density{tag}"), || Ok(grid_json(&signed_density(g, 1e-6)?)));
        c.attempt(&format!("dori{tag}"), || Ok(grid_json(&dori(g, 1e-6)?)));
        c.attempt(&format!("smooth{tag}"), || Ok(grid_json(&smooth(g, 0.6, 2)?)));
        for axis in 0..3 {
            c.attempt(&format!("planar average{tag}{axis}"), || Ok(vec2s(&planar_average(g, axis)?)));
        }
        c.attempt(&format!("macroscopic{tag}"), || Ok(vec2s(&macroscopic_average(g, 2, 3)?)));
    }
    c.attempt("energy density S", || Ok(grids(&energy_density(s, 1e-6)?)));
    c.attempt("rdg units", || Ok(grid_json(&reduced_density_gradient(l, 1e-6)?)));
    c.attempt("smooth bad", || Ok(grid_json(&smooth(p, 0.5, 60)?)));
    c.attempt("line profile D", || Ok(vec2s(&line_profile(d, d.position(1, 1, 1), d.position(9, 8, 6), 17)?)));
    c.attempt("line profile P", || Ok(vec2s(&line_profile(p, DVec3::new(-1.0, -2.0, 0.5), DVec3::new(7.0, 9.0, 4.0), 25)?)));
    c.attempt("macroscopic bad", || Ok(vec2s(&macroscopic_average(p, 2, 2)?)));
    c.attempt("section and contours", || {
        let slice = section(d, d.position(1, 1, 3), DVec3::new(2.0, 0.1, 0.0), DVec3::new(0.05, 1.5, 0.2), 14, 11)?;
        let mut o = Json::object();
        o.set("section", grid_json(&slice));
        o.set("contours", Json::Array(contours(&slice, 0.05)?.into_iter().map(vec3).collect()));
        Ok(o)
    });
    c.attempt("section P", || Ok(grid_json(&section(p, DVec3::new(0.3, 0.2, 0.1), DVec3::new(8.0, 1.0, 0.0), DVec3::new(-1.0, 7.0, 2.0), 9, 7)?)));
    c.attempt("section bad", || Ok(grid_json(&section(p, DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0), 9, 7)?)));

    // Fourier.
    let mut all = Vec::new();
    for h in -p.shape[0] / 2..=(p.shape[0] - 1) / 2 {
        for k in -p.shape[1] / 2..=(p.shape[1] - 1) / 2 {
            for l in -p.shape[2] / 2..=(p.shape[2] - 1) / 2 {
                all.push(IVec3::new(h, k, l));
            }
        }
    }
    let hkl = |list: &[[i32; 3]]| list.iter().map(|v| IVec3::from_array(*v)).collect::<Vec<_>>();
    c.attempt("structure factors", || Ok(reflections(&structure_factors(p, &hkl(&[[0, 0, 0], [1, 0, 0], [-1, 2, 1], [3, -3, 2]]))?)));
    c.attempt("structure factors range", || Ok(reflections(&structure_factors(p, &hkl(&[[4, 0, 0]]))?)));
    c.attempt("structure factors finite", || Ok(reflections(&structure_factors(d, &hkl(&[[0, 0, 0]]))?)));
    c.attempt("fourier synthesis", || Ok(grid_json(&fourier_synthesis(p, &structure_factors(p, &all)?)?)));
    c.attempt("fourier synthesis partial", || Ok(grid_json(&fourier_synthesis(p, &structure_factors(p, &hkl(&[[1, 0, 0]]))?)?)));
    c.attempt("fourier synthesis duplicate", || Ok(grid_json(&fourier_synthesis(p, &structure_factors(p, &hkl(&[[1, 0, 0], [1, 0, 0]]))?)?)));
    c.attempt("patterson", || Ok(grid_json(&patterson(p)?)));
    c.attempt("patterson odd", || {
        let mut odd = Grid { shape: [5, 6, 3], values: vec![0.0; 90], ..p.clone() };
        for (i, v) in odd.values.iter_mut().enumerate() {
            *v = (0.37 * i as f64).sin() + 1.0;
        }
        Ok(grid_json(&patterson(&odd)?))
    });
    c.attempt("atomic structure factors", || {
        let indices = hkl(&[[1, 0, 0], [1, 1, 0], [2, -1, 1]]);
        let mut f = Vec::new();
        for r in 0..indices.len() {
            for j in 0..chg_sites.len() {
                f.push(Complex::new(29.0 - 3.0 * r as f64 - 10.0 * j as f64, 0.5 * r as f64));
            }
        }
        Ok(reflections(&atomic_structure_factors(p, &chg_sites, &[1.0, 0.8], &[0.5, 1.2], &indices, &f)?))
    });
    let ewald_json = |r: EwaldResult| {
        let mut o = Json::object();
        o.set("potentials", numbers(&r.potentials));
        o.set("energy", r.energy);
        o
    };
    c.attempt("ewald", || {
        let cell = DMat3::from_diagonal(DVec3::splat(5.64));
        let sites: Vec<DVec3> = [[0.0, 0.0, 0.0], [2.82, 2.82, 0.0], [2.82, 0.0, 2.82], [0.0, 2.82, 2.82], [2.82, 0.0, 0.0], [0.0, 2.82, 0.0], [0.0, 0.0, 2.82], [2.82, 2.82, 2.82]]
            .iter()
            .map(|v| DVec3::from_array(*v))
            .collect();
        Ok(ewald_json(ewald(&cell, &sites, &[1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0], 0.6, 9.0, 4.0)?))
    });
    c.attempt("ewald skew", || {
        let cell = DMat3::from_cols(DVec3::new(4.0, 0.0, 0.0), DVec3::new(1.3, 3.7, 0.0), DVec3::new(0.6, -0.5, 4.4));
        let sites = [DVec3::new(0.1, 0.2, 0.3), DVec3::new(2.4, 1.9, 2.0), DVec3::new(1.0, 3.0, 3.5)];
        Ok(ewald_json(ewald(&cell, &sites, &[2.0, -1.5, -0.5], 0.8, 7.0, 5.0)?))
    });
    c.attempt("ewald charged", || Ok(Json::from(ewald(&DMat3::from_diagonal(DVec3::splat(3.0)), &[DVec3::ZERO], &[1.0], 0.5, 5.0, 5.0)?.energy)));

    // Charge analysis.
    let p2 = scale(p, 0.9).unwrap();
    c.attempt("density difference", || Ok(grid_json(&density_difference(p, &p2, 0.5)?)));
    c.attempt("density difference units", || Ok(grid_json(&density_difference(p, l, 1.0)?)));
    let mask_a = threshold_mask(p, 0.5, 3.0).unwrap();
    let mask_b = threshold_mask(s, -0.1, 0.1).unwrap();
    c.attempt("threshold mask", || Ok(grid_json(&mask_a)));
    for op in ["union", "intersection", "difference", "xor", "nand"] {
        c.attempt(&format!("boolean {op}"), || Ok(grid_json(&boolean_mask(&mask_a, &mask_b, op)?)));
    }
    c.attempt("boolean unmasked", || Ok(grid_json(&boolean_mask(p, &mask_b, "union")?)));
    c.attempt("apply mask", || Ok(grid_json(&apply_mask(p, &mask_a)?)));
    c.attempt("split", || Ok(grids(&split_density(s)?)));
    c.attempt("charge summary", || Ok(numbers(&charge_summary(s)?)));
    c.attempt("cumulative P", || Ok(vec2s(&cumulative_charge(p, 2)?)));
    c.attempt("cumulative D", || Ok(vec2s(&cumulative_charge(d, 0)?)));

    // Bader.
    for which in 0..2 {
        let (g, sites) = if which == 1 { (p, &chg_sites) } else { (d, &cube_sites) };
        c.attempt(&format!("bader {}", if which == 1 { "P" } else { "D" }), || {
            let partition = bader_on_grid(g)?;
            let mut o = Json::object();
            o.set("basin", Json::Array(partition.basin.iter().map(|&b| Json::from(b)).collect()));
            o.set(
                "basins",
                Json::Array(
                    partition
                        .basins
                        .iter()
                        .map(|b| Json::Array(vec![b.charge.into(), b.volume.into(), b.maximum.x.into(), b.maximum.y.into(), b.maximum.z.into()]))
                        .collect(),
                ),
            );
            o.set("populations", numbers(&bader_populations(g, sites)?));
            Ok(o)
        });
    }
    c.attempt("bader negative", || Ok(Json::from(bader_on_grid(s)?.basins.len() as i32)));

    // Topology.
    c.attempt("betti curve", || {
        let rows = betti_curve(d, &[0.001, 0.01, 0.05, 0.2, 1.0, 5.0])?;
        Ok(Json::Array(rows.iter().map(|r| Json::Array(vec![r.x.into(), r.y.into(), r.z.into(), r.w.into()])).collect()))
    });
    c.attempt("betti shell", || {
        let mut g = Grid { shape: [9, 9, 9], values: vec![0.0; 729], ..Default::default() };
        for z in 0..9 {
            for y in 0..9 {
                for x in 0..9 {
                    let r = ((x as f64 - 4.0).powi(2) + (y as f64 - 4.0).powi(2) + (z as f64 - 4.0).powi(2)).sqrt();
                    let i = g.index(x, y, z);
                    g.values[i] = if r > 1.5 && r < 3.5 { 1.0 } else { 0.0 };
                }
            }
        }
        let b = betti_numbers(&g, 0.5)?;
        Ok(Json::Array(vec![b.betti0.into(), b.betti1.into(), b.betti2.into()]))
    });
    c.attempt("betti periodic", || Ok(Json::from(betti_numbers(p, 1.0)?.betti0)));

    // LOBSTER.
    for name in ["COHPCAR.lobster", "COHPCAR.spin", "COHPCAR.noheaders", "COHPCAR.badheader", "COHPCAR.nolabel", "COHPCAR.short", "missing"] {
        c.attempt(&format!("cohpcar {name}"), || {
            let data = read_cohpcar(&path(name))?;
            let mut o = Json::object();
            o.set("spin", data.spin_polarized);
            o.set("fermi", data.fermi_energy);
            o.set("energies", numbers(&data.energies));
            o.set("bonds", Json::Array(data.bonds.iter().map(|b| Json::Array(vec![b.atom1.into(), b.atom2.into(), b.length.into()])).collect()));
            let mut curves = Vec::new();
            for spin in 0..data.cohp.len() {
                curves.push(numbers(&data.average_cohp[spin]));
                curves.push(numbers(&data.average_icohp[spin]));
                for b in 0..data.cohp[spin].len() {
                    curves.push(numbers(&data.cohp[spin][b]));
                    curves.push(numbers(&data.icohp[spin][b]));
                }
            }
            o.set("curves", Json::Array(curves));
            Ok(o)
        });
    }
    for name in ["ICOHPLIST.lobster", "ICOHPLIST.spin", "ICOHPLIST.old", "ICOHPLIST.bad", "missing"] {
        c.attempt(&format!("icohplist {name}"), || {
            let list = read_icohplist(&path(name))?;
            let mut o = Json::object();
            o.set("spin", list.spin_polarized);
            o.set(
                "entries",
                Json::Array(
                    list.entries
                        .iter()
                        .map(|e| Json::Array(vec![e.atom1.into(), e.atom2.into(), e.length.into(), e.num_bonds.into(), e.icohp[0].into(), e.icohp[1].into()]))
                        .collect(),
                ),
            );
            Ok(o)
        });
    }

    // Wannier90.
    for name in ["short_hr.dat", "missing_hr.dat"] {
        c.attempt(&format!("wannier {name}"), || Ok(Json::from(read_wannier_hr(&path(name))?.num_wann)));
    }
    let hr = read_wannier_hr(&path("haldane_hr.dat")).unwrap();
    let non_hermitian = read_wannier_hr(&path("nonhermitian_hr.dat")).unwrap();
    c.attempt("wannier read", || {
        let mut o = Json::object();
        o.set("num_wann", hr.num_wann);
        let mut blocks = Vec::new();
        for i in 0..hr.vectors.len() {
            let mut row = vec![hr.vectors[i].x.into(), hr.vectors[i].y.into(), hr.vectors[i].z.into(), hr.degeneracies[i].into()];
            for v in &hr.hamiltonians[i] {
                row.push(v.re.into());
                row.push(v.im.into());
            }
            blocks.push(Json::Array(row));
        }
        o.set("blocks", Json::Array(blocks));
        Ok(o)
    });
    c.attempt("bloch", || Ok(Json::Array(bloch_hamiltonian(&hr, DVec3::new(0.13, 0.27, 0.0))?.iter().flat_map(|v| [Json::from(v.re), Json::from(v.im)]).collect())));
    let path_k = [DVec3::ZERO, DVec3::new(0.5, 0.0, 0.0), DVec3::new(1.0 / 3.0, 1.0 / 3.0, 0.0), DVec3::new(0.2, 0.7, 0.1)];
    c.attempt("bands", || Ok(numbers(&wannier_bands(&hr, &path_k)?)));
    c.attempt("bands non-Hermitian", || Ok(numbers(&wannier_bands(&non_hermitian, &[DVec3::new(0.1, 0.2, 0.3)])?)));
    c.attempt("berry", || Ok(numbers(&berry_curvature(&hr, DVec3::new(0.21, 0.37, 0.0), 0, 1)?)));
    c.attempt("berry plane", || Ok(numbers(&berry_curvature(&hr, DVec3::new(0.21, 0.37, 0.0), 1, 1)?)));
    for band in [0, 1] {
        c.attempt(&format!("chern {band}"), || Ok(Json::from(chern_number(&hr, band, 0, 1, 12, 0.0)?)));
    }
    c.attempt("chern band", || Ok(Json::from(chern_number(&hr, 2, 0, 1, 12, 0.0)?)));

    // Comparison.
    let golden = Json::parse(&std::fs::read_to_string("rust/tests/golden/electronic.json").unwrap()).unwrap();
    let expected = golden.items().unwrap();
    assert_eq!(c.0.len(), expected.len(), "case count");
    let mut differences = Vec::new();
    for (got, want) in c.0.iter().zip(expected) {
        let name = want.at("name").unwrap().string().unwrap();
        let relative = ["bands", "berry", "chern 0", "chern 1"].contains(&name);
        diff(name, got, want, relative, &mut differences);
    }
    for (_, _, name) in saves {
        let cpp = std::fs::read(format!("rust/tests/golden/electronic_saved/{name}")).ok();
        let rust = std::fs::read(format!("{saved_rust}{name}")).ok();
        if cpp != rust {
            differences.push(format!("saved {name} differs"));
        }
    }
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.iter().take(40).cloned().collect::<Vec<_>>().join("\n"));
}
