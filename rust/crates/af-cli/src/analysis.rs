//! `--analyze MODE` and `--convert` (the C++ cli/AnalysisCLI): structure
//! analyses written as CSV, the cell sculptor, drag preparation and format
//! conversion.

use std::collections::{BTreeMap, BTreeSet};

use af_algorithms::adf::{compute_adf, AdfCentreMode, AdfParams};
use af_algorithms::cell_sculptor::{csc_apply_slabs, csc_build_supercell, CellSlabPlane};
use af_algorithms::cna::{compute_cna, CnaParams};
use af_algorithms::interstitial_voids::{detect_interstitial_voids, VoidParams};
use af_algorithms::nano_pattern::{build_pattern, detect_pattern, CrystalPattern};
use af_algorithms::nye::{compute_nye_tensor, compute_vitek_map, fit_dislocation_from_nye};
use af_algorithms::rdf::{compute_rdf, RdfParams};
use af_algorithms::sqs_drag::prepare_drag;
use af_algorithms::sro::analyze_sro;
use af_core::format::stod;
use af_core::structure::Structure;
use af_core::suggest;
use af_io::loader::{load_structure_from_file, save_structure};
use glam::DVec3;

use crate::args::{detect_format, err, g17, out, text_mode, write_text_file, Args, IStream};

const HELP: &str = concat!(
    "Analysis: --analyze MODE --input FILE --output FILE [options]\n",
    "cna: --cutoff R (0 uses --scale 1.18), --no-pbc\n",
    "rdf: --minimum 0 --cutoff 8 --bins 200 --smooth 0 --centre-z 0 --neighbor-z 0 --raw --no-pbc\n",
    "adf: --cutoff 3.5 --bins 180 --smooth 0 --centre SYMBOL --neighbor-a SYMBOL --neighbor-b SYMBOL --raw --no-pbc\n",
    "sro: --shells 3 --tolerance 0.1\n",
    "interstitial: --resolution 14 --limit 0 --clearance 0.6 --separation 0.5\n",
    "nye: --reference FILE --cutoff 3.0 --no-pbc  (dislocation density tensor;\n",
    "     --input is the deformed structure, --reference the perfect one, same\n",
    "     atom count/order in both)\n",
    "vitek: --reference FILE --line \"u v w\" --burgers M --cutoff 0 --no-pbc\n",
    "     (differential-displacement map; --line is the dislocation line/Burgers\n",
    "     direction, --burgers its magnitude in Angstrom -- 0 disables screw-\n",
    "     component wrapping; --cutoff 0 auto-detects the neighbor shell)\n",
    "pattern-init: --cutoff C --no-pbc  (--input is the perfect reference\n",
    "     structure; writes a dx,dy,dz pattern file to --output)\n",
    "pattern-detect: --pattern FILE --cutoff C --angle-threshold 10 --no-pbc\n",
    "     (--input is the possibly-defective structure, --pattern a file\n",
    "     written by pattern-init)\n",
    "drag: --final FILE --zeta Z --clip-displacement --no-pbc  (prepares a\n",
    "     constrained-minimization/migration-barrier calculation; --input is\n",
    "     the initial configuration, --final the final one, same atom count/\n",
    "     order; writes the interpolated structure to --output and the\n",
    "     per-atom reaction-coordinate direction to --output + '.constraint.csv')\n",
    "fit-dislocation: --reference FILE --line \"u v w\" --cutoff 3.0 --area 1.0\n",
    "     --no-pbc  (recovers a dislocation's position/Burgers vector from a\n",
    "     measured field, by Nye-tensor centroid/integral; --input is the\n",
    "     deformed structure, --reference the perfect one; --area is the\n",
    "     cross-sectional area per atom (A^2) used to scale the Burgers vector\n",
    "     estimate)\n",
    "sculpt: --slabs \"h k l lower upper periodic;...\" --nx 1 --ny 1 --nz 1\n",
    "Lengths are Angstrom. Periodic sculpt uses integer start/count bounds.\n",
    "Conversion: --convert --input FILE --output FILE [--format FORMAT]\n",
    "CSV is used for analysis results; sculpt writes VASP.\n",
);

pub fn run(argv: &Args) -> i32 {
    if argv.has("--help") {
        out(HELP);
        return 0;
    }
    match analyse(argv) {
        Ok(()) => 0,
        Err(e) => {
            err(&format!("Error: {e}\n"));
            1
        }
    }
}

fn f(v: f32) -> String {
    g17(f64::from(v))
}

fn load(path: &str) -> Result<Structure, String> {
    load_structure_from_file(path)
}

/// Lines as std::getline reads a text-mode file (no empty last line).
fn text_lines(text: &str) -> Vec<String> {
    let text = text.replace("\r\n", "\n");
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    if lines.last().is_some_and(|l| l.is_empty()) {
        lines.pop();
    }
    lines
}

fn analyse(argv: &Args) -> Result<(), String> {
    let mut args: BTreeMap<String, String> = BTreeMap::new();
    let a = &argv.0;
    let mut i = 0;
    while i < a.len() {
        let key = a[i].clone();
        if matches!(key.as_str(), "--convert" | "--no-pbc" | "--raw" | "--clip-displacement") {
            if args.insert(key.clone(), "1".into()).is_some() {
                return Err(format!("Repeated option {key}"));
            }
        } else {
            if i + 1 == a.len() {
                return Err(format!("Missing value for {key}"));
            }
            i += 1;
            if args.contains_key(&key) {
                return Err(format!("Repeated option {key}"));
            }
            args.insert(key, a[i].clone());
        }
        i += 1;
    }
    let number = |key: &str, fallback: f64| -> Result<f64, String> {
        let Some(text) = args.get(key) else { return Ok(fallback) };
        match stod(text) {
            Ok((value, used)) if used > 0 && used == text.len() && value.is_finite() => Ok(value),
            _ => Err(format!("{key} needs a finite number, not \"{text}\"")),
        }
    };
    let integer = |key: &str, fallback: i32, low: i32, high: i32| -> Result<i32, String> {
        let value = number(key, f64::from(fallback))?;
        if value != value.floor() || value < f64::from(low) || value > f64::from(high) {
            return Err(format!("{key} must be a whole number from {low} to {high}"));
        }
        Ok(value as i32)
    };
    if !args.contains_key("--input") || !args.contains_key("--output") {
        return Err("--input FILE and --output FILE are required".into());
    }
    let options: BTreeMap<&str, BTreeSet<&str>> = [
        ("convert", vec!["--convert", "--format"]),
        ("cna", vec!["--cutoff", "--scale", "--no-pbc"]),
        ("rdf", vec!["--cutoff", "--minimum", "--bins", "--smooth", "--centre-z", "--neighbor-z", "--raw", "--no-pbc"]),
        ("adf", vec!["--cutoff", "--bins", "--smooth", "--centre", "--neighbor-a", "--neighbor-b", "--raw", "--no-pbc"]),
        ("sro", vec!["--shells", "--tolerance"]),
        ("interstitial", vec!["--resolution", "--limit", "--clearance", "--separation"]),
        ("nye", vec!["--reference", "--cutoff", "--no-pbc"]),
        ("vitek", vec!["--reference", "--line", "--burgers", "--cutoff", "--no-pbc"]),
        ("pattern-init", vec!["--cutoff", "--no-pbc"]),
        ("pattern-detect", vec!["--pattern", "--cutoff", "--angle-threshold", "--no-pbc"]),
        ("drag", vec!["--final", "--zeta", "--clip-displacement", "--no-pbc"]),
        ("fit-dislocation", vec!["--reference", "--line", "--cutoff", "--area", "--no-pbc"]),
        ("sculpt", vec!["--slabs", "--nx", "--ny", "--nz"]),
    ]
    .into_iter()
    .map(|(k, v)| (k, v.into_iter().collect()))
    .collect();
    let selected: String = if args.contains_key("--convert") {
        "convert".into()
    } else {
        args.get("--analyze").cloned().ok_or("map::at")?
    };
    let Some(allowed) = options.get(selected.as_str()) else {
        return Err(format!("Unknown analysis: {selected}"));
    };
    for key in args.keys() {
        if key != "--input" && key != "--output" && !(key == "--analyze" && selected != "convert") && !allowed.contains(key.as_str()) {
            let names: Vec<&str> = allowed.iter().copied().filter(|o| *o != "--convert").collect();
            return Err(format!(
                "Unsupported option for {selected}: {key}{}{}",
                suggest::did_you_mean(key, &names),
                if names.is_empty() { String::new() } else { format!("; it takes {}", suggest::list(&names)) }
            ));
        }
    }
    let structure = load(&args["--input"])?;
    if !structure.load_note.is_empty() {
        err(&format!("Note: {}\n", structure.load_note));
    }
    let output = args["--output"].clone();
    if args.contains_key("--convert") {
        let mut format = args.get("--format").cloned().unwrap_or_else(|| detect_format(&output));
        if format == "pwi" {
            format = "pwscf".into();
        }
        if format == "extxyz" {
            let mut text = format!("{}\n", structure.atoms.len());
            if structure.has_unit_cell {
                text += "Lattice=\"";
                for row in &structure.cell_vectors {
                    for &v in row {
                        text += &format!("{} ", g17(v));
                    }
                }
                text += "\" ";
            }
            text += "Properties=species:S:1:pos:R:3\n";
            for atom in &structure.atoms {
                text += &format!("{} {} {} {}\n", atom.symbol, g17(atom.x), g17(atom.y), g17(atom.z));
            }
            return write_text_file(&output, &text).map_err(|_| "Cannot write extended XYZ".to_string());
        }
        return save_structure(&structure, &output, &format).map_err(|e| if e.is_empty() { "Conversion failed".into() } else { e });
    }
    let mode = args["--analyze"].clone();
    let mut text = String::new();
    if mode == "sculpt" {
        let mut slabs = Vec::new();
        let source = args.get("--slabs").cloned().unwrap_or_default();
        // std::getline(..., ';'): no empty piece after a final separator.
        let mut rows: Vec<&str> = source.split(';').collect();
        if rows.last() == Some(&"") {
            rows.pop();
        }
        for row in rows {
            let mut v = IStream::new(row);
            let parsed = (|| Some((v.int()?, v.int()?, v.int()?, v.float()?, v.float()?, v.int()?)))();
            let bad = "Slabs use h k l lower upper periodic, separated by semicolons";
            let Some((h, k, l, d1, d2, periodic)) = parsed else { return Err(bad.into()) };
            if (h == 0 && k == 0 && l == 0)
                || !d1.is_finite()
                || !d2.is_finite()
                || !(-1000..=1000).contains(&h)
                || !(-1000..=1000).contains(&k)
                || !(-1000..=1000).contains(&l)
                || (periodic != 0 && periodic != 1)
            {
                return Err(bad.into());
            }
            let mut slab = CellSlabPlane { h, k, l, d1, d2, use_periodic: periodic != 0, ..Default::default() };
            if slab.use_periodic {
                if !structure.has_unit_cell || d1 != d1.floor() || d2 != d2.floor() || d1.abs() > 10000.0 || d2 < 1.0 || d2 > 10000.0 {
                    return Err("Periodic slabs need a cell, integer start and positive period count".into());
                }
                slab.start_plane = d1 as i32;
                slab.n_periods = d2 as i32;
            } else if d2 < d1 {
                return Err("Slab upper bound is below its lower bound".into());
            }
            slabs.push(slab);
        }
        if slabs.is_empty() {
            return Err("At least one slab is required".into());
        }
        let (nx, ny, nz) = (integer("--nx", 1, 1, 100)?, integer("--ny", 1, 1, 100)?, integer("--nz", 1, 1, 100)?);
        if (nx as u64) * (ny as u64) * (nz as u64) * (structure.atoms.len() as u64) > 10000000 {
            return Err("Sculpted supercell exceeds 10 million atoms".into());
        }
        let supercell = csc_build_supercell(&structure, nx, ny, nz);
        let result = csc_apply_slabs(&supercell, &slabs, &structure);
        return save_structure(&result, &output, "vasp").map_err(|_| "Cannot save sculpted structure".to_string());
    }
    let periodic = !args.contains_key("--no-pbc");
    if mode == "drag" {
        let Some(final_path) = args.get("--final") else { return Err("drag requires --final FILE (the final configuration)".into()) };
        let final_structure = load(final_path)?;
        let zeta = number("--zeta", 0.5)?;
        let result = prepare_drag(&structure, &final_structure, zeta, args.contains_key("--clip-displacement"), periodic);
        if !result.success {
            return Err(result.message);
        }
        let mut fmt = std::path::Path::new(&output).extension().map(|e| e.to_string_lossy().into_owned()).unwrap_or_default();
        if fmt.is_empty() {
            fmt = "vasp".into();
        }
        save_structure(&result.interpolated, &output, &fmt).map_err(|_| "Cannot save interpolated structure".to_string())?;
        let mut csv = String::from("index,dx,dy,dz\n");
        for (i, d) in result.constraint_direction.iter().enumerate() {
            csv += &format!("{i},{},{},{}\n", g17(d.x), g17(d.y), g17(d.z));
        }
        write_text_file(&format!("{output}.constraint.csv"), &csv).map_err(|_| "Cannot write constraint file".to_string())?;
        out(&format!("{}\n", result.message));
        return Ok(());
    }
    let parse_line = |key: &str| -> Result<DVec3, String> {
        let mut s = IStream::new(&args[key]);
        match (s.double(), s.double(), s.double()) {
            (Some(x), Some(y), Some(z)) => Ok(DVec3::new(x, y, z)),
            _ => Err("Cannot parse --line \"u v w\"".into()),
        }
    };
    match mode.as_str() {
        "cna" => {
            let p = CnaParams { use_pbc: periodic, cutoff_radius: number("--cutoff", 0.0)? as f32, cutoff_scale: number("--scale", 1.18)? as f32 };
            let result = compute_cna(&structure, &p);
            if !result.valid {
                return Err(result.message);
            }
            text += "index,symbol,coordination,common,bonds,chain,environment\n";
            for row in &result.atom_rows {
                let s = &row.dominant_signature;
                text += &format!("{},{},{},{},{},{},{}\n", row.index, row.symbol, row.coordination, s.common, s.bonds, s.chain, row.environment);
            }
        }
        "rdf" => {
            let p = RdfParams {
                use_pbc: periodic,
                normalize: !args.contains_key("--raw"),
                r_min: number("--minimum", 0.0)? as f32,
                r_max: number("--cutoff", 8.0)? as f32,
                bin_count: integer("--bins", 200, 1, 100000)?,
                ref_atomic_number_filter: integer("--centre-z", 0, 0, 118)?,
                target_atomic_number_filter: integer("--neighbor-z", 0, 0, 118)?,
                smoothing_passes: integer("--smooth", 0, 0, 1000)?,
                ..Default::default()
            };
            let result = compute_rdf(&structure, &p);
            if !result.valid {
                return Err(result.message);
            }
            text += "radius_A,g,raw_count,coordination\n";
            for row in &result.bins {
                text += &format!("{},{},{},{}\n", f(row.r_center), f(row.g), f(row.raw_count), f(row.cumulative));
            }
        }
        "adf" => {
            let mut p = AdfParams {
                use_pbc: periodic,
                normalize: !args.contains_key("--raw"),
                r_cutoff: number("--cutoff", 3.5)? as f32,
                bin_count: integer("--bins", 180, 1, 100000)?,
                smooth_passes: integer("--smooth", 0, 0, 1000)?,
                ..Default::default()
            };
            if let Some(centre) = args.get("--centre") {
                p.centre_mode = AdfCentreMode::ByElement;
                p.centre_symbol = centre.clone();
            }
            if args.contains_key("--neighbor-a") || args.contains_key("--neighbor-b") {
                if p.centre_symbol.is_empty() {
                    return Err("Pair selection requires --centre".into());
                }
                p.centre_mode = AdfCentreMode::ByPair;
                p.neigh_symbol1 = args.get("--neighbor-a").cloned().unwrap_or_default();
                p.neigh_symbol2 = args.get("--neighbor-b").cloned().unwrap_or_default();
            }
            let result = compute_adf(&structure, &p);
            if !result.valid {
                return Err(result.message);
            }
            text += "angle_deg,raw_count,value\n";
            for row in &result.bins {
                text += &format!("{},{},{}\n", f(row.angle_deg), f(row.count), f(row.value));
            }
        }
        "sro" => {
            let tolerance = number("--tolerance", 0.1)?;
            if tolerance <= 0.0 {
                return Err("Tolerance must be positive".into());
            }
            let result = analyze_sro(&structure, integer("--shells", 3, 1, 100)?, tolerance);
            if !result.success {
                return Err(result.message);
            }
            text += "method,shell,distance_A,centre,neighbor,alpha\n";
            for row in &result.warren_cowley.entries {
                text += &format!("Warren-Cowley,{},{},{},{},{}\n", row.shell, g17(row.shell_distance), row.center_element, row.neighbor_element, g17(row.alpha));
            }
            for row in &result.rao_curtin.entries {
                text += &format!("Rao-Curtin,{},{},{},{},{}\n", row.shell, g17(row.shell_distance), row.elem_i, row.elem_j, g17(row.alpha));
            }
        }
        "interstitial" => {
            let p = VoidParams {
                grid_resolution: integer("--resolution", 14, 2, 100)?,
                max_voids: integer("--limit", 0, 0, 100000)?,
                min_clearance: number("--clearance", 0.6)? as f32,
                min_separation: number("--separation", 0.5)? as f32,
                ..Default::default()
            };
            if p.min_clearance < 0.0 || p.min_separation <= 0.0 {
                return Err("Invalid clearance/separation".into());
            }
            let result = detect_interstitial_voids(&structure, &p);
            if !result.success {
                return Err(result.message);
            }
            text += "x_A,y_A,z_A,clearance_A,coordination,volume_A3,kind\n";
            for row in &result.regions {
                text += &format!(
                    "{},{},{},{},{},{},{}\n",
                    f(row.position.x),
                    f(row.position.y),
                    f(row.position.z),
                    f(row.clearance),
                    row.coordination,
                    f(row.volume),
                    row.kind as i32
                );
            }
        }
        "nye" => {
            let Some(reference) = args.get("--reference") else { return Err("nye requires --reference FILE (the undeformed structure)".into()) };
            let reference = load(reference)?;
            let cutoff = number("--cutoff", 3.0)? as f32;
            if cutoff.is_nan() || cutoff <= 0.0 {
                return Err("--cutoff must be positive".into());
            }
            let result = compute_nye_tensor(&reference, &structure, cutoff, periodic);
            if !result.success {
                return Err(result.message);
            }
            text += "index,symbol,alpha_xx,alpha_xy,alpha_xz,alpha_yx,alpha_yy,alpha_yz,alpha_zx,alpha_zy,alpha_zz,norm\n";
            for (i, alpha) in result.alpha.iter().enumerate() {
                text += &format!("{i},{}", structure.atoms.get(i).map(|a| a.symbol.as_str()).unwrap_or(""));
                for row in alpha {
                    for &v in row {
                        text += &format!(",{}", g17(v));
                    }
                }
                text += &format!(",{}\n", g17(result.alpha_norm[i]));
            }
        }
        "fit-dislocation" => {
            let Some(reference) = args.get("--reference") else { return Err("fit-dislocation requires --reference FILE (the undeformed structure)".into()) };
            if !args.contains_key("--line") {
                return Err("fit-dislocation requires --line \"u v w\" (Cartesian assumed line direction)".into());
            }
            let reference = load(reference)?;
            let line = parse_line("--line")?;
            let cutoff = number("--cutoff", 3.0)? as f32;
            if cutoff.is_nan() || cutoff <= 0.0 {
                return Err("--cutoff must be positive".into());
            }
            let area = number("--area", 1.0)?;
            let nye = compute_nye_tensor(&reference, &structure, cutoff, periodic);
            if !nye.success {
                return Err(nye.message);
            }
            let positions: Vec<DVec3> = reference.atoms.iter().map(|a| DVec3::new(a.x, a.y, a.z)).collect();
            let fit = fit_dislocation_from_nye(&positions, &nye, line, area);
            if !fit.success {
                return Err(fit.message);
            }
            text += "quantity,x,y,z\n";
            let (p, b) = (fit.line_position, fit.burgers_vector);
            text += &format!("line_position,{},{},{}\n", g17(p.x), g17(p.y), g17(p.z));
            text += &format!("burgers_vector,{},{},{}\n", g17(b.x), g17(b.y), g17(b.z));
        }
        "vitek" => {
            let Some(reference) = args.get("--reference") else { return Err("vitek requires --reference FILE (the undeformed structure)".into()) };
            if !args.contains_key("--line") {
                return Err("vitek requires --line \"u v w\" (Cartesian dislocation line/Burgers direction)".into());
            }
            let reference = load(reference)?;
            let line = parse_line("--line")?;
            let burgers = number("--burgers", 0.0)?;
            let cutoff = number("--cutoff", 0.0)? as f32;
            let result = compute_vitek_map(&reference, &structure, line, burgers, cutoff, periodic);
            if !result.success {
                return Err(result.message);
            }
            text += "index_i,index_j,xi_A,yi_A,zi_A,xj_A,yj_A,zj_A,screw_A,edge_1_A,edge_2_A\n";
            for p in &result.pairs {
                let (a, b) = (p.position_i, p.position_j);
                text += &format!(
                    "{},{},{},{},{},{},{},{},{},{},{}\n",
                    p.index_i,
                    p.index_j,
                    g17(a.x),
                    g17(a.y),
                    g17(a.z),
                    g17(b.x),
                    g17(b.y),
                    g17(b.z),
                    g17(p.screw_component),
                    g17(p.edge_component.x),
                    g17(p.edge_component.y)
                );
            }
        }
        "pattern-init" => {
            let cutoff = number("--cutoff", 0.0)? as f32;
            if cutoff.is_nan() || cutoff <= 0.0 {
                return Err("--cutoff must be positive".into());
            }
            let result = build_pattern(&structure, cutoff, periodic);
            if !result.success {
                return Err(result.message);
            }
            text += "dx,dy,dz\n";
            for d in &result.directions {
                text += &format!("{},{},{}\n", g17(d.x), g17(d.y), g17(d.z));
            }
        }
        "pattern-detect" => {
            let Some(pattern_path) = args.get("--pattern") else { return Err("pattern-detect requires --pattern FILE (written by pattern-init)".into()) };
            let cutoff = number("--cutoff", 0.0)? as f32;
            if cutoff.is_nan() || cutoff <= 0.0 {
                return Err("--cutoff must be positive".into());
            }
            let angle_threshold = number("--angle-threshold", 10.0)?;
            let content = std::fs::read(pattern_path).map_err(|_| format!("Cannot read pattern file '{pattern_path}'"))?;
            let mut pattern = CrystalPattern::default();
            for row in text_lines(&String::from_utf8_lossy(&content)).iter().skip(1) {
                if row.is_empty() {
                    continue;
                }
                let mut s = IStream::new(row);
                let parsed = (|| Some((s.double()?, s.char()?, s.double()?, s.char()?, s.double()?)))();
                let Some((dx, _, dy, _, dz)) = parsed else { return Err(format!("Malformed pattern file row: {row}")) };
                pattern.directions.push(DVec3::new(dx, dy, dz));
            }
            pattern.success = !pattern.directions.is_empty();
            if !pattern.success {
                return Err("Pattern file contains no directions".into());
            }
            let result = detect_pattern(&structure, &pattern, cutoff, angle_threshold, periodic);
            if !result.success {
                return Err(result.message);
            }
            text += "index,matched,neighbor_count,max_angle_deviation_deg\n";
            for row in &result.rows {
                text += &format!("{},{},{},{}\n", row.index, i32::from(row.matched), row.neighbor_count, g17(row.max_angle_deviation_deg));
            }
        }
        _ => return Err("Analysis must be cna, rdf, adf, sro, interstitial, nye, vitek, pattern-init or pattern-detect".into()),
    }
    std::fs::write(&output, text_mode(&text)).map_err(|_| "Cannot write analysis output".to_string())
}
