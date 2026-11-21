//! Crystal builders and editing tools (the C++ cli/BuildCrystalCLI): bulk,
//! custom mesh fill, solid solution, primitive cell, surface, SQS, vacancy and strain.

use std::collections::BTreeMap;

use af_algorithms::bulk_crystal::{
    apply_element_to_atom, apply_system_constraints, build_bulk_crystal, crystal_system_label, validate_parameters, CrystalSystem, LatticeParameters,
    SPACE_GROUP_RANGES,
};
use af_algorithms::csl::reduce_to_primitive;
use af_algorithms::mesh::{parse_obj_mesh, parse_stl_mesh};
use af_algorithms::nanocrystal::{build_nanocrystal, NanoGenerationMode, NanoParams, NanoShape};
use af_algorithms::solid_solution::{build_substitutional_solid_solution, ElementFraction, SssParams};
use af_algorithms::sqs_drag::{build_sqs, SqsParams};
use af_algorithms::strain_vacancy::{apply_strain, build_vacancies, engineering_strain, VacancyParams};
use af_algorithms::surface::{build_surface, SurfaceParams};
use af_core::elements::default_colours;
use af_core::format::{parse_prefix_f32, stod};
use af_core::structure::{AtomSite, Structure};
use glam::Vec3;

use crate::args::{detect_format, err, g, load, out, save, symbol_z, Args, IStream};

const PI: f64 = 3.14159265358979323846;

fn saved(structure: &Structure, path: &str) -> Result<i32, String> {
    if !save(structure, path, &detect_format(path)) {
        err(&format!("Error: failed to save structure to '{path}'\n"));
        return Ok(1);
    }
    out(&format!("Saved to: {path}\n"));
    Ok(0)
}

/// A fractional coordinate as a decimal or a fraction ("0.25", "1/3", "-2/3").
fn parse_fraction(text: &str) -> Option<f64> {
    match text.find('/') {
        None => match stod(text) {
            Ok((v, used)) if used == text.len() && v.is_finite() => Some(v),
            _ => None,
        },
        Some(slash) => {
            let (numerator, used) = stod(&text[..slash]).ok()?;
            if used != slash {
                return None;
            }
            let rest = &text[slash + 1..];
            let (denominator, used) = stod(rest).ok()?;
            if used != rest.len() || denominator == 0.0 {
                return None;
            }
            let v = numerator / denominator;
            v.is_finite().then_some(v)
        }
    }
}

pub fn run_bulk(args: &Args) -> Result<i32, String> {
    let system_text = args.find("--system");
    let mut system = CrystalSystem::Cubic;
    if let Some(s) = system_text {
        system = match s.to_ascii_lowercase().as_str() {
            "triclinic" => CrystalSystem::Triclinic,
            "monoclinic" => CrystalSystem::Monoclinic,
            "orthorhombic" => CrystalSystem::Orthorhombic,
            "tetragonal" => CrystalSystem::Tetragonal,
            "trigonal" => CrystalSystem::Trigonal,
            "hexagonal" => CrystalSystem::Hexagonal,
            "cubic" => CrystalSystem::Cubic,
            _ => {
                err(&format!("Error: unknown crystal system '{s}'\n"));
                return Ok(1);
            }
        };
    }
    let sg = args.int("--spacegroup", 225)?;
    if !(1..=230).contains(&sg) {
        err("Error: --spacegroup must be a number from 1 to 230\n");
        return Ok(1);
    }
    // The space group fixes the crystal system; --system may only confirm it.
    for range in &SPACE_GROUP_RANGES {
        if sg < range.first || sg > range.last {
            continue;
        }
        if let Some(given) = system_text {
            if range.system != system {
                let label = range.label.to_ascii_lowercase();
                err(&format!("Error: space group {sg} is {label}, not {given}; leave out --system or give --system {label}\n"));
                return Ok(1);
            }
        }
        system = range.system;
    }
    let mut lp = LatticeParameters::default();
    lp.a = args.double("--a", 4.0)?;
    lp.b = args.double("--b", lp.a)?;
    lp.c = args.double("--c", lp.a)?;
    lp.alpha = args.double("--alpha", 90.0)?;
    lp.beta = args.double("--beta", 90.0)?;
    lp.gamma = args.double("--gamma", 90.0)?;
    // A value the crystal system fixes is refused rather than overwritten.
    let given = lp;
    apply_system_constraints(system, &mut lp);
    let checks = [("--b", given.b, lp.b), ("--c", given.c, lp.c), ("--alpha", given.alpha, lp.alpha), ("--beta", given.beta, lp.beta), ("--gamma", given.gamma, lp.gamma)];
    for (flag, typed, fixed) in checks {
        if args.find(flag).is_none() || (typed - fixed).abs() <= 1e-6 * f64::max(1.0, fixed.abs()) {
            continue;
        }
        let label = crystal_system_label(system).to_ascii_lowercase();
        let mut message = format!("Error: a {label} cell has {} = {}, not {}; leave out {flag}.", &flag[2..], g(fixed), g(typed));
        if system == CrystalSystem::Trigonal && given.alpha == given.beta && given.beta == given.gamma && given.alpha != 90.0 {
            // Rhombohedral axes to the hexagonal ones the builder uses.
            let alpha = given.alpha * PI / 180.0;
            message += &format!(
                " Trigonal groups are built on hexagonal axes: rhombohedral a = {}, alpha = {} is --a {} --c {} (with the sites in hexagonal coordinates).",
                g(given.a),
                g(given.alpha),
                g(2.0 * given.a * (alpha / 2.0).sin()),
                g(given.a * (3.0 * (1.0 + 2.0 * alpha.cos())).sqrt())
            );
        }
        err(&format!("{message}\n"));
        return Ok(1);
    }
    if let Err(e) = validate_parameters(system, &lp) {
        err(&format!("Error: invalid lattice parameters \u{2013} {e}\n"));
        return Ok(1);
    }
    let colours = default_colours();
    let mut asymmetric = Vec::new();
    for text in args.find_all("--atom") {
        let mut s = IStream::new(&text);
        let parsed = (|| {
            let (sym, tx, ty, tz) = (s.word()?, s.word()?, s.word()?, s.word()?);
            if s.word().is_some() {
                return None;
            }
            Some((sym, parse_fraction(&tx)?, parse_fraction(&ty)?, parse_fraction(&tz)?))
        })();
        let Some((sym, fx, fy, fz)) = parsed else {
            err(&format!("Error: cannot parse --atom value '{text}'.  Expected: \"SYMBOL fx fy fz\" (decimals or fractions such as 1/3)\n"));
            return Ok(1);
        };
        let z = symbol_z(&sym);
        if z <= 0 {
            err(&format!("Error: unknown element symbol '{sym}'\n"));
            return Ok(1);
        }
        let mut site = AtomSite { x: fx, y: fy, z: fz, ..Default::default() };
        apply_element_to_atom(&mut site, z, &colours);
        asymmetric.push(site);
    }
    if asymmetric.is_empty() {
        err("Error: no --atom arguments specified.  At least one asymmetric unit atom is required.\n");
        return Ok(1);
    }
    let Some(output) = args.find("--output") else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let mut structure = Structure::default();
    let result = build_bulk_crystal(&mut structure, system, sg, &lp, &asymmetric, &colours, args.int("--origin", 1)?);
    if !result.success {
        err(&format!("Error building bulk crystal: {}\n", result.message));
        return Ok(1);
    }
    out(&format!(
        "Built bulk crystal: {} atoms, space group {} ({})\n",
        result.generated_atoms, result.space_group_number, result.space_group_symbol
    ));
    saved(&structure, output)
}

pub fn run_custom(args: &Args) -> Result<i32, String> {
    let (input, mesh_path, output) = (args.find("--input"), args.find("--mesh"), args.find("--output"));
    let Some(input) = input else {
        err("Error: --input <file> is required for --build custom\n");
        return Ok(1);
    };
    let Some(mesh_path) = mesh_path else {
        err("Error: --mesh <file> (OBJ or STL) is required for --build custom\n");
        return Ok(1);
    };
    let Some(output) = output else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let reference = match load(input) {
        Ok(s) => s,
        Err(e) => {
            err(&format!("Error loading reference structure: {e}\n"));
            return Ok(1);
        }
    };
    let ext = mesh_path.rfind('.').map(|d| mesh_path[d + 1..].to_ascii_lowercase()).unwrap_or_default();
    let mesh = match ext.as_str() {
        "obj" => parse_obj_mesh(std::path::Path::new(mesh_path)),
        "stl" => parse_stl_mesh(std::path::Path::new(mesh_path)),
        _ => {
            err(&format!("Error: unsupported mesh format '.{ext}'.  Use OBJ or STL.\n"));
            return Ok(1);
        }
    };
    let mesh = match mesh {
        Ok(m) => m,
        Err(e) => {
            err(&format!("Error loading mesh: {e}\n"));
            return Ok(1);
        }
    };
    if mesh.vertices.is_empty() || mesh.indices.is_empty() {
        err("Error: mesh file contains no geometry\n");
        return Ok(1);
    }
    let (mut mn, mut mx) = (mesh.vertices[0], mesh.vertices[0]);
    for &v in &mesh.vertices {
        mn = mn.min(v);
        mx = mx.max(v);
    }
    let half: Vec3 = (mx - mn) * 0.5;
    let mut params = NanoParams {
        generation_mode: NanoGenerationMode::Shape,
        shape: NanoShape::MeshModel,
        model_scale: args.double("--scale", 1.0)? as f32,
        model_hx: half.x,
        model_hy: half.y,
        model_hz: half.z,
        vacuum_padding: args.double("--vacuum", 5.0)? as f32,
        set_output_cell: !args.has("--no-cell"),
        auto_center_from_atoms: true,
        ..Default::default()
    };
    if let Some(center) = args.find("--center") {
        let mut s = IStream::new(center);
        let (Some(x), Some(y), Some(z)) = (s.float(), s.float(), s.float()) else {
            err("Error: --center needs three numbers \"x y z\"\n");
            return Ok(1);
        };
        params.auto_center_from_atoms = false;
        (params.cx, params.cy, params.cz) = (x, y, z);
    }
    replication(args, &mut params)?;
    if let Some(miller) = args.find("--miller") {
        params.apply_crystal_orientation = true;
        params.use_miller_orientation = true;
        let mut s = IStream::new(miller);
        // A failed read stores 0 (C++11 streams).
        params.miller_h = s.int().unwrap_or(0);
        params.miller_k = s.int().unwrap_or(0);
        params.miller_l = s.int().unwrap_or(0);
    } else {
        let (rx, ry, rz) = (args.double("--rotx", 0.0)? as f32, args.double("--roty", 0.0)? as f32, args.double("--rotz", 0.0)? as f32);
        if rx != 0.0 || ry != 0.0 || rz != 0.0 {
            params.apply_crystal_orientation = true;
            params.use_miller_orientation = false;
            (params.orient_x_deg, params.orient_y_deg, params.orient_z_deg) = (rx, ry, rz);
        }
    }
    let colours = default_colours();
    let mut structure = Structure::default();
    let result = build_nanocrystal(&mut structure, &reference, &params, &colours, &mesh.vertices, &mesh.indices);
    if !result.success {
        err(&format!("Error building custom structure: {}\n", result.message));
        return Ok(1);
    }
    out(&format!("Built custom mesh-fill: {} atoms\n", result.output_atoms));
    saved(&structure, output)
}

/// --repa/--repb/--repc: manual replication (5 for the unset ones), else automatic.
pub fn replication(args: &Args, params: &mut NanoParams) -> Result<(), String> {
    let (a, b, c) = (args.int("--repa", 0)?, args.int("--repb", 0)?, args.int("--repc", 0)?);
    if a > 0 || b > 0 || c > 0 {
        params.auto_replicate = false;
        params.rep_a = if a > 0 { a } else { 5 };
        params.rep_b = if b > 0 { b } else { 5 };
        params.rep_c = if c > 0 { c } else { 5 };
    } else {
        params.auto_replicate = true;
    }
    Ok(())
}

pub fn run_sss(args: &Args) -> Result<i32, String> {
    let Some(input) = args.find("--input") else {
        err("Error: --input <file> is required for --build sss\n");
        return Ok(1);
    };
    let Some(output) = args.find("--output") else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let base = match load(input) {
        Ok(s) => s,
        Err(e) => {
            err(&format!("Error loading input structure: {e}\n"));
            return Ok(1);
        }
    };
    let Some(fractions) = args.find("--frac") else {
        err("Error: --frac \"SYM=frac,...\" is required (e.g. --frac \"Cu=0.7,Zn=0.3\")\n");
        return Ok(1);
    };
    let mut params = SssParams { seed: args.int("--seed", 12345)? as u32, ..Default::default() };
    let mut tokens: Vec<&str> = fractions.split(',').collect();
    if tokens.last() == Some(&"") {
        tokens.pop();
    }
    for token in tokens {
        let token = token.trim_matches([' ', '\t']);
        if token.is_empty() {
            continue;
        }
        let Some(eq) = token.find('=') else {
            err(&format!("Error: cannot parse --frac token '{token}'.  Expected format: SYM=fraction\n"));
            return Ok(1);
        };
        let (sym, text) = (&token[..eq], &token[eq + 1..]);
        let fraction = match parse_prefix_f32(text) {
            Some((v, used)) if used == text.len() && v.is_finite() => v,
            _ => {
                err(&format!("Error: cannot parse fraction '{text}' in --frac\n"));
                return Ok(1);
            }
        };
        let z = symbol_z(sym);
        if z <= 0 {
            err(&format!("Error: unknown element symbol '{sym}' in --frac\n"));
            return Ok(1);
        }
        params.composition.push(ElementFraction { atomic_number: z, fraction });
    }
    if params.composition.is_empty() {
        err("Error: --frac string parsed no valid entries\n");
        return Ok(1);
    }
    // Fractions, not the desktop's percentages: they must sum to 1.
    let mut total = 0.0f32;
    for entry in &params.composition {
        total += entry.fraction;
    }
    if (total - 1.0).abs() > 1e-3 {
        err(&format!("Error: --frac fractions sum to {}; they must sum to 1 (e.g. \"Cu=0.7,Zn=0.3\")\n", g(f64::from(total))));
        return Ok(1);
    }
    let result = build_substitutional_solid_solution(&base, &params);
    if !result.success {
        err(&format!("Error building solid solution: {}\n", result.message));
        return Ok(1);
    }
    let mut text = format!("Built solid solution: {} atoms\n", result.output.atoms.len());
    for &(z, count) in &result.element_counts {
        text += &format!("  {}: {count} atoms\n", af_core::elements::symbol(z));
    }
    out(&text);
    saved(&result.output, output)
}

/// --input and --output of the tools that throw on errors.
fn input_output(args: &Args) -> Result<(Structure, String), String> {
    let (Some(input), Some(output)) = (args.find("--input"), args.find("--output")) else {
        return Err("--input and --output are required".into());
    };
    Ok((load(input)?, output.to_string()))
}

pub fn run_primitive(args: &Args) -> Result<i32, String> {
    let (mut structure, output) = input_output(args)?;
    let symprec = args.double("--symprec", 1e-3)?;
    let before = structure.atoms.len();
    if !reduce_to_primitive(&mut structure, symprec) {
        return Err("Could not reduce to a primitive cell (needs spglib support and a valid unit cell).".into());
    }
    if !save(&structure, &output, &detect_format(&output)) {
        return Err("Failed to save primitive structure".into());
    }
    out(&format!("Reduced from {before} to {} atoms.\n", structure.atoms.len()));
    Ok(0)
}

pub fn run_surface(args: &Args) -> Result<i32, String> {
    let (source, output) = input_output(args)?;
    let params = SurfaceParams {
        h: args.int("--h", 1)?,
        k: args.int("--k", 1)?,
        l: args.int("--l", 1)?,
        layers: args.int("--layers", 4)?,
        vacuum: args.double("--vacuum", 15.0)?,
        nmax: args.int("--nmax", 8)?,
        primitive_input: !args.has("--no-primitive"),
        primitive_symprec: args.double("--symprec", 1e-3)?,
        ..Default::default()
    };
    let result = build_surface(&source, &params);
    if !result.success {
        return Err(result.message);
    }
    if !save(&result.structure, &output, &detect_format(&output)) {
        return Err("Failed to save surface".into());
    }
    out(&format!("{}\n", result.message));
    Ok(0)
}

pub fn run_sqs(args: &Args) -> Result<i32, String> {
    let (source, output) = input_output(args)?;
    let mut composition = BTreeMap::new();
    for spec in args.find_all("--element") {
        let mut s = IStream::new(&spec);
        let (Some(symbol), Some(fraction)) = (s.word(), s.double()) else {
            return Err(format!("Cannot parse --element value '{spec}'. Expected: \"SYMBOL FRACTION\""));
        };
        composition.insert(symbol, fraction);
    }
    let params = SqsParams {
        composition,
        shells: args.int("--shells", 2)?,
        shell_tolerance: args.double("--shell-tolerance", 0.2)?,
        steps: args.int("--steps", 3000)?,
        start_temperature: args.double("--start-temp", 1.0)?,
        end_temperature: args.double("--end-temp", 0.02)?,
        seed: args.int("--seed", 1)? as u32,
    };
    let result = build_sqs(&source, &params);
    if !result.success {
        return Err(result.message);
    }
    if !save(&result.structure, &output, &detect_format(&output)) {
        return Err("Failed to save SQS structure".into());
    }
    out(&format!("{}\n", result.message));
    Ok(0)
}

pub fn run_vacancy(args: &Args) -> Result<i32, String> {
    let (source, output) = input_output(args)?;
    let params = VacancyParams {
        element: args.find("--element").unwrap_or("").to_string(),
        target_percentage: args.double("--percent", 0.0)?,
        target_count: args.int("--count", 0)?,
        min_separation: args.double("--min-separation", 0.0)?,
        seed: args.int("--seed", 1)? as u32,
    };
    let result = build_vacancies(&source, &params);
    if !result.success {
        return Err(result.message);
    }
    if !save(&result.structure, &output, &detect_format(&output)) {
        return Err("Failed to save vacancy structure".into());
    }
    out(&format!("{}\n", result.message));
    Ok(0)
}

pub fn run_strain(args: &Args) -> Result<i32, String> {
    let (source, output) = input_output(args)?;
    let mut params = engineering_strain(
        args.double("--exx", 0.0)?,
        args.double("--eyy", 0.0)?,
        args.double("--ezz", 0.0)?,
        args.double("--exy", 0.0)?,
        args.double("--exz", 0.0)?,
        args.double("--eyz", 0.0)?,
    );
    if let Some(matrix) = args.find("--matrix") {
        let mut s = IStream::new(matrix);
        let mut f = [0.0; 9];
        for value in f.iter_mut() {
            *value = s.double().ok_or("--matrix needs 9 numbers")?;
        }
        for i in 0..3 {
            for j in 0..3 {
                params.f[i][j] = f[i * 3 + j];
            }
        }
    }
    let result = apply_strain(&source, &params);
    if !result.success {
        return Err(result.message);
    }
    if !save(&result.structure, &output, &detect_format(&output)) {
        return Err("Failed to save strained structure".into());
    }
    out(&format!("{}\n", result.message));
    Ok(0)
}
