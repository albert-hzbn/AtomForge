//! Polycrystals, nanocrystals, amorphous packing, nanowires and core-shell
//! particles (the C++ cli/BuildNanoCLI).

use af_algorithms::amorphous::{build_amorphous_structure, AmorphousBoxMode, AmorphousElementSpec, AmorphousPairDist, AmorphousParams};
use af_algorithms::nano_pattern::{apply_core_shell, build_nanowire, CoreShellParams, NanowireParams};
use af_algorithms::nanocrystal::{build_nanocrystal, shape_label, NanoGenerationMode, NanoParams, NanoShape, WulffPlaneInput};
use af_algorithms::polycrystal::{build_polycrystal, GrainOrientation, GrainOrientationMode, PolyParams};
use af_core::elements::{default_colours, literature_covalent_radii};
use af_core::structure::Structure;

use crate::args::{detect_format, err, g, load, out, parse_vec3, save, symbol_z, Args, IStream};
use crate::build_crystal::replication;

fn saved(structure: &Structure, path: &str) -> Result<i32, String> {
    if !save(structure, path, &detect_format(path)) {
        err(&format!("Error: failed to save structure to '{path}'\n"));
        return Ok(1);
    }
    out(&format!("Saved to: {path}\n"));
    Ok(0)
}

fn reference(args: &Args, mode: &str) -> Result<Option<Structure>, String> {
    let Some(input) = args.find("--input") else {
        err(&format!("Error: --input <file> is required for --build {mode}\n"));
        return Ok(None);
    };
    let structure = match load(input) {
        Ok(s) => s,
        Err(e) => {
            err(&format!("Error loading input structure: {e}\n"));
            return Ok(None);
        }
    };
    if !structure.has_unit_cell {
        err("Error: input structure has no unit cell\n");
        return Ok(None);
    }
    Ok(Some(structure))
}

pub fn run_poly(args: &Args) -> Result<i32, String> {
    let Some(reference) = reference(args, "poly")? else { return Ok(1) };
    let mut params = PolyParams {
        size_x: args.double("--sizex", 50.0)? as f32,
        size_y: args.double("--sizey", 50.0)? as f32,
        size_z: args.double("--sizez", 50.0)? as f32,
        num_grains: args.int("--grains", 8)?,
        seed: args.int("--seed", 42)?,
        ..Default::default()
    };
    let euler = args.find_all("--euler");
    if !euler.is_empty() {
        for (i, text) in euler.iter().enumerate() {
            let mut s = IStream::new(text);
            let (Some(phi1), Some(phi), Some(phi2)) = (s.float(), s.float(), s.float()) else {
                err(&format!("Error: cannot parse --euler value '{text}'.  Expected: \"phi1 Phi phi2\"\n"));
                return Ok(1);
            };
            params.specified_orientations.push(GrainOrientation { grain_index: i as i32, phi1, phi, phi2 });
        }
        params.orientation_mode =
            if euler.len() as i32 >= params.num_grains { GrainOrientationMode::AllSpecified } else { GrainOrientationMode::PartialSpecified };
    }
    let grain_euler = args.find_all("--grain-euler");
    for text in &grain_euler {
        let mut s = IStream::new(text);
        let parsed = (|| Some((s.int()?, s.float()?, s.float()?, s.float()?)))();
        let Some((grain, phi1, phi, phi2)) = parsed.filter(|p| p.0 >= 1) else {
            err(&format!("Error: cannot parse --grain-euler value '{text}'.  Expected: \"N phi1 Phi phi2\" with N >= 1\n"));
            return Ok(1);
        };
        params.specified_orientations.push(GrainOrientation { grain_index: grain - 1, phi1, phi, phi2 });
    }
    if !grain_euler.is_empty() {
        params.orientation_mode = GrainOrientationMode::PartialSpecified;
    }
    let Some(output) = args.find("--output") else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let colours = default_colours();
    let mut structure = Structure::default();
    let result = build_polycrystal(&mut structure, &reference, &params, &colours);
    if !result.success {
        err(&format!("Error building polycrystal: {}\n", result.message));
        return Ok(1);
    }
    out(&format!("Built polycrystal: {} atoms, {} grains\n", result.output_atoms, result.num_grains));
    saved(&structure, output)
}

pub fn run_nano(args: &Args) -> Result<i32, String> {
    let Some(reference) = reference(args, "nano")? else { return Ok(1) };
    let Some(output) = args.find("--output") else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let mut params = NanoParams { generation_mode: NanoGenerationMode::Shape, ..Default::default() };
    if let Some(shape) = args.find("--shape") {
        match shape.to_ascii_lowercase().as_str() {
            "sphere" => params.shape = NanoShape::Sphere,
            "ellipsoid" => params.shape = NanoShape::Ellipsoid,
            "box" => params.shape = NanoShape::Box,
            "cylinder" => params.shape = NanoShape::Cylinder,
            "octahedron" => params.shape = NanoShape::Octahedron,
            "truncated-octahedron" => params.shape = NanoShape::TruncatedOctahedron,
            "cuboctahedron" => params.shape = NanoShape::Cuboctahedron,
            "wulff" => params.generation_mode = NanoGenerationMode::WulffConstruction,
            _ => {
                err(&format!("Error: unknown shape '{shape}'\n"));
                return Ok(1);
            }
        }
    }
    let facets = args.find_all("--facet");
    if params.generation_mode == NanoGenerationMode::WulffConstruction {
        if facets.is_empty() {
            return Err("Wulff construction requires at least one --facet".into());
        }
        for text in &facets {
            let mut s = IStream::new(text);
            let parsed = (|| Some((s.int()?, s.int()?, s.int()?, s.float()?)))();
            let extra = s.word().is_some();
            match parsed {
                Some((h, k, l, energy)) if !extra && energy.is_finite() && energy > 0.0 && !(h == 0 && k == 0 && l == 0) => {
                    params.wulff_planes.push(WulffPlaneInput { h, k, l, surface_energy: energy });
                }
                _ => return Err("Invalid --facet: expected nonzero h k l and positive finite energy".into()),
            }
        }
        params.wulff_max_radius = args.double("--radius", 20.0)? as f32;
    } else if !facets.is_empty() {
        return Err("--facet requires --shape wulff".into());
    }
    // Sizes are lengths: a zero or negative one is a mistake, not a shape.
    for size in ["--radius", "--rx", "--ry", "--rz", "--hx", "--hy", "--hz", "--trunc", "--cylradius", "--cylheight"] {
        if args.find(size).is_some() && !(args.double(size, 1.0)? > 0.0) {
            return Err(format!("{size} must be a positive length in Angstrom"));
        }
    }
    let radius = args.double("--radius", 15.0)? as f32;
    params.sphere_radius = radius;
    params.oct_radius = radius;
    params.trunc_oct_radius = args.double("--radius", 18.0)? as f32;
    params.cubo_radius = radius;
    params.ellip_rx = args.double("--rx", 15.0)? as f32;
    params.ellip_ry = args.double("--ry", 12.0)? as f32;
    params.ellip_rz = args.double("--rz", 10.0)? as f32;
    params.box_hx = args.double("--hx", 15.0)? as f32;
    params.box_hy = args.double("--hy", 15.0)? as f32;
    params.box_hz = args.double("--hz", 15.0)? as f32;
    params.trunc_oct_trunc = args.double("--trunc", 12.0)? as f32;
    params.cyl_radius = args.double("--cylradius", 12.0)? as f32;
    params.cyl_height = args.double("--cylheight", 30.0)? as f32;
    params.cyl_axis = args.int("--cylaxis", 2)?;
    params.vacuum_padding = args.double("--vacuum", 5.0)? as f32;
    params.set_output_cell = !args.has("--no-cell");
    params.auto_center_from_atoms = true;
    if let Some(center) = args.find("--center") {
        let Some(c) = parse_vec3(center) else {
            err(&format!("Error: cannot parse --center value '{center}'.  Expected: \"x y z\"\n"));
            return Ok(1);
        };
        params.auto_center_from_atoms = false;
        (params.cx, params.cy, params.cz) = (c[0], c[1], c[2]);
    }
    replication(args, &mut params)?;
    let colours = default_colours();
    let mut structure = Structure::default();
    let result = build_nanocrystal(&mut structure, &reference, &params, &colours, &[], &[]);
    if !result.success {
        err(&format!("Error building nanocrystal: {}\n", result.message));
        return Ok(1);
    }
    let shape = if params.generation_mode == NanoGenerationMode::WulffConstruction { "Wulff" } else { shape_label(params.shape) };
    out(&format!("Built nanocrystal: {} atoms (shape: {shape})\n", result.output_atoms));
    saved(&structure, output)
}

pub fn run_amorphous(args: &Args) -> Result<i32, String> {
    let Some(output) = args.find("--output") else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let specs = args.find_all("--element");
    if specs.is_empty() {
        err("Error: at least one --element \"SYMBOL N\" is required\n");
        return Ok(1);
    }
    let mut params = AmorphousParams::default();
    for text in &specs {
        let mut s = IStream::new(text);
        let (Some(sym), Some(count)) = (s.word(), s.int()) else {
            err(&format!("Error: cannot parse --element value '{text}'.  Expected: \"SYMBOL N\"\n"));
            return Ok(1);
        };
        let z = symbol_z(&sym);
        if z <= 0 {
            err(&format!("Error: unknown element symbol '{sym}'\n"));
            return Ok(1);
        }
        params.elements.push(AmorphousElementSpec { atomic_number: z, count });
    }
    // Box mode: manual if any of boxa/boxb/boxc are given, else auto density.
    if args.find("--boxa").is_some() || args.find("--boxb").is_some() || args.find("--boxc").is_some() {
        params.box_mode = AmorphousBoxMode::Manual;
        params.box_a = args.double("--boxa", 20.0)? as f32;
        params.box_b = args.double("--boxb", 20.0)? as f32;
        params.box_c = args.double("--boxc", 20.0)? as f32;
    } else {
        params.box_mode = AmorphousBoxMode::AutoFromDensity;
        params.target_density = args.double("--density", 2.0)? as f32;
    }
    params.cell_scale_factor = args.double("--scale", 1.0)? as f32;
    params.covalent_tolerance = args.double("--covtol", 0.75)? as f32;
    params.seed = args.int("--seed", 42)? as u32;
    params.max_attempts = args.int("--attempts", 1000)?;
    for text in args.find_all("--mindist") {
        let mut s = IStream::new(&text);
        let (Some(s1), Some(s2), Some(dist)) = (s.word(), s.word(), s.float()) else {
            err(&format!("Error: cannot parse --mindist value '{text}'.  Expected: \"SYMBOL1 SYMBOL2 dist\"\n"));
            return Ok(1);
        };
        let (z1, z2) = (symbol_z(&s1), symbol_z(&s2));
        if z1 <= 0 || z2 <= 0 {
            err(&format!("Error: unknown element symbol in --mindist '{text}'\n"));
            return Ok(1);
        }
        params.pair_distances.push(AmorphousPairDist { z1: z1.min(z2), z2: z1.max(z2), min_dist: dist });
    }
    params.periodic = !args.has("--no-periodic");
    let result = build_amorphous_structure(&params, &literature_covalent_radii(), &default_colours());
    if !result.success {
        err(&format!("Error building amorphous structure: {}\n", result.message));
        return Ok(1);
    }
    out(&format!("Built amorphous structure: {} atoms (density: {} g/cm3)\n", result.placed_atoms, g(f64::from(result.actual_density))));
    saved(&result.output, output)
}

pub fn run_nanowire(args: &Args) -> Result<i32, String> {
    let (Some(input), Some(output)) = (args.find("--input"), args.find("--output")) else {
        return Err("--input and --output are required".into());
    };
    let source = load(input)?;
    let params = NanowireParams {
        axis: args.int("--axis", 2)?,
        radius: args.double("--radius", 10.0)?,
        sides: args.int("--sides", 0)?,
        vacuum: args.double("--vacuum", 10.0)?,
        axis_repeats: args.int("--axis-repeats", 1)?,
    };
    let result = build_nanowire(&source, &params);
    if !result.success {
        return Err(result.message);
    }
    if !save(&result.structure, output, &detect_format(output)) {
        return Err("Failed to save nanowire".into());
    }
    out(&format!("{}\n", result.message));
    Ok(0)
}

pub fn run_core_shell(args: &Args) -> Result<i32, String> {
    let (Some(input), Some(output)) = (args.find("--input"), args.find("--output")) else {
        return Err("--input and --output are required".into());
    };
    let source = load(input)?;
    let mut params = CoreShellParams { core_radius: args.double("--core-radius", 5.0)?, ..Default::default() };
    if let Some(e) = args.find("--core-element") {
        params.core_element = e.to_string();
    }
    if let Some(e) = args.find("--shell-element") {
        params.shell_element = e.to_string();
    }
    if let Some(c) = args.find("--center") {
        let mut s = IStream::new(c);
        let (Some(x), Some(y), Some(z)) = (s.double(), s.double(), s.double()) else {
            return Err("Cannot parse --center \"x y z\"".into());
        };
        params.center = [x, y, z];
        params.use_centroid = false;
    }
    let result = apply_core_shell(&source, &params);
    if !result.success {
        return Err(result.message);
    }
    if !save(&result.structure, output, &detect_format(output)) {
        return Err("Failed to save core-shell structure".into());
    }
    out(&format!("{}\n", result.message));
    Ok(0)
}
