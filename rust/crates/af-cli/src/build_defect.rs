//! Defect and interface builders (the C++ cli/BuildDefectCLI): dislocations,
//! CSL grain boundaries, interfaces and stacking-fault sequences.

use af_algorithms::csl::{build_csl_bicrystal, classify_boundary_type, compute_gb_info};
use af_algorithms::dislocation::{build_dislocation, DislocationCharacter, DislocationParams, DislocationShape, ElasticSymmetry};
use af_algorithms::interface::{build_interface, find_matches, prepare_slab, repeat_interface_xy, suggested_match, Assembly, SearchParams, Side, StrainOn};
use af_algorithms::stacking_fault::{build_stacking_fault_sequence, StackingFaultCellMode, StackingFaultParams, StackingFaultPlane};
use glam::{IVec3, Vec2, Vec3};

use crate::args::{detect_format, err, g, load, out, parse_ivec3, parse_vec3, save, Args, IStream};

fn saved(structure: &af_core::structure::Structure, path: &str) -> Result<i32, String> {
    if !save(structure, path, &detect_format(path)) {
        err(&format!("Error: failed to save structure to '{path}'\n"));
        return Ok(1);
    }
    out(&format!("Saved to: {path}\n"));
    Ok(0)
}

/// Pieces of std::getline(..., ';') (no empty piece after a final separator).
fn semicolon_pieces(text: &str) -> Vec<&str> {
    let mut pieces: Vec<&str> = text.split(';').collect();
    if pieces.last() == Some(&"") {
        pieces.pop();
    }
    pieces
}

pub fn run_dislocation(args: &Args) -> Result<i32, String> {
    let (input, output) = (args.find("--input"), args.find("--output"));
    let Some(input) = input else {
        err("Error: --input <file> is required for --build dislocation\n");
        return Ok(1);
    };
    let Some(output) = output else {
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
    let mut p = DislocationParams { auto_directions: !args.has("--manual-vectors"), ..Default::default() };
    if let Some(character) = args.find("--character") {
        p.character = match character.to_ascii_lowercase().as_str() {
            "edge" => DislocationCharacter::Edge,
            "screw" => DislocationCharacter::Screw,
            "mixed" => DislocationCharacter::Mixed,
            _ => {
                err(&format!("Error: unsupported --character value '{character}'\n"));
                return Ok(1);
            }
        };
    }
    if let Some(shape) = args.find("--shape") {
        p.shape = match shape.to_ascii_lowercase().as_str() {
            "halfplane" => DislocationShape::HalfPlane,
            "cylinder" => DislocationShape::Cylinder,
            "sphere" => DislocationShape::Sphere,
            "ellipsoid" => DislocationShape::Ellipsoid,
            "freeform" => DislocationShape::Freeform2D,
            _ => {
                err(&format!("Error: unsupported --shape value '{shape}'\n"));
                return Ok(1);
            }
        };
    }
    for (flag, what, target) in [("--plane", "h k l", &mut p.plane_hkl), ("--burgers", "u v w", &mut p.burgers_uvw), ("--line", "u v w", &mut p.line_uvw)] {
        if let Some(text) = args.find(flag) {
            match parse_ivec3(text) {
                Some(v) => *target = IVec3::from_array(v),
                None => {
                    err(&format!("Error: could not parse {flag} \"{what}\"\n"));
                    return Ok(1);
                }
            }
        }
    }
    if let Some(text) = args.find("--line-frac") {
        p.use_fractional_line_point = true;
        match parse_vec3(text) {
            Some(v) => p.line_point_fractional = Vec3::from_array(v),
            None => {
                err("Error: could not parse --line-frac \"fx fy fz\"\n");
                return Ok(1);
            }
        }
    }
    if args.has("--line-center") {
        p.use_fractional_line_point = false;
    }
    if let Some(text) = args.find("--line-offset") {
        match parse_vec3(text) {
            Some(v) => p.line_point_cartesian_offset = Vec3::from_array(v),
            None => {
                err("Error: could not parse --line-offset \"x y z\"\n");
                return Ok(1);
            }
        }
    }
    p.burgers_scale = args.double("--bscale", f64::from(p.burgers_scale))? as f32;
    p.burgers_override_magnitude = args.double("--bmag", f64::from(p.burgers_override_magnitude))? as f32;
    p.mixed_character_angle_deg = args.double("--mixed-angle", f64::from(p.mixed_character_angle_deg))? as f32;
    p.poisson_ratio = args.double("--nu", f64::from(p.poisson_ratio))? as f32;
    p.core_radius = args.double("--core", f64::from(p.core_radius))? as f32;
    p.cutoff_radius = args.double("--cutoff", f64::from(p.cutoff_radius))? as f32;
    p.line_half_length = args.double("--line-half", f64::from(p.line_half_length))? as f32;
    p.cylinder_radius = args.double("--cyl-radius", f64::from(p.cylinder_radius))? as f32;
    p.sphere_radius = args.double("--sphere-radius", f64::from(p.sphere_radius))? as f32;
    if let Some(text) = args.find("--ellipsoid") {
        match parse_vec3(text) {
            Some(v) => p.ellipsoid_radii = Vec3::from_array(v),
            None => {
                err("Error: could not parse --ellipsoid \"rx ry rz\"\n");
                return Ok(1);
            }
        }
    }
    if let Some(poly) = args.find("--poly2d") {
        p.freeform_points.clear();
        for token in semicolon_pieces(poly) {
            let mut s = IStream::new(token);
            let (Some(x), Some(y)) = (s.float(), s.float()) else {
                err(&format!("Error: could not parse polygon point '{token}' in --poly2d\n"));
                return Ok(1);
            };
            p.freeform_points.push(Vec2::new(x, y));
        }
        if p.freeform_points.len() < 3 {
            err("Error: --poly2d needs at least 3 points\n");
            return Ok(1);
        }
    }
    if args.has("--dipole") {
        p.dipole = true;
        if let Some(offset) = args.find("--dipole-offset") {
            let mut s = IStream::new(offset);
            let (Some(dx), Some(dy)) = (s.float(), s.float()) else {
                err("Error: could not parse --dipole-offset \"dx dy\"\n");
                return Ok(1);
            };
            p.dipole_offset = Vec2::new(dx, dy);
        }
    }
    if args.has("--anisotropic") {
        p.anisotropic_elasticity = true;
        if let Some(symmetry) = args.find("--elastic-symmetry") {
            p.elastic_symmetry = match symmetry.to_ascii_lowercase().as_str() {
                "cubic" => ElasticSymmetry::Cubic,
                "hexagonal" => ElasticSymmetry::Hexagonal,
                _ => {
                    err("Error: --elastic-symmetry must be cubic or hexagonal\n");
                    return Ok(1);
                }
            };
        }
        p.elastic_c11 = args.double("--elastic-c11", 0.0)?;
        p.elastic_c12 = args.double("--elastic-c12", 0.0)?;
        p.elastic_c44 = args.double("--elastic-c44", 0.0)?;
        p.elastic_c13 = args.double("--elastic-c13", 0.0)?;
        p.elastic_c33 = args.double("--elastic-c33", 0.0)?;
        p.elastic_noise_amplitude = args.double("--elastic-noise", p.elastic_noise_amplitude)?;
        if p.elastic_c11 <= 0.0 || p.elastic_c12 < 0.0 || p.elastic_c44 <= 0.0 {
            err("Error: --anisotropic requires positive --elastic-c11/--elastic-c44 (and --elastic-c12 >= 0), in GPa\n");
            return Ok(1);
        }
        if p.elastic_symmetry == ElasticSymmetry::Hexagonal && (p.elastic_c13 <= 0.0 || p.elastic_c33 <= 0.0) {
            err("Error: --elastic-symmetry hexagonal also requires positive --elastic-c13/--elastic-c33, in GPa\n");
            return Ok(1);
        }
    }
    let result = build_dislocation(&base, &p);
    if !result.success {
        err(&format!("Error building dislocation: {}\n", result.message));
        return Ok(1);
    }
    out(&format!("{}\nShifted atoms: {}\nValidation: {}\n", result.message, result.shifted_atom_count, result.validation.message));
    saved(&result.output, output)
}

pub fn run_gb(args: &Args) -> Result<i32, String> {
    let Some(input) = args.find("--input") else {
        err("Error: --input <file> is required for --build gb\n");
        return Ok(1);
    };
    let structure = match load(input) {
        Ok(s) => s,
        Err(e) => {
            err(&format!("Error loading input structure: {e}\n"));
            return Ok(1);
        }
    };
    if !structure.has_unit_cell {
        err("Error: input structure has no unit cell\n");
        return Ok(1);
    }
    let mut axis = [0, 0, 1];
    if let Some(text) = args.find("--axis") {
        match parse_ivec3(text) {
            Some(v) => axis = v,
            None => {
                err(&format!("Error: cannot parse --axis value '{text}'.  Expected: \"u v w\" (three integers)\n"));
                return Ok(1);
            }
        }
    }
    let sigma_target = args.int("--sigma", 0)?;
    // An explicit sigma names a coincidence site lattice: an odd number >= 3.
    if args.find("--sigma").is_some() && (sigma_target < 3 || sigma_target % 2 == 0) {
        err("Error: --sigma must be an odd number of at least 3 (omit it to search up to --sigmamax)\n");
        return Ok(1);
    }
    let sigma_max = args.int("--sigmamax", 100)?;
    let plane_index = args.int("--plane", 0)?;
    let (uc_a, uc_b) = (args.int("--uca", 1)?, args.int("--ucb", 1)?);
    let vacuum = args.double("--vacuum", 0.0)? as f32;
    let gap = args.double("--gap", 0.0)? as f32;
    let overlap = args.double("--overlap", 0.0)? as f32;
    let conventional = args.has("--conventional");
    let Some(output) = args.find("--output") else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let candidates = compute_gb_info(&axis, sigma_max);
    if candidates.is_empty() {
        err(&format!("Error: no \u{3a3} candidates found for axis [{} {} {}] with max \u{3a3}={sigma_max}\n", axis[0], axis[1], axis[2]));
        return Ok(1);
    }
    let mut selected = 0;
    if sigma_target > 0 {
        match candidates.iter().position(|c| c.sigma == sigma_target) {
            Some(i) => selected = i,
            None => {
                let list: String = candidates.iter().map(|c| format!(" {}", c.sigma)).collect();
                err(&format!(
                    "Error: \u{3a3}={sigma_target} not found for axis [{} {} {}]\nAvailable \u{3a3} values:{list}\n",
                    axis[0], axis[1], axis[2]
                ));
                return Ok(1);
            }
        }
    }
    if !(0..=2).contains(&plane_index) {
        err("Error: --plane must be 0, 1, or 2\n");
        return Ok(1);
    }
    let sel = &candidates[selected];
    let plane = sel.plane[plane_index as usize];
    out(&format!(
        "Using \u{3a3}={} (\u{3b8}={}\u{b0}), axis [{} {} {}], plane ({} {} {})\n",
        sel.sigma,
        g(f64::from(sel.theta_deg)),
        axis[0],
        axis[1],
        axis[2],
        plane[0],
        plane[1],
        plane[2]
    ));
    let bicrystal = build_csl_bicrystal(&structure, &axis, sel.sigma, sigma_max, plane_index as usize, uc_a, uc_b, vacuum, gap, overlap, conventional)?;
    out(&format!(
        "Built GB bicrystal: {} atoms (removed {} overlapping)\nBoundary type: {}\n",
        bicrystal.structure.atoms.len(),
        bicrystal.removed,
        classify_boundary_type(&axis, &plane)
    ));
    saved(&bicrystal.structure, output)
}

pub fn run_interface(args: &Args) -> Result<i32, String> {
    let (Some(path_a), Some(path_b)) = (args.find("--layerA"), args.find("--layerB")) else {
        err("Error: --layerA <file> and --layerB <file> are required\n");
        return Ok(1);
    };
    let Some(output) = args.find("--output") else {
        err("Error: --output <file> is required\n");
        return Ok(1);
    };
    let mut side_a = Side::default();
    let mut side_b = Side::default();
    match load(path_a) {
        Ok(s) => side_a.crystal = s,
        Err(e) => {
            err(&format!("Error loading layer A: {e}\n"));
            return Ok(1);
        }
    }
    match load(path_b) {
        Ok(s) => side_b.crystal = s,
        Err(e) => {
            err(&format!("Error loading layer B: {e}\n"));
            return Ok(1);
        }
    }
    for (flag, side) in [("--hklA", &mut side_a), ("--hklB", &mut side_b)] {
        let text = args.find(flag);
        side.cleave = text.is_some();
        if let Some(text) = text {
            let hkl = parse_vec3(text).filter(|v| v.iter().all(|x| *x == x.round())).ok_or(format!("{flag} needs three whole numbers, e.g. \"1 1 1\""))?;
            (side.h, side.k, side.l) = (hkl[0] as i32, hkl[1] as i32, hkl[2] as i32);
        }
    }
    let (layers_a, layers_b) = (args.int("--layersA", 4)?, args.int("--layersB", 4)?);
    side_a.layers = layers_a;
    side_b.layers = layers_b;
    let max_strain = args.double("--max-strain", 0.03)?;
    let max_area = args.double("--max-area", 150.0)?;
    let nmax = args.int("--nmax", 6)?;
    let search = SearchParams {
        max_strain,
        max_area,
        max_repeats: args.int("--repeats", nmax)?,
        min_cell_angle: args.double("--min-angle", 20.0)?,
        angle_min: args.double("--angle-min", 0.0)?,
        angle_max: args.double("--angle-max", 180.0)?,
        angle_step: args.double("--angle-step", 1.0)?,
        ..Default::default()
    };
    let mut assembly = Assembly { gap: args.double("--gap", 2.0)?, vacuum: args.double("--vacuum", 0.0)?, ..Default::default() };
    assembly.strain_on = match args.find("--strain-on").unwrap_or("B") {
        "B" | "b" => StrainOn::B,
        "A" | "a" => StrainOn::A,
        "both" => StrainOn::Both,
        _ => return Err("--strain-on must be B, A or both".into()),
    };
    if let Some(text) = args.find("--shift") {
        let mut s = IStream::new(text);
        let (Some(x), Some(y)) = (s.double(), s.double()) else {
            return Err("--shift needs two fractions, e.g. \"0.5 0\"".into());
        };
        assembly.shift = [x, y];
    }
    let (repx, repy) = (args.int("--repx", 1)?, args.int("--repy", 1)?);
    if layers_a < 1
        || layers_b < 1
        || layers_a > 1000
        || layers_b > 1000
        || search.max_strain <= 0.0
        || search.max_area <= 0.0
        || search.max_repeats < 1
        || search.max_repeats > 30
        || search.angle_step <= 0.0
        || search.angle_max < search.angle_min
        || assembly.gap < 0.0
        || assembly.vacuum < 0.0
        || !(1..=100).contains(&repx)
        || !(1..=100).contains(&repy)
    {
        return Err("Invalid interface layers, search limits, gap, vacuum or repeats".into());
    }
    let slab_a = match prepare_slab(&side_a) {
        Ok(s) => s,
        Err(e) => {
            err(&format!("Error: layer A {e}\n"));
            return Ok(1);
        }
    };
    let slab_b = match prepare_slab(&side_b) {
        Ok(s) => s,
        Err(e) => {
            err(&format!("Error: layer B {e}\n"));
            return Ok(1);
        }
    };
    let matches = find_matches(&slab_a, &slab_b, &search);
    if matches.is_empty() {
        err(&format!(
            "Error: no matches within {}% strain and {} A^2; raise --max-strain, --max-area or --repeats\n",
            g(search.max_strain * 100.0),
            g(search.max_area)
        ));
        return Ok(1);
    }
    let pick = args.int("--pick", suggested_match(&matches))?;
    if pick < 0 || pick as usize >= matches.len() {
        err(&format!("Error: --pick {pick} is out of range ({} matches)\n", matches.len()));
        return Ok(1);
    }
    let list = args.has("--list");
    let mut text = format!(
        "{} matches (fewest atoms first; * Pareto-optimal, > built):\n   #   atoms  strain%   area A2  angle  theta   cells A/B\n",
        matches.len()
    );
    for (i, m) in matches.iter().enumerate() {
        let i = i as i32;
        if !list && i >= 12 && i != pick {
            continue;
        }
        text += &format!(
            "{}{}{:>3} {:>7} {:>8.3} {:>9.2} {:>6.1} {:>6.1}   {}/{}\n",
            if i == pick { '>' } else { ' ' },
            if m.pareto { '*' } else { ' ' },
            i,
            m.atoms,
            m.mean_strain * 100.0,
            m.area,
            m.cell_angle,
            m.theta_deg,
            m.cells_a,
            m.cells_b
        );
    }
    if !list && matches.len() > 12 {
        text += "  ... (--list prints all)\n";
    }
    out(&text);
    let mut result = build_interface(&slab_a, &slab_b, &matches[pick as usize], &assembly);
    if repx > 1 || repy > 1 {
        result = repeat_interface_xy(&result, repx, repy);
    }
    out(&format!("Built interface: {} atoms\n", result.atoms.len()));
    saved(&result, output)
}

pub fn run_stacking_fault(args: &Args) -> Result<i32, String> {
    let (Some(input), Some(output)) = (args.find("--input"), args.find("--output")) else {
        return Err("--input and --output are required".into());
    };
    let source = load(input)?;
    let plane = args.int("--plane", 0)?;
    let plane = match plane {
        0 => StackingFaultPlane::Auto,
        1 => StackingFaultPlane::Fcc111,
        2 => StackingFaultPlane::HcpBasal,
        3 => StackingFaultPlane::HcpPrismatic,
        4 => StackingFaultPlane::HcpPyramidal,
        5 => StackingFaultPlane::Bcc110,
        6 => StackingFaultPlane::Bcc112,
        _ => return Err("Plane must be 0..6".into()),
    };
    let mut params = StackingFaultParams {
        plane,
        layer_count: args.int("--layers", 9)?,
        interval: args.double("--interval", 0.1)? as f32,
        max_displacement_factor: args.double("--maximum", 2.0)? as f32,
        ..Default::default()
    };
    if args.has("--orthogonal") {
        params.cell_mode = StackingFaultCellMode::OrthogonalCell;
    }
    let result = build_stacking_fault_sequence(&source, &params);
    if !result.success {
        return Err(result.message);
    }
    let frame = args.int("--frame", 0)?;
    if frame < 0 || frame as usize >= result.sequence.len() {
        return Err("Frame out of range".into());
    }
    if !save(&result.sequence[frame as usize].structure, output, &detect_format(output)) {
        return Err("Failed to save stacking fault".into());
    }
    if let Some(directory) = args.find("--sequence") {
        std::fs::create_dir_all(directory).map_err(|e| e.to_string())?;
        for (i, item) in result.sequence.iter().enumerate() {
            let path = std::path::Path::new(directory).join(format!("frame-{i}.vasp"));
            if !save(&item.structure, &path.to_string_lossy(), "vasp") {
                return Err("Failed to save sequence frame".into());
            }
        }
    }
    out(&format!("{}\nFrames: {}\n", result.message, result.sequence.len()));
    Ok(0)
}
