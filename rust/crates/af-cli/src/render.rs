//! `--render` (the C++ cli/RenderCLI): headless snapshots of a structure as
//! PNG (optionally a turntable series), or an animated GIF of a trajectory or
//! a turntable, with the display settings' look.

use std::path::{Path, PathBuf};

use af_core::display::{default_look, load_defaults, saved_defaults_path, Look as DisplayLook};
use af_core::occlusion::darken_occluded;
use af_core::structure::Structure;
use af_render::scene::{apply_default_view, build_frame_view, build_instance_data, BondStyle, OrbitCamera, Scene};
use af_render::structure::{DrawOptions, Look, StructureRenderer};
use af_science::data::TrajectoryStream;
use af_science::gif::GifWriter;

use crate::args::{err, load, out, symbol_z, unknown_option, Args, IStream};
use crate::help;

/// The values of every --flag that takes `words` words: one quoted argument
/// ("Cu 0.9 0.5 0.2") or the separate words the help shows (Cu 0.9 0.5 0.2).
fn multi_word_args(args: &Args, flag: &str, words: usize) -> Vec<String> {
    let a = &args.0;
    let mut values = Vec::new();
    for i in 0..a.len() {
        if a[i] != flag || i + 1 >= a.len() {
            continue;
        }
        let mut value = a[i + 1].clone();
        if !value.contains([' ', '\t']) {
            let mut k = 2;
            while k <= words && i + k < a.len() && !a[i + k].starts_with("--") {
                value += " ";
                value += &a[i + k];
                k += 1;
            }
        }
        values.push(value);
    }
    values
}

fn apply_colour_overrides(specs: &[String], colours: &mut [[f32; 3]]) -> Result<(), String> {
    for spec in specs {
        let mut s = IStream::new(spec);
        let (Some(symbol), Some(r), Some(g), Some(b)) = (s.word(), s.float(), s.float(), s.float()) else {
            return Err(format!("Cannot parse --color value '{spec}'. Expected: \"SYMBOL r g b\""));
        };
        let z = symbol_z(&symbol);
        if z <= 0 {
            return Err(format!("Unknown element symbol in --color: '{symbol}'"));
        }
        colours[z as usize] = [r, g, b];
    }
    Ok(())
}

fn apply_radius_overrides(specs: &[String], radii: &mut [f32]) -> Result<(), String> {
    for spec in specs {
        let mut s = IStream::new(spec);
        let parsed = (s.word(), s.float());
        let (Some(symbol), Some(radius)) = parsed else {
            return Err(format!("Cannot parse --radius value '{spec}'. Expected: \"SYMBOL radiusAngstrom\""));
        };
        if !radius.is_finite() || radius <= 0.0 {
            return Err(format!("Cannot parse --radius value '{spec}'. Expected: \"SYMBOL radiusAngstrom\""));
        }
        let z = symbol_z(&symbol);
        if z <= 0 {
            return Err(format!("Unknown element symbol in --radius: '{symbol}'"));
        }
        radii[z as usize] = radius;
    }
    Ok(())
}

/// Standard PNG chunk CRC-32 (polynomial 0xEDB88320, reflected).
fn png_crc32(data: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &byte in data {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = if crc & 1 != 0 { (crc >> 1) ^ 0xEDB8_8320 } else { crc >> 1 };
        }
    }
    crc ^ 0xFFFF_FFFF
}

/// Splices a pHYs chunk (pixels per metre) right after the PNG's IHDR chunk.
fn embed_png_dpi(path: &str, dpi: f64) -> Result<(), String> {
    let mut bytes = std::fs::read(path).unwrap_or_default();
    const SIGNATURE: [u8; 8] = [0x89, b'P', b'N', b'G', b'\r', b'\n', 0x1A, b'\n'];
    let after_ihdr = 8 + 25;
    if bytes.len() < after_ihdr || bytes[..8] != SIGNATURE || &bytes[12..16] != b"IHDR" {
        return Err("Unexpected PNG layout; cannot embed DPI metadata.".into());
    }
    let ppm = (dpi / 0.0254).round() as u32;
    let mut type_and_data = b"pHYs".to_vec();
    type_and_data.extend_from_slice(&ppm.to_be_bytes());
    type_and_data.extend_from_slice(&ppm.to_be_bytes());
    type_and_data.push(1); // unit: metre
    let mut chunk = 9u32.to_be_bytes().to_vec();
    chunk.extend_from_slice(&type_and_data);
    chunk.extend_from_slice(&png_crc32(&type_and_data).to_be_bytes());
    bytes.splice(after_ihdr..after_ihdr, chunk);
    std::fs::write(path, bytes).map_err(|_| format!("Failed to rewrite '{path}' with DPI metadata."))
}

fn write_png(path: &str, width: u32, height: u32, rgba: &[u8]) -> Result<(), String> {
    let file = std::fs::File::create(path).map_err(|_| "Failed to write PNG file.".to_string())?;
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().map_err(|_| "Failed to write PNG file.".to_string())?;
    writer.write_image_data(rgba).map_err(|_| "Failed to write PNG file.".to_string())?;
    writer.finish().map_err(|_| "Failed to write PNG file.".to_string())
}

pub fn run(args: &Args) -> i32 {
    if args.has("--help") || args.has("-h") {
        out(help::RENDER);
        return 0;
    }
    let problem = unknown_option(&args.0, help::RENDER, "--render", "AtomForge --render --help");
    if !problem.is_empty() {
        err(&format!("Error: {problem}\n"));
        return 1;
    }
    match render(args) {
        Ok(code) => code,
        Err(e) => {
            err(&format!("Error: {e}\n"));
            1
        }
    }
}

fn render(args: &Args) -> Result<i32, String> {
    if let Some(types) = args.find("--types") {
        af_io::trajectory::set_lammps_type_elements(types)?;
    }
    let (Some(input), Some(output)) = (args.find("--input"), args.find("--output")) else {
        return Err("--input FILE and --output FILE.png are required".into());
    };
    let width = args.int("--width", 1600)?;
    let height = args.int("--height", 1200)?;
    if width <= 0 || height <= 0 {
        return Err("--width/--height must be positive".into());
    }
    let yaw = args.double("--yaw", 0.0)? as f32;
    let pitch = args.double("--pitch", 0.0)? as f32;
    let roll = args.double("--roll", 0.0)? as f32;
    let explicit_distance = args.find("--distance").is_some();
    let distance = args.double("--distance", 0.0)? as f32;
    let orthographic = args.has("--orthographic");
    let show_bonds = !args.has("--no-bonds");
    let show_box = !args.has("--no-box");
    let radius_scale = args.double("--radius-scale", 1.0)? as f32;
    if radius_scale.is_nan() || radius_scale <= 0.0 {
        return Err("--radius-scale must be positive".into());
    }
    let has_dpi = args.find("--dpi").is_some();
    let dpi = args.double("--dpi", 0.0)?;
    if has_dpi && dpi.is_nan() || has_dpi && dpi <= 0.0 {
        return Err("--dpi must be positive".into());
    }
    // The look: display settings and element tables (factory defaults without --display).
    let mut look: DisplayLook = default_look();
    if let Some(display_path) = args.find("--display") {
        let path = if display_path == "saved" { saved_defaults_path() } else { PathBuf::from(display_path) };
        if !path.exists() {
            return Err(format!("--display file not found: {}", path.display()));
        }
        look = load_defaults(&path);
    }
    let display = look.display;
    let mut background = if display.custom_background { [display.background[0], display.background[1], display.background[2], 1.0] } else { [1.0; 4] };
    if let Some(bg) = multi_word_args(args, "--background", 3).last() {
        let mut s = IStream::new(bg);
        let (Some(r), Some(g), Some(b)) = (s.float(), s.float(), s.float()) else {
            return Err("Cannot parse --background \"r g b\"".into());
        };
        background = [r, g, b, 1.0];
    }
    let lower = output.to_ascii_lowercase();
    let gif_output = lower.len() > 4 && lower.ends_with(".gif");
    let fps = args.double("--fps", 12.0)?;
    if fps.is_nan() || fps <= 0.0 || fps > 100.0 {
        return Err("--fps must be between 0 and 100".into());
    }
    let every = args.int("--every", 1)?;
    if every <= 0 {
        return Err("--every must be positive".into());
    }
    // Trajectory frames are read one at a time while rendering.
    let mut structures: Vec<Structure> = Vec::new();
    let mut trajectory: Option<TrajectoryStream> = None;
    let mut trajectory_error = String::new();
    let mut open_trajectory = |structures: &mut Vec<Structure>, trajectory: &mut Option<TrajectoryStream>| match TrajectoryStream::new(Path::new(input), None)
        .and_then(|t| Ok((t.frame(0)?.structure, t)))
    {
        Ok((first, stream)) => {
            structures.push(first);
            *trajectory = Some(stream);
        }
        Err(e) => {
            trajectory_error = e;
            *trajectory = None;
            structures.clear();
        }
    };
    // A GIF reads trajectory formats first; anything else is one structure.
    if gif_output {
        open_trajectory(&mut structures, &mut trajectory);
    }
    if structures.is_empty() {
        match load(input) {
            Ok(s) => structures.push(s),
            Err(load_error) => {
                // A still image of a trajectory (a LAMMPS dump, say) shows its first frame.
                if !gif_output {
                    open_trajectory(&mut structures, &mut trajectory);
                }
                if structures.is_empty() {
                    let large = af_io::loader::is_large_data_file(input);
                    return Err(format!(
                        "Error loading input structure: {}",
                        if trajectory_error.is_empty() || !large { load_error } else { trajectory_error }
                    ));
                }
                trajectory = None; // one frame: not an animation
            }
        }
    }
    if structures.is_empty() || structures[0].atoms.is_empty() {
        return Err("Error loading input structure: no atoms found".into());
    }
    let trajectory_gif = gif_output && trajectory.as_ref().is_some_and(|t| t.len() > 1);
    let has_frames = args.find("--frames").is_some();
    let frame_count = if trajectory_gif {
        trajectory.as_ref().map_or(1, |t| t.len().div_ceil(every as usize)) as i32
    } else {
        args.int("--frames", if gif_output { 36 } else { 1 })?
    };
    if frame_count <= 0 || (trajectory_gif && has_frames) {
        return Err(if trajectory_gif { "--frames does not apply to trajectory GIFs; use --every".into() } else { "--frames must be positive".into() });
    }
    let yaw_step = args.double("--yaw-step", if gif_output && !trajectory_gif { 360.0 / f64::from(frame_count) } else { 0.0 })? as f32;
    let mut colours = look.colours.clone();
    let mut radii = look.radii.clone();
    apply_colour_overrides(&multi_word_args(args, "--color", 4), &mut colours)?;
    apply_radius_overrides(&multi_word_args(args, "--radius", 2), &mut radii)?;
    // The drawn radius: the element radius times both size factors (bonds
    // are found from the element radii alone, as in the desktop).
    let atom_scale = radius_scale * display.atom_scale;
    let apply_colours = |target: &mut Structure| {
        for atom in target.atoms.iter_mut() {
            if atom.atomic_number >= 0 && (atom.atomic_number as usize) < colours.len() {
                [atom.r, atom.g, atom.b] = colours[atom.atomic_number as usize];
            }
        }
    };
    for s in structures.iter_mut() {
        apply_colours(s);
    }
    let mut renderer = StructureRenderer::new().map_err(|e| format!("cannot create a graphics device ({e}); --render needs a display/GPU driver."))?;
    let shininess = vec![32.0f32; 119];
    let bond_style = BondStyle::from_display(&display);
    let make_scene = |source: &Structure| {
        let mut data = build_instance_data(source, false, &[[1, 0, 0], [0, 1, 0], [0, 0, 1]], &radii, &shininess);
        for scale in data.scales.iter_mut() {
            *scale *= atom_scale;
        }
        darken_occluded(&data.positions, &data.scales, &mut data.colors, display.ambient_occlusion);
        Scene::new(data, &bond_style)
    };
    let mut scene = make_scene(&structures[0]);
    renderer.upload(&scene);
    let mut camera = OrbitCamera { yaw, pitch, roll, ..Default::default() };
    if explicit_distance {
        camera.distance = distance;
    } else {
        apply_default_view(&mut camera, &scene, width, height, true);
    }
    let mut gif = if gif_output { Some(GifWriter::new(Path::new(output), width, height, 2.max((100.0 / fps).round() as i32))?) } else { None };
    let options = DrawOptions {
        background,
        show_bonds,
        show_atoms: true,
        show_box,
        box_colour: if display.custom_cell_colour { Some(display.cell_colour) } else { None },
        selection_colour: display.selection_colour,
    };
    let mut text = String::new();
    for frame in 0..frame_count {
        if trajectory_gif && frame > 0 {
            let stream = trajectory.as_ref().expect("trajectory");
            let mut next = stream.frame(frame as usize * every as usize)?.structure;
            apply_colours(&mut next);
            scene = make_scene(&next);
            renderer.upload(&scene);
        }
        let view = build_frame_view(&camera, &scene, orthographic, width, height, &display);
        let look = Look::new(&display, &view, [background[0], background[1], background[2]]);
        let pixels = renderer.render(&view, &look, &options).map_err(|e| format!("Error rendering frame: {e}"))?;
        if let Some(g) = gif.as_mut() {
            g.add_frame(&pixels, 4)?;
            camera.yaw += yaw_step;
            continue;
        }
        let path = if frame_count > 1 {
            let suffix = format!("-{frame:03}");
            match output.rfind('.') {
                Some(dot) => format!("{}{suffix}{}", &output[..dot], &output[dot..]),
                None => format!("{output}{suffix}"),
            }
        } else {
            output.to_string()
        };
        write_png(&path, width as u32, height as u32, &pixels).map_err(|e| format!("Error writing image: {e}"))?;
        if has_dpi {
            embed_png_dpi(&path, dpi)?;
        }
        out(&format!("Saved to: {path}\n"));
        camera.yaw += yaw_step;
    }
    if let Some(mut g) = gif {
        g.close()?;
        text += &format!("Saved {} frames to: {output}\n", g.frames());
    }
    out(&text);
    Ok(0)
}
