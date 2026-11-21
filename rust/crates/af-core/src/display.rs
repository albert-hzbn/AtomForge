//! How the 3D view draws a structure (the C++ model/DisplaySettings): light,
//! surface finish, sizes, colours, depth cues, outlines and bond rules, with
//! the "key=value" text form used in projects and the saved defaults.

use std::path::{Path, PathBuf};

use crate::elements;
use crate::format::{format_g, parse_prefix_f32, read_f32_words};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplaySettings {
    pub ambient: f32,
    pub shadow_strength: f32,
    pub contrast: f32,
    pub saturation: f32,
    pub light_azimuth: f32,
    pub light_elevation: f32,
    pub light_follows_camera: bool,
    pub specular_intensity: f32,
    pub shininess_scale: f32,
    pub shininess_floor: f32,
    pub atom_scale: f32,
    pub bond_scale: f32,
    pub uniform_bond_colour: bool,
    pub bond_colour: [f32; 3],
    pub custom_background: bool,
    pub background: [f32; 3],
    pub custom_cell_colour: bool,
    pub cell_colour: [f32; 3],
    pub cell_line_width: f32,
    pub depth_cueing: f32,
    pub ambient_occlusion: f32,
    pub outline_width: f32,
    pub outline_colour: [f32; 3],
    pub selection_colour: [f32; 3],
    pub bond_tolerance: f32,
    pub metal_bonds: bool,
    pub bond_atom_limit: f32,
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            ambient: 0.40,
            shadow_strength: 0.20,
            contrast: 1.25,
            saturation: 1.55,
            light_azimuth: 45.0,
            light_elevation: 35.26439,
            light_follows_camera: false,
            specular_intensity: 0.65,
            shininess_scale: 1.5,
            shininess_floor: 32.0,
            atom_scale: 1.0,
            bond_scale: 1.0,
            uniform_bond_colour: false,
            bond_colour: [0.62; 3],
            custom_background: false,
            background: [1.0; 3],
            custom_cell_colour: false,
            cell_colour: [0.15; 3],
            cell_line_width: 3.5,
            depth_cueing: 0.0,
            ambient_occlusion: 0.0,
            outline_width: 0.0,
            outline_colour: [0.0; 3],
            selection_colour: [1.0, 0.66, 0.12],
            bond_tolerance: 1.18,
            metal_bonds: false,
            bond_atom_limit: 500000.0,
        }
    }
}

/// The allowed range of a number (values read from files are clamped to it).
#[derive(Clone, Copy, Debug)]
pub struct Range {
    pub low: f32,
    pub high: f32,
}
pub const AMBIENT: Range = Range { low: 0.0, high: 1.0 };
pub const SHADOW_STRENGTH: Range = Range { low: 0.0, high: 1.0 };
pub const CONTRAST: Range = Range { low: 0.5, high: 2.5 };
pub const SATURATION: Range = Range { low: 0.0, high: 3.0 };
pub const AZIMUTH: Range = Range { low: 0.0, high: 360.0 };
pub const ELEVATION: Range = Range { low: -85.0, high: 85.0 };
pub const SPECULAR_INTENSITY: Range = Range { low: 0.0, high: 2.0 };
pub const SHININESS_SCALE: Range = Range { low: 0.1, high: 5.0 };
pub const SHININESS_FLOOR: Range = Range { low: 1.0, high: 256.0 };
pub const ATOM_SCALE: Range = Range { low: 0.2, high: 3.0 };
pub const BOND_SCALE: Range = Range { low: 0.2, high: 4.0 };
pub const CELL_LINE_WIDTH: Range = Range { low: 1.0, high: 8.0 };
pub const DEPTH_CUEING: Range = Range { low: 0.0, high: 1.0 };
pub const AMBIENT_OCCLUSION: Range = Range { low: 0.0, high: 1.0 };
pub const OUTLINE_WIDTH: Range = Range { low: 0.0, high: 0.6 };
pub const BOND_TOLERANCE: Range = Range { low: 0.8, high: 1.6 };
pub const BOND_ATOM_LIMIT: Range = Range { low: 0.0, high: 5.0e6 };

const DEGREES: f32 = 3.14159265358979 / 180.0;

fn clamp_to(value: f32, range: Range) -> f32 {
    value.clamp(range.low, range.high)
}

fn near(a: f32, b: f32) -> bool {
    (a - b).abs() <= 1e-4 * 1.0f32.max(b.abs())
}

fn colour_text(c: [f32; 3]) -> String {
    format!("{} {} {}", format_g(c[0] as f64, 9), format_g(c[1] as f64, 9), format_g(c[2] as f64, 9))
}

fn read_colour(text: &str, colour: &mut [f32; 3]) -> bool {
    let Some(values) = read_f32_words(text, 3) else { return false };
    if values.iter().any(|v| !v.is_finite()) {
        return false;
    }
    *colour = [values[0].clamp(0.0, 1.0), values[1].clamp(0.0, 1.0), values[2].clamp(0.0, 1.0)];
    true
}

fn read_number(text: &str, value: &mut f32, range: Range) -> bool {
    let Some((number, _)) = parse_prefix_f32(text) else { return false };
    if !number.is_finite() {
        return false;
    }
    *value = clamp_to(number, range);
    true
}

fn read_switch(text: &str, value: &mut bool) -> bool {
    match text {
        "true" | "1" | "yes" => *value = true,
        "false" | "0" | "no" => *value = false,
        _ => return false,
    }
    true
}

fn trimmed(text: &str) -> &str {
    text.trim_matches(|c| c == ' ' || c == '\t' || c == '\r')
}

/// Unit vector towards the light: in the scene's frame (y up) for a fixed
/// light; in the camera's frame for one following the camera.
pub fn light_direction(s: &DisplaySettings) -> [f32; 3] {
    let azimuth = s.light_azimuth * DEGREES;
    let elevation = clamp_to(s.light_elevation, ELEVATION) * DEGREES;
    let horizontal = elevation.cos();
    if s.light_follows_camera {
        [horizontal * azimuth.sin(), elevation.sin(), horizontal * azimuth.cos()]
    } else {
        [horizontal * azimuth.cos(), elevation.sin(), horizontal * azimuth.sin()]
    }
}

/// A named look: lighting, surface finish, background and cell colours.
#[derive(Clone, Copy, Debug)]
pub struct Preset {
    pub name: &'static str,
    pub description: &'static str,
    pub ambient: f32,
    pub shadow_strength: f32,
    pub contrast: f32,
    pub saturation: f32,
    pub specular_intensity: f32,
    pub shininess_scale: f32,
    pub shininess_floor: f32,
    pub white_background: bool,
}

const fn preset(
    name: &'static str,
    description: &'static str,
    v: [f32; 7],
    white_background: bool,
) -> Preset {
    Preset {
        name,
        description,
        ambient: v[0],
        shadow_strength: v[1],
        contrast: v[2],
        saturation: v[3],
        specular_intensity: v[4],
        shininess_scale: v[5],
        shininess_floor: v[6],
        white_background,
    }
}

pub const PRESETS: [Preset; 7] = [
    preset("Default", "Balanced light with soft shadows and moderate highlights.", [0.40, 0.20, 1.25, 1.55, 0.65, 1.5, 32.0], false),
    preset("Soft", "Even light, faint shadows and gentle highlights.", [0.55, 0.10, 1.05, 1.25, 0.35, 1.0, 16.0], false),
    preset("Glossy", "Polished spheres with small, bright highlights.", [0.35, 0.30, 1.30, 1.60, 1.25, 3.0, 96.0], false),
    preset("Matte", "Chalky surfaces without highlights.", [0.45, 0.20, 1.15, 1.30, 0.05, 0.6, 8.0], false),
    preset("Dramatic", "Low ambient light and deep shadows for strong depth.", [0.15, 0.75, 1.60, 1.70, 0.90, 2.0, 48.0], false),
    preset("Flat", "No shadows or highlights: plain colours, as in a schematic.", [0.75, 0.0, 1.0, 1.30, 0.0, 1.0, 32.0], false),
    preset("Publication", "White background, black cell and soft light, for figures.", [0.50, 0.15, 1.10, 1.35, 0.50, 1.5, 32.0], true),
];

pub fn apply_preset(s: &mut DisplaySettings, p: &Preset) {
    s.ambient = p.ambient;
    s.shadow_strength = p.shadow_strength;
    s.contrast = p.contrast;
    s.saturation = p.saturation;
    s.specular_intensity = p.specular_intensity;
    s.shininess_scale = p.shininess_scale;
    s.shininess_floor = p.shininess_floor;
    s.custom_background = p.white_background;
    s.background = [1.0; 3];
    s.custom_cell_colour = p.white_background;
    s.cell_colour = [0.0; 3];
}

/// The preset the settings match, or -1 for a customised look.
pub fn matching_preset(s: &DisplaySettings) -> i32 {
    for (i, p) in PRESETS.iter().enumerate() {
        let look = near(s.ambient, p.ambient)
            && near(s.shadow_strength, p.shadow_strength)
            && near(s.contrast, p.contrast)
            && near(s.saturation, p.saturation)
            && near(s.specular_intensity, p.specular_intensity)
            && near(s.shininess_scale, p.shininess_scale)
            && near(s.shininess_floor, p.shininess_floor);
        let colours = if p.white_background {
            s.custom_background && s.background == [1.0; 3] && s.custom_cell_colour && s.cell_colour == [0.0; 3]
        } else {
            !s.custom_background && !s.custom_cell_colour
        };
        if look && colours {
            return i as i32;
        }
    }
    -1
}

/// Text form, one "key=value" per line.
pub fn to_text(s: &DisplaySettings) -> String {
    let g = |v: f32| format_g(v as f64, 9);
    let switch = |b: bool| if b { "true" } else { "false" };
    let mut out = String::new();
    let mut line = |key: &str, value: String| {
        out.push_str(key);
        out.push('=');
        out.push_str(&value);
        out.push('\n');
    };
    line("ambient", g(s.ambient));
    line("shadow_strength", g(s.shadow_strength));
    line("contrast", g(s.contrast));
    line("saturation", g(s.saturation));
    line("light_azimuth", g(s.light_azimuth));
    line("light_elevation", g(s.light_elevation));
    line("light_follows_camera", switch(s.light_follows_camera).into());
    line("specular_intensity", g(s.specular_intensity));
    line("shininess_scale", g(s.shininess_scale));
    line("shininess_floor", g(s.shininess_floor));
    line("atom_scale", g(s.atom_scale));
    line("bond_scale", g(s.bond_scale));
    line("bond_colour", if s.uniform_bond_colour { colour_text(s.bond_colour) } else { "atoms".into() });
    line("background", if s.custom_background { colour_text(s.background) } else { "theme".into() });
    line("cell_colour", if s.custom_cell_colour { colour_text(s.cell_colour) } else { "theme".into() });
    line("cell_line_width", g(s.cell_line_width));
    line("depth_cueing", g(s.depth_cueing));
    line("ambient_occlusion", g(s.ambient_occlusion));
    line("outline_width", g(s.outline_width));
    line("outline_colour", colour_text(s.outline_colour));
    line("selection_colour", colour_text(s.selection_colour));
    line("bond_tolerance", g(s.bond_tolerance));
    line("metal_bonds", switch(s.metal_bonds).into());
    line("bond_atom_limit", g(s.bond_atom_limit));
    out
}

/// Settings from text, keeping `base` for missing or unreadable keys and
/// clamping numbers to their ranges.
pub fn from_text(text: &str, base: &DisplaySettings) -> DisplaySettings {
    let mut s = *base;
    for line in text.split('\n') {
        let Some(equals) = line.find('=') else { continue };
        let key = trimmed(&line[..equals]);
        let value = trimmed(&line[equals + 1..]);
        match key {
            "ambient" => _ = read_number(value, &mut s.ambient, AMBIENT),
            "shadow_strength" => _ = read_number(value, &mut s.shadow_strength, SHADOW_STRENGTH),
            "contrast" => _ = read_number(value, &mut s.contrast, CONTRAST),
            "saturation" => _ = read_number(value, &mut s.saturation, SATURATION),
            "light_azimuth" => {
                let mut azimuth = 0.0;
                if read_number(value, &mut azimuth, Range { low: -1e6, high: 1e6 }) {
                    s.light_azimuth = ((azimuth % 360.0) + 360.0) % 360.0;
                }
            }
            "light_elevation" => _ = read_number(value, &mut s.light_elevation, ELEVATION),
            "light_follows_camera" => _ = read_switch(value, &mut s.light_follows_camera),
            "specular_intensity" => _ = read_number(value, &mut s.specular_intensity, SPECULAR_INTENSITY),
            "shininess_scale" => _ = read_number(value, &mut s.shininess_scale, SHININESS_SCALE),
            "shininess_floor" => _ = read_number(value, &mut s.shininess_floor, SHININESS_FLOOR),
            "atom_scale" => _ = read_number(value, &mut s.atom_scale, ATOM_SCALE),
            "bond_scale" => _ = read_number(value, &mut s.bond_scale, BOND_SCALE),
            "bond_colour" => {
                if value == "atoms" {
                    s.uniform_bond_colour = false;
                } else if read_colour(value, &mut s.bond_colour) {
                    s.uniform_bond_colour = true;
                }
            }
            "background" => {
                if value == "theme" {
                    s.custom_background = false;
                } else if read_colour(value, &mut s.background) {
                    s.custom_background = true;
                }
            }
            "cell_colour" => {
                if value == "theme" {
                    s.custom_cell_colour = false;
                } else if read_colour(value, &mut s.cell_colour) {
                    s.custom_cell_colour = true;
                }
            }
            "cell_line_width" => _ = read_number(value, &mut s.cell_line_width, CELL_LINE_WIDTH),
            "depth_cueing" => _ = read_number(value, &mut s.depth_cueing, DEPTH_CUEING),
            "ambient_occlusion" => _ = read_number(value, &mut s.ambient_occlusion, AMBIENT_OCCLUSION),
            "outline_width" => _ = read_number(value, &mut s.outline_width, OUTLINE_WIDTH),
            "outline_colour" => _ = read_colour(value, &mut s.outline_colour),
            "selection_colour" => _ = read_colour(value, &mut s.selection_colour),
            "bond_tolerance" => _ = read_number(value, &mut s.bond_tolerance, BOND_TOLERANCE),
            "metal_bonds" => _ = read_switch(value, &mut s.metal_bonds),
            "bond_atom_limit" => {
                if read_number(value, &mut s.bond_atom_limit, BOND_ATOM_LIMIT) {
                    s.bond_atom_limit = s.bond_atom_limit.round();
                }
            }
            _ => {}
        }
    }
    s
}

/// Element colours and radii that differ from the built-in tables.
pub fn element_text(colours: &[[f32; 3]], radii: &[f32]) -> String {
    let default_colours = elements::default_colours();
    let default_radii = elements::literature_covalent_radii();
    let mut out = String::new();
    for z in 1..119usize {
        if z < colours.len() && colours[z] != default_colours[z] {
            out.push_str(&format!("colour.{}={}\n", elements::symbol(z as i32), colour_text(colours[z])));
        }
        if z < radii.len() && radii[z] != default_radii[z] {
            out.push_str(&format!("radius.{}={}\n", elements::symbol(z as i32), format_g(radii[z] as f64, 9)));
        }
    }
    out
}

pub fn read_element_text(text: &str, colours: &mut [[f32; 3]], radii: &mut [f32]) {
    for line in text.split('\n') {
        let Some(equals) = line.find('=') else { continue };
        let key = trimmed(&line[..equals]);
        let value = trimmed(&line[equals + 1..]);
        let Some(dot) = key.find('.') else { continue };
        let kind = &key[..dot];
        let z = elements::atomic_number(&key[dot + 1..]);
        if z == 0 {
            continue;
        }
        let k = z as usize;
        if kind == "colour" && k < colours.len() {
            read_colour(value, &mut colours[k]);
        } else if kind == "radius" && k < radii.len() {
            read_number(value, &mut radii[k], Range { low: 0.05, high: 5.0 });
        }
    }
}

/// The whole look of a window: display settings with element colours and radii.
#[derive(Clone, Debug, PartialEq)]
pub struct Look {
    pub display: DisplaySettings,
    pub colours: Vec<[f32; 3]>,
    pub radii: Vec<f32>,
}

pub fn default_look() -> Look {
    Look { display: DisplaySettings::default(), colours: elements::default_colours(), radii: elements::literature_covalent_radii() }
}

pub fn look_text(look: &Look) -> String {
    to_text(&look.display) + &element_text(&look.colours, &look.radii)
}

pub fn read_look(text: &str) -> Look {
    let mut look = default_look();
    look.display = from_text(text, &DisplaySettings::default());
    read_element_text(text, &mut look.colours, &mut look.radii);
    look
}

/// The user's saved default (~/.atomforge_display.txt).
pub fn saved_defaults_path() -> PathBuf {
    let home = if cfg!(windows) { std::env::var_os("USERPROFILE") } else { std::env::var_os("HOME") };
    let base = home.filter(|h| !h.is_empty()).map(PathBuf::from).unwrap_or_else(|| PathBuf::from("."));
    base.join(".atomforge_display.txt")
}

pub fn save_defaults(look: &Look, path: &Path) -> bool {
    let text = format!("# AtomForge display settings (Settings > Display Settings > Save as My Default)\n{}", look_text(look));
    std::fs::write(path, text).is_ok()
}

pub fn load_defaults(path: &Path) -> Look {
    match std::fs::read_to_string(path) {
        Ok(text) => read_look(&text),
        Err(_) => default_look(),
    }
}
