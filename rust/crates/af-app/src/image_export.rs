//! File > Export Image (the C++ app/ImageExport): PNG, JPEG (quality 95) and
//! SVG (atoms as circles, bonds as lines, painted back to front), at a
//! multiple of the view's size, with or without the background, and with the
//! orientation gizmo drawn into the pixels by the C++ software rasteriser.

use af_core::structure::Structure;
use af_render::scene::Scene;
use glam::{Mat3, Mat4, Vec3, Vec4};

/// The file formats of Export Image.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Png,
    Jpeg,
    Svg,
}

/// The format a file name asks for (PNG unless it ends in .jpg/.jpeg or .svg).
pub fn format_of(path: &str) -> Format {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".jpg") || lower.ends_with(".jpeg") {
        Format::Jpeg
    } else if lower.ends_with(".svg") {
        Format::Svg
    } else {
        Format::Png
    }
}

/// The camera of an export (ImageExportView): size in pixels and the glm
/// matrices of the frame.
#[derive(Clone, Copy, Debug)]
pub struct ExportView {
    pub width: i32,
    pub height: i32,
    pub projection: Mat4,
    pub view: Mat4,
}

fn clamp_to_byte(value: f32) -> i32 {
    (value.clamp(0.0, 1.0) * 255.0).round() as i32
}

/// "%.3f".
pub fn svg_float(value: f32) -> String {
    format!("{value:.3}")
}

/// "rgb(r,g,b)" with each channel rounded to a byte.
pub fn svg_rgb(c: [f32; 3]) -> String {
    format!("rgb({},{},{})", clamp_to_byte(c[0]), clamp_to_byte(c[1]), clamp_to_byte(c[2]))
}

fn is_orthographic(projection: &Mat4) -> bool {
    (projection.w_axis.w - 1.0).abs() <= 1e-4
}

fn pixels_per_world_unit(view: &ExportView, view_depth: f32) -> f32 {
    let scale = view.projection.y_axis.y.abs();
    if is_orthographic(&view.projection) {
        return 0.5 * view.height as f32 * scale;
    }
    0.5 * view.height as f32 * scale / view_depth.max(1e-4)
}

/// (screen x, screen y, NDC depth, view depth), or None when clipped.
fn project_point(world: Vec3, view: &ExportView) -> Option<(f32, f32, f32, f32)> {
    let clip = view.projection * view.view * Vec4::new(world.x, world.y, world.z, 1.0);
    if clip.w <= 1e-6 {
        return None;
    }
    let ndc = clip.truncate() / clip.w;
    if ndc.z < -1.0 || ndc.z > 1.0 {
        return None;
    }
    let x = (ndc.x * 0.5 + 0.5) * view.width as f32;
    let y = (1.0 - (ndc.y * 0.5 + 0.5)) * view.height as f32;
    let view_pos = view.view * Vec4::new(world.x, world.y, world.z, 1.0);
    Some((x, y, ndc.z, -view_pos.z))
}

/// The SVG of a scene (writeSvgImage).
pub fn svg_text(view: &ExportView, background: [f32; 4], include_background: bool, show_bonds: bool, show_atoms: bool, scene: &Scene) -> String {
    struct Bond {
        depth: f32,
        x1: f32,
        y1: f32,
        x2: f32,
        y2: f32,
        stroke_width: f32,
        colour: [f32; 3],
    }
    struct Atom {
        depth: f32,
        x: f32,
        y: f32,
        radius: f32,
        colour: [f32; 3],
    }
    let mut bonds: Vec<Bond> = Vec::new();
    let mut atoms: Vec<Atom> = Vec::new();
    if show_bonds && show_atoms {
        let b = &scene.bonds;
        let count = b.starts.len().min(b.ends.len()).min(b.colors_a.len().min(b.colors_b.len()));
        for i in 0..count {
            let Some((x1, y1, _, _)) = project_point(b.starts[i], view) else { continue };
            let Some((x2, y2, _, _)) = project_point(b.ends[i], view) else { continue };
            let mid = 0.5 * (b.starts[i] + b.ends[i]);
            let Some((_, _, z_mid, depth_mid)) = project_point(mid, view) else { continue };
            let radius = b.radii.get(i).copied().unwrap_or(0.08);
            let stroke_width = (radius * pixels_per_world_unit(view, depth_mid) * 2.0).max(1.0);
            let (a, c) = (b.colors_a[i], b.colors_b[i]);
            let colour = [0.5 * (a[0] + c[0]), 0.5 * (a[1] + c[1]), 0.5 * (a[2] + c[2])];
            bonds.push(Bond { depth: z_mid, x1, y1, x2, y2, stroke_width, colour });
        }
        af_core::cxx::sort_by(&mut bonds, |a, b| a.depth > b.depth);
    }
    if show_atoms {
        let d = &scene.instances;
        let count = d.positions.len().min(d.scales.len()).min(d.colors.len());
        for i in 0..count {
            let Some((x, y, z, depth)) = project_point(d.positions[i], view) else { continue };
            let radius = (d.scales[i] * pixels_per_world_unit(view, depth)).max(1.0);
            atoms.push(Atom { depth: z, x, y, radius, colour: d.colors[i] });
        }
        af_core::cxx::sort_by(&mut atoms, |a, b| a.depth > b.depth);
    }
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">\n",
        w = view.width,
        h = view.height
    ));
    if include_background {
        out.push_str(&format!(
            "  <rect x=\"0\" y=\"0\" width=\"{}\" height=\"{}\" fill=\"{}\" />\n",
            view.width,
            view.height,
            svg_rgb([background[0], background[1], background[2]])
        ));
    }
    for b in &bonds {
        out.push_str(&format!(
            "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"{}\" stroke-linecap=\"round\" />\n",
            svg_float(b.x1),
            svg_float(b.y1),
            svg_float(b.x2),
            svg_float(b.y2),
            svg_rgb(b.colour),
            svg_float(b.stroke_width)
        ));
    }
    for a in &atoms {
        let stroke = [(a.colour[0] * 0.65).clamp(0.0, 1.0), (a.colour[1] * 0.65).clamp(0.0, 1.0), (a.colour[2] * 0.65).clamp(0.0, 1.0)];
        out.push_str(&format!(
            "  <circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"1\" />\n",
            svg_float(a.x),
            svg_float(a.y),
            svg_float(a.radius),
            svg_rgb(a.colour),
            svg_rgb(stroke)
        ));
    }
    out.push_str("</svg>\n");
    out
}

pub fn write_svg(path: &str, view: &ExportView, background: [f32; 4], include_background: bool, show_bonds: bool, show_atoms: bool, scene: &Scene) -> Result<(), String> {
    std::fs::write(path, svg_text(view, background, include_background, show_bonds, show_atoms, scene)).map_err(|_| "Failed to open SVG output file.".to_string())
}

// ---------------------------------------------------------------------------
// Software gizmo rasteriser (drawGizmoOnPixels)
// ---------------------------------------------------------------------------

fn blend(pixels: &mut [u8], width: i32, height: i32, px: i32, py: i32, r: u8, g: u8, b: u8, alpha: u8) {
    if px < 0 || py < 0 || px >= width || py >= height {
        return;
    }
    let i = (py as usize * width as usize + px as usize) * 4;
    let fa = alpha as f32 / 255.0;
    let fb = 1.0 - fa;
    pixels[i] = (pixels[i] as f32 * fb + r as f32 * fa) as u8;
    pixels[i + 1] = (pixels[i + 1] as f32 * fb + g as f32 * fa) as u8;
    pixels[i + 2] = (pixels[i + 2] as f32 * fb + b as f32 * fa) as u8;
    pixels[i + 3] = 255;
}

#[allow(clippy::too_many_arguments)]
fn soft_circle(pixels: &mut [u8], width: i32, height: i32, cx: f32, cy: f32, radius: f32, r: u8, g: u8, b: u8, a: u8) {
    let x0 = 0.max((cx - radius - 1.5) as i32);
    let x1 = (width - 1).min((cx + radius + 1.5) as i32);
    let y0 = 0.max((cy - radius - 1.5) as i32);
    let y1 = (height - 1).min((cy + radius + 1.5) as i32);
    for y in y0..=y1 {
        for x in x0..=x1 {
            let (dx, dy) = (x as f32 - cx, y as f32 - cy);
            let dist = (dx * dx + dy * dy).sqrt();
            if dist <= radius + 0.5 {
                let k = (radius + 0.5 - dist).min(1.0);
                blend(pixels, width, height, x, y, r, g, b, (a as f32 * k) as u8);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn soft_line(pixels: &mut [u8], width: i32, height: i32, x0: f32, y0: f32, x1: f32, y1: f32, thickness: f32, r: u8, g: u8, b: u8, a: u8) {
    let (dx, dy) = (x1 - x0, y1 - y0);
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-3 {
        return;
    }
    let (nx, ny) = (-dy / len, dx / len);
    let half = thickness * 0.5;
    let bx0 = 0.max((x0.min(x1) - half - 1.5) as i32);
    let bx1 = (width - 1).min((x0.max(x1) + half + 1.5) as i32);
    let by0 = 0.max((y0.min(y1) - half - 1.5) as i32);
    let by1 = (height - 1).min((y0.max(y1) + half + 1.5) as i32);
    let inv_len2 = 1.0 / (len * len);
    for y in by0..=by1 {
        for x in bx0..=bx1 {
            let (px, py) = (x as f32 - x0, y as f32 - y0);
            let t = (px * dx + py * dy) * inv_len2;
            let perp = (px * nx + py * ny).abs();
            if (-0.01..=1.01).contains(&t) && perp <= half + 0.5 {
                let k = (half + 0.5 - perp).min(1.0);
                blend(pixels, width, height, x, y, r, g, b, (a as f32 * k) as u8);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn soft_ring(pixels: &mut [u8], width: i32, height: i32, cx: f32, cy: f32, radius: f32, ring: f32, r: u8, g: u8, b: u8, a: u8) {
    let x0 = 0.max((cx - radius - ring - 1.5) as i32);
    let x1 = (width - 1).min((cx + radius + ring + 1.5) as i32);
    let y0 = 0.max((cy - radius - ring - 1.5) as i32);
    let y1 = (height - 1).min((cy + radius + ring + 1.5) as i32);
    for y in y0..=y1 {
        for x in x0..=x1 {
            let (dx, dy) = (x as f32 - cx, y as f32 - cy);
            let from_ring = ((dx * dx + dy * dy).sqrt() - radius).abs();
            if from_ring <= ring * 0.5 + 0.5 {
                let k = (ring * 0.5 + 0.5 - from_ring).min(1.0);
                blend(pixels, width, height, x, y, r, g, b, (a as f32 * k) as u8);
            }
        }
    }
}

/// The orientation gizmo at the bottom left of the image (drawGizmoOnPixels).
pub fn draw_gizmo(pixels: &mut [u8], width: i32, height: i32, view: &Mat4, light: bool, scale: i32, structure: &Structure) {
    // 5x7 bitmap glyphs (rows, most significant bit at the left).
    const X: [u8; 7] = [0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001];
    const Y: [u8; 7] = [0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100, 0b00100];
    const Z: [u8; 7] = [0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b11111];
    const A: [u8; 7] = [0b00000, 0b01110, 0b00001, 0b01111, 0b10001, 0b10001, 0b01111];
    const B: [u8; 7] = [0b10000, 0b10000, 0b11110, 0b10001, 0b10001, 0b10001, 0b11110];
    const C: [u8; 7] = [0b00000, 0b01110, 0b10000, 0b10000, 0b10000, 0b10001, 0b01110];
    let lattice = structure.has_unit_cell;
    let glyphs: [&[u8; 7]; 3] = if lattice { [&A, &B, &C] } else { [&X, &Y, &Z] };
    let draw_glyph = |pixels: &mut [u8], ox: f32, oy: f32, glyph: &[u8; 7], r: u8, g: u8, b: u8, ps: f32| {
        let cell = (ps * 1.5).max(1.0);
        for (row, bits) in glyph.iter().enumerate() {
            for col in 0..5 {
                if bits & (1 << (4 - col)) != 0 {
                    let px = ox + col as f32 * cell;
                    let py = oy + row as f32 * cell;
                    let (ix0, iy0) = (px as i32, py as i32);
                    let ix1 = ((px + cell) as i32).min(ix0 + cell as i32 + 1);
                    let iy1 = ((py + cell) as i32).min(iy0 + cell as i32 + 1);
                    for y in iy0..=iy1 {
                        for x in ix0..=ix1 {
                            blend(pixels, width, height, x, y, r, g, b, 255);
                        }
                    }
                }
            }
        }
    };
    let s = scale.max(1) as f32;
    let bg_radius = 53.0 * s;
    let (cx, cy) = (73.0 * s, height as f32 - 73.0 * s);
    let proj_scale = 53.0 * s;
    let axis_length = 0.77f32;
    // Backdrop: a disk and two rings.
    let (bg, bg_a) = if light { ([235u8, 238, 242], 200u8) } else { ([18u8, 22, 28], 190u8) };
    soft_circle(pixels, width, height, cx, cy, bg_radius, bg[0], bg[1], bg[2], bg_a);
    let ring = if light { 0u8 } else { 255u8 };
    soft_ring(pixels, width, height, cx, cy, bg_radius, s, ring, ring, ring, 55);
    soft_ring(pixels, width, height, cx, cy, bg_radius - 7.0 * s, s, ring, ring, ring, 20);
    let colours: [[u8; 3]; 3] = if light { [[200, 50, 50], [40, 160, 50], [40, 110, 220]] } else { [[235, 92, 92], [110, 220, 120], [110, 175, 255]] };
    let outline: [u8; 3] = if light { [30, 30, 40] } else { [245, 245, 250] };
    let world: [Vec3; 3] = if lattice {
        let v = |k: usize| Vec3::new(structure.cell_vectors[k][0] as f32, structure.cell_vectors[k][1] as f32, structure.cell_vectors[k][2] as f32).normalize();
        [v(0), v(1), v(2)]
    } else {
        [Vec3::X, Vec3::Y, Vec3::Z]
    };
    let rotation = Mat3::from_mat4(*view);
    let mut axes: Vec<(f32, f32, f32, usize)> = (0..3)
        .map(|i| {
            let dir = (rotation * world[i]).normalize();
            (cx + dir.x * axis_length * proj_scale, cy - dir.y * axis_length * proj_scale, -dir.z, i)
        })
        .collect();
    // Painter's order: the deepest first.
    af_core::cxx::sort_by(&mut axes, |a, b| a.2 < b.2);
    let glyph_scale = (s * 1.6).max(1.0);
    let (glyph_w, glyph_h) = (5.0 * glyph_scale * 1.5, 7.0 * glyph_scale * 1.5);
    for &(sx, sy, depth, i) in &axes {
        let c = colours[i];
        let depth01 = ((depth + 1.0) * 0.5).clamp(0.0, 1.0);
        let thickness = (1.8 + 1.2 * depth01) * s;
        let dot = (2.5 + 2.0 * depth01) * s;
        soft_line(pixels, width, height, cx, cy, sx, sy, thickness, c[0], c[1], c[2], 255);
        soft_circle(pixels, width, height, sx, sy, dot, c[0], c[1], c[2], 255);
        soft_ring(pixels, width, height, sx, sy, dot + s, s, outline[0], outline[1], outline[2], 200);
        draw_glyph(pixels, sx + 5.5 * s - glyph_w * 0.5, sy + 4.0 * s - glyph_h * 0.5, glyphs[i], c[0], c[1], c[2], glyph_scale);
    }
    let centre: [u8; 3] = if light { [20, 20, 25] } else { [240, 240, 245] };
    soft_circle(pixels, width, height, cx, cy, 3.2 * s, centre[0], centre[1], centre[2], 240);
}

/// Writes rendered RGBA pixels (rows top first) as PNG or JPEG (writeRasterImage).
pub fn write_raster(path: &str, format: Format, rgba: &[u8], width: u32, height: u32, include_background: bool) -> Result<(), String> {
    if format == Format::Png {
        let write = || -> Result<(), String> {
            let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
            let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), width, height);
            encoder.set_color(png::ColorType::Rgba);
            encoder.set_depth(png::BitDepth::Eight);
            let mut writer = encoder.write_header().map_err(|e| e.to_string())?;
            writer.write_image_data(rgba).map_err(|e| e.to_string())
        };
        return write().map_err(|_| "Failed to write PNG file.".to_string());
    }
    // JPEG has no alpha: without the background, pixels are composited over white.
    let mut rgb = Vec::with_capacity(width as usize * height as usize * 3);
    for p in rgba.chunks_exact(4) {
        if include_background {
            rgb.extend_from_slice(&p[..3]);
        } else {
            let a = p[3] as f32 / 255.0;
            for &c in &p[..3] {
                rgb.push((c as f32 * a + 255.0 * (1.0 - a)).round() as u8);
            }
        }
    }
    let write = || -> Result<(), String> {
        let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
        let mut encoder = image::codecs::jpeg::JpegEncoder::new_with_quality(std::io::BufWriter::new(file), 95);
        encoder.encode(&rgb, width, height, image::ExtendedColorType::Rgb8).map_err(|e| e.to_string())
    };
    write().map_err(|_| "Failed to write JPEG file.".to_string())
}

/// The width and height of an exported PNG, JPEG or SVG (UI script checks).
pub fn image_size(path: &str) -> Result<(u32, u32), String> {
    match format_of(path) {
        Format::Png => {
            let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
            let reader = png::Decoder::new(std::io::BufReader::new(file)).read_info().map_err(|e| e.to_string())?;
            let info = reader.info();
            Ok((info.width, info.height))
        }
        Format::Jpeg => {
            let reader = image::ImageReader::open(path).map_err(|e| e.to_string())?.with_guessed_format().map_err(|e| e.to_string())?;
            reader.into_dimensions().map_err(|e| e.to_string())
        }
        Format::Svg => {
            let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
            let attribute = |name: &str| -> Option<u32> {
                let at = text.find(&format!(" {name}=\""))? + name.len() + 3;
                text[at..].split('"').next()?.parse().ok()
            };
            Ok((attribute("width").ok_or("no width")?, attribute("height").ok_or("no height")?))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_render::scene::{BondStyle, InstanceData};

    #[test]
    fn formats_and_numbers() {
        assert_eq!(format_of("a.PNG"), Format::Png);
        assert_eq!(format_of("a.jpeg"), Format::Jpeg);
        assert_eq!(format_of("x/y.svg"), Format::Svg);
        assert_eq!(svg_float(1.23456), "1.235");
        assert_eq!(svg_float(-0.0001), "-0.000");
        assert_eq!(svg_rgb([1.2, 0.5, -1.0]), "rgb(255,128,0)");
    }

    /// Two atoms and their bond seen down -z by an orthographic camera:
    /// the C++ writer's output, worked through by hand.
    #[test]
    fn svg_of_a_tiny_scene() {
        let instances = InstanceData {
            positions: vec![Vec3::new(-1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, -1.0)],
            colors: vec![[1.0, 0.0, 0.0], [0.0, 0.5, 1.0]],
            scales: vec![0.5, 0.25],
            bond_radii: vec![0.1, 0.1],
            shininess: vec![32.0; 2],
            atomic_numbers: vec![8, 1],
            box_lines: Vec::new(),
            atom_indices: vec![0, 1],
            orbit_center: Vec3::ZERO,
        };
        let mut scene = Scene::new(instances, &BondStyle::default());
        scene.bonds.starts = vec![Vec3::new(-1.0, 0.0, 0.0)];
        scene.bonds.ends = vec![Vec3::new(1.0, 0.0, -1.0)];
        scene.bonds.colors_a = vec![[1.0, 0.0, 0.0]];
        scene.bonds.colors_b = vec![[0.0, 0.5, 1.0]];
        scene.bonds.radii = vec![0.1];
        // glm::ortho(-2, 2, -1, 1, 0.1, 10) and a camera at z = 5 looking down -z.
        let projection = Mat4::orthographic_rh_gl(-2.0, 2.0, -1.0, 1.0, 0.1, 10.0);
        let view = Mat4::look_at_rh(Vec3::new(0.0, 0.0, 5.0), Vec3::ZERO, Vec3::Y);
        let v = ExportView { width: 200, height: 100, projection, view };
        let svg = svg_text(&v, [0.8, 0.8, 0.8, 1.0], true, true, true, &scene);
        let expected = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"200\" height=\"100\" viewBox=\"0 0 200 100\">\n  \
<rect x=\"0\" y=\"0\" width=\"200\" height=\"100\" fill=\"rgb(204,204,204)\" />\n  \
<line x1=\"50.000\" y1=\"50.000\" x2=\"150.000\" y2=\"50.000\" stroke=\"rgb(128,64,128)\" stroke-width=\"10.000\" stroke-linecap=\"round\" />\n  \
<circle cx=\"150.000\" cy=\"50.000\" r=\"12.500\" fill=\"rgb(0,128,255)\" stroke=\"rgb(0,83,166)\" stroke-width=\"1\" />\n  \
<circle cx=\"50.000\" cy=\"50.000\" r=\"25.000\" fill=\"rgb(255,0,0)\" stroke=\"rgb(166,0,0)\" stroke-width=\"1\" />\n\
</svg>\n";
        assert_eq!(svg, expected);
        // Without the background there is no rect.
        assert!(!svg_text(&v, [0.9; 4], false, true, true, &scene).contains("<rect"));
    }

    /// AF_SVG_FORMAT_REFERENCE=FILE: lines "HEX %.3f HEX lround" written by
    /// the C++ svgFloat / clampToByte, compared with these functions.
    #[test]
    fn formatting_matches_cpp_reference() {
        let Ok(path) = std::env::var("AF_SVG_FORMAT_REFERENCE") else { return };
        let hex = |t: &str| -> f32 {
            let (sign, t) = if let Some(r) = t.strip_prefix('-') { (-1.0f64, r) } else { (1.0, t) };
            let t = t.trim_start_matches("0x");
            let (mantissa, exponent) = t.split_once('p').unwrap();
            let (int, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
            let mut m = u64::from_str_radix(int, 16).unwrap() as f64;
            for (k, c) in frac.chars().enumerate() {
                m += c.to_digit(16).unwrap() as f64 / 16f64.powi(k as i32 + 1);
            }
            (sign * m * 2f64.powi(exponent.parse::<i32>().unwrap())) as f32
        };
        let text = std::fs::read_to_string(path).unwrap();
        let mut n = 0;
        for line in text.lines() {
            let w: Vec<&str> = line.split(' ').collect();
            let (v, c) = (hex(w[0]), hex(w[2]));
            assert_eq!(svg_float(v), w[1], "{line}");
            assert_eq!(clamp_to_byte(c).to_string(), w[3], "{line}");
            n += 1;
        }
        assert!(n > 0);
    }

    #[test]
    fn gizmo_marks_the_corner() {
        let (w, h) = (300, 300);
        let mut pixels = vec![0u8; (w * h * 4) as usize];
        draw_gizmo(&mut pixels, w, h, &Mat4::IDENTITY, true, 1, &Structure::default());
        // The backdrop at the gizmo's centre, nothing far from it.
        let at = |x: i32, y: i32| pixels[((y * w + x) * 4 + 3) as usize];
        assert_eq!(at(73, h - 73), 255);
        assert_eq!(at(250, 20), 0);
    }
}
