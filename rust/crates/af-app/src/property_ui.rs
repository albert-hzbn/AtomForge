//! What the C++ draws for the colouring modes: the Atom Property Colouring
//! panel (FileBrowser::drawAtomPropertySettings), the property legend at the
//! bottom right (drawAtomPropertyLegend: a viridis bar, or swatches for coded
//! values such as "Structure type (0 other, 1 fcc, 2 hcp)"), and the IPF
//! triangle of crystal-orientation colouring (drawIPFTriangleLegend).

use af_science::atom_properties::{finite_range, viridis, PropertyDisplay};
use eframe::egui;

use crate::theme;

/// The codes a property name lists in parentheses: "Structure type (0 other,
/// 1 fcc)" gives ("Structure type", [(0, "other"), (1, "fcc")]); none when
/// it lists none (propertyCodes).
pub fn property_codes(name: &str) -> (String, Vec<(i64, String)>) {
    let none = (name.to_string(), Vec::new());
    let (Some(open), Some(close)) = (name.rfind('('), name.rfind(')')) else { return none };
    if close < open {
        return none;
    }
    let mut codes = Vec::new();
    for item in name[open + 1..close].split(',') {
        let item = item.trim_start();
        if item.is_empty() {
            return none;
        }
        let digits = item.char_indices().take_while(|(i, c)| c.is_ascii_digit() || (*i == 0 && (*c == '-' || *c == '+'))).last().map(|(i, c)| i + c.len_utf8()).unwrap_or(0);
        if digits == 0 || !item[digits..].starts_with(' ') {
            return none;
        }
        let Ok(code) = item[..digits].parse::<i64>() else { return none };
        let label = item[digits + 1..].trim_end();
        if label.is_empty() {
            return none;
        }
        codes.push((code, label.to_string()));
    }
    (name[..open].trim_end().to_string(), codes)
}

fn colour(c: [f32; 3]) -> egui::Color32 {
    egui::Color32::from_rgb((c[0].clamp(0.0, 1.0) * 255.0) as u8, (c[1].clamp(0.0, 1.0) * 255.0) as u8, (c[2].clamp(0.0, 1.0) * 255.0) as u8)
}

fn card(painter: &egui::Painter, rect: egui::Rect, dark: bool) {
    let p = theme::palette(dark);
    let scale = 18.0 / 16.0;
    painter.rect_filled(rect, 8.0 * scale, p.popup.gamma_multiply(0.9));
    painter.rect_stroke(rect, 8.0 * scale, egui::Stroke::new(1.0, p.border), egui::StrokeKind::Inside);
}

/// The legend card at the bottom right of the view (above the status bar).
pub fn draw_property_legend(painter: &egui::Painter, view: egui::Rect, name: &str, low: f64, high: f64, values: &[f64], dark: bool) {
    if view.width() < 200.0 || view.height() < 200.0 {
        return;
    }
    let p = theme::palette(dark);
    let scale = 18.0 / 16.0;
    let font = egui::FontId::proportional(theme::BASE_FONT);
    let (margin, pad) = (16.0 * scale, 10.0 * scale);
    let wrap = 240.0 * scale;
    let (title, codes) = property_codes(name);
    let mut counts = vec![0usize; codes.len()];
    let mut coded = !codes.is_empty();
    let (mut valid, mut only, mut constant) = (0usize, f64::NAN, true);
    for &v in values {
        if !v.is_finite() {
            continue;
        }
        if valid == 0 {
            only = v;
        } else if v != only {
            constant = false;
        }
        valid += 1;
        if coded {
            match codes.iter().position(|c| c.0 as f64 == v) {
                Some(k) => counts[k] += 1,
                None => coded = false,
            }
        }
    }
    let bottom_right = egui::pos2(view.right() - margin, view.bottom() - margin);
    let span = high - low;
    let colour_of = |v: f64| colour(viridis(if span > 0.0 { (v - low) / span } else { 0.5 }));
    if valid > 0 && (coded || constant) {
        // Coded values or a single value: swatches, not a bar.
        let rows: Vec<(String, egui::Color32)> = if coded {
            codes.iter().zip(&counts).filter(|(_, n)| **n > 0).map(|((code, label), n)| (format!("{label}  {n}"), colour_of(*code as f64))).collect()
        } else {
            vec![(format!("All atoms: {}", af_core::format::format_g(only, 4)), colour_of(only))]
        };
        let (swatch, gap) = (12.0 * scale, 8.0 * scale);
        let title_galley = painter.layout(if coded { title } else { name.to_string() }, font.clone(), p.text, wrap);
        let line = painter.ctx().fonts_mut(|f| f.row_height(&font)) + 4.0;
        let row_width = rows.iter().map(|r| painter.layout_no_wrap(r.0.clone(), font.clone(), p.text).size().x).fold(0.0f32, f32::max);
        let width = title_galley.size().x.max(swatch + gap + row_width) + 2.0 * pad;
        let height = pad * 2.0 + title_galley.size().y + 6.0 * scale + line * rows.len() as f32;
        let rect = egui::Rect::from_min_max(bottom_right - egui::vec2(width, height), bottom_right);
        card(painter, rect, dark);
        painter.galley(rect.min + egui::vec2(pad, pad), title_galley.clone(), p.text);
        let mut y = rect.top() + pad + title_galley.size().y + 6.0 * scale;
        for (label, c) in rows {
            let centre = y + (line - 4.0) * 0.5;
            painter.rect_filled(egui::Rect::from_min_size(egui::pos2(rect.left() + pad, centre - swatch * 0.5), egui::vec2(swatch, swatch)), 3.0 * scale, c);
            painter.text(egui::pos2(rect.left() + pad + swatch + gap, centre), egui::Align2::LEFT_CENTER, label, font.clone(), p.text);
            y += line;
        }
        return;
    }
    // The property name, then the viridis bar with its high and low values beside it.
    let (bar_width, bar_height) = (16.0 * scale, 150.0 * scale);
    let title_galley = painter.layout(name.to_string(), font.clone(), p.text, wrap);
    let (high_text, low_text) = (af_core::format::format_g(high, 4), af_core::format::format_g(low, 4));
    let value_width = [&high_text, &low_text].iter().map(|t| painter.layout_no_wrap(t.to_string(), font.clone(), p.text).size().x).fold(0.0f32, f32::max);
    let width = title_galley.size().x.max(bar_width + 8.0 * scale + value_width) + 2.0 * pad;
    let height = pad * 2.0 + title_galley.size().y + 8.0 * scale + bar_height;
    let rect = egui::Rect::from_min_max(bottom_right - egui::vec2(width, height), bottom_right);
    card(painter, rect, dark);
    painter.galley(rect.min + egui::vec2(pad, pad), title_galley.clone(), p.text);
    let (x0, y1) = (rect.left() + pad, rect.top() + pad + title_galley.size().y + 8.0 * scale);
    let steps = 48;
    let mut mesh = egui::Mesh::default();
    for i in 0..steps {
        let top = colour(viridis(1.0 - i as f64 / steps as f64));
        let bottom = colour(viridis(1.0 - (i + 1) as f64 / steps as f64));
        let ya = y1 + bar_height * i as f32 / steps as f32;
        let yb = y1 + bar_height * (i + 1) as f32 / steps as f32;
        let base = mesh.vertices.len() as u32;
        for (x, y, c) in [(x0, ya, top), (x0 + bar_width, ya, top), (x0 + bar_width, yb, bottom), (x0, yb, bottom)] {
            mesh.colored_vertex(egui::pos2(x, y), c);
        }
        mesh.add_triangle(base, base + 1, base + 2);
        mesh.add_triangle(base, base + 2, base + 3);
    }
    painter.add(egui::Shape::mesh(mesh));
    painter.rect_stroke(egui::Rect::from_min_size(egui::pos2(x0, y1), egui::vec2(bar_width, bar_height)), 0.0, egui::Stroke::new(1.0, p.text_muted), egui::StrokeKind::Inside);
    painter.text(egui::pos2(x0 + bar_width + 8.0 * scale, y1), egui::Align2::LEFT_TOP, high_text, font.clone(), p.text_muted);
    painter.text(egui::pos2(x0 + bar_width + 8.0 * scale, y1 + bar_height), egui::Align2::LEFT_BOTTOM, low_text, font, p.text_muted);
}

/// The IPF colour of a direction in the standard triangle, h >= k >= l >= 0 (ipfDirToColor).
fn ipf_colour(h: f32, k: f32, l: f32) -> egui::Color32 {
    let (r, g, b) = (h - k, k - l, l * 1.7320508);
    let m = r.max(g).max(b).max(1e-6);
    egui::Color32::from_rgb((r / m * 255.0) as u8, (g / m * 255.0) as u8, (b / m * 255.0) as u8)
}

/// The cubic IPF-Z triangle at the bottom right ([001] red, [011] green, [111] blue).
pub fn draw_ipf_legend(painter: &egui::Painter, view: egui::Rect, dark: bool) {
    if view.width() < 200.0 || view.height() < 200.0 {
        return;
    }
    let light = !dark;
    let font = egui::FontId::proportional(theme::BASE_FONT);
    let (tri, margin, label_gap, pad_x, pad_y) = (140.0f32, 20.0f32, 6.0f32, 3.0f32, 1.0f32);
    let size = |t: &str| painter.layout_no_wrap(t.to_string(), font.clone(), egui::Color32::WHITE).size();
    let (title, l001, l011, l111) = ("IPF-Z (cubic)", "[001]", "[011]", "[111]");
    let (title_size, s001, s011, s111) = (size(title), size(l001), size(l011), size(l111));
    let overhang = ((title_size.x - tri) * 0.5).max(0.0);
    let reserve = (s011.x + pad_x * 2.0).max(overhang + 4.0);
    let ox = view.right() - margin - reserve - tri;
    let oy = view.bottom() - margin;
    let v0 = egui::pos2(ox, oy);
    let v1 = egui::pos2(ox + tri, oy);
    let v2 = egui::pos2(ox + tri * 0.58, oy - tri * 0.82);
    let bary_to_dir = |b0: f32, b1: f32, b2: f32| {
        let mut h = b2 * (1.0 / 1.7320508);
        let mut k = b1 * (1.0 / 1.4142136) + b2 * (1.0 / 1.7320508);
        let mut l = b0 + b1 * (1.0 / 1.4142136) + b2 * (1.0 / 1.7320508);
        if l < k {
            std::mem::swap(&mut l, &mut k);
        }
        if l < h {
            std::mem::swap(&mut l, &mut h);
        }
        if k < h {
            std::mem::swap(&mut k, &mut h);
        }
        let len = (h * h + k * k + l * l).sqrt();
        if len > 1e-10 {
            (h / len, k / len, l / len)
        } else {
            (h, k, l)
        }
    };
    let at = |b0: f32, b1: f32, b2: f32| egui::pos2(b0 * v0.x + b1 * v1.x + b2 * v2.x, b0 * v0.y + b1 * v1.y + b2 * v2.y);
    let mut mesh = egui::Mesh::default();
    let mut triangle = |a: egui::Pos2, b: egui::Pos2, c: egui::Pos2, col: egui::Color32| {
        let base = mesh.vertices.len() as u32;
        for q in [a, b, c] {
            mesh.colored_vertex(q, col);
        }
        mesh.add_triangle(base, base + 1, base + 2);
    };
    let steps = 30;
    for i in 0..steps {
        for j in 0..=i {
            let u0 = 1.0 - i as f32 / steps as f32;
            let u1 = j as f32 / steps as f32;
            let u2 = 1.0 - u0 - u1;
            if u2 < -0.001 {
                continue;
            }
            let u0b = 1.0 - (i + 1) as f32 / steps as f32;
            let u1b = j as f32 / steps as f32;
            let u2b = 1.0 - u0b - u1b;
            let u1c = (j + 1) as f32 / steps as f32;
            let u2c = 1.0 - u0b - u1c;
            if u2b >= -0.001 && u2c >= -0.001 {
                let bc0 = (u0 + u0b + u0b) / 3.0;
                let bc1 = (u1 + u1b + u1c) / 3.0;
                let (h, k, l) = bary_to_dir(bc0, bc1, (1.0 - bc0 - bc1).max(0.0));
                triangle(at(u0, u1, u2.max(0.0)), at(u0b, u1b, u2b.max(0.0)), at(u0b, u1c, u2c.max(0.0)), ipf_colour(l, k, h));
            }
            if j < i {
                let u1d = (j + 1) as f32 / steps as f32;
                let u2d = 1.0 - u0 - u1d;
                if u2d >= -0.001 && u2c >= -0.001 {
                    let bc0 = (u0 + u0 + u0b) / 3.0;
                    let bc1 = (u1 + u1d + u1c) / 3.0;
                    let (h, k, l) = bary_to_dir(bc0, bc1, (1.0 - bc0 - bc1).max(0.0));
                    triangle(at(u0, u1, u2.max(0.0)), at(u0, u1d, u2d.max(0.0)), at(u0b, u1c, u2c.max(0.0)), ipf_colour(l, k, h));
                }
            }
        }
    }
    painter.add(egui::Shape::mesh(mesh));
    let border = if light { egui::Color32::from_rgb(60, 60, 60) } else { egui::Color32::from_rgb(200, 200, 200) };
    painter.add(egui::Shape::closed_line(vec![v0, v1, v2], egui::Stroke::new(2.0, border)));
    let label_bg = if light { egui::Color32::from_rgba_unmultiplied(240, 240, 240, 210) } else { egui::Color32::from_rgba_unmultiplied(30, 30, 30, 200) };
    let label_fg = if light { egui::Color32::from_rgb(30, 30, 30) } else { egui::Color32::from_rgb(220, 220, 220) };
    let pos111 = egui::pos2(v2.x - s111.x * 0.5, v2.y - s111.y - label_gap);
    let title_pos = egui::pos2(ox + tri * 0.5 - title_size.x * 0.5, pos111.y - title_size.y - label_gap);
    painter.rect_filled(egui::Rect::from_min_max(title_pos - egui::vec2(4.0, 2.0), title_pos + title_size + egui::vec2(4.0, 2.0)), 3.0, label_bg);
    painter.text(title_pos, egui::Align2::LEFT_TOP, title, font.clone(), label_fg);
    let label = |pos: egui::Pos2, s: egui::Vec2, text: &str, c: egui::Color32| {
        painter.rect_filled(egui::Rect::from_min_max(pos - egui::vec2(pad_x, pad_y), pos + s + egui::vec2(pad_x, pad_y)), 2.0, label_bg);
        painter.text(pos, egui::Align2::LEFT_TOP, text, font.clone(), c);
    };
    label(egui::pos2(v0.x, v0.y - s001.y - label_gap), s001, l001, egui::Color32::from_rgb(255, 80, 80));
    label(egui::pos2(v1.x - s011.x, v1.y - s011.y - label_gap), s011, l011, egui::Color32::from_rgb(80, 255, 80));
    label(pos111, s111, l111, egui::Color32::from_rgb(100, 100, 255));
}

/// What the Atom Property Colouring panel asked for.
#[derive(Default)]
pub struct PanelEvents {
    pub changed: bool,
    pub show_elements: bool,
}

/// The Atom Property Colouring panel at the top left of the view (drawAtomPropertySettings).
pub fn property_panel(ctx: &egui::Context, name: &str, values: &[f64], display: &mut PropertyDisplay, top: f32) -> PanelEvents {
    let mut events = PanelEvents::default();
    let range = finite_range(values);
    let valid = values.iter().filter(|v| v.is_finite()).count();
    theme::window("Atom Property Colouring").id(egui::Id::unique("atom property colouring")).default_pos(egui::pos2(18.0, top)).pivot(egui::Align2::LEFT_TOP).default_size([360.0, 250.0]).show(ctx, |ui| {
        let p = theme::palette(ui.visuals().dark_mode);
        ui.label(if name.is_empty() { "Per-atom property" } else { name });
        ui.label(
            egui::RichText::new(format!("{valid} valid of {} atoms; data range {} to {}", values.len(), af_core::format::format_g(range.low, 4), af_core::format::format_g(range.high, 4)))
                .color(p.text_muted),
        );
        events.changed |= ui.checkbox(&mut display.auto_range, "Automatic colour range").changed();
        if display.auto_range {
            display.range = range;
        }
        ui.add_enabled_ui(!display.auto_range, |ui| {
            for (label, v) in [("Minimum", &mut display.range.low), ("Maximum", &mut display.range.high)] {
                ui.horizontal(|ui| {
                    events.changed |= ui.add(egui::DragValue::new(v).speed(0.01).max_decimals(6)).changed();
                    ui.label(label);
                });
            }
        });
        events.changed |= ui.checkbox(&mut display.hide_outside, "Hide atoms outside the range").changed();
        events.changed |= ui.checkbox(&mut display.hide_invalid, "Hide atoms without a value").changed();
        if ui.button("Show Element Colours").clicked() {
            events.show_elements = true;
        }
    });
    events
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes() {
        let (title, codes) = property_codes("Structure type (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)");
        assert_eq!(title, "Structure type");
        assert_eq!(codes, vec![(0, "other".into()), (1, "fcc".into()), (2, "hcp".into()), (3, "bcc".into()), (4, "ico".into())]);
        assert!(property_codes("Coordination").1.is_empty());
        assert!(property_codes("Strain (percent)").1.is_empty());
        assert!(property_codes("x (1 a, b)").1.is_empty());
    }
}
