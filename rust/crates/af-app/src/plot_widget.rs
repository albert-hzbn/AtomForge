//! Numeric x/y plots of scientific results (the C++ ui/SciencePlot and
//! ui/UnitText): a true x axis with nice ticks, a legend, an optional log y
//! axis, reference lines, labelled markers and a nearest-point readout.

use af_science::plots::PlotSpec;
use eframe::egui;

/// Units for display: "e/A^3" -> "e/Å³", "cutoff (A)" -> "cutoff (Å)",
/// "(Angstrom)" -> "(Å)", "(degrees)" -> "(°)", "2-theta" -> "2θ".
pub fn unit_display(text: &str) -> String {
    let t = text.as_bytes();
    let alnum = |c: u8| c.is_ascii_alphanumeric();
    let powers: [(&str, &str); 5] = [("^-3", "⁻³"), ("^-2", "⁻²"), ("^-1", "⁻¹"), ("^3", "³"), ("^2", "²")];
    let mut out: Vec<u8> = Vec::with_capacity(t.len() + 8);
    let mut i = 0;
    while i < t.len() {
        let lone = i == 0 || !alnum(t[i - 1]);
        // The word: "Angstrom", "Angstroms", "angstrom".
        for word in ["Angstroms", "angstroms", "Angstrom", "angstrom"] {
            let w = word.as_bytes();
            if lone && t[i..].starts_with(w) && (i + w.len() == t.len() || !alnum(t[i + w.len()])) {
                out.extend_from_slice("Å".as_bytes());
                i += w.len();
                break;
            }
        }
        if i >= t.len() {
            break;
        }
        // An Angstrom sign: a lone "A" followed by a power.
        if t[i] == b'A' && lone && i + 1 < t.len() && t[i + 1] == b'^' {
            if let Some((ascii, shown)) = powers.iter().find(|(a, _)| t[i + 1..].starts_with(a.as_bytes())) {
                out.extend_from_slice("Å".as_bytes());
                out.extend_from_slice(shown.as_bytes());
                i += ascii.len() + 1;
                continue;
            }
        }
        // A unit in brackets or per Angstrom: "(A)", "(A, ...", "1/A".
        if t[i] == b'A' && i > 0 && (t[i - 1] == b'(' || t[i - 1] == b'/') && (i + 1 == t.len() || matches!(t[i + 1], b')' | b',' | b' ')) {
            out.extend_from_slice("Å".as_bytes());
            i += 1;
            continue;
        }
        out.push(t[i]);
        i += 1;
    }
    let mut out = String::from_utf8_lossy(&out).into_owned();
    for (from, to) in [("(degrees)", "(°)"), ("(degrees,", "(°,"), ("2-theta", "2θ"), ("Å squared", "Å²"), ("Å cubed", "Å³")] {
        let mut at = out.find(from);
        while let Some(a) = at {
            out.replace_range(a..a + from.len(), to);
            at = out[a + to.len()..].find(from).map(|k| k + a + to.len());
        }
    }
    out
}

/// Round tick values covering [low, high] with about `target` steps (1, 2, 5 x 10^n).
pub fn nice_ticks(low: f64, high: f64, target: i32) -> Vec<f64> {
    let mut ticks = Vec::new();
    if !(high > low) {
        return ticks;
    }
    let raw = (high - low) / target.max(1) as f64;
    let magnitude = 10f64.powf(raw.log10().floor());
    let mut step = magnitude;
    for m in [1.0, 2.0, 5.0, 10.0] {
        if raw <= m * magnitude {
            step = m * magnitude;
            break;
        }
    }
    let mut t = (low / step).ceil() * step;
    while t <= high + step * 1e-9 {
        ticks.push(if t.abs() < step * 1e-9 { 0.0 } else { t });
        t += step;
    }
    ticks
}

/// ImColor::HSV: hue, saturation and value straight to sRGB (ImGui::ColorConvertHSVtoRGB).
fn hsv(h: f32, s: f32, v: f32) -> egui::Color32 {
    if s == 0.0 {
        return crate::theme::srgb(v, v, v);
    }
    let h = (h % 1.0) / (60.0 / 360.0);
    let i = h as i32;
    let f = h - i as f32;
    let (p, q, t) = (v * (1.0 - s), v * (1.0 - s * f), v * (1.0 - s * (1.0 - f)));
    let (r, g, b) = match i {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    };
    crate::theme::srgb(r, g, b)
}

/// The colour of series `s` (golden-ratio hues).
pub fn series_colour(s: usize, light: bool) -> egui::Color32 {
    hsv((0.58 + s as f32 * 0.61803).fract(), 0.7, if light { 0.6 } else { 0.95 })
}

/// Draws a plot `height` points tall across the available width, in a
/// bordered child as the C++ uiPlot::drawSciencePlot.
pub fn science_plot(ui: &mut egui::Ui, spec: &PlotSpec, height: f32) {
    let p = crate::theme::palette(ui.visuals().dark_mode);
    egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| {
        ui.set_width(ui.available_width());
        plot_canvas(ui, spec, (height - 18.0).max(80.0));
    });
}

fn plot_canvas(ui: &mut egui::Ui, spec: &PlotSpec, height: f32) {
    let light = !ui.visuals().dark_mode;
    let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width().max(80.0), height.max(80.0)), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let rgb = egui::Color32::from_rgb;
    let (background, frame, grid, text) = if light {
        (rgb(248, 249, 251), rgb(150, 155, 170), rgb(215, 218, 226), rgb(45, 45, 60))
    } else {
        (rgb(20, 24, 31), rgb(80, 90, 110), rgb(45, 54, 68), rgb(205, 205, 210))
    };
    painter.rect_filled(rect, 0.0, background);
    // The interface's text font, as ImGui draws plot text.
    let font = egui::TextStyle::Body.resolve(ui.style());
    let line_height = ui.fonts_mut(|f| f.row_height(&font));
    let transform = |y: f64| if spec.log_y { if y > 0.0 { y.log10() } else { f64::NAN } } else { y };
    let (mut x_min, mut x_max, mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for s in &spec.series {
        for (&x, &y) in s.x.iter().zip(&s.y) {
            let y = transform(y);
            if x.is_finite() && y.is_finite() {
                x_min = x_min.min(x);
                x_max = x_max.max(x);
                y_min = y_min.min(y);
                y_max = y_max.max(y);
            }
        }
    }
    for &h in &spec.horizontal {
        let y = transform(h);
        if y.is_finite() {
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
    }
    if !(x_max >= x_min) || !(y_max >= y_min) {
        painter.text(rect.min + egui::vec2(10.0, 10.0), egui::Align2::LEFT_TOP, "No data to plot.", font, text);
        return;
    }
    if x_max == x_min {
        x_min -= 0.5;
        x_max += 0.5;
    }
    if y_max == y_min {
        let pad = 1e-12f64.max(y_min.abs() * 0.05);
        y_min -= pad;
        y_max += pad;
    }
    let pad = (y_max - y_min) * 0.06;
    y_min -= pad;
    y_max += pad;
    let g = egui::Rect::from_min_max(egui::pos2(rect.min.x + 64.0, rect.min.y + line_height * 1.6), egui::pos2(rect.max.x - 14.0, rect.max.y - line_height * 2.4));
    let px = |x: f64| g.min.x + ((x - x_min) / (x_max - x_min)) as f32 * g.width();
    let py = |y: f64| g.max.y - ((y - y_min) / (y_max - y_min)) as f32 * g.height();
    let g4 = |v: f64| af_core::format::format_g(v, 4);
    for t in nice_ticks(x_min, x_max, 6) {
        painter.line_segment([egui::pos2(px(t), g.min.y), egui::pos2(px(t), g.max.y)], egui::Stroke::new(1.0, grid));
        painter.text(egui::pos2(px(t), g.max.y + 3.0), egui::Align2::CENTER_TOP, g4(t), font.clone(), text);
    }
    for t in nice_ticks(y_min, y_max, 5) {
        painter.line_segment([egui::pos2(g.min.x, py(t)), egui::pos2(g.max.x, py(t))], egui::Stroke::new(1.0, grid));
        painter.text(egui::pos2(g.min.x - 4.0, py(t)), egui::Align2::RIGHT_CENTER, g4(if spec.log_y { 10f64.powf(t) } else { t }), font.clone(), text);
    }
    painter.rect_stroke(g, 0.0, egui::Stroke::new(1.0, frame), egui::StrokeKind::Middle);
    let title = format!("{}   [{}]", spec.title, unit_display(&spec.y_label));
    let title_galley = painter.layout(title.clone(), font.clone(), text, g.width());
    let title_response = ui.interact(egui::Rect::from_min_size(egui::pos2(g.min.x, rect.min.y), egui::vec2(g.width(), line_height + 2.0)), response.id.with("title"), egui::Sense::hover());
    if title_galley.rows.len() > 1 {
        title_response.on_hover_text(&title);
    }
    painter.with_clip_rect(egui::Rect::from_min_max(egui::pos2(g.min.x, rect.min.y), egui::pos2(g.max.x, rect.min.y + line_height + 2.0))).galley(egui::pos2(g.min.x, rect.min.y + 2.0), title_galley, text);
    painter.text(egui::pos2(g.center().x, g.max.y + line_height + 4.0), egui::Align2::CENTER_TOP, unit_display(&spec.x_label), font.clone(), text);
    let clipped = painter.with_clip_rect(g);
    for &h in &spec.horizontal {
        let y = transform(h);
        if y.is_finite() {
            clipped.line_segment([egui::pos2(g.min.x, py(y)), egui::pos2(g.max.x, py(y))], egui::Stroke::new(1.0, frame));
        }
    }
    // Marker labels on the first of three rows where they do not overlap an earlier one.
    let mut row_end = [f32::MIN; 3];
    for m in &spec.markers {
        let x = px(m.x);
        clipped.line_segment([egui::pos2(x, g.min.y), egui::pos2(x, g.max.y)], egui::Stroke::new(1.0, frame));
        let width = clipped.layout_no_wrap(m.label.clone(), font.clone(), text).size().x;
        for (row, end) in row_end.iter_mut().enumerate() {
            if x + 2.0 <= *end + 4.0 {
                continue;
            }
            clipped.text(egui::pos2(x + 2.0, g.min.y + 2.0 + row as f32 * line_height), egui::Align2::LEFT_TOP, &m.label, font.clone(), text);
            *end = x + 2.0 + width;
            break;
        }
    }
    let mouse = response.hover_pos();
    let mut best: Option<(f32, f64, f64, &str)> = None;
    for (s, series) in spec.series.iter().enumerate() {
        let colour = series_colour(s, light);
        let mut points: Vec<egui::Pos2> = Vec::new();
        let flush = |points: &mut Vec<egui::Pos2>| {
            if points.len() > 1 && !series.points {
                clipped.add(egui::Shape::line(std::mem::take(points), egui::Stroke::new(1.6, colour)));
            }
            points.clear();
        };
        for (&x, &y0) in series.x.iter().zip(&series.y) {
            let y = transform(y0);
            if !x.is_finite() || !y.is_finite() {
                flush(&mut points);
                continue;
            }
            let p = egui::pos2(px(x), py(y));
            points.push(p);
            if series.points {
                clipped.circle_filled(p, 3.0, colour);
            }
            if let Some(m) = mouse {
                let d = (p - m).length_sq();
                if d < 144.0 && best.is_none_or(|b| d < b.0) {
                    best = Some((d, x, y0, &series.name));
                }
            }
        }
        flush(&mut points);
    }
    if spec.series.len() > 1 {
        let mut y = g.min.y + 4.0;
        for (s, series) in spec.series.iter().enumerate() {
            let width = painter.layout_no_wrap(series.name.clone(), font.clone(), text).size().x;
            painter.rect_filled(egui::Rect::from_min_max(egui::pos2(g.max.x - width - 22.0, y + 4.0), egui::pos2(g.max.x - width - 12.0, y + line_height - 4.0)), 0.0, series_colour(s, light));
            painter.text(egui::pos2(g.max.x - width - 6.0, y), egui::Align2::LEFT_TOP, &series.name, font.clone(), text);
            y += line_height;
        }
    }
    if let Some((_, x, y, name)) = best {
        let g6 = |v: f64| af_core::format::format_g(v, 6);
        response.on_hover_text(format!("{name}\n{} = {}\n{} = {}", unit_display(&spec.x_label), g6(x), unit_display(&spec.y_label), g6(y)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn units_as_the_cpp_shows_them() {
        assert_eq!(unit_display("e/A^3"), "e/Å³");
        assert_eq!(unit_display("(e/A^3)^2*A^3"), "(e/Å³)^2*Å³");
        assert_eq!(unit_display("A^-1"), "Å⁻¹");
        assert_eq!(unit_display("DATA^3 and MSD (A^2)"), "DATA^3 and MSD (Å²)");
        assert_eq!(unit_display("Cutoff (A)"), "Cutoff (Å)");
        assert_eq!(unit_display("Neighbor cutoff (Angstrom, 0 = automatic)"), "Neighbor cutoff (Å, 0 = automatic)");
        assert_eq!(unit_display("k (1/A)"), "k (1/Å)");
        assert_eq!(unit_display("Site A, Alpha (A1)"), "Site A, Alpha (A1)");
        assert_eq!(unit_display("2-theta range (degrees)"), "2θ range (°)");
        assert_eq!(unit_display("Debye-Waller B (Angstrom squared)"), "Debye-Waller B (Å²)");
    }

    #[test]
    fn ticks() {
        assert_eq!(nice_ticks(0.0, 10.0, 5), vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
        // 0.61 / 6 rounds up to steps of 0.2.
        assert_eq!(nice_ticks(-0.3, 0.31, 6), vec![-0.2, 0.0, 0.2]);
        assert!(nice_ticks(1.0, 1.0, 5).is_empty());
    }
}
