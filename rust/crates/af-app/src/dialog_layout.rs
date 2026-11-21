//! Dialog building blocks shared by the builder and edit dialogs, after the
//! C++ ui/DialogLayout.h and the ImGui widgets the C++ dialogs use: labelled
//! form rows, input fields with − / + step buttons (ImGui InputInt/InputFloat),
//! callout notes, help markers, framed panels and right-aligned footers.

use eframe::egui;

use crate::theme::{self, FRAME_HEIGHT};

/// The kind of a callout (dialogLayout::Note).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Note {
    Info,
    Success,
    Warning,
    Error,
}

pub fn note_colour(ui: &egui::Ui, kind: Note) -> egui::Color32 {
    let dark = ui.visuals().dark_mode;
    match kind {
        Note::Info => theme::palette(dark).accent,
        Note::Success => {
            if dark {
                egui::Color32::from_rgb(110, 200, 130)
            } else {
                egui::Color32::from_rgb(30, 130, 60)
            }
        }
        Note::Warning => ui.visuals().warn_fg_color,
        Note::Error => ui.visuals().error_fg_color,
    }
}

/// A callout: a tinted, rounded box with a coloured left edge and wrapped text (dialogLayout::note).
pub fn note(ui: &mut egui::Ui, kind: Note, text: &str) {
    let colour = note_colour(ui, kind);
    let tint = egui::Color32::from_rgba_unmultiplied(colour.r(), colour.g(), colour.b(), if ui.visuals().dark_mode { 41 } else { 23 });
    let icon = match kind {
        Note::Info => theme::icon::INFO,
        Note::Success => theme::icon::CHECK,
        Note::Warning => theme::icon::WARNING,
        Note::Error => theme::icon::ERROR,
    };
    let response = egui::Frame::new().fill(tint).corner_radius(6).inner_margin(egui::Margin { left: 12, right: 9, top: 9, bottom: 9 }).show(ui, |ui| {
        ui.set_width(ui.available_width());
        ui.horizontal_top(|ui| {
            ui.label(egui::RichText::new(icon).color(colour));
            ui.add(egui::Label::new(text).wrap());
        });
    });
    let r = response.response.rect;
    ui.painter().rect_filled(egui::Rect::from_min_max(r.min, egui::pos2(r.min.x + 3.0, r.max.y)), egui::CornerRadius { nw: 6, sw: 6, ne: 0, se: 0 }, colour);
}


/// Muted text (ImGui::TextDisabled).
pub fn muted(ui: &mut egui::Ui, text: impl Into<String>) -> egui::Response {
    let muted = theme::palette(ui.visuals().dark_mode).text_muted;
    ui.label(egui::RichText::new(text.into()).color(muted))
}

/// A semibold section heading over a thin rule (dialogLayout::section).
pub fn section(ui: &mut egui::Ui, title: &str) {
    theme::section(ui, title);
}


/// A form row: the label at the left (ImGui text before a field), the widget after it.
pub fn row<R>(ui: &mut egui::Ui, label: &str, label_width: f32, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    ui.horizontal(|ui| {
        let (rect, _) = ui.allocate_exact_size(egui::vec2(label_width, FRAME_HEIGHT), egui::Sense::hover());
        let font = egui::TextStyle::Body.resolve(ui.style());
        ui.painter().text(egui::pos2(rect.left(), rect.center().y), egui::Align2::LEFT_CENTER, label, font, ui.visuals().text_color());
        add(ui)
    })
    .inner
}

/// A form row with a tooltip on its label.
pub fn row_hint<R>(ui: &mut egui::Ui, label: &str, hint: &str, label_width: f32, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    ui.horizontal(|ui| {
        let (rect, response) = ui.allocate_exact_size(egui::vec2(label_width, FRAME_HEIGHT), egui::Sense::hover());
        let font = egui::TextStyle::Body.resolve(ui.style());
        ui.painter().text(egui::pos2(rect.left(), rect.center().y), egui::Align2::LEFT_CENTER, label, font, ui.visuals().text_color());
        if !hint.is_empty() {
            response.on_hover_text(hint);
        }
        add(ui)
    })
    .inner
}

/// A text field (ImGui InputText): white, framed, left-aligned.
pub fn text_field(ui: &mut egui::Ui, text: &mut String, width: f32, hint: &str) -> egui::Response {
    ui.add_sized(egui::vec2(width, FRAME_HEIGHT), egui::TextEdit::singleline(text).hint_text(hint).align(egui::Align2::LEFT_CENTER).margin(egui::vec2(8.0, 4.0)))
}

/// A number typed in a field (ImGui InputScalar without step buttons): the
/// text is kept while the field has focus, and the value follows it when it parses.
pub fn number_field(ui: &mut egui::Ui, id_salt: impl std::hash::Hash + std::fmt::Debug, value: &mut f64, decimals: usize, width: f32) -> egui::Response {
    let id = ui.make_persistent_id(id_salt);
    let shown = format!("{:.*}", decimals, value);
    let mut text = ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| shown.clone());
    let response = ui.add_sized(egui::vec2(width, FRAME_HEIGHT), egui::TextEdit::singleline(&mut text).id(id.with("edit")).align(egui::Align2::LEFT_CENTER).margin(egui::vec2(8.0, 4.0)));
    let mut changed = false;
    if response.has_focus() {
        if response.changed() {
            if let Ok(v) = text.trim().parse::<f64>() {
                if v.is_finite() && v != *value {
                    *value = v;
                    changed = true;
                }
            }
        }
        ui.data_mut(|d| d.insert_temp(id, text));
    } else {
        ui.data_mut(|d| d.remove::<String>(id));
    }
    let mut response = response;
    if changed {
        response.mark_changed();
    }
    response
}

/// A float field with − and + buttons (ImGui InputFloat with a step).
pub fn float_input(ui: &mut egui::Ui, id_salt: impl std::hash::Hash + std::fmt::Debug, value: &mut f64, step: f64, decimals: usize, width: f32, range: std::ops::RangeInclusive<f64>) -> bool {
    let field_width = (width - 2.0 * (FRAME_HEIGHT + ui.spacing().item_spacing.x)).max(40.0);
    let mut changed = number_field(ui, id_salt, value, decimals, field_width).changed();
    if ui.add_sized(egui::vec2(FRAME_HEIGHT, FRAME_HEIGHT), egui::Button::new("-")).clicked() {
        *value -= step;
        changed = true;
    }
    if ui.add_sized(egui::vec2(FRAME_HEIGHT, FRAME_HEIGHT), egui::Button::new("+")).clicked() {
        *value += step;
        changed = true;
    }
    if changed {
        *value = value.clamp(*range.start(), *range.end());
    }
    changed
}


/// An integer field with − and + buttons (ImGui InputInt).
pub fn int_input(ui: &mut egui::Ui, id_salt: impl std::hash::Hash + std::fmt::Debug, value: &mut i32, width: f32, range: std::ops::RangeInclusive<i32>) -> bool {
    let mut v = *value as f64;
    let changed = float_input(ui, id_salt, &mut v, 1.0, 0, width, *range.start() as f64..=*range.end() as f64);
    if changed {
        *value = v.round() as i32;
    }
    changed
}

/// Integer fields without buttons side by side (ImGui InputInt3), filling `width`.
pub fn int_fields(ui: &mut egui::Ui, id_salt: impl std::hash::Hash + std::fmt::Debug + Copy, values: &mut [i32], width: f32, range: std::ops::RangeInclusive<i32>) -> bool {
    let n = values.len().max(1) as f32;
    let each = ((width - (n - 1.0) * ui.spacing().item_spacing.x) / n).max(30.0);
    let mut changed = false;
    for (k, v) in values.iter_mut().enumerate() {
        let mut f = *v as f64;
        if number_field(ui, (id_salt, k), &mut f, 0, each).changed() {
            *v = (f.round() as i32).clamp(*range.start(), *range.end());
            changed = true;
        }
    }
    changed
}

/// Float fields without buttons side by side (ImGui InputFloat3), filling `width`.
pub fn float_fields(ui: &mut egui::Ui, id_salt: impl std::hash::Hash + std::fmt::Debug + Copy, values: &mut [f64], decimals: usize, width: f32) -> bool {
    let n = values.len().max(1) as f32;
    let each = ((width - (n - 1.0) * ui.spacing().item_spacing.x) / n).max(30.0);
    let mut changed = false;
    for (k, v) in values.iter_mut().enumerate() {
        changed |= number_field(ui, (id_salt, k), v, decimals, each).changed();
    }
    changed
}

/// A bordered, rounded panel (ImGui child window with a border).
pub fn panel<R>(ui: &mut egui::Ui, size: egui::Vec2, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    let p = theme::palette(ui.visuals().dark_mode);
    egui::Frame::new()
        .fill(p.surface)
        .stroke(egui::Stroke::new(1.0, p.border))
        .corner_radius(8)
        .inner_margin(egui::Margin::same(10))
        .show(ui, |ui| {
            let inner = size - egui::vec2(22.0, 22.0);
            ui.set_min_size(inner);
            ui.set_max_width(inner.x);
            ui.vertical(|ui| add(ui)).inner
        })
        .inner
}

/// The minimum width of a footer button (5.5 font sizes, as the C++).
pub fn footer_button_width() -> f32 {
    (theme::BASE_FONT * 18.0 / 13.5 * 5.5).max(90.0)
}

/// The common footer: a rule, then [secondary] [primary] at the right edge
/// (dialogLayout::footer). Returns 1 for the primary action, 2 for the secondary, else 0.
pub fn footer(ui: &mut egui::Ui, primary: &str, secondary: Option<&str>, primary_enabled: bool) -> i32 {
    ui.add_space(2.0);
    ui.separator();
    action_bar(ui, "", Note::Info, primary, secondary, primary_enabled)
}

/// A footer row without the rule: the last action's outcome at the left
/// (one line, muted for Info), then [secondary] [primary] at the right (dialogLayout::actionBar).
pub fn action_bar(ui: &mut egui::Ui, status: &str, kind: Note, primary: &str, secondary: Option<&str>, primary_enabled: bool) -> i32 {
    let mut result = 0;
    let colour = if kind == Note::Info { theme::palette(ui.visuals().dark_mode).text_muted } else { note_colour(ui, kind) };
    ui.horizontal(|ui| {
        let w = footer_button_width();
        let buttons = if secondary.is_some() { 2.0 * w + ui.spacing().item_spacing.x } else { w };
        let room = (ui.available_width() - buttons - 2.0 * ui.spacing().item_spacing.x).max(0.0);
        ui.allocate_ui_with_layout(egui::vec2(room, FRAME_HEIGHT), egui::Layout::left_to_right(egui::Align::Center), |ui| {
            ui.set_min_width(room);
            if !status.is_empty() {
                ui.add(egui::Label::new(egui::RichText::new(status).color(colour)).truncate()).on_hover_text(status);
            }
        });
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui.add_enabled_ui(primary_enabled, |ui| theme::primary_button(ui, primary, Some(egui::vec2(w, FRAME_HEIGHT)))).inner.clicked() {
                result = 1;
            }
            if let Some(s) = secondary {
                if theme::secondary_button(ui, s, Some(egui::vec2(w, FRAME_HEIGHT))).clicked() {
                    result = 2;
                }
            }
        });
    });
    result
}

/// A table header row on the button colour (ImGui TableHeadersRow).
pub fn header_cell(ui: &mut egui::Ui, text: &str, width: f32) {
    let p = theme::palette(ui.visuals().dark_mode);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, FRAME_HEIGHT - 4.0), egui::Sense::hover());
    ui.painter().rect_filled(rect, 0.0, p.button);
    ui.painter().text(egui::pos2(rect.left() + 6.0, rect.center().y), egui::Align2::LEFT_CENTER, text, egui::TextStyle::Body.resolve(ui.style()), p.text);
}

// --- Fork 2 (pipeline, palette, file pickers, measurements): ImGui widgets. ---

/// themeStatusError.
pub fn status_error(ui: &egui::Ui) -> egui::Color32 {
    if ui.visuals().dark_mode {
        egui::Color32::from_rgb(255, 115, 107)
    } else {
        egui::Color32::from_rgb(191, 38, 38)
    }
}

/// A destructive action (delete, discard), its text in red (dialogLayout::dangerButton).
pub fn danger_button(ui: &mut egui::Ui, text: &str) -> egui::Response {
    let red = status_error(ui);
    ui.button(egui::RichText::new(text).color(red))
}

/// The primary action in the accent colour (dialogLayout::primaryButton), sized to its text.
pub fn accent_button(ui: &mut egui::Ui, text: &str, enabled: bool) -> egui::Response {
    let p = theme::palette(ui.visuals().dark_mode);
    let button = egui::Button::new(egui::RichText::new(text).color(egui::Color32::WHITE)).fill(p.accent).stroke(egui::Stroke::new(1.0, p.accent_strong));
    ui.add_enabled(enabled, button)
}

fn paint_check(ui: &egui::Ui, rect: egui::Rect, on: bool, hovered: bool, enabled: bool) {
    let p = theme::palette(ui.visuals().dark_mode);
    let fill = if hovered && enabled { p.frame_hover } else { p.frame };
    ui.painter().rect(rect, 6.0, fill, egui::Stroke::new(1.0, p.border), egui::StrokeKind::Inside);
    if on {
        let colour = if enabled { p.accent } else { p.text_muted };
        let r = rect.shrink(rect.width() * 0.25);
        let points = vec![
            egui::pos2(r.left(), r.center().y + r.height() * 0.04),
            egui::pos2(r.left() + r.width() * 0.36, r.bottom() - r.height() * 0.12),
            egui::pos2(r.right(), r.top() + r.height() * 0.12),
        ];
        ui.painter().add(egui::Shape::line(points, egui::Stroke::new(2.6, colour)));
    }
}

/// An ImGui check box without a label: a frame-sized square with an accent check mark.
pub fn check_box(ui: &mut egui::Ui, value: &mut bool) -> egui::Response {
    let sense = if ui.is_enabled() { egui::Sense::click() } else { egui::Sense::hover() };
    let (rect, mut response) = ui.allocate_exact_size(egui::vec2(FRAME_HEIGHT, FRAME_HEIGHT), sense);
    if response.clicked() {
        *value = !*value;
        response.mark_changed();
    }
    paint_check(ui, rect, *value, response.hovered(), ui.is_enabled());
    response
}

/// An ImGui check box with its label after it (the label toggles it too).
pub fn labelled_check(ui: &mut egui::Ui, value: &mut bool, label: &str) -> egui::Response {
    let font = egui::TextStyle::Body.resolve(ui.style());
    let colour = if ui.is_enabled() { ui.visuals().text_color() } else { ui.visuals().weak_text_color() };
    let galley = ui.painter().layout_no_wrap(label.to_string(), font, colour);
    let gap = 7.0;
    let sense = if ui.is_enabled() { egui::Sense::click() } else { egui::Sense::hover() };
    let (rect, mut response) = ui.allocate_exact_size(egui::vec2(FRAME_HEIGHT + gap + galley.size().x, FRAME_HEIGHT), sense);
    if response.clicked() {
        *value = !*value;
        response.mark_changed();
    }
    let square = egui::Rect::from_min_size(rect.min, egui::vec2(FRAME_HEIGHT, FRAME_HEIGHT));
    paint_check(ui, square, *value, response.hovered(), ui.is_enabled());
    ui.painter().galley(egui::pos2(square.right() + gap, rect.center().y - galley.size().y / 2.0), galley, colour);
    response
}

fn square_button(ui: &mut egui::Ui, enabled: bool) -> egui::Response {
    ui.add_enabled(enabled, egui::Button::new("").min_size(egui::vec2(FRAME_HEIGHT, FRAME_HEIGHT)))
}

/// ImGui::ArrowButton: a frame-sized square button with a filled triangle.
pub fn arrow_button(ui: &mut egui::Ui, up: bool, enabled: bool) -> egui::Response {
    let response = square_button(ui, enabled);
    let c = response.rect.center();
    let h = FRAME_HEIGHT * 0.2;
    let colour = if enabled { ui.visuals().text_color() } else { ui.visuals().weak_text_color() };
    let points = if up {
        vec![egui::pos2(c.x, c.y - h), egui::pos2(c.x + h, c.y + h * 0.85), egui::pos2(c.x - h, c.y + h * 0.85)]
    } else {
        vec![egui::pos2(c.x - h, c.y - h * 0.85), egui::pos2(c.x + h, c.y - h * 0.85), egui::pos2(c.x, c.y + h)]
    };
    ui.painter().add(egui::Shape::convex_polygon(points, colour, egui::Stroke::NONE));
    response
}

/// A frame-sized square button with a drawn cross (the same in every font).
pub fn cross_button(ui: &mut egui::Ui, enabled: bool) -> egui::Response {
    let response = square_button(ui, enabled);
    let c = response.rect.center();
    let r = FRAME_HEIGHT * 0.2;
    let stroke = egui::Stroke::new(1.5, if enabled { ui.visuals().text_color() } else { ui.visuals().weak_text_color() });
    ui.painter().line_segment([c - egui::vec2(r, r), c + egui::vec2(r, r)], stroke);
    ui.painter().line_segment([c + egui::vec2(-r, r), c + egui::vec2(r, -r)], stroke);
    response
}

/// ImGui::SeparatorText: the text, then a rule to the right edge.
pub fn separator_text(ui: &mut egui::Ui, text: &str) {
    let font = egui::TextStyle::Body.resolve(ui.style());
    let colour = ui.visuals().text_color();
    let galley = ui.painter().layout_no_wrap(text.to_string(), font, colour);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), galley.size().y + 8.0), egui::Sense::hover());
    let y = rect.center().y;
    let start = rect.left() + galley.size().x + 8.0;
    ui.painter().galley(egui::pos2(rect.left(), y - galley.size().y / 2.0), galley, colour);
    if start < rect.right() {
        ui.painter().hline(start..=rect.right(), y, egui::Stroke::new(1.0, theme::palette(ui.visuals().dark_mode).separator));
    }
}

/// A bordered child panel on the surface colour with the window padding (ImGui child with borders).
pub fn child_frame<R>(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui) -> R) -> egui::InnerResponse<R> {
    let p = theme::palette(ui.visuals().dark_mode);
    egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(16, 14)).show(ui, |ui| {
        ui.set_width(ui.available_width());
        add(ui)
    })
}

/// A form row as the C++ two-column tables: the label in the left `fraction`
/// of the width, the field filling the rest.
pub fn form_row<R>(ui: &mut egui::Ui, label: &str, fraction: f32, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    let width = ui.available_width();
    ui.horizontal(|ui| {
        let column = (width * fraction - ui.spacing().item_spacing.x).max(40.0);
        let font = egui::TextStyle::Body.resolve(ui.style());
        let galley = ui.painter().layout(label.to_string(), font, ui.visuals().text_color(), column);
        let (rect, _) = ui.allocate_exact_size(egui::vec2(column, FRAME_HEIGHT.max(galley.size().y)), egui::Sense::hover());
        ui.painter().galley(egui::pos2(rect.left(), rect.center().y - galley.size().y / 2.0), galley, ui.visuals().text_color());
        ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), FRAME_HEIGHT), egui::Layout::left_to_right(egui::Align::Center), |ui| {
            ui.spacing_mut().text_edit_width = ui.available_width();
            add(ui)
        })
        .inner
    })
    .inner
}

// ---------------------------------------------------------------------------
// Widgets of the analysis dialogs and the scientific-tool window (fork 3):
// ImGui SeparatorText, SliderInt/SliderFloat, CollapsingHeader, tables with
// headers and borders, and centred empty states.
// ---------------------------------------------------------------------------

/// An ImGui slider: a framed field with a grab marker at the value's place
/// and the value centred in it; dragging or clicking sets the value.
fn slider_frame(ui: &mut egui::Ui, width: f32, t: f64, text: &str) -> (egui::Response, egui::Rect) {
    let p = theme::palette(ui.visuals().dark_mode);
    let (rect, response) = ui.allocate_exact_size(egui::vec2(width, FRAME_HEIGHT), egui::Sense::click_and_drag());
    let fill = if response.dragged() { p.frame_active } else if response.hovered() { p.frame_hover } else { p.frame };
    ui.painter().rect(rect, 6.0, fill, egui::Stroke::new(1.0, p.border), egui::StrokeKind::Inside);
    let grab_w = 12.0f32;
    let inner = rect.shrink(2.0);
    let x = inner.left() + (inner.width() - grab_w) * t.clamp(0.0, 1.0) as f32;
    let grab = egui::Rect::from_min_size(egui::pos2(x, inner.top() + 1.0), egui::vec2(grab_w, inner.height() - 2.0));
    let a = p.accent;
    let alpha = if response.dragged() { 255 } else { 166 };
    ui.painter().rect_filled(grab, 5.0, egui::Color32::from_rgba_unmultiplied(a.r(), a.g(), a.b(), alpha));
    let font = egui::TextStyle::Body.resolve(ui.style());
    let colour = if ui.is_enabled() { p.text } else { p.text_muted };
    ui.painter().text(rect.center(), egui::Align2::CENTER_CENTER, text, font, colour);
    (response, inner)
}

fn slider_fraction(response: &egui::Response, inner: egui::Rect) -> Option<f64> {
    if response.dragged() || response.clicked() || response.drag_started() {
        let pos = response.interact_pointer_pos()?;
        let grab_w = 12.0f32;
        let t = ((pos.x - inner.left() - grab_w / 2.0) / (inner.width() - grab_w).max(1.0)).clamp(0.0, 1.0);
        return Some(t as f64);
    }
    None
}

/// ImGui::SliderInt with a format such as "%d" or "%d pass".
pub fn slider_i32(ui: &mut egui::Ui, value: &mut i32, range: std::ops::RangeInclusive<i32>, width: f32, text: impl Fn(i32) -> String) -> egui::Response {
    let (lo, hi) = (*range.start(), *range.end());
    let t = if hi > lo { (*value - lo) as f64 / (hi - lo) as f64 } else { 0.0 };
    let (mut response, inner) = slider_frame(ui, width, t, &text(*value));
    if let Some(f) = slider_fraction(&response, inner) {
        let v = (lo as f64 + f * (hi - lo) as f64).round() as i32;
        if v != *value {
            *value = v;
            response.mark_changed();
        }
    }
    response
}

/// ImGui::SliderFloat with `decimals` digits (and an optional suffix in `text`).
pub fn slider_f64(ui: &mut egui::Ui, value: &mut f64, range: std::ops::RangeInclusive<f64>, width: f32, text: impl Fn(f64) -> String) -> egui::Response {
    let (lo, hi) = (*range.start(), *range.end());
    let t = if hi > lo { (*value - lo) / (hi - lo) } else { 0.0 };
    let (mut response, inner) = slider_frame(ui, width, t, &text(*value));
    if let Some(f) = slider_fraction(&response, inner) {
        let v = lo + f * (hi - lo);
        if v != *value {
            *value = v;
            response.mark_changed();
        }
    }
    response
}

pub fn slider_f32(ui: &mut egui::Ui, value: &mut f32, range: std::ops::RangeInclusive<f32>, width: f32, text: impl Fn(f32) -> String) -> egui::Response {
    let mut v = *value as f64;
    let response = slider_f64(ui, &mut v, *range.start() as f64..=*range.end() as f64, width, |x| text(x as f32));
    if response.changed() {
        *value = v as f32;
    }
    response
}

/// ImGui::CollapsingHeader: a full-width bar on the header colour with an
/// arrow; the body below it when open.
pub fn collapsing<R>(ui: &mut egui::Ui, id_salt: impl std::hash::Hash + std::fmt::Debug, title: &str, default_open: bool, body: impl FnOnce(&mut egui::Ui) -> R) -> Option<R> {
    let p = theme::palette(ui.visuals().dark_mode);
    let id = ui.make_persistent_id(id_salt);
    let mut state = egui::collapsing_header::CollapsingState::load_with_default_open(ui.ctx(), id, default_open);
    let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width(), FRAME_HEIGHT - 2.0), egui::Sense::click());
    if response.clicked() {
        state.toggle(ui);
    }
    let a = p.accent;
    let alpha = match (response.is_pointer_button_down_on(), response.hovered(), ui.visuals().dark_mode) {
        (true, _, false) => 82,
        (true, _, true) => 117,
        (false, true, false) => 56,
        (false, true, true) => 87,
        (false, false, false) => 36,
        (false, false, true) => 61,
    };
    ui.painter().rect_filled(rect, 6.0, egui::Color32::from_rgba_unmultiplied(a.r(), a.g(), a.b(), alpha));
    // The arrow: down when open, right when closed.
    let c = egui::pos2(rect.left() + 14.0, rect.center().y);
    let s = 5.0;
    let points = if state.is_open() { vec![c + egui::vec2(-s, -s * 0.6), c + egui::vec2(s, -s * 0.6), c + egui::vec2(0.0, s * 0.7)] } else { vec![c + egui::vec2(-s * 0.6, -s), c + egui::vec2(-s * 0.6, s), c + egui::vec2(s * 0.7, 0.0)] };
    ui.painter().add(egui::Shape::convex_polygon(points, p.text, egui::Stroke::NONE));
    ui.painter().text(egui::pos2(rect.left() + 30.0, rect.center().y), egui::Align2::LEFT_CENTER, title, egui::TextStyle::Body.resolve(ui.style()), p.text);
    state.store(ui.ctx());
    if state.is_open() {
        Some(body(ui))
    } else {
        None
    }
}

/// A table (ImGui BeginTable with TableHeadersRow, BordersOuter/Inner and
/// RowBg). `widths`: above 1 a fixed width in points, else a fraction of the
/// width; the last column takes the rest. Blank headers draw no header row
/// (for the rows of a scrolled table under its header). `cell` draws a cell.
pub fn table(ui: &mut egui::Ui, headers: &[&str], widths: &[f32], rows: usize, mut cell: impl FnMut(&mut egui::Ui, usize, usize)) {
    let p = theme::palette(ui.visuals().dark_mode);
    let total = ui.available_width();
    let row_h = FRAME_HEIGHT - 6.0;
    let n = headers.len();
    let header = headers.iter().any(|h| !h.is_empty());
    let mut xs = Vec::with_capacity(n + 1);
    let mut x = 0.0f32;
    xs.push(0.0);
    for k in 0..n {
        let w = if k + 1 == n {
            (total - x).max(0.0)
        } else {
            let w = widths.get(k).copied().unwrap_or(1.0 / n as f32);
            if w > 1.0 { w } else { w * total }
        };
        x += w;
        xs.push(x.min(total));
    }
    let header_rows = if header { 1 } else { 0 };
    let height = row_h * (rows + header_rows) as f32;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(total, height), egui::Sense::hover());
    let painter = ui.painter().clone();
    let font = egui::TextStyle::Body.resolve(ui.style());
    if header {
        painter.rect_filled(egui::Rect::from_min_size(rect.min, egui::vec2(total, row_h)), 0.0, p.button);
        for (k, h) in headers.iter().enumerate() {
            painter.text(egui::pos2(rect.left() + xs[k] + 6.0, rect.top() + row_h / 2.0), egui::Align2::LEFT_CENTER, *h, font.clone(), p.text);
        }
    }
    let alt = if ui.visuals().dark_mode { egui::Color32::from_rgba_unmultiplied(255, 255, 255, 8) } else { egui::Color32::from_rgba_unmultiplied(0, 0, 0, 6) };
    for r in 0..rows {
        let top = rect.top() + row_h * (r + header_rows) as f32;
        if r % 2 == 1 {
            painter.rect_filled(egui::Rect::from_min_size(egui::pos2(rect.left(), top), egui::vec2(total, row_h)), 0.0, alt);
        }
        for k in 0..n {
            if xs[k + 1] - xs[k] < 2.0 {
                continue;
            }
            let cell_rect = egui::Rect::from_min_max(egui::pos2(rect.left() + xs[k] + 6.0, top), egui::pos2(rect.left() + xs[k + 1] - 4.0, top + row_h));
            let mut child = ui.new_child(egui::UiBuilder::new().max_rect(cell_rect).layout(egui::Layout::left_to_right(egui::Align::Center)));
            child.spacing_mut().item_spacing.x = 4.0;
            cell(&mut child, r, k);
        }
        if r > 0 || header {
            painter.hline(rect.x_range(), top, egui::Stroke::new(1.0, p.separator));
        }
    }
    // Borders: outer, and between the columns.
    painter.rect_stroke(rect, 0.0, egui::Stroke::new(1.0, p.border), egui::StrokeKind::Inside);
    for &x in &xs[1..n] {
        if x < total - 1.0 {
            painter.vline(rect.left() + x, rect.y_range(), egui::Stroke::new(1.0, p.separator));
        }
    }
}

/// The width of ImGui's default items: 65 % of the window's content width.
pub fn item_width(ui: &egui::Ui) -> f32 {
    (ui.available_width() * 0.65).max(80.0)
}

/// dialogLayout::emptyState: muted text centred in the space left.
pub fn empty_state(ui: &mut egui::Ui, text: &str) {
    let p = theme::palette(ui.visuals().dark_mode);
    let rect = ui.available_rect_before_wrap();
    let font = egui::TextStyle::Body.resolve(ui.style());
    let galley = ui.painter().layout(text.to_string(), font, p.text_muted, (rect.width() * 0.8).max(1.0));
    let size = galley.size();
    let pos = egui::pos2(rect.left() + ((rect.width() - size.x) * 0.5).max(0.0), rect.top() + ((rect.height() - size.y) * 0.45).max(0.0));
    ui.painter().galley(pos, galley, p.text_muted);
    ui.allocate_rect(rect, egui::Sense::hover());
}

/// ImGui::RadioButton: a frame-height circle on the frame colour, an accent
/// dot inside when selected, the label after it (clicking either selects it).
pub fn radio(ui: &mut egui::Ui, selected: bool, label: &str) -> egui::Response {
    let p = theme::palette(ui.visuals().dark_mode);
    let font = egui::TextStyle::Body.resolve(ui.style());
    let colour = if ui.is_enabled() { p.text } else { p.text_muted };
    let galley = ui.painter().layout_no_wrap(label.to_string(), font, colour);
    let size = FRAME_HEIGHT - 6.0;
    let width = size + ui.spacing().item_spacing.x + galley.size().x;
    let (rect, response) = ui.allocate_exact_size(egui::vec2(width, FRAME_HEIGHT), egui::Sense::click());
    let centre = egui::pos2(rect.left() + size / 2.0, rect.center().y);
    let fill = if response.hovered() && ui.is_enabled() { p.frame_hover } else { p.frame };
    ui.painter().circle(centre, size / 2.0, fill, egui::Stroke::new(1.0, p.border));
    if selected {
        ui.painter().circle_filled(centre, size / 2.0 - 4.0, if ui.is_enabled() { p.accent } else { p.text_muted });
    }
    ui.painter().galley(egui::pos2(rect.left() + size + ui.spacing().item_spacing.x, rect.center().y - galley.size().y / 2.0), galley, colour);
    response
}

/// The primary button across the available width, its label centred
/// (dialogLayout::primaryButton with ImVec2(-1, height)).
pub fn wide_primary(ui: &mut egui::Ui, text: &str, height: f32, enabled: bool) -> egui::Response {
    let width = ui.available_width();
    ui.allocate_ui_with_layout(egui::vec2(width, height), egui::Layout::top_down_justified(egui::Align::Center), |ui| ui.add_enabled_ui(enabled, |ui| theme::primary_button(ui, text, Some(egui::vec2(width, height)))).inner).inner
}

/// An ImGui tab bar: tabs on the button colour, the selected one on the
/// surface colour with an accent line along its top, over a rule. True when
/// the selection changed.
pub fn tab_bar(ui: &mut egui::Ui, selected: &mut usize, names: &[&str]) -> bool {
    let p = theme::palette(ui.visuals().dark_mode);
    let font = egui::TextStyle::Body.resolve(ui.style());
    let mut changed = false;
    let row = ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 4.0;
        for (i, name) in names.iter().enumerate() {
            let galley = ui.painter().layout_no_wrap(name.to_string(), font.clone(), p.text);
            let size = egui::vec2(galley.size().x + 20.0, FRAME_HEIGHT - 2.0);
            let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click());
            let on = *selected == i;
            let fill = if on { p.surface } else if response.hovered() { p.button_hover } else { p.button };
            ui.painter().rect_filled(rect, egui::CornerRadius { nw: 6, ne: 6, sw: 0, se: 0 }, fill);
            if on {
                ui.painter().line_segment([rect.left_top() + egui::vec2(3.0, 1.0), rect.right_top() + egui::vec2(-3.0, 1.0)], egui::Stroke::new(2.0, p.accent));
            }
            ui.painter().galley(egui::pos2(rect.left() + 10.0, rect.center().y - galley.size().y / 2.0), galley, p.text);
            if response.clicked() && !on {
                *selected = i;
                changed = true;
            }
        }
    });
    let r = row.response.rect;
    ui.painter().hline(ui.max_rect().x_range(), r.bottom(), egui::Stroke::new(1.0, p.surface));
    changed
}

/// ImGui::SmallButton: a compact button (still padded, so its label does not touch the frame).
pub fn small_button(ui: &mut egui::Ui, text: &str, enabled: bool) -> egui::Response {
    ui.scope(|ui| {
        ui.spacing_mut().button_padding = egui::vec2(8.0, 3.0);
        ui.add_enabled(enabled, egui::Button::new(text))
    })
    .inner
}

/// ImGui::ColorEdit3: red, green and blue (0-255) in framed drag fields
/// with coloured edges, a swatch that opens a colour picker, and the label
/// after them. `width` covers the fields and the swatch.
pub fn colour_edit3(ui: &mut egui::Ui, id_salt: impl std::hash::Hash + std::fmt::Debug, colour: &mut [f32; 3], width: f32, label: &str) -> bool {
    let p = theme::palette(ui.visuals().dark_mode);
    let mut changed = false;
    ui.push_id(id_salt, |ui| {
        ui.horizontal(|ui| {
            let spacing = ui.spacing().item_spacing.x;
            let each = ((width - FRAME_HEIGHT - 3.0 * spacing) / 3.0).max(36.0);
            let edges = [egui::Color32::from_rgb(220, 60, 60), egui::Color32::from_rgb(60, 180, 70), egui::Color32::from_rgb(70, 110, 230)];
            for k in 0..3 {
                let mut v = (colour[k].clamp(0.0, 1.0) * 255.0).round() as i32;
                let r = ui.scope(|ui| {
                    let w = &mut ui.visuals_mut().widgets;
                    w.inactive.weak_bg_fill = p.frame;
                    w.hovered.weak_bg_fill = p.frame_hover;
                    w.active.weak_bg_fill = p.frame_active;
                    ui.add_sized([each, FRAME_HEIGHT], egui::DragValue::new(&mut v).range(0..=255).speed(1.0))
                });
                let rect = r.inner.rect;
                ui.painter().line_segment([rect.left_top() + egui::vec2(1.5, 3.0), rect.left_bottom() + egui::vec2(1.5, -3.0)], egui::Stroke::new(3.0, edges[k]));
                if r.inner.changed() {
                    colour[k] = v as f32 / 255.0;
                    changed = true;
                }
            }
            // The colours are sRGB values (not linear), as the C++ keeps them.
            let mut srgb = [0u8; 3];
            for k in 0..3 {
                srgb[k] = (colour[k].clamp(0.0, 1.0) * 255.0).round() as u8;
            }
            ui.scope(|ui| {
                ui.spacing_mut().interact_size = egui::vec2(FRAME_HEIGHT, FRAME_HEIGHT);
                if egui::color_picker::color_edit_button_srgb(ui, &mut srgb).changed() {
                    *colour = [srgb[0] as f32 / 255.0, srgb[1] as f32 / 255.0, srgb[2] as f32 / 255.0];
                    changed = true;
                }
            });
            if !label.is_empty() {
                ui.label(label);
            }
        });
    });
    changed
}

/// The height of a dialog body that fills a resizable window: the room left
/// minus what the dialog drew below the body last frame (measured by
/// `FitHeight::finish`), so the window keeps its size instead of growing
/// frame by frame when a guessed footer height is too small.
pub struct FitHeight {
    id: egui::Id,
    top: f32,
    pub height: f32,
}

pub fn fit_height(ui: &egui::Ui, salt: &str, guess: f32, min: f32) -> FitHeight {
    let id = ui.scope_id().with(("fit height", salt));
    let below = ui.ctx().data(|d| d.get_temp::<f32>(id)).unwrap_or(guess);
    FitHeight { id, top: ui.cursor().top(), height: (ui.available_height() - below).max(min) }
}

impl FitHeight {
    /// Call after the last widget of the window.
    pub fn finish(&self, ui: &egui::Ui) {
        let below = (ui.min_rect().bottom() - self.top - self.height).max(0.0);
        ui.ctx().data_mut(|d| d.insert_temp(self.id, below));
    }
}
