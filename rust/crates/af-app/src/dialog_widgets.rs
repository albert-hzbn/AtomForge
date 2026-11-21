//! ImGui-shaped controls for the dialogs ported from the C++ (the tool,
//! overlay and electronic windows): framed fields with the value centred and
//! the label at their right (DragInt, DragFloat), sliders with a grab block
//! (SliderFloat, SliderInt), frame-height check boxes, colour swatches
//! (ColorEdit3 without inputs), callouts and right-aligned footers
//! (ui/DialogLayout.h).

use eframe::egui::{self, Color32, CornerRadius, Stroke};

use crate::theme::{self, FRAME_HEIGHT};

fn palette(ui: &egui::Ui) -> theme::Palette {
    theme::palette(ui.visuals().dark_mode)
}

fn label_after(ui: &mut egui::Ui, label: &str) {
    if !label.is_empty() {
        ui.label(label);
    }
}

/// ImGui::DragInt `width` wide, the label at its right.
pub fn drag_int(ui: &mut egui::Ui, value: &mut i32, width: f32, speed: f64, range: std::ops::RangeInclusive<i32>, label: &str) -> egui::Response {
    ui.horizontal(|ui| {
        field_style(ui);
        let r = ui.add_sized([width, FRAME_HEIGHT], egui::DragValue::new(value).speed(speed).range(range));
        label_after(ui, label);
        r
    })
    .inner
}

/// ImGui::DragFloat with a "%.Nf" format.
pub fn drag_f32(ui: &mut egui::Ui, value: &mut f32, width: f32, speed: f64, range: std::ops::RangeInclusive<f32>, decimals: usize, label: &str) -> egui::Response {
    ui.horizontal(|ui| {
        field_style(ui);
        let r = ui.add_sized([width, FRAME_HEIGHT], egui::DragValue::new(value).speed(speed).range(range).fixed_decimals(decimals));
        label_after(ui, label);
        r
    })
    .inner
}

/// ImGui::SliderFloat / SliderInt: a framed track with the grab as an
/// accent block and the value centred over it.
pub fn slider(ui: &mut egui::Ui, value: &mut f64, range: std::ops::RangeInclusive<f64>, width: f32, text: &str, label: &str, integer: bool) -> egui::Response {
    ui.horizontal(|ui| {
        let p = palette(ui);
        let (rect, mut response) = ui.allocate_exact_size(egui::vec2(width, FRAME_HEIGHT), egui::Sense::click_and_drag());
        let (lo, hi) = (*range.start(), *range.end());
        let track = rect.shrink(2.0);
        let grab_w = if integer { ((track.width()) / ((hi - lo).abs() + 1.0) as f32).clamp(8.0, track.width() * 0.5) } else { 8.0 };
        let usable = (track.width() - grab_w).max(1.0);
        if ui.is_enabled() && (response.dragged() || response.clicked()) {
            if let Some(pos) = response.interact_pointer_pos() {
                let t = ((pos.x - track.left() - grab_w * 0.5) / usable).clamp(0.0, 1.0) as f64;
                let mut v = lo + t * (hi - lo);
                if integer {
                    v = v.round();
                }
                if v != *value {
                    *value = v;
                    response.mark_changed();
                }
            }
        }
        let t = if hi > lo { ((*value - lo) / (hi - lo)).clamp(0.0, 1.0) as f32 } else { 0.0 };
        let fill = if response.hovered() || response.dragged() { p.frame_hover } else { p.frame };
        let fade = |c: Color32| if ui.is_enabled() { c } else { c.gamma_multiply(0.5) };
        let painter = ui.painter();
        painter.rect(rect, CornerRadius::same(6), fade(fill), Stroke::new(1.0, fade(p.border)), egui::StrokeKind::Inside);
        let grab = egui::Rect::from_min_size(egui::pos2(track.left() + t * usable, track.top()), egui::vec2(grab_w, track.height()));
        let grab_colour = if response.dragged() { p.accent_strong } else { Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), 166) };
        painter.rect_filled(grab, CornerRadius::same(5), fade(grab_colour));
        painter.text(rect.center(), egui::Align2::CENTER_CENTER, text, egui::TextStyle::Body.resolve(ui.style()), fade(p.text));
        label_after(ui, label);
        response
    })
    .inner
}

pub fn slider_f32(ui: &mut egui::Ui, value: &mut f32, range: std::ops::RangeInclusive<f32>, width: f32, decimals: usize, label: &str) -> egui::Response {
    let mut v = *value as f64;
    let r = slider(ui, &mut v, (*range.start() as f64)..=(*range.end() as f64), width, &format!("{:.*}", decimals, *value), label, false);
    *value = v as f32;
    r
}

pub fn slider_i32(ui: &mut egui::Ui, value: &mut i32, range: std::ops::RangeInclusive<i32>, width: f32, label: &str) -> egui::Response {
    let mut v = *value as f64;
    let r = slider(ui, &mut v, (*range.start() as f64)..=(*range.end() as f64), width, &value.to_string(), label, true);
    *value = v as i32;
    r
}

pub fn slider_usize(ui: &mut egui::Ui, value: &mut usize, range: std::ops::RangeInclusive<usize>, width: f32, label: &str) -> egui::Response {
    let mut v = *value as f64;
    let r = slider(ui, &mut v, (*range.start() as f64)..=(*range.end() as f64), width, &value.to_string(), label, true);
    *value = v.max(0.0) as usize;
    r
}

/// ImGui::Checkbox: a frame-height square with the check mark in the accent colour.
pub fn checkbox(ui: &mut egui::Ui, value: &mut bool, label: impl AsRef<str>) -> egui::Response {
    let label = label.as_ref();
    let p = palette(ui);
    let font = egui::TextStyle::Body.resolve(ui.style());
    let side = FRAME_HEIGHT - 4.0;
    // A label longer than the room left wraps instead of widening the window.
    let wrap = (ui.available_width() - side - 7.0).max(60.0);
    let galley = ui.painter().layout(label.to_string(), font, p.text, wrap);
    let width = side + if label.is_empty() { 0.0 } else { 7.0 + galley.size().x };
    let (rect, mut response) = ui.allocate_exact_size(egui::vec2(width, FRAME_HEIGHT.max(galley.size().y + 4.0)), egui::Sense::click());
    if response.clicked() {
        *value = !*value;
        response.mark_changed();
    }
    let enabled = ui.is_enabled();
    let fade = |c: Color32| if enabled { c } else { c.gamma_multiply(0.5) };
    let b = egui::Rect::from_min_size(egui::pos2(rect.left(), rect.top() + (FRAME_HEIGHT - side) / 2.0), egui::vec2(side, side));
    let fill = if response.hovered() { p.frame_hover } else { p.frame };
    ui.painter().rect(b, CornerRadius::same(6), fade(fill), Stroke::new(1.0, fade(p.border)), egui::StrokeKind::Inside);
    if *value {
        let c = b.center();
        let s = side * 0.26;
        ui.painter().line(vec![c + egui::vec2(-s * 1.1, 0.0), c + egui::vec2(-s * 0.3, s * 0.8), c + egui::vec2(s * 1.2, -s * 0.9)], Stroke::new(2.2, fade(p.accent)));
    }
    if !label.is_empty() {
        ui.painter().galley(egui::pos2(b.right() + 7.0, rect.center().y - galley.size().y / 2.0), galley, fade(p.text));
    }
    response
}

/// ImGui::ColorEdit3 without inputs: a frame-height swatch (a picker on
/// click), the label at its right.
pub fn colour_edit(ui: &mut egui::Ui, colour: &mut [f32; 3], label: &str) -> egui::Response {
    ui.horizontal(|ui| {
        ui.spacing_mut().interact_size = egui::vec2(FRAME_HEIGHT - 4.0, FRAME_HEIGHT - 4.0);
        let mut rgb = *colour;
        let r = ui.color_edit_button_rgb(&mut rgb);
        if r.changed() {
            *colour = rgb;
        }
        label_after(ui, label);
        r
    })
    .inner
}

/// A single-line text field `width` wide with a hint and the label at its right.
pub fn input_text(ui: &mut egui::Ui, text: &mut String, width: f32, hint: &str, label: &str) -> egui::Response {
    ui.horizontal(|ui| {
        let r = ui.add_sized([width, FRAME_HEIGHT], egui::TextEdit::singleline(text).hint_text(hint).align(egui::Align2::LEFT_CENTER));
        label_after(ui, label);
        r
    })
    .inner
}

/// Muted text (ImGui::TextDisabled).
pub fn disabled(ui: &mut egui::Ui, text: impl Into<String>) -> egui::Response {
    let c = palette(ui).text_muted;
    ui.label(egui::RichText::new(text.into()).color(c))
}

/// ImGui::SeparatorText: the text, then a rule to the right edge on its line.
pub fn separator_text(ui: &mut egui::Ui, text: &str) {
    ui.add_space(2.0);
    let p = palette(ui);
    let response = ui.label(text);
    let y = response.rect.center().y;
    let x0 = response.rect.right() + 8.0;
    let x1 = ui.max_rect().right();
    if x1 > x0 {
        ui.painter().hline(x0..=x1, y, Stroke::new(1.0, p.separator));
    }
}

/// ImGui::Separator across the available width.
pub fn rule(ui: &mut egui::Ui) {
    let p = palette(ui);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), 1.0), egui::Sense::hover());
    ui.painter().hline(rect.x_range(), rect.center().y, Stroke::new(1.0, p.separator));
}

/// A section heading (dialogLayout::section): semibold title over a thin rule.
pub fn section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(4.0);
    ui.label(egui::RichText::new(title).font(theme::heading(1.0)));
    rule(ui);
    ui.add_space(2.0);
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Note {
    Info,
    Success,
    Warning,
    Error,
}

pub fn note_colour(ui: &egui::Ui, kind: Note) -> Color32 {
    let dark = ui.visuals().dark_mode;
    match kind {
        Note::Info => palette(ui).accent,
        Note::Success => {
            if dark {
                Color32::from_rgb(110, 200, 130)
            } else {
                Color32::from_rgb(30, 130, 60)
            }
        }
        Note::Warning => ui.visuals().warn_fg_color,
        Note::Error => ui.visuals().error_fg_color,
    }
}

/// A callout (dialogLayout::note): tinted rounded box, coloured left edge,
/// icon and wrapped text.
pub fn note(ui: &mut egui::Ui, kind: Note, text: &str) {
    let colour = note_colour(ui, kind);
    let p = palette(ui);
    let width = ui.available_width();
    let pad = 9.0;
    let icon_w = 22.0;
    let font = egui::TextStyle::Body.resolve(ui.style());
    let galley = ui.painter().layout(text.to_string(), font.clone(), p.text, (width - 2.0 * pad - icon_w).max(40.0));
    let height = galley.size().y + 2.0 * pad;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, height), egui::Sense::hover());
    let tint = Color32::from_rgba_unmultiplied(colour.r(), colour.g(), colour.b(), if ui.visuals().dark_mode { 41 } else { 23 });
    ui.painter().rect_filled(rect, CornerRadius::same(6), tint);
    ui.painter().rect_filled(egui::Rect::from_min_size(rect.min, egui::vec2(3.0, height)), CornerRadius { nw: 6, sw: 6, ne: 0, se: 0 }, colour);
    let glyph = match kind {
        Note::Warning => theme::icon::WARNING,
        Note::Error => theme::icon::ERROR,
        Note::Success => theme::icon::CHECK,
        Note::Info => theme::icon::INFO,
    };
    let row = galley.rows.first().map(|r| r.height()).unwrap_or(16.0);
    ui.painter().text(egui::pos2(rect.left() + pad + icon_w * 0.4, rect.top() + pad + row * 0.5), egui::Align2::CENTER_CENTER, glyph, font, colour);
    ui.painter().galley(egui::pos2(rect.left() + pad + icon_w, rect.top() + pad), galley, p.text);
}

/// A destructive action (dialogLayout::dangerButton): the label in red.
pub fn danger_button(ui: &mut egui::Ui, label: &str) -> egui::Response {
    let red = note_colour(ui, Note::Error);
    ui.add(egui::Button::new(egui::RichText::new(label).color(red)))
}

/// The primary action in the accent colour (dialogLayout::primaryButton).
pub fn primary_button(ui: &mut egui::Ui, label: &str, size: egui::Vec2) -> egui::Response {
    let p = palette(ui);
    let border = Color32::from_rgb((p.accent.r() as f32 * 0.9) as u8, (p.accent.g() as f32 * 0.9) as u8, (p.accent.b() as f32 * 0.9) as u8);
    let button = egui::Button::new(egui::RichText::new(label).color(Color32::WHITE)).fill(p.accent).stroke(Stroke::new(1.0, border));
    // Sized buttons centre their label, as ImGui buttons do.
    if size.x > 0.0 { ui.add_sized([size.x, size.y.max(FRAME_HEIGHT)], button) } else { ui.add(button) }
}

fn footer_width(ui: &egui::Ui, label: &str) -> f32 {
    let font = egui::TextStyle::Button.resolve(ui.style());
    let w = ui.painter().layout_no_wrap(label.to_string(), font, Color32::WHITE).size().x + 2.0 * ui.spacing().button_padding.x;
    w.max(18.0 * 5.5)
}

/// The common footer (dialogLayout::footer): a rule, then [secondary]
/// [primary] at the right edge. 1: primary pressed, 2: secondary, else 0.
pub fn footer(ui: &mut egui::Ui, primary: &str, secondary: Option<&str>, primary_enabled: bool) -> i32 {
    ui.add_space(2.0);
    rule(ui);
    let mut result = 0;
    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
        let pw = footer_width(ui, primary);
        if ui.add_enabled_ui(primary_enabled, |ui| primary_button(ui, primary, egui::vec2(pw, FRAME_HEIGHT))).inner.clicked() {
            result = 1;
        }
        if let Some(s) = secondary {
            let sw = footer_width(ui, s);
            if ui.add(egui::Button::new(s).min_size(egui::vec2(sw, FRAME_HEIGHT))).clicked() {
                result = 2;
            }
        }
    });
    result
}

/// A placeholder centred in a results area (dialogLayout::emptyState).
pub fn empty_state(ui: &mut egui::Ui, text: &str, height: f32) {
    let c = palette(ui).text_muted;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), height), egui::Sense::hover());
    let font = egui::TextStyle::Body.resolve(ui.style());
    let galley = ui.painter().layout(text.to_string(), font, c, rect.width() * 0.8);
    let pos = egui::pos2(rect.center().x - galley.size().x / 2.0, rect.top() + (rect.height() - galley.size().y) * 0.45);
    ui.painter().galley(pos, galley, c);
}

/// A bordered child region (ImGui::BeginChild with a border).
pub fn child<R>(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    let p = palette(ui);
    egui::Frame::new()
        .stroke(Stroke::new(1.0, p.border))
        .corner_radius(8)
        .inner_margin(egui::Margin::same(8))
        .show(ui, |ui| {
            ui.set_width(ui.available_width());
            add(ui)
        })
        .inner
}

/// A full-width button (responsive::button with the available width).
pub fn wide_button(ui: &mut egui::Ui, label: &str) -> egui::Response {
    ui.add_sized([ui.available_width(), FRAME_HEIGHT], egui::Button::new(label))
}

/// ImGui::InputInt: the field, then − and + buttons, then the label; `width`
/// covers the field and the buttons.
pub fn input_int(ui: &mut egui::Ui, value: &mut i64, width: f32, range: std::ops::RangeInclusive<i64>, label: &str) -> egui::Response {
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 4.0;
        let field = (width - 2.0 * (FRAME_HEIGHT + 4.0)).max(30.0);
        field_style(ui);
        let mut r = ui.add_sized([field, FRAME_HEIGHT], egui::DragValue::new(value).speed(0.0).range(range.clone()));
        let before = *value;
        if ui.add_sized([FRAME_HEIGHT, FRAME_HEIGHT], egui::Button::new("-")).clicked() {
            *value = (*value - 1).clamp(*range.start(), *range.end());
        }
        if ui.add_sized([FRAME_HEIGHT, FRAME_HEIGHT], egui::Button::new("+")).clicked() {
            *value = (*value + 1).clamp(*range.start(), *range.end());
        }
        if *value != before {
            r.mark_changed();
        }
        ui.spacing_mut().item_spacing.x = 8.0;
        label_after(ui, label);
        r
    })
    .inner
}

/// ImGui::InputFloat with a step: the field, − and + buttons, the label.
pub fn input_float(ui: &mut egui::Ui, value: &mut f32, step: f32, decimals: usize, width: f32, range: std::ops::RangeInclusive<f32>, label: &str) -> egui::Response {
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 4.0;
        let field = (width - 2.0 * (FRAME_HEIGHT + 4.0)).max(30.0);
        field_style(ui);
        let mut r = ui.add_sized([field, FRAME_HEIGHT], egui::DragValue::new(value).speed(0.0).range(range.clone()).fixed_decimals(decimals));
        let before = *value;
        if ui.add_sized([FRAME_HEIGHT, FRAME_HEIGHT], egui::Button::new("-")).clicked() {
            *value = (*value - step).clamp(*range.start(), *range.end());
        }
        if ui.add_sized([FRAME_HEIGHT, FRAME_HEIGHT], egui::Button::new("+")).clicked() {
            *value = (*value + step).clamp(*range.start(), *range.end());
        }
        if *value != before {
            r.mark_changed();
        }
        ui.spacing_mut().item_spacing.x = 8.0;
        label_after(ui, label);
        r
    })
    .inner
}

/// Drag fields drawn on the field colour (ImGui frames), not the button colour.
pub fn field_style(ui: &mut egui::Ui) {
    let p = palette(ui);
    let v = &mut ui.visuals_mut().widgets;
    v.inactive.weak_bg_fill = p.frame;
    v.hovered.weak_bg_fill = p.frame_hover;
    v.active.weak_bg_fill = p.frame_active;
}

/// ImGui::CollapsingHeader: a full-width bar in the header colour with a
/// triangle at the left (the same calls as egui's CollapsingHeader).
pub struct CollapsingHeader {
    label: String,
    default_open: bool,
}

impl CollapsingHeader {
    pub fn new(label: impl Into<String>) -> Self {
        Self { label: label.into(), default_open: false }
    }

    pub fn default_open(mut self, open: bool) -> Self {
        self.default_open = open;
        self
    }

    pub fn show<R>(self, ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui) -> R) -> Option<R> {
        let p = palette(ui);
        let dark = ui.visuals().dark_mode;
        let id = ui.make_persistent_id(("collapsing header", self.label.as_str()));
        let mut open = ui.data_mut(|d| *d.get_persisted_mut_or(id, self.default_open));
        let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width(), FRAME_HEIGHT - 4.0), egui::Sense::click());
        if response.clicked() {
            open = !open;
            ui.data_mut(|d| d.insert_persisted(id, open));
        }
        let a = if response.hovered() { if dark { 87 } else { 56 } } else if dark { 61 } else { 36 };
        ui.painter().rect_filled(rect, CornerRadius::same(6), Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), a));
        let c = egui::pos2(rect.left() + 14.0, rect.center().y);
        let s = 4.5;
        let triangle = if open { vec![c + egui::vec2(-s, -s * 0.6), c + egui::vec2(s, -s * 0.6), c + egui::vec2(0.0, s * 0.8)] } else { vec![c + egui::vec2(-s * 0.6, -s), c + egui::vec2(s * 0.8, 0.0), c + egui::vec2(-s * 0.6, s)] };
        ui.painter().add(egui::Shape::convex_polygon(triangle, p.text, Stroke::NONE));
        ui.painter().text(egui::pos2(rect.left() + 28.0, rect.center().y), egui::Align2::LEFT_CENTER, &self.label, egui::TextStyle::Body.resolve(ui.style()), p.text);
        if open {
            Some(ui.push_id(id, |ui| add(ui)).inner)
        } else {
            None
        }
    }
}

/// ImGui::RadioButton: a frame-height circle, the accent dot when selected.
pub fn radio(ui: &mut egui::Ui, selected: bool, label: impl AsRef<str>) -> egui::Response {
    let p = palette(ui);
    let font = egui::TextStyle::Body.resolve(ui.style());
    let side = FRAME_HEIGHT - 4.0;
    // A label longer than the room left wraps instead of widening the window.
    let wrap = (ui.available_width() - side - 7.0).max(60.0);
    let galley = ui.painter().layout(label.as_ref().to_string(), font, p.text, wrap);
    let (rect, response) = ui.allocate_exact_size(egui::vec2(side + 7.0 + galley.size().x, FRAME_HEIGHT.max(galley.size().y + 4.0)), egui::Sense::click());
    let fade = |c: Color32| if ui.is_enabled() { c } else { c.gamma_multiply(0.5) };
    let centre = egui::pos2(rect.left() + side / 2.0, rect.top() + FRAME_HEIGHT / 2.0);
    let fill = if response.hovered() { p.frame_hover } else { p.frame };
    ui.painter().circle(centre, side / 2.0, fade(fill), Stroke::new(1.0, fade(p.border)));
    if selected {
        ui.painter().circle_filled(centre, side * 0.28, fade(p.accent));
    }
    ui.painter().galley(egui::pos2(rect.left() + side + 7.0, rect.center().y - galley.size().y / 2.0), galley, fade(p.text));
    response
}

/// A radio button setting `value` to `choice` when clicked.
pub fn radio_value<T: PartialEq>(ui: &mut egui::Ui, value: &mut T, choice: T, label: impl AsRef<str>) -> egui::Response {
    let mut r = radio(ui, *value == choice, label);
    if r.clicked() && *value != choice {
        *value = choice;
        r.mark_changed();
    }
    r
}

/// ImGui::SliderFloat/SliderInt with egui's Slider calls (new, fixed_decimals,
/// text, suffix): the framed slider, `slider_width` wide, the text after it.
pub struct Slider<'a> {
    get_set: Box<dyn FnMut(Option<f64>) -> f64 + 'a>,
    range: std::ops::RangeInclusive<f64>,
    decimals: Option<usize>,
    integer: bool,
    text: String,
    suffix: String,
}

impl<'a> Slider<'a> {
    pub fn new<N: egui::emath::Numeric>(value: &'a mut N, range: std::ops::RangeInclusive<N>) -> Self {
        let integer = N::INTEGRAL;
        Self {
            get_set: Box::new(move |v: Option<f64>| {
                if let Some(v) = v {
                    *value = N::from_f64(v);
                }
                value.to_f64()
            }),
            range: range.start().to_f64()..=range.end().to_f64(),
            decimals: None,
            integer,
            text: String::new(),
            suffix: String::new(),
        }
    }

    pub fn fixed_decimals(mut self, decimals: usize) -> Self {
        self.decimals = Some(decimals);
        self
    }

    pub fn text(mut self, text: impl Into<String>) -> Self {
        self.text = text.into();
        self
    }

    pub fn suffix(mut self, suffix: impl Into<String>) -> Self {
        self.suffix = suffix.into();
        self
    }
}

impl egui::Widget for Slider<'_> {
    fn ui(mut self, ui: &mut egui::Ui) -> egui::Response {
        let mut v = (self.get_set)(None);
        let shown = if self.integer { format!("{}{}", v.round() as i64, self.suffix) } else { format!("{:.*}{}", self.decimals.unwrap_or(3), v, self.suffix) };
        let width = ui.spacing().slider_width.max(120.0);
        let r = slider(ui, &mut v, self.range.clone(), width, &shown, &self.text, self.integer);
        if r.changed() {
            if let Some(d) = self.decimals.filter(|_| !self.integer) {
                let k = 10f64.powi(d as i32);
                v = (v * k).round() / k;
            }
            (self.get_set)(Some(v));
        }
        r
    }
}

/// ImGui::Selectable `width` wide: the label on the header colour when
/// selected or hovered.
pub fn selectable(ui: &mut egui::Ui, selected: bool, label: &str, width: f32) -> egui::Response {
    let p = palette(ui);
    let dark = ui.visuals().dark_mode;
    let height = crate::theme::BASE_FONT * 1.33 + 6.0;
    let (rect, response) = ui.allocate_exact_size(egui::vec2(width, height), egui::Sense::click());
    let a = if selected { if dark { 117 } else { 82 } } else if response.hovered() { if dark { 87 } else { 56 } } else { 0 };
    if a > 0 {
        ui.painter().rect_filled(rect, CornerRadius::same(4), Color32::from_rgba_unmultiplied(p.accent.r(), p.accent.g(), p.accent.b(), a));
    }
    ui.painter().text(egui::pos2(rect.left() + 4.0, rect.center().y), egui::Align2::LEFT_CENTER, label, egui::TextStyle::Body.resolve(ui.style()), p.text);
    response
}
