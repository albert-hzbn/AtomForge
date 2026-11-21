//! The look of the interface, as the C++ ui/ImGuiSetup: the platform UI font
//! at 18 px with a semibold heading face, the interface icons (a Codicons
//! subset, res/AtomForgeIcons.ttf) merged into the text font, compact softly
//! rounded geometry, and the light and dark palettes (neutral surfaces, one
//! accent, quiet borders). The light theme is the default, as in the C++.

use eframe::egui::{self, Color32, CornerRadius, FontFamily, FontId, Margin, Stroke, TextStyle};

/// The C++ kBaseFontSizePixels (18 px from ascent to descent, as ImGui sizes
/// fonts), as egui sizes them (by the em: 18 x 2048 / 2724 for Segoe UI).
pub const BASE_FONT: f32 = 13.5;
/// The C++ frame height: 18 px text and 6 px padding above and below.
pub const FRAME_HEIGHT: f32 = 30.0;

/// The interface icons (the C++ Icon namespace; Codicons glyph names after //).
#[allow(dead_code)]
pub mod icon {
    pub const OPEN: &str = "\u{EAF7}"; // folder-opened
    pub const SAVE: &str = "\u{EB4B}"; // save
    pub const UNDO: &str = "\u{EAE2}"; // discard
    pub const REDO: &str = "\u{EBB0}"; // redo
    pub const HOME: &str = "\u{EB06}"; // home
    pub const RULER: &str = "\u{EA99}"; // arrow-both
    pub const ROTATE: &str = "\u{EB37}"; // refresh
    pub const INFO: &str = "\u{EA74}"; // info
    pub const WARNING: &str = "\u{EA6C}"; // warning
    pub const ERROR: &str = "\u{EA87}"; // error
    pub const CHECK: &str = "\u{EBA4}"; // pass
    pub const ADD: &str = "\u{EA60}"; // add
    pub const DELETE: &str = "\u{EA81}"; // trash
    pub const SETTINGS: &str = "\u{EB51}"; // settings-gear
    pub const BUILD: &str = "\u{EA73}"; // edit
    pub const CUBE: &str = "\u{EB29}"; // package
    pub const HELP: &str = "\u{EB32}"; // question
    pub const CAMERA: &str = "\u{EADA}"; // device-camera
    pub const FOLDER: &str = "\u{EA83}"; // folder
    pub const DOCUMENT: &str = "\u{EA7B}"; // symbol-file
    pub const PLAY: &str = "\u{EB2C}"; // play
    pub const PIPELINE: &str = "\u{EBD2}"; // layers
    pub const SEARCH: &str = "\u{EA6D}"; // search
    pub const PAUSE: &str = "\u{EAD1}"; // debug-pause
    pub const PREVIOUS: &str = "\u{EB8E}"; // debug-reverse-continue
    pub const NEXT: &str = "\u{EACF}"; // debug-continue
}

/// "<icon>  <label>" (uiFonts::withIcon).
pub fn with_icon(icon: &str, label: &str) -> String {
    format!("{icon}  {label}")
}

/// The heading face (semibold when the platform has one).
pub fn heading_family() -> FontFamily {
    FontFamily::Name("heading".into())
}

pub fn heading(size_scale: f32) -> FontId {
    FontId::new(BASE_FONT * size_scale, heading_family())
}

/// The platform UI font (regular and semibold) and the icon font.
pub fn install_fonts(ctx: &egui::Context) {
    let (text, bold): (&[&str], &[&str]) = if cfg!(windows) {
        (&["C:/Windows/Fonts/segoeui.ttf", "C:/Windows/Fonts/arial.ttf"], &["C:/Windows/Fonts/seguisb.ttf", "C:/Windows/Fonts/segoeuib.ttf", "C:/Windows/Fonts/arialbd.ttf"])
    } else {
        (
            &["/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "/usr/share/fonts/truetype/freefont/FreeSans.ttf"],
            &["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf"],
        )
    };
    let mut fonts = egui::FontDefinitions::default();
    let mut proportional: Vec<String> = Vec::new();
    if let Some(bytes) = text.iter().find_map(|p| std::fs::read(p).ok()) {
        fonts.font_data.insert("system".into(), std::sync::Arc::new(egui::FontData::from_owned(bytes)));
        proportional.push("system".into());
    }
    // Icons share the text line: drawn on a 16-unit grid, slightly larger than
    // the text and nudged down to sit on its centre line (the C++ GlyphOffset).
    let icons = egui::FontData::from_static(include_bytes!("../res/AtomForgeIcons.ttf")).tweak(egui::FontTweak { scale: 19.0 / BASE_FONT, y_offset_factor: 0.16, ..Default::default() });
    fonts.font_data.insert("icons".into(), std::sync::Arc::new(icons));
    proportional.push("icons".into());
    let defaults = fonts.families.get(&FontFamily::Proportional).cloned().unwrap_or_default();
    let mut family = proportional.clone();
    family.extend(defaults.iter().cloned());
    fonts.families.insert(FontFamily::Proportional, family.clone());
    // Headings: the semibold face, else the text font.
    let mut heading = Vec::new();
    if let Some(bytes) = bold.iter().find_map(|p| std::fs::read(p).ok()) {
        fonts.font_data.insert("system-bold".into(), std::sync::Arc::new(egui::FontData::from_owned(bytes)));
        heading.push("system-bold".to_string());
    }
    heading.extend(family);
    fonts.families.insert(heading_family(), heading);
    // Monospace text keeps its font; the system and icon fonts fill in missing symbols.
    if let Some(mono) = fonts.families.get_mut(&FontFamily::Monospace) {
        mono.extend(proportional);
    }
    ctx.set_fonts(fonts);
}

fn rgb(r: f32, g: f32, b: f32) -> Color32 {
    Color32::from_rgb((r * 255.0).round() as u8, (g * 255.0).round() as u8, (b * 255.0).round() as u8)
}

fn alpha(c: Color32, a: f32) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), (a * 255.0).round() as u8)
}

/// The palette of a theme (the C++ Palette; every colour is kept, as in the C++).
#[allow(dead_code)]
#[derive(Clone, Copy)]
pub struct Palette {
    pub text: Color32,
    pub text_muted: Color32,
    pub window: Color32,
    pub surface: Color32,
    pub popup: Color32,
    pub border: Color32,
    pub frame: Color32,
    pub frame_hover: Color32,
    pub frame_active: Color32,
    pub title: Color32,
    pub title_active: Color32,
    pub menu_bar: Color32,
    pub button: Color32,
    pub button_hover: Color32,
    pub button_active: Color32,
    pub accent: Color32,
    pub accent_strong: Color32,
    pub separator: Color32,
    pub dim: Color32,
}

pub fn light_palette() -> Palette {
    Palette {
        text: rgb(0.11, 0.12, 0.14),
        text_muted: rgb(0.40, 0.43, 0.48),
        window: rgb(0.965, 0.969, 0.976),
        surface: rgb(1.0, 1.0, 1.0),
        popup: rgb(1.0, 1.0, 1.0),
        border: rgb(0.800, 0.820, 0.850),
        frame: rgb(1.0, 1.0, 1.0),
        frame_hover: rgb(0.945, 0.955, 0.972),
        frame_active: rgb(0.905, 0.930, 0.965),
        title: rgb(0.930, 0.937, 0.950),
        title_active: rgb(0.900, 0.912, 0.932),
        menu_bar: rgb(0.975, 0.978, 0.984),
        button: rgb(0.935, 0.945, 0.960),
        button_hover: rgb(0.885, 0.905, 0.935),
        button_active: rgb(0.830, 0.865, 0.915),
        accent: rgb(0.12, 0.44, 0.75),
        accent_strong: rgb(0.08, 0.35, 0.64),
        separator: rgb(0.865, 0.880, 0.900),
        dim: Color32::from_rgba_unmultiplied(26, 31, 41, 77),
    }
}

pub fn dark_palette() -> Palette {
    Palette {
        text: rgb(0.90, 0.92, 0.95),
        text_muted: rgb(0.56, 0.60, 0.66),
        window: rgb(0.105, 0.115, 0.135),
        surface: rgb(0.125, 0.137, 0.160),
        popup: rgb(0.135, 0.147, 0.172),
        border: rgb(0.235, 0.260, 0.300),
        frame: rgb(0.155, 0.170, 0.198),
        frame_hover: rgb(0.185, 0.203, 0.237),
        frame_active: rgb(0.215, 0.240, 0.282),
        title: rgb(0.115, 0.126, 0.148),
        title_active: rgb(0.140, 0.153, 0.180),
        menu_bar: rgb(0.120, 0.131, 0.154),
        button: rgb(0.180, 0.197, 0.230),
        button_hover: rgb(0.225, 0.248, 0.290),
        button_active: rgb(0.265, 0.295, 0.345),
        accent: rgb(0.33, 0.63, 0.93),
        accent_strong: rgb(0.45, 0.72, 0.98),
        separator: rgb(0.215, 0.238, 0.275),
        dim: Color32::from_rgba_unmultiplied(5, 8, 13, 140),
    }
}

impl Palette {
    pub fn window_fill(&self) -> Color32 {
        self.window
    }
}

pub fn palette(dark: bool) -> Palette {
    if dark {
        dark_palette()
    } else {
        light_palette()
    }
}

/// The accent colour of the current theme (dialogLayout::accent).
#[allow(dead_code)]
pub fn accent(ui: &egui::Ui) -> Color32 {
    palette(ui.visuals().dark_mode).accent
}

/// Muted text (ImGui TextDisabled).
#[allow(dead_code)]
pub fn muted(ui: &egui::Ui) -> Color32 {
    palette(ui.visuals().dark_mode).text_muted
}

fn style_for(dark: bool, style: &mut egui::Style) {
    let p = palette(dark);
    // Text: the platform font at 18 px (headings semibold).
    style.text_styles = [
        (TextStyle::Small, FontId::proportional(BASE_FONT * 0.82)),
        (TextStyle::Body, FontId::proportional(BASE_FONT)),
        (TextStyle::Button, FontId::proportional(BASE_FONT)),
        (TextStyle::Heading, FontId::new(BASE_FONT * 1.2, heading_family())),
        (TextStyle::Monospace, FontId::monospace(BASE_FONT * 0.85)),
    ]
    .into();
    // Geometry shared by both themes: compact, evenly spaced, softly rounded.
    let s = &mut style.spacing;
    s.item_spacing = egui::vec2(8.0, 8.0);
    s.button_padding = egui::vec2(12.0, 6.0);
    s.window_margin = Margin::symmetric(16, 14);
    s.menu_margin = Margin::symmetric(8, 6);
    s.indent = 18.0;
    s.interact_size = egui::vec2(40.0, FRAME_HEIGHT);
    // ImGui check boxes and radio buttons are as tall as a frame.
    s.icon_width = 22.0;
    s.icon_width_inner = 12.0;
    s.icon_spacing = 7.0;
    s.combo_height = 400.0;
    s.scroll.bar_width = 12.0;
    s.slider_rail_height = 6.0;

    let v = &mut style.visuals;
    *v = if dark { egui::Visuals::dark() } else { egui::Visuals::light() };
    v.override_text_color = None;
    v.weak_text_color = Some(p.text_muted);
    v.hyperlink_color = p.accent;
    v.faint_bg_color = if dark { Color32::from_rgba_unmultiplied(255, 255, 255, 8) } else { Color32::from_rgba_unmultiplied(0, 0, 0, 6) };
    v.extreme_bg_color = p.frame;
    v.text_edit_bg_color = Some(p.frame);
    v.code_bg_color = p.button;
    v.window_fill = p.window;
    v.panel_fill = p.window;
    v.window_stroke = Stroke::new(1.0, p.border);
    v.window_corner_radius = CornerRadius::same(10);
    v.menu_corner_radius = CornerRadius::same(8);
    v.window_shadow = egui::Shadow { offset: [0, 4], blur: 16, spread: 0, color: Color32::from_black_alpha(if dark { 90 } else { 28 }) };
    v.popup_shadow = egui::Shadow { offset: [0, 3], blur: 10, spread: 0, color: Color32::from_black_alpha(if dark { 80 } else { 24 }) };
    v.selection.bg_fill = alpha(p.accent, if dark { 0.34 } else { 0.22 });
    v.selection.stroke = Stroke::new(1.0, p.accent);
    v.error_fg_color = if dark { rgb(0.95, 0.45, 0.42) } else { rgb(0.78, 0.18, 0.16) };
    v.warn_fg_color = if dark { rgb(0.95, 0.74, 0.36) } else { rgb(0.70, 0.45, 0.05) };
    v.striped = false;
    v.slider_trailing_fill = true;
    v.indent_has_left_vline = false;

    let radius = CornerRadius::same(6);
    let w = &mut v.widgets;
    // Labels and separators.
    w.noninteractive.bg_fill = p.window;
    w.noninteractive.weak_bg_fill = p.surface;
    w.noninteractive.bg_stroke = Stroke::new(1.0, p.separator);
    w.noninteractive.fg_stroke = Stroke::new(1.0, p.text);
    w.noninteractive.corner_radius = radius;
    // Buttons (weak_bg_fill), and the frames of fields and checkboxes (bg_fill).
    w.inactive.weak_bg_fill = p.button;
    w.inactive.bg_fill = p.frame;
    w.inactive.bg_stroke = Stroke::new(1.0, p.border);
    w.inactive.fg_stroke = Stroke::new(1.0, p.text);
    w.inactive.corner_radius = radius;
    w.inactive.expansion = 0.0;
    w.hovered.weak_bg_fill = p.button_hover;
    w.hovered.bg_fill = p.frame_hover;
    w.hovered.bg_stroke = Stroke::new(1.0, alpha(p.accent, 0.55));
    w.hovered.fg_stroke = Stroke::new(1.5, p.text);
    w.hovered.corner_radius = radius;
    w.hovered.expansion = 0.0;
    w.active.weak_bg_fill = p.button_active;
    w.active.bg_fill = p.frame_active;
    w.active.bg_stroke = Stroke::new(1.0, p.accent);
    w.active.fg_stroke = Stroke::new(2.0, p.text);
    w.active.corner_radius = radius;
    w.active.expansion = 0.0;
    w.open.weak_bg_fill = p.button_active;
    w.open.bg_fill = p.frame_active;
    w.open.bg_stroke = Stroke::new(1.0, p.border);
    w.open.fg_stroke = Stroke::new(1.0, p.text);
    w.open.corner_radius = radius;
}

/// Applies both themes and makes the light one current (the C++ default).
pub fn apply(ctx: &egui::Context) {
    ctx.style_mut_of(egui::Theme::Light, |s| style_for(false, s));
    ctx.style_mut_of(egui::Theme::Dark, |s| style_for(true, s));
    ctx.set_theme(egui::Theme::Light);
}

/// The primary (accent-filled) button of a dialog (dialogLayout::primaryButton).
pub fn primary_button(ui: &mut egui::Ui, text: impl Into<String>, size: Option<egui::Vec2>) -> egui::Response {
    let p = palette(ui.visuals().dark_mode);
    let label = egui::RichText::new(text.into()).color(Color32::WHITE);
    let mut button = egui::Button::new(label).fill(p.accent).stroke(Stroke::new(1.0, p.accent_strong));
    if let Some(s) = size {
        button = button.min_size(s);
    }
    ui.add(button)
}

/// A secondary button (dialogLayout::secondaryButton).
pub fn secondary_button(ui: &mut egui::Ui, text: impl Into<String>, size: Option<egui::Vec2>) -> egui::Response {
    let mut button = egui::Button::new(text.into());
    if let Some(s) = size {
        button = button.min_size(s);
    }
    ui.add(button)
}

/// A section heading with a rule under it (dialogLayout::section).
pub fn section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(6.0);
    ui.label(egui::RichText::new(title).font(heading(1.0)));
    let r = ui.available_rect_before_wrap();
    ui.painter().hline(r.x_range(), r.top(), Stroke::new(1.0, palette(ui.visuals().dark_mode).separator));
    ui.add_space(6.0);
}

/// The style of menus (ImGui menus): compact rows, no frames on the items.
pub fn menu_style(style: &mut egui::Style) {
    egui::containers::menu::menu_style(style);
    style.spacing.interact_size.y = 24.0;
    style.spacing.item_spacing.y = 0.0;
    style.spacing.button_padding = egui::vec2(8.0, 2.0);
    style.spacing.menu_margin = Margin::symmetric(6, 6);
}

/// A menu item that toggles a setting, with a check mark at the right when
/// on (ImGui::MenuItem with `selected`).
pub fn menu_toggle(ui: &mut egui::Ui, enabled: bool, value: &mut bool, label: &str) -> egui::Response {
    let mut response = ui.add_enabled(enabled, egui::Button::new(label).shortcut_text("    "));
    check_mark(ui, &response, *value);
    if response.clicked() {
        *value = !*value;
        response.mark_changed();
    }
    response
}

/// A menu item for one of several choices, checked when it is the current one.
pub fn menu_radio(ui: &mut egui::Ui, enabled: bool, selected: bool, label: impl Into<egui::WidgetText>) -> egui::Response {
    let response = ui.add_enabled(enabled, egui::Button::new(label).shortcut_text("    "));
    check_mark(ui, &response, selected);
    response
}

/// ImGui::RenderCheckMark at the right of a menu item (the text fonts have
/// no check mark glyph).
fn check_mark(ui: &egui::Ui, response: &egui::Response, on: bool) {
    if !on {
        return;
    }
    let size = ui.text_style_height(&egui::TextStyle::Body) * 0.62;
    let r = response.rect;
    let x = r.right() - ui.spacing().button_padding.x - size;
    let y = r.center().y - size * 0.5;
    let colour = ui.style().interact(response).text_color();
    let points = vec![egui::pos2(x, y + size * 0.55), egui::pos2(x + size * 0.36, y + size * 0.9), egui::pos2(x + size, y + size * 0.1)];
    ui.painter().add(egui::Shape::line(points, egui::Stroke::new((size * 0.16).max(1.5), colour)));
}

thread_local! {
    static SCREEN: std::cell::Cell<egui::Rect> = const { std::cell::Cell::new(egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1280.0, 800.0))) };
}

thread_local! {
    static DARK: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// The window's area and theme this frame (dialogs open centred in it, with the theme's title bar).
pub fn set_screen(rect: egui::Rect, dark: bool) {
    SCREEN.with(|s| s.set(rect));
    DARK.with(|d| d.set(dark));
}

/// A dialog window as the C++ ones: opened centred on the screen, a compact
/// title bar with the title at the left in the text font, no collapse arrow.
pub fn window(title: impl Into<String>) -> egui::Window<'static> {
    let title: String = title.into();
    let mut atoms = egui::Atoms::new(egui::RichText::new(title).font(FontId::proportional(BASE_FONT)));
    // egui centres the title between two growing spaces; a wide spacer after it,
    // shrunk to the room that is left, leaves them nothing, so the title sits at
    // the left as in ImGui (without widening the window).
    let mut spacer = egui::Atom::from(egui::RichText::new(" ".repeat(600)).color(Color32::TRANSPARENT));
    spacer.shrink = true;
    atoms.push_right(spacer);
    let centre = SCREEN.with(|s| s.get()).center();
    egui::Window::new(atoms)
        .collapsible(false)
        .pivot(egui::Align2::CENTER_CENTER)
        .default_pos(centre)
        .title_frame(egui::Frame::new().fill(palette(DARK.with(|d| d.get())).title_active).inner_margin(Margin { left: 10, right: 6, top: 4, bottom: 4 }))
}

thread_local! {
    static MODAL_SHOWN: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Modal dialogs drawn this frame (title, window id).
    static MODALS: std::cell::RefCell<Vec<(String, egui::Id)>> = const { std::cell::RefCell::new(Vec::new()) };
    /// Minimised dialogs: title and where the window was.
    static MINIMISED: std::cell::RefCell<Vec<(String, egui::Pos2)>> = const { std::cell::RefCell::new(Vec::new()) };
    /// A dialog restored from the bar: put back where it was on its next frame.
    static RESTORE: std::cell::RefCell<Option<(String, egui::Pos2)>> = const { std::cell::RefCell::new(None) };
}

/// The dialogs that can be minimised (responsive::minimizable): builders and analyses.
const MINIMISABLE: [&str; 17] = [
    "Build Bulk Crystal", "Amorphous Structure Builder", "Build Nanocrystal", "Build Polycrystal",
    "CSL Grain Boundary Builder", "Interface Builder", "Stacking Faults", "Custom Structure",
    "Substitutional Solid Solution", "Insert Dislocation", "Add Interstitial Atoms", "Merge Structures",
    "Cell Sculptor", "Radial Distribution Function", "Angular Distribution Function",
    "Common Neighbour Analysis", "Short-Range Order Analysis",
];

/// A modal dialog (responsive::beginModal): the view behind is dimmed and
/// takes no input while it is open (see `modal_backdrop`). Builders and
/// analyses also get a minimise button: hidden, settings kept, restored from
/// the bar above the status bar.
pub fn modal_window(title: impl Into<String>) -> egui::Window<'static> {
    let title: String = title.into();
    let minimisable = MINIMISABLE.contains(&title.as_str());
    let id = egui::Id::unique(("af dialog", title.clone()));
    let mut w = window(title.clone()).order(egui::Order::Foreground);
    if minimisable {
        w = w.id(id);
        if MINIMISED.with(|m| m.borrow().iter().any(|(t, _)| *t == title)) {
            // Out of sight while minimised (its contents keep running and keep their state).
            return w.current_pos(egui::pos2(-100_000.0, -100_000.0));
        }
        if let Some(pos) = RESTORE.with(|r| r.borrow_mut().take_if(|(t, _)| *t == title).map(|(_, p)| p)) {
            w = w.current_pos(pos).pivot(egui::Align2::LEFT_TOP);
        }
        MODALS.with(|m| m.borrow_mut().push((title, id)));
    }
    MODAL_SHOWN.with(|m| m.set(true));
    w
}

/// Dims the window behind open modal dialogs and absorbs its clicks, draws
/// the minimise buttons on their title bars and the bar of minimised dialogs;
/// call once per frame after every dialog was drawn.
pub fn modal_backdrop(ctx: &egui::Context) {
    let dark = ctx.theme() == egui::Theme::Dark;
    let p = palette(dark);
    let modals = MODALS.with(|m| std::mem::take(&mut *m.borrow_mut()));
    // Minimise buttons, left of each dialog's close cross (titleBarMinimizeButton).
    let heading = ctx.fonts_mut(|f| f.row_height(&FontId::new(BASE_FONT * 1.2, heading_family())));
    for (title, id) in &modals {
        let Some(rect) = ctx.memory(|m| m.area_rect(*id)) else { continue };
        let size = heading;
        let close_left = rect.right() - 6.0 - size;
        let bar_top = rect.top() + 1.0;
        let button = egui::Rect::from_min_size(egui::pos2(close_left - 4.0 - size, bar_top + 4.0), egui::vec2(size, size));
        let response = egui::Area::new(egui::Id::unique(("minimise", title.clone()))).order(egui::Order::Foreground).fixed_pos(button.min).show(ctx, |ui| {
            let (r, response) = ui.allocate_exact_size(button.size(), egui::Sense::click());
            if response.hovered() {
                ui.painter().rect_filled(r, 6.0, if response.is_pointer_button_down_on() { p.button_active } else { p.button_hover });
            }
            let inset = r.width() * 0.25;
            let y = r.bottom() - r.height() * 0.32;
            ui.painter().line_segment([egui::pos2(r.left() + inset, y), egui::pos2(r.right() - inset, y)], egui::Stroke::new((r.width() * 0.09).max(1.0), p.text));
            response.on_hover_text("Minimize: hide the dialog and keep its settings; restore it from the bar above the status bar")
        });
        ctx.move_to_top(response.response.layer_id);
        if response.inner.clicked() {
            MINIMISED.with(|m| m.borrow_mut().push((title.clone(), rect.min)));
        }
    }
    // The bar of minimised dialogs above the status bar (drawMinimizedDialogs).
    let minimised = MINIMISED.with(|m| m.borrow().clone());
    if !minimised.is_empty() {
        let screen = ctx.content_rect();
        let bottom = screen.bottom() - FRAME_HEIGHT - 10.0;
        egui::Area::new(egui::Id::unique("minimised dialogs")).order(egui::Order::Foreground).pivot(egui::Align2::CENTER_BOTTOM).fixed_pos(egui::pos2(screen.center().x, bottom)).show(ctx, |ui| {
            egui::Frame::new().fill(p.popup.gamma_multiply(0.92)).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(10).inner_margin(egui::Margin::symmetric(12, 8)).show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new("Minimized:").color(p.text_muted));
                    for (title, pos) in &minimised {
                        if ui.button(title).on_hover_text(format!("Restore {title}")).clicked() {
                            MINIMISED.with(|m| m.borrow_mut().retain(|(t, _)| t != title));
                            RESTORE.with(|r| *r.borrow_mut() = Some((title.clone(), *pos)));
                        }
                    }
                });
            });
        });
    }
    if !MODAL_SHOWN.with(|m| m.replace(false)) {
        return;
    }
    let dim = p.dim;
    let screen = ctx.content_rect();
    let response = egui::Area::new(egui::Id::unique("modal backdrop")).order(egui::Order::Middle).fixed_pos(screen.min).interactable(true).show(ctx, |ui| {
        let (rect, response) = ui.allocate_exact_size(screen.size(), egui::Sense::click_and_drag());
        ui.painter().rect_filled(rect, 0.0, dim);
        response
    });
    ctx.move_to_top(response.response.layer_id);
}

/// Whether a dialog is minimised (its window is out of sight); UI scripts.
pub fn minimise(title: &str, minimised: bool) {
    MINIMISED.with(|m| {
        let mut m = m.borrow_mut();
        m.retain(|(t, _)| t != title);
        if minimised {
            m.push((title.to_string(), egui::pos2(100.0, 100.0)));
        }
    });
}

/// A colour from 0..1 sRGB components, as ImGui (and the C++) take them;
/// egui's `Rgba` would read them as linear light and lighten them.
pub fn srgb(r: f32, g: f32, b: f32) -> Color32 {
    srgba(r, g, b, 1.0)
}

/// `srgb` with straight (unmultiplied) alpha.
pub fn srgba(r: f32, g: f32, b: f32, a: f32) -> Color32 {
    let c = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
    Color32::from_rgba_unmultiplied(c(r), c(g), c(b), c(a))
}
