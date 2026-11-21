//! Settings > Display Settings and Settings > Atomic Sizes (the C++
//! ui/EditMenuDialogs drawDisplaySettings and drawRadiusEditor): looks,
//! background and cell colours, depth cues, outlines, lighting, atom and
//! bond sizes and rules, element radii and colours, saved as the default.

use af_core::display::{self as d, apply_preset, matching_preset, save_defaults, saved_defaults_path, DisplaySettings, Look, Range, PRESETS};
use af_core::elements;
use af_core::structure::Structure;
use eframe::egui;

use crate::dialog_layout::{self as dl, Note};
use crate::view::AtomStyle;

/// What a frame of the window changed.
#[derive(Default)]
pub struct Changes {
    /// Sizes, bonds or element tables: the scene is built again.
    pub rebuild: bool,
    /// Element colours: the atoms are coloured again.
    pub recolour: bool,
    pub style: Option<AtomStyle>,
}

pub struct DisplaySettingsWindow {
    pub open: bool,
    pub sizes_open: bool,
    tab: usize,
    colour_element: i32,
    radius_element: i32,
    status: String,
    rebuild_pending: bool,
}

impl Default for DisplaySettingsWindow {
    fn default() -> Self {
        Self { open: false, sizes_open: false, tab: 0, colour_element: 6, radius_element: 6, status: String::new(), rebuild_pending: false }
    }
}

const TABS: [&str; 4] = ["Style", "Lighting", "Atoms and Bonds", "Element Colours"];
/// The C++ labelColumn() and sliderWidth(): 10.5 and 17 font sizes (18 px).
const LABEL_COLUMN: f32 = 18.0 * 10.5;
const SLIDER_WIDTH: f32 = 18.0 * 17.0;
/// The width of the look and colour-scheme buttons (9 font sizes).
const PRESET_WIDTH: f32 = 18.0 * 9.0;

fn differs(value: f32, reference: f32) -> bool {
    (value - reference).abs() > 1e-4 * 1.0f32.max(reference.abs())
}

/// A label in the left column, then the field (dialogLayout::fieldLabel).
fn field_label(ui: &mut egui::Ui, label: &str) -> egui::Response {
    let (rect, response) = ui.allocate_exact_size(egui::vec2(LABEL_COLUMN, crate::theme::FRAME_HEIGHT), egui::Sense::hover());
    let font = egui::TextStyle::Body.resolve(ui.style());
    ui.painter().text(egui::pos2(rect.left(), rect.center().y), egui::Align2::LEFT_CENTER, label, font, ui.visuals().text_color());
    response
}

/// One setting: its name, a slider, a Reset button (enabled while the value
/// differs from the default) and its explanation on hover (settingSlider).
#[allow(clippy::too_many_arguments)]
fn setting(ui: &mut egui::Ui, label: &str, value: &mut f32, range: Range, factory: f32, decimals: usize, suffix: &str, help: &str) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        field_label(ui, label).on_hover_text(help);
        changed |= dl::slider_f32(ui, value, range.low..=range.high, SLIDER_WIDTH, |v| format!("{v:.decimals$}{suffix}")).on_hover_text(help).changed();
        if dl::small_button(ui, "Reset", differs(*value, factory)).on_hover_text(format!("Back to the default, {}.", af_core::format::format_g(factory as f64, 6))).clicked() {
            *value = factory;
            changed = true;
        }
    });
    changed
}

/// "Theme colour" or a chosen one, with a few common swatches (colourChoice).
fn colour_choice(ui: &mut egui::Ui, id: &str, custom: &mut bool, colour: &mut [f32; 3]) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        changed |= crate::dialog_widgets::radio_value(ui, custom, false, "Theme colour").changed();
        changed |= crate::dialog_widgets::radio_value(ui, custom, true, "Custom colour").changed();
    });
    if *custom {
        ui.horizontal(|ui| {
            changed |= dl::colour_edit3(ui, (id, "colour"), colour, SLIDER_WIDTH, "");
            let swatches: [(&str, [f32; 3]); 6] =
                [("White", [1.0; 3]), ("Light grey", [0.82; 3]), ("Mid grey", [0.5; 3]), ("Dark grey", [0.18; 3]), ("Black", [0.0; 3]), ("Navy", [0.06, 0.09, 0.2])];
            for (name, c) in swatches {
                let side = crate::theme::FRAME_HEIGHT;
                let (rect, response) = ui.allocate_exact_size(egui::vec2(side, side), egui::Sense::click());
                ui.painter().rect_filled(rect, 4.0, srgb(c));
                ui.painter().rect_stroke(rect, 4.0, egui::Stroke::new(1.0, crate::theme::palette(ui.visuals().dark_mode).border), egui::StrokeKind::Inside);
                if response.on_hover_text(name).clicked() {
                    *colour = c;
                    changed = true;
                }
            }
        });
    }
    changed
}

/// A stored colour (sRGB components 0-1) for painting.
fn srgb(c: [f32; 3]) -> egui::Color32 {
    let byte = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
    egui::Color32::from_rgb(byte(c[0]), byte(c[1]), byte(c[2]))
}

/// A look or colour-scheme button: the current one as the primary button.
fn choice_button(ui: &mut egui::Ui, name: &str, selected: bool) -> egui::Response {
    let size = Some(egui::vec2(PRESET_WIDTH, crate::theme::FRAME_HEIGHT));
    if selected {
        crate::theme::primary_button(ui, name, size)
    } else {
        crate::theme::secondary_button(ui, name, size)
    }
}

/// Row and column of an element in the periodic table (lanthanides and
/// actinides in rows 8 and 9).
pub fn table_position(z: i32) -> (usize, usize) {
    let z = z as usize;
    match z {
        1 => (0, 0),
        2 => (0, 17),
        3..=4 => (1, z - 3),
        5..=10 => (1, z + 7),
        11..=12 => (2, z - 11),
        13..=18 => (2, z - 1),
        19..=36 => (3, z - 19),
        37..=54 => (4, z - 37),
        55..=56 => (5, z - 55),
        57..=71 => (7, z - 57 + 2),
        72..=86 => (5, z - 69),
        87..=88 => (6, z - 87),
        89..=103 => (8, z - 89 + 2),
        _ => (6, z - 101),
    }
}

/// The C++ periodic table's category of an element (kCatColors index).
fn category(z: i32) -> usize {
    match z {
        1 => 0,
        2 | 10 | 18 | 36 | 54 | 86 | 118 => 10,
        3 | 11 | 19 | 37 | 55 | 87 => 1,
        4 | 12 | 20 | 38 | 56 | 88 => 2,
        57..=71 => 3,
        89..=103 => 4,
        9 | 17 | 35 | 53 | 85 | 117 => 9,
        6 | 7 | 8 | 15 | 16 | 34 => 8,
        5 | 14 | 32 | 33 | 51 | 52 => 7,
        13 | 31 | 49 | 50 | 81 | 82 | 83 | 84 | 113 | 114 | 115 | 116 => 6,
        _ => 5,
    }
}

/// H, alkali and alkaline-earth metals, lanthanides, actinides, transition
/// and post-transition metals, metalloids, nonmetals, halogens, noble gases.
const CATEGORY_COLOURS: [[f32; 3]; 11] = [
    [0.78, 0.78, 0.20],
    [0.92, 0.40, 0.25],
    [0.92, 0.76, 0.30],
    [1.00, 0.72, 0.77],
    [1.00, 0.78, 0.60],
    [0.72, 0.72, 0.95],
    [0.55, 0.85, 0.80],
    [0.72, 0.90, 0.55],
    [0.60, 0.90, 0.58],
    [0.55, 0.90, 0.70],
    [0.50, 0.78, 0.95],
];

fn brighter(c: [f32; 3], d: f32) -> [f32; 3] {
    [(c[0] + d).min(1.0), (c[1] + d).min(1.0), (c[2] + d).min(1.0)]
}

/// The periodic table in a framed panel, coloured by element category
/// (drawPeriodicTableInlineSelector); clicking an element selects it. The
/// colours argument is kept for callers; the table shows categories, as the C++.
pub fn periodic_table(ui: &mut egui::Ui, selected: &mut i32, _colours: Option<&[[f32; 3]]>) -> bool {
    let mut changed = false;
    let p = crate::theme::palette(ui.visuals().dark_mode);
    let (cell_w, cell_h) = (18.0f32 * 2.8, 18.0f32 * 2.1);
    let (pad_x, pad_y, f_gap) = (6.0f32, 10.0f32, 10.0f32);
    let row_y = |row: usize| if row <= 6 { pad_y + row as f32 * cell_h } else { pad_y + 7.0 * cell_h + f_gap + (row - 7) as f32 * cell_h };
    let size = egui::vec2(pad_x * 2.0 + 18.0 * cell_w + 32.0, pad_y + 7.0 * cell_h + f_gap + 2.0 * cell_h + 10.0 + 28.0);
    let frame = egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(16, 14));
    frame.show(ui, |ui| {
        let (rect, _) = ui.allocate_exact_size(size - egui::vec2(32.0, 28.0), egui::Sense::hover());
        let font = egui::FontId::proportional(crate::theme::BASE_FONT);
        let paint = |c: [f32; 3]| egui::Color32::from_rgb((c[0] * 255.0) as u8, (c[1] * 255.0) as u8, (c[2] * 255.0) as u8);
        for z in 1..=118 {
            let (row, col) = table_position(z);
            let min = egui::pos2(rect.min.x + pad_x + col as f32 * cell_w, rect.min.y + row_y(row));
            let r = egui::Rect::from_min_size(min, egui::vec2(cell_w - 1.0, cell_h - 1.0));
            let response = ui.interact(r, ui.scope_id().with(("element", z)), egui::Sense::click());
            let on = *selected == z;
            let base = CATEGORY_COLOURS[category(z)];
            let bg = if on { brighter(base, 0.35) } else if response.hovered() { brighter(base, 0.15) } else { base };
            ui.painter().rect_filled(r, 3.0, paint(bg));
            if on {
                ui.painter().rect_stroke(r, 3.0, egui::Stroke::new(2.0, egui::Color32::from_rgba_unmultiplied(0, 0, 100, 200)), egui::StrokeKind::Inside);
            } else {
                ui.painter().rect_stroke(r, 3.0, egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(0, 0, 0, 55)), egui::StrokeKind::Inside);
            }
            ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, elements::symbol(z), font.clone(), egui::Color32::BLACK);
            if response.on_hover_text(format!("{z}  {}\n{}", elements::symbol(z), elements::name(z))).clicked() {
                *selected = z;
                changed = true;
            }
        }
        // The places of the f-block rows, tinted like them.
        for (row, label, cat) in [(5usize, "57-71", 3usize), (6, "89-103", 4)] {
            let min = egui::pos2(rect.min.x + pad_x + 2.0 * cell_w, rect.min.y + row_y(row));
            let r = egui::Rect::from_min_size(min, egui::vec2(cell_w - 1.0, cell_h - 1.0));
            let c = CATEGORY_COLOURS[cat];
            ui.painter().rect_filled(r, 3.0, egui::Color32::from_rgba_unmultiplied((c[0] * 255.0) as u8, (c[1] * 255.0) as u8, (c[2] * 255.0) as u8, 140));
            ui.painter().rect_stroke(r, 3.0, egui::Stroke::new(1.0, egui::Color32::from_rgba_unmultiplied(0, 0, 0, 55)), egui::StrokeKind::Inside);
            ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, label, egui::FontId::proportional(crate::theme::BASE_FONT * 0.8), egui::Color32::from_rgb(40, 40, 40));
        }
    });
    changed
}

impl DisplaySettingsWindow {
    pub fn select_tab(&mut self, name: &str) -> bool {
        let key = name.to_lowercase();
        let found = TABS.iter().position(|t| t.to_lowercase().starts_with(&key) || (key == "colors" && *t == "Element Colours") || (key == "atoms" && *t == "Atoms and Bonds"));
        if let Some(i) = found {
            self.tab = i;
        }
        found.is_some()
    }

    /// The radius of the picked element (drawRadiusEditor).
    fn radius_editor(&mut self, ui: &mut egui::Ui, radii: &mut [f32], changes: &mut Changes) {
        dl::muted(ui, "Pick an element and set the radius its atoms are drawn with. Defaults: covalent radii of Cordero et al., Dalton Trans. 2008.");
        ui.add_space(4.0);
        periodic_table(ui, &mut self.radius_element, None);
        let z = self.radius_element;
        if (1..=118).contains(&z) {
            dl::section(ui, &format!("{} (Z = {z})", elements::symbol(z)));
            ui.horizontal(|ui| {
                let p = crate::theme::palette(ui.visuals().dark_mode);
                let r = ui.scope(|ui| {
                    let w = &mut ui.visuals_mut().widgets;
                    w.inactive.weak_bg_fill = p.frame;
                    w.hovered.weak_bg_fill = p.frame_hover;
                    w.active.weak_bg_fill = p.frame_active;
                    ui.add_sized([240.0, crate::theme::FRAME_HEIGHT], egui::DragValue::new(&mut radii[z as usize]).range(0.2..=3.5).speed(0.01).fixed_decimals(2).suffix(" Å"))
                });
                if r.inner.changed() {
                    changes.rebuild = true;
                }
                ui.label("Radius");
                if ui.button("Reset Element").clicked() {
                    radii[z as usize] = elements::literature_covalent_radii()[z as usize];
                    changes.rebuild = true;
                }
            });
        } else {
            dl::muted(ui, "No element selected.");
        }
    }

    /// Settings > Atomic Sizes.
    pub fn show_sizes(&mut self, ctx: &egui::Context, radii: &mut [f32]) -> Changes {
        let mut changes = Changes::default();
        if !self.sizes_open {
            return changes;
        }
        let mut open = true;
        crate::theme::modal_window("Atomic Sizes").open(&mut open).resizable(false).show(ctx, |ui| {
            self.radius_editor(ui, radii, &mut changes);
            match dl::footer(ui, "Close", Some("Reset All Elements"), true) {
                2 => {
                    radii.copy_from_slice(&elements::literature_covalent_radii()[..radii.len()]);
                    changes.rebuild = true;
                }
                1 => self.sizes_open = false,
                _ => {}
            }
        });
        self.sizes_open = open && self.sizes_open;
        changes
    }

    /// Settings > Display Settings. `colours` and `radii` are the open
    /// structure's element tables; `style` the View > Atom Style.
    pub fn show(&mut self, ctx: &egui::Context, display: &mut DisplaySettings, colours: &mut [[f32; 3]], radii: &mut [f32], style: AtomStyle, structure: &Structure) -> Changes {
        let mut changes = Changes::default();
        if !self.open {
            return changes;
        }
        let factory = DisplaySettings::default();
        let mut rebuild = false;
        let mut open = true;
        let max_height = ctx.content_rect().height() * 0.72;
        crate::theme::modal_window("Display Settings").open(&mut open).resizable(false).show(ctx, |ui| {
            dl::tab_bar(ui, &mut self.tab, &TABS);
            egui::ScrollArea::vertical().max_height(max_height).auto_shrink([true, true]).show(ui, |ui| match self.tab {
                0 => rebuild |= self.style_tab(ui, display, &factory),
                1 => self.lighting_tab(ui, display, &factory),
                2 => rebuild |= self.atoms_tab(ui, display, &factory, radii, style, structure, &mut changes),
                _ => self.colours_tab(ui, colours, &mut changes),
            });
            // Footer: keep this look for new windows, back to the factory look, close.
            ui.add_space(4.0);
            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("Save as My Default").on_hover_text("Keep every setting of this dialog, element colours and radii included,\nfor new windows and the next start.").clicked() {
                    let look = Look { display: *display, colours: colours.to_vec(), radii: radii.to_vec() };
                    let path = saved_defaults_path();
                    self.status = if save_defaults(&look, &path) { "Saved: new windows and the next start use this look.".into() } else { format!("Could not write {}.", path.display()) };
                }
                if ui.button("Restore Factory Look").on_hover_text("Every setting of this dialog, element colours and radii included, back to its default.").clicked() {
                    let look = d::default_look();
                    *display = look.display;
                    colours.copy_from_slice(&look.colours[..colours.len()]);
                    radii.copy_from_slice(&look.radii[..radii.len()]);
                    rebuild = true;
                    changes.recolour = true;
                    self.status = "Factory look restored for this window (Save as My Default keeps it).".into();
                }
            });
            if !self.status.is_empty() {
                dl::muted(ui, &self.status);
            }
            // Large structures rebuild when the slider is released.
            if rebuild {
                self.rebuild_pending = true;
            }
            if self.rebuild_pending && (structure.atoms.len() < 20000 || !ui.ctx().egui_is_using_pointer()) {
                self.rebuild_pending = false;
                changes.rebuild = true;
            }
            if dl::footer(ui, "Close", None, true) == 1 {
                self.status.clear();
                self.open = false;
            }
        });
        if !open {
            self.status.clear();
        }
        self.open = open && self.open;
        changes
    }

    fn style_tab(&mut self, ui: &mut egui::Ui, s: &mut DisplaySettings, factory: &DisplaySettings) -> bool {
        let mut rebuild = false;
        dl::section(ui, "Look");
        dl::muted(ui, "Start from a ready-made look, then fine-tune it in the other tabs.");
        let current = matching_preset(s);
        for row in PRESETS.chunks(4).enumerate() {
            ui.horizontal(|ui| {
                for (k, preset) in row.1.iter().enumerate() {
                    let i = row.0 * 4 + k;
                    if choice_button(ui, preset.name, i as i32 == current).on_hover_text(preset.description).clicked() {
                        apply_preset(s, preset);
                        rebuild = true;
                    }
                }
            });
        }
        if current >= 0 {
            let p = &PRESETS[current as usize];
            dl::muted(ui, format!("Current look: {}. {}", p.name, p.description));
        } else {
            dl::muted(ui, "Current look: custom. Pick a look above to start again from it.");
        }
        dl::section(ui, "Background");
        colour_choice(ui, "background", &mut s.custom_background, &mut s.background);
        dl::section(ui, "Cell Outline");
        if colour_choice(ui, "cell", &mut s.custom_cell_colour, &mut s.cell_colour) {
            rebuild = true;
        }
        setting(ui, "Line width", &mut s.cell_line_width, d::CELL_LINE_WIDTH, factory.cell_line_width, 1, " px", "Thickness of the cell edges in pixels, where the graphics driver draws wide lines.");
        dl::muted(ui, "Background and cell are used in exported images too.");
        dl::section(ui, "Depth Cues");
        setting(ui, "Depth cueing", &mut s.depth_cueing, d::DEPTH_CUEING, factory.depth_cueing, 2, "", "Far atoms fade into the background, so near and far are easy to tell apart: 0 for none, 1 for full fading at the back.");
        rebuild |= setting(ui, "Ambient occlusion", &mut s.ambient_occlusion, d::AMBIENT_OCCLUSION, factory.ambient_occlusion, 2, "", "Atoms darken where neighbours hide the sky from them, so crevices and buried atoms look deeper: 0 for none.\nComputed per atom when the structure changes (up to 2 million atoms).");
        dl::section(ui, "Outlines");
        setting(ui, "Outline width", &mut s.outline_width, d::OUTLINE_WIDTH, factory.outline_width, 2, "", "Dark edges around atoms and bonds, as in a drawing: 0 for none, 0.2-0.3 for a thin line.");
        if s.outline_width > 0.0 {
            dl::colour_edit3(ui, "outline colour", &mut s.outline_colour, SLIDER_WIDTH, "Outline colour");
        }
        dl::section(ui, "Selection");
        ui.horizontal(|ui| {
            rebuild |= dl::colour_edit3(ui, "selection colour", &mut s.selection_colour, SLIDER_WIDTH, "Selection ring colour");
            if ui.button("Reset").clicked() {
                s.selection_colour = factory.selection_colour;
                rebuild = true;
            }
        });
        rebuild
    }

    fn lighting_tab(&mut self, ui: &mut egui::Ui, s: &mut DisplaySettings, f: &DisplaySettings) {
        dl::section(ui, "Light Direction");
        crate::dialog_widgets::checkbox(ui, &mut s.light_follows_camera, "Light follows the camera")
            .on_hover_text("On: the light turns with the view, so the side facing you stays lit as you rotate.\nOff: the light is fixed in the scene, like a lamp in a room.");
        let azimuth_help = if s.light_follows_camera {
            "Where the light comes from around you: 0° from the camera, 90° from the right, 270° from the left."
        } else {
            "Where the light comes from around the vertical (y) axis: 0° along +x, 90° along +z."
        };
        setting(ui, "Horizontal angle", &mut s.light_azimuth, d::AZIMUTH, f.light_azimuth, 0, "°", azimuth_help);
        setting(ui, "Height", &mut s.light_elevation, d::ELEVATION, f.light_elevation, 0, "°", "How high the light stands: 0° level with the scene, 85° from almost straight above, negative from below.");
        dl::section(ui, "Brightness and Shadows");
        setting(ui, "Shadow-side light", &mut s.ambient, d::AMBIENT, f.ambient, 2, "", "Light reaching the sides of atoms facing away from the lamp. Raise it to see into dense regions.");
        setting(ui, "Shadow darkness", &mut s.shadow_strength, d::SHADOW_STRENGTH, f.shadow_strength, 2, "", "How dark the shadows atoms cast on each other are: 0 for none, 1 for black.");
        setting(ui, "Contrast", &mut s.contrast, d::CONTRAST, f.contrast, 2, "", "Stretch between dark and light: 1 for none, more for a punchier picture.");
        setting(ui, "Colour vividness", &mut s.saturation, d::SATURATION, f.saturation, 2, "", "Colour saturation: 0 for greys, 1 for the element colours as they are, more for richer colours.");
        dl::section(ui, "Surface Finish");
        setting(ui, "Highlight strength", &mut s.specular_intensity, d::SPECULAR_INTENSITY, f.specular_intensity, 2, "", "Brightness of the shiny spot on atoms and bonds: 0 for a matte look.");
        setting(ui, "Highlight sharpness", &mut s.shininess_scale, d::SHININESS_SCALE, f.shininess_scale, 2, "", "Multiplies each element's shininess: larger values give smaller, sharper highlights.");
        setting(ui, "Least sharpness", &mut s.shininess_floor, d::SHININESS_FLOOR, f.shininess_floor, 0, "", "The broadest a highlight may spread (the lowest specular exponent).");
        ui.add_space(4.0);
        if ui.button("Reset Lighting").on_hover_text("Every value on this tab back to its default.").clicked() {
            s.light_azimuth = f.light_azimuth;
            s.light_elevation = f.light_elevation;
            s.light_follows_camera = f.light_follows_camera;
            s.ambient = f.ambient;
            s.shadow_strength = f.shadow_strength;
            s.contrast = f.contrast;
            s.saturation = f.saturation;
            s.specular_intensity = f.specular_intensity;
            s.shininess_scale = f.shininess_scale;
            s.shininess_floor = f.shininess_floor;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn atoms_tab(&mut self, ui: &mut egui::Ui, s: &mut DisplaySettings, f: &DisplaySettings, radii: &mut [f32], style: AtomStyle, structure: &Structure, changes: &mut Changes) -> bool {
        let mut rebuild = false;
        dl::section(ui, "Atom Style");
        ui.horizontal(|ui| {
            for k in AtomStyle::ALL {
                let enabled = k != AtomStyle::Polyhedral || !structure.atoms.is_empty();
                if ui.add_enabled(enabled, egui::RadioButton::new(style == k, k.label())).on_hover_text(k.hint()).clicked() && style != k {
                    changes.style = Some(k);
                }
            }
        });
        dl::section(ui, "Sizes");
        rebuild |= setting(ui, "Atom size", &mut s.atom_scale, d::ATOM_SCALE, f.atom_scale, 2, "×", "Multiplies every atom's radius (on top of View > Atom Style). Smaller atoms show bonds and the inside of a structure.");
        rebuild |= setting(ui, "Bond thickness", &mut s.bond_scale, d::BOND_SCALE, f.bond_scale, 2, "×", "Multiplies the radius of the bond cylinders.");
        if !matches!(changes.style.unwrap_or(style), AtomStyle::BallAndStick | AtomStyle::Polyhedral) {
            dl::note(ui, Note::Info, "Bonds are drawn in the Ball and Stick and Polyhedral styles; choose one above to see them.");
        }
        dl::section(ui, "Bond Colour");
        ui.horizontal(|ui| {
            rebuild |= crate::dialog_widgets::radio_value(ui, &mut s.uniform_bond_colour, false, "Half of each atom's colour").changed();
            rebuild |= crate::dialog_widgets::radio_value(ui, &mut s.uniform_bond_colour, true, "One colour").changed();
        });
        if s.uniform_bond_colour {
            rebuild |= dl::colour_edit3(ui, "bond colour", &mut s.bond_colour, SLIDER_WIDTH, "Bond colour");
        }
        dl::section(ui, "Bond Rules");
        rebuild |= setting(ui, "Bond tolerance", &mut s.bond_tolerance, d::BOND_TOLERANCE, f.bond_tolerance, 2, "", "Atoms closer than this times the sum of their element radii are bonded. Raise it to bond stretched or distorted structures.");
        rebuild |= crate::dialog_widgets::checkbox(ui, &mut s.metal_bonds, "Bonds between metals").on_hover_text("Off: bonds between two metal atoms are left out, so metallic crystals show no bond network.").changed();
        ui.horizontal(|ui| {
            let hint = "No bonds are searched for structures with at least this many atoms (bond search takes time and memory).";
            field_label(ui, "Bonds up to").on_hover_text(hint);
            // ImGui::InputInt with steps of 100000.
            let mut limit = s.bond_atom_limit as f64;
            let field = SLIDER_WIDTH - 2.0 * (crate::theme::FRAME_HEIGHT + ui.spacing().item_spacing.x);
            let mut edited = dl::number_field(ui, "bond limit", &mut limit, 0, field).on_hover_text(hint).changed();
            for (label, step) in [("-", -100000.0), ("+", 100000.0)] {
                if ui.add_sized([crate::theme::FRAME_HEIGHT, crate::theme::FRAME_HEIGHT], egui::Button::new(label)).clicked() {
                    limit += step;
                    edited = true;
                }
            }
            ui.label("atoms");
            if edited {
                s.bond_atom_limit = (limit as f32).clamp(d::BOND_ATOM_LIMIT.low, d::BOND_ATOM_LIMIT.high);
                rebuild = true;
            }
        });
        if structure.atoms.len() as f32 >= s.bond_atom_limit {
            ui.colored_label(ui.visuals().warn_fg_color, format!("This structure has {} atoms: raise the limit to draw its bonds.", structure.atoms.len()));
        }
        ui.add_space(4.0);
        dl::collapsing(ui, "element radii", "Element Radii", false, |ui| {
            self.radius_editor(ui, radii, changes);
            if ui.button("Reset All Radii").clicked() {
                radii.copy_from_slice(&elements::literature_covalent_radii()[..radii.len()]);
                changes.rebuild = true;
            }
        });
        rebuild
    }

    fn colours_tab(&mut self, ui: &mut egui::Ui, colours: &mut [[f32; 3]], changes: &mut Changes) {
        ui.add_space(4.0);
        dl::muted(ui, "The colour each element is drawn with (View > Colour Atoms By > Element).");
        dl::section(ui, "Colour Scheme");
        let scheme = elements::matching_scheme(colours);
        ui.horizontal(|ui| {
            for (i, sc) in elements::colour_schemes().iter().enumerate() {
                if choice_button(ui, sc.name, i as i32 == scheme).on_hover_text(sc.description).clicked() {
                    let table = elements::scheme_colours(i as i32);
                    colours.copy_from_slice(&table[..colours.len()]);
                    changes.recolour = true;
                }
            }
        });
        if scheme < 0 {
            dl::muted(ui, "Current colours: a scheme with your own changes. Pick a scheme to replace every colour.");
        } else {
            dl::muted(ui, format!("Current colours: {}. Change single elements below.", elements::colour_schemes()[scheme as usize].name));
        }
        dl::section(ui, "One Element");
        periodic_table(ui, &mut self.colour_element, Some(colours));
        ui.add_space(4.0);
        ui.separator();
        ui.add_space(4.0);
        let z = self.colour_element;
        if (1..=118).contains(&z) {
            ui.label(format!("{} (Z = {z})", elements::symbol(z)));
            ui.horizontal(|ui| {
                if dl::colour_edit3(ui, "element colour", &mut colours[z as usize], SLIDER_WIDTH, "Colour") {
                    changes.recolour = true;
                }
                if ui.button("Reset Element").clicked() {
                    colours[z as usize] = elements::default_colour(z);
                    changes.recolour = true;
                }
            });
        } else {
            dl::muted(ui, "No element selected.");
        }
        ui.add_space(4.0);
        if ui.button("Reset All Colours").clicked() {
            let table = elements::default_colours();
            colours.copy_from_slice(&table[..colours.len()]);
            changes.recolour = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn periodic_table_positions_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for z in 1..=118 {
            let p = table_position(z);
            assert!(p.1 < 18 && p.0 < 9, "{z} at {p:?}");
            assert!(seen.insert(p), "{z} shares {p:?}");
        }
        assert_eq!(table_position(26), (3, 7));
        assert_eq!(table_position(79), (5, 10));
        assert_eq!(table_position(118), (6, 17));
        assert_eq!(table_position(57), (7, 2));
        assert_eq!(table_position(103), (8, 16));
    }
}
