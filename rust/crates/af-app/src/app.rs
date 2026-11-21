//! The AtomForge desktop (the C++ app/EditorApplication, ui/FileBrowser's main
//! menu bar and camera/Camera): tabs of open structures, the menus, the 3D
//! viewport drawn by af-render on the UI's own GPU device, orbit/pan/zoom,
//! picking and selection, undo and redo, opening, saving and image export.

use std::path::PathBuf;

use af_core::display::{default_look, load_defaults, saved_defaults_path, DisplaySettings, Look as SavedLook};
use af_core::picking::{pick_atom, pick_ray_dir, pick_ray_origin};
use af_core::structure::Structure;
use af_render::colouring::ColourMode;
use af_render::scene::{apply_default_view, build_frame_view, radians, MAX_DISTANCE, MIN_DISTANCE};
use af_render::structure::{DrawOptions, Look, StructureRenderer};
use af_science::catalog::catalog;
use eframe::egui;
use glam::Vec3;

use crate::display_settings::DisplaySettingsWindow;
use crate::atom_menu::{ElementPicker, MenuAction, Pick};
use crate::document::Document;
use crate::measure_ui::{Kind, MeasureCards};
use crate::structure_info::StructureInfoWindow;
use crate::file_browser::{FileBrowser, Purpose};
use crate::menu_table::MENU_LOCATIONS;
use crate::operation::{needs_structure, OperationDialog, Outcome};
use crate::palette::{shortcuts_window, Action, Command, Palette};
use crate::pipeline_ui::PipelinePanel;
use crate::tools::ToolWindow;
use crate::ui_script::{split, UiScript};
use crate::view::{has_property, AtomStyle, LabelContent, ViewOptions};

const DARK_BACKGROUND: [f32; 4] = [0.09, 0.11, 0.15, 1.0];
const LIGHT_BACKGROUND: [f32; 4] = [0.92, 0.93, 0.94, 1.0];
/// Settings > Interface Size.
const INTERFACE_SIZES: [(&str, f32); 5] = [("Automatic (display scaling)", 1.0), ("Compact (85%)", 0.85), ("Larger (125%)", 1.25), ("Large (150%)", 1.5), ("Extra Large (200%)", 2.0)];
const STRUCTURE_EXTENSIONS: &[&str] =
    &["cif", "xyz", "extxyz", "vasp", "poscar", "contcar", "pdb", "mol", "sdf", "mol2", "gjf", "xsf", "pwi", "in", "lmp", "data", "cell", "gen"];

/// Edit > Box Select Mode / Lasso Select Mode: right drags select.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SelectMode {
    Off,
    Box,
    Lasso,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Pending {
    None,
    CloseTab(usize),
    Quit,
}

pub struct AtomForgeApp {
    renderer: StructureRenderer,
    render_state: eframe::egui_wgpu::RenderState,
    texture: Option<egui::TextureId>,
    /// Large files being converted to atom clouds.
    cloud_conversion: crate::large_data::CloudConversion,
    /// Frame timing for the Large Dataset panel.
    frame_clock: std::time::Instant,
    frame_ms: f32,
    /// wait-cloud's deadline.
    cloud_wait: Option<std::time::Instant>,
    tabs: Vec<Document>,
    active: usize,
    display: DisplaySettings,
    /// Element colours and radii of new tabs (the saved default look).
    look: SavedLook,
    view: ViewOptions,
    settings: DisplaySettingsWindow,
    measure: MeasureCards,
    info: StructureInfoWindow,
    palette: Palette,
    picker: ElementPicker,
    pipeline: PipelinePanel,
    bulk: crate::bulk_dialog::BulkDialog,
    dislocation: crate::dislocation_dialog::DislocationDialog,
    /// The dedicated builder, edit and analysis dialogs (structure_dialog::all).
    dialogs: Vec<Box<dyn crate::structure_dialog::StructureDialog>>,
    rdf: crate::rdf_dialog::RdfDialog,
    adf: crate::adf_dialog::AdfDialog,
    cna: crate::cna_dialog::CnaDialog,
    sro: crate::sro_dialog::SroDialog,
    trajectory: crate::trajectory_dialog::TrajectoryDialog,
    overlays: crate::overlays::Overlays,
    polyhedra: crate::polyhedra_ui::PolyhedraView,
    shortcuts_open: bool,
    /// File > Autosave project every minute, and when it last ran.
    autosave: bool,
    last_autosave: std::time::Instant,
    interface_scale: f32,
    dark: bool,
    status: String,
    browser: FileBrowser,
    operation: OperationDialog,
    tools: ToolWindow,
    recent: Vec<String>,
    pending: Pending,
    about: bool,
    /// UI tour: menus held open for a screenshot (menu, then submenu).
    tour_menus: Vec<String>,
    /// The dialogs' 3D previews, drawn by the structure renderer.
    previews: crate::preview3d::Previews,
    /// The selection's context menu held open at this point (manual screenshots).
    tour_context_menu: Option<egui::Pos2>,
    /// Manual screenshots: the folder and each image's crop bounds.
    capture: Capture,
    /// Save in the close prompt: close this tab once it is saved.
    close_after_save: Option<u64>,
    /// What the docked pipeline panel asked for this frame.
    docked_requests: Option<crate::pipeline_ui::Requests>,
    /// Settings > Renderer: the saved choice, and the backend drawing now.
    renderer_choice: &'static str,
    renderer_in_use: String,
    /// Pointer travel of the current press (a click is under 4 pixels).
    drag_travel: f32,
    /// Where the pointer rests over the view, and since when (the hover tooltip).
    hover: Option<(egui::Pos2, f64)>,
    /// UI tour: a pointer position for the hover tooltip (screen points).
    tour_hover: Option<egui::Pos2>,
    select_mode: SelectMode,
    /// G: the selected atoms follow the mouse until a click or Enter (Esc cancels).
    grab: Option<af_render::grab::Grab>,
    /// The toolbar's rotation angle (degrees).
    rotate_angle: f32,
    /// The width of the toolbar's right-hand group in the last frame.
    /// The last message shown as a notice, and since when.
    notice: (String, std::time::Instant),
    /// Notification cards: message, error, when shown.
    notifications: Vec<(String, bool, std::time::Instant)>,
    /// The box (start, end) or lasso being dragged, in viewport pixels.
    region: Vec<[f32; 2]>,
    /// ATOMFORGE_SCREENSHOT: save the window to this PNG after a few frames and quit (UI checks).
    screenshot: Option<String>,
    frames: u32,
    script: Option<UiScript>,
    /// The viewport's size (pixels) and top-left corner (points) in the last frame.
    viewport_size: (i32, i32),
    viewport_origin: egui::Pos2,
}

/// A UI script failed (the process exits with code 3).
pub static SCRIPT_FAILED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn preferences_path() -> PathBuf {
    let home = std::env::var("USERPROFILE").or_else(|_| std::env::var("HOME")).unwrap_or_else(|_| ".".into());
    PathBuf::from(home).join(".atomforge_recent_files.txt")
}

impl AtomForgeApp {
    pub fn new(cc: &eframe::CreationContext<'_>, startup: Vec<String>) -> Self {
        let state = cc.wgpu_render_state.clone().expect("wgpu renderer");
        let renderer_in_use = format!("{:?}", state.adapter.get_info().backend);
        crate::theme::install_fonts(&cc.egui_ctx);
        crate::theme::apply(&cc.egui_ctx);
        // UI scripts take screenshots: windows appear at once, without fading in.
        if std::env::var_os("ATOMFORGE_UI_SCRIPT").is_some() || std::env::var_os("ATOMFORGE_SCREENSHOT").is_some() {
            cc.egui_ctx.all_styles_mut(|s| s.animation_time = 0.0);
        }
        let renderer = StructureRenderer::with_device(state.device.clone(), state.queue.clone(), state.adapter.get_info().name);
        let saved = saved_defaults_path();
        let look = if saved.exists() { load_defaults(&saved) } else { default_look() };
        let recent = std::fs::read_to_string(preferences_path()).map(|t| t.lines().map(str::to_string).filter(|l| !l.is_empty()).collect()).unwrap_or_default();
        let mut app = Self {
            renderer,
            render_state: state,
            texture: None,
            cloud_conversion: crate::large_data::CloudConversion::default(),
            frame_clock: std::time::Instant::now(),
            frame_ms: 0.0,
            cloud_wait: None,
            tabs: Vec::new(),
            active: 0,
            display: look.display,
            look,
            view: ViewOptions::default(),
            settings: DisplaySettingsWindow::default(),
            measure: MeasureCards::default(),
            info: StructureInfoWindow::default(),
            palette: Palette::default(),
            picker: ElementPicker::default(),
            pipeline: PipelinePanel::default(),
            bulk: crate::bulk_dialog::BulkDialog::default(),
            dislocation: crate::dislocation_dialog::DislocationDialog::default(),
            dialogs: crate::structure_dialog::all(),
            rdf: crate::rdf_dialog::RdfDialog::default(),
            adf: crate::adf_dialog::AdfDialog::default(),
            cna: crate::cna_dialog::CnaDialog::default(),
            sro: crate::sro_dialog::SroDialog::default(),
            trajectory: crate::trajectory_dialog::TrajectoryDialog::default(),
            overlays: crate::overlays::Overlays::default(),
            polyhedra: crate::polyhedra_ui::PolyhedraView::default(),
            shortcuts_open: false,
            // Scripted runs never write the recovery file.
            autosave: std::env::var_os("ATOMFORGE_UI_SCRIPT").is_none() && std::env::var_os("ATOMFORGE_SCREENSHOT").is_none(),
            last_autosave: std::time::Instant::now(),
            interface_scale: 1.0,
            dark: true,
            status: String::new(),
            browser: FileBrowser::default(),
            operation: OperationDialog::default(),
            tools: ToolWindow::default(),
            recent,
            pending: Pending::None,
            about: false,
            tour_menus: Vec::new(),
            previews: Default::default(),
            tour_context_menu: None,
            capture: Capture::default(),
            close_after_save: None,
            docked_requests: None,
            renderer_choice: crate::renderer_pref::saved(),
            renderer_in_use,
            drag_travel: 0.0,
            hover: None,
            tour_hover: None,
            select_mode: SelectMode::Off,
            grab: None,
            rotate_angle: 5.0,
            notice: (String::new(), std::time::Instant::now()),
            notifications: Vec::new(),
            region: Vec::new(),
            screenshot: std::env::var("ATOMFORGE_SCREENSHOT").ok(),
            frames: 0,
            script: UiScript::from_env(),
            viewport_size: (1, 1),
            viewport_origin: egui::Pos2::ZERO,
        };
        for path in startup {
            app.open_path(&path);
        }
        app
    }

    fn doc(&self) -> Option<&Document> {
        self.tabs.get(self.active)
    }

    fn doc_mut(&mut self) -> Option<&mut Document> {
        self.tabs.get_mut(self.active)
    }

    fn structure(&self) -> Structure {
        self.doc().map(|d| d.structure.clone()).unwrap_or_default()
    }

    fn remember(&mut self, path: &str) {
        self.recent.retain(|p| p != path);
        self.recent.insert(0, path.to_string());
        self.recent.truncate(5);
        let _ = std::fs::write(preferences_path(), self.recent.join("\n"));
    }

    pub fn open_path(&mut self, path: &str) {
        // Large datasets open as atom clouds (converted in the background first).
        let p = std::path::Path::new(path);
        if crate::large_data::opens_as_cloud(p) {
            if af_cloud::is_cloud_file(p) {
                self.open_cloud(path);
            } else {
                self.cloud_conversion.start(p);
                self.status = format!("Converting {} to an atom cloud", p.file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default());
            }
            return;
        }
        match af_io::loader::load_structure_from_file(path) {
            Ok(structure) => {
                let n = structure.atoms.len();
                let note = structure.load_note.clone();
                let name = std::path::Path::new(path).file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string());
                self.status = format!("Opened {name} ({n} {}).{}", if n == 1 { "atom" } else { "atoms" }, if note.is_empty() { String::new() } else { format!(" {note}") });
                // An empty untitled tab is replaced.
                if self.doc().is_some_and(|d| d.path.is_none() && d.structure.atoms.is_empty()) {
                    self.tabs.remove(self.active);
                }
                self.tabs.push(Document::new(structure, Some(path.to_string()), &self.look));
                self.active = self.tabs.len() - 1;
                self.remember(path);
            }
            Err(e) => self.status = format!("Cannot open {path}: {e}"),
        }
    }

    /// Opens an atom cloud in the active tab when it is empty, else in a new tab.
    fn open_cloud(&mut self, path: &str) -> bool {
        let mut cloud = af_render::cloud::CloudRenderer::with_device(self.render_state.device.clone(), self.render_state.queue.clone());
        if let Err(e) = cloud.open(std::path::Path::new(path)) {
            self.status = format!("Cannot open {path}: {e}");
            return false;
        }
        let total = cloud.stats().total_atoms;
        let mut doc = Document::new(Structure::default(), Some(path.to_string()), &self.look);
        doc.scene = crate::large_data::cloud_scene(&cloud);
        doc.scene_dirty = false;
        doc.pending_default_view = false;
        doc.pending_cloud_fit = true;
        doc.cloud = Some(Box::new(cloud));
        let reuse = self.doc().is_some_and(|d| d.path.is_none() && d.structure.atoms.is_empty() && d.cloud.is_none());
        if reuse {
            self.tabs[self.active] = doc;
        } else {
            self.tabs.push(doc);
            self.active = self.tabs.len() - 1;
        }
        self.remember(path);
        self.status = format!("Large dataset: {total} atoms (view only)");
        true
    }

    /// Shows a rendered texture in the viewport (the structure's or a cloud's).
    fn show_texture(&mut self, ui: &mut egui::Ui, rect: egui::Rect, view: &eframe::wgpu::TextureView) {
        let mut egui_renderer = self.render_state.renderer.write();
        let id = match self.texture {
            Some(id) => {
                egui_renderer.update_egui_texture_from_wgpu_texture(&self.render_state.device, view, eframe::wgpu::FilterMode::Linear, id);
                id
            }
            None => {
                let id = egui_renderer.register_native_texture(&self.render_state.device, view, eframe::wgpu::FilterMode::Linear);
                self.texture = Some(id);
                id
            }
        };
        drop(egui_renderer);
        ui.painter().image(id, rect, egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)), egui::Color32::WHITE);
    }

    /// A cloud tab's frame view: the structure camera, with the cloud's clip planes.
    fn cloud_frame(&self, doc: &Document, width: i32, height: i32) -> af_render::scene::FrameView {
        let mut frame = build_frame_view(&doc.camera, &doc.scene, self.view.orthographic, width, height, &self.display);
        if let Some(cloud) = &doc.cloud {
            crate::large_data::adjust_cloud_projection(&mut frame, &doc.camera, cloud.centre(), cloud.radius(), self.view.orthographic);
        }
        frame
    }

    fn open_structure(&mut self, structure: Structure, title: &str) {
        let mut doc = Document::new(structure, None, &self.look);
        doc.title = title.to_string();
        self.tabs.push(doc);
        self.active = self.tabs.len() - 1;
    }

    fn save_to(&mut self, path: &str) {
        let format = af_cli::args::detect_format(path);
        let Some(doc) = self.doc_mut() else { return };
        match af_io::loader::save_structure(&doc.structure, path, &format) {
            Ok(()) => {
                doc.mark_saved(path);
                self.status = format!("Saved {path}");
                self.remember(path);
            }
            Err(e) => self.status = format!("Cannot save {path}: {e}"),
        }
    }

    fn save(&mut self) {
        match self.doc().and_then(|d| d.path.clone()) {
            Some(path) => self.save_to(&path),
            None => self.browser.show_for(Purpose::SaveAs, "Save structure as", "structure.cif", STRUCTURE_EXTENSIONS),
        }
    }

    fn export_image(&mut self, path: &str) {
        // File > Export Image (handleImageExportIfRequested / exportStructureImage):
        // the view at 1-8 times its size, PNG, JPEG or SVG by the file's
        // extension, with or without the background and the orientation gizmo.
        use crate::image_export::{self as export, Format};
        let Some(doc) = self.tabs.get(self.active) else { return };
        let format = export::format_of(path);
        let scale = self.browser.export_scale.clamp(1, 8);
        let (w, h) = (self.viewport_size.0.max(1) * scale, self.viewport_size.1.max(1) * scale);
        let (include_background, include_gizmo) = (self.browser.export_background, self.browser.export_gizmo);
        let background = self.background();
        let clear = if include_background { background } else { [0.0; 4] };
        let options = self.draw_options();
        let result: Result<(), String> = if doc.cloud.is_some() {
            let frame = self.cloud_frame(doc, w, h);
            if format == Format::Svg {
                let v = export::ExportView { width: w, height: h, projection: frame.projection, view: frame.view };
                export::write_svg(path, &v, background, include_background, options.show_bonds, options.show_atoms, &doc.scene)
            } else {
                // The cloud at the export resolution: redrawn while its atoms stream in (at most 30 s).
                let orthographic = self.view.orthographic;
                let light = clear[0] + clear[1] + clear[2] > 1.5;
                let show_box = self.view.show_box;
                let structure = doc.structure.clone();
                let view = frame.view;
                match self.tabs.get_mut(self.active).and_then(|d| d.cloud.as_mut()) {
                    None => return,
                    Some(cloud) => {
                        cloud.settings.show_box = show_box;
                        cloud.capture(&frame.projection, &frame.view, w as u32, h as u32, orthographic, light, clear, std::time::Duration::from_secs(30)).and_then(|mut rgba| {
                            if include_gizmo {
                                export::draw_gizmo(&mut rgba, w, h, &view, background[0] > 0.5, scale, &structure);
                            }
                            export::write_raster(path, format, &rgba, w as u32, h as u32, include_background)
                        })
                    }
                }
            }
        } else {
            let frame = build_frame_view(&doc.camera, &doc.scene, self.view.orthographic, w, h, &self.display);
            if format == Format::Svg {
                let v = export::ExportView { width: w, height: h, projection: frame.projection, view: frame.view };
                export::write_svg(path, &v, background, include_background, options.show_bonds, options.show_atoms, &doc.scene)
            } else {
                let look = Look::new(&self.display, &frame, [background[0], background[1], background[2]]);
                let options = DrawOptions { background: clear, ..options };
                let structure = doc.structure.clone();
                self.renderer.upload(&doc.scene);
                let rendered = self.renderer.render(&frame, &look, &options);
                if let Some(doc) = self.tabs.get_mut(self.active) {
                    doc.upload_needed = true;
                }
                rendered.and_then(|mut rgba| {
                    if include_gizmo {
                        export::draw_gizmo(&mut rgba, w, h, &frame.view, background[0] > 0.5, scale, &structure);
                    }
                    export::write_raster(path, format, &rgba, w as u32, h as u32, include_background)
                })
            }
        };
        self.status = match result {
            Ok(()) => format!("Image exported: {path}"),
            Err(e) => format!("Failed to export image: {e}"),
        };
    }

    fn background(&self) -> [f32; 4] {
        if self.display.custom_background {
            [self.display.background[0], self.display.background[1], self.display.background[2], 1.0]
        } else if self.dark {
            DARK_BACKGROUND
        } else {
            LIGHT_BACKGROUND
        }
    }

    fn draw_options(&self) -> DrawOptions {
        DrawOptions {
            background: self.background(),
            show_bonds: self.view.bonds_drawn(),
            show_atoms: self.view.show_atoms,
            show_box: self.view.show_box && self.doc().is_some_and(|d| d.structure.has_unit_cell),
            box_colour: if self.display.custom_cell_colour { Some(self.display.cell_colour) } else { None },
            selection_colour: self.display.selection_colour,
        }
    }

    fn close_tab(&mut self, index: usize, force: bool) {
        if index >= self.tabs.len() {
            return;
        }
        if !force && self.tabs[index].is_modified() {
            self.pending = Pending::CloseTab(index);
            return;
        }
        self.tabs.remove(index);
        if self.active >= self.tabs.len() {
            self.active = self.tabs.len().saturating_sub(1);
        }
        for doc in self.tabs.iter_mut() {
            doc.upload_needed = true;
        }
    }

    fn apply_outcome(&mut self, outcome: Outcome) {
        let outcome = match outcome {
            Outcome::UpdateStep(edit) => return self.update_step(edit),
            other => other,
        };
        let step = self.operation.step().to_string();
        let Some(doc) = self.doc_mut() else {
            if let Outcome::Structure(s) = outcome {
                self.open_structure(s, "Untitled");
            }
            return;
        };
        match outcome {
            Outcome::Structure(s) => {
                if doc.structure.atoms.is_empty() && doc.path.is_none() {
                    doc.title = "Untitled".into();
                }
                doc.replace_structure(s);
                doc.apply_element_colours();
                if step.starts_with("build-") {
                    doc.pending_default_view = true;
                }
            }
            Outcome::Selection(atoms) => {
                let set: std::collections::HashSet<usize> = atoms.into_iter().collect();
                doc.selection = doc.scene.instances.atom_indices.iter().enumerate().filter(|(_, &a)| a >= 0 && set.contains(&(a as usize))).map(|(i, _)| i).collect();
                doc.scene_dirty = true;
            }
            Outcome::UpdateStep(_) => {}
        }
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        use egui::{Key, KeyboardShortcut, Modifiers};
        let pressed = |m: Modifiers, k: Key| ctx.input_mut(|i| i.consume_shortcut(&KeyboardShortcut::new(m, k)));
        if pressed(Modifiers::COMMAND | Modifiers::SHIFT, Key::S) && self.doc().is_some() {
            self.browser.show_for(Purpose::SaveAs, "Save structure as", "structure.cif", STRUCTURE_EXTENSIONS);
        }
        if pressed(Modifiers::COMMAND, Key::O) {
            self.browser.show_for(Purpose::Open, "Open structure", "", STRUCTURE_EXTENSIONS);
        }
        if pressed(Modifiers::COMMAND | Modifiers::ALT, Key::S) && self.doc().is_some() {
            self.run_action(ctx, Action::ExportImage);
        }
        if pressed(Modifiers::COMMAND, Key::S) && self.doc().is_some() {
            self.save();
        }
        if !ctx.input(|i| i.modifiers.shift) && pressed(Modifiers::COMMAND, Key::Z) {
            if let Some(d) = self.doc_mut() {
                d.undo();
            }
        }
        if pressed(Modifiers::COMMAND, Key::Y) || pressed(Modifiers::COMMAND | Modifiers::SHIFT, Key::Z) {
            if let Some(d) = self.doc_mut() {
                d.redo();
            }
        }
        if pressed(Modifiers::COMMAND, Key::K) || pressed(Modifiers::COMMAND | Modifiers::SHIFT, Key::P) {
            self.palette.show_palette();
        }
        if pressed(Modifiers::COMMAND, Key::W) {
            self.close_tab(self.active, false);
        }
        if ctx.egui_wants_keyboard_input() || self.palette.open {
            return;
        }
        if let Some(grab) = self.grab.as_mut() {
            for (key, axis) in [(Key::X, af_render::grab::Axis::X), (Key::Y, af_render::grab::Axis::Y), (Key::Z, af_render::grab::Axis::Z)] {
                if pressed(Modifiers::NONE, key) {
                    grab.axis = grab.axis.toggled(axis);
                }
            }
            if pressed(Modifiers::NONE, Key::Enter) {
                self.grab_end(true);
            } else if pressed(Modifiers::NONE, Key::Escape) {
                self.grab_end(false);
            }
            return;
        }
        if self.doc().is_some_and(|d| !d.selection.is_empty()) && pressed(Modifiers::NONE, Key::G) {
            let mouse = ctx.input(|i| i.pointer.latest_pos()).unwrap_or_default();
            let origin = self.viewport_origin;
            let scale = ctx.pixels_per_point();
            self.grab_start([(mouse.x - origin.x) * scale, (mouse.y - origin.y) * scale]);
        }
        if self.doc().is_some_and(|d| !d.selection.is_empty()) {
            if pressed(Modifiers::NONE, Key::S) {
                self.picker.open(Pick::Substitute);
            }
            if pressed(Modifiers::NONE, Key::I) && self.doc().is_some_and(|d| d.selection.len() >= 2) {
                self.picker.open(Pick::InsertMidpoint);
            }
        }
        if pressed(Modifiers::NONE, Key::F1) {
            self.shortcuts_open = true;
        }
        if pressed(Modifiers::COMMAND, Key::A) {
            self.run_action(ctx, Action::SelectAll);
        }
        if pressed(Modifiers::COMMAND, Key::D) {
            self.run_action(ctx, Action::ClearSelection);
        }
        if self.doc().is_some_and(|d| !d.structure.atoms.is_empty()) {
            for (k, key) in [Key::Num1, Key::Num2, Key::Num3].into_iter().enumerate() {
                if pressed(Modifiers::SHIFT, key) {
                    self.look_along_lattice(k);
                }
                if pressed(Modifiers::NONE, key) {
                    self.look_along_axis(k);
                }
            }
        }
        if pressed(Modifiers::NONE, Key::Delete) {
            if let Some(d) = self.doc_mut() {
                let n = d.delete_selected();
                if n > 0 {
                    self.status = format!("Deleted {n} {}", if n == 1 { "atom" } else { "atoms" });
                }
            }
        }
        if pressed(Modifiers::NONE, Key::R) {
            if let Some(d) = self.doc_mut() {
                d.pending_default_view = true;
            }
        }
        if pressed(Modifiers::NONE, Key::Escape) {
            if let Some(d) = self.doc_mut() {
                d.selection.clear();
                d.scene_dirty = true;
            }
        }
    }

    fn step_items(&mut self, ui: &mut egui::Ui, menu: &str, submenu: &str) {
        let has = self.doc().is_some_and(|d| !d.structure.atoms.is_empty());
        for &(step, m, sub, label, _) in MENU_LOCATIONS {
            if m != menu || sub != submenu {
                continue;
            }
            let Some(t) = af_pipeline::pipeline::find_modifier_type(step) else { continue };
            let enabled = has || !needs_structure(step);
            let help = t.help.split('\n').next().unwrap_or("").to_string();
            if ui.add_enabled(enabled, egui::Button::new(label)).on_hover_text(format!("{help}\nAlso a pipeline step: {step}")).clicked() {
                self.open_step_dialog(step);
                ui.close();
            }
        }
    }

    /// Measures the selected atoms (in pick order) and shows the result card.
    fn measure_selection(&mut self, kind: Kind) {
        let Some(d) = self.tabs.get(self.active) else { return };
        let m = match kind {
            Kind::Distance => af_render::measure::distance(&d.selection, &d.scene, &d.structure),
            Kind::Angle => af_render::measure::angle(&d.selection, &d.scene, &d.structure),
            Kind::Dihedral => af_render::measure::dihedral(&d.selection, &d.scene, &d.structure),
        };
        let revision = d.history.current_id();
        self.measure.set(kind, m, revision);
    }

    fn atom_info(&mut self) {
        let Some(d) = self.tabs.get(self.active) else { return };
        self.measure.atom_info = Some(af_render::measure::atom_info(&d.selection, &d.scene, &d.structure, &d.element_radii));
    }

    /// Every command of the menus, for Find a Command.
    fn commands(&self) -> Vec<Command> {
        let has = self.doc().is_some_and(|d| !d.structure.atoms.is_empty());
        let cell = self.doc().is_some_and(|d| d.structure.has_unit_cell);
        let (undo, redo) = match self.doc() {
            Some(d) => (d.history.undo_description(), d.history.redo_description()),
            None => (String::new(), String::new()),
        };
        let mut c = vec![
            Command::new("File", "Open", Action::Open, true),
            Command::new("File", "Save", Action::Save, has),
            Command::new("File", "Save As", Action::SaveAs, has),
            Command::new("File", "Export Image", Action::ExportImage, has),
            Command::new("File", "Open Project", Action::Toggle("open-project"), true),
            Command::new("File", "Save Project As", Action::Toggle("save-project"), !self.tabs.is_empty()),
            Command::new("File", "Close", Action::Close, has),
            Command::new("File", "Quit", Action::Quit, true),
            Command::new("Edit", if undo.is_empty() { "Undo".to_string() } else { format!("Undo {undo}") }, Action::Undo, !undo.is_empty()),
            Command::new("Edit", if redo.is_empty() { "Redo".to_string() } else { format!("Redo {redo}") }, Action::Redo, !redo.is_empty()),
            Command::new("Edit", "Select All Atoms", Action::SelectAll, has),
            Command::new("Edit", "Clear Selection", Action::ClearSelection, has),
            Command::new("Edit", "Structure Pipeline", Action::Toggle("pipeline"), has),
            Command::new("Analysis > Structure and Defects", "Common Neighbour Analysis", Action::Toggle("cna"), has),
            Command::new("Analysis > Structure and Defects", "Radial Distribution Function", Action::Toggle("rdf"), has),
            Command::new("Analysis > Structure and Defects", "Angular Distribution Function", Action::Toggle("adf"), has),
            Command::new("Analysis > Structure and Defects", "Short-Range Order (SRO)", Action::Toggle("sro"), has),
            Command::new("Analysis > Trajectories and Transport", "Trajectory Playback", Action::Toggle("trajectory"), true),
            Command::new("Edit", "Box Select Mode", Action::Toggle("box-select"), true),
            Command::new("Edit", "Lasso Select Mode", Action::Toggle("lasso-select"), true),
        ];
        for &(step, menu, submenu, label, _) in MENU_LOCATIONS {
            if af_pipeline::pipeline::find_modifier_type(step).is_none() {
                continue;
            }
            let path = if submenu.is_empty() { menu.to_string() } else { format!("{menu} > {submenu}") };
            c.push(Command::new(&path, label, Action::Step(step), has || !needs_structure(step)));
        }
        for d in &self.dialogs {
            let step = d.steps()[0];
            if !MENU_LOCATIONS.iter().any(|m| m.0 == step) {
                let path = if step == "dislocation-analysis" { "Analysis > Structure and Defects" } else { "Analysis > Electronic Structure" };
                c.push(Command::new(path, d.title(), Action::Step(step), true));
            }
        }
        for tool in catalog() {
            c.push(Command::new(&format!("Analysis > {}", tool.category), tool.title, Action::Tool(tool.id), true));
        }
        for (label, which) in [("Atoms", "atoms"), ("Bonds", "bonds"), ("Atom Labels", "labels"), ("Element Legend", "legend"), ("Bounding Box", "box")] {
            c.push(Command::new("View > Show", label, Action::Toggle(which), which != "box" || cell));
        }
        c.push(Command::new("View > Projection", "Perspective", Action::Projection(false), true));
        c.push(Command::new("View > Projection", "Orthographic", Action::Projection(true), true));
        c.push(Command::new("View", "Reset View", Action::ResetView, true));
        for (k, axis) in ["x", "y", "z"].iter().enumerate() {
            c.push(Command::new("View", format!("Look Along {axis}"), Action::ViewAxis(k), has));
        }
        for (k, axis) in ["a", "b", "c"].iter().enumerate() {
            c.push(Command::new("View", format!("Look Along Lattice Vector {axis}"), Action::ViewLattice(k), has && cell));
        }
        for mode in [ColourMode::Element, ColourMode::CrystalOrientation, ColourMode::GrainBoundary, ColourMode::AtomProperty] {
            let enabled = mode != ColourMode::AtomProperty || self.doc().is_some_and(|d| has_property(&d.structure));
            let label = match mode {
                ColourMode::Element => "Element",
                ColourMode::CrystalOrientation => "Crystal Orientation",
                ColourMode::GrainBoundary => "Grain Boundaries",
                ColourMode::AtomProperty => "Atom Property",
            };
            c.push(Command::new("View > Colour Atoms By", label, Action::ColourMode(mode), enabled));
        }
        for style in AtomStyle::ALL {
            c.push(Command::new("View > Atom Style", style.label(), Action::Style(style), style != AtomStyle::Polyhedral || has));
        }
        c.push(Command::new("View", "Lattice Planes Settings", Action::Toggle("planes-settings"), cell));
        c.push(Command::new("View", "Miller Directions Settings", Action::Toggle("directions-settings"), cell));
        c.push(Command::new("View", "Structure Info", Action::StructureInfo, has));
        c.push(Command::new("View", "Atom Info", Action::AtomInfo, has));
        c.push(Command::new("View", "Measure Distance", Action::Measure(Kind::Distance), has));
        c.push(Command::new("View", "Measure Angle", Action::Measure(Kind::Angle), has));
        c.push(Command::new("View", "Measure Dihedral", Action::Measure(Kind::Dihedral), has));
        c.push(Command::new("View > Theme", "Light theme", Action::Theme(false), true));
        c.push(Command::new("View > Theme", "Dark theme", Action::Theme(true), true));
        c.push(Command::new("View > Theme", "Match the system theme", Action::SystemTheme, true));
        for (label, scale) in INTERFACE_SIZES {
            c.push(Command::new("Settings", format!("Interface Size: {label}"), Action::InterfaceScale(scale), true));
        }
        c.push(Command::new("Settings", "Atomic Sizes", Action::AtomicSizes, true));
        c.push(Command::new("Settings", "Display Settings", Action::DisplaySettings, true));
        c.push(Command::new("Settings", "Polyhedral Settings", Action::Toggle("polyhedral-settings"), true));
        c.push(Command::new("View", "Polyhedral Viewer", Action::Toggle("polyhedral-viewer"), has));
        c.push(Command::new("Help", "Manual", Action::Manual, true));
        c.push(Command::new("Help", "Keyboard Shortcuts", Action::Shortcuts, true));
        c.push(Command::new("Help", "About", Action::About, true));
        c
    }

    fn run_action(&mut self, ctx: &egui::Context, action: Action) {
        match action {
            Action::Open => self.browser.show_for(Purpose::Open, "Open structure", "", STRUCTURE_EXTENSIONS),
            Action::Save => self.save(),
            Action::SaveAs => self.browser.show_for(Purpose::SaveAs, "Save structure as", "structure.cif", STRUCTURE_EXTENSIONS),
            Action::ExportImage => {
                self.browser.view_size = self.viewport_size;
                self.browser.show_for(Purpose::ExportImage, "Export image", "structure.png", &["png", "jpg", "jpeg", "svg"]);
            }
            Action::Close => self.close_tab(self.active, false),
            Action::Quit => self.pending = Pending::Quit,
            Action::Undo => {
                if let Some(d) = self.doc_mut() {
                    d.undo();
                }
            }
            Action::Redo => {
                if let Some(d) = self.doc_mut() {
                    d.redo();
                }
            }
            Action::Step(step) => self.open_step_dialog(step),
            Action::Tool(id) => self.tools.open_tool(id),
            Action::ColourMode(mode) => self.set_colour_mode(mode),
            Action::Style(style) => self.set_style(style),
            Action::Toggle("pipeline") => self.pipeline.open = true,
            Action::Toggle("rdf") => self.rdf.open = true,
            Action::Toggle("cna") => self.cna.open = true,
            Action::Toggle("sro") => self.sro.open = true,
            Action::Toggle("trajectory") => self.trajectory.open = true,
            Action::Toggle("planes-settings") => self.overlays.planes_open = true,
            Action::Toggle("open-project") => self.browser.show_for(Purpose::OpenProject, "Open AtomForge project", "project.afproject", &["afproject"]),
            Action::Toggle("save-project") => self.browser.show_for(Purpose::SaveProject, "Save AtomForge project", "project.afproject", &["afproject"]),
            Action::Toggle("polyhedral-settings") => self.polyhedra.settings_open = true,
            Action::Toggle("polyhedral-viewer") => self.polyhedra.viewer = !self.polyhedra.viewer,
            Action::Toggle("directions-settings") => self.overlays.directions_open = true,
            Action::Toggle("adf") => self.adf.open = true,
            Action::Toggle("box-select") => self.select_mode = if self.select_mode == SelectMode::Box { SelectMode::Off } else { SelectMode::Box },
            Action::Toggle("lasso-select") => self.select_mode = if self.select_mode == SelectMode::Lasso { SelectMode::Off } else { SelectMode::Lasso },
            Action::Toggle(which) => {
                let flag = match which {
                    "atoms" => &mut self.view.show_atoms,
                    "bonds" => &mut self.view.show_bonds,
                    "labels" => &mut self.view.show_labels,
                    "legend" => &mut self.view.show_legend,
                    _ => &mut self.view.show_box,
                };
                *flag = !*flag;
            }
            Action::Projection(orthographic) => self.set_projection(orthographic),
            Action::ResetView => {
                if let Some(d) = self.doc_mut() {
                    d.pending_default_view = true;
                }
            }
            Action::ViewAxis(k) => self.look_along_axis(k),
            Action::ViewLattice(k) => self.look_along_lattice(k),
            Action::StructureInfo => self.info.open = true,
            Action::AtomInfo => self.atom_info(),
            Action::Measure(kind) => self.measure_selection(kind),
            Action::Theme(dark) => ctx.set_theme(if dark { egui::Theme::Dark } else { egui::Theme::Light }),
            Action::SystemTheme => ctx.set_theme(egui::ThemePreference::System),
            Action::Manual => self.open_manual(),
            Action::InterfaceScale(scale) => {
                self.interface_scale = scale;
                ctx.set_zoom_factor(scale);
            }
            Action::AtomicSizes => self.settings.sizes_open = true,
            Action::DisplaySettings => self.settings.open = true,
            Action::SelectAll => {
                if let Some(d) = self.doc_mut() {
                    d.selection = (0..d.scene.instances.positions.len()).collect();
                    d.scene_dirty = true;
                }
            }
            Action::ClearSelection => {
                if let Some(d) = self.doc_mut() {
                    d.selection.clear();
                    d.scene_dirty = true;
                }
            }
            Action::Shortcuts => self.shortcuts_open = true,
            Action::About => self.about = true,
        }
    }

    /// 1, 2, 3: look along x, y, z (handleAxisViewRequest).
    fn look_along_axis(&mut self, k: usize) {
        let Some(d) = self.doc_mut() else { return };
        let (yaw, pitch) = [(90.0, 0.0), (0.0, 90.0), (0.0, 0.0)][k.min(2)];
        d.camera.yaw = yaw;
        d.camera.pitch = pitch;
        d.camera.roll = 0.0;
    }

    /// Shift+1, 2, 3: look along the lattice vectors a, b, c.
    fn look_along_lattice(&mut self, k: usize) {
        let Some(d) = self.doc_mut() else { return };
        if !d.structure.has_unit_cell || d.structure.atoms.is_empty() {
            return;
        }
        let v = d.structure.cell_vectors[k.min(2)];
        let (vx, vy, vz) = (v[0] as f32, v[1] as f32, v[2] as f32);
        let len = (vx * vx + vy * vy + vz * vz).sqrt();
        if len < 1e-6 {
            return;
        }
        let (nx, ny, nz) = (vx / len, vy / len, vz / len);
        d.camera.pitch = ny.clamp(-1.0, 1.0).asin() * 57.295_779_5;
        d.camera.yaw = nx.atan2(nz) * 57.295_779_5;
    }

    /// Every tab as a project tab (captureProject).
    fn capture_project(&self) -> Vec<af_electronic::project::Workspace> {
        let extra: Option<af_electronic::project::Workspace> = self.dialogs.iter().find_map(|d| d.workspace());
        self.tabs
            .iter()
            .enumerate()
            .map(|(k, d)| {
                let mut w = crate::project_io::capture(d, &self.view, &self.display, !self.dark);
                w.pipeline = self.pipeline.snapshot(&d.pipeline);
                w.has_pipeline_input = d.pipeline.active();
                w.pipeline_input = d.pipeline.input().cloned().unwrap_or_default();
                if k == self.active {
                    // The scientific-tool window goes with the active tab, as the C++ per-tab one.
                    w.science = self.tools.snapshot();
                    if let Some(e) = &extra {
                        w.volume = e.volume.clone();
                        w.reference = e.reference.clone();
                        w.surface = e.surface.clone();
                        w.table = e.table.clone();
                        w.columns = e.columns;
                        w.heading = e.heading.clone();
                        w.source_path = e.source_path.clone();
                        w.reference_path = e.reference_path.clone();
                        w.history = e.history.clone();
                        w.slice_settings = e.slice_settings.clone();
                        w.settings = e.settings.clone();
                    }
                }
                w
            })
            .collect()
    }

    /// A calculation is running in a dialog (projects wait for it).
    fn dialogs_busy(&self) -> bool {
        self.dialogs.iter().any(|d| d.busy())
    }

    fn save_project(&mut self, path: &str) {
        if self.dialogs_busy() {
            self.status = "Wait for or cancel calculations before saving the project".into();
            return;
        }
        let tabs = self.capture_project();
        self.status = match af_electronic::project::save_project(&tabs, path, af_electronic::project::FORMAT) {
            Ok(()) => {
                for d in self.tabs.iter_mut() {
                    d.mark_saved_in_place();
                }
                format!("Project saved: {path}")
            }
            Err(e) => format!("Cannot save the project: {e}"),
        };
    }

    /// Opens a project's tabs after the open ones (existing work stays open).
    pub fn open_project(&mut self, path: &str) {
        match af_electronic::project::load_project(path) {
            Ok(tabs) => {
                let n = tabs.len();
                let first = self.tabs.len();
                let mut restored_workspace = false;
                let mut restored_science = false;
                for (k, w) in tabs.into_iter().enumerate() {
                    if k == 0 {
                        crate::project_io::restore_view(&w, &mut self.view);
                    }
                    if !restored_science && !w.science.is_empty() {
                        restored_science = true;
                        self.tools.restore(&w.science);
                    }
                    // The electronic workspace: from the first tab that saved one.
                    if !restored_workspace && (!w.volume.fields.is_empty() || !w.settings.is_empty()) {
                        restored_workspace = true;
                        for d in self.dialogs.iter_mut().filter(|d| d.workspace().is_some()) {
                            if let Err(e) = d.restore_workspace(&w) {
                                self.status = e;
                            }
                        }
                    }
                    let pipeline_text = w.pipeline.clone();
                    let pipeline_input = w.has_pipeline_input.then(|| w.pipeline_input.clone());
                    let (mut doc, display) = crate::project_io::restore(w, &self.look);
                    self.pipeline.restore(&mut doc.pipeline, &pipeline_text, pipeline_input.as_ref());
                    if k == 0 {
                        if let Some(d) = display {
                            self.display = d;
                        }
                    }
                    self.tabs.push(doc);
                }
                // An empty untitled tab is replaced.
                if first > 0 && self.tabs[first - 1].path.is_none() && self.tabs[first - 1].structure.atoms.is_empty() && self.tabs[first - 1].project.is_none() {
                    self.tabs.remove(first - 1);
                }
                self.active = self.tabs.len() - n;
                for d in self.tabs.iter_mut() {
                    d.upload_needed = true;
                    d.invalidate_colours();
                }
                self.status = format!("Opened project {path} ({n} {})", if n == 1 { "tab" } else { "tabs" });
            }
            Err(e) => self.status = format!("Cannot open the project {path}: {e}"),
        }
    }

    /// Saves the recovery project once a minute in the background.
    fn autosave(&mut self) {
        if !self.autosave || self.last_autosave.elapsed().as_secs() < 60 {
            return;
        }
        self.last_autosave = std::time::Instant::now();
        if self.dialogs_busy() {
            return;
        }
        let tabs = self.capture_project();
        if !tabs.iter().any(|w| !w.structure.atoms.is_empty() || !w.volume.fields.is_empty()) {
            return;
        }
        let Some(path) = crate::project_io::autosave_path() else { return };
        std::thread::spawn(move || {
            if let Err(e) = af_electronic::project::save_project(&tabs, &path.to_string_lossy(), af_electronic::project::FORMAT) {
                eprintln!("Autosave failed: {e}");
            }
        });
    }

    /// Builds the active tab's scene now if it is out of date (UI scripts).
    fn refresh_scene(&mut self) {
        let (display, view) = (self.display, self.view);
        if let Some(d) = self.doc_mut() {
            if d.scene_dirty {
                d.rebuild_scene(&display, &view);
            }
        }
    }

    /// Help > Manual: the PDF manual in the system's viewer.
    fn open_manual(&mut self) {
        let Some(path) = find_manual() else {
            self.status = "Manual PDF not found. Restore AtomForge-manual.pdf from the application package.".into();
            return;
        };
        let text = path.to_string_lossy().into_owned();
        let opened = if cfg!(windows) {
            std::process::Command::new("cmd").args(["/C", "start", "", &text]).spawn()
        } else if cfg!(target_os = "macos") {
            std::process::Command::new("open").arg(&text).spawn()
        } else {
            std::process::Command::new("xdg-open").arg(&text).spawn()
        };
        self.status = match opened {
            Ok(_) => format!("Opened the manual ({text})"),
            Err(_) => "Unable to open the manual. Install or select a default PDF viewer.".into(),
        };
    }

    /// The toolbar under the menu bar (FileBrowser::drawToolbar): groups of
    /// controls with muted captions and thin rules between them, icons on the
    /// file, history, view and measure buttons, and the commands, pipeline,
    /// projection and theme at the right edge.
    fn toolbar(&mut self, ui: &mut egui::Ui) {
        use crate::theme::icon;
        let has = self.doc().is_some_and(|d| !d.structure.atoms.is_empty());
        let cell = has && self.doc().is_some_and(|d| d.structure.has_unit_cell);
        let (undo, redo, can_undo, can_redo) = match self.doc() {
            Some(d) => (d.history.undo_description(), d.history.redo_description(), d.history.can_undo(), d.history.can_redo()),
            None => (String::new(), String::new(), false, false),
        };
        let ctx = ui.ctx().clone();
        let palette = crate::theme::palette(ui.visuals().dark_mode);
        ui.spacing_mut().item_spacing = egui::vec2(5.0, 5.0);
        ui.spacing_mut().button_padding = egui::vec2(9.0, 4.0);
        let unit = crate::theme::FRAME_HEIGHT - 4.0;
        // A button with a tooltip (shown also when disabled).
        let tool = |ui: &mut egui::Ui, text: String, tip: &str, enabled: bool| -> bool { ui.add_enabled(enabled, egui::Button::new(text)).on_hover_text(tip).on_disabled_hover_text(tip).clicked() };
        // A mode button: highlighted while the mode is on.
        let toggle = |ui: &mut egui::Ui, text: &str, tip: &str, on: bool, enabled: bool| -> bool {
            let mut b = egui::Button::new(text);
            if on {
                b = b.fill(egui::Color32::from_rgba_unmultiplied(palette.accent.r(), palette.accent.g(), palette.accent.b(), if ui.visuals().dark_mode { 117 } else { 82 })).stroke(egui::Stroke::new(1.0, palette.accent));
            }
            ui.add_enabled(enabled, b).on_hover_text(tip).on_disabled_hover_text(tip).clicked()
        };
        // The next group on this row when it fits, else on a new row, with a thin rule between groups.
        let next_group = |ui: &mut egui::Ui, width: f32| {
            if ui.available_width() < width {
                ui.end_row();
                return;
            }
            let (rect, _) = ui.allocate_exact_size(egui::vec2(5.0, unit), egui::Sense::hover());
            let x = rect.left() + 2.0;
            ui.painter().vline(x, (rect.top() + unit * 0.18)..=(rect.top() + unit * 0.82), egui::Stroke::new(1.0, palette.separator));
        };
        let caption = |ui: &mut egui::Ui, text: &str| {
            ui.label(egui::RichText::new(text).color(palette.text_muted));
        };
        // Right-aligned: commands, the structure pipeline, the projection and
        // the theme, laid out from the right edge first; the groups before
        // them wrap in the width that is left.
        let full = ui.available_rect_before_wrap();
        let row = egui::Rect::from_min_size(full.min, egui::vec2(full.width(), unit));
        let right_rect;
        {
            let ui = &mut ui.new_child(egui::UiBuilder::new().max_rect(row).layout(egui::Layout::right_to_left(egui::Align::Center)));
            ui.spacing_mut().item_spacing = egui::vec2(5.0, 5.0);
            ui.spacing_mut().button_padding = egui::vec2(9.0, 4.0);
            let light = !self.dark;
            let theme = ui.button(if light { "Dark theme" } else { "Light theme" }).on_hover_text("Switch the colour theme (Settings > Theme can also follow the system)");
            if theme.clicked() {
                self.run_action(&ctx, Action::Theme(light));
            }
            let projection = if self.view.orthographic { "Orthographic" } else { "Perspective" };
            if tool(ui, projection.into(), "Switch between perspective and orthographic projection", true) {
                self.set_projection(!self.view.orthographic);
            }
            let pipeline_ok = has || self.doc().is_some_and(|d| d.pipeline.active());
            if tool(ui, crate::theme::with_icon(icon::PIPELINE, "Pipeline"), "Edit > Structure pipeline: non-destructive editing steps", pipeline_ok) {
                self.pipeline.open = true;
            }
            if tool(ui, crate::theme::with_icon(icon::SEARCH, "Commands"), "Find any command by name (Ctrl+K)", true) {
                self.palette.show_palette();
            }
            right_rect = ui.min_rect();
        }
        let left_width = (right_rect.left() - full.left() - 12.0).max(200.0);
        ui.allocate_ui_with_layout(egui::vec2(left_width, unit), egui::Layout::left_to_right(egui::Align::Center).with_main_wrap(true), |ui| {
            if tool(ui, crate::theme::with_icon(icon::OPEN, "Open"), "Open a structure file (Ctrl+O)", true) {
                self.run_action(&ctx, Action::Open);
            }
            if tool(ui, crate::theme::with_icon(icon::SAVE, "Save"), "Save the structure (Ctrl+S)", has) {
                self.run_action(&ctx, Action::Save);
            }
            next_group(ui, unit * 4.0);
            let tip = |what: &str, d: &str, keys: &str| if d.is_empty() { format!("{what} ({keys})") } else { format!("{what} {d} ({keys})") };
            if tool(ui, icon::UNDO.into(), &tip("Undo", &undo, "Ctrl+Z"), can_undo) {
                self.run_action(&ctx, Action::Undo);
            }
            if tool(ui, icon::REDO.into(), &tip("Redo", &redo, "Ctrl+Y"), can_redo) {
                self.run_action(&ctx, Action::Redo);
            }
            next_group(ui, unit * 9.5);
            caption(ui, "View");
            for (k, axis) in ["X", "Y", "Z"].iter().enumerate() {
                if tool(ui, axis.to_string(), &format!("Look along {} ({})", axis.to_lowercase(), k + 1), true) {
                    self.look_along_axis(k);
                }
            }
            ui.add_space(4.0);
            for (k, axis) in ["a", "b", "c"].iter().enumerate() {
                if tool(ui, axis.to_string(), &format!("Look along the lattice vector {axis} (Shift+{})", k + 1), cell) {
                    self.look_along_lattice(k);
                }
            }
            ui.add_space(4.0);
            if tool(ui, icon::HOME.into(), "Reset the view (R)", true) {
                self.run_action(&ctx, Action::ResetView);
            }
            next_group(ui, unit * 8.0);
            caption(ui, "Rotate");
            ui.add_sized(
                egui::vec2(62.0, unit),
                egui::DragValue::new(&mut self.rotate_angle).range(-360.0..=360.0).speed(1.0).custom_formatter(|v, _| format!("{v:.1}°")).custom_parser(|t| t.trim().trim_end_matches('°').trim().parse::<f64>().ok()),
            )
            .on_hover_text("Rotation angle (degrees); drag or double-click to type");
            for (k, axis) in ["X", "Y", "Z"].iter().enumerate() {
                if tool(ui, axis.to_string(), &format!("Rotate the crystal about {}", axis.to_lowercase()), has) {
                    self.rotate_crystal(k, self.rotate_angle as f64);
                }
            }
            next_group(ui, unit * 8.0);
            caption(ui, "Measure");
            if tool(ui, crate::theme::with_icon(icon::RULER, "Distance"), "Measure the distance between two atoms (pick them in the view)", has) {
                self.measure_selection(Kind::Distance);
            }
            if tool(ui, "Angle".into(), "Measure the angle between three atoms (pick them in the view)", has) {
                self.measure_selection(Kind::Angle);
            }
            next_group(ui, unit * 7.0);
            caption(ui, "Select");
            if toggle(ui, "Box", "Box selection: drag a rectangle in the view to select atoms (Ctrl adds)", self.select_mode == SelectMode::Box, has) {
                self.run_action(&ctx, Action::Toggle("box-select"));
            }
            if toggle(ui, "Lasso", "Lasso selection: draw around atoms in the view to select them (Ctrl adds)", self.select_mode == SelectMode::Lasso, has) {
                self.run_action(&ctx, Action::Toggle("lasso-select"));
            }
        });
    }

    /// The toolbar's Rotate: the crystal turns by `angle` about an axis (the camera orbits).
    fn rotate_crystal(&mut self, axis: usize, angle: f64) {
        if let Some(d) = self.doc_mut() {
            af_render::scene::rotate_crystal_view(&mut d.camera, axis, angle);
        }
    }

    /// The status bar (FileBrowser::drawStatusBar): the file, atom count,
    /// selection and cell, then the projection and colouring (both clickable).
    fn status_bar(&mut self, ui: &mut egui::Ui) {
        let accent = crate::theme::palette(ui.visuals().dark_mode).accent;
        ui.horizontal(|ui| {
            let mut first = true;
            let mut field = |ui: &mut egui::Ui, text: String, colour: Option<egui::Color32>| {
                if !first {
                    ui.weak("•");
                }
                first = false;
                match colour {
                    Some(c) => ui.colored_label(c, text),
                    None => ui.weak(text),
                };
            };
            match self.doc() {
                None => field(ui, "No structure open".into(), None),
                Some(d) if d.cloud.is_some() => {
                    let c = d.cloud.as_ref().unwrap();
                    field(ui, c.path().file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default(), Some(ui.visuals().text_color()));
                    field(ui, format!("{} atoms (large dataset, view only)", crate::large_data::count_text(c.stats().total_atoms)), None);
                }
                Some(d) if d.structure.atoms.is_empty() => field(ui, "No structure open".into(), None),
                Some(d) => {
                    let s = &d.structure;
                    let name = d.path.as_ref().map(|p| std::path::Path::new(p).file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default()).unwrap_or_else(|| "Untitled".into());
                    field(ui, name, Some(ui.visuals().text_color()));
                    field(ui, format!("{} atoms", s.atoms.len()), None);
                    let selected = d.selected_atoms();
                    if selected.len() == 1 {
                        let a = &s.atoms[selected[0]];
                        field(ui, format!("Selected {} #{} at ({:.3}, {:.3}, {:.3}) Å", a.symbol, selected[0], a.x, a.y, a.z), Some(accent));
                    } else if !selected.is_empty() {
                        let mut composition = String::new();
                        if selected.len() <= 200_000 {
                            let mut counts: Vec<(String, usize)> = Vec::new();
                            for &i in &selected {
                                let symbol = &s.atoms[i].symbol;
                                match counts.iter_mut().find(|c| &c.0 == symbol) {
                                    Some(c) => c.1 += 1,
                                    None => counts.push((symbol.clone(), 1)),
                                }
                            }
                            af_core::cxx::sort_by(&mut counts, |a, b| a.1 > b.1);
                            for (k, (symbol, n)) in counts.iter().take(3).enumerate() {
                                composition += &format!("{}{symbol} {n}", if k == 0 { " (" } else { ", " });
                            }
                            if !composition.is_empty() {
                                composition += if counts.len() > 3 { ", ...)" } else { ")" };
                            }
                        }
                        field(ui, format!("{} selected{composition}", selected.len()), Some(accent));
                    }
                    if s.has_unit_cell {
                        let l = |r: usize| (s.cell_vectors[r][0].powi(2) + s.cell_vectors[r][1].powi(2) + s.cell_vectors[r][2].powi(2)).sqrt();
                        field(ui, format!("Cell {:.3} × {:.3} × {:.3} Å", l(0), l(1), l(2)), None);
                    } else {
                        field(ui, "No periodic cell".into(), None);
                    }
                }
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let colour = format!("Colour: {}", self.view.colour_mode.label());
                let has = self.doc().is_some_and(|d| !d.structure.atoms.is_empty());
                let response = ui.add(egui::Label::new(egui::RichText::new(colour).weak()).sense(egui::Sense::click())).on_hover_text("Choose how atoms are coloured");
                egui::Popup::menu(&response).show(|ui| {
                    for (mode, label) in [(ColourMode::Element, "Element"), (ColourMode::CrystalOrientation, "Crystal Orientation"), (ColourMode::GrainBoundary, "Grain Boundaries")] {
                        if ui.add_enabled(mode == ColourMode::Element || has, egui::RadioButton::new(self.view.colour_mode == mode, label)).clicked() {
                            self.set_colour_mode(mode);
                            ui.close();
                        }
                    }
                });
                ui.weak("•");
                let projection = if self.view.orthographic { "Orthographic" } else { "Perspective" };
                if ui.add(egui::Label::new(egui::RichText::new(projection).weak()).sense(egui::Sense::click())).on_hover_text("Switch between perspective and orthographic projection").clicked() {
                    self.set_projection(!self.view.orthographic);
                }
            });
        });
    }

    /// Messages as notification cards (FileBrowser::drawNotifications):
    /// stacked upwards from above the status bar at the right, in the popup
    /// colours with a green (or red, for errors) accent bar; each lives 3.5 s
    /// and fades in and out.
    fn notices(&mut self, ctx: &egui::Context) {
        const LIFETIME: f32 = 3.5;
        if self.status != self.notice.0 {
            self.notice = (self.status.clone(), std::time::Instant::now());
            // Commands run from Find a Command are not announced (as in the C++).
            if !self.status.is_empty() && !self.status.starts_with("Ran ") {
                let lower = self.status.to_lowercase();
                let error = lower.starts_with("error") || lower.contains("cannot ") || lower.contains("failed") || lower.contains("could not") || lower.starts_with("no ");
                self.notifications.push((self.status.clone(), error, std::time::Instant::now()));
            }
        }
        self.notifications.retain(|n| n.2.elapsed().as_secs_f32() < LIFETIME);
        if self.notifications.is_empty() {
            return;
        }
        ctx.request_repaint_after(std::time::Duration::from_millis(30));
        let dark = ctx.theme() == egui::Theme::Dark;
        let p = crate::theme::palette(dark);
        let (good, bad) = if dark {
            (crate::theme::srgb(0.22, 0.90, 0.35), crate::theme::srgb(1.0, 0.45, 0.42))
        } else {
            (crate::theme::srgb(0.06, 0.50, 0.12), crate::theme::srgb(0.75, 0.15, 0.15))
        };
        let screen = ctx.content_rect();
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::unique("notifications")));
        let font = egui::FontId::proportional(crate::theme::BASE_FONT);
        let (margin, pad_x, pad_y, accent_w, gap) = (16.0, 14.0, 10.0, 4.0, 8.0);
        let mut bottom = screen.bottom() - crate::theme::FRAME_HEIGHT - margin;
        for (text, error, born) in self.notifications.iter().rev() {
            let age = born.elapsed().as_secs_f32();
            let alpha = (age / 0.15).min((LIFETIME - age) / 0.4).clamp(0.0, 1.0);
            let galley = painter.layout(text.clone(), font.clone(), p.text.gamma_multiply(alpha), 380.0);
            let size = egui::vec2(galley.size().x + accent_w + pad_x * 2.0, galley.size().y + pad_y * 2.0);
            let rect = egui::Rect::from_min_size(egui::pos2(screen.right() - margin - size.x, bottom - size.y), size);
            painter.rect_filled(rect.translate(egui::vec2(2.0, 3.0)), 8.0, egui::Color32::from_black_alpha((30.0 * alpha) as u8));
            painter.rect_filled(rect, 8.0, p.popup.gamma_multiply(alpha));
            let accent: egui::Color32 = if *error { bad } else { good }.gamma_multiply(alpha);
            painter.rect_filled(egui::Rect::from_min_size(rect.min, egui::vec2(accent_w, size.y)), egui::CornerRadius { nw: 8, sw: 8, ne: 0, se: 0 }, accent);
            painter.rect_stroke(rect, 8.0, egui::Stroke::new(1.0, p.border.gamma_multiply(alpha)), egui::StrokeKind::Inside);
            painter.galley(rect.min + egui::vec2(accent_w + pad_x, pad_y), galley, p.text);
            bottom = rect.top() - gap;
        }
    }

    /// The Structure Pipeline panel of the active tab and what it asks for.
    fn pipeline_panel(&mut self, ctx: &egui::Context) {
        let Some(d) = self.tabs.get_mut(self.active) else { return };
        let mut requests = self.pipeline.update(&mut d.pipeline, false);
        let shown = match self.docked_requests.take() {
            Some(r) if self.pipeline.docked => r,
            _ => self.pipeline.show(ctx, &mut d.pipeline, &d.structure),
        };
        if shown.output.is_some() || shown.bake.is_some() {
            requests.output = shown.output;
            requests.bake = shown.bake;
        }
        requests.show_property |= shown.show_property;
        requests.save_file = shown.save_file;
        requests.load_file = shown.load_file;
        requests.edit_in_dialog = self.pipeline.edit_request.take();
        if let Some(s) = requests.output {
            // The pipeline's output: shown, not an undo step.
            if s.atoms.len() != d.structure.atoms.len() {
                d.selection.clear();
            }
            d.structure = s;
            d.invalidate_colours();
        }
        if let Some(s) = requests.bake {
            d.replace_structure(s);
            d.invalidate_colours();
            self.status = "Baked the pipeline's output into the structure".into();
        }
        if requests.show_property {
            self.set_colour_mode(ColourMode::AtomProperty);
        }
        if let Some(i) = requests.edit_in_dialog {
            self.edit_step_in_dialog(i);
        }
        if requests.save_file {
            self.browser.show_for(Purpose::SavePipeline, "Save pipeline", "pipeline.json", &["json"]);
        }
        if requests.load_file {
            self.browser.show_for(Purpose::LoadPipeline, "Load pipeline (JSON or text)", "pipeline.json", &["json", "txt"]);
        }
    }

    /// The rectangle a manual screenshot shows (TutorialCapture's targets):
    /// the open dialogs, the menus with the menu bar, the topmost popup, the
    /// tool windows, the 3D view, or the whole window.
    fn capture_rect(&self, ctx: &egui::Context, target: &str) -> egui::Rect {
        let screen = ctx.content_rect();
        if target == "viewport" {
            let k = ctx.pixels_per_point();
            return egui::Rect::from_min_size(self.viewport_origin, egui::vec2(self.viewport_size.0 as f32 / k, self.viewport_size.1 as f32 / k));
        }
        if target == "full" {
            return screen;
        }
        let skip = [egui::Id::unique("modal backdrop"), egui::Id::unique("minimised dialogs"), egui::Id::unique("welcome panel"), egui::Id::unique("atom hover")];
        let layers: Vec<(egui::LayerId, egui::Rect)> = ctx.memory(|m| {
            m.layer_ids()
                .filter(|l| m.areas().visible_last_frame(l) || m.areas().is_visible(l))
                .filter(|l| !skip.contains(&l.id))
                .filter_map(|l| m.area_rect(l.id).map(|r| (l, r)))
                .filter(|(_, r)| r.width() > 40.0 && r.height() > 30.0 && r.intersects(screen))
                .collect()
        });
        let union = |orders: &[egui::Order]| {
            layers.iter().filter(|(l, _)| orders.contains(&l.order)).fold(egui::Rect::NOTHING, |a, (_, r)| a.union(*r))
        };
        let r = match target {
            "menus" => {
                let r = union(&[egui::Order::Foreground]);
                if r.is_positive() {
                    egui::Rect::from_min_max(egui::pos2((r.min.x - 40.0).max(0.0), 0.0), r.max)
                } else {
                    r
                }
            }
            "top" => layers.iter().rev().find(|(l, _)| l.order == egui::Order::Foreground).map(|(_, r)| *r).unwrap_or(egui::Rect::NOTHING),
            "tools" => union(&[egui::Order::Middle]),
            _ => union(&[egui::Order::Foreground, egui::Order::Middle]),
        };
        if r.is_positive() { r } else { screen }
    }

    /// Opens a menu for a UI tour screenshot (FileBrowser::forceTourMenu).
    fn force_tour_menu(&mut self, ui: &egui::Ui, name: &str, response: &egui::Response) {
        if !self.tour_menus.iter().any(|m| m == name) {
            return;
        }
        if egui::containers::menu::is_in_menu(ui) {
            // A submenu: the open item of its parent menu, shown last frame.
            use egui::containers::menu::{MenuState, SubMenu};
            let sub = SubMenu::id_from_widget_id(response.id);
            MenuState::mark_shown(ui.ctx(), sub);
            MenuState::from_ui(ui, |state, _| state.open_item = Some(sub));
            ui.ctx().request_repaint();
            return;
        }
        let id = egui::Popup::default_response_id(response);
        if !egui::Popup::is_id_open(ui.ctx(), id) {
            egui::Popup::open_id(ui.ctx(), id);
        }
    }

    /// Edit in Dialog (requestDialog / openStepDialog): the step's dedicated
    /// dialog, else the operation dialog, filled from the step and the
    /// structure entering it.
    fn edit_step_in_dialog(&mut self, index: usize) {
        let Some(d) = self.tabs.get_mut(self.active) else { return };
        if index >= d.pipeline.size() {
            return;
        }
        let m = d.pipeline.modifier(index).clone();
        if crate::menu_table::menu_path_of(&m.type_id).is_empty() {
            return;
        }
        let mut input = Structure::default();
        if d.pipeline.active() {
            let _ = d.pipeline.evaluate();
            let from_input = index == 0 || d.pipeline.stages().len() < index;
            input = if from_input { d.pipeline.input().cloned().unwrap_or_default() } else { d.pipeline.stage_output(index - 1).map(|s| s.structure.clone()).unwrap_or_default() };
        }
        let edit = crate::structure_dialog::StepEdit { tab: d.id, index, step: m.type_id.clone(), parameters: m.parameters.clone(), input };
        let step = m.type_id.as_str();
        let taken = match step {
            "build-bulk" => self.bulk.edit_step(edit.clone()),
            "build-dislocation" => self.dislocation.edit_step(edit.clone()),
            "relax" | "nvt-dynamics" | "npt-dynamics" | "standardize-cell" => self.tools.edit_step(edit.clone()),
            _ => match self.dialogs.iter_mut().find(|d| d.steps().contains(&step)) {
                Some(d) => d.edit_step(edit.clone()),
                None => false,
            },
        };
        if !taken {
            self.operation.edit_step(edit);
        }
    }

    /// Update Step (PipelineDialog::updateStep): the step's new parameters,
    /// if it is still where the dialog found it.
    fn update_step(&mut self, edit: crate::structure_dialog::StepEdit) {
        // The tab the step came from (the active one may have changed meanwhile).
        let Some(k) = self.tabs.iter().position(|d| d.id == edit.tab) else {
            self.pipeline.message = "The tab of the step edited in the dialog is closed".into();
            return;
        };
        let d = &mut self.tabs[k];
        if edit.index >= d.pipeline.size() || d.pipeline.modifier(edit.index).type_id != edit.step {
            self.pipeline.message = format!("The step edited in the dialog is no longer at position {}", edit.index + 1);
            return;
        }
        if let Ok(members) = edit.parameters.members() {
            for (name, value) in members.clone() {
                let _ = d.pipeline.set_parameter(edit.index, &name, value);
            }
        }
        if k == self.active {
            self.pipeline.select_step(edit.index);
        }
        self.pipeline.message.clear();
        self.status = format!("Updated pipeline step {}", edit.index + 1);
    }

    /// Colours the atoms by a tool's per-atom values (View > Colour Atoms By > Atom Property).
    fn colour_by_property(&mut self, p: af_science::atom_properties::AtomProperty) {
        let Some(d) = self.doc_mut() else { return };
        if p.values.len() != d.structure.atoms.len() {
            return;
        }
        d.structure.atom_property = p.values;
        d.structure.atom_property_name = p.name.clone();
        d.invalidate_colours();
        self.view.colour_mode = ColourMode::Element;
        self.set_colour_mode(ColourMode::AtomProperty);
        self.status = format!("Atoms coloured by {}", p.name);
    }

    /// A menu operation: its dedicated dialog where there is one, else the generic one.
    fn open_step_dialog(&mut self, step: &str) {
        match step {
            "build-bulk" => self.bulk.open_from_menu(),
            "build-dislocation" => self.dislocation.open_from_menu(),
            _ if self.dialogs.iter().any(|d| d.steps().contains(&step)) => {
                if let Some(d) = self.dialogs.iter_mut().find(|d| d.steps().contains(&step)) {
                    d.open(step);
                }
            }
            _ => {
                self.operation.open_step(step);
            }
        }
    }

    /// Dropped files go to an open dialog that takes them (processDroppedFiles:
    /// Custom Structure and Merge take them all first), else open in tabs.
    fn drop_files(&mut self, paths: Vec<String>) {
        if paths.is_empty() {
            return;
        }
        if let Some(d) = self.dialogs.iter_mut().find(|d| d.is_open() && d.takes_all_drops()) {
            for p in &paths {
                d.drop_file(p);
            }
            return;
        }
        for path in paths {
            let taken = self.dialogs.iter_mut().filter(|d| d.is_open()).any(|d| d.drop_file(&path));
            if !taken {
                self.open_path(&path);
            }
        }
    }

    /// A dedicated dialog's result.
    fn take_outcome(&mut self, o: crate::structure_dialog::Outcome) {
        use crate::structure_dialog::Outcome;
        match o {
            Outcome::Built(s, what) => {
                let n = s.atoms.len();
                self.take_built(s, &what);
                self.status = format!("{what}: {n} atoms");
            }
            Outcome::Edited(s, what) => {
                match self.doc_mut() {
                    Some(d) => {
                        d.replace_structure(s);
                        d.invalidate_colours();
                    }
                    None => self.open_structure(s, "Untitled"),
                }
                self.status = what;
            }
            Outcome::Select(atoms, what) => {
                if let Some(d) = self.doc_mut() {
                    let wanted: std::collections::HashSet<usize> = atoms.into_iter().collect();
                    d.selection = d.scene.instances.atom_indices.iter().enumerate().filter(|(_, &a)| a >= 0 && wanted.contains(&(a as usize))).map(|(i, _)| i).collect();
                    d.scene_dirty = true;
                }
                self.status = what;
            }
            Outcome::Status(text) => self.status = text,
            Outcome::UpdateStep(edit) => self.update_step(edit),
            Outcome::AddSteps(steps, what) => {
                let Some(d) = self.tabs.get_mut(self.active) else {
                    self.status = "Open a structure first: the pipeline belongs to a tab".into();
                    return;
                };
                self.pipeline.add_steps(&mut d.pipeline, &d.structure.clone(), steps);
                self.status = what;
            }
        }
    }

    /// Replace Main Scene from the dislocation builder: an undoable edit, with
    /// the loop shown (View > Show > Dislocation Lines).
    fn take_dislocation(&mut self, s: Structure) {
        let n = s.atoms.len();
        if s.dislocation_loop_points.len() >= 3 {
            self.overlays.show_dislocation_lines = true;
        }
        match self.doc_mut() {
            Some(d) => {
                d.replace_structure(s);
                d.invalidate_colours();
            }
            None => self.open_structure(s, "Untitled"),
        }
        self.status = format!("Inserted a dislocation: {n} atoms");
    }

    /// A builder's structure: replaces the open one (or opens a tab), framed anew.
    fn take_built(&mut self, s: Structure, what: &str) {
        let n = s.atoms.len();
        match self.doc_mut() {
            Some(d) => {
                if d.structure.atoms.is_empty() && d.path.is_none() {
                    d.title = "Untitled".into();
                }
                d.replace_structure(s);
                d.invalidate_colours();
                d.pending_default_view = true;
            }
            None => self.open_structure(s, "Untitled"),
        }
        self.status = format!("{what}: {n} atoms");
    }

    /// Trajectory playback: the frame to show, a finished load, file pickers
    /// and one GIF frame per UI frame while exporting.
    fn trajectory_frame(&mut self, ctx: &egui::Context) {
        let ev = self.trajectory.show(ctx);
        if ev.open_file {
            self.browser.show_for(Purpose::OpenTrajectory, "Open Trajectory (XYZ/extXYZ, XDATCAR, POSCAR or LAMMPS dump)", "trajectory.xyz", &["xyz", "extxyz", "lammpstrj", "dump", "vasp", "poscar", "xdatcar"]);
        }
        if ev.export_gif {
            self.browser.show_for(Purpose::ExportGif, "Export trajectory animation as GIF", "trajectory.gif", &["gif"]);
        }
        if let Some(s) = ev.show {
            self.show_trajectory_frame(s, ev.loaded, ev.colour_by_property);
        }
        if self.trajectory.exporting() {
            self.gif_step();
            ctx.request_repaint();
        }
    }

    /// A trajectory frame as the active structure (not an undo step); the first
    /// frame of a newly loaded trajectory opens a tab when none is open.
    fn show_trajectory_frame(&mut self, s: Structure, loaded: bool, property: bool) {
        if self.tabs.is_empty() {
            let title = std::path::Path::new(&self.trajectory.loaded).file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_else(|| "Trajectory".into());
            self.open_structure(s, &title);
        } else if let Some(d) = self.doc_mut() {
            if d.structure.atoms.len() != s.atoms.len() {
                d.selection.clear();
            }
            d.structure = s;
            d.invalidate_colours();
        }
        if let Some(d) = self.doc_mut() {
            if loaded {
                d.pending_default_view = true;
            }
        }
        if property {
            self.set_colour_mode(ColourMode::AtomProperty);
        }
    }

    /// Renders the next GIF frame with the current camera and adds it.
    fn gif_step(&mut self) {
        let Some(index) = self.trajectory.gif_frame() else { return };
        let s = match self.trajectory.frame_structure(index) {
            Ok(s) => s,
            Err(e) => {
                let restore = self.trajectory.finish_gif(format!("GIF export failed: {e}"));
                if let Ok(s) = self.trajectory.frame_structure(restore) {
                    self.show_trajectory_frame(s, false, false);
                }
                return;
            }
        };
        let property = self.trajectory.analyse;
        self.show_trajectory_frame(s, false, property);
        let (display, view) = (self.display, self.view);
        let (w, h) = (self.viewport_size.0.max(16), self.viewport_size.1.max(16));
        let background = self.background();
        let options = self.draw_options();
        let Some(d) = self.tabs.get_mut(self.active) else { return };
        d.rebuild_scene(&display, &view);
        let frame = build_frame_view(&d.camera, &d.scene, view.orthographic, w, h, &display);
        let look = Look::new(&display, &frame, [background[0], background[1], background[2]]);
        self.renderer.upload(&d.scene);
        d.upload_needed = true;
        let done = match self.renderer.render(&frame, &look, &options) {
            Ok(rgba) => self.trajectory.add_gif_frame(&rgba, w as usize, h as usize),
            Err(e) => Some(self.trajectory.finish_gif(format!("GIF export failed: {e}"))),
        };
        if let Some(restore) = done {
            if let Ok(s) = self.trajectory.frame_structure(restore) {
                self.show_trajectory_frame(s, false, property);
            }
            self.status = self.trajectory.gif_message.clone();
        }
    }

    /// The Atom Property Colouring panel while atoms are coloured by a property.
    fn property_colouring_panel(&mut self, ctx: &egui::Context) {
        if self.view.colour_mode != ColourMode::AtomProperty {
            return;
        }
        let Some(d) = self.tabs.get(self.active) else { return };
        if d.structure.atoms.is_empty() || d.structure.atom_property.len() != d.structure.atoms.len() {
            return;
        }
        let (name, values) = (d.structure.atom_property_name.clone(), d.structure.atom_property.clone());
        let top = crate::theme::FRAME_HEIGHT * 4.5;
        let events = crate::property_ui::property_panel(ctx, &name, &values, &mut self.view.property, top);
        if events.show_elements {
            self.set_colour_mode(ColourMode::Element);
        } else if events.changed {
            if let Some(d) = self.doc_mut() {
                d.invalidate_colours();
            }
        }
    }

    fn set_colour_mode(&mut self, mode: ColourMode) {
        if self.view.colour_mode == mode {
            return;
        }
        self.view.colour_mode = mode;
        for d in self.tabs.iter_mut() {
            d.invalidate_colours();
        }
    }

    fn set_style(&mut self, style: AtomStyle) {
        if self.view.style == style {
            return;
        }
        self.view.style = style;
        for d in self.tabs.iter_mut() {
            d.scene_dirty = true;
        }
    }

    /// A new projection frames the structure again (as the C++ does).
    fn set_projection(&mut self, orthographic: bool) {
        if self.view.orthographic != orthographic {
            self.view.orthographic = orthographic;
            if let Some(d) = self.doc_mut() {
                d.pending_default_view = true;
            }
        }
    }

    fn theme_menu(&mut self, ui: &mut egui::Ui) {
        let tour_submenu_0 = ui.menu_button("Theme", |ui| {
            let current = ui.ctx().options(|o| o.theme_preference);
            for (label, preference) in [("Light", egui::ThemePreference::Light), ("Dark", egui::ThemePreference::Dark), ("Match the system", egui::ThemePreference::System)] {
                if crate::theme::menu_radio(ui, true, current == preference, label).clicked() {
                    ui.ctx().set_theme(preference);
                    ui.close();
                }
            }
        });
        self.force_tour_menu(ui, "Theme", &tour_submenu_0.response);
    }

    /// The Display Settings and Atomic Sizes windows edit the display and the
    /// open structure's element tables.
    fn settings_windows(&mut self, ctx: &egui::Context) {
        if !self.settings.open && !self.settings.sizes_open {
            return;
        }
        let style = self.view.style;
        let mut display = self.display;
        let (mut colours, mut radii, structure) = match self.doc() {
            Some(d) => (d.element_colours.clone(), d.element_radii.clone(), d.structure.clone()),
            None => (self.look.colours.clone(), self.look.radii.clone(), Structure::default()),
        };
        let mut changes = self.settings.show(ctx, &mut display, &mut colours, &mut radii, style, &structure);
        let sizes = self.settings.show_sizes(ctx, &mut radii);
        changes.rebuild |= sizes.rebuild;
        if let Some(s) = changes.style {
            self.set_style(s);
        }
        self.display = display;
        match self.tabs.get_mut(self.active) {
            Some(d) => {
                if changes.recolour {
                    d.element_colours = colours;
                    d.invalidate_colours();
                }
                if changes.rebuild {
                    d.element_radii = radii;
                }
            }
            None => {
                self.look.colours = colours;
                self.look.radii = radii;
            }
        }
        // Sizes, bonds, occlusion and selection colour live in the scene.
        if changes.rebuild {
            for d in self.tabs.iter_mut() {
                d.scene_dirty = true;
            }
        }
    }

    /// Labels, the element legend and the axes over the 3D picture.
    fn overlays(&mut self, ui: &mut egui::Ui, rect: egui::Rect, size: (i32, i32)) {
        let light = !self.dark;
        let orthographic = self.view.orthographic;
        let display = self.display;
        let mut legend_click = None;
        if let Some(doc) = self.tabs.get(self.active) {
            if doc.structure.atoms.is_empty() {
                return;
            }
            let frame = build_frame_view(&doc.camera, &doc.scene, orthographic, size.0, size.1, &display);
            let painter = ui.painter_at(rect);
            let scale = ui.ctx().pixels_per_point();
            if self.view.show_labels {
                crate::view::draw_labels(&painter, rect, scale, &doc.scene, &doc.structure, &frame, self.view.label_content, light);
            }
            crate::view::draw_axes(&painter, rect, &doc.structure, &frame, light);
            self.measure.draw(&painter, rect, scale, &doc.scene, &frame);
            self.trajectory.draw_lines(&painter, rect, scale, &frame);
            self.overlays.draw(&painter, rect, scale, &doc.structure, &frame);
            self.tools.draw_overlay(&painter, rect, scale, &frame);
            // Cells are per drawn instance (periodic images included), as the selection.
            self.overlays.draw_cells(&painter, rect, scale, &doc.structure, &doc.selection, &frame);
            if self.polyhedra.viewer || self.view.style == AtomStyle::Polyhedral {
                self.polyhedra.draw(&painter, rect, scale, &doc.structure, &doc.element_colours, &frame);
            }
            if let Some(grab) = &self.grab {
                draw_grab(&painter, rect, scale, grab, &doc.selection, &doc.scene, &doc.structure, &frame, light);
            }
            if self.region.len() >= 2 {
                let points: Vec<egui::Pos2> = self.region.iter().map(|p| rect.min + egui::vec2(p[0] / scale, p[1] / scale)).collect();
                let fill = egui::Color32::from_rgba_unmultiplied(80, 160, 255, 45);
                let line = egui::Stroke::new(1.5, egui::Color32::from_rgba_unmultiplied(80, 160, 255, 220));
                if self.select_mode == SelectMode::Box {
                    let r = egui::Rect::from_two_pos(points[0], points[points.len() - 1]);
                    painter.rect(r, 0.0, fill, line, egui::StrokeKind::Middle);
                } else {
                    let polygon = &points[1..];
                    painter.add(egui::Shape::line(polygon.to_vec(), line));
                    if let (Some(&a), Some(&b)) = (polygon.first(), polygon.last()) {
                        painter.line_segment([b, a], egui::Stroke::new(1.5, egui::Color32::from_rgba_unmultiplied(80, 160, 255, 110)));
                    }
                }
            }
            if self.view.show_legend && doc.effective_colour_mode(self.view.colour_mode) == ColourMode::Element {
                legend_click = crate::view::draw_legend(ui, rect, &doc.structure);
            }
            // The colouring's legend at the bottom right (the property bar or codes, or the IPF triangle).
            let legend_painter = ui.painter_at(rect);
            match doc.effective_colour_mode(self.view.colour_mode) {
                ColourMode::AtomProperty if doc.structure.atom_property.len() == doc.structure.atoms.len() => {
                    let values = &doc.structure.atom_property;
                    let r = if self.view.property.auto_range { af_science::atom_properties::finite_range(values) } else { self.view.property.range };
                    crate::property_ui::draw_property_legend(&legend_painter, rect, &doc.structure.atom_property_name, r.low, r.high, values, self.dark);
                }
                ColourMode::CrystalOrientation => crate::property_ui::draw_ipf_legend(&legend_painter, rect, self.dark),
                _ => {}
            }
        }
        // A legend row selects that element's atoms (Ctrl adds them).
        if let Some((symbol, additive)) = legend_click {
            if let Some(d) = self.doc_mut() {
                let mut chosen: std::collections::HashSet<usize> = if additive { d.selected_atoms().into_iter().collect() } else { Default::default() };
                chosen.extend(d.structure.atoms.iter().enumerate().filter(|(_, a)| a.symbol == symbol).map(|(i, _)| i));
                d.selection = d.scene.instances.atom_indices.iter().enumerate().filter(|(_, &a)| a >= 0 && chosen.contains(&(a as usize))).map(|(i, _)| i).collect();
                d.scene_dirty = true;
                let n = chosen.len();
                self.status = format!("Selected {n} {}", if n == 1 { "atom" } else { "atoms" });
            }
        }
    }

    fn menus(&mut self, ui: &mut egui::Ui) {
        let has = self.doc().is_some_and(|d| !d.structure.atoms.is_empty());
        egui::MenuBar::new().config(egui::containers::menu::MenuConfig::default().style(crate::theme::menu_style)).ui(ui, |ui| {
            let tour_menu_45324 = ui.menu_button("File", |ui| {
                if ui.add(egui::Button::new("Open").shortcut_text("Ctrl+O")).clicked() {
                    self.browser.show_for(Purpose::Open, "Open structure", "", STRUCTURE_EXTENSIONS);
                    ui.close();
                }
                ui.add_enabled_ui(!self.recent.is_empty(), |ui| {
                    let tour_menu_90121 = ui.menu_button("Open Recent", |ui| {
                        for path in self.recent.clone() {
                            let exists = std::path::Path::new(&path).exists();
                            if ui.add_enabled(exists, egui::Button::new(&path)).clicked() {
                                self.open_path(&path);
                                ui.close();
                            }
                        }
                    });
                    self.force_tour_menu(ui, "Open Recent", &tour_menu_90121.response);
                });
                if ui.add_enabled(has, egui::Button::new("Close").shortcut_text("Ctrl+W")).clicked() {
                    self.close_tab(self.active, false);
                    ui.close();
                }
                if ui.add_enabled(has, egui::Button::new("Save").shortcut_text("Ctrl+S")).clicked() {
                    self.save();
                    ui.close();
                }
                if ui.add_enabled(has, egui::Button::new("Save As").shortcut_text("Ctrl+Shift+S")).clicked() {
                    self.browser.show_for(Purpose::SaveAs, "Save structure as", "structure.cif", STRUCTURE_EXTENSIONS);
                    ui.close();
                }
                if ui.add_enabled(has, egui::Button::new("Export Image")).clicked() {
                    self.browser.view_size = self.viewport_size;
                    self.browser.show_for(Purpose::ExportImage, "Export image", "structure.png", &["png", "jpg", "jpeg", "svg"]);
                    ui.close();
                }
                ui.separator();
                if ui.button("Open Project").clicked() {
                    self.browser.show_for(Purpose::OpenProject, "Open AtomForge project", "project.afproject", &["afproject"]);
                    ui.close();
                }
                if ui.add_enabled(!self.tabs.is_empty(), egui::Button::new("Save Project As")).clicked() {
                    self.browser.show_for(Purpose::SaveProject, "Save AtomForge project", "project.afproject", &["afproject"]);
                    ui.close();
                }
                crate::theme::menu_toggle(ui, true, &mut self.autosave, "Autosave project every minute");
                let recovery = crate::project_io::autosave_path().filter(|p| p.exists());
                if ui.add_enabled(recovery.is_some(), egui::Button::new("Recover Last Autosave")).clicked() {
                    if let Some(path) = recovery {
                        self.open_project(&path.to_string_lossy());
                    }
                    ui.close();
                }
                ui.separator();
                if ui.button("Quit").clicked() {
                    self.pending = Pending::Quit;
                    ui.close();
                }
            });
            self.force_tour_menu(ui, "File", &tour_menu_45324.response);
            let tour_menu_68489 = ui.menu_button("Edit", |ui| {
                let (undo, redo, can_undo, can_redo) = match self.doc() {
                    Some(d) => (d.history.undo_description(), d.history.redo_description(), d.history.can_undo(), d.history.can_redo()),
                    None => (String::new(), String::new(), false, false),
                };
                let undo_label = if undo.is_empty() { "Undo".to_string() } else { format!("Undo {undo}") };
                let redo_label = if redo.is_empty() { "Redo".to_string() } else { format!("Redo {redo}") };
                if ui.add_enabled(can_undo, egui::Button::new(undo_label).shortcut_text("Ctrl+Z")).clicked() {
                    if let Some(d) = self.doc_mut() {
                        d.undo();
                    }
                    ui.close();
                }
                if ui.add_enabled(can_redo, egui::Button::new(redo_label).shortcut_text("Ctrl+Y")).clicked() {
                    if let Some(d) = self.doc_mut() {
                        d.redo();
                    }
                    ui.close();
                }
                let steps = self.doc().map(|d| d.history.len()).unwrap_or(0);
                ui.add_enabled_ui(steps > 1, |ui| {
                    let tour_menu_74779 = ui.menu_button("History", |ui| {
                        let mut jump = None;
                        if let Some(d) = self.doc() {
                            for i in (0..d.history.len()).rev() {
                                let text = if i == 0 { "The opened structure".to_string() } else { capitalised(d.history.description(i)) };
                                if ui.selectable_label(i == d.history.index(), text).clicked() {
                                    jump = Some(i);
                                }
                            }
                        }
                        if let Some(i) = jump {
                            if let Some(d) = self.doc_mut() {
                                d.jump_to(i);
                            }
                            ui.close();
                        }
                    });
                    self.force_tour_menu(ui, "History", &tour_menu_74779.response);
                });
                ui.separator();
                if ui.add_enabled(has || self.doc().is_some_and(|d| d.pipeline.active()), egui::Button::new("Structure Pipeline")).on_hover_text("Non-destructive editing steps on a copy of the structure").clicked() {
                    self.pipeline.open = true;
                    ui.close();
                }
                let mut boxes = self.select_mode == SelectMode::Box;
                if crate::theme::menu_toggle(ui, true, &mut boxes, "Box Select Mode").on_hover_text("Right drag selects the atoms in a rectangle (Ctrl adds to the selection).").changed() {
                    self.select_mode = if boxes { SelectMode::Box } else { SelectMode::Off };
                }
                let mut lasso = self.select_mode == SelectMode::Lasso;
                if crate::theme::menu_toggle(ui, true, &mut lasso, "Lasso Select Mode").on_hover_text("Right drag draws a lasso around the atoms to select (Ctrl adds to the selection).").changed() {
                    self.select_mode = if lasso { SelectMode::Lasso } else { SelectMode::Off };
                }
                ui.separator();
                self.step_items(ui, "Edit", "");
                for submenu in ["Structure Operations", "Geometry", "Cell", "Chemistry and Defects", "Select Atoms"] {
                    ui.menu_button(submenu, |ui| self.step_items(ui, "Edit", submenu));
                }
                ui.separator();
                if ui.add_enabled(self.doc().is_some_and(|d| !d.selection.is_empty()), egui::Button::new("Delete Selected Atoms").shortcut_text("Del")).clicked() {
                    if let Some(d) = self.doc_mut() {
                        d.delete_selected();
                    }
                    ui.close();
                }
            });
            self.force_tour_menu(ui, "Edit", &tour_menu_68489.response);
            let tour_menu_66712 = ui.menu_button("Build", |ui| {
                self.step_items(ui, "Build", "");
                ui.separator();
                for submenu in ["Low-Dimensional", "Clusters and Molecules"] {
                    ui.menu_button(submenu, |ui| self.step_items(ui, "Build", submenu));
                }
            });
            self.force_tour_menu(ui, "Build", &tour_menu_66712.response);
            let tour_menu_2801 = ui.menu_button("View", |ui| {
                let property = self.doc().is_some_and(|d| has_property(&d.structure));
                let tour_submenu_1 = ui.menu_button("Colour Atoms By", |ui| {
                    for mode in [ColourMode::Element, ColourMode::CrystalOrientation, ColourMode::GrainBoundary, ColourMode::AtomProperty] {
                        let enabled = mode != ColourMode::AtomProperty || property;
                        let label = match mode {
                            ColourMode::Element => "Element",
                            ColourMode::CrystalOrientation => "Crystal Orientation",
                            ColourMode::GrainBoundary => "Grain Boundaries",
                            ColourMode::AtomProperty => "Atom Property",
                        };
                        let response = crate::theme::menu_radio(ui, enabled, self.view.colour_mode == mode, label);
                        let response = if enabled { response } else { response.on_disabled_hover_text("Available after a per-atom analysis (Analysis menu or a compute-property step).") };
                        if response.clicked() {
                            self.set_colour_mode(mode);
                            ui.close();
                        }
                    }
                });
                self.force_tour_menu(ui, "Colour Atoms By", &tour_submenu_1.response);
                let tour_submenu_2 = ui.menu_button("Atom Style", |ui| {
                    for style in AtomStyle::ALL {
                        let enabled = style != AtomStyle::Polyhedral || has;
                        let response = crate::theme::menu_radio(ui, enabled, self.view.style == style, style.label()).on_hover_text(style.hint());
                        if response.clicked() {
                            self.set_style(style);
                            ui.close();
                        }
                    }
                });
                self.force_tour_menu(ui, "Atom Style", &tour_submenu_2.response);
                ui.menu_button("Show", |ui| {
                    crate::theme::menu_toggle(ui, has, &mut self.view.show_atoms, "Atoms");
                    ui.add_enabled(!matches!(self.view.style, AtomStyle::BallAndStick | AtomStyle::Polyhedral), egui::Checkbox::new(&mut self.view.show_bonds, "Bonds"))
                        .on_disabled_hover_text("The Ball and Stick and Polyhedral styles always draw bonds.");
                    crate::theme::menu_toggle(ui, true, &mut self.view.show_labels, "Atom Labels");
                    ui.menu_button("Label Atoms With", |ui| {
                        for content in LabelContent::ALL {
                            let enabled = matches!(content, LabelContent::Element | LabelContent::Index) || property;
                            if crate::theme::menu_radio(ui, enabled, self.view.label_content == content, content.label()).clicked() {
                                self.view.label_content = content;
                                self.view.show_labels = true;
                                ui.close();
                            }
                        }
                    });
                    crate::theme::menu_toggle(ui, true, &mut self.view.show_legend, "Element Legend");
                    let cell = self.doc().is_some_and(|d| d.structure.has_unit_cell);
                    crate::theme::menu_toggle(ui, cell, &mut self.view.show_box, "Bounding Box");
                    ui.separator();
                    crate::theme::menu_toggle(ui, cell, &mut self.overlays.show_planes, "Lattice Planes");
                    crate::theme::menu_toggle(ui, cell, &mut self.overlays.show_directions, "Miller Directions");
                    crate::theme::menu_toggle(ui, has, &mut self.overlays.show_voronoi, "Voronoi Volumes");
                    crate::theme::menu_toggle(ui, has, &mut self.overlays.show_dislocation_lines, "Dislocation Lines");
                });
                ui.menu_button("Projection", |ui| {
                    if crate::theme::menu_radio(ui, true, !self.view.orthographic, "Perspective").clicked() {
                        self.set_projection(false);
                        ui.close();
                    }
                    if crate::theme::menu_radio(ui, true, self.view.orthographic, "Orthographic").clicked() {
                        self.set_projection(true);
                        ui.close();
                    }
                });
                if ui.add(egui::Button::new("Reset View").shortcut_text("R")).clicked() {
                    if let Some(d) = self.doc_mut() {
                        d.pending_default_view = true;
                    }
                    ui.close();
                }
                ui.separator();
                let cell = self.doc().is_some_and(|d| d.structure.has_unit_cell);
                if ui.add_enabled(cell, egui::Button::new("Lattice Planes Settings")).clicked() {
                    self.overlays.planes_open = true;
                    ui.close();
                }
                if ui.add_enabled(cell, egui::Button::new("Miller Directions Settings")).clicked() {
                    self.overlays.directions_open = true;
                    ui.close();
                }
                crate::theme::menu_toggle(ui, has, &mut self.polyhedra.viewer, "Polyhedral Viewer");
                ui.separator();
                if ui.add_enabled(has, egui::Button::new("Structure Info")).clicked() {
                    self.info.open = true;
                    ui.close();
                }
                if ui.add_enabled(has, egui::Button::new("Atom Info")).clicked() {
                    self.atom_info();
                    ui.close();
                }
                for (kind, label) in [(Kind::Distance, "Measure Distance"), (Kind::Angle, "Measure Angle"), (Kind::Dihedral, "Measure Dihedral")] {
                    let response = ui.add_enabled(has, egui::Button::new(label));
                    let response = if kind == Kind::Dihedral { response.on_hover_text("The torsion angle of four selected atoms, in the order they were picked") } else { response };
                    if response.clicked() {
                        self.measure_selection(kind);
                        ui.close();
                    }
                }
                ui.separator();
                self.theme_menu(ui);
                if ui.button("Dock Structure Pipeline at the Right").clicked() {
                    self.pipeline.open = true;
                    self.pipeline.docked = true;
                    ui.close();
                }
                if ui.button("Reset Window Layout").on_hover_text("Undock every tool window.").clicked() {
                    self.pipeline.docked = false;
                    ui.close();
                }
            });
            self.force_tour_menu(ui, "View", &tour_menu_2801.response);
            let tour_menu_2421 = ui.menu_button("Analysis", |ui| {
                // The C++ order and names; Simulation after a rule.
                const CATEGORIES: [(&str, &str); 7] = [
                    ("Structure and Defects", "Structure and defects"),
                    ("Property Data", "Property data"),
                    ("Trajectories and Transport", "Trajectories and transport"),
                    ("Electronic Structure", "Electronic structure"),
                    ("Reciprocal Space", "Reciprocal space"),
                    ("Thermal and Mechanical", "Thermal and mechanical"),
                    ("Simulation", "Simulation"),
                ];
                for (label, category) in CATEGORIES {
                    if category == "Simulation" {
                        ui.separator();
                    }
                    ui.menu_button(label, |ui| {
                        // The dedicated analysis dialogs first (as in the C++ menu).
                        if category == "Structure and defects" {
                            if ui.add_enabled(has, egui::Button::new("Common Neighbour Analysis")).clicked() {
                                self.cna.open = true;
                                ui.close();
                            }
                            if ui.add_enabled(has, egui::Button::new("Radial Distribution Function")).clicked() {
                                self.rdf.open = true;
                                ui.close();
                            }
                            if ui.add_enabled(has, egui::Button::new("Angular Distribution Function")).clicked() {
                                self.adf.open = true;
                                ui.close();
                            }
                            if ui.add_enabled(has, egui::Button::new("Short-Range Order (SRO)")).clicked() {
                                self.sro.open = true;
                                ui.close();
                            }
                            if ui.button("Dislocation Post-Processing").clicked() {
                                self.open_step_dialog("dislocation-analysis");
                                ui.close();
                            }
                            ui.separator();
                        }
                        if category == "Electronic structure" {
                            // Dialogs with no pipeline step (post-processing of calculations).
                            let mut chosen = None;
                            for d in &self.dialogs {
                                let step = d.steps()[0];
                                if !MENU_LOCATIONS.iter().any(|m| m.0 == step) && step != "dislocation-analysis" && ui.button(d.title()).clicked() {
                                    chosen = Some(step);
                                    ui.close();
                                }
                            }
                            if let Some(step) = chosen {
                                self.open_step_dialog(step);
                            }
                            ui.separator();
                        }
                        if category == "Trajectories and transport" {
                            if ui.button("Trajectory Playback").clicked() {
                                self.trajectory.open = true;
                                ui.close();
                            }
                            ui.separator();
                        }
                        for tool in catalog().iter().filter(|t| t.category == category) {
                            if ui.button(tool.title).on_hover_text(tool.help.split('\n').next().unwrap_or("")).clicked() {
                                self.tools.open_tool(tool.id);
                                ui.close();
                            }
                        }
                        for &(step, m, sub, label, _) in MENU_LOCATIONS {
                            if m == "Analysis" && sub == category && ui.add_enabled(has, egui::Button::new(label)).clicked() {
                                self.open_step_dialog(step);
                                ui.close();
                            }
                        }
                    });
                }
            });
            self.force_tour_menu(ui, "Analysis", &tour_menu_2421.response);
            let tour_menu_16330 = ui.menu_button("Settings", |ui| {
                let tour_submenu_3 = ui.menu_button("Interface Size", |ui| {
                    for (label, scale) in INTERFACE_SIZES {
                        if crate::theme::menu_radio(ui, true, self.interface_scale == scale, label).clicked() {
                            self.interface_scale = scale;
                            ui.ctx().set_zoom_factor(scale);
                            ui.close();
                        }
                    }
                });
                self.force_tour_menu(ui, "Interface Size", &tour_submenu_3.response);
                ui.separator();
                if ui.button("Atomic Sizes").clicked() {
                    self.settings.sizes_open = true;
                    ui.close();
                }
                if ui.button("Display Settings").clicked() {
                    self.settings.open = true;
                    ui.close();
                }
                if ui.button("Polyhedral Settings").clicked() {
                    self.polyhedra.settings_open = true;
                    ui.close();
                }
                ui.separator();
                self.theme_menu(ui);
                ui.menu_button("Renderer", |ui| {
                    // The choice applies at the next start; the one in use is shown below it.
                    for (name, label) in crate::renderer_pref::CHOICES {
                        if !crate::renderer_pref::available(name) {
                            continue;
                        }
                        if crate::theme::menu_radio(ui, true, self.renderer_choice == name, label).clicked() {
                            self.renderer_choice = name;
                            self.status = match crate::renderer_pref::save(name) {
                                Ok(()) => format!("AtomForge will draw with {label} from its next start."),
                                Err(e) => format!("Cannot save the renderer choice: {e}"),
                            };
                            ui.close();
                        }
                    }
                    ui.separator();
                    ui.weak(format!("In use: {}", self.renderer_in_use));
                });
            });
            self.force_tour_menu(ui, "Settings", &tour_menu_16330.response);
            let tour_menu_89338 = ui.menu_button("Help", |ui| {
                if ui.button("Manual").clicked() {
                    self.open_manual();
                    ui.close();
                }
                if ui.add(egui::Button::new("Find a Command").shortcut_text("Ctrl+K")).clicked() {
                    self.palette.show_palette();
                    ui.close();
                }
                if ui.add(egui::Button::new("Keyboard Shortcuts").shortcut_text("F1")).clicked() {
                    self.shortcuts_open = true;
                    ui.close();
                }
                if ui.button("About AtomForge").clicked() {
                    self.about = true;
                    ui.close();
                }
            });
            self.force_tour_menu(ui, "Help", &tour_menu_89338.response);
        });
    }

    fn tabs_bar(&mut self, ui: &mut egui::Ui) {
        // ImGui tab items: the selected tab on the surface colour with an
        // accent line along its top, the others on the button colour; an
        // unsaved document shows a dot, which turns into the close cross on hover.
        let p = crate::theme::palette(ui.visuals().dark_mode);
        let font = egui::TextStyle::Body.resolve(ui.style());
        let mut close = None;
        let mut select = None;
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 4.0;
            for (i, doc) in self.tabs.iter().enumerate() {
                let selected = i == self.active;
                let galley = ui.painter().layout_no_wrap(doc.title.clone(), font.clone(), p.text);
                let marker = 16.0;
                let size = egui::vec2(galley.size().x + 20.0 + marker + 6.0, crate::theme::FRAME_HEIGHT);
                let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click());
                let hovered = response.hovered();
                let fill = if selected { p.surface } else if hovered { p.button_hover } else { p.button };
                ui.painter().rect_filled(rect, egui::CornerRadius { nw: 6, ne: 6, sw: 0, se: 0 }, fill);
                if selected {
                    ui.painter().line_segment([rect.left_top() + egui::vec2(3.0, 1.0), rect.right_top() + egui::vec2(-3.0, 1.0)], egui::Stroke::new(2.0, p.accent));
                }
                ui.painter().galley(egui::pos2(rect.left() + 10.0, rect.center().y - galley.size().y / 2.0), galley, p.text);
                let marker_rect = egui::Rect::from_center_size(egui::pos2(rect.right() - 10.0 - marker / 2.0 + 4.0, rect.center().y), egui::vec2(marker, marker));
                let marker_response = ui.interact(marker_rect, ui.scope_id().with(("tab close", i)), egui::Sense::click());
                let show_cross = marker_response.hovered() || (hovered && !doc.is_modified()) || (selected && !doc.is_modified());
                if marker_response.hovered() {
                    ui.painter().rect_filled(marker_rect, 3.0, p.button_active);
                }
                if show_cross {
                    let c = marker_rect.center();
                    let d = 3.5;
                    let stroke = egui::Stroke::new(1.2, p.text);
                    ui.painter().line_segment([c - egui::vec2(d, d), c + egui::vec2(d, d)], stroke);
                    ui.painter().line_segment([c + egui::vec2(-d, d), c + egui::vec2(d, -d)], stroke);
                } else if doc.is_modified() {
                    ui.painter().circle_filled(marker_rect.center(), 3.0, p.text);
                }
                if marker_response.clicked() {
                    close = Some(i);
                } else if response.clicked() && !selected {
                    select = Some(i);
                }
                response.on_hover_text(doc.path.clone().unwrap_or_else(|| doc.title.clone()));
            }
        });
        if let Some(i) = select {
            self.active = i;
        }
        if let Some(i) = close {
            self.close_tab(i, false);
        }
    }

    /// Camera controls of the C++ Camera: left drag orbits, right drag pans,
    /// middle drag and the wheel zoom; a click (under 4 pixels of travel) picks.
    fn viewport_input(&mut self, response: &egui::Response, ctx: &egui::Context, size: (i32, i32), origin: egui::Pos2) {
        let scale = ctx.pixels_per_point();
        if self.grab.is_some() {
            if let Some(pos) = ctx.input(|i| i.pointer.latest_pos()) {
                self.grab_move([(pos.x - origin.x) * scale, (pos.y - origin.y) * scale], size);
            }
            if ctx.input(|i| i.pointer.button_clicked(egui::PointerButton::Primary)) || response.drag_started_by(egui::PointerButton::Primary) {
                self.grab_end(true);
            } else if ctx.input(|i| i.pointer.button_pressed(egui::PointerButton::Secondary)) {
                self.grab_end(false);
            }
            return;
        }
        let delta = response.drag_delta() * scale;
        if response.drag_started() {
            self.drag_travel = 0.0;
        }
        self.drag_travel += delta.x.abs() + delta.y.abs();
        let Some(doc) = self.tabs.get_mut(self.active) else { return };
        let camera = &mut doc.camera;
        if response.dragged_by(egui::PointerButton::Primary) {
            camera.yaw -= delta.x * 0.3;
            camera.pitch += delta.y * 0.3;
        }
        if self.select_mode != SelectMode::Off && (response.dragged_by(egui::PointerButton::Secondary) || response.drag_stopped_by(egui::PointerButton::Secondary)) {
            if let Some(pos) = ctx.input(|i| i.pointer.latest_pos()) {
                let p = [(pos.x - origin.x) * scale, (pos.y - origin.y) * scale];
                if response.drag_started_by(egui::PointerButton::Secondary) || self.region.is_empty() {
                    self.region = vec![p, p];
                } else if self.select_mode == SelectMode::Box {
                    self.region[1] = p;
                } else {
                    let last = *self.region.last().unwrap();
                    if (p[0] - last[0]).powi(2) + (p[1] - last[1]).powi(2) >= 9.0 {
                        self.region.push(p);
                    }
                }
            }
            if response.drag_stopped_by(egui::PointerButton::Secondary) {
                let additive = ctx.input(|i| i.modifiers.command);
                let region = std::mem::take(&mut self.region);
                self.select_region(&region, size, additive);
            }
            return;
        }
        let Some(doc) = self.tabs.get_mut(self.active) else { return };
        let camera = &mut doc.camera;
        if response.dragged_by(egui::PointerButton::Secondary) {
            let (yaw, pitch) = (radians(camera.yaw), radians(camera.pitch));
            let forward = Vec3::new(pitch.cos() * yaw.sin(), pitch.sin(), pitch.cos() * yaw.cos());
            let right = forward.cross(Vec3::Y).normalize();
            let up = right.cross(forward).normalize();
            let pan = 0.01 * (camera.distance / 30.0).max(1.0);
            camera.pan_offset += right * (delta.x * pan) + up * (delta.y * pan);
        }
        if response.dragged_by(egui::PointerButton::Middle) {
            camera.distance = (camera.distance + delta.y * 0.5 * 0.004 * camera.distance).clamp(MIN_DISTANCE, MAX_DISTANCE);
        }
        if response.hovered() {
            let scroll = ctx.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                // egui's points per wheel notch (50) to GLFW's notches.
                doc.scroll_velocity += scroll / 50.0 * 0.5;
            }
        }
        // Exponential decay of the accumulated zoom (time constant 1/12 s).
        if doc.scroll_velocity.abs() >= 1e-5 {
            let dt = ctx.input(|i| i.stable_dt).min(0.1);
            let retain = (-12.0 * dt).exp();
            let consumed = doc.scroll_velocity * (1.0 - retain);
            doc.scroll_velocity *= retain;
            if doc.scroll_velocity.abs() < 0.001 {
                doc.scroll_velocity = 0.0;
            }
            camera.distance = (camera.distance * (1.0 - consumed * 0.15)).clamp(MIN_DISTANCE, MAX_DISTANCE);
            ctx.request_repaint();
        }
        if response.clicked() && self.drag_travel < 4.0 {
            if let Some(pos) = response.interact_pointer_pos() {
                let (mx, my) = (f64::from((pos.x - origin.x) * scale), f64::from((pos.y - origin.y) * scale));
                // Ctrl+click adds or removes an atom (as in the C++).
                let shift = ctx.input(|i| i.modifiers.command);
                self.pick_at(mx, my, size, shift);
            }
        }
    }

    /// Starts grab mode at the mouse position (viewport pixels).
    fn grab_start(&mut self, mouse: [f32; 2]) {
        self.refresh_scene();
        let Some(d) = self.tabs.get(self.active) else { return };
        if d.selection.is_empty() {
            return;
        }
        self.grab = Some(af_render::grab::Grab::start(&d.scene, &d.structure, &d.selection, mouse));
        self.status = "Grab: move the mouse; X, Y or Z locks an axis; click or Enter confirms, Esc cancels".into();
    }

    /// The grabbed atoms follow the mouse (viewport pixels).
    fn grab_move(&mut self, mouse: [f32; 2], size: (i32, i32)) {
        let (display, orthographic) = (self.display, self.view.orthographic);
        let Some(grab) = self.grab.as_ref() else { return };
        let Some(d) = self.tabs.get_mut(self.active) else { return };
        let frame = build_frame_view(&d.camera, &d.scene, orthographic, size.0, size.1, &display);
        let delta = grab.delta(&frame, mouse);
        let selection = d.selection.clone();
        grab.apply(&mut d.scene, &mut d.structure, &selection, delta);
        d.upload_needed = true;
    }

    /// Confirms (one undoable step) or cancels the grab.
    fn grab_end(&mut self, confirm: bool) {
        let Some(grab) = self.grab.take() else { return };
        let Some(d) = self.tabs.get_mut(self.active) else { return };
        if confirm {
            grab.confirm(&mut d.structure);
            let selection = d.selection.clone();
            d.commit();
            d.selection = selection;
            self.status = "Moved the selected atoms".into();
        } else {
            grab.cancel(&mut d.scene, &mut d.structure);
            self.status = "Move cancelled".into();
        }
        d.scene_dirty = true;
        d.upload_needed = true;
    }

    /// Selects the instances inside a box (two corners) or lasso, in viewport pixels.
    fn select_region(&mut self, region: &[[f32; 2]], size: (i32, i32), additive: bool) {
        let (display, view, mode) = (self.display, self.view, self.select_mode);
        let Some(doc) = self.tabs.get_mut(self.active) else { return };
        let frame = build_frame_view(&doc.camera, &doc.scene, view.orthographic, size.0, size.1, &display);
        match mode {
            SelectMode::Box if region.len() >= 2 => af_render::region::select_box(&mut doc.selection, &doc.scene, &frame, region[0], region[region.len() - 1], additive),
            SelectMode::Lasso => {
                // The first point is repeated when the drag starts.
                let polygon: Vec<[f32; 2]> = region.iter().skip(1).copied().collect();
                af_render::region::select_lasso(&mut doc.selection, &doc.scene, &frame, &polygon, additive);
            }
            _ => return,
        }
        doc.scene_dirty = true;
        let n = doc.selected_atoms().len();
        self.status = format!("{n} {} selected", if n == 1 { "atom" } else { "atoms" });
    }

    /// Runs an item of the selected atoms' menu.
    fn atom_action(&mut self, action: MenuAction) {
        match action {
            MenuAction::Pick(pick) => self.picker.open(pick),
            MenuAction::Distance => self.measure_selection(Kind::Distance),
            MenuAction::Angle => self.measure_selection(Kind::Angle),
            MenuAction::Dihedral => self.measure_selection(Kind::Dihedral),
            MenuAction::AtomInfo => self.atom_info(),
            MenuAction::Delete => {
                if let Some(d) = self.doc_mut() {
                    let n = d.delete_selected();
                    self.status = format!("Deleted {n} {}", if n == 1 { "atom" } else { "atoms" });
                }
            }
            MenuAction::Deselect => {
                if let Some(d) = self.doc_mut() {
                    d.selection.clear();
                    d.scene_dirty = true;
                }
            }
        }
    }

    /// Substitutes the selected atoms by, or inserts at their midpoint, an atom of element `z`.
    fn apply_pick(&mut self, pick: Pick, z: i32) {
        let Some(d) = self.tabs.get_mut(self.active) else { return };
        let mut structure = d.structure.clone();
        let colours = d.element_colours.clone();
        let symbol = af_core::elements::symbol(z);
        match pick {
            Pick::Substitute => {
                let n = crate::atom_menu::substitute(&mut structure, &d.scene, &d.selection, z, &colours);
                if n == 0 {
                    return;
                }
                self.status = format!("Substituted {n} {} by {symbol}", if n == 1 { "atom" } else { "atoms" });
            }
            Pick::InsertMidpoint => {
                if !crate::atom_menu::insert_midpoint(&mut structure, &d.scene, &d.selection, z, &colours) {
                    self.status = "Select two atoms or more to insert an atom between them.".into();
                    return;
                }
                self.status = format!("Inserted {symbol} at the midpoint of the selection");
            }
        }
        let selection = d.selection.clone();
        d.structure = structure;
        d.commit();
        // The substituted atoms stay selected.
        if pick == Pick::Substitute {
            d.selection = selection;
        }
        d.invalidate_colours();
    }

    /// The atom under a resting mouse: its element, index and position (drawAtomHoverTooltip).
    fn hover_tooltip(&mut self, ui: &egui::Ui, response: &egui::Response, size: (i32, i32), origin: egui::Pos2) {
        let ctx = ui.ctx();
        let Some(pos) = self.tour_hover.or(response.hover_pos()) else {
            self.hover = None;
            return;
        };
        let now = ctx.input(|i| i.time);
        match self.hover {
            Some((p, _)) if p == pos => {}
            _ => {
                self.hover = Some((pos, now));
                ctx.request_repaint_after(std::time::Duration::from_millis(450));
                return;
            }
        }
        let since = self.hover.map(|h| h.1).unwrap_or(now);
        // Scripted runs (tour screenshots) do not wait for the pointer to rest.
        let since = if self.tour_hover.is_some() { now - 1.0 } else { since };
        if now - since < 0.4 || ctx.input(|i| i.pointer.any_down()) || self.palette.open {
            if now - since < 0.4 {
                ctx.request_repaint_after(std::time::Duration::from_millis(100));
            }
            return;
        }
        let scale = ctx.pixels_per_point();
        let (display, orthographic) = (self.display, self.view.orthographic);
        let Some(doc) = self.tabs.get(self.active) else { return };
        if doc.structure.atoms.is_empty() {
            return;
        }
        let frame = build_frame_view(&doc.camera, &doc.scene, orthographic, size.0, size.1, &display);
        let (mx, my) = (f64::from((pos.x - origin.x) * scale), f64::from((pos.y - origin.y) * scale));
        let Some(instance) = pick_instance(&frame, &doc.scene, orthographic, mx, my) else { return };
        let Some(&index) = doc.scene.instances.atom_indices.get(instance) else { return };
        let Some(atom) = (index >= 0).then(|| doc.structure.atoms.get(index as usize)).flatten() else { return };
        let p = doc.scene.instances.positions[instance];
        let text = format!("{} #{index}", atom.symbol);
        let position = format!("({:.3}, {:.3}, {:.3}) Å", p.x, p.y, p.z);
        // The C++ tooltip: the atom and its position, beside the pointer.
        egui::Area::new(egui::Id::unique("atom hover")).order(egui::Order::Tooltip).fixed_pos(pos + egui::vec2(16.0, 16.0)).interactable(false).show(ctx, |ui| {
            egui::Frame::popup(ui.style()).show(ui, |ui| {
                ui.label(text);
                ui.weak(position);
            });
        });
    }

    /// A click in the viewport at pixel (mx, my): selects the atom under it
    /// (shift toggles it in the selection; a click on nothing clears it).
    fn pick_at(&mut self, mx: f64, my: f64, size: (i32, i32), shift: bool) {
        let display = self.display;
        let orthographic = self.view.orthographic;
        let Some(doc) = self.tabs.get_mut(self.active) else { return };
        let frame = build_frame_view(&doc.camera, &doc.scene, orthographic, size.0, size.1, &display);
        let hit = pick_instance(&frame, &doc.scene, orthographic, mx, my);
        match hit {
            Some(i) => {
                if shift {
                    if let Some(k) = doc.selection.iter().position(|&s| s == i) {
                        doc.selection.remove(k);
                    } else {
                        doc.selection.push(i);
                    }
                } else {
                    doc.selection = vec![i];
                }
                let atom = doc.scene.instances.atom_indices[i] as usize;
                if let Some(a) = doc.structure.atoms.get(atom) {
                    self.status = format!("Atom {} ({}) at {:.4} {:.4} {:.4}; {} selected", atom + 1, a.symbol, a.x, a.y, a.z, doc.selection.len());
                }
            }
            None if !shift => doc.selection.clear(),
            None => {}
        }
        doc.scene_dirty = true;
    }

    fn viewport(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let available = ui.available_size();
        let scale = ctx.pixels_per_point();
        let size = (((available.x * scale) as i32).max(1), ((available.y * scale) as i32).max(1));
        let (rect, response) = ui.allocate_exact_size(available, egui::Sense::click_and_drag());
        self.viewport_size = size;
        self.viewport_origin = rect.min;
        let empty = self.doc().is_none_or(|d| d.structure.atoms.is_empty() && d.cloud.is_none());
        if empty {
            // The welcome panel over the view's background (drawWelcome).
            let b = self.background();
            // The background in sRGB, as the scene draws it.
            let srgb = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
            ui.painter().rect_filled(rect, 0.0, egui::Color32::from_rgb(srgb(b[0]), srgb(b[1]), srgb(b[2])));
            // Measurement cards belong to a structure.
            self.measure.distance = None;
            self.measure.angle = None;
            self.measure.dihedral = None;
            let recent = self.recent.clone();
            match crate::welcome::show(&ctx, rect, &recent) {
                Some(crate::welcome::WelcomeAction::OpenStructure) => self.run_action(&ctx, Action::Open),
                Some(crate::welcome::WelcomeAction::BuildCrystal) => self.bulk.open_from_menu(),
                Some(crate::welcome::WelcomeAction::OpenProject) => self.browser.show_for(Purpose::OpenProject, "Open AtomForge project", "project.afproject", &["afproject"]),
                Some(crate::welcome::WelcomeAction::ReadManual) => self.open_manual(),
                Some(crate::welcome::WelcomeAction::Open(path)) => self.open_path(&path),
                None => {}
            }
            if self.tabs.is_empty() {
                return;
            }
        }
        self.viewport_input(&response, &ctx, size, rect.min);
        self.hover_tooltip(ui, &response, size, rect.min);
        let selected = self.doc().map(|d| d.selection.len()).unwrap_or(0);
        if selected > 0 && self.select_mode == SelectMode::Off {
            let mut action = None;
            if let Some(pos) = self.tour_context_menu {
                // Scripts have no pointer to open the menu at: the same menu, held open.
                egui::Area::new(egui::Id::unique("tour context menu")).order(egui::Order::Foreground).fixed_pos(pos).show(&ctx, |ui| {
                    egui::Frame::menu(ui.style()).show(ui, |ui| {
                        crate::theme::menu_style(ui.style_mut());
                        ui.set_min_width(220.0);
                        crate::atom_menu::menu(ui, selected);
                    });
                });
            }
            response.context_menu(|ui| action = crate::atom_menu::menu(ui, selected));
            if let Some(a) = action {
                self.atom_action(a);
            }
        }
        let display = self.display;
        let view = self.view;
        let options = self.draw_options();
        let background = self.background();
        if self.tabs.get(self.active).is_some_and(|d| d.cloud.is_some()) {
            let doc = &mut self.tabs[self.active];
            if doc.scene_dirty {
                doc.rebuild_scene(&display, &view);
            }
            if doc.pending_cloud_fit {
                let radius = doc.cloud.as_ref().map(|c| c.radius()).unwrap_or(1.0);
                crate::large_data::fit_camera_to_cloud(&mut doc.camera, radius, size.0, size.1);
                doc.pending_cloud_fit = false;
            }
            let frame = self.cloud_frame(&self.tabs[self.active], size.0, size.1);
            let doc = &mut self.tabs[self.active];
            let cloud = doc.cloud.as_mut().unwrap();
            cloud.settings.show_box = view.show_box;
            let rendered = cloud.render_texture(&frame.projection, &frame.view, size.0 as u32, size.1 as u32, view.orthographic, !self.dark, background);
            if !cloud.stats().complete {
                ctx.request_repaint();
            }
            if let Some(texture) = rendered {
                self.show_texture(ui, rect, &texture);
                self.overlays(ui, rect, size);
            }
            return;
        }
        let Some(doc) = self.tabs.get_mut(self.active) else { return };
        if doc.scene_dirty {
            doc.rebuild_scene(&display, &view);
        }
        if doc.pending_default_view {
            apply_default_view(&mut doc.camera, &doc.scene, size.0, size.1, true);
            doc.pending_default_view = false;
        }
        if doc.upload_needed {
            self.renderer.upload(&doc.scene);
            doc.upload_needed = false;
            for (i, other) in self.tabs.iter_mut().enumerate() {
                if i != self.active {
                    other.upload_needed = true;
                }
            }
        }
        let Some(doc) = self.tabs.get(self.active) else { return };
        let frame = build_frame_view(&doc.camera, &doc.scene, self.view.orthographic, size.0, size.1, &display);
        let look = Look::new(&display, &frame, [background[0], background[1], background[2]]);
        match self.renderer.render_texture(&frame, &look, &options) {
            Ok(view) => {
                self.show_texture(ui, rect, &view);
                self.overlays(ui, rect, size);
            }
            Err(e) => self.status = e,
        }
    }

    fn dialogs(&mut self, ctx: &egui::Context) {
        // A finished large-data conversion opens its cloud.
        match self.cloud_conversion.finished() {
            Some(Ok(cloud)) => {
                self.open_cloud(&cloud.to_string_lossy());
            }
            Some(Err(e)) if e != "Cancelled" => self.status = format!("Cannot convert the large dataset: {e}"),
            _ => {}
        }
        self.cloud_conversion.show(ctx);
        let frame_ms = self.frame_ms;
        if let Some(cloud) = self.tabs.get_mut(self.active).and_then(|d| d.cloud.as_mut()) {
            crate::large_data::cloud_panel(ctx, cloud, frame_ms);
        }
        if let Some(path) = self.browser.show(ctx) {
            let text = path.to_string_lossy().into_owned();
            match self.browser.purpose {
                Purpose::Open => self.open_path(&text),
                Purpose::SaveAs => self.save_to(&text),
                Purpose::ExportImage => self.export_image(&text),
                Purpose::OpenProject => self.open_project(&text),
                Purpose::OpenTrajectory => self.trajectory.load(&text),
                Purpose::ExportGif => self.trajectory.start_gif(&text),
                Purpose::SavePipeline => {
                    if let Some(d) = self.tabs.get(self.active) {
                        self.pipeline.save_file(&d.pipeline, &text);
                    }
                }
                Purpose::LoadPipeline => {
                    if let Some(d) = self.tabs.get_mut(self.active) {
                        if !d.pipeline.active() && !d.structure.atoms.is_empty() {
                            d.pipeline.set_input(&d.structure);
                        }
                        self.pipeline.load_file(&mut d.pipeline, &text);
                    }
                }
                Purpose::SaveProject => self.save_project(&text),
            }
        }
        let structure = self.structure();
        let view_selection = self.doc().map(|d| d.selected_atoms()).unwrap_or_default();
        if let Some(outcome) = self.operation.show(ctx, &structure, &view_selection) {
            self.apply_outcome(outcome);
        }
        if let Some(steps) = self.operation.added.take() {
            if let Some(d) = self.tabs.get_mut(self.active) {
                let s = d.structure.clone();
                self.pipeline.add_steps(&mut d.pipeline, &s, steps);
                self.status = "Added to the structure pipeline".into();
            }
        }
        let trajectory = self.trajectory.loaded.clone();
        self.tools.set_trajectory_source(&trajectory);
        let events = self.tools.show(ctx, &structure);
        if let Some(edit) = self.tools.updated.take() {
            self.update_step(edit);
        }
        for (s, title) in events.opened {
            self.open_structure(s, &title);
        }
        if let Some(p) = events.colour {
            self.colour_by_property(p);
        }
        self.settings_windows(ctx);
        self.property_colouring_panel(ctx);
        self.pipeline_panel(ctx);
        let colours = self.doc().map(|d| d.element_colours.clone()).unwrap_or_else(|| self.look.colours.clone());
        if let Some(s) = self.bulk.show(ctx, &colours) {
            self.take_built(s, "Built a bulk crystal");
        }
        if let Some(edit) = self.bulk.updated.take() {
            self.update_step(edit);
        }
        let structure = self.structure();
        if let Some(s) = self.dislocation.show(ctx, &structure) {
            self.take_dislocation(s);
        }
        if let Some(edit) = self.dislocation.updated.take() {
            self.update_step(edit);
        }
        let mut outcomes = Vec::new();
        for d in self.dialogs.iter_mut() {
            let input = crate::structure_dialog::Input { source: &structure, colours: &colours };
            if let Some(o) = d.show(ctx, &input) {
                outcomes.push(o);
            }
        }
        for o in outcomes {
            self.take_outcome(o);
        }
        let structure = self.structure();
        self.rdf.show(ctx, &structure);
        self.adf.show(ctx, &structure);
        self.cna.show(ctx, &structure);
        self.sro.show(ctx, &structure);
        self.trajectory_frame(ctx);
        self.overlays.show(ctx, &structure);
        self.polyhedra.show_settings(ctx);
        if let Some(d) = self.tabs.get(self.active) {
            self.measure.structure_changed(d.history.current_id());
            let structure = d.structure.clone();
            self.info.show(ctx, &structure, crate::structure_info::fingerprint(&structure));
        } else {
            self.info.open = false;
        }
        // As in the C++ docking: a panel docked at the right covers the
        // result cards (they come back when it is undocked or closed).
        if !(self.pipeline.open && self.pipeline.docked && !self.tabs.is_empty()) {
            self.measure.show(ctx, f32::INFINITY);
        }
        let colours = self.doc().map(|d| d.element_colours.clone()).unwrap_or_else(|| self.look.colours.clone());
        if let Some((pick, z)) = self.picker.show(ctx, &colours) {
            self.apply_pick(pick, z);
        }
        let commands = if self.palette.open { self.commands() } else { Vec::new() };
        if let Some(action) = self.palette.show(ctx, &commands) {
            self.run_action(ctx, action);
        }
        if self.shortcuts_open {
            shortcuts_window(ctx, &mut self.shortcuts_open);
        }
        if self.about {
            // Help > About (the C++ About popup): the icon, version and a
            // scrolling list of credits, with Manual and Close at the foot.
            let mut open = true;
            let mut manual = false;
            let renderer = self.renderer.adapter_name.clone();
            crate::theme::modal_window("About").open(&mut open).resizable(false).fixed_size([828.0, 560.0]).show(ctx, |ui| {
                let p = crate::theme::palette(ui.visuals().dark_mode);
                ui.horizontal(|ui| {
                    let (rect, _) = ui.allocate_exact_size(egui::vec2(60.0, 52.0), egui::Sense::hover());
                    let o = rect.min;
                    let (p0, p1, p2) = (o.x + 6.0, o.x + 26.0, o.x + 46.0);
                    let (q0, q1, q2) = (o.y + 6.0, o.y + 26.0, o.y + 46.0);
                    let box_col = egui::Color32::from_rgba_unmultiplied(46, 94, 119, 166);
                    let stroke = egui::Stroke::new(3.0, box_col);
                    let painter = ui.painter();
                    painter.line_segment([egui::pos2(p0, q0), egui::pos2(p2, q0)], stroke);
                    painter.line_segment([egui::pos2(p0, q2), egui::pos2(p2, q2)], stroke);
                    painter.line_segment([egui::pos2(p0, q0), egui::pos2(p0, q2)], stroke);
                    painter.line_segment([egui::pos2(p2, q0), egui::pos2(p2, q2)], stroke);
                    let node = egui::Color32::from_rgb(63, 110, 134);
                    for (x, y) in [(p0, q0), (p2, q0), (p0, q2), (p2, q2)] {
                        painter.circle_filled(egui::pos2(x, y), 8.5, node);
                    }
                    painter.circle_filled(egui::pos2(p1, q1), 10.2, egui::Color32::from_rgb(18, 63, 87));
                    ui.vertical(|ui| {
                        ui.label(egui::RichText::new("AtomForge").font(crate::theme::heading(1.0)));
                        ui.label(
                            egui::RichText::new(format!("Version {} - Atomistic structure builder, viewer, and editor with periodic-cell and finite-geometry workflows", env!("CARGO_PKG_VERSION")))
                                .color(p.text_muted),
                        );
                    });
                });
                ui.separator();
                let footer = crate::theme::FRAME_HEIGHT + 20.0;
                egui::ScrollArea::vertical().max_height(ui.available_height() - footer).auto_shrink([false, false]).show(ui, |ui| {
                    let bullet = |ui: &mut egui::Ui, text: &str| {
                        ui.horizontal_wrapped(|ui| {
                            ui.add_space(10.0);
                            ui.label("•");
                            ui.label(text);
                        });
                    };
                    crate::theme::section(ui, "Creator");
                    bullet(ui, "Albert Linda");
                    crate::theme::section(ui, "Core Libraries");
                    bullet(ui, "egui and eframe: immediate-mode GUI framework (by Emil Ernerfeldt)");
                    bullet(ui, &format!("wgpu: GPU rendering through Vulkan, Direct3D 12, Metal or OpenGL (in use: {renderer})"));
                    bullet(ui, "winit: window and input management");
                    bullet(ui, "glam: mathematics library for transformations");
                    crate::theme::section(ui, "Chemistry Libraries");
                    bullet(ui, "AtomForge itself reads and writes CIF, POSCAR, XYZ and extended XYZ, PDB, SDF, MOL, MOL2, PWI, GJF and LAMMPS data");
                    bullet(ui, "spglib: crystallographic space group and symmetry operations");
                    crate::theme::section(ui, "Data & References");
                    bullet(ui, "Periodic table element data (atomic numbers, masses, default colours)");
                    bullet(ui, "Cordero et al. (Dalton Trans. 2008): covalent radii, the default atom sizes");
                    bullet(ui, "Standard crystallographic conventions for cell vectors and positions");
                    crate::theme::section(ui, "Algorithm References");
                    bullet(ui, "CSL Grain Boundary Builder: Jianli Cheng, Jian Luo, and Kesong Yang, Aimsgb: An Algorithm and Open-Source Python Library to Generate Periodic Grain Boundary Structures, Comput. Mater. Sci. 155, 92-103 (2018).");
                    bullet(ui, "Nanocrystal Builder (Wulff construction): Georgios D. Barmparis, Zbigniew Lodziana, Nuria Lopez, and Ioannis N. Remediakis, Nanoparticle shapes by using Wulff constructions and first-principles calculations, Beilstein J. Nanotechnol. 6, 361-368 (2015).");
                    bullet(ui, "Dislocation Builder and Post-Processing: native reimplementations of the physical models and algorithms used by the BABEL dislocation toolkit (Douin, Clouet et al.). Isotropic and anisotropic (Stroh sextic formalism) displacement fields: A. N. Stroh, Philos. Mag. 3, 625 (1958); J. Math. Phys. 41, 77 (1962); D. M. Barnett and J. Lothe, Phys. Norv. 7, 13 (1973); J. P. Hirth and J. Lothe, Theory of Dislocations, 2nd ed. (1982), ch. 13. Nye tensor: J. F. Nye, Acta Metall. 1, 153 (1953); Hartley and Mishin, Acta Mater. 53, 1313 (2005), with the index-order correction documented in BABEL's own source (attributed there to David Olmsted). Vitek differential-displacement map: V. Vitek, Philos. Mag. 18, 773 (1968).");
                    crate::theme::section(ui, "Features");
                    for f in [
                        "Interactive 3D visualization with real-time rendering",
                        "Full undo/redo support for all editing operations",
                        "Multi-format structure loading and exporting",
                        "Crystallographic analysis: space groups, symmetry, cell metrics",
                        "Bulk crystal, CSL grain boundary, nanocrystal, polycrystal, substitutional solid solution, and custom mesh-filled structure builders",
                        "Wulff-construction nanocrystal builder with symmetry-expanded facet preview for periodic reference crystals, including non-cubic systems",
                        "Interactive Merge Structures workflow with per-structure 3D gizmo transforms",
                        "Customisable atom colours, sizes, and materials",
                        "Editable Wulff facet colours with dedicated 3D preview and high-contrast face boundaries",
                        "Crystal-orientation IPF colouring with in-view legend and saved sidecar metadata",
                        "Distance and angle measurement tools",
                        "Light and dark UI themes with matching overlay rendering",
                        "High-resolution image export with adjustable scale and HiDPI-aware interface scaling",
                        "Supercell transformation with periodic boundary visualization",
                        "Headless CLI builder workflows via --build modes (bulk, gb, poly, nano, amorphous, sss, custom)",
                    ] {
                        bullet(ui, f);
                    }
                });
                ui.separator();
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if crate::theme::primary_button(ui, "Close", Some(egui::vec2(90.0, crate::theme::FRAME_HEIGHT))).clicked() {
                        self.about = false;
                    }
                    if crate::theme::secondary_button(ui, "Manual", Some(egui::vec2(90.0, crate::theme::FRAME_HEIGHT))).clicked() {
                        manual = true;
                    }
                });
            });
            if !open {
                self.about = false;
            }
            if manual {
                self.open_manual();
            }
        }
        match self.pending {
            Pending::None => {}
            Pending::CloseTab(i) if i >= self.tabs.len() => self.pending = Pending::None,
            _ => {
                let modified: Vec<String> = match self.pending {
                    Pending::CloseTab(i) => vec![self.tabs[i].title.clone()],
                    _ => self.tabs.iter().filter(|d| d.is_modified()).map(|d| d.title.clone()).collect(),
                };
                if modified.is_empty() {
                    if self.pending == Pending::Quit {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                    self.pending = Pending::None;
                    return;
                }
                // The C++ prompts: closing one tab (Don't Save / Cancel / Save)
                // or quitting with unsaved tabs (Quit Without Saving / Cancel).
                let closing = match self.pending {
                    Pending::CloseTab(i) => Some(i),
                    _ => None,
                };
                let title = if closing.is_some() { "Unsaved Changes##close" } else { "Unsaved Changes##quit" };
                let mut choice = 0;
                crate::theme::modal_window(title.split("##").next().unwrap_or(title)).id(egui::Id::unique(title)).resizable(false).show(ctx, |ui| {
                    // A fixed width: the right-aligned buttons would otherwise
                    // stretch the window across the screen.
                    ui.set_width(if closing.is_some() { 460.0 } else { 520.0 });
                    let p = crate::theme::palette(ui.visuals().dark_mode);
                    match closing {
                        Some(_) => {
                            ui.label(format!("Save the changes to {} before closing it?", modified[0]));
                            ui.label(egui::RichText::new("Changes that are not saved are lost.").color(p.text_muted));
                        }
                        None => {
                            let n = modified.len();
                            ui.label(format!("{} unsaved changes: {}", if n == 1 { "One structure has".to_string() } else { format!("{n} structures have") }, modified.join(", ")));
                            ui.label(egui::RichText::new("Quit anyway and lose them, or cancel and save them first.").color(p.text_muted));
                        }
                    }
                    ui.add_space(4.0);
                    ui.separator();
                    let width = |ui: &egui::Ui, label: &str| {
                        let text = ui.painter().layout_no_wrap(label.to_string(), egui::TextStyle::Button.resolve(ui.style()), p.text).size().x;
                        egui::vec2((crate::theme::BASE_FONT * 5.5 * 18.0 / 13.5).max(text + 20.0), crate::theme::FRAME_HEIGHT)
                    };
                    ui.allocate_ui_with_layout(egui::vec2(ui.available_width(), crate::theme::FRAME_HEIGHT), egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if closing.is_some() && crate::theme::primary_button(ui, "Save", Some(width(ui, "Save"))).clicked() {
                            choice = 3;
                        }
                        if ui.add(egui::Button::new("Cancel").min_size(width(ui, "Cancel"))).clicked() || ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                            choice = 2;
                        }
                        let danger = if closing.is_some() { "Don't Save" } else { "Quit Without Saving" };
                        let red = if ui.visuals().dark_mode { egui::Color32::from_rgb(255, 115, 107) } else { egui::Color32::from_rgb(191, 38, 38) };
                        if ui.add(egui::Button::new(egui::RichText::new(danger).color(red)).min_size(width(ui, danger))).clicked() {
                            choice = 1;
                        }
                    });
                });
                match choice {
                    1 => {
                        match closing {
                            Some(i) => self.close_tab(i, true),
                            None => ctx.send_viewport_cmd(egui::ViewportCommand::Close),
                        }
                        self.pending = Pending::None;
                    }
                    2 => self.pending = Pending::None,
                    3 => {
                        // Save (Save As when the tab has no file yet), then close it.
                        if let Some(i) = closing {
                            self.active = i;
                            self.close_after_save = Some(self.tabs[i].id);
                            self.run_action(ctx, Action::Save);
                        }
                        self.pending = Pending::None;
                    }
                    _ => {}
                }
            }
        }
    }
}

impl eframe::App for AtomForgeApp {
    fn ui(&mut self, root: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = root.ctx().clone();
        self.dark = ctx.theme() == egui::Theme::Dark;
        crate::theme::set_screen(ctx.content_rect(), ctx.theme() == egui::Theme::Dark);
        if let Some(id) = self.close_after_save {
            match self.tabs.iter().position(|d| d.id == id) {
                Some(i) if !self.tabs[i].is_modified() => {
                    self.close_after_save = None;
                    self.close_tab(i, true);
                }
                Some(_) if !self.browser.open => self.close_after_save = None,
                None => self.close_after_save = None,
                _ => {}
            }
        }
        let now = std::time::Instant::now();
        self.frame_ms = (now - self.frame_clock).as_secs_f32() * 1000.0;
        self.frame_clock = now;
        self.run_script(&ctx);
        // Files dropped on the window open in new tabs.
        let dropped: Vec<String> = ctx.input(|i| i.raw.dropped_files.iter().map(|f| f.path().to_string_lossy().into_owned()).filter(|p| !p.is_empty()).collect());
        self.drop_files(dropped);
        self.shortcuts(&ctx);
        let title = match self.doc() {
            Some(d) => format!("AtomForge - {}{}", d.title, if d.is_modified() { " *" } else { "" }),
            None => "AtomForge".into(),
        };
        ctx.send_viewport_cmd(egui::ViewportCommand::Title(title));
        let bars = crate::theme::palette(self.dark);
        let bar = |margin: egui::Margin| egui::Frame::new().fill(bars.menu_bar).inner_margin(margin);
        egui::Panel::top("menu").frame(bar(egui::Margin::symmetric(10, 2)).stroke(egui::Stroke::new(1.0, bars.separator))).show(root, |ui| self.menus(ui));
        egui::Panel::top("toolbar").frame(bar(egui::Margin::symmetric(10, 6))).show(root, |ui| self.toolbar(ui));
        if !self.tabs.is_empty() {
            egui::Panel::top("tabs").frame(egui::Frame::new().fill(bars.menu_bar).inner_margin(egui::Margin { left: 4, right: 4, top: 4, bottom: 0 }).stroke(egui::Stroke::new(1.0, bars.separator))).show(root, |ui| self.tabs_bar(ui));
        }
        egui::Panel::bottom("status").frame(bar(egui::Margin::symmetric(12, 4)).stroke(egui::Stroke::new(1.0, bars.separator))).show(root, |ui| self.status_bar(ui));
        if self.pipeline.open && self.pipeline.docked && !self.tabs.is_empty() {
            let width = ctx.content_rect().width() * 0.32;
            egui::Panel::right("pipeline dock").resizable(true).default_size(width).frame(egui::Frame::new().fill(bars.window_fill()).inner_margin(egui::Margin::symmetric(12, 10)).stroke(egui::Stroke::new(1.0, bars.border))).show(root, |ui| {
                let d = &mut self.tabs[self.active];
                let shown = self.pipeline.show_docked(ui, &mut d.pipeline, &d.structure);
                self.docked_requests = Some(shown);
            });
        }
        egui::CentralPanel::default().frame(egui::Frame::NONE).show(root, |ui| self.viewport(ui));
        self.dialogs(&ctx);
        self.notices(&ctx);
        crate::theme::modal_backdrop(&ctx);
        // The previews the dialogs laid out this frame.
        let (colours, radii) = match self.doc() {
            Some(d) => (d.element_colours.clone(), d.element_radii.clone()),
            None => (self.look.colours.clone(), self.look.radii.clone()),
        };
        self.previews.render(&self.render_state, &self.display, &colours, &radii);
        self.autosave();
        self.capture(&ctx);
    }
}

impl AtomForgeApp {
    /// The screenshot hook: after 8 frames save the window and quit.
    fn capture(&mut self, ctx: &egui::Context) {
        let Some(path) = self.screenshot.clone() else { return };
        self.frames += 1;
        if self.frames == 8 {
            save_screenshot(ctx, path, true);
        }
        ctx.request_repaint();
    }

    /// UI script: pipeline open|start|text TEXT|add TYPE|set N NAME TEXT|enable N|disable N|move A B|remove N|replace N TYPE|update|bake|discard|save PATH|load PATH|colour on|off.
    fn script_pipeline(&mut self, rest: &str) -> Result<(), String> {
        let (what, value) = split(rest);
        let Some(d) = self.tabs.get_mut(self.active) else { return Err("no structure".into()) };
        let index = |v: &str| v.split_whitespace().next().and_then(|n| n.parse::<usize>().ok()).filter(|&n| n >= 1).map(|n| n - 1).ok_or(format!("pipeline {what}: a step number (from 1)"));
        let editor = &mut d.pipeline;
        let panel = &mut self.pipeline;
        match what {
            "open" => panel.open = true,
            "start" => editor.set_input(&d.structure),
            "text" => {
                if !editor.active() {
                    editor.set_input(&d.structure);
                }
                panel.apply_text(editor, value)?;
            }
            "add" => panel.add_step(editor, &d.structure, value.trim(), None),
            "set" => {
                let i = index(value)?;
                let (_, rest) = split(value);
                let (name, text) = split(rest);
                let t = af_pipeline::pipeline::find_modifier_type(&editor.modifier(i).type_id).ok_or("unknown step")?;
                let p = t.find_parameter(name).ok_or(format!("no parameter {name}"))?;
                editor.set_parameter(i, name, af_pipeline::pipeline::value_for(p, text)?)?;
            }
            "enable" | "disable" => editor.set_enabled(index(value)?, what == "enable")?,
            "move" => {
                let numbers: Vec<usize> = value.split_whitespace().filter_map(|w| w.parse().ok()).collect();
                if numbers.len() != 2 || numbers[0] == 0 || numbers[1] == 0 {
                    return Err("pipeline move FROM TO".into());
                }
                panel.move_step(editor, numbers[0] - 1, numbers[1] - 1);
            }
            "remove" => panel.remove_step(editor, index(value)?),
            "replace" => {
                let i = index(value)?;
                editor.replace(i, split(value).1.trim())?;
            }
            "duplicate" => {
                editor.duplicate(index(value)?)?;
            }
            "colour" => panel.colour_selection = value.trim() == "on",
            "dock" => {
                panel.open = true;
                panel.docked = true;
            }
            "undock" => panel.docked = false,
            "edit-dialog" => panel.edit_request = Some(index(value)?),
            "update" => {}
            "bake" | "discard" | "save" | "load" => {}
            other => return Err(format!("unknown pipeline command {other}")),
        }
        if what == "edit-dialog" {
            let i = self.pipeline.edit_request.take().unwrap_or(0);
            self.edit_step_in_dialog(i);
            return Ok(());
        }
        match what {
            "bake" => {
                let s = d.pipeline.bake();
                d.replace_structure(s);
                d.invalidate_colours();
            }
            "discard" => {
                let input = d.pipeline.input().cloned().ok_or("the pipeline is not running")?;
                d.pipeline.reset();
                d.structure = input;
                d.invalidate_colours();
            }
            "save" => self.pipeline.save_file(&d.pipeline, value.trim()),
            "load" => {
                if !d.pipeline.active() {
                    d.pipeline.set_input(&d.structure);
                }
                self.pipeline.load_file(&mut d.pipeline, value.trim());
            }
            // Showing the panel changes nothing: a structure edited outside the
            // pipeline stays, with the panel's warning (as when opened by hand).
            "open" | "dock" | "undock" => {}
            _ => {
                // The output as after a frame (update follows every edit).
                let requests = self.pipeline.update(&mut d.pipeline, true);
                if let Some(s) = requests.output {
                    d.structure = s;
                    d.selection.clear();
                    d.invalidate_colours();
                }
                if requests.show_property {
                    self.set_colour_mode(ColourMode::AtomProperty);
                }
            }
        }
        if !self.pipeline.message.is_empty() {
            return Err(self.pipeline.message.clone());
        }
        Ok(())
    }

    /// UI script: view style|colour|show|hide|labels|projection VALUE.
    fn script_view(&mut self, what: &str, value: &str) -> Result<(), String> {
        let key = value.to_lowercase().replace(['-', '_'], " ");
        match what {
            "style" => {
                let style = AtomStyle::ALL.into_iter().find(|s| s.label().to_lowercase() == key).ok_or(format!("no atom style {value}"))?;
                self.set_style(style);
            }
            "colour" | "color" => {
                let mode = match key.as_str() {
                    "element" => ColourMode::Element,
                    "crystal orientation" | "orientation" => ColourMode::CrystalOrientation,
                    "grain boundaries" | "grain boundary" => ColourMode::GrainBoundary,
                    "atom property" | "property" => ColourMode::AtomProperty,
                    _ => return Err(format!("no colour mode {value}")),
                };
                self.set_colour_mode(mode);
            }
            "show" | "hide" => {
                let on = what == "show";
                match key.as_str() {
                    "atoms" => self.view.show_atoms = on,
                    "bonds" => self.view.show_bonds = on,
                    "labels" => self.view.show_labels = on,
                    "legend" => self.view.show_legend = on,
                    "box" => self.view.show_box = on,
                    _ => return Err(format!("nothing called {value} to {what}")),
                }
            }
            "labels" => {
                self.view.label_content = LabelContent::ALL.into_iter().find(|c| c.label().to_lowercase().starts_with(&key)).ok_or(format!("no label content {value}"))?;
                self.view.show_labels = true;
            }
            "projection" => self.set_projection(key == "orthographic"),
            _ => return Err(format!("unknown view setting {what}")),
        }
        Ok(())
    }

    /// UI script: settings open|sizes|tab NAME|close|atom-size X|preset NAME|scheme NAME.
    fn script_settings(&mut self, what: &str, value: &str) -> Result<(), String> {
        match what {
            "open" => self.settings.open = true,
            "sizes" => self.settings.sizes_open = true,
            "close" => {
                self.settings.open = false;
                self.settings.sizes_open = false;
            }
            "tab" => {
                self.settings.open = true;
                if !self.settings.select_tab(value) {
                    return Err(format!("no settings tab {value}"));
                }
            }
            "atom-size" => {
                self.display.atom_scale = value.trim().parse::<f32>().map_err(|e| e.to_string())?.clamp(af_core::display::ATOM_SCALE.low, af_core::display::ATOM_SCALE.high);
                for d in self.tabs.iter_mut() {
                    d.scene_dirty = true;
                }
            }
            "preset" => {
                let preset = af_core::display::PRESETS.iter().find(|p| p.name.eq_ignore_ascii_case(value.trim())).ok_or(format!("no look {value}"))?;
                af_core::display::apply_preset(&mut self.display, preset);
                for d in self.tabs.iter_mut() {
                    d.scene_dirty = true;
                }
            }
            "scheme" => {
                let schemes = af_core::elements::colour_schemes();
                let i = schemes.iter().position(|s| s.name.eq_ignore_ascii_case(value.trim())).ok_or(format!("no colour scheme {value}"))?;
                let table = af_core::elements::scheme_colours(i as i32);
                if let Some(d) = self.doc_mut() {
                    let n = d.element_colours.len();
                    d.element_colours.copy_from_slice(&table[..n]);
                    d.invalidate_colours();
                }
            }
            _ => return Err(format!("unknown settings command {what}")),
        }
        Ok(())
    }

    fn fail(&mut self, message: String) {
        eprintln!("UI script: {message}");
        SCRIPT_FAILED.store(true, std::sync::atomic::Ordering::Relaxed);
        if let Some(script) = self.script.as_mut() {
            script.failed = Some(message);
        }
    }

    /// Runs the next line of the UI script (ATOMFORGE_UI_SCRIPT).
    fn run_script(&mut self, ctx: &egui::Context) {
        let Some(script) = self.script.as_mut() else { return };
        ctx.request_repaint();
        if script.failed.is_some() {
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            return;
        }
        let Some(line) = script.next_command() else {
            if script.finished() && script.wait_frames == 0 {
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
            return;
        };
        let result = self.script_command(ctx, &line);
        if let Err(e) = result {
            self.fail(format!("{line}: {e}"));
        }
    }

    fn script_command(&mut self, ctx: &egui::Context, line: &str) -> Result<(), String> {
        let (command, rest) = split(line);
        let words: Vec<&str> = rest.split_whitespace().collect();
        let number = |i: usize| words.get(i).and_then(|w| w.parse::<f64>().ok());
        if command == "try" {
            // try COMMAND: runs a command whose failure is part of the scene (tour screenshots).
            let _ = self.script_command(ctx, rest);
            return Ok(());
        }
        match command {
            "open" => {
                let before = self.tabs.len();
                // Scripts convert large files on this thread.
                let p = std::path::Path::new(rest);
                if crate::large_data::opens_as_cloud(p) && !af_cloud::is_cloud_file(p) {
                    let cloud = crate::large_data::CloudConversion::convert_now(p)?;
                    if !self.open_cloud(&cloud.to_string_lossy()) {
                        return Err(self.status.clone());
                    }
                    return Ok(());
                }
                self.open_path(rest);
                if self.tabs.len() == before {
                    return Err(self.status.clone());
                }
            }
            "step" => {
                if !self.operation.open_step(rest) {
                    return Err(format!("no step {rest}"));
                }
            }
            "set" => {
                let (name, text) = split(rest);
                self.operation.set_parameter_text(name, text)?;
            }
            "target" => {
                let (which, condition) = split(rest);
                let target = match which {
                    "all" => crate::operation::Target::AllAtoms,
                    "view" => crate::operation::Target::ViewSelection,
                    _ => crate::operation::Target::Condition,
                };
                self.operation.set_target(target, condition);
            }
            "apply" => {
                let structure = self.structure();
                let selection = self.doc().map(|d| d.selected_atoms()).unwrap_or_default();
                let outcome = self.operation.apply_now(&structure, &selection)?;
                self.apply_outcome(outcome);
            }
            "update-step" => {
                // Update Step in the operation dialog editing a pipeline step.
                let edit = self.operation.update_step().ok_or("the operation dialog is not editing a pipeline step")?;
                self.update_step(edit);
            }
            "tool" => {
                let id = af_science::catalog::find_tool(rest).ok_or_else(|| format!("no tool {rest}"))?.id;
                self.tools.open_tool(id);
            }
            "toolset" => {
                let (name, text) = split(rest);
                self.tools.set_text(name, text);
            }
            "toolwin" => { let (c, a) = split(rest); let s = self.structure(); self.tools.script(c, a, &s)?; }
            "run" => {
                let structure = self.structure();
                self.tools.run_blocking(&structure)?;
            }
            "pick" => {
                let (x, y) = (number(0).ok_or("pick X Y")?, number(1).ok_or("pick X Y")?);
                let shift = words.get(2) == Some(&"shift");
                self.pick_at(x, y, self.viewport_size, shift);
            }
            "select-all" => {
                if let Some(d) = self.doc_mut() {
                    d.selection = (0..d.scene.instances.positions.len()).collect();
                    d.scene_dirty = true;
                }
            }
            "clear-selection" => {
                if let Some(d) = self.doc_mut() {
                    d.selection.clear();
                    d.scene_dirty = true;
                }
            }
            "delete" => {
                if let Some(d) = self.doc_mut() {
                    d.delete_selected();
                }
            }
            "undo" => {
                if let Some(d) = self.doc_mut() {
                    d.undo();
                }
            }
            "redo" => {
                if let Some(d) = self.doc_mut() {
                    d.redo();
                }
            }
            "reset-view" => {
                if let Some(d) = self.doc_mut() {
                    d.pending_default_view = true;
                }
            }
            "orbit" => {
                let (dx, dy) = (number(0).unwrap_or(0.0) as f32, number(1).unwrap_or(0.0) as f32);
                if let Some(d) = self.doc_mut() {
                    d.camera.yaw -= dx * 0.3;
                    d.camera.pitch += dy * 0.3;
                }
            }
            "export-options" => {
                // export-options SCALE [background on|off] [gizmo on|off]
                self.browser.export_scale = number(0).unwrap_or(1.0).clamp(1.0, 8.0) as i32;
                let w: Vec<&str> = rest.split_whitespace().collect();
                for pair in w[1..].chunks(2) {
                    let on = pair.get(1).is_some_and(|v| *v == "on");
                    match pair[0] {
                        "background" => self.browser.export_background = on,
                        "gizmo" => self.browser.export_gizmo = on,
                        other => return Err(format!("export-options: unknown option {other}")),
                    }
                }
            }
            "export-image" => {
                // export-image PATH: File > Export Image with the current options.
                self.export_image(rest.trim());
                if !self.status.starts_with("Image exported") {
                    return Err(self.status.clone());
                }
            }
            "expect-file" => {
                // expect-file PATH [W H | view N]: the file exists (and is an image of
                // this size, or of N times the view's size).
                let path = words.first().ok_or("expect-file PATH [W H | view N]")?;
                let meta = std::fs::metadata(path).map_err(|e| format!("{path}: {e}"))?;
                if meta.len() == 0 {
                    return Err(format!("{path} is empty"));
                }
                let size = if words.get(1) == Some(&"view") {
                    let n = number(2).unwrap_or(1.0);
                    Some((self.viewport_size.0 as f64 * n, self.viewport_size.1 as f64 * n))
                } else {
                    number(1).zip(number(2))
                };
                if let Some((w, h)) = size {
                    let (fw, fh) = crate::image_export::image_size(path)?;
                    if fw as f64 != w || fh as f64 != h {
                        return Err(format!("{path} is {fw} x {fh}, expected {w} x {h}"));
                    }
                }
            }
            "save" => {
                self.save_to(rest);
                if !self.status.starts_with("Saved") {
                    return Err(self.status.clone());
                }
            }
            "tab" => {
                let n = number(0).ok_or("tab N")? as usize;
                if n == 0 || n > self.tabs.len() {
                    return Err(format!("no tab {n}"));
                }
                self.active = n - 1;
                self.tabs[self.active].upload_needed = true;
            }
            "expect-atoms" => {
                let want = number(0).ok_or("expect-atoms N")? as usize;
                let got = self.doc().map(|d| d.structure.atoms.len()).unwrap_or(0);
                if got != want {
                    return Err(format!("expected {want} atoms, the structure has {got}"));
                }
            }
            "expect-selected" => {
                let want = number(0).ok_or("expect-selected N")? as usize;
                let got = self.doc().map(|d| d.selected_atoms().len()).unwrap_or(0);
                if got != want {
                    return Err(format!("expected {want} selected atoms, {got} are selected"));
                }
            }
            "expect-tabs" => {
                let want = number(0).ok_or("expect-tabs N")? as usize;
                if self.tabs.len() != want {
                    return Err(format!("expected {want} tabs, {} are open", self.tabs.len()));
                }
            }
            "capture-folder" => {
                // capture-folder DIR: where `capture` writes images and the bounds file.
                self.capture.folder = rest.trim().to_string();
                self.capture.bounds.clear();
            }
            "capture" => {
                // capture FILE popups|menus|top|viewport|tools|full: a manual screenshot
                // and the rectangle of what it shows (the C++ TutorialCapture).
                let (file, target) = split(rest);
                let r = self.capture_rect(ctx, target.trim());
                let k = ctx.pixels_per_point();
                let screen = ctx.content_rect();
                let pad = 6.0;
                let r = egui::Rect::from_min_max(
                    egui::pos2((r.min.x - pad).max(0.0), (r.min.y - pad).max(0.0)),
                    egui::pos2((r.max.x + pad).min(screen.max.x), (r.max.y + pad).min(screen.max.y)),
                );
                self.capture.bounds.push((file.to_string(), [(r.min.x * k) as i32, (r.min.y * k) as i32, (r.max.x * k) as i32, (r.max.y * k) as i32]));
                let path = std::path::Path::new(&self.capture.folder).join(file).to_string_lossy().into_owned();
                save_screenshot(ctx, path, false);
                if let Some(s) = self.script.as_mut() {
                    s.wait_frames = 4;
                }
            }
            "capture-done" => {
                // capture-done NAME: writes FOLDER/bounds-NAME.json.
                let mut text = String::from("{");
                for (i, (file, b)) in self.capture.bounds.iter().enumerate() {
                    text += &format!("{}\n  \"{}\": [{}, {}, {}, {}]", if i > 0 { "," } else { "" }, file, b[0], b[1], b[2], b[3]);
                }
                text += "\n}\n";
                let path = std::path::Path::new(&self.capture.folder).join(format!("bounds-{}.json", rest.trim()));
                std::fs::write(path, text).map_err(|e| e.to_string())?;
            }
            "screenshot" => {
                save_screenshot(ctx, rest.to_string(), false);
                if let Some(s) = self.script.as_mut() {
                    s.wait_frames = 4;
                }
            }
            "expect-cloud" => {
                // expect-cloud ATOMS: the active tab is an atom cloud of ATOMS atoms.
                let want = number(0).ok_or("expect-cloud ATOMS")? as u64;
                let cloud = self.doc().and_then(|d| d.cloud.as_ref()).ok_or("the active tab is not an atom cloud")?;
                if cloud.stats().total_atoms != want {
                    return Err(format!("the cloud has {} atoms", cloud.stats().total_atoms));
                }
            }
            "expect-cloud-drawn" => {
                // expect-cloud-drawn MIN [MAX]: atoms drawn in the last frame.
                let cloud = self.doc().and_then(|d| d.cloud.as_ref()).ok_or("the active tab is not an atom cloud")?;
                let drawn = cloud.stats().drawn_atoms as f64;
                if drawn < number(0).unwrap_or(1.0) || number(1).is_some_and(|m| drawn > m) {
                    return Err(format!("{drawn} atoms drawn ({:?})", cloud.stats()));
                }
            }
            "cloud-setting" => {
                // cloud-setting NAME VALUE (the Large Dataset panel).
                let (name, value) = split(rest);
                let cloud = self.tabs.get_mut(self.active).and_then(|d| d.cloud.as_mut()).ok_or("the active tab is not an atom cloud")?;
                crate::large_data::set_cloud_setting(cloud, name, value)?;
            }
            "wait-cloud" => {
                // wait-cloud [SECONDS]: until every visible chunk has streamed in.
                let complete = self.doc().and_then(|d| d.cloud.as_ref()).ok_or("the active tab is not an atom cloud")?.stats().complete;
                if complete {
                    self.cloud_wait = None;
                } else {
                    let deadline = *self.cloud_wait.get_or_insert_with(|| std::time::Instant::now() + std::time::Duration::from_secs_f64(number(0).unwrap_or(60.0)));
                    if std::time::Instant::now() > deadline {
                        self.cloud_wait = None;
                        return Err("the cloud did not finish streaming in time".into());
                    }
                    if let Some(s) = self.script.as_mut() {
                        s.repeat_last();
                        s.wait_frames = 1;
                    }
                }
            }
            "cloud-generate" => {
                // cloud-generate PATH LATTICE A NX NY NZ ELEMENT: writes a crystal cloud (--cloud generate).
                let w: Vec<&str> = rest.split_whitespace().collect();
                if w.len() != 7 {
                    return Err("cloud-generate PATH LATTICE A NX NY NZ ELEMENT".into());
                }
                let n = |k: usize| w[k].parse::<u64>().map_err(|_| format!("{} is not a cell count", w[k]));
                let a = w[2].parse::<f64>().map_err(|_| "A is a number".to_string())?;
                if let Some(parent) = std::path::Path::new(w[0]).parent() {
                    let _ = std::fs::create_dir_all(parent);
                }
                af_cloud::generate_crystal(std::path::Path::new(w[0]), w[1], a, 0.0, [n(3)?, n(4)?, n(5)?], w[6], "", None, af_cloud::BuildOptions::default())?;
            }
            "copy" => {
                // copy SRC DST: copies a file (test inputs into rust/target).
                let (from, to) = split(rest);
                std::fs::copy(from, to).map_err(|e| format!("cannot copy {from} to {to}: {e}"))?;
            }
            "remove" => {
                // remove PATH: deletes a file if it exists (fresh test outputs).
                let _ = std::fs::remove_file(rest);
            }
            "wait" => {
                let n = number(0).unwrap_or(1.0) as u32;
                if let Some(s) = self.script.as_mut() {
                    s.wait_frames = n;
                }
            }
            "view" => {
                let (what, value) = split(rest);
                self.script_view(what, value)?;
            }
            "settings" => {
                let (what, value) = split(rest);
                self.script_settings(what, value)?;
            }
            "expect-bonds" => {
                let want = number(0).ok_or("expect-bonds N")? as usize;
                let view = self.view;
                let display = self.display;
                let Some(d) = self.doc_mut() else { return Err("no structure".into()) };
                if d.scene_dirty {
                    d.rebuild_scene(&display, &view);
                }
                let got = if view.bonds_drawn() { d.scene.bonds.starts.len() } else { 0 };
                if got != want {
                    return Err(format!("expected {want} bonds drawn, {got} are"));
                }
            }
            "expect-colour" => {
                // expect-colour ATOM R G B: the atom's colour (0-1, within 0.01).
                let atom = number(0).ok_or("expect-colour ATOM R G B")? as usize;
                let want = [number(1).unwrap_or(-1.0), number(2).unwrap_or(-1.0), number(3).unwrap_or(-1.0)];
                let view = self.view;
                let display = self.display;
                let Some(d) = self.doc_mut() else { return Err("no structure".into()) };
                if d.scene_dirty {
                    d.rebuild_scene(&display, &view);
                }
                let a = d.structure.atoms.get(atom).ok_or(format!("no atom {atom}"))?;
                let got = [a.r as f64, a.g as f64, a.b as f64];
                if (0..3).any(|k| (got[k] - want[k]).abs() > 0.01) {
                    return Err(format!("atom {atom} has colour {got:?}, expected {want:?}"));
                }
            }
            "expect-coloured" => {
                // expect-coloured R G B N: N atoms have this colour (within 0.01).
                let want = [number(0).unwrap_or(-1.0), number(1).unwrap_or(-1.0), number(2).unwrap_or(-1.0)];
                let count = number(3).ok_or("expect-coloured R G B N")? as usize;
                let view = self.view;
                let display = self.display;
                let Some(d) = self.doc_mut() else { return Err("no structure".into()) };
                if d.scene_dirty {
                    d.rebuild_scene(&display, &view);
                }
                let got = d.structure.atoms.iter().filter(|a| ((a.r as f64 - want[0]).abs() <= 0.01) && ((a.g as f64 - want[1]).abs() <= 0.01) && ((a.b as f64 - want[2]).abs() <= 0.01)).count();
                if got != count {
                    return Err(format!("expected {count} atoms coloured {want:?}, {got} are"));
                }
            }
            "measure" => {
                let view = self.view;
                let display = self.display;
                if let Some(d) = self.doc_mut() {
                    if d.scene_dirty {
                        d.rebuild_scene(&display, &view);
                    }
                }
                match rest.trim() {
                    "distance" => self.measure_selection(Kind::Distance),
                    "angle" => self.measure_selection(Kind::Angle),
                    "dihedral" => self.measure_selection(Kind::Dihedral),
                    "atom-info" => self.atom_info(),
                    other => return Err(format!("no measurement {other}")),
                }
            }
            "expect-measure" => {
                // expect-measure distance|angle|dihedral|atom-info TEXT: the value (or the Atom Info text) contains TEXT.
                let (what, text) = split(rest);
                let got = match what {
                    "distance" => self.measure.get(Kind::Distance).map(|m| format!("{} {}", m.value, m.message)),
                    "angle" => self.measure.get(Kind::Angle).map(|m| format!("{} {}", m.value, m.message)),
                    "dihedral" => self.measure.get(Kind::Dihedral).map(|m| format!("{} {}", m.value, m.message)),
                    "atom-info" => self.measure.atom_info.clone(),
                    other => return Err(format!("no measurement {other}")),
                }
                .ok_or(format!("no {what} result is shown"))?;
                if !got.contains(text) {
                    return Err(format!("{what} shows {got:?}, without {text:?}"));
                }
            }
            "structure-info" => {
                // Opens the window and waits (up to 10 s) for the background symmetry.
                self.info.open = true;
                let Some(d) = self.tabs.get(self.active) else { return Err("no structure".into()) };
                let (structure, revision) = (d.structure.clone(), crate::structure_info::fingerprint(&d.structure));
                let start = std::time::Instant::now();
                while self.info.data(revision).is_none() && start.elapsed().as_secs() < 10 {
                    self.info.update(&structure, revision);
                    std::thread::sleep(std::time::Duration::from_millis(10));
                }
                if self.info.data(revision).is_none() {
                    return Err("structure info did not finish".into());
                }
            }
            "expect-info" => {
                // expect-info LABEL = TEXT: the Structure Info row LABEL contains TEXT.
                let (label, text) = rest.split_once('=').map(|(l, t)| (l.trim(), t.trim())).ok_or("expect-info LABEL = TEXT")?;
                let Some(d) = self.tabs.get(self.active) else { return Err("no structure".into()) };
                let data = self.info.data(crate::structure_info::fingerprint(&d.structure));
                let rows: Vec<(String, String)> = crate::structure_info::summary_rows(&d.structure, data).into_iter().chain(crate::structure_info::symmetry_rows(data)).collect();
                let value = rows.iter().find(|(l, _)| l == label).map(|(_, v)| v.clone()).ok_or(format!("no Structure Info row {label}"))?;
                if !value.contains(text) {
                    return Err(format!("{label} is {value:?}, without {text:?}"));
                }
            }
            "pick-atom" => {
                // pick-atom N [shift]: picks atom N's first instance (as a click on it).
                let atom = number(0).ok_or("pick-atom N")? as i32;
                let shift = words.get(1) == Some(&"shift");
                let view = self.view;
                let display = self.display;
                let Some(d) = self.doc_mut() else { return Err("no structure".into()) };
                if d.scene_dirty {
                    d.rebuild_scene(&display, &view);
                }
                let i = d.scene.instances.atom_indices.iter().position(|&a| a == atom).ok_or(format!("no atom {atom}"))?;
                if shift {
                    d.selection.push(i);
                } else {
                    d.selection = vec![i];
                }
                d.scene_dirty = true;
            }
            "command" => {
                // command QUERY: runs the best enabled match of Find a Command.
                let commands = self.commands();
                let best = crate::palette::ranked(&commands, rest).into_iter().find(|c| c.enabled).ok_or(format!("no command matches {rest}"))?;
                let action = best.action.clone();
                self.status = format!("Ran {} > {}", best.path, best.label);
                self.run_action(ctx, action);
            }
            "expect-camera" => {
                // expect-camera YAW PITCH: the camera angles (degrees, within 0.01).
                let (yaw, pitch) = (number(0).ok_or("expect-camera YAW PITCH")?, number(1).ok_or("expect-camera YAW PITCH")?);
                let d = self.doc().ok_or("no structure")?;
                // An optional third number checks the roll.
                let roll_off = number(2).is_some_and(|r| (d.camera.roll as f64 - r).abs() > 0.01);
                if (d.camera.yaw as f64 - yaw).abs() > 0.01 || (d.camera.pitch as f64 - pitch).abs() > 0.01 || roll_off {
                    return Err(format!("camera yaw {} pitch {} roll {}, expected {}", d.camera.yaw, d.camera.pitch, d.camera.roll, rest));
                }
            }
            "expect-status" => {
                if !self.status.contains(rest) {
                    return Err(format!("status is {:?}, without {rest:?}", self.status));
                }
            }
            "save-project" => {
                self.save_project(rest);
                if !self.status.starts_with("Project saved") {
                    return Err(self.status.clone());
                }
            }
            "open-project" => {
                let before = self.tabs.len();
                self.open_project(rest);
                if self.tabs.len() == before {
                    return Err(self.status.clone());
                }
            }
            "close-all" => {
                self.tabs.clear();
                self.active = 0;
            }
            "box" => {
                // box X0 Y0 X1 Y1 [add]: a right drag in Box Select Mode (viewport pixels).
                let corners: Vec<f32> = words.iter().take(4).filter_map(|w| w.parse().ok()).collect();
                if corners.len() != 4 {
                    return Err("box X0 Y0 X1 Y1 [add]".into());
                }
                self.select_mode = SelectMode::Box;
                self.refresh_scene();
                self.select_region(&[[corners[0], corners[1]], [corners[2], corners[3]]], self.viewport_size, words.get(4) == Some(&"add"));
            }
            "lasso" => {
                // lasso X Y X Y ... [add]: a right drag in Lasso Select Mode through these points.
                let numbers: Vec<f32> = words.iter().filter_map(|w| w.parse().ok()).collect();
                let mut region: Vec<[f32; 2]> = numbers.chunks_exact(2).map(|c| [c[0], c[1]]).collect();
                if region.len() < 3 {
                    return Err("lasso needs 3 points or more".into());
                }
                region.insert(0, region[0]);
                self.select_mode = SelectMode::Lasso;
                self.refresh_scene();
                self.select_region(&region, self.viewport_size, words.last() == Some(&"add"));
            }
            "rotate" => {
                // rotate x|y|z ANGLE: the toolbar's Rotate.
                let axis = ["x", "y", "z"].iter().position(|a| words.first() == Some(a)).ok_or("rotate x|y|z ANGLE")?;
                self.rotate_crystal(axis, number(1).ok_or("rotate x|y|z ANGLE")?);
            }
            "substitute" | "insert" => {
                // substitute ELEMENT / insert ELEMENT: the menu's Substitute Atom and Insert Atom at Midpoint.
                let z = af_core::elements::atomic_number(rest.trim());
                if z <= 0 {
                    return Err(format!("no element {rest}"));
                }
                self.refresh_scene();
                self.apply_pick(if command == "substitute" { Pick::Substitute } else { Pick::InsertMidpoint }, z);
            }
            "grab" => {
                // grab start | move DX DY | axis x|y|z | confirm | cancel (pixels from where it started).
                match words.first().copied() {
                    Some("start") => {
                        self.grab_start([0.0, 0.0]);
                        if self.grab.is_none() {
                            return Err("nothing selected to grab".into());
                        }
                    }
                    Some("move") => self.grab_move([number(1).unwrap_or(0.0) as f32, number(2).unwrap_or(0.0) as f32], self.viewport_size),
                    Some("axis") => {
                        let axis = match words.get(1).copied() {
                            Some("x") => af_render::grab::Axis::X,
                            Some("y") => af_render::grab::Axis::Y,
                            Some("z") => af_render::grab::Axis::Z,
                            _ => return Err("grab axis x|y|z".into()),
                        };
                        if let Some(g) = self.grab.as_mut() {
                            g.axis = g.axis.toggled(axis);
                        }
                    }
                    Some("confirm") => self.grab_end(true),
                    Some("cancel") => self.grab_end(false),
                    _ => return Err("grab start|move DX DY|axis A|confirm|cancel".into()),
                }
            }
            "expect-moved" => {
                // expect-moved ATOM X Y Z: a checked coordinate differs ("*" unchecked).
                let atom = number(0).ok_or("expect-moved ATOM X Y Z")? as usize;
                let d = self.doc().ok_or("no structure")?;
                let a = d.structure.atoms.get(atom).ok_or(format!("no atom {atom}"))?;
                let want = [number(1), number(2), number(3)];
                let moved = [a.x, a.y, a.z].iter().zip(want).any(|(&got, w)| w.is_some_and(|w| (got - w).abs() > 1e-3));
                if !moved {
                    return Err(format!("atom {atom} did not move from ({}, {}, {})", a.x, a.y, a.z));
                }
            }
            "expect-position" => {
                // expect-position ATOM X Y Z (within 1e-3 A).
                let atom = number(0).ok_or("expect-position ATOM X Y Z")? as usize;
                let d = self.doc().ok_or("no structure")?;
                let a = d.structure.atoms.get(atom).ok_or(format!("no atom {atom}"))?;
                // "*" leaves a coordinate unchecked.
                let want = [number(1), number(2), number(3)];
                let off = |got: f64, want: Option<f64>| want.is_some_and(|w| (got - w).abs() > 1e-3);
                if off(a.x, want[0]) || off(a.y, want[1]) || off(a.z, want[2]) {
                    return Err(format!("atom {atom} is at ({}, {}, {}), expected {want:?}", a.x, a.y, a.z));
                }
            }
            "pipeline" => self.script_pipeline(rest)?,
            "expect-stage" => {
                // expect-stage N ATOMS SELECTED (1-based step; after the last update).
                let n = number(0).ok_or("expect-stage N ATOMS SELECTED")? as usize;
                let d = self.doc().ok_or("no structure")?;
                let stage = d.pipeline.stages().get(n.wrapping_sub(1)).ok_or(format!("no evaluated step {n}"))?;
                if !stage.error.is_empty() {
                    return Err(format!("step {n} failed: {}", stage.error));
                }
                let (atoms, selected) = (number(1).unwrap_or(-1.0) as usize, number(2).unwrap_or(-1.0) as usize);
                if stage.atoms != atoms || stage.selected != selected {
                    return Err(format!("step {n}: {} atoms, {} selected; expected {atoms}, {selected}", stage.atoms, stage.selected));
                }
            }
            "expect-plots" => {
                // expect-plots N [TITLE]: the tool's last run has N plots (one titled TITLE).
                let n = number(0).ok_or("expect-plots N [TITLE]")? as usize;
                let titles = self.tools.plot_titles();
                let title = rest.split_once(' ').map(|(_, t)| t.trim()).unwrap_or("");
                if titles.len() != n || (!title.is_empty() && !titles.iter().any(|t| t == title)) {
                    return Err(format!("the tool has plots {titles:?}"));
                }
            }
            "expect-property" => {
                // expect-property VALUE N: N atoms have this per-atom value.
                let (value, n) = (number(0).ok_or("expect-property VALUE N")?, number(1).ok_or("expect-property VALUE N")? as usize);
                let d = self.doc().ok_or("no structure")?;
                let got = d.structure.atom_property.iter().filter(|&&v| v == value).count();
                if got != n {
                    return Err(format!("{got} atoms have {} = {value}, expected {n}", d.structure.atom_property_name));
                }
            }
            "colour-by-tool" => {
                // colour-by-tool [N]: colour the atoms by the tool's per-atom values N (from 1).
                let k = number(0).unwrap_or(1.0) as usize;
                let p = self.tools.properties().get(k.saturating_sub(1)).cloned().ok_or("the tool gave no such per-atom values")?;
                self.colour_by_property(p);
            }
            "bulk" => {
                // bulk open | set WHAT VALUE | build
                let (what, value) = split(rest);
                let colours = self.doc().map(|d| d.element_colours.clone()).unwrap_or_else(|| self.look.colours.clone());
                match what {
                    "open" => self.bulk.open = true,
                    // bulk systems: the crystal-system list held open.
                    "systems" => self.bulk.show_systems = true,
                    "set" => {
                        let (name, text) = split(value);
                        self.bulk.set(name, text, &colours)?;
                    }
                    "update-step" => {
                        let edit = self.bulk.update_step().ok_or("the bulk dialog is not editing a pipeline step")?;
                        self.update_step(edit);
                    }
                    "build" => match self.bulk.build(&colours) {
                        Some(s) => self.take_built(s, "Built a bulk crystal"),
                        None => return Err(self.bulk.result().message.clone()),
                    },
                    other => return Err(format!("unknown bulk command {other}")),
                }
            }
            "dislocation" => {
                // dislocation open | set WHAT VALUE | detect | generate | replace
                let (what, value) = split(rest);
                let view = self.structure();
                let source = self.dislocation.source(&view).clone();
                match what {
                    "open" => self.dislocation.open = true,
                    "update-step" => {
                        let edit = self.dislocation.update_step().ok_or("the dislocation dialog is not editing a pipeline step")?;
                        self.update_step(edit);
                    }
                    "expect-options" => {
                        let options = self.dislocation.step_options();
                        if !options.contains(value.trim()) {
                            return Err(format!("the step options are {options}"));
                        }
                    }
                    "set" => {
                        let (name, text) = split(value);
                        self.dislocation.set(name, text)?;
                    }
                    "detect" => {
                        let d = self.dislocation.detect_now(&source);
                        if !d.success {
                            return Err(d.message.clone());
                        }
                    }
                    "generate" => {
                        let r = self.dislocation.generate_now(&source);
                        if !r.success {
                            return Err(r.message.clone());
                        }
                    }
                    "replace" => {
                        let r = self.dislocation.result.as_ref().filter(|r| r.success).ok_or("generate a dislocation first")?;
                        let s = crate::dislocation_dialog::result_structure(r);
                        self.take_dislocation(s);
                    }
                    other => return Err(format!("unknown dislocation command {other}")),
                }
            }
            "dialog" => {
                // dialog STEP COMMAND ARGS: a dedicated dialog's script command.
                let (step, tail) = split(rest);
                let (command, args) = split(tail);
                let structure = self.structure();
                let colours = self.doc().map(|d| d.element_colours.clone()).unwrap_or_else(|| self.look.colours.clone());
                let input = crate::structure_dialog::Input { source: &structure, colours: &colours };
                let d = self.dialogs.iter_mut().find(|d| d.steps().contains(&step)).ok_or(format!("no dialog for {step}"))?;
                if command == "open" {
                    d.open(step);
                } else if let Some(o) = d.script(command, args, &input)? {
                    self.take_outcome(o);
                }
            }
            "minimise" | "restore" => {
                // minimise TITLE | restore TITLE: a builder or analysis dialog (tour screenshots).
                crate::theme::minimise(rest.trim(), command == "minimise");
            }
            "theme" => {
                // theme light|dark
                ctx.set_theme(if rest.trim() == "dark" { egui::Theme::Dark } else { egui::Theme::Light });
            }
            "menu" => {
                // menu NAME [> SUBMENU] | menu off: holds menus open for a screenshot.
                self.tour_menus = if rest.trim() == "off" { Vec::new() } else { rest.split('>').map(|m| m.trim().to_string()).filter(|m| !m.is_empty()).collect() };
            }
            "hover" => {
                // hover FX FY: the pointer rests at this fraction of the view (tooltip screenshots); hover off.
                self.tour_hover = match (number(0), number(1)) {
                    (Some(x), Some(y)) => {
                        let s = self.viewport_size;
                        let k = ctx.pixels_per_point();
                        Some(self.viewport_origin + egui::vec2(x as f32 * s.0 as f32 / k, y as f32 * s.1 as f32 / k))
                    }
                    _ => None,
                };
            }
            "element-picker" => {
                // element-picker open substitute|insert | type LETTERS
                let (c, a) = split(rest);
                self.picker.script(c, a)?;
            }
            "palette" => {
                // palette [QUERY]: Find a Command, with this query typed.
                self.palette.show_palette();
                self.palette.set_query(rest.trim());
            }
            "notify" => {
                // notify TEXT: a notice, as an operation's message; notify alone clears them.
                self.status = rest.trim().to_string();
                if self.status.is_empty() {
                    self.notifications.clear();
                }
            }
            "select-atoms" => {
                // select-atoms I J ...: selects these atoms (0-based).
                let wanted: std::collections::HashSet<usize> = rest.split_whitespace().filter_map(|w| w.parse().ok()).collect();
                let d = self.doc_mut().ok_or("no structure")?;
                d.selection = d.scene.instances.atom_indices.iter().enumerate().filter(|(_, &a)| a >= 0 && wanted.contains(&(a as usize))).map(|(i, _)| i).collect();
                d.scene_dirty = true;
            }
            "browser-name" => {
                // browser-name NAME: the file name typed in the open file browser.
                self.browser.name = rest.trim().to_string();
            }
            "context-menu" => {
                // context-menu: the selected atoms' right-click menu (at the view's centre).
                let k = ctx.pixels_per_point();
                let s = self.viewport_size;
                self.tour_context_menu = Some(self.viewport_origin + egui::vec2(s.0 as f32 / k * 0.5, s.1 as f32 / k * 0.45));
            }
            "close-windows" => {
                // Closes what a tour step opened.
                self.tour_menus.clear();
                egui::Popup::close_all(ctx);
                self.measure.atom_info = None;
                self.tour_context_menu = None;
                self.operation.open = false;
                self.tools.open = false;
                self.pipeline.open = false;
                self.about = false;
                self.trajectory.open = false;
                self.bulk.open = false;
                self.bulk.show_systems = false;
                self.dislocation.open = false;
                self.cna.open = false;
                self.rdf.open = false;
                self.adf.open = false;
                self.sro.open = false;
                self.info.open = false;
                self.settings.open = false;
                self.settings.sizes_open = false;
                self.overlays.planes_open = false;
                self.overlays.directions_open = false;
                self.polyhedra.settings_open = false;
                self.shortcuts_open = false;
                self.pending = Pending::None;
                self.picker = ElementPicker::default();
                self.tour_hover = None;
                self.palette.open = false;
                self.browser.open = false;
                for d in self.dialogs.iter_mut() {
                    d.close();
                }
            }
            "drop" => {
                // drop PATH: as if the file were dropped on the window.
                self.drop_files(vec![rest.trim().to_string()]);
            }
            "expect-shifted" => {
                // expect-shifted MIN: the last generated dislocation shifted at least MIN atoms.
                let r = self.dislocation.result.as_ref().ok_or("no dislocation generated")?;
                if (r.shifted_atom_count as f64) < number(0).unwrap_or(1.0) {
                    return Err(format!("{} atoms shifted", r.shifted_atom_count));
                }
            }
            "expect-loop" => {
                // expect-loop MIN: the structure's dislocation loop has at least MIN points (and is shown).
                let d = self.doc().ok_or("no structure")?;
                let n = d.structure.dislocation_loop_points.len();
                if (n as f64) < number(0).unwrap_or(3.0) || !self.overlays.show_dislocation_lines {
                    return Err(format!("loop has {n} points, shown {}", self.overlays.show_dislocation_lines));
                }
            }
            "rdf" => {
                // rdf open | pair REF TARGET | rmax R | bins N | pbc on|off | normalize on|off | distortion on|off | run
                let (what, value) = split(rest);
                let on = value.trim() == "on";
                let structure = self.structure();
                let p = &mut self.rdf.params;
                match what {
                    "open" => self.rdf.open = true,
                    "pair" => {
                        let (a, b) = split(value);
                        self.rdf.set_pair(a, b)?;
                    }
                    "rmax" => p.r_max = value.trim().parse().map_err(|e| format!("rmax: {e}"))?,
                    "bins" => p.bin_count = value.trim().parse().map_err(|e| format!("bins: {e}"))?,
                    "pbc" => p.use_pbc = on,
                    "normalize" => p.normalize = on,
                    "distortion" => p.enable_distortion_analysis = on,
                    "run" => {
                        let r = self.rdf.run_now(&structure);
                        if !r.valid {
                            return Err(r.message.clone());
                        }
                    }
                    other => return Err(format!("unknown rdf command {other}")),
                }
            }
            "adf" => {
                // adf open | cutoff R | bins N | smooth N | run
                let (what, value) = split(rest);
                let structure = self.structure();
                let p = &mut self.adf.params;
                let parse = |v: &str| v.trim().parse::<f64>().map_err(|e| format!("adf {what}: {e}"));
                match what {
                    "open" => self.adf.open = true,
                    "cutoff" => p.r_cutoff = parse(value)? as f32,
                    "bins" => p.bin_count = parse(value)? as i32,
                    "smooth" => p.smooth_passes = parse(value)? as i32,
                    "run" => {
                        let r = self.adf.run_now(&structure);
                        if !r.valid {
                            return Err(r.message.clone());
                        }
                    }
                    other => return Err(format!("unknown adf command {other}")),
                }
            }
            "cna" => {
                // cna open | cutoff R | scale S | run
                let (what, value) = split(rest);
                let structure = self.structure();
                let parse = |v: &str| v.trim().parse::<f32>().map_err(|e| format!("cna {what}: {e}"));
                match what {
                    "open" => self.cna.open = true,
                    "cutoff" => self.cna.params.cutoff_radius = parse(value)?,
                    "scale" => self.cna.params.cutoff_scale = parse(value)?,
                    "run" => {
                        let r = self.cna.run_now(&structure);
                        if !r.valid {
                            return Err(r.message.clone());
                        }
                    }
                    other => return Err(format!("unknown cna command {other}")),
                }
            }
            "sro" => {
                // sro open | shells N | tolerance T | run
                let (what, value) = split(rest);
                let structure = self.structure();
                let parse = |v: &str| v.trim().parse::<f64>().map_err(|e| format!("sro {what}: {e}"));
                match what {
                    "open" => self.sro.open = true,
                    "tab" => self.sro.tab = (parse(value)? as usize).min(3),
                    "shells" => self.sro.shells = parse(value)? as i32,
                    "tolerance" => self.sro.tolerance = parse(value)?,
                    "run" => {
                        let r = self.sro.run_now(&structure);
                        if !r.success {
                            return Err(r.message.clone());
                        }
                    }
                    other => return Err(format!("unknown sro command {other}")),
                }
            }
            "expect-sro" => {
                // expect-sro SHELL CENTRE NEIGHBOUR ALPHA: Warren-Cowley alpha (within 0.01), shell from 1.
                let r = self.sro.result.as_ref().ok_or("SRO has not run")?;
                let (shell, centre, neighbour) = (number(0).ok_or("expect-sro SHELL I J ALPHA")? as i32 - 1, words.get(1).copied().unwrap_or(""), words.get(2).copied().unwrap_or(""));
                let want = number(3).ok_or("expect-sro SHELL I J ALPHA")?;
                let e = r.warren_cowley.entries.iter().find(|e| e.shell == shell && e.center_element == centre && e.neighbor_element == neighbour).ok_or("no such SRO entry")?;
                if (e.alpha - want).abs() > 0.01 {
                    return Err(format!("alpha {centre}-{neighbour} in shell {} is {}, expected {want}", shell + 1, e.alpha));
                }
            }
            "trajectory" => {
                // trajectory open PATH | frame N (from 1) | analyse on|off | analysis N | lines on|off | gif PATH
                let (what, value) = split(rest);
                match what {
                    "open" => {
                        let s = self.trajectory.load_now(value.trim())?;
                        self.trajectory.open = true;
                        self.show_trajectory_frame(s, true, false);
                    }
                    "frame" => {
                        let n = value.trim().parse::<usize>().map_err(|e| format!("frame: {e}"))?;
                        if n == 0 || n > self.trajectory.count() {
                            return Err(format!("no frame {n} of {}", self.trajectory.count()));
                        }
                        self.trajectory.playback.frame = n - 1;
                        let s = self.trajectory.frame_structure(n - 1)?;
                        let property = self.trajectory.analyse;
                        self.show_trajectory_frame(s, false, property);
                    }
                    "analyse" => self.trajectory.analyse = value.trim() == "on",
                    "analysis" => self.trajectory.analysis = value.trim().parse::<usize>().map_err(|e| e.to_string())?.saturating_sub(1).min(3),
                    "lines" => {
                        // lines on|off|whole (whole: every frame's path, not up to the shown one).
                        self.trajectory.lines.visible = matches!(value.trim(), "on" | "whole");
                        self.trajectory.lines.whole_path = value.trim() == "whole";
                        if self.trajectory.lines.visible {
                            self.trajectory.trace_lines();
                        }
                    }
                    "gif" => {
                        self.trajectory.start_gif(value.trim());
                        while self.trajectory.exporting() {
                            self.gif_step();
                        }
                        if !self.trajectory.gif_message.starts_with("Saved") {
                            return Err(self.trajectory.gif_message.clone());
                        }
                    }
                    other => return Err(format!("unknown trajectory command {other}")),
                }
            }
            "plane" => {
                // plane H K L [OFFSET]: adds a lattice plane; plane clear removes them.
                if rest.trim() == "clear" {
                    self.overlays.planes.clear();
                    return Ok(());
                }
                let v: Vec<f64> = words.iter().filter_map(|w| w.parse().ok()).collect();
                if v.len() < 3 {
                    return Err("plane H K L [OFFSET]".into());
                }
                let p = crate::overlays::LatticePlane { h: v[0] as i32, k: v[1] as i32, l: v[2] as i32, offset: v.get(3).copied().unwrap_or(1.0) as f32, ..Default::default() };
                self.overlays.add_plane(p)?;
            }
            "direction" => {
                // direction U V W [LENGTH]: adds a Miller direction arrow; direction clear removes them.
                if rest.trim() == "clear" {
                    self.overlays.directions.clear();
                    return Ok(());
                }
                let v: Vec<f64> = words.iter().filter_map(|w| w.parse().ok()).collect();
                if v.len() < 3 {
                    return Err("direction U V W [LENGTH]".into());
                }
                let d = crate::overlays::MillerDirection { u: v[0] as i32, v: v[1] as i32, w: v[2] as i32, length: v.get(3).copied().unwrap_or(5.0) as f32, ..Default::default() };
                self.overlays.add_direction(d)?;
            }
            "voronoi" => {
                // voronoi on|off; expect-voronoi N checks the cell count.
                self.overlays.show_voronoi = rest.trim() == "on";
            }
            "expect-voronoi" => {
                let n = number(0).ok_or("expect-voronoi CELLS [VOLUME]")? as usize;
                let s = self.structure();
                let (diagram, volumes) = self.overlays.voronoi(&s);
                let total: f32 = volumes.iter().sum();
                if diagram.cells.len() != n || number(1).is_some_and(|v| (total as f64 - v).abs() > 0.01 * v) {
                    return Err(format!("{} Voronoi cells, total volume {total}", diagram.cells.len()));
                }
            }
            "expect-polyhedra" => {
                // expect-polyhedra N [FACES]: N polyhedra (each with FACES triangles).
                let d = self.doc().ok_or("no structure")?;
                let (structure, colours) = (d.structure.clone(), d.element_colours.clone());
                let p = self.polyhedra.polyhedra(&structure, &colours);
                let n = number(0).ok_or("expect-polyhedra N [FACES]")? as usize;
                if p.len() != n || number(1).is_some_and(|f| p.iter().any(|x| x.triangles.len() != f as usize)) {
                    return Err(format!("{} polyhedra with {:?} faces", p.len(), p.iter().map(|x| x.triangles.len()).collect::<Vec<_>>()));
                }
            }
            "polyhedra-centres" => {
                self.polyhedra.settings.centre_element_filter = true;
                self.polyhedra.settings.centre_elements = af_render::polyhedra::parse_elements(rest);
            }
            "expect-frames" => {
                let n = number(0).ok_or("expect-frames N")? as usize;
                if self.trajectory.count() != n {
                    return Err(format!("the trajectory has {} frames", self.trajectory.count()));
                }
            }
            "expect-cna" => {
                // expect-cna ENVIRONMENT N: N atoms in this environment.
                let r = self.cna.result.as_ref().ok_or("CNA has not run")?;
                let n = words.last().and_then(|w| w.parse::<i32>().ok()).ok_or("expect-cna ENVIRONMENT N")?;
                let name = rest.rsplit_once(' ').map(|(e, _)| e.trim()).unwrap_or("");
                let got = r.environment_counts.get(name).copied().unwrap_or(0);
                if got != n {
                    return Err(format!("{got} atoms are {name}; environments {:?}", r.environment_counts));
                }
            }
            "expect-adf-peak" => {
                // expect-adf-peak ANGLE: a detected peak within 1.5 degrees.
                let r = self.adf.result.as_ref().ok_or("ADF has not run")?;
                let want = number(0).ok_or("expect-adf-peak ANGLE")?;
                if !r.peaks.iter().any(|p| (p.angle_deg as f64 - want).abs() <= 1.5) {
                    let peaks: Vec<f32> = r.peaks.iter().map(|p| p.angle_deg).collect();
                    return Err(format!("no ADF peak near {want}; peaks at {peaks:?}"));
                }
            }
            "expect-rdf" => {
                // expect-rdf peak R | coordination C (within 0.01).
                let r = self.rdf.result.as_ref().ok_or("RDF has not run")?;
                let want = number(1).ok_or("expect-rdf peak|coordination VALUE")?;
                let got = match words.first().copied() {
                    Some("peak") => r.first_peak_r as f64,
                    Some("coordination") => r.shell_coordination as f64,
                    _ => return Err("expect-rdf peak|coordination VALUE".into()),
                };
                if (got - want).abs() > 0.01 {
                    return Err(format!("RDF {} is {got}, expected {want}", words[0]));
                }
            }
            "quit" => ctx.send_viewport_cmd(egui::ViewportCommand::Close),
            "fail" => return Err(rest.to_string()),
            other => return Err(format!("unknown command {other}")),
        }
        Ok(())
    }
}

/// Saves the next frame of the window as PNG (and closes the window after it).
fn save_screenshot(ctx: &egui::Context, path: String, close: bool) {
    let ctx2 = ctx.clone();
    ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::ScreenshotCallback::new(move |image| {
        let [w, h] = image.size;
        let bytes: Vec<u8> = image.pixels.iter().flat_map(|c| c.to_array()).collect();
        if let Some(parent) = std::path::Path::new(&path).parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".jpg") || lower.ends_with(".jpeg") {
            // JPEG at quality 92, as the C++ tutorial capture wrote them.
            let rgb: Vec<u8> = bytes.chunks_exact(4).flat_map(|p| [p[0], p[1], p[2]]).collect();
            if let Ok(file) = std::fs::File::create(&path) {
                let mut encoder = image::codecs::jpeg::JpegEncoder::new_with_quality(std::io::BufWriter::new(file), 92);
                let _ = encoder.encode(&rgb, w as u32, h as u32, image::ExtendedColorType::Rgb8);
            }
            if close {
                ctx2.send_viewport_cmd(egui::ViewportCommand::Close);
            }
            ctx2.request_repaint();
            return;
        }
        if let Ok(file) = std::fs::File::create(&path) {
            let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), w as u32, h as u32);
            encoder.set_color(png::ColorType::Rgba);
            encoder.set_depth(png::BitDepth::Eight);
            if let Ok(mut writer) = encoder.write_header() {
                let _ = writer.write_image_data(&bytes);
            }
        }
        if close {
            ctx2.send_viewport_cmd(egui::ViewportCommand::Close);
        }
        ctx2.request_repaint();
    })));
}


/// Manual screenshots (`capture-folder`, `capture`, `capture-done`).
#[derive(Default)]
struct Capture {
    folder: String,
    bounds: Vec<(String, [i32; 4])>,
}

/// The first letter in upper case (the History menu's labels).
fn capitalised(text: &str) -> String {
    let mut chars = text.chars();
    match chars.next() {
        Some(c) => c.to_uppercase().collect::<String>() + chars.as_str(),
        None => String::new(),
    }
}

/// The manual beside the program, in a package, or (development) in the sources (findManualPdf).
fn find_manual() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    let directory = exe.parent()?.to_path_buf();
    let parent = directory.parent().map(|p| p.to_path_buf()).unwrap_or_default();
    let candidates = [
        parent.join("docs/manual/AtomForge-manual.pdf"),
        parent.parent().map(|p| p.join("docs/manual/AtomForge-manual.pdf")).unwrap_or_default(),
        directory.join("AtomForge-manual.pdf"),
        parent.join("AtomForge-manual.pdf"),
        parent.join("share/doc/atomforge/AtomForge-manual.pdf"),
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../docs/manual/AtomForge-manual.pdf"),
    ];
    candidates.into_iter().find(|c| c.is_file())
}

/// The instance a pixel's ray hits first (handlePendingAtomPick): from the
/// camera in perspective, from the near plane in orthographic projection.
fn pick_instance(frame: &af_render::scene::FrameView, scene: &af_render::scene::Scene, orthographic: bool, mx: f64, my: f64) -> Option<usize> {
    let origin = if orthographic { pick_ray_origin(mx, my, frame.width, frame.height, &frame.projection, &frame.view) } else { frame.camera_position };
    let dir = pick_ray_dir(mx, my, frame.width, frame.height, &frame.projection, &frame.view);
    pick_atom(origin, dir, &scene.instances.positions, &scene.instances.scales, 1.0)
}

/// The grab overlay (drawGrabOverlay): the mode and its keys at the top, an
/// arrow from each selected atom's start to where it is, coordinates beside
/// the first eight.
#[allow(clippy::too_many_arguments)]
fn draw_grab(painter: &egui::Painter, rect: egui::Rect, scale: f32, grab: &af_render::grab::Grab, selection: &[usize], scene: &af_render::scene::Scene, structure: &Structure, frame: &af_render::scene::FrameView, light: bool) {
    use af_render::grab::Axis;
    let rgba = egui::Color32::from_rgba_unmultiplied;
    let (axis_colour, arrow) = match grab.axis {
        Axis::X => (rgba(255, 80, 80, 220), (255, 80, 80)),
        Axis::Y => (rgba(80, 255, 80, 220), (80, 255, 80)),
        Axis::Z => (rgba(80, 130, 255, 220), (80, 130, 255)),
        Axis::Free => (if light { rgba(30, 30, 30, 220) } else { rgba(255, 255, 255, 220) }, (255, 200, 50)),
    };
    let header = format!("Grab Mode{}  |  LMB/Enter: confirm  |  Esc/RMB: cancel  |  X/Y/Z: constrain", grab.axis.label());
    let font = egui::FontId::proportional(14.0);
    let galley = painter.layout_no_wrap(header, font.clone(), axis_colour);
    let at = egui::pos2(rect.center().x - galley.size().x / 2.0, rect.min.y + 8.0);
    painter.rect_filled(egui::Rect::from_min_size(at - egui::vec2(8.0, 4.0), galley.size() + egui::vec2(16.0, 8.0)), 4.0, if light { rgba(240, 240, 240, 210) } else { rgba(20, 20, 20, 200) });
    painter.galley(at, galley, axis_colour);
    let to_screen = |p: glam::Vec3| crate::view::project(p, frame).map(|(x, y)| rect.min + egui::vec2(x / scale, y / scale));
    for (k, &instance) in selection.iter().enumerate() {
        let (Some(&now), Some(&start)) = (scene.instances.positions.get(instance), grab.original_positions.get(k)) else { continue };
        let detailed = k < 8;
        let alpha = if detailed { 200 } else { 140 };
        let (Some(a), Some(b)) = (to_screen(start), to_screen(now)) else { continue };
        if (b - a).length() > 3.0 {
            let colour = rgba(arrow.0, arrow.1, arrow.2, alpha);
            painter.circle_stroke(a, if detailed { 5.0 } else { 4.0 }, egui::Stroke::new(1.5, if light { rgba(100, 100, 100, 150) } else { rgba(180, 180, 180, 150) }));
            painter.line_segment([a, b], egui::Stroke::new(if detailed { 2.5 } else { 1.5 }, colour));
            let n = (b - a).normalized();
            let side = egui::vec2(-n.y, n.x);
            let size = if detailed { 12.0 } else { 9.0 };
            painter.add(egui::Shape::convex_polygon(vec![b, b - n * size + side * size * 0.4, b - n * size - side * size * 0.4], colour, egui::Stroke::NONE));
        }
        if detailed {
            let mut text = format!("Cart: ({:.4}, {:.4}, {:.4})", now.x, now.y, now.z);
            if let Some((_, inv)) = af_core::glm::cell_matrices(structure) {
                let o = glam::Vec3::new(structure.cell_offset[0] as f32, structure.cell_offset[1] as f32, structure.cell_offset[2] as f32);
                let f = af_core::glm::mul(&inv, now - o);
                text += &format!("\nFrac: ({:.4}, {:.4}, {:.4})", f.x, f.y, f.z);
            }
            let galley = painter.layout_no_wrap(text, egui::FontId::monospace(12.0), if light { rgba(120, 100, 20, 255) } else { rgba(255, 255, 180, 255) });
            let at = b + egui::vec2(12.0, -10.0);
            painter.rect_filled(egui::Rect::from_min_size(at - egui::vec2(4.0, 4.0), galley.size() + egui::vec2(12.0, 8.0)), 3.0, if light { rgba(240, 240, 240, 210) } else { rgba(10, 10, 10, 200) });
            painter.galley(at, galley, egui::Color32::WHITE);
        }
    }
    if selection.len() > 8 {
        let text = format!("... and {} more atoms", selection.len() - 8);
        painter.text(egui::pos2(rect.center().x, rect.max.y - 30.0), egui::Align2::CENTER_TOP, text, font, if light { rgba(60, 60, 60, 200) } else { rgba(200, 200, 200, 200) });
    }
}

