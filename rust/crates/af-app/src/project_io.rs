//! File > Open Project, Save Project As, autosave and Recover Last Autosave
//! (the C++ EditorApplication captureProject / restore and FileBrowser
//! workspace / restoreWorkspace): each tab's structure, camera, title, look
//! and view settings; the parts this app does not edit yet (volumetric data,
//! tables, tool and pipeline state) are kept as loaded and saved again.

use af_core::display::{look_text, read_look, DisplaySettings, Look};
use af_electronic::project::Workspace;
use af_render::colouring::ColourMode;
use glam::Vec3;

use crate::document::Document;
use crate::view::{AtomStyle, LabelContent, ViewOptions};

/// Where autosaves go: recovery.afproject in %LOCALAPPDATA%/AtomForge
/// (Windows) or the XDG state directory, the C++ desktop's file (File >
/// Recover Last Autosave reads projects of either).
pub fn autosave_path() -> Option<std::path::PathBuf> {
    let root = if cfg!(windows) {
        std::env::var_os("LOCALAPPDATA").map(std::path::PathBuf::from)
    } else {
        std::env::var_os("XDG_STATE_HOME").map(std::path::PathBuf::from).or_else(|| std::env::var_os("HOME").map(|h| std::path::PathBuf::from(h).join(".local").join("state")))
    }?;
    let dir = root.join("AtomForge");
    std::fs::create_dir_all(&dir).ok()?;
    Some(dir.join("recovery.afproject"))
}

fn colour_mode_index(m: ColourMode) -> f64 {
    match m {
        ColourMode::Element => 0.0,
        ColourMode::CrystalOrientation => 1.0,
        ColourMode::GrainBoundary => 2.0,
        ColourMode::AtomProperty => 3.0,
    }
}

fn style_index(s: AtomStyle) -> f64 {
    AtomStyle::ALL.iter().position(|&x| x == s).unwrap_or(0) as f64
}

/// A tab as a project tab (FileBrowser::workspace with the camera, title and look).
pub fn capture(doc: &Document, view: &ViewOptions, display: &DisplaySettings, light_theme: bool) -> Workspace {
    let mut w = doc.project.clone().unwrap_or_default();
    w.structure = doc.structure.clone();
    let s = &mut w.settings;
    let flag = |b: bool| if b { 1.0 } else { 0.0 };
    s.insert("desktop.useLightTheme".into(), flag(light_theme));
    s.insert("desktop.showElementLabels".into(), flag(view.show_labels));
    s.insert("desktop.atomLabelContent".into(), LabelContent::ALL.iter().position(|&c| c == view.label_content).unwrap_or(0) as f64);
    s.insert("desktop.showBonds".into(), flag(view.show_bonds));
    s.insert("desktop.showAtoms".into(), flag(view.show_atoms));
    s.insert("desktop.showBoundingBox".into(), flag(view.show_box));
    s.insert("desktop.viewMode".into(), flag(view.orthographic));
    s.insert("desktop.atomColorMode".into(), colour_mode_index(view.colour_mode));
    s.insert("desktop.atomDisplayMode".into(), style_index(view.style));
    s.insert("desktop.property.autoRange".into(), flag(view.property.auto_range));
    s.insert("desktop.property.low".into(), view.property.range.low);
    s.insert("desktop.property.high".into(), view.property.range.high);
    s.insert("desktop.property.hideOutside".into(), flag(view.property.hide_outside));
    s.insert("desktop.property.hideInvalid".into(), flag(view.property.hide_invalid));
    let c = &doc.camera;
    w.camera = vec![c.yaw as f64, c.pitch as f64, c.roll as f64, c.distance as f64, c.pan_offset.x as f64, c.pan_offset.y as f64, c.pan_offset.z as f64];
    w.title = doc.title.clone();
    w.display = look_text(&Look { display: *display, colours: doc.element_colours.clone(), radii: doc.element_radii.clone() });
    w
}

/// The view settings a project tab saved (restoreWorkspace), over `view`.
pub fn restore_view(w: &Workspace, view: &mut ViewOptions) {
    let read = |key: &str, fallback: f64| w.settings.get(key).copied().unwrap_or(fallback);
    let on = |key: &str, current: bool| read(key, if current { 1.0 } else { 0.0 }) != 0.0;
    view.show_labels = on("desktop.showElementLabels", view.show_labels);
    let content = read("desktop.atomLabelContent", 0.0).clamp(0.0, 3.0) as usize;
    view.label_content = LabelContent::ALL[content];
    view.show_bonds = on("desktop.showBonds", view.show_bonds);
    view.show_atoms = on("desktop.showAtoms", view.show_atoms);
    view.show_box = on("desktop.showBoundingBox", view.show_box);
    view.orthographic = read("desktop.viewMode", 0.0).clamp(0.0, 1.0) as i32 == 1;
    view.colour_mode = match read("desktop.atomColorMode", 0.0).clamp(0.0, 3.0) as i32 {
        1 => ColourMode::CrystalOrientation,
        2 => ColourMode::GrainBoundary,
        3 if w.structure.atom_property.len() == w.structure.atoms.len() => ColourMode::AtomProperty,
        _ => ColourMode::Element,
    };
    view.property.auto_range = on("desktop.property.autoRange", view.property.auto_range);
    view.property.range.low = read("desktop.property.low", view.property.range.low);
    view.property.range.high = read("desktop.property.high", view.property.range.high);
    view.property.hide_outside = on("desktop.property.hideOutside", view.property.hide_outside);
    view.property.hide_invalid = on("desktop.property.hideInvalid", view.property.hide_invalid);
    view.style = AtomStyle::ALL[read("desktop.atomDisplayMode", 0.0).clamp(0.0, 3.0) as usize];
}

/// A project tab as a document: its look (or `fallback`), camera and title;
/// the history starts at the restored state.
pub fn restore(w: Workspace, fallback: &Look) -> (Document, Option<DisplaySettings>) {
    let look = if w.display.is_empty() { None } else { Some(read_look(&w.display)) };
    let mut doc = Document::new(w.structure.clone(), None, look.as_ref().unwrap_or(fallback));
    doc.title = if w.title.is_empty() { "Restored project".into() } else { w.title.clone() };
    if w.camera.len() == 7 {
        let c = &w.camera;
        doc.camera.yaw = c[0] as f32;
        doc.camera.pitch = c[1] as f32;
        doc.camera.roll = c[2] as f32;
        doc.camera.distance = c[3].max(0.01) as f32;
        doc.camera.pan_offset = Vec3::new(c[4] as f32, c[5] as f32, c[6] as f32);
        doc.pending_default_view = false;
    }
    doc.reset_history();
    let display = look.map(|l| l.display);
    doc.project = Some(w);
    (doc, display)
}
