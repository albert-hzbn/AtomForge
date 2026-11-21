//! The AtomForge desktop application (egui/eframe on wgpu): structures in
//! tabs, the menus of the C++ desktop, and the 3D viewport, which af-render
//! draws on the UI's own GPU device and egui shows as a texture (no copy
//! through the CPU). Paths on the command line open at start.

mod adf_dialog;
mod analysis_panels;
mod dislocation_analysis_dialog;
mod electronic_dialog;
mod electronic_slice;
mod lobster_dialog;
mod wannier_dialog;
mod amorphous_dialog;
mod app;
mod stacking_fault_dialog;
mod interface_dialog;
mod csl_dialog;
mod atom_menu;
mod builder_preview;
mod bulk_dialog;
mod cna_dialog;
mod custom_dialog;
mod dialog_layout;
mod dialog_widgets;
mod dislocation_dialog;
mod display_settings;
mod edit_structure_dialog;
mod document;
mod file_browser;
mod image_export;
mod interstitial_dialog;
mod large_data;
mod measure_ui;
mod merge_dialog;
mod nano_dialog;
mod operation;
mod overlays;
mod palette;
mod pipeline_ui;
mod plot_widget;
mod params_form;
mod poly_dialog;
mod polyhedra_ui;
mod project_io;
mod property_ui;
mod rdf_dialog;
mod sculptor_dialog;
mod shape_dialog;
mod sro_dialog;
mod sss_dialog;
mod structure_dialog;
mod renderer_pref;
mod structure_info;
mod theme;
mod welcome;
mod tools;
mod trajectory_dialog;
mod transform_dialog;
mod ui_script;
mod undo;
mod view;

/// Where every pipeline step appears in the menus (shared with the command line).
mod menu_table {
    pub use af_cli::menu::*;
}

use eframe::egui;

/// A path with * or ? becomes the files it matches (sorted); without a
/// match, or without wildcards, it is kept as given (expandWildcard).
fn expand_wildcard(pattern: &str, out: &mut Vec<String>) {
    if !pattern.contains(['*', '?']) {
        out.push(pattern.to_string());
        return;
    }
    let path = std::path::Path::new(pattern);
    let directory = path.parent().filter(|p| !p.as_os_str().is_empty());
    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let listing = std::fs::read_dir(directory.unwrap_or(std::path::Path::new(".")));
    let mut matches: Vec<String> = match listing {
        Ok(entries) => entries
            .filter_map(|e| e.ok())
            .filter(|e| wildcard_match(&name, &e.file_name().to_string_lossy()))
            .map(|e| match directory {
                Some(d) => d.join(e.file_name()).to_string_lossy().into_owned(),
                None => e.file_name().to_string_lossy().into_owned(),
            })
            .collect(),
        Err(_) => Vec::new(),
    };
    if matches.is_empty() {
        out.push(pattern.to_string());
        return;
    }
    matches.sort();
    out.extend(matches);
}

/// * matches any run of characters, ? one character (case-insensitive on Windows).
fn wildcard_match(pattern: &str, name: &str) -> bool {
    let fold = |s: &str| -> Vec<char> { if cfg!(windows) { s.to_lowercase().chars().collect() } else { s.chars().collect() } };
    let (p, n) = (fold(pattern), fold(name));
    let (mut pi, mut ni, mut star, mut mark) = (0usize, 0usize, None::<usize>, 0usize);
    while ni < n.len() {
        if pi < p.len() && (p[pi] == '?' || p[pi] == n[ni]) {
            pi += 1;
            ni += 1;
        } else if pi < p.len() && p[pi] == '*' {
            star = Some(pi);
            pi += 1;
            mark = ni;
        } else if let Some(s) = star {
            pi = s + 1;
            mark += 1;
            ni = mark;
        } else {
            return false;
        }
    }
    p[pi..].iter().all(|&c| c == '*')
}

/// The window's icon: the largest PNG image of assets/icon/atomforge.ico.
fn window_icon(viewport: egui::ViewportBuilder) -> egui::ViewportBuilder {
    const ICO: &[u8] = include_bytes!("../../../../assets/icon/atomforge.ico");
    let entry = |i: usize| -> Option<(u32, &[u8])> {
        let e = ICO.get(6 + 16 * i..22 + 16 * i)?;
        let width = if e[0] == 0 { 256 } else { e[0] as u32 };
        let size = u32::from_le_bytes(e[8..12].try_into().ok()?) as usize;
        let offset = u32::from_le_bytes(e[12..16].try_into().ok()?) as usize;
        Some((width, ICO.get(offset..offset + size)?))
    };
    let count = ICO.get(4..6).map(|b| u16::from_le_bytes([b[0], b[1]]) as usize).unwrap_or(0);
    let Some((_, bytes)) = (0..count).filter_map(entry).filter(|(_, b)| b.starts_with(&[0x89, b'P', b'N', b'G'])).max_by_key(|(w, _)| *w) else { return viewport };
    let mut decoder = png::Decoder::new(std::io::Cursor::new(bytes));
    decoder.set_transformations(png::Transformations::EXPAND | png::Transformations::ALPHA);
    let Ok(mut reader) = decoder.read_info() else { return viewport };
    let mut buffer = vec![0; reader.output_buffer_size().unwrap_or(0)];
    let Ok(info) = reader.next_frame(&mut buffer) else { return viewport };
    if info.color_type != png::ColorType::Rgba || info.bit_depth != png::BitDepth::Eight {
        return viewport;
    }
    buffer.truncate(info.buffer_size());
    viewport.with_icon(std::sync::Arc::new(egui::IconData { rgba: buffer, width: info.width, height: info.height }))
}

fn main() -> eframe::Result {
    // Build and Edit operations are pipeline steps, as in the C++ desktop.
    af_cli::build_modifiers::register();
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    // One program, as the C++ AtomForge: command-line modes run without a window.
    let cli = af_cli::args::Args(arguments.clone());
    if af_cli::is_cli_mode(&cli) {
        std::process::exit(af_cli::run_cli(&cli));
    }
    let mut startup = Vec::new();
    let mut command_line_renderer = None;
    let mut i = 0;
    while i < arguments.len() {
        let a = &arguments[i];
        if a == "--renderer" && i + 1 < arguments.len() {
            command_line_renderer = Some(arguments[i + 1].clone());
            i += 2;
            continue;
        }
        if let Some(name) = a.strip_prefix("--renderer=") {
            command_line_renderer = Some(name.to_string());
            i += 1;
            continue;
        }
        expand_wildcard(a, &mut startup);
        i += 1;
    }
    // --renderer, else ATOMFORGE_RENDERER, else Settings > Renderer.
    let backends = renderer_pref::backends(renderer_pref::requested(command_line_renderer.as_deref()));
    let mut wgpu_options = eframe::egui_wgpu::WgpuConfiguration::default();
    if let (Some(b), eframe::egui_wgpu::WgpuSetup::CreateNew(setup)) = (backends, &mut wgpu_options.wgpu_setup) {
        setup.instance_descriptor.backends = b;
    }
    let options = eframe::NativeOptions {
        viewport: window_icon(egui::ViewportBuilder::default().with_inner_size([1280.0, 800.0]).with_maximized(true).with_title("AtomForge").with_drag_and_drop(true)),
        renderer: eframe::Renderer::Wgpu,
        wgpu_options,
        ..Default::default()
    };
    let result = eframe::run_native("AtomForge", options, Box::new(move |cc| Ok(Box::new(app::AtomForgeApp::new(cc, startup)))));
    if app::SCRIPT_FAILED.load(std::sync::atomic::Ordering::Relaxed) {
        std::process::exit(3);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wildcards() {
        assert!(wildcard_match("*.cif", "cu_fcc.cif"));
        assert!(wildcard_match("cu_???.cif", "cu_fcc.cif"));
        assert!(!wildcard_match("*.xyz", "cu_fcc.cif"));
        assert!(wildcard_match("*", ""));
        let mut out = Vec::new();
        expand_wildcard("no-such-dir/*.cif", &mut out);
        assert_eq!(out, vec!["no-such-dir/*.cif"]);
    }
}
