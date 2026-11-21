//! The desktop's large-data path (the C++ app/LargeData): atom clouds
//! (af_cloud) shown with the cloud renderer instead of the editable structure
//! renderer. Huge LAMMPS dumps and XYZ files are converted to a cloud in the
//! background on opening; the "Large Dataset" panel shows statistics and the
//! level-of-detail settings.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use af_render::cloud::CloudRenderer;
use af_render::scene::{perspective, radians, FrameView, OrbitCamera, Scene, MAX_DISTANCE, MIN_DISTANCE};
use eframe::egui;
use glam::Vec3;

/// Structures above this many atoms open as atom clouds.
pub const LARGE_DATA_ATOMS: u64 = 5_000_000;

fn extension(path: &Path) -> String {
    path.extension().map(|e| format!(".{}", e.to_string_lossy().to_ascii_lowercase())).unwrap_or_default()
}

pub fn is_lammps_dump(path: &Path) -> bool {
    let e = extension(path);
    if e == ".dump" || e == ".lammpstrj" {
        return true;
    }
    use std::io::BufRead;
    let Ok(file) = std::fs::File::open(path) else { return false };
    let mut first = String::new();
    let _ = std::io::BufReader::new(file).read_line(&mut first);
    first.starts_with("ITEM: TIMESTEP")
}

/// Whether a file opens through the large-data path: an .afcloud, a LAMMPS
/// dump, or an XYZ file of more than LARGE_DATA_ATOMS atoms (opensAsCloud).
pub fn opens_as_cloud(path: &Path) -> bool {
    if !path.is_file() {
        return false;
    }
    if af_cloud::is_cloud_file(path) || is_lammps_dump(path) {
        return true;
    }
    let e = extension(path);
    (e == ".xyz" || e == ".extxyz") && af_cloud::peek_atom_count(path) > LARGE_DATA_ATOMS
}

/// "2.50 billion", "12.3 million" or the number.
pub fn count_text(n: u64) -> String {
    if n >= 1_000_000_000 {
        format!("{:.2} billion", n as f64 / 1e9)
    } else if n >= 1_000_000 {
        format!("{:.1} million", n as f64 / 1e6)
    } else {
        n.to_string()
    }
}

// ---------------------------------------------------------------- conversion

#[derive(Default)]
struct Shared {
    progress: f64,
    message: String,
    error: String,
}

/// Background conversion of a large file to an atom cloud (CloudConversion).
#[derive(Default)]
pub struct CloudConversion {
    input: PathBuf,
    output: PathBuf,
    thread: Option<JoinHandle<()>>,
    finished: Arc<AtomicBool>,
    cancel: Arc<AtomicBool>,
    shared: Arc<Mutex<Shared>>,
}

/// Converts `input` to `output` (through a partial file, removed on failure).
fn convert(input: &Path, output: &Path, progress: af_cloud::Progress) -> Result<(), String> {
    let temporary = output.with_extension("afcloud-part");
    let built = if is_lammps_dump(input) {
        af_cloud::build_from_lammps_dump(input, &temporary, &[], 0, progress, af_cloud::BuildOptions::default())
    } else {
        af_cloud::build_from_xyz(input, &temporary, 0, progress, af_cloud::BuildOptions::default())
    };
    if let Err(e) = built {
        let _ = std::fs::remove_file(&temporary);
        return Err(e);
    }
    std::fs::rename(&temporary, output).map_err(|e| format!("Cannot write {}: {e}", output.display()))
}

/// An up-to-date cloud built from `input` before.
fn reusable(input: &Path, output: &Path) -> bool {
    let modified = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
    af_cloud::is_cloud_file(output) && modified(output) >= modified(input) && af_cloud::CloudFile::open(output).is_ok()
}

impl CloudConversion {
    /// Where the cloud of `input` is kept: next to it when that folder is
    /// writable, else in the user's cache folder.
    pub fn cache_path(input: &Path) -> PathBuf {
        let mut beside = input.as_os_str().to_owned();
        beside.push(".afcloud");
        let beside = PathBuf::from(beside);
        let mut probe = beside.as_os_str().to_owned();
        probe.push(".probe");
        if std::fs::File::create(&probe).is_ok() {
            let _ = std::fs::remove_file(&probe);
            return beside;
        }
        let root = std::env::var_os("LOCALAPPDATA").map(PathBuf::from).unwrap_or_else(std::env::temp_dir).join("AtomForge").join("clouds");
        let _ = std::fs::create_dir_all(&root);
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        input.to_string_lossy().hash(&mut h);
        root.join(format!("{}-{}.afcloud", input.file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default(), h.finish()))
    }

    /// Converts on this thread (UI scripts); returns the cloud.
    pub fn convert_now(input: &Path) -> Result<PathBuf, String> {
        let output = Self::cache_path(input);
        if !reusable(input, &output) {
            convert(input, &output, None)?;
        }
        Ok(output)
    }

    /// Starts converting `input` (an existing up-to-date cloud is reused at
    /// once); a conversion still running is cancelled.
    pub fn start(&mut self, input: &Path) {
        self.cancel.store(true, Ordering::Release);
        if let Some(t) = self.thread.take() {
            let _ = t.join();
        }
        self.input = input.to_path_buf();
        self.output = Self::cache_path(input);
        self.finished = Arc::new(AtomicBool::new(false));
        self.cancel = Arc::new(AtomicBool::new(false));
        self.shared = Arc::new(Mutex::new(Shared { message: "Preparing".into(), ..Default::default() }));
        if reusable(input, &self.output) {
            self.finished.store(true, Ordering::Release);
            self.thread = Some(std::thread::spawn(|| {}));
            return;
        }
        let (input, output) = (self.input.clone(), self.output.clone());
        let (finished, cancel, shared) = (self.finished.clone(), self.cancel.clone(), self.shared.clone());
        self.thread = Some(std::thread::spawn(move || {
            let report = |fraction: f64, message: &str| -> bool {
                let mut s = shared.lock().unwrap();
                if fraction >= 0.0 {
                    s.progress = fraction;
                }
                s.message = message.to_string();
                !cancel.load(Ordering::Acquire)
            };
            if let Err(e) = convert(&input, &output, Some(&report)) {
                shared.lock().unwrap().error = e;
            }
            finished.store(true, Ordering::Release);
        }));
    }

    pub fn running(&self) -> bool {
        self.thread.is_some() && !self.finished.load(Ordering::Acquire)
    }

    /// When the conversion ended (once): the cloud to open, or why not.
    pub fn finished(&mut self) -> Option<Result<PathBuf, String>> {
        if !self.finished.load(Ordering::Acquire) {
            return None;
        }
        let t = self.thread.take()?;
        let _ = t.join();
        let error = std::mem::take(&mut self.shared.lock().unwrap().error);
        Some(if error.is_empty() { Ok(self.output.clone()) } else { Err(error) })
    }

    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Release);
    }

    /// Progress window of a running conversion (drawConversionProgress).
    pub fn show(&self, ctx: &egui::Context) {
        if !self.running() {
            return;
        }
        let (progress, message) = {
            let s = self.shared.lock().unwrap();
            (s.progress, s.message.clone())
        };
        ctx.request_repaint_after(std::time::Duration::from_millis(100));
        crate::theme::window("Preparing a large dataset")
            .id(egui::Id::unique("cloud conversion"))
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .fixed_size([480.0, 0.0])
            .show(ctx, |ui| {
                ui.label(format!(
                    "{} is too large to edit; AtomForge converts it to an atom cloud (once) and shows it with the large-data renderer.",
                    self.input.file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default()
                ));
                ui.add(egui::ProgressBar::new(progress as f32).show_percentage());
                ui.weak(message);
                if ui.button("Cancel").clicked() {
                    self.cancel();
                }
            });
    }
}

impl Drop for CloudConversion {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Release);
        if let Some(t) = self.thread.take() {
            let _ = t.join();
        }
    }
}

// ---------------------------------------------------------------- view

/// The scene of a cloud tab: orbit centre and outline, so that the camera's
/// framing, near/far planes and picking rays see the cloud's extent.
pub fn cloud_scene(cloud: &CloudRenderer) -> Scene {
    let mut scene = Scene { orbit_center: cloud.centre(), ..Default::default() };
    if let Some(info) = cloud.info() {
        scene.instances.box_lines = af_render::cloud::outline(info).into_iter().map(Vec3::from_array).collect();
    }
    scene
}

/// Fits the camera to a cloud (iso view; fitCameraToCloud).
pub fn fit_camera_to_cloud(camera: &mut OrbitCamera, radius: f32, width: i32, height: i32) {
    camera.yaw = 45.0;
    camera.pitch = 35.2643897;
    camera.roll = 0.0;
    camera.pan_offset = Vec3::ZERO;
    let aspect = if height > 0 { width as f32 / height as f32 } else { 1.0 };
    let vertical = radians(45.0);
    let horizontal = 2.0 * ((vertical * 0.5).tan() * aspect).atan();
    let half = radians(10.0).max(0.5 * vertical.min(horizontal));
    camera.distance = (radius / half.sin() * 1.1).clamp(MIN_DISTANCE, MAX_DISTANCE);
}

/// Perspective clip planes for a cloud: a sensible near plane inside huge
/// scenes (adjustCloudProjection). Orthographic views keep their projection.
pub fn adjust_cloud_projection(frame: &mut FrameView, camera: &OrbitCamera, centre: Vec3, radius: f32, orthographic: bool) {
    if orthographic || frame.height <= 0 {
        return;
    }
    let aspect = frame.width as f32 / frame.height as f32;
    let to_centre = (frame.camera_position - centre).length();
    // Near plane: well in front of the nearest atoms, never so close that depth precision is lost.
    let nearest = to_centre - radius;
    let near = if nearest > 0.0 { 0.5f32.max(nearest * 0.5) } else { 0.5f32.max(camera.distance * 0.002) };
    let far = (near * 2.0).max(to_centre + radius * 1.05);
    frame.projection = perspective(radians(45.0), aspect, near, far);
}

/// The "Large Dataset" panel: statistics and level-of-detail settings (drawCloudPanel).
pub fn cloud_panel(ctx: &egui::Context, cloud: &mut CloudRenderer, frame_milliseconds: f32) {
    let stats = *cloud.stats();
    let Some(info) = cloud.info().cloned() else { return };
    let screen = ctx.content_rect();
    crate::theme::window("Large Dataset")
        .id(egui::Id::unique("large dataset"))
        .collapsible(false)
        .default_pos(egui::pos2(12.0, screen.max.y - 380.0))
        .default_width(430.0)
        .show(ctx, |ui| {
            ui.label(cloud.path().file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default());
            if !info.source.is_empty() {
                ui.weak(&info.source);
            }
            ui.label(format!("{} atoms in {} chunks", count_text(stats.total_atoms), stats.chunks));
            ui.weak(format!("Drawn: {} atoms, {} chunks visible, {} hidden", count_text(stats.drawn_atoms), stats.visible_chunks, stats.occluded_chunks));
            ui.weak(format!("Frame {frame_milliseconds:.1} ms, {:.0} MB resident", stats.resident_mb));
            if !stats.complete {
                ui.colored_label(egui::Color32::from_rgb(230, 170, 60), format!("Streaming ({} loads pending)", stats.loading_chunks));
            } else {
                ui.colored_label(egui::Color32::from_rgb(70, 170, 90), "Up to date");
            }
            ui.separator();
            ui.strong("Rendering");
            let s = &mut cloud.settings;
            let mut budget = (s.point_budget / 1e6) as f32;
            if ui.add(egui::Slider::new(&mut budget, 1.0..=500.0).logarithmic(true).fixed_decimals(0).text("Atoms per frame (M)")).on_hover_text("The most atoms drawn in one frame; lower it for a faster view.").changed() {
                s.point_budget = budget as f64 * 1e6;
            }
            let mut memory = s.memory_budget_mb as f32;
            if ui
                .add(egui::Slider::new(&mut memory, 256.0..=16384.0).logarithmic(true).fixed_decimals(0).text("GPU memory (MB)"))
                .on_hover_text("Graphics memory for loaded chunks; the least recently seen are released beyond it.")
                .changed()
            {
                s.memory_budget_mb = memory as f64;
            }
            crate::dialog_widgets::checkbox(ui, &mut s.level_of_detail, "Level of detail").on_hover_text("Draw distant regions with a random subsample of their atoms (as many as the screen resolves).");
            ui.add_enabled(s.level_of_detail, egui::Slider::new(&mut s.density, 0.5..=20.0).logarithmic(true).fixed_decimals(1).text("Detail"))
                .on_hover_text("Atoms drawn per screen area at a distance; higher is denser and slower.");
            ui.add(egui::Slider::new(&mut s.radius_scale, 0.2..=3.0).fixed_decimals(2).text("Atom size")).on_hover_text("Scales the drawn atom radius.");
            const COLOURS: [&str; 4] = ["Species", "Height along x", "Height along y", "Height along z"];
            egui::ComboBox::from_label("Colour").selected_text(COLOURS[cloud.colour_axis.clamp(0, 3) as usize]).show_ui(ui, |ui| {
                for (k, name) in COLOURS.iter().enumerate() {
                    ui.selectable_value(&mut cloud.colour_axis, k as i32, *name);
                }
            });
            let s = &mut cloud.settings;
            ui.horizontal(|ui| {
                crate::dialog_widgets::checkbox(ui, &mut s.show_box, "Show box");
                crate::dialog_widgets::checkbox(ui, &mut s.occlusion_culling, "Occlusion culling").on_hover_text("Skip chunks hidden behind nearer ones.");
            });
            ui.add_space(4.0);
            ui.weak("Large datasets are view-only; File > Open a smaller file to edit.");
        });
}

/// Sets a panel setting by name (UI scripts): budget (M atoms), memory (MB),
/// lod on|off, detail, size, colour species|x|y|z, box on|off, occlusion on|off.
pub fn set_cloud_setting(cloud: &mut CloudRenderer, name: &str, value: &str) -> Result<(), String> {
    let on = value.trim() == "on";
    let number = || value.trim().parse::<f64>().map_err(|_| format!("cloud-setting {name} wants a number"));
    let s = &mut cloud.settings;
    match name {
        "budget" => s.point_budget = number()? * 1e6,
        "memory" => s.memory_budget_mb = number()?,
        "lod" => s.level_of_detail = on,
        "detail" => s.density = number()? as f32,
        "size" => s.radius_scale = number()? as f32,
        "box" => s.show_box = on,
        "occlusion" => s.occlusion_culling = on,
        "colour" => {
            cloud.colour_axis = match value.trim() {
                "species" => 0,
                "x" => 1,
                "y" => 2,
                "z" => 3,
                other => return Err(format!("no colour {other} (species, x, y or z)")),
            }
        }
        other => return Err(format!("no cloud setting {other}")),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_and_dumps() {
        assert_eq!(count_text(2_500_000_000), "2.50 billion");
        assert_eq!(count_text(12_345_678), "12.3 million");
        assert_eq!(count_text(999_999), "999999");
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/cloud");
        assert!(is_lammps_dump(&root.join("elements.dump")));
        assert!(opens_as_cloud(&root.join("elements.dump")));
        assert!(!opens_as_cloud(&root.join("missing.xyz")));
        let cu = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../docs/manual/examples/cu_fcc.cif");
        assert!(!opens_as_cloud(&cu));
    }

    #[test]
    fn fits_the_camera() {
        let mut c = OrbitCamera::default();
        fit_camera_to_cloud(&mut c, 100.0, 800, 600);
        // 1.1 r / sin(22.5 degrees) for a wide view.
        assert!((c.distance - 110.0 / radians(22.5).sin()).abs() < 1e-3, "{}", c.distance);
        assert_eq!((c.yaw, c.roll), (45.0, 0.0));
    }
}
