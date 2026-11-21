//! Analysis > Trajectories and Transport > Trajectory Playback (the C++
//! ui/TrajectoryDialog and util/TrajectoryPlayback): frames indexed in the
//! background and read on demand, transport controls and a frame slider,
//! a per-frame analysis coloured onto the atoms, unwrapped trajectory lines
//! and GIF export of the view.

use std::path::Path;
use std::sync::mpsc::{channel, Receiver};
use std::sync::Arc;

use af_core::structure::Structure;
use af_render::scene::FrameView;
use af_science::data::TrajectoryStream;
use eframe::egui;

/// Playback state (TrajectoryPlayback, kept free of the UI for tests).
#[derive(Clone, Copy, Debug, Default)]
pub struct Playback {
    pub playing: bool,
    pub frame: usize,
    pub previous: f64,
}

impl Playback {
    /// Applies this UI frame's events; true when `frame` should be shown.
    pub fn tick(&mut self, toggle_play: bool, jump_to_first: bool, slider_changed: bool, frames: usize, fps: f32, now: f64) -> bool {
        let mut changed = slider_changed;
        if toggle_play {
            self.playing = !self.playing;
            self.previous = now;
            changed = true;
        }
        if jump_to_first {
            self.frame = 0;
            self.playing = false;
            changed = true;
        }
        if self.playing && frames > 0 && now - self.previous >= 1.0 / fps as f64 {
            self.frame = (self.frame + 1) % frames;
            self.previous = now;
            changed = true;
        }
        changed
    }
}

/// Per-atom analyses that can colour every frame (catalog tool id, label).
pub const FRAME_ANALYSES: [(&str, &str); 4] = [
    ("structure-type", "Structure type (Ackland-Jones)"),
    ("diamond-structure", "Diamond structure (cubic, hexagonal)"),
    ("centrosymmetry", "Centrosymmetry parameter"),
    ("bond-order", "Steinhardt q6 bond order"),
];

/// What the dialog asks the app to do this frame.
#[derive(Default)]
pub struct Events {
    /// Show this frame (its structure, analysed when asked).
    pub show: Option<Structure>,
    /// A trajectory finished loading (its first frame is in `show`).
    pub loaded: bool,
    pub open_file: bool,
    pub export_gif: bool,
    /// Switch the view to atom-property colouring.
    pub colour_by_property: bool,
}

/// The unwrapped atom paths drawn over the view.
#[derive(Default)]
pub struct Lines {
    pub visible: bool,
    pub whole_path: bool,
    pub paths: Vec<Vec<[f32; 3]>>,
    /// The trajectory frame of each path point.
    pub frames: Vec<usize>,
}

struct GifExport {
    path: String,
    writer: Option<af_science::gif::GifWriter>,
    next: usize,
    restore: usize,
}

pub struct TrajectoryDialog {
    pub open: bool,
    pub stream: Option<Arc<TrajectoryStream>>,
    loading: Option<(String, Receiver<Result<TrajectoryStream, String>>)>,
    pub loaded: String,
    pub error: String,
    pub playback: Playback,
    pub fps: f32,
    lammps_types: String,
    pub analyse: bool,
    pub analysis: usize,
    pub cutoff: f32,
    analysis_message: String,
    pub line_stride: usize,
    pub lines: Lines,
    line_message: String,
    pub gif_stride: usize,
    gif: Option<GifExport>,
    pub gif_message: String,
    shown_atoms: usize,
}

impl Default for TrajectoryDialog {
    fn default() -> Self {
        Self {
            open: false,
            stream: None,
            loading: None,
            loaded: String::new(),
            error: String::new(),
            playback: Playback::default(),
            fps: 12.0,
            lammps_types: String::new(),
            analyse: false,
            analysis: 0,
            cutoff: 3.0,
            analysis_message: String::new(),
            line_stride: 1,
            lines: Lines::default(),
            line_message: String::new(),
            gif_stride: 1,
            gif: None,
            gif_message: String::new(),
            shown_atoms: 0,
        }
    }
}

fn reader(path: &Path) -> Result<Structure, String> {
    af_io::loader::load_structure_from_file(&path.to_string_lossy())
}

impl TrajectoryDialog {
    pub fn count(&self) -> usize {
        self.stream.as_ref().map(|s| s.len()).unwrap_or(0)
    }

    /// Indexes a trajectory file in the background (Open Trajectory).
    pub fn load(&mut self, path: &str) {
        self.playback.playing = false;
        self.error.clear();
        let (tx, rx) = channel();
        let p = path.to_string();
        std::thread::spawn(move || {
            let _ = tx.send(TrajectoryStream::new(Path::new(&p), Some(&reader)));
        });
        self.loading = Some((path.to_string(), rx));
    }

    /// Indexes on this thread (UI scripts).
    pub fn load_now(&mut self, path: &str) -> Result<Structure, String> {
        let stream = TrajectoryStream::new(Path::new(path), Some(&reader))?;
        self.stream = Some(Arc::new(stream));
        self.loaded = path.to_string();
        self.playback = Playback::default();
        if self.lines.visible {
            self.trace_lines();
        }
        self.frame_structure(0)
    }

    /// Frame `index` with the per-frame analysis when enabled.
    pub fn frame_structure(&mut self, index: usize) -> Result<Structure, String> {
        let stream = self.stream.as_ref().ok_or("No trajectory is loaded")?;
        let mut s = stream.frame(index)?.structure;
        self.shown_atoms = s.atoms.len();
        if self.analyse {
            let start = std::time::Instant::now();
            let (tool, label) = FRAME_ANALYSES[self.analysis.min(FRAME_ANALYSES.len() - 1)];
            match af_science::tools::analyse_per_atom(&s, tool, self.cutoff as f64, "") {
                Ok(p) => {
                    s.atom_property = p.values;
                    s.atom_property_name = p.name;
                    self.analysis_message = format!("{label}: {:.0} ms per frame", start.elapsed().as_secs_f64() * 1000.0);
                }
                Err(e) => {
                    self.analyse = false;
                    self.analysis_message = format!("Per-frame analysis stopped: {e}");
                }
            }
        }
        Ok(s)
    }

    /// Unwrapped paths of every Nth atom (at most 400 atoms, 4000 frames).
    pub fn trace_lines(&mut self) {
        self.lines.paths.clear();
        self.lines.frames.clear();
        self.line_message.clear();
        let Some(stream) = self.stream.clone().filter(|s| !s.is_empty()) else { return };
        let result = stream.frame(0).and_then(|f| {
            let atoms: Vec<usize> = (0..f.structure.atoms.len()).step_by(self.line_stride.max(1)).take(400).collect();
            let mut frames = Vec::new();
            let paths = af_science::trajectory_tools::unwrapped_paths(&stream, &atoms, stream.len().div_ceil(4000).max(1), &mut frames)?;
            Ok((atoms.len(), paths, frames))
        });
        match result {
            Ok((n, paths, frames)) => {
                self.lines.paths = paths.iter().map(|p| p.iter().map(|x| [x[0] as f32, x[1] as f32, x[2] as f32]).collect()).collect();
                self.line_message = format!("{n} atoms traced over {} frames", frames.len());
                self.lines.frames = frames;
            }
            Err(e) => {
                self.lines.visible = false;
                self.line_message = format!("Cannot trace: {e}");
            }
        }
    }

    pub fn exporting(&self) -> bool {
        self.gif.is_some()
    }

    /// Starts writing a GIF of every gif_stride-th frame.
    pub fn start_gif(&mut self, path: &str) {
        let mut path = path.to_string();
        if Path::new(&path).extension().is_none() {
            path += ".gif";
        }
        self.playback.playing = false;
        self.gif = Some(GifExport { path, writer: None, next: 0, restore: self.playback.frame });
        self.gif_message.clear();
    }

    /// The GIF frame to render next, if exporting.
    pub fn gif_frame(&self) -> Option<usize> {
        self.gif.as_ref().map(|g| g.next)
    }

    /// Adds the rendered view of the current GIF frame (RGBA, rows top to
    /// bottom); returns the frame to show afterwards when the export ends.
    pub fn add_gif_frame(&mut self, rgba: &[u8], width: usize, height: usize) -> Option<usize> {
        let count = self.count();
        let fps = self.fps;
        let stride = self.gif_stride.max(1);
        let g = self.gif.as_mut()?;
        let result = (|| -> Result<(), String> {
            // Box-downsample large views so the GIF stays a sensible size.
            let factor = width.max(height).div_ceil(960).max(1);
            let (w, h) = ((width / factor).max(1), (height / factor).max(1));
            let mut small = vec![0u8; w * h * 3];
            for y in 0..h {
                for x in 0..w {
                    for c in 0..3 {
                        let mut sum = 0usize;
                        for dy in 0..factor {
                            for dx in 0..factor {
                                sum += rgba[((y * factor + dy) * width + x * factor + dx) * 4 + c] as usize;
                            }
                        }
                        small[(y * w + x) * 3 + c] = (sum / (factor * factor)) as u8;
                    }
                }
            }
            if g.writer.is_none() {
                g.writer = Some(af_science::gif::GifWriter::new(Path::new(&g.path), w as i32, h as i32, ((100.0 / fps).round() as i32).max(2))?);
            }
            g.writer.as_mut().unwrap().add_frame(&small, 3)
        })();
        match result {
            Err(e) => Some(self.finish_gif(format!("GIF export failed: {e}"))),
            Ok(()) => {
                g.next += stride;
                if g.next >= count {
                    let n = g.writer.as_ref().map(|w| w.frames()).unwrap_or(0);
                    let path = g.path.clone();
                    Some(self.finish_gif(format!("Saved {n} frames to {path}")))
                } else {
                    None
                }
            }
        }
    }

    /// Ends the export; returns the frame to show again.
    pub fn finish_gif(&mut self, message: String) -> usize {
        let restore = match self.gif.take() {
            Some(mut g) => {
                if let Some(w) = g.writer.as_mut() {
                    let _ = w.close();
                }
                g.restore
            }
            None => self.playback.frame,
        };
        self.gif_message = message;
        self.playback.frame = restore.min(self.count().saturating_sub(1));
        self.playback.frame
    }

    pub fn show(&mut self, ctx: &egui::Context) -> Events {
        let mut ev = Events::default();
        let mut changed = false;
        if let Some((path, rx)) = &self.loading {
            match rx.try_recv() {
                Ok(Ok(stream)) => {
                    self.loaded = path.clone();
                    self.stream = Some(Arc::new(stream));
                    self.playback = Playback::default();
                    self.loading = None;
                    changed = true;
                    ev.loaded = true;
                    if self.lines.visible {
                        self.trace_lines();
                    }
                }
                Ok(Err(e)) => {
                    self.error = e;
                    self.loading = None;
                }
                Err(_) => ctx.request_repaint_after(std::time::Duration::from_millis(100)),
            }
        }
        let count = self.count();
        if !self.open {
            if self.playback.playing {
                self.playback.playing = false;
            }
        } else {
            let exporting = self.exporting();
            let mut open = self.open;
            crate::theme::window("Trajectory Playback").id(egui::Id::unique("trajectory playback")).open(&mut open).resizable(false).default_width(560.0).show(ctx, |ui| {
                use crate::dialog_widgets as w;
                use crate::theme::icon;
                ui.set_width(528.0);
                let full = 528.0;
                ui.add_enabled_ui(self.loading.is_none() && !exporting, |ui| {
                    ui.horizontal(|ui| {
                        ev.open_file |= ui.button("Open Trajectory").clicked();
                        let edit = ui.add_sized([220.0, crate::theme::FRAME_HEIGHT], egui::TextEdit::singleline(&mut self.lammps_types).hint_text("LAMMPS types, e.g. Cu Ni").align(egui::Align2::LEFT_CENTER));
                        let edit = edit.on_hover_text("Elements of LAMMPS atom types 1, 2, ... for a dump without an element column\n(the default dump atom); not needed with dump_modify ... element.");
                        if edit.changed() {
                            match af_io::trajectory::set_lammps_type_elements(&self.lammps_types) {
                                Ok(()) => self.error.clear(),
                                Err(e) => self.error = e,
                            }
                        }
                    });
                });
                if self.loading.is_some() {
                    ui.label("Indexing frames...");
                }
                if count > 0 {
                    ui.add_enabled_ui(!exporting, |ui| {
                        let focused = ui.ui_contains_pointer() && !ui.ctx().egui_wants_keyboard_input();
                        let (left, right, space) = ui.input(|i| (i.key_pressed(egui::Key::ArrowLeft), i.key_pressed(egui::Key::ArrowRight), i.key_pressed(egui::Key::Space)));
                        let mut step = 0i64;
                        let mut first = false;
                        let mut last = false;
                        let mut toggle = false;
                        // Transport: first, previous, play/pause, next, last.
                        ui.horizontal(|ui| {
                            first = ui.button(crate::theme::with_icon(icon::PREVIOUS, "First")).on_hover_text("First frame").clicked();
                            if ui.button("<").on_hover_text("Previous frame (Left)").clicked() || (focused && left) {
                                step = -1;
                            }
                            let play = if self.playback.playing { crate::theme::with_icon(icon::PAUSE, "Pause") } else { crate::theme::with_icon(icon::PLAY, "Play") };
                            toggle = w::primary_button(ui, &play, egui::Vec2::ZERO).on_hover_text("Play or pause (Space)").clicked() || (focused && space);
                            if ui.button(">").on_hover_text("Next frame (Right)").clicked() || (focused && right) {
                                step = 1;
                            }
                            last = ui.button(crate::theme::with_icon(icon::NEXT, "Last")).on_hover_text("Last frame").clicked();
                        });
                        let mut stepped = false;
                        if step != 0 || last {
                            self.playback.playing = false;
                            self.playback.frame = if last { count - 1 } else { (self.playback.frame as i64 + step).clamp(0, count as i64 - 1) as usize };
                            stepped = true;
                        }
                        let mut shown = self.playback.frame as f64 + 1.0;
                        let text = format!("Frame {} of {count}", self.playback.frame + 1);
                        let dragged = w::slider(ui, &mut shown, 1.0..=count as f64, full, &text, "", true).changed();
                        if dragged {
                            self.playback.frame = (shown as usize).saturating_sub(1).min(count - 1);
                        }
                        ui.label("Frames per second");
                        let mut fps = self.fps as f64;
                        let fps_text = format!("{:.0}", self.fps);
                        if w::slider(ui, &mut fps, 1.0..=60.0, full, &fps_text, "", false).changed() {
                            self.fps = fps as _;
                        }
                        let now = ui.input(|i| i.time);
                        changed |= self.playback.tick(toggle, first, dragged || stepped, count, self.fps, now);
                        if self.playback.playing {
                            ui.ctx().request_repaint();
                        }
                        let streamed = self.stream.as_ref().is_some_and(|s| s.streamed());
                        ui.label(format!("{count} frames{}, {} atoms in current frame", if streamed { " (read on demand)" } else { "" }, self.shown_atoms));
                        // Per-frame analysis: recomputed and shown as atom colours for every frame.
                        w::separator_text(ui, "Colour each frame");
                        if w::checkbox(ui, &mut self.analyse, "Analyse every displayed frame").changed() {
                            changed = true;
                            self.analysis_message.clear();
                        }
                        egui::ComboBox::from_id_salt("frame analysis").selected_text(FRAME_ANALYSES[self.analysis].1).width(320.0).show_ui(ui, |ui| {
                            for (k, (_, label)) in FRAME_ANALYSES.iter().enumerate() {
                                if ui.selectable_value(&mut self.analysis, k, *label).changed() {
                                    changed |= self.analyse;
                                }
                            }
                        });
                        let id = FRAME_ANALYSES[self.analysis].0;
                        if id != "structure-type" && id != "diamond-structure" {
                            if w::input_float(ui, &mut self.cutoff, 0.1, 2, 200.0, 0.5..=50.0, "Neighbour cutoff (Å)").changed() {
                                changed |= self.analyse;
                            }
                        }
                        if !self.analysis_message.is_empty() {
                            w::disabled(ui, &self.analysis_message);
                        }
                        // Trajectory lines: unwrapped atom paths drawn over the view.
                        w::separator_text(ui, "Trajectory lines");
                        let mut stride = self.line_stride as i64;
                        if w::input_int(ui, &mut stride, 200.0, 1..=100000, "Trace every Nth atom").changed() {
                            self.line_stride = stride.max(1) as usize;
                        }
                        ui.horizontal(|ui| {
                            if w::checkbox(ui, &mut self.lines.visible, "Show trajectory lines").changed() && self.lines.visible {
                                self.trace_lines();
                            }
                            w::checkbox(ui, &mut self.lines.whole_path, "Whole path");
                        });
                        if !self.line_message.is_empty() {
                            w::disabled(ui, &self.line_message);
                        }
                        w::separator_text(ui, "Animation");
                        ui.label("GIF: use every Nth frame");
                        let mut gif = self.gif_stride as f64;
                        let gif_text = self.gif_stride.to_string();
                        if w::slider(ui, &mut gif, 1.0..=((count / 2).max(1)) as f64, full, &gif_text, "", true).changed() {
                            self.gif_stride = (gif as usize).max(1);
                        }
                        ev.export_gif |= w::primary_button(ui, "Export GIF", egui::Vec2::ZERO).clicked();
                    });
                    if let Some(g) = &self.gif {
                        ui.label(format!("Writing GIF: frame {} of {count}", g.next + 1));
                        if ui.button("Stop Export").clicked() {
                            let n = g.writer.as_ref().map(|w| w.frames()).unwrap_or(0);
                            let message = format!("GIF export stopped; {n} frames written to {}", g.path);
                            self.finish_gif(message);
                            changed = true;
                        }
                    }
                    if !self.gif_message.is_empty() {
                        ui.add(egui::Label::new(&self.gif_message).wrap());
                    }
                }
                ui.add(egui::Label::new("Playback replaces the active structure. Frames are read from the file as they are shown, so long trajectories need little memory. Loaded trajectories can be analysed directly with the trajectory and transport tools. GIF export renders every chosen frame with the current camera and frame rate, at most 960 pixels wide.").wrap());
                if !self.error.is_empty() {
                    ui.colored_label(egui::Color32::from_rgb(220, 80, 60), &self.error);
                }
            });
            self.open = open;
        }
        if changed && count > 0 {
            match self.frame_structure(self.playback.frame) {
                Ok(s) => {
                    ev.colour_by_property = self.analyse;
                    ev.show = Some(s);
                }
                Err(e) => self.error = e,
            }
        }
        ev
    }

    /// The trajectory lines over the view, up to the shown frame (or whole).
    pub fn draw_lines(&self, painter: &egui::Painter, rect: egui::Rect, scale: f32, frame: &FrameView) {
        if !self.lines.visible || self.lines.paths.is_empty() {
            return;
        }
        let shown = self.lines.frames.iter().take_while(|&&f| self.lines.whole_path || f <= self.playback.frame).count();
        for (k, path) in self.lines.paths.iter().enumerate() {
            let points: Vec<egui::Pos2> = path
                .iter()
                .take(shown.max(1))
                .filter_map(|p| crate::view::project(glam::Vec3::from_array(*p), frame))
                .map(|(x, y)| rect.min + egui::vec2(x / scale, y / scale))
                .collect();
            let colour = crate::plot_widget::series_colour(k, false);
            painter.add(egui::Shape::line(points, egui::Stroke::new(1.5, colour)));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn playback_ticks_like_the_cpp() {
        let mut p = Playback::default();
        assert!(!p.tick(false, false, false, 10, 10.0, 0.0), "nothing happens");
        assert!(p.tick(true, false, false, 10, 10.0, 1.0) && p.playing, "play");
        assert!(!p.tick(false, false, false, 10, 10.0, 1.05), "not yet");
        assert!(p.tick(false, false, false, 10, 10.0, 1.11) && p.frame == 1, "next frame after 1/fps");
        p.frame = 9;
        assert!(p.tick(false, false, false, 10, 10.0, 1.3) && p.frame == 0, "wraps");
        assert!(p.tick(false, true, false, 10, 10.0, 1.31) && !p.playing && p.frame == 0, "first frame stops");
        assert!(p.tick(false, false, true, 10, 10.0, 2.0), "slider");
    }
}
