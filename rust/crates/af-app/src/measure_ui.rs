//! The measurement result cards beside the view and the Atom Info window
//! (the C++ drawMeasurementPopups and drawMeasurementOverlays): the value
//! large, the atoms below, Copy and Close; dashed lines between the measured
//! atoms while a card is open.

use af_render::measure::Measurement;
use af_render::scene::{FrameView, Scene};
use eframe::egui;

use crate::view::project;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Distance,
    Angle,
    Dihedral,
}

impl Kind {
    pub fn title(self) -> &'static str {
        match self {
            Kind::Distance => "Distance",
            Kind::Angle => "Angle",
            Kind::Dihedral => "Dihedral",
        }
    }
}

#[derive(Default)]
pub struct MeasureCards {
    pub distance: Option<Measurement>,
    pub angle: Option<Measurement>,
    pub dihedral: Option<Measurement>,
    pub atom_info: Option<String>,
    /// The structure revision measured; lines are dropped when it changes.
    revision: u64,
}

impl MeasureCards {
    pub fn set(&mut self, kind: Kind, m: Measurement, revision: u64) {
        self.revision = revision;
        match kind {
            Kind::Distance => self.distance = Some(m),
            Kind::Angle => self.angle = Some(m),
            Kind::Dihedral => self.dihedral = Some(m),
        }
    }

    pub fn get(&self, kind: Kind) -> Option<&Measurement> {
        match kind {
            Kind::Distance => self.distance.as_ref(),
            Kind::Angle => self.angle.as_ref(),
            Kind::Dihedral => self.dihedral.as_ref(),
        }
    }

    /// The structure changed: the instances measured may be gone.
    pub fn structure_changed(&mut self, revision: u64) {
        if revision != self.revision {
            for m in [&mut self.distance, &mut self.angle].into_iter().flatten() {
                m.lines.clear();
            }
        }
    }

    /// A result card beside the view (drawResultCard: not modal, so the
    /// measured atoms stay in sight): the value large, the atoms below, Copy
    /// and Close. Cards stack downwards from `top`; returns the next top.
    fn card(ctx: &egui::Context, kind: Kind, slot: &mut Option<Measurement>, top: f32) -> f32 {
        let Some(m) = slot.as_ref() else { return top };
        let mut open = true;
        let mut close = false;
        let screen = ctx.content_rect();
        let margin = 16.0;
        let response = crate::theme::window(kind.title())
            .id(egui::Id::unique(("measurement", kind.title())))
            .open(&mut open)
            .resizable(false)
            .fixed_pos(egui::pos2(screen.max.x - margin, top))
            .pivot(egui::Align2::RIGHT_TOP)
            .resizable(false)
            .show(ctx, |ui| {
                // As narrow as the C++ card: 280 wide, wider only for a long text.
                let text_width = ui.painter().layout_no_wrap(m.message.clone(), egui::TextStyle::Body.resolve(ui.style()), egui::Color32::WHITE).size().x;
                ui.set_width(text_width.clamp(248.0, 388.0));
                if m.value.is_empty() {
                    crate::dialog_layout::note(ui, crate::dialog_layout::Note::Warning, &m.message);
                } else {
                    ui.label(egui::RichText::new(&m.value).font(crate::theme::heading(1.5)));
                    let muted = crate::theme::palette(ui.visuals().dark_mode).text_muted;
                    ui.label(egui::RichText::new(&m.message).color(muted));
                }
                match crate::dialog_layout::footer(ui, "Close", (!m.value.is_empty()).then_some("Copy"), true) {
                    // The number alone, for pasting into a spreadsheet or report.
                    2 => ui.ctx().copy_text(m.value.split(' ').next().unwrap_or("").to_string()),
                    1 => close = true,
                    _ => {}
                }
                if ui.ctx().memory(|mem| mem.focused().is_none()) && ui.ui_contains_pointer() && ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                    close = true;
                }
            });
        if !open || close {
            *slot = None;
        }
        response.map(|r| r.response.rect.max.y + margin * 0.5).unwrap_or(top)
    }

    pub fn show(&mut self, ctx: &egui::Context) {
        let mut top = ctx.content_rect().top() + crate::theme::FRAME_HEIGHT * 4.5 + 16.0;
        top = Self::card(ctx, Kind::Distance, &mut self.distance, top);
        top = Self::card(ctx, Kind::Angle, &mut self.angle, top);
        Self::card(ctx, Kind::Dihedral, &mut self.dihedral, top);
        if let Some(text) = self.atom_info.clone() {
            let mut open = true;
            let mut close = false;
            crate::theme::modal_window("Atom Info").open(&mut open).resizable(false).min_width(360.0).max_width(520.0).show(ctx, |ui| {
                ui.label(&text);
                match crate::dialog_layout::footer(ui, "Close", Some("Copy"), true) {
                    1 => close = true,
                    2 => ui.ctx().copy_text(text.clone()),
                    _ => {}
                }
            });
            if !open || close {
                self.atom_info = None;
            }
        }
    }

    /// Dashed lines between the measured atoms (and the angle's arc).
    pub fn draw(&self, painter: &egui::Painter, rect: egui::Rect, scale: f32, scene: &Scene, frame: &FrameView) {
        let point = |i: usize| -> Option<egui::Pos2> {
            let p = *scene.instances.positions.get(i)?;
            let (x, y) = project(p, frame)?;
            Some(rect.min + egui::vec2(x / scale, y / scale))
        };
        let dashed = |a: egui::Pos2, b: egui::Pos2, colour: egui::Color32, dash: f32, gap: f32| {
            let d = b - a;
            let length = d.length();
            if length < 1.0 {
                return;
            }
            let u = d / length;
            let mut t = 0.0;
            while t < length {
                let t2 = (t + dash).min(length);
                painter.line_segment([a + u * t, a + u * t2], egui::Stroke::new(2.0, colour));
                t += dash + gap;
            }
        };
        if let Some(m) = &self.distance {
            if let [a, b] = m.lines[..] {
                if let (Some(pa), Some(pb)) = (point(a), point(b)) {
                    dashed(pa, pb, egui::Color32::from_rgba_unmultiplied(255, 255, 80, 230), 7.0, 5.0);
                }
            }
        }
        if let Some(m) = &self.angle {
            if let [i0, i1, i2] = m.lines[..] {
                let colour = egui::Color32::from_rgba_unmultiplied(80, 220, 255, 230);
                let (p0, p1, p2) = (point(i0), point(i1), point(i2));
                if let (Some(a), Some(v)) = (p0, p1) {
                    dashed(a, v, colour, 7.0, 5.0);
                }
                if let (Some(c), Some(v)) = (p2, p1) {
                    dashed(c, v, colour, 7.0, 5.0);
                }
                if let (Some(a), Some(v), Some(c)) = (p0, p1, p2) {
                    let (d0, d2) = (a - v, c - v);
                    let r = d0.length().min(d2.length()).min(30.0) * 0.35;
                    if r > 4.0 {
                        let a0 = d0.y.atan2(d0.x);
                        let mut da = d2.y.atan2(d2.x) - a0;
                        while da > std::f32::consts::PI {
                            da -= 2.0 * std::f32::consts::PI;
                        }
                        while da < -std::f32::consts::PI {
                            da += 2.0 * std::f32::consts::PI;
                        }
                        let points: Vec<egui::Pos2> = (0..=20).map(|s| {
                            let t = a0 + da * s as f32 / 20.0;
                            v + egui::vec2(r * t.cos(), r * t.sin())
                        }).collect();
                        painter.add(egui::Shape::line(points, egui::Stroke::new(1.5, colour)));
                    }
                }
            }
        }
    }
}
