//! View > Show > Lattice Planes and Miller Directions (the C++
//! ui/LatticePlaneOverlay, ui/MillerDirectionOverlay and their settings
//! windows): (hkl) planes clipped to the unit cell as translucent polygons,
//! and [uvw] arrows from the cell origin.

use af_core::glm;
use af_core::structure::Structure;
use af_render::scene::FrameView;
use eframe::egui;
use glam::{Vec3, Vec4};

#[derive(Clone, Debug, PartialEq)]
pub struct LatticePlane {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    /// The plane h x + k y + l z = offset (fractional coordinates).
    pub offset: f32,
    pub opacity: f32,
    pub colour: [f32; 3],
    pub visible: bool,
}

impl Default for LatticePlane {
    fn default() -> Self {
        Self { h: 1, k: 0, l: 0, offset: 1.0, opacity: 0.20, colour: [0.95, 0.62, 0.20], visible: true }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MillerDirection {
    pub u: i32,
    pub v: i32,
    pub w: i32,
    /// Display length (Å).
    pub length: f32,
    pub colour: [f32; 3],
    pub visible: bool,
}

impl Default for MillerDirection {
    fn default() -> Self {
        Self { u: 1, v: 0, w: 0, length: 5.0, colour: [0.25, 0.75, 1.0], visible: true }
    }
}

const CUBE: [Vec3; 8] = [
    Vec3::new(0.0, 0.0, 0.0),
    Vec3::new(1.0, 0.0, 0.0),
    Vec3::new(0.0, 1.0, 0.0),
    Vec3::new(1.0, 1.0, 0.0),
    Vec3::new(0.0, 0.0, 1.0),
    Vec3::new(1.0, 0.0, 1.0),
    Vec3::new(0.0, 1.0, 1.0),
    Vec3::new(1.0, 1.0, 1.0),
];
const EDGES: [(usize, usize); 12] = [(0, 1), (1, 3), (3, 2), (2, 0), (4, 5), (5, 7), (7, 6), (6, 4), (0, 4), (1, 5), (2, 6), (3, 7)];

fn add_unique(points: &mut Vec<Vec3>, p: Vec3) {
    if !points.iter().any(|q| glm::dot(*q - p, *q - p) <= 1e-8) {
        points.push(p);
    }
}

/// The plane's polygon inside the unit cube (fractional coordinates), in
/// order around its centre (computePlanePolygonFractional); empty when it
/// misses the cell.
pub fn plane_polygon(h: i32, k: i32, l: i32, offset: f32) -> Vec<Vec3> {
    let n = Vec3::new(h as f32, k as f32, l as f32);
    let eps = 1e-6;
    let mut points = Vec::new();
    for (a, b) in EDGES {
        let (p0, p1) = (CUBE[a], CUBE[b]);
        let (f0, f1) = (glm::dot(n, p0) - offset, glm::dot(n, p1) - offset);
        if f0.abs() <= eps && f1.abs() <= eps {
            add_unique(&mut points, p0);
            add_unique(&mut points, p1);
            continue;
        }
        let denom = f0 - f1;
        if denom.abs() <= eps {
            continue;
        }
        let t = f0 / denom;
        if t < -eps || t > 1.0 + eps {
            continue;
        }
        add_unique(&mut points, p0 + (p1 - p0) * t.clamp(0.0, 1.0));
    }
    if points.len() < 3 {
        return Vec::new();
    }
    let centre = points.iter().fold(Vec3::ZERO, |a, &p| a + p) / points.len() as f32;
    let normal = glm::normalize(n);
    let reference = if normal.z.abs() < 0.8 { Vec3::Z } else { Vec3::Y };
    let u = glm::normalize(glm::cross(reference, normal));
    let v = glm::cross(normal, u);
    let angle = |p: &Vec3| {
        let d = *p - centre;
        glm::dot(d, v).atan2(glm::dot(d, u))
    };
    af_core::cxx::sort_by(&mut points, |a, b| angle(a) < angle(b));
    points
}

/// The part of a world polygon in front of the near plane (view depth).
fn clip_near(points: &[Vec3], frame: &FrameView, near: f32) -> Vec<Vec3> {
    let z = |p: Vec3| (frame.view * Vec4::new(p.x, p.y, p.z, 1.0)).z;
    let mut out = Vec::new();
    for i in 0..points.len() {
        let (a, b) = (points[i], points[(i + 1) % points.len()]);
        let (za, zb) = (z(a), z(b));
        let (ia, ib) = (za <= -near, zb <= -near);
        if ia {
            out.push(a);
        }
        if ia != ib && (zb - za).abs() > 1e-6 {
            out.push(a + (b - a) * ((-near - za) / (zb - za)).clamp(0.0, 1.0));
        }
    }
    if out.len() >= 3 {
        out
    } else {
        Vec::new()
    }
}

fn colour(c: [f32; 3], alpha: f32) -> egui::Color32 {
    crate::theme::srgba(c[0], c[1], c[2], alpha).into()
}

fn cell_frame(s: &Structure) -> Option<(glam::Mat3, Vec3)> {
    let (cell, _) = glm::cell_matrices(s)?;
    Some((cell, Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32)))
}

#[derive(Default)]
pub struct Overlays {
    pub show_planes: bool,
    pub planes: Vec<LatticePlane>,
    pub show_directions: bool,
    pub directions: Vec<MillerDirection>,
    pub planes_open: bool,
    pub directions_open: bool,
    new_plane: LatticePlane,
    new_direction: MillerDirection,
    /// View > Show > Voronoi Volumes, and the diagram (kept per structure).
    pub show_voronoi: bool,
    voronoi: Option<(u64, af_algorithms::voronoi::VoronoiDiagram, Vec<f32>)>,
    /// View > Show > Dislocation Lines (the dislocation builder's loop).
    pub show_dislocation_lines: bool,
}

/// The volume of a convex cell from its faces (divergence theorem).
pub fn cell_volume(cell: &af_algorithms::voronoi::VoronoiCell) -> f32 {
    let mut volume = 0.0f32;
    for face in &cell.faces {
        if face.vertices.len() < 3 {
            continue;
        }
        let v0 = face.vertices[0];
        for i in 1..face.vertices.len() - 1 {
            volume += glm::dot(v0, glm::cross(face.vertices[i], face.vertices[i + 1]));
        }
    }
    volume.abs() / 6.0
}

fn structure_key(s: &Structure) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    for a in &s.atoms {
        a.atomic_number.hash(&mut h);
        for v in [a.x, a.y, a.z] {
            v.to_bits().hash(&mut h);
        }
    }
    for row in &s.cell_vectors {
        for v in row {
            v.to_bits().hash(&mut h);
        }
    }
    s.has_unit_cell.hash(&mut h);
    h.finish()
}

impl Overlays {
    pub fn add_plane(&mut self, p: LatticePlane) -> Result<(), String> {
        if p.h == 0 && p.k == 0 && p.l == 0 {
            return Err("h, k and l cannot all be zero.".into());
        }
        self.planes.push(p);
        self.show_planes = true;
        Ok(())
    }

    pub fn add_direction(&mut self, d: MillerDirection) -> Result<(), String> {
        if d.u == 0 && d.v == 0 && d.w == 0 {
            return Err("u, v and w cannot all be zero.".into());
        }
        self.directions.push(d);
        self.show_directions = true;
        Ok(())
    }

    /// The Voronoi cells of a structure and their volumes (computed once per structure).
    pub fn voronoi(&mut self, s: &Structure) -> (&af_algorithms::voronoi::VoronoiDiagram, &[f32]) {
        let key = structure_key(s);
        if self.voronoi.as_ref().is_none_or(|v| v.0 != key) {
            let diagram = af_algorithms::voronoi::compute_voronoi(s);
            let volumes = diagram.cells.iter().map(cell_volume).collect();
            self.voronoi = Some((key, diagram, volumes));
        }
        let v = self.voronoi.as_ref().unwrap();
        (&v.1, &v.2)
    }

    /// Voronoi cells (light blue small to dark blue large; selected atoms red)
    /// and the dislocation loop (drawVoronoiOverlay, drawLineLoop).
    pub fn draw_cells(&mut self, painter: &egui::Painter, rect: egui::Rect, scale: f32, s: &Structure, selected: &[usize], frame: &FrameView) {
        let screen = |p: Vec3| crate::view::project(p, frame).map(|(x, y)| rect.min + egui::vec2(x / scale, y / scale));
        if self.show_voronoi && !s.atoms.is_empty() {
            let (diagram, volumes) = self.voronoi(s);
            let (lo, hi) = volumes.iter().fold((f32::MAX, 0.0f32), |(a, b), &v| (a.min(v), b.max(v)));
            let range = if hi - lo > 1e-10 { hi - lo } else { 1.0 };
            for (ci, cell) in diagram.cells.iter().enumerate() {
                let t = (volumes[ci] - lo) / range;
                let (fill, edge) = if selected.contains(&ci) {
                    (crate::theme::srgba(1.0, 0.2, 0.2, 0.15), crate::theme::srgba(1.0, 0.2, 0.2, 0.85))
                } else {
                    (
                        crate::theme::srgba(0.6 * (1.0 - t), 0.8 * (1.0 - t) + 0.1 * t, (1.0 - t) + 0.6 * t, 0.08 + 0.12 * t),
                        crate::theme::srgba(0.4 * (1.0 - t), 0.6 * (1.0 - t) + 0.1 * t, (1.0 - t) + 0.7 * t, 0.75),
                    )
                };
                for face in &cell.faces {
                    if face.vertices.len() < 3 {
                        continue;
                    }
                    let clipped = clip_near(&face.vertices, frame, 0.03);
                    let points: Option<Vec<egui::Pos2>> = clipped.into_iter().map(screen).collect();
                    let Some(points) = points.filter(|p| p.len() >= 3) else { continue };
                    painter.add(egui::Shape::convex_polygon(points.clone(), fill, egui::Stroke::NONE));
                    painter.add(egui::Shape::closed_line(points, egui::Stroke::new(1.5, edge)));
                }
            }
        }
        if self.show_dislocation_lines && s.dislocation_loop_points.len() >= 2 {
            let points: Vec<egui::Pos2> = s.dislocation_loop_points.iter().filter_map(|p| screen(Vec3::from_array(*p))).collect();
            painter.add(egui::Shape::closed_line(points, egui::Stroke::new(3.5, egui::Color32::from_rgb(255, 64, 26))));
        }
    }

    /// Draws the planes and arrows over the view.
    pub fn draw(&self, painter: &egui::Painter, rect: egui::Rect, scale: f32, s: &Structure, frame: &FrameView) {
        let Some((cell, origin)) = cell_frame(s) else { return };
        let screen = |p: Vec3| crate::view::project(p, frame).map(|(x, y)| rect.min + egui::vec2(x / scale, y / scale));
        if self.show_planes {
            for plane in self.planes.iter().filter(|p| p.visible && (p.h, p.k, p.l) != (0, 0, 0)) {
                let world: Vec<Vec3> = plane_polygon(plane.h, plane.k, plane.l, plane.offset).into_iter().map(|f| origin + glm::mul(&cell, f)).collect();
                let points: Vec<egui::Pos2> = clip_near(&world, frame, 0.03).into_iter().filter_map(screen).collect();
                if points.len() < 3 {
                    continue;
                }
                painter.add(egui::Shape::convex_polygon(points.clone(), colour(plane.colour, plane.opacity.clamp(0.0, 1.0)), egui::Stroke::NONE));
                painter.add(egui::Shape::closed_line(points, egui::Stroke::new(2.0, colour(plane.colour, 0.95))));
            }
        }
        if self.show_directions {
            for d in self.directions.iter().filter(|d| d.visible && (d.u, d.v, d.w) != (0, 0, 0)) {
                let along = glm::mul(&cell, Vec3::new(d.u as f32, d.v as f32, d.w as f32));
                let len = glm::length(along);
                if len < 1e-6 {
                    continue;
                }
                let (mut p0, mut p1) = (origin, origin + along / len * d.length);
                let z = |p: Vec3| (frame.view * Vec4::new(p.x, p.y, p.z, 1.0)).z;
                let (z0, z1) = (z(p0), z(p1));
                let (in0, in1) = (z0 <= -0.03, z1 <= -0.03);
                if !in0 && !in1 {
                    continue;
                }
                if in0 != in1 && (z0 - z1).abs() > 1e-6 {
                    let mid = p0 + (p1 - p0) * ((z0 + 0.03) / (z0 - z1)).clamp(0.0, 1.0);
                    if in0 {
                        p1 = mid;
                    } else {
                        p0 = mid;
                    }
                }
                let (Some(a), Some(b)) = (screen(p0), screen(p1)) else { continue };
                let c = colour(d.colour, 1.0);
                painter.line_segment([a, b], egui::Stroke::new(2.5, c));
                let n = (b - a).normalized();
                let side = egui::vec2(-n.y, n.x);
                painter.add(egui::Shape::convex_polygon(vec![b, b - n * 10.0 + side * 4.5, b - n * 10.0 - side * 4.5], c, egui::Stroke::NONE));
                let label = format!("[{}{}{}]", d.u, d.v, d.w);
                let at = b + egui::vec2(6.0, -4.0);
                painter.text(at + egui::vec2(1.0, 1.0), egui::Align2::LEFT_TOP, &label, egui::FontId::proportional(13.0), egui::Color32::from_black_alpha(180));
                painter.text(at, egui::Align2::LEFT_TOP, label, egui::FontId::proportional(13.0), c);
            }
        }
    }

    /// The Lattice Planes and Miller Directions windows (the C++ popups of
    /// FileBrowser): 780 and 680 px wide, the new item's fields, the list in
    /// a bordered panel and a Close footer.
    pub fn show(&mut self, ctx: &egui::Context, s: &Structure) {
        use crate::dialog_widgets as w;
        let cell = s.has_unit_cell;
        if self.planes_open {
            let mut open = true;
            let mut close = false;
            crate::theme::modal_window("Lattice Planes").id(egui::Id::unique("lattice planes")).open(&mut open).resizable(false).min_width(748.0).max_width(748.0).show(ctx, |ui| {
                ui.set_width(748.0);
                if !cell {
                    w::note(ui, w::Note::Info, "Lattice planes need a unit cell, and this structure has none.");
                    close = w::footer(ui, "Close", None, true) == 1;
                    return;
                }
                w::checkbox(ui, &mut self.show_planes, "Show lattice planes");
                w::disabled(ui, "Planes update in the view as you edit them.");
                w::section(ui, "New Plane");
                let p = &mut self.new_plane;
                ui.horizontal(|ui| {
                    w::drag_int(ui, &mut p.h, 85.0, 0.2, -50..=50, "h");
                    w::drag_int(ui, &mut p.k, 85.0, 0.2, -50..=50, "k");
                    w::drag_int(ui, &mut p.l, 85.0, 0.2, -50..=50, "l");
                });
                w::drag_f32(ui, &mut p.offset, 210.0, 0.01, -50.0..=50.0, 3, "Offset n");
                ui.horizontal(|ui| {
                    w::slider_f32(ui, &mut p.opacity, 0.0..=1.0, 210.0, 2, "Opacity");
                    w::colour_edit(ui, &mut p.colour, "Colour");
                });
                let invalid = (p.h, p.k, p.l) == (0, 0, 0);
                if invalid {
                    w::note(ui, w::Note::Warning, "h, k and l cannot all be zero.");
                }
                if ui.add_enabled(!invalid, egui::Button::new("Add Plane").min_size(egui::vec2(140.0, crate::theme::FRAME_HEIGHT))).clicked() {
                    let new = self.new_plane.clone();
                    let _ = self.add_plane(new);
                }
                w::section(ui, "Planes");
                if self.planes.is_empty() {
                    w::disabled(ui, "No lattice planes added.");
                } else {
                    let mut delete = None;
                    w::child(ui, |ui| {
                        egui::ScrollArea::vertical().max_height(204.0).show(ui, |ui| {
                            for (i, plane) in self.planes.iter_mut().enumerate() {
                                ui.push_id(i, |ui| {
                                    ui.horizontal(|ui| {
                                        w::checkbox(ui, &mut plane.visible, "");
                                        ui.label(format!("Plane {}", i + 1));
                                    });
                                    ui.horizontal(|ui| {
                                        w::drag_int(ui, &mut plane.h, 80.0, 0.2, -50..=50, "h");
                                        w::drag_int(ui, &mut plane.k, 80.0, 0.2, -50..=50, "k");
                                        w::drag_int(ui, &mut plane.l, 80.0, 0.2, -50..=50, "l");
                                        w::drag_f32(ui, &mut plane.offset, 120.0, 0.01, -50.0..=50.0, 3, "n");
                                        w::slider_f32(ui, &mut plane.opacity, 0.0..=1.0, 95.0, 2, "opacity");
                                        if (plane.h, plane.k, plane.l) == (0, 0, 0) {
                                            plane.h = 1;
                                        }
                                        w::colour_edit(ui, &mut plane.colour, "");
                                        if w::danger_button(ui, "Delete").clicked() {
                                            delete = Some(i);
                                        }
                                    });
                                    w::rule(ui);
                                });
                            }
                        });
                    });
                    if let Some(i) = delete {
                        self.planes.remove(i);
                    }
                }
                close = w::footer(ui, "Close", None, true) == 1;
            });
            self.planes_open = open && !close;
        }
        if self.directions_open {
            let mut open = true;
            let mut close = false;
            crate::theme::modal_window("Miller Directions").id(egui::Id::unique("miller directions")).open(&mut open).resizable(false).min_width(648.0).max_width(648.0).show(ctx, |ui| {
                ui.set_width(648.0);
                if !cell {
                    w::note(ui, w::Note::Info, "Miller directions need a unit cell, and this structure has none.");
                    close = w::footer(ui, "Close", None, true) == 1;
                    return;
                }
                w::checkbox(ui, &mut self.show_directions, "Show Miller directions");
                w::disabled(ui, "Arrows are drawn from the cell origin along the [uvw] direction.");
                w::section(ui, "New Direction");
                let d = &mut self.new_direction;
                ui.horizontal(|ui| {
                    w::drag_int(ui, &mut d.u, 85.0, 0.2, -50..=50, "u");
                    w::drag_int(ui, &mut d.v, 85.0, 0.2, -50..=50, "v");
                    w::drag_int(ui, &mut d.w, 85.0, 0.2, -50..=50, "w");
                });
                ui.horizontal(|ui| {
                    w::drag_f32(ui, &mut d.length, 210.0, 0.1, 0.5..=100.0, 1, "Length (Å)");
                    w::colour_edit(ui, &mut d.colour, "Colour");
                });
                let invalid = (d.u, d.v, d.w) == (0, 0, 0);
                if invalid {
                    w::note(ui, w::Note::Warning, "u, v and w cannot all be zero.");
                }
                if ui.add_enabled(!invalid, egui::Button::new("Add Direction").min_size(egui::vec2(140.0, crate::theme::FRAME_HEIGHT))).clicked() {
                    let new = self.new_direction.clone();
                    let _ = self.add_direction(new);
                }
                w::section(ui, "Directions");
                if self.directions.is_empty() {
                    w::disabled(ui, "No Miller directions added.");
                } else {
                    let mut delete = None;
                    w::child(ui, |ui| {
                        egui::ScrollArea::vertical().max_height(184.0).show(ui, |ui| {
                            for (i, d) in self.directions.iter_mut().enumerate() {
                                ui.push_id(i, |ui| {
                                    ui.horizontal(|ui| {
                                        w::checkbox(ui, &mut d.visible, "");
                                        ui.label(format!("Direction {}", i + 1));
                                    });
                                    ui.horizontal(|ui| {
                                        w::drag_int(ui, &mut d.u, 80.0, 0.2, -50..=50, "u");
                                        w::drag_int(ui, &mut d.v, 80.0, 0.2, -50..=50, "v");
                                        w::drag_int(ui, &mut d.w, 80.0, 0.2, -50..=50, "w");
                                        w::drag_f32(ui, &mut d.length, 100.0, 0.1, 0.5..=100.0, 1, "Å");
                                        if (d.u, d.v, d.w) == (0, 0, 0) {
                                            d.u = 1;
                                        }
                                        w::colour_edit(ui, &mut d.colour, "");
                                        if w::danger_button(ui, "Delete").clicked() {
                                            delete = Some(i);
                                        }
                                    });
                                    w::rule(ui);
                                });
                            }
                        });
                    });
                    if let Some(i) = delete {
                        self.directions.remove(i);
                    }
                }
                close = w::footer(ui, "Close", None, true) == 1;
            });
            self.directions_open = open && !close;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planes_cut_the_cell() {
        assert_eq!(plane_polygon(1, 0, 0, 0.5).len(), 4, "a square");
        assert_eq!(plane_polygon(1, 1, 1, 1.0).len(), 3, "a triangle");
        assert_eq!(plane_polygon(1, 1, 1, 1.5).len(), 6, "a hexagon");
        assert!(plane_polygon(1, 0, 0, 2.0).is_empty(), "outside the cell");
        assert_eq!(plane_polygon(1, 0, 0, 1.0).len(), 4, "a face");
        // In order around the centre: consecutive corners share an edge of the square.
        let sq = plane_polygon(0, 0, 1, 0.25);
        for i in 0..4 {
            let d = sq[(i + 1) % 4] - sq[i];
            assert!((d.length() - 1.0).abs() < 1e-5, "{sq:?}");
        }
    }
}
