//! View > Atom Style > Polyhedral and View > Polyhedral Viewer (the C++
//! drawPolyhedralOverlay and Settings > Polyhedral Settings): translucent
//! coordination polyhedra over the view, faces of all polyhedra sorted back
//! to front, ligands outside the cell drawn as spheres; the geometry is
//! kept until the structure or the settings change.

use af_core::structure::Structure;
use af_render::polyhedra::{build_polyhedra, parse_elements, parse_indices, Polyhedron, PolyhedraSettings};
use af_render::scene::FrameView;
use eframe::egui;
use glam::{Vec3, Vec4};

#[derive(Default)]
pub struct PolyhedraView {
    pub settings: PolyhedraSettings,
    /// View > Polyhedral Viewer (the Polyhedral atom style also shows them).
    pub viewer: bool,
    pub settings_open: bool,
    index_text: String,
    centre_text: String,
    ligand_text: String,
    cache: Option<(u64, Vec<Polyhedron>)>,
}

fn key(s: &Structure, settings: &PolyhedraSettings, colours: &[[f32; 3]]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    format!("{settings:?}").hash(&mut h);
    for a in &s.atoms {
        a.atomic_number.hash(&mut h);
        for v in [a.x, a.y, a.z] {
            v.to_bits().hash(&mut h);
        }
    }
    s.has_unit_cell.hash(&mut h);
    for row in &s.cell_vectors {
        for v in row {
            v.to_bits().hash(&mut h);
        }
    }
    for c in colours {
        for v in c {
            v.to_bits().hash(&mut h);
        }
    }
    h.finish()
}

fn colour(c: [f32; 3], alpha: f32) -> egui::Color32 {
    crate::theme::srgba(c[0].clamp(0.0, 1.0), c[1].clamp(0.0, 1.0), c[2].clamp(0.0, 1.0), alpha).into()
}

impl PolyhedraView {
    /// The polyhedra of a structure (rebuilt only when it or the settings change).
    pub fn polyhedra(&mut self, s: &Structure, colours: &[[f32; 3]]) -> &[Polyhedron] {
        let k = key(s, &self.settings, colours);
        if self.cache.as_ref().is_none_or(|c| c.0 != k) {
            self.cache = Some((k, build_polyhedra(s, &self.settings, colours)));
        }
        &self.cache.as_ref().unwrap().1
    }

    pub fn draw(&mut self, painter: &egui::Painter, rect: egui::Rect, scale: f32, s: &Structure, colours: &[[f32; 3]], frame: &FrameView) {
        let settings = self.settings.clone();
        let polyhedra = self.polyhedra(s, colours).to_vec();
        if polyhedra.is_empty() {
            return;
        }
        let near = 0.03;
        let view_z = |p: Vec3| (frame.view * Vec4::new(p.x, p.y, p.z, 1.0)).z;
        let screen = |p: Vec3| crate::view::project(p, frame).map(|(x, y)| rect.min + egui::vec2(x / scale, y / scale));
        // Ligands outside the cell: spheres, so they look like the atoms in it.
        let inverse = af_core::glm::cell_matrices(s).map(|(_, inv)| inv);
        let origin = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
        let mut outside: Vec<(Vec3, [f32; 3])> = Vec::new();
        if let Some(inv) = &inverse {
            for p in &polyhedra {
                for (q, c) in p.ligands.iter().zip(&p.ligand_colours) {
                    let f = af_core::glm::mul(inv, *q - origin);
                    let out = [f.x, f.y, f.z].iter().any(|&v| !(-1e-4..=1.0 + 1e-4).contains(&v));
                    if out && !outside.iter().any(|(o, _)| (*o - *q).length_squared() <= 1e-6) {
                        outside.push((*q, *c));
                    }
                }
            }
        }
        let projection_y = frame.projection.y_axis.y;
        for (q, c) in outside {
            let z = view_z(q);
            if z >= -near {
                continue;
            }
            let radius = 0.65 * projection_y * (frame.height as f32 * 0.5) / -z / scale;
            if radius < 1.0 {
                continue;
            }
            let Some(at) = screen(q) else { continue };
            let base = crate::theme::srgb(c[0], c[1], c[2]);
            painter.circle_filled(at, radius + 1.0, colour([c[0] * 0.45, c[1] * 0.45, c[2] * 0.45], 1.0));
            painter.circle_filled(at, radius, base);
            painter.circle_filled(at - egui::vec2(radius * 0.28, radius * 0.28), radius * 0.38, colour([c[0] * 0.55 + 0.45, c[1] * 0.55 + 0.45, c[2] * 0.55 + 0.45], 0.7));
        }
        // Every face of every polyhedron, back to front.
        let mut faces: Vec<(f32, Vec<egui::Pos2>, egui::Color32, egui::Color32)> = Vec::new();
        for p in &polyhedra {
            let c = p.centre_atomic_number;
            let centre = if c >= 1 && (c as usize) < colours.len() { colours[c as usize] } else { [0.8; 3] };
            let fill = colour(centre, settings.face_opacity.clamp(0.0, 1.0));
            let edge = colour([centre[0] * 0.55, centre[1] * 0.55, centre[2] * 0.55], settings.edge_opacity.clamp(0.0, 1.0));
            for t in &p.triangles {
                let tri = [p.ligands[t[0]], p.ligands[t[1]], p.ligands[t[2]]];
                // Near-plane clipping of the triangle.
                let mut clipped = Vec::new();
                for i in 0..3 {
                    let (a, b) = (tri[i], tri[(i + 1) % 3]);
                    let (za, zb) = (view_z(a), view_z(b));
                    if za <= -near {
                        clipped.push(a);
                    }
                    if (za <= -near) != (zb <= -near) && (zb - za).abs() > 1e-6 {
                        clipped.push(a + (b - a) * ((-near - za) / (zb - za)).clamp(0.0, 1.0));
                    }
                }
                if clipped.len() < 3 {
                    continue;
                }
                let depth = clipped.iter().map(|&v| view_z(v)).sum::<f32>() / clipped.len() as f32;
                let points: Option<Vec<egui::Pos2>> = clipped.into_iter().map(screen).collect();
                if let Some(points) = points {
                    faces.push((depth, points, fill, edge));
                }
            }
        }
        af_core::cxx::sort_by(&mut faces, |a, b| a.0 < b.0);
        for (_, points, fill, edge) in faces {
            painter.add(egui::Shape::convex_polygon(points.clone(), fill, egui::Stroke::NONE));
            if settings.show_edges {
                painter.add(egui::Shape::closed_line(points, egui::Stroke::new(1.2, edge)));
            }
        }
    }

    /// Settings > Polyhedral Settings.
    pub fn show_settings(&mut self, ctx: &egui::Context) {
        if !self.settings_open {
            return;
        }
        use crate::dialog_widgets as w;
        let mut open = true;
        let mut close = false;
        crate::theme::modal_window("Polyhedral Settings").id(egui::Id::unique("polyhedral settings")).open(&mut open).resizable(false).min_width(588.0).max_width(588.0).show(ctx, |ui| {
            ui.set_width(588.0);
            let s = &mut self.settings;
            w::section(ui, "Polyhedra");
            w::slider_usize(ui, &mut s.max_centres, 1..=2000, 240.0, "Max displayed centers");
            w::slider_usize(ui, &mut s.max_neighbours, 4..=64, 240.0, "Nearest neighbors");
            if w::checkbox(ui, &mut s.centre_index_filter, "Specify center atom indices").changed() && s.centre_index_filter {
                s.centre_indices = parse_indices(&self.index_text);
            }
            if s.centre_index_filter {
                if w::input_text(ui, &mut self.index_text, 420.0, "", "Center atom IDs").changed() {
                    s.centre_indices = parse_indices(&self.index_text);
                }
                w::disabled(ui, "1-based IDs, comma-separated, ranges allowed (e.g. 1,4,10-20)");
                w::disabled(ui, format!("Parsed centers: {}", s.centre_indices.len()));
            }
            w::section(ui, "Appearance");
            w::checkbox(ui, &mut s.show_edges, "Show polyhedral edges");
            w::slider_f32(ui, &mut s.face_opacity, 0.0..=1.0, 220.0, 2, "Face opacity");
            w::slider_f32(ui, &mut s.edge_opacity, 0.0..=1.0, 220.0, 2, "Edge opacity");
            w::section(ui, "Elements");
            if w::checkbox(ui, &mut s.centre_element_filter, "Filter center elements").changed() && s.centre_element_filter {
                s.centre_elements = parse_elements(&self.centre_text);
            }
            if s.centre_element_filter {
                if w::input_text(ui, &mut self.centre_text, 340.0, "", "Centers").changed() {
                    s.centre_elements = parse_elements(&self.centre_text);
                }
                w::disabled(ui, "Comma-separated symbols, e.g. Ti,Zr,Fe");
            }
            ui.add_space(4.0);
            if w::checkbox(ui, &mut s.ligand_element_filter, "Filter ligand elements").changed() && s.ligand_element_filter {
                s.ligand_elements = parse_elements(&self.ligand_text);
            }
            if s.ligand_element_filter {
                if w::input_text(ui, &mut self.ligand_text, 340.0, "", "Ligands").changed() {
                    s.ligand_elements = parse_elements(&self.ligand_text);
                }
                w::disabled(ui, "Comma-separated symbols, e.g. O,N,S,F");
            }
            match w::footer(ui, "Close", Some("Reset Defaults"), true) {
                1 => close = true,
                2 => {
                    self.settings = PolyhedraSettings::default();
                    self.index_text.clear();
                    self.centre_text.clear();
                    self.ligand_text.clear();
                }
                _ => {}
            }
        });
        self.settings_open = open && !close;
    }
}
