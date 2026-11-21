//! Isosurfaces by marching tetrahedra (six per grid cube) and mesh files
//! (the C++ electronic/Isosurface).

use std::path::Path;

use af_core::cxx::sort_by;
use af_core::format::format_g;
use af_core::glm::normalize_d;
use glam::DVec3;

use crate::volume::{Grid, Mesh};

fn mix(a: DVec3, b: DVec3, t: f64) -> DVec3 {
    a * (1.0 - t) + b * t
}

fn length(v: DVec3) -> f64 {
    v.dot(v).sqrt()
}

/// The surface where the field crosses `level`; vertex scalars from `color`
/// (sampled) or the level.
pub fn isosurface(grid: &Grid, level: f64, color: Option<&Grid>) -> Result<Mesh, String> {
    grid.validate()?;
    if !level.is_finite() {
        return Err("Isovalue must be finite".into());
    }
    if let Some(c) = color {
        c.validate()?;
    }
    let mut out = Mesh::default();
    let offsets = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
    let tetrahedra = [[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]];
    let edges = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
    let open = if grid.periodic { 0 } else { 1 };
    for z in 0..grid.shape[2] - open {
        for y in 0..grid.shape[1] - open {
            for x in 0..grid.shape[0] - open {
                let mut p = [DVec3::ZERO; 8];
                let mut values = [0.0; 8];
                for i in 0..8 {
                    let q = [x + offsets[i][0], y + offsets[i][1], z + offsets[i][2]];
                    p[i] = grid.position(q[0], q[1], q[2]);
                    values[i] = grid.values[grid.index(q[0] % grid.shape[0], q[1] % grid.shape[1], q[2] % grid.shape[2])];
                }
                for tetra in &tetrahedra {
                    let mut hits: Vec<DVec3> = Vec::new();
                    let (mut inside, mut outside) = (DVec3::ZERO, DVec3::ZERO);
                    let (mut ni, mut no) = (0, 0);
                    for &i in tetra {
                        if values[i] < level {
                            inside += p[i];
                            ni += 1;
                        } else {
                            outside += p[i];
                            no += 1;
                        }
                    }
                    if ni == 0 || no == 0 {
                        continue;
                    }
                    for edge in &edges {
                        let (a, b) = (tetra[edge[0]], tetra[edge[1]]);
                        if (values[a] < level) != (values[b] < level) {
                            let hit = mix(p[a], p[b], (level - values[a]) / (values[b] - values[a]));
                            if !hits.iter().any(|v| length(hit - *v) < 1e-12) {
                                hits.push(hit);
                            }
                        }
                    }
                    if hits.len() < 3 {
                        continue;
                    }
                    let normal = normalize_d(outside / no as f64 - inside / ni as f64);
                    let mut center = DVec3::ZERO;
                    for v in &hits {
                        center += *v;
                    }
                    center /= hits.len() as f64;
                    let u = normalize_d(hits[0] - center);
                    let v = normalize_d(normal.cross(u));
                    let angle = |a: &DVec3| (*a - center).dot(v).atan2((*a - center).dot(u));
                    sort_by(&mut hits, |a, b| angle(a) < angle(b));
                    for i in 1..hits.len() - 1 {
                        if length((hits[i] - hits[0]).cross(hits[i + 1] - hits[0])) < 1e-16 {
                            continue;
                        }
                        for vertex in [hits[0], hits[i], hits[i + 1]] {
                            out.vertices.push(vertex);
                            out.colors.push(match color {
                                Some(c) => c.sample(vertex)?,
                                None => level,
                            });
                        }
                    }
                    if out.vertices.len() > 12_000_000 {
                        return Err("Isosurface exceeds four million triangles; reduce resolution".into());
                    }
                }
            }
        }
    }
    Ok(out)
}

fn g17(x: f64) -> String {
    format_g(x, 17)
}

/// Writes a mesh as OBJ (scalars in vertex comments) or, for .ply, PLY.
pub fn save_mesh(mesh: &Mesh, path: &Path) -> Result<(), String> {
    if mesh.vertices.len() % 3 != 0 || mesh.colors.len() != mesh.vertices.len() {
        return Err("Mesh requires triangles and one scalar per vertex".into());
    }
    let mut text = String::new();
    if path.extension().is_some_and(|e| e == "ply") {
        text += &format!(
            "ply\nformat ascii 1.0\nelement vertex {}\nproperty double x\nproperty double y\nproperty double z\nproperty double scalar\nelement face {}\nproperty list uchar int vertex_indices\nend_header\n",
            mesh.vertices.len(),
            mesh.vertices.len() / 3
        );
        for (p, c) in mesh.vertices.iter().zip(&mesh.colors) {
            text += &format!("{} {} {} {}\n", g17(p.x), g17(p.y), g17(p.z), g17(*c));
        }
        for i in (0..mesh.vertices.len()).step_by(3) {
            text += &format!("3 {} {} {}\n", i, i + 1, i + 2);
        }
    } else {
        text += "# AtomForge isosurface; scalar colors in vertex comments\n";
        for (p, c) in mesh.vertices.iter().zip(&mesh.colors) {
            text += &format!("v {} {} {} # {}\n", g17(p.x), g17(p.y), g17(p.z), g17(*c));
        }
        let mut i = 0;
        while i + 2 < mesh.vertices.len() {
            text += &format!("f {} {} {}\n", i + 1, i + 2, i + 3);
            i += 3;
        }
    }
    std::fs::write(path, text).map_err(|_| "Cannot write OBJ mesh".to_string())
}
