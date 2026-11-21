//! Voronoi cells by bisector clipping (the C++ algorithms/VoronoiComputation):
//! each atom's cell starts as a box and is clipped by the bisector planes of
//! its neighbours, nearest first; periodic cells of atoms on the cell faces
//! are replicated as the renderer shows their images.

use af_core::cxx::sort_by;
use af_core::glm::{self, cross, dot, length, mul, normalize};
use af_core::neighbours;
use af_core::Structure;
use glam::{DMat3, Vec3};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VoronoiFace {
    /// Ordered convex polygon vertices in world space.
    pub vertices: Vec<Vec3>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VoronoiCell {
    pub faces: Vec<VoronoiFace>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VoronoiDiagram {
    pub cells: Vec<VoronoiCell>,
}

fn max_f(a: f32, b: f32) -> f32 {
    if a < b {
        b
    } else {
        a
    }
}
fn min_f(a: f32, b: f32) -> f32 {
    if b < a {
        b
    } else {
        a
    }
}

/// Orders polygon vertices around their centroid into a convex loop.
fn sort_face_vertices(face: &mut VoronoiFace) {
    if face.vertices.len() < 3 {
        return;
    }
    let mut center = Vec3::ZERO;
    for &v in &face.vertices {
        center += v;
    }
    center /= face.vertices.len() as f32;
    let v = &face.vertices;
    let normal = normalize(cross(v[1] - v[0], v[2] - v[0]));
    let reference = normalize(v[0] - center);
    let binormal = cross(normal, reference);
    sort_by(&mut face.vertices, |a, b| {
        let (da, db) = (*a - center, *b - center);
        dot(da, binormal).atan2(dot(da, reference)) < dot(db, binormal).atan2(dot(db, reference))
    });
}

fn add_unique(points: &mut Vec<Vec3>, p: Vec3) {
    if points.iter().any(|&q| dot(p - q, p - q) <= 1e-8) {
        return;
    }
    points.push(p);
}

/// The part of a convex polygon where dot(v - point, normal) <= 0.
fn clip_polygon(poly: &[Vec3], point: Vec3, normal: Vec3) -> Vec<Vec3> {
    if poly.len() < 3 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(poly.len() + 2);
    for i in 0..poly.len() {
        let curr = poly[i];
        let next = poly[(i + 1) % poly.len()];
        let d_curr = dot(curr - point, normal);
        let d_next = dot(next - point, normal);
        let curr_inside = d_curr <= 1e-6;
        let next_inside = d_next <= 1e-6;
        if curr_inside {
            out.push(curr);
        }
        if curr_inside != next_inside {
            let denom = d_curr - d_next;
            if denom.abs() > 1e-10 {
                let t = max_f(0.0, min_f(1.0, d_curr / denom));
                out.push(curr + t * (next - curr));
            }
        }
    }
    if out.len() >= 3 {
        out
    } else {
        Vec::new()
    }
}

/// A convex polyhedron clipped by a half-space, with the new cap face.
fn clip_polyhedron(poly: &[VoronoiFace], point: Vec3, normal: Vec3) -> Vec<VoronoiFace> {
    let mut result: Vec<VoronoiFace> = poly
        .iter()
        .map(|f| clip_polygon(&f.vertices, point, normal))
        .filter(|c| c.len() >= 3)
        .map(|vertices| VoronoiFace { vertices })
        .collect();
    let mut cut = Vec::new();
    for face in &result {
        for &v in &face.vertices {
            if dot(v - point, normal).abs() < 1e-4 {
                add_unique(&mut cut, v);
            }
        }
    }
    if cut.len() >= 3 {
        let mut cap = VoronoiFace { vertices: cut };
        sort_face_vertices(&mut cap);
        result.push(cap);
    }
    result
}

fn make_box(lo: Vec3, hi: Vec3) -> Vec<VoronoiFace> {
    let v = |x: f32, y: f32, z: f32| Vec3::new(x, y, z);
    let face = |vertices: [Vec3; 4]| VoronoiFace { vertices: vertices.to_vec() };
    vec![
        face([v(lo.x, lo.y, lo.z), v(hi.x, lo.y, lo.z), v(hi.x, lo.y, hi.z), v(lo.x, lo.y, hi.z)]),
        face([v(lo.x, hi.y, lo.z), v(lo.x, hi.y, hi.z), v(hi.x, hi.y, hi.z), v(hi.x, hi.y, lo.z)]),
        face([v(lo.x, lo.y, lo.z), v(lo.x, hi.y, lo.z), v(hi.x, hi.y, lo.z), v(hi.x, lo.y, lo.z)]),
        face([v(lo.x, lo.y, hi.z), v(hi.x, lo.y, hi.z), v(hi.x, hi.y, hi.z), v(lo.x, hi.y, hi.z)]),
        face([v(lo.x, lo.y, lo.z), v(lo.x, lo.y, hi.z), v(lo.x, hi.y, hi.z), v(lo.x, hi.y, lo.z)]),
        face([v(hi.x, lo.y, lo.z), v(hi.x, hi.y, lo.z), v(hi.x, hi.y, hi.z), v(hi.x, lo.y, hi.z)]),
    ]
}

#[derive(Clone, Copy)]
struct Neighbor {
    dist2: f32,
    pos: Vec3,
}

/// Voronoi cells of all atoms (periodic with a unit cell, else within a
/// padded bounding box).
pub fn compute_voronoi(structure: &Structure) -> VoronoiDiagram {
    let mut diagram = VoronoiDiagram::default();
    if structure.atoms.is_empty() {
        return diagram;
    }
    let use_pbc = structure.has_unit_cell;
    let (mut cell, mut inv, mut origin) = (glam::Mat3::IDENTITY, glam::Mat3::IDENTITY, Vec3::ZERO);
    if use_pbc {
        let Some((c, i)) = glm::cell_matrices(structure) else { return diagram };
        cell = c;
        inv = i;
        origin = Vec3::new(structure.cell_offset[0] as f32, structure.cell_offset[1] as f32, structure.cell_offset[2] as f32);
    }
    let n = structure.atoms.len();
    let mut positions: Vec<Vec3> = structure.atoms.iter().map(|a| Vec3::new(a.x as f32, a.y as f32, a.z as f32)).collect();
    // Wrap as the renderer does: into [0, 1), near-0 and near-1 snapped to 0.
    const SNAP: f32 = 1e-4;
    if use_pbc {
        for p in positions.iter_mut() {
            let mut f = mul(&inv, *p - origin);
            for k in 0..3 {
                f[k] -= f[k].floor();
                if f[k].abs() <= SNAP {
                    f[k] = 0.0;
                }
                if (1.0 - f[k]).abs() <= SNAP {
                    f[k] = 0.0;
                }
            }
            *p = origin + mul(&cell, f);
        }
    }
    let mut bbox_lo = Vec3::splat(1e30);
    let mut bbox_hi = Vec3::splat(-1e30);
    if !use_pbc {
        for &p in &positions {
            bbox_lo = glm::min(bbox_lo, p);
            bbox_hi = glm::max(bbox_hi, p);
        }
        let extent = bbox_hi - bbox_lo;
        let pad = max_f(max_f(extent.x, extent.y), extent.z) * 0.5 + 2.0;
        bbox_lo -= Vec3::splat(pad);
        bbox_hi += Vec3::splat(pad);
    }
    // The cell radius is at most the longest cell vector.
    let max_cell_len = if use_pbc {
        // std::max({a, b, c}): the first largest.
        let (a, b, c) = (length(cell.x_axis), length(cell.y_axis), length(cell.z_axis));
        let mut best = a;
        for x in [b, c] {
            if best < x {
                best = x;
            }
        }
        best
    } else {
        0.0
    };

    let full_neighbours = |i: usize| -> Vec<Neighbor> {
        let pi = positions[i];
        let mut list = Vec::new();
        if use_pbc {
            for (j, &pj0) in positions.iter().enumerate() {
                for da in -1..=1 {
                    for db in -1..=1 {
                        for dc in -1..=1 {
                            if j == i && da == 0 && db == 0 && dc == 0 {
                                continue;
                            }
                            let pj = pj0 + mul(&cell, Vec3::new(da as f32, db as f32, dc as f32));
                            let diff = pj - pi;
                            let d2 = dot(diff, diff);
                            if d2 < 1e-6 {
                                continue; // coincident atoms
                            }
                            list.push(Neighbor { dist2: d2, pos: pj });
                        }
                    }
                }
            }
        } else {
            for (j, &pj) in positions.iter().enumerate() {
                if j != i {
                    let diff = pj - pi;
                    list.push(Neighbor { dist2: dot(diff, diff), pos: pj });
                }
            }
        }
        list
    };
    // The cell of atom i clipped by the bisectors of `neighbors`, with the
    // squared distance to its farthest vertex.
    let clip = |i: usize, mut neighbors: Vec<Neighbor>| -> (Vec<VoronoiFace>, f32) {
        let pi = positions[i];
        sort_by(&mut neighbors, |a, b| a.dist2 < b.dist2);
        let mut poly = if use_pbc { make_box(pi - Vec3::splat(max_cell_len), pi + Vec3::splat(max_cell_len)) } else { make_box(bbox_lo, bbox_hi) };
        let farthest = |poly: &[VoronoiFace]| {
            let mut r2 = 0.0f32;
            for face in poly {
                for &v in &face.vertices {
                    r2 = max_f(r2, dot(v - pi, v - pi));
                }
            }
            r2
        };
        let mut rmax2 = farthest(&poly);
        for nb in &neighbors {
            if poly.is_empty() {
                break;
            }
            // The bisector lies at dist/2: beyond rmax no neighbour can clip.
            if nb.dist2 > 4.0 * rmax2 {
                break;
            }
            let len = nb.dist2.sqrt();
            let midpoint = 0.5 * (pi + nb.pos);
            let normal = (nb.pos - pi) / len;
            poly = clip_polyhedron(&poly, midpoint, normal);
            rmax2 = farthest(&poly);
        }
        (poly, rmax2)
    };

    // Candidate neighbours within a few atomic spacings; a cell whose
    // farthest vertex lies within half the cutoff is complete.
    let volume = if use_pbc {
        let c = DMat3::from_cols(cell.x_axis.as_dvec3(), cell.y_axis.as_dvec3(), cell.z_axis.as_dvec3());
        glm::determinant_d(&c).abs()
    } else {
        let mut lo = Vec3::splat(1e30);
        let mut hi = Vec3::splat(-1e30);
        for &p in &positions {
            lo = glm::min(lo, p);
            hi = glm::max(hi, p);
        }
        let e = glm::max(hi - lo, Vec3::ONE);
        e.x as f64 * e.y as f64 * e.z as f64
    };
    let mut search_cutoff = 3.0 * ((if volume < 1e-9 { 1e-9 } else { volume }) / n as f64).cbrt();
    let mut candidates = Vec::new();
    let mut start = vec![0usize; n + 1];
    match neighbours::within(structure, search_cutoff, use_pbc) {
        Ok(pairs) => {
            for p in &pairs {
                start[p.i + 1] += 1;
            }
            for i in 0..n {
                start[i + 1] += start[i];
            }
            candidates = pairs;
        }
        Err(_) => search_cutoff = 0.0, // a cell too small to search: every atom takes the full list
    }
    diagram.cells = vec![VoronoiCell::default(); n];
    for i in 0..n {
        let pi = positions[i];
        let mut result = None;
        if search_cutoff > 0.0 {
            let nearby: Vec<Neighbor> = candidates[start[i]..start[i + 1]]
                .iter()
                .filter_map(|pair| {
                    let d2 = (pair.r * pair.r) as f32;
                    (d2 >= 1e-6).then(|| Neighbor { dist2: d2, pos: pi + Vec3::new(pair.d.x as f32, pair.d.y as f32, pair.d.z as f32) })
                })
                .collect();
            let (poly, rmax2) = clip(i, nearby);
            if !poly.is_empty() && 4.0 * (rmax2 as f64) < search_cutoff * search_cutoff * 0.999 {
                result = Some(poly);
            }
        }
        let mut faces = result.unwrap_or_else(|| clip(i, full_neighbours(i)).0);
        for face in faces.iter_mut() {
            sort_face_vertices(face);
        }
        diagram.cells[i] = VoronoiCell { faces };
    }
    // Copies of the cells of atoms on the cell faces, as the renderer shows
    // their periodic images.
    if use_pbc {
        let tol = if structure.pbc_boundary_tol > 0.0 { structure.pbc_boundary_tol } else { SNAP };
        let safe_tol = min_f(tol, 0.5);
        let mut replicated = Vec::new();
        for i in 0..n {
            let base = &diagram.cells[i];
            if base.faces.is_empty() {
                continue;
            }
            let frac = mul(&inv, positions[i] - origin);
            let options = |x: f32| if x.abs() <= safe_tol { vec![0, 1] } else { vec![0] };
            for a in options(frac.x) {
                for b in options(frac.y) {
                    for c in options(frac.z) {
                        if a == 0 && b == 0 && c == 0 {
                            continue;
                        }
                        let shift = mul(&cell, Vec3::new(a as f32, b as f32, c as f32));
                        replicated.push(VoronoiCell {
                            faces: base.faces.iter().map(|f| VoronoiFace { vertices: f.vertices.iter().map(|&v| v + shift).collect() }).collect(),
                        });
                    }
                }
            }
        }
        diagram.cells.extend(replicated);
    }
    diagram
}

