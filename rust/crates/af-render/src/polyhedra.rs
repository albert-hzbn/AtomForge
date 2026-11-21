//! Coordination polyhedra (the C++ ui/PolyhedralOverlay geometry): centre
//! atoms (anion-like elements left out unless filtered), their first shell
//! of the dominant ligand element (periodic images included) and the
//! convex hull's triangles.

use std::collections::{BTreeMap, BTreeSet};

use af_core::glm;
use af_core::Structure;
use glam::{Mat3, Vec3};

pub const MIN_NEIGHBOURS: usize = 4;
const FACET_EPSILON: f32 = 1e-4;

#[derive(Clone, Debug, PartialEq)]
pub struct PolyhedraSettings {
    pub centre_index_filter: bool,
    /// 0-based atom indices.
    pub centre_indices: Vec<usize>,
    pub max_centres: usize,
    pub max_neighbours: usize,
    pub show_edges: bool,
    pub face_opacity: f32,
    pub edge_opacity: f32,
    pub centre_element_filter: bool,
    pub ligand_element_filter: bool,
    pub centre_elements: BTreeSet<i32>,
    pub ligand_elements: BTreeSet<i32>,
}

impl Default for PolyhedraSettings {
    fn default() -> Self {
        Self {
            centre_index_filter: false,
            centre_indices: Vec::new(),
            max_centres: 350,
            max_neighbours: 24,
            show_edges: true,
            face_opacity: 0.08,
            edge_opacity: 0.82,
            centre_element_filter: false,
            ligand_element_filter: false,
            centre_elements: BTreeSet::new(),
            ligand_elements: BTreeSet::new(),
        }
    }
}

/// Elements that are almost always ligands or anions, never centres.
const DEFAULT_EXCLUDE: [i32; 17] = [1, 2, 7, 8, 9, 10, 16, 17, 18, 34, 35, 36, 52, 53, 54, 85, 86];

/// The centre atoms (gatherCenters).
pub fn centres(s: &Structure, settings: &PolyhedraSettings) -> Vec<usize> {
    let allowed = |i: usize| {
        let Some(a) = s.atoms.get(i) else { return false };
        let z = a.atomic_number;
        if !(1..=118).contains(&z) {
            return false;
        }
        if settings.centre_element_filter {
            settings.centre_elements.contains(&z)
        } else {
            !DEFAULT_EXCLUDE.contains(&z)
        }
    };
    if settings.centre_index_filter && !settings.centre_indices.is_empty() {
        return settings.centre_indices.iter().copied().filter(|&i| allowed(i)).collect();
    }
    (0..s.atoms.len()).filter(|&i| allowed(i)).take(settings.max_centres.max(1)).collect()
}

/// One polyhedron: the ligand positions and the hull's outward triangles.
#[derive(Clone, Debug, PartialEq)]
pub struct Polyhedron {
    pub centre: Vec3,
    pub centre_atomic_number: i32,
    pub ligands: Vec<Vec3>,
    pub ligand_colours: Vec<[f32; 3]>,
    pub triangles: Vec<[usize; 3]>,
}

fn colour_of(s: &Structure, i: usize, colours: &[[f32; 3]]) -> [f32; 3] {
    let a = &s.atoms[i];
    let z = a.atomic_number;
    if z >= 1 && (z as usize) < colours.len() {
        colours[z as usize]
    } else {
        [a.r, a.g, a.b]
    }
}

/// The first shell of the nearest ligand element around `centre`
/// (buildNeighborCloud): within 1.2 times the nearest distance.
pub fn ligand_shell(s: &Structure, centre: Vec3, centre_z: i32, cell: Option<&Mat3>, settings: &PolyhedraSettings, colours: &[[f32; 3]]) -> Vec<(Vec3, [f32; 3])> {
    let keep = settings.max_neighbours.max(MIN_NEIGHBOURS);
    let mut by_element: BTreeMap<i32, Vec<(f32, Vec3, [f32; 3])>> = BTreeMap::new();
    let shifts: Vec<Vec3> = match cell {
        Some(c) => {
            let mut v = Vec::new();
            for sx in -1..=1 {
                for sy in -1..=1 {
                    for sz in -1..=1 {
                        v.push(c.x_axis * sx as f32 + c.y_axis * sy as f32 + c.z_axis * sz as f32);
                    }
                }
            }
            v
        }
        None => vec![Vec3::ZERO],
    };
    for (j, a) in s.atoms.iter().enumerate() {
        if a.atomic_number == centre_z {
            continue;
        }
        if settings.ligand_element_filter && !settings.ligand_elements.contains(&a.atomic_number) {
            continue;
        }
        let p = Vec3::new(a.x as f32, a.y as f32, a.z as f32);
        for &shift in &shifts {
            let q = p + shift;
            let d = glm::length(q - centre);
            if d <= 1e-5 {
                continue;
            }
            by_element.entry(a.atomic_number).or_default().push((d, q, colour_of(s, j, colours)));
        }
    }
    let mut candidates: Vec<(f32, Vec3, [f32; 3], i32)> = Vec::new();
    for (z, mut list) in by_element {
        list.sort_by(|a, b| a.0.total_cmp(&b.0));
        list.truncate(keep);
        candidates.extend(list.into_iter().map(|(d, q, c)| (d, q, c, z)));
    }
    if candidates.len() < MIN_NEIGHBOURS {
        return Vec::new();
    }
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (nearest, ligand) = (candidates[0].0, candidates[0].3);
    let cutoff = nearest * 1.2;
    let mut shell: Vec<(Vec3, [f32; 3])> = Vec::new();
    for (d, q, c, z) in candidates {
        if z != ligand {
            continue;
        }
        if d > cutoff || shell.len() >= keep {
            break;
        }
        if !shell.iter().any(|(p, _)| glm::dot(*p - q, *p - q) <= 1e-8) {
            shell.push((q, c));
        }
    }
    shell
}

/// The outward triangles of the convex hull of `points` (buildConvexHullTriangles).
pub fn hull_triangles(centre: Vec3, points: &[Vec3]) -> Vec<[usize; 3]> {
    let n = points.len();
    let mut triangles = Vec::new();
    if n < 4 {
        return triangles;
    }
    let mut seen = BTreeSet::new();
    for a in 0..n - 2 {
        for b in a + 1..n - 1 {
            for c in b + 1..n {
                let (pa, pb, pc) = (points[a] - centre, points[b] - centre, points[c] - centre);
                let mut normal = glm::cross(pb - pa, pc - pa);
                let area = glm::length(normal);
                if area < 1e-6 {
                    continue;
                }
                normal /= area;
                let (mut pos, mut neg) = (false, false);
                for (d, p) in points.iter().enumerate() {
                    if d == a || d == b || d == c {
                        continue;
                    }
                    let side = glm::dot(normal, *p - points[a]);
                    if side > FACET_EPSILON {
                        pos = true;
                    } else if side < -FACET_EPSILON {
                        neg = true;
                    }
                    if pos && neg {
                        break;
                    }
                }
                if pos && neg {
                    continue;
                }
                let (ib, ic) = if glm::dot(normal, (pa + pb + pc) * (1.0 / 3.0)) < 0.0 { (c, b) } else { (b, c) };
                if seen.insert((a, b, c)) {
                    triangles.push([a, ib, ic]);
                }
            }
        }
    }
    triangles
}

/// Every polyhedron of a structure, centres' boundary images included.
pub fn build_polyhedra(s: &Structure, settings: &PolyhedraSettings, colours: &[[f32; 3]]) -> Vec<Polyhedron> {
    if s.atoms.len() < 4 {
        return Vec::new();
    }
    let matrices = glm::cell_matrices(s);
    let origin = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
    let mut entries: Vec<(Vec3, i32)> = Vec::new();
    for i in centres(s, settings) {
        let a = &s.atoms[i];
        let p = Vec3::new(a.x as f32, a.y as f32, a.z as f32);
        entries.push((p, a.atomic_number));
        if let Some((cell, inv)) = &matrices {
            let f = glm::mul(inv, p - origin);
            let tol = 0.002;
            for sx in -1..=1 {
                for sy in -1..=1 {
                    for sz in -1..=1 {
                        if (sx, sy, sz) == (0, 0, 0) {
                            continue;
                        }
                        let g = f + Vec3::new(sx as f32, sy as f32, sz as f32);
                        if [g.x, g.y, g.z].iter().all(|&v| (-tol..=1.0 + tol).contains(&v)) {
                            entries.push((origin + glm::mul(cell, g), a.atomic_number));
                        }
                    }
                }
            }
        }
    }
    let cell = matrices.map(|(c, _)| c);
    entries
        .into_iter()
        .filter_map(|(centre, z)| {
            let shell = ligand_shell(s, centre, z, cell.as_ref(), settings, colours);
            if shell.len() < MIN_NEIGHBOURS {
                return None;
            }
            let ligands: Vec<Vec3> = shell.iter().map(|x| x.0).collect();
            let triangles = hull_triangles(centre, &ligands);
            (!triangles.is_empty()).then(|| Polyhedron { centre, centre_atomic_number: z, ligand_colours: shell.iter().map(|x| x.1).collect(), ligands, triangles })
        })
        .collect()
}

/// "1,4,10-20" (1-based) as 0-based indices.
pub fn parse_indices(text: &str) -> Vec<usize> {
    let mut out = Vec::new();
    for part in text.split(',').map(str::trim).filter(|p| !p.is_empty()) {
        match part.split_once('-') {
            Some((a, b)) => {
                if let (Ok(a), Ok(b)) = (a.trim().parse::<usize>(), b.trim().parse::<usize>()) {
                    for i in a.min(b)..=a.max(b) {
                        if i >= 1 {
                            out.push(i - 1);
                        }
                    }
                }
            }
            None => {
                if let Ok(i) = part.parse::<usize>() {
                    if i >= 1 {
                        out.push(i - 1);
                    }
                }
            }
        }
    }
    out
}

/// "Ti, Zr,Fe" as atomic numbers.
pub fn parse_elements(text: &str) -> BTreeSet<i32> {
    text.split([',', ' ']).map(str::trim).filter(|p| !p.is_empty()).map(af_core::elements::atomic_number).filter(|&z| z > 0).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::structure::AtomSite;

    fn atom(symbol: &str, p: [f64; 3]) -> AtomSite {
        AtomSite { symbol: symbol.into(), atomic_number: af_core::elements::atomic_number(symbol), x: p[0], y: p[1], z: p[2], ..Default::default() }
    }

    #[test]
    fn perovskite_octahedron() {
        // SrTiO3: Ti at the body centre, O on the face centres, Sr at the corner.
        let a = 3.905;
        let mut s = Structure { has_unit_cell: true, cell_vectors: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], ..Default::default() };
        s.atoms.push(atom("Sr", [0.0, 0.0, 0.0]));
        s.atoms.push(atom("Ti", [a / 2.0; 3]));
        for f in [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]] {
            s.atoms.push(atom("O", [f[0] * a, f[1] * a, f[2] * a]));
        }
        let settings = PolyhedraSettings { centre_element_filter: true, centre_elements: parse_elements("Ti"), ..Default::default() };
        let p = build_polyhedra(&s, &settings, &af_core::elements::default_colours());
        assert_eq!(p.len(), 1, "one Ti inside the cell");
        assert_eq!(p[0].ligands.len(), 6, "six O");
        assert_eq!(p[0].triangles.len(), 8, "an octahedron");
        // Outward triangles: normals point away from the centre.
        for t in &p[0].triangles {
            let (a, b, c) = (p[0].ligands[t[0]], p[0].ligands[t[1]], p[0].ligands[t[2]]);
            let n = (b - a).cross(c - a);
            assert!(n.dot((a + b + c) / 3.0 - p[0].centre) > 0.0);
        }
        // By default O is never a centre; Sr and Ti are.
        assert_eq!(centres(&s, &PolyhedraSettings::default()), vec![0, 1]);
    }

    #[test]
    fn index_and_element_lists() {
        assert_eq!(parse_indices("1,4,10-12, 3"), vec![0, 3, 9, 10, 11, 2]);
        assert_eq!(parse_elements("Ti, Zr,Fe"), [22, 26, 40].into_iter().collect());
    }
}
