//! Voronoi tessellation, VoroTop and adaptive common-neighbour analysis (the
//! C++ science/VoronoiAnalysis).

use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use af_core::cxx::sort_by;
use af_core::cxx_random::{Mt19937, Normal};
use af_core::json::Json;

use crate::analysis::{ji, jv, R};
use crate::core::*;
use crate::data::{read_text_file, Parameters};
use crate::local_environment::{add_classification, read_atom_set, AtomSet};

fn fmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

/// The cube root of the volume per atom (cell, or padded bounding box).
pub(crate) fn mean_atomic_spacing(a: &AtomSet) -> Result<f64, String> {
    if a.positions.is_empty() {
        return Err("No atoms".into());
    }
    let mut volume = 1.0;
    if a.has_cell {
        volume = cell_volume(&a.cell);
    } else {
        let (mut low, mut high) = (a.positions[0], a.positions[0]);
        for x in &a.positions {
            for d in 0..3 {
                if x[d] < low[d] {
                    low[d] = x[d];
                }
                if high[d] < x[d] {
                    high[d] = x[d];
                }
            }
        }
        for d in 0..3 {
            volume *= high[d] - low[d] + 2.0;
        }
    }
    Ok((volume / a.positions.len() as f64).cbrt())
}

/// Each atom's neighbour vectors within the cutoff, nearest first.
pub(crate) fn sorted_neighbor_vectors(a: &AtomSet, cutoff: f64) -> Result<Vec<Vec<Vec3>>, String> {
    let mut lists = vec![Vec::new(); a.positions.len()];
    for link in neighbor_list(&a.positions, &a.cell, a.pbc, cutoff)? {
        if norm(link.vector) < 1e-10 {
            return Err("Coincident atoms make local environment analysis undefined".into());
        }
        lists[link.i].push(link.vector);
    }
    for list in lists.iter_mut() {
        sort_by(list, |x, y| dot(*x, *x) < dot(*y, *y));
    }
    Ok(lists)
}

// ---------------------------------------------------------------- Voronoi cells
#[derive(Clone)]
struct Face {
    neighbor: i32,
    points: Vec<Vec3>,
}

fn polygon_area(pts: &[Vec3]) -> f64 {
    let mut sum = [0.0; 3];
    for k in 1..pts.len().saturating_sub(1) {
        sum = add(sum, cross(sub(pts[k], pts[0]), sub(pts[k + 1], pts[0])));
    }
    0.5 * norm(sum)
}

/// Cuts the convex cell with the half-space n.x <= d; false when the plane misses it.
fn clip(cell: &mut Vec<Face>, n: Vec3, d: f64, neighbor: i32, eps: f64) -> bool {
    let cuts = cell.iter().any(|f| f.points.iter().any(|p| dot(n, *p) - d > eps));
    if !cuts {
        return false;
    }
    let mut out = Vec::new();
    let mut cap: Vec<Vec3> = Vec::new();
    for face in cell.iter() {
        let mut clipped = Face { neighbor: face.neighbor, points: Vec::new() };
        let m = face.points.len();
        for k in 0..m {
            let current = face.points[k];
            let next = face.points[(k + 1) % m];
            let dc = dot(n, current) - d;
            let dn = dot(n, next) - d;
            let inside = dc <= eps;
            let next_inside = dn <= eps;
            if inside {
                clipped.points.push(current);
                if dc.abs() <= eps {
                    cap.push(current);
                }
            }
            if inside != next_inside && (dc - dn).abs() > 0.0 {
                let crossing = add(current, scale(sub(next, current), dc / (dc - dn)));
                clipped.points.push(crossing);
                cap.push(crossing);
            }
        }
        let mut unique: Vec<Vec3> = Vec::new();
        for point in &clipped.points {
            if unique.is_empty() || norm(sub(*point, *unique.last().unwrap())) > eps {
                unique.push(*point);
            }
        }
        while unique.len() > 1 && norm(sub(unique[0], *unique.last().unwrap())) <= eps {
            unique.pop();
        }
        if unique.len() >= 3 {
            clipped.points = unique;
            out.push(clipped);
        }
    }
    let mut cap_unique: Vec<Vec3> = Vec::new();
    for point in &cap {
        if !cap_unique.iter().any(|q| norm(sub(*point, *q)) <= eps) {
            cap_unique.push(*point);
        }
    }
    if cap_unique.len() >= 3 {
        let mut centre = [0.0; 3];
        for point in &cap_unique {
            centre = add(centre, scale(*point, 1.0 / cap_unique.len() as f64));
        }
        let helper = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let u = scale(cross(n, helper), 1.0 / norm(cross(n, helper)));
        let w = cross(n, u);
        let angle = |x: &Vec3| {
            let dx = sub(*x, centre);
            dot(dx, w).atan2(dot(dx, u))
        };
        sort_by(&mut cap_unique, |x, y| angle(x) < angle(y));
        out.push(Face { neighbor, points: cap_unique });
    }
    *cell = out;
    true
}

fn cube(h: f64) -> Vec<Face> {
    let v = |x: f64, y: f64, z: f64| [x * h, y * h, z * h];
    let f = |points: Vec<Vec3>| Face { neighbor: -1, points };
    vec![
        f(vec![v(-1.0, -1.0, -1.0), v(-1.0, 1.0, -1.0), v(1.0, 1.0, -1.0), v(1.0, -1.0, -1.0)]),
        f(vec![v(-1.0, -1.0, 1.0), v(1.0, -1.0, 1.0), v(1.0, 1.0, 1.0), v(-1.0, 1.0, 1.0)]),
        f(vec![v(-1.0, -1.0, -1.0), v(1.0, -1.0, -1.0), v(1.0, -1.0, 1.0), v(-1.0, -1.0, 1.0)]),
        f(vec![v(-1.0, 1.0, -1.0), v(-1.0, 1.0, 1.0), v(1.0, 1.0, 1.0), v(1.0, 1.0, -1.0)]),
        f(vec![v(-1.0, -1.0, -1.0), v(-1.0, -1.0, 1.0), v(-1.0, 1.0, 1.0), v(-1.0, 1.0, -1.0)]),
        f(vec![v(1.0, -1.0, -1.0), v(1.0, 1.0, -1.0), v(1.0, 1.0, 1.0), v(1.0, -1.0, 1.0)]),
    ]
}

#[derive(Clone, Default)]
struct CellResult {
    faces: Vec<Face>,
    radius: f64,
    open: bool,
}

fn voronoi_cell(neighbors: &[Vec3], half: f64, eps: f64) -> CellResult {
    let mut faces = cube(half);
    let mut radius = 3.0f64.sqrt() * half;
    for (k, nb) in neighbors.iter().enumerate() {
        let r = norm(*nb);
        if r > 2.0 * radius {
            break;
        }
        clip(&mut faces, scale(*nb, 1.0 / r), r / 2.0, k as i32, eps);
        if k % 8 == 7 || k + 1 == neighbors.len() {
            radius = 0.0;
            for face in &faces {
                for point in &face.points {
                    radius = fmax(radius, norm(*point));
                }
            }
        }
    }
    let mut result = CellResult { faces, radius: 0.0, open: false };
    for face in &result.faces {
        if face.neighbor < 0 {
            result.open = true;
        }
        for point in &face.points {
            result.radius = fmax(result.radius, norm(*point));
        }
    }
    result
}

fn all_cells(a: &AtomSet, cutoff: &mut f64) -> Result<Vec<CellResult>, String> {
    let n = a.positions.len();
    let periodic = a.pbc[0] && a.pbc[1] && a.pbc[2];
    if *cutoff <= 0.0 {
        *cutoff = 2.5 * mean_atomic_spacing(a)?;
    }
    let eps = 1e-9 * mean_atomic_spacing(a)?;
    let mut cells = vec![CellResult::default(); n];
    let mut done = vec![false; n];
    for attempt in 0..6 {
        let neighbors = sorted_neighbor_vectors(a, *cutoff)?;
        let mut pending = 0;
        for i in 0..n {
            if done[i] {
                continue;
            }
            cells[i] = voronoi_cell(&neighbors[i], *cutoff, eps);
            let complete = !cells[i].open && 2.0 * cells[i].radius <= *cutoff;
            if complete || (cells[i].open && !periodic && attempt > 0) {
                done[i] = true;
            } else {
                pending += 1;
            }
        }
        if pending == 0 {
            break;
        }
        *cutoff *= 1.6;
    }
    if !a.pbc.iter().any(|&b| b) {
        let spacing = mean_atomic_spacing(a)?;
        for i in 0..n {
            if !cells[i].open {
                continue;
            }
            let mut others: Vec<Vec3> = (0..n).filter(|&j| j != i).map(|j| sub(a.positions[j], a.positions[i])).collect();
            sort_by(&mut others, |x, y| dot(*x, *x) < dot(*y, *y));
            cells[i] = voronoi_cell(&others, 1e6 * spacing, eps);
        }
    }
    Ok(cells)
}

// ---------------------------------------------------------------- topology
fn topology_of(cell: &CellResult, eps: f64, vertex_count: &mut i32) -> Vec<Vec<i32>> {
    let mut vertices: Vec<Vec3> = Vec::new();
    let mut faces = Vec::new();
    for face in &cell.faces {
        let mut ids: Vec<i32> = Vec::new();
        for point in &face.points {
            let id = match vertices.iter().position(|v| norm(sub(*v, *point)) <= eps) {
                Some(k) => k as i32,
                None => {
                    vertices.push(*point);
                    vertices.len() as i32 - 1
                }
            };
            if ids.last() != Some(&id) {
                ids.push(id);
            }
        }
        while ids.len() > 1 && ids[0] == *ids.last().unwrap() {
            ids.pop();
        }
        if ids.len() >= 3 {
            faces.push(ids);
        }
    }
    *vertex_count = vertices.len() as i32;
    faces
}

/// The canonical Weinberg code of a polyhedral graph (empty when not closed).
fn weinberg_code(faces: &[Vec<i32>], n: i32) -> Vec<i32> {
    if n < 4 || faces.len() < 4 {
        return Vec::new();
    }
    let nu = n as usize;
    let mut best: Vec<i32> = Vec::new();
    for orientation in 0..2 {
        let mut turn = vec![-1i32; nu * nu];
        let mut edges: Vec<(i32, i32)> = Vec::new();
        for original in faces {
            let mut face = original.clone();
            if orientation == 1 {
                face.reverse();
            }
            let m = face.len();
            for k in 0..m {
                let (u, v, w) = (face[k], face[(k + 1) % m], face[(k + 2) % m]);
                let slot = &mut turn[u as usize * nu + v as usize];
                if *slot >= 0 {
                    return Vec::new();
                }
                *slot = w;
                edges.push((u, v));
            }
        }
        if edges.iter().any(|&(u, v)| turn[v as usize * nu + u as usize] < 0) {
            return Vec::new();
        }
        let total = edges.len();
        let mut used = vec![false; nu * nu];
        let mut label = vec![0i32; nu];
        let mut code: Vec<i32> = Vec::new();
        for &(a, b) in &edges {
            used.iter_mut().for_each(|x| *x = false);
            label.iter_mut().for_each(|x| *x = 0);
            code.clear();
            let mut next = 1;
            label[a as usize] = next;
            code.push(next);
            let (mut u, mut v) = (a, b);
            let mut count = 0usize;
            let mut worse = false;
            loop {
                used[u as usize * nu + v as usize] = true;
                count += 1;
                let w;
                if label[v as usize] == 0 {
                    next += 1;
                    label[v as usize] = next;
                    code.push(next);
                    w = turn[u as usize * nu + v as usize];
                } else {
                    code.push(label[v as usize]);
                    if count == total {
                        break;
                    }
                    if !used[v as usize * nu + u as usize] {
                        w = u;
                    } else {
                        let mut ww = turn[u as usize * nu + v as usize];
                        let mut guard = 0;
                        while used[v as usize * nu + ww as usize] && {
                            let ok = guard < n;
                            guard += 1;
                            ok
                        } {
                            ww = turn[ww as usize * nu + v as usize];
                        }
                        if used[v as usize * nu + ww as usize] {
                            break;
                        }
                        w = ww;
                    }
                }
                if !best.is_empty() {
                    let k = code.len() - 1;
                    if k < best.len() && code[k] != best[k] && code[k] > best[k] && code[..k] == best[..k] {
                        worse = true;
                        break;
                    }
                }
                u = v;
                v = w;
            }
            if !worse && count == total && (best.is_empty() || code < best) {
                best = code.clone();
            }
        }
    }
    best
}

fn code_key(code: &[i32]) -> String {
    code.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(",")
}

fn cell_codes(a: &AtomSet) -> Result<Vec<Vec<i32>>, String> {
    let mut cutoff = 0.0;
    let cells = all_cells(a, &mut cutoff)?;
    let eps = 1e-6 * mean_atomic_spacing(a)?;
    Ok(cells
        .iter()
        .map(|cell| {
            if cell.open {
                return Vec::new();
            }
            let mut n = 0;
            let faces = topology_of(cell, eps, &mut n);
            weinberg_code(&faces, n)
        })
        .collect())
}

#[derive(Clone)]
struct Filter {
    names: Vec<String>,
    types: BTreeMap<String, i32>,
}

/// Topologies of thermally perturbed fcc, bcc and hcp crystals.
fn built_in_filter(noise: f64, samples: i32) -> Result<Filter, String> {
    static CACHE: OnceLock<Mutex<HashMap<(i64, i32), Filter>>> = OnceLock::new();
    let key = ((noise * 1e6).round() as i64, samples);
    if let Some(found) = CACHE.get_or_init(Default::default).lock().unwrap().get(&key) {
        return Ok(found.clone());
    }
    let mut filter = Filter { names: ["other", "fcc", "bcc", "hcp"].iter().map(|s| s.to_string()).collect(), types: BTreeMap::new() };
    let mut seen: BTreeMap<String, [i64; 4]> = BTreeMap::new();
    let mut sampled = [0i64; 4];
    let mut random = Mt19937::new(12345);
    let mut gauss = Normal::new(0.0, 1.0);
    for kind in 1..=3usize {
        let (unit, basis, repeat, nearest): (Mat3, Vec<Vec3>, [i32; 3], f64) = match kind {
            1 => (identity(), vec![[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]], [4, 4, 4], 0.5f64.sqrt()),
            2 => (identity(), vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], [5, 5, 5], 0.75f64.sqrt()),
            _ => (
                [[1.0, 0.0, 0.0], [-0.5, 3.0f64.sqrt() / 2.0, 0.0], [0.0, 0.0, (8.0f64 / 3.0).sqrt()]],
                vec![[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.5]],
                [5, 5, 3],
                1.0,
            ),
        };
        let mut crystal = AtomSet { positions: Vec::new(), cell: [[0.0; 3]; 3], pbc: [true; 3], has_cell: true };
        for i in 0..repeat[0] {
            for j in 0..repeat[1] {
                for k in 0..repeat[2] {
                    for b in &basis {
                        crystal.positions.push(row_times([b[0] + i as f64, b[1] + j as f64, b[2] + k as f64], &unit));
                    }
                }
            }
        }
        for r in 0..3 {
            crystal.cell[r] = scale(unit[r], repeat[r] as f64);
        }
        let ideal = crystal.positions.clone();
        for sample in 0..samples {
            let level = noise * (0.5 + sample as f64 / (samples - 1) as f64);
            for i in 0..ideal.len() {
                for d in 0..3 {
                    crystal.positions[i][d] = ideal[i][d] + level * nearest * gauss.sample(&mut random);
                }
            }
            for code in cell_codes(&crystal)? {
                if code.is_empty() {
                    continue;
                }
                seen.entry(code_key(&code)).or_insert([0; 4])[kind] += 1;
                sampled[kind] += 1;
            }
        }
    }
    for (k, counts) in &seen {
        let (mut best, mut best_frequency) = (0, 0.0);
        for kind in 1..=3usize {
            let frequency = counts[kind] as f64 / 1i64.max(sampled[kind]) as f64;
            if frequency > best_frequency {
                best_frequency = frequency;
                best = kind as i32;
            }
        }
        if best != 0 {
            filter.types.insert(k.clone(), best);
        }
    }
    CACHE.get().unwrap().lock().unwrap().insert(key, filter.clone());
    Ok(filter)
}

/// A VoroTop filter file: "* id name" lines name the types, "id (codes)" lines list codes.
fn read_filter(path: &Path) -> Result<Filter, String> {
    let text = read_text_file(path)?;
    let mut filter = Filter { names: vec!["other".into()], types: BTreeMap::new() };
    let mut names: BTreeMap<i32, String> = BTreeMap::new();
    let mut codes: Vec<(i32, String)> = Vec::new();
    for line in text.split('\n') {
        let Some(start) = line.find(|c: char| c != ' ' && c != '\t' && c != '\r') else { continue };
        let rest = &line[start..];
        if rest.starts_with('#') {
            continue;
        }
        if let Some(stripped) = rest.strip_prefix('*') {
            let mut words = stripped.split_ascii_whitespace();
            if let (Some(id), Some(name)) = (words.next(), words.next()) {
                if let Ok(id) = id.parse::<i32>() {
                    names.insert(id, name.to_string());
                }
            }
            continue;
        }
        let (Some(open), Some(close)) = (line.find('('), line.find(')')) else { continue };
        if close < open {
            continue;
        }
        let id = atoi(&line[start..open]);
        let inner: String = line[open + 1..close].chars().filter(|c| !c.is_ascii_whitespace()).collect();
        codes.push((id, inner));
    }
    if codes.is_empty() {
        return Err(format!("No codes in the VoroTop filter {}", path.to_string_lossy()));
    }
    let mut index: BTreeMap<i32, i32> = BTreeMap::new();
    for (id, name) in &names {
        index.insert(*id, filter.names.len() as i32);
        filter.names.push(name.clone());
    }
    for (id, code) in codes {
        if !index.contains_key(&id) {
            index.insert(id, filter.names.len() as i32);
            filter.names.push(format!("type {id}"));
        }
        filter.types.insert(code, index[&id]);
    }
    Ok(filter)
}

/// C atoi: leading whitespace, sign, digits.
fn atoi(text: &str) -> i32 {
    let t = text.trim_start();
    let (negative, digits) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let mut value: i64 = 0;
    for b in digits.bytes() {
        if !b.is_ascii_digit() {
            break;
        }
        value = (value * 10 + (b - b'0') as i64).min(i32::MAX as i64 + 1);
    }
    (if negative { -value } else { value }) as i32
}

// ---------------------------------------------------------------- CNA
struct Signature {
    common: i32,
    bonds: i32,
    chain: i32,
}

fn signature(v: &[Vec3], j: usize, cutoff: f64) -> Signature {
    let common: Vec<usize> = (0..v.len()).filter(|&k| k != j && norm(sub(v[k], v[j])) < cutoff).collect();
    let mut bonds = Vec::new();
    for a in 0..common.len() {
        for b in a + 1..common.len() {
            if norm(sub(v[common[a]], v[common[b]])) < cutoff {
                bonds.push((a, b));
            }
        }
    }
    let mut parent: Vec<usize> = (0..common.len()).collect();
    let find = |parent: &mut Vec<usize>, mut x: usize| {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    };
    for &(a, b) in &bonds {
        let rb = find(&mut parent, b);
        let ra = find(&mut parent, a);
        parent[ra] = rb;
    }
    let mut per_cluster: BTreeMap<usize, i32> = BTreeMap::new();
    for &(a, _) in &bonds {
        let root = find(&mut parent, a);
        *per_cluster.entry(root).or_insert(0) += 1;
    }
    let chain = per_cluster.values().fold(0, |m, &c| m.max(c));
    Signature { common: common.len() as i32, bonds: bonds.len() as i32, chain }
}

fn cna_type(v: &[Vec3], cutoff: f64) -> i32 {
    let (mut n421, mut n422, mut n444, mut n555, mut n666) = (0, 0, 0, 0, 0);
    for j in 0..v.len() {
        let s = signature(v, j, cutoff);
        match (s.common, s.bonds, s.chain) {
            (4, 2, 1) => n421 += 1,
            (4, 2, 2) => n422 += 1,
            (4, 4, 4) => n444 += 1,
            (5, 5, 5) => n555 += 1,
            (6, 6, 6) => n666 += 1,
            _ => {}
        }
    }
    if v.len() == 12 && n421 == 12 {
        1
    } else if v.len() == 12 && n421 == 6 && n422 == 6 {
        2
    } else if v.len() == 12 && n555 == 12 {
        4
    } else if v.len() == 14 && n444 == 6 && n666 == 8 {
        3
    } else {
        0
    }
}

pub fn voronoi_analysis(p: &Parameters) -> R {
    let a = read_atom_set(p, "positions")?;
    let n = a.positions.len();
    if n == 0 {
        return Err("No atoms".into());
    }
    let area_threshold = p.number_or("face_area_threshold_A2", 0.0)?;
    let relative_threshold = p.number_or("relative_face_area_threshold", 0.0)?;
    let edge_threshold = p.number_or("edge_length_threshold_A", 0.0)?;
    if area_threshold < 0.0 || relative_threshold < 0.0 || edge_threshold < 0.0 {
        return Err("Thresholds must be nonnegative".into());
    }
    let mut cutoff = p.number_or("cutoff_A", 0.0)?;
    let cells = all_cells(&a, &mut cutoff)?;
    let (mut volume, mut max_order) = (vec![f64::NAN; n], vec![f64::NAN; n]);
    let mut coordination = vec![0i32; n];
    let mut open = vec![false; n];
    let mut indices = Vec::new();
    let mut index_counts: BTreeMap<String, i64> = BTreeMap::new();
    let mut total = 0.0;
    for i in 0..n {
        let cell = &cells[i];
        open[i] = cell.open;
        let mut index = vec![0i32; 10];
        let mut surface = 0.0;
        for face in &cell.faces {
            surface += polygon_area(&face.points);
        }
        let mut v = 0.0;
        let mut largest = 0;
        for face in &cell.faces {
            let area = polygon_area(&face.points);
            let mut signed6 = 0.0;
            for k in 1..face.points.len().saturating_sub(1) {
                signed6 += dot(face.points[0], cross(face.points[k], face.points[k + 1]));
            }
            v += signed6.abs() / 6.0;
            if face.neighbor < 0 || area <= area_threshold || area <= relative_threshold * surface {
                continue;
            }
            coordination[i] += 1;
            let m = face.points.len();
            let edges = (0..m).filter(|&k| norm(sub(face.points[(k + 1) % m], face.points[k])) >= edge_threshold).count() as i32;
            if (3..=12).contains(&edges) {
                index[(edges - 3) as usize] += 1;
            }
            largest = largest.max(edges);
        }
        if !cell.open {
            volume[i] = v;
            total += v;
        }
        max_order[i] = largest as f64;
        indices.push(ji(&index));
        *index_counts.entry(format!("<{},{},{},{}>", index[0], index[1], index[2], index[3])).or_insert(0) += 1;
    }
    let mut ranked: Vec<(i64, String)> = index_counts.into_iter().map(|(k, c)| (c, k)).collect();
    ranked.sort_by(|x, y| y.0.cmp(&x.0).then_with(|| x.1.cmp(&y.1)));
    let common: Vec<Json> = ranked
        .iter()
        .take(20)
        .map(|(count, key)| {
            let mut item = Json::object();
            item.set("index", key.as_str());
            item.set("atoms", *count as f64);
            item
        })
        .collect();
    let mut result = Json::object();
    result.set("atomic_volume_A3", jv(&volume));
    result.set("coordination", ji(&coordination));
    result.set("voronoi_index", Json::Array(indices));
    result.set("index_columns", "faces with 3, 4, ..., 12 edges");
    result.set("max_face_order", jv(&max_order));
    result.set("open_cell", Json::Array(open.iter().map(|&b| Json::Bool(b)).collect()));
    result.set("most_common_indices", Json::Array(common));
    result.set("total_volume_A3", total);
    if a.has_cell {
        result.set("cell_volume_A3", cell_volume(&a.cell));
    }
    result.set("search_radius_A", cutoff);
    Ok(result)
}

pub fn voro_top(p: &Parameters) -> R {
    let a = read_atom_set(p, "positions")?;
    let n = a.positions.len();
    if n == 0 {
        return Err("No atoms".into());
    }
    let noise = p.number_or("noise", 0.06)?;
    if !(noise > 0.0 && noise < 0.3) {
        return Err("noise is a fraction of the nearest-neighbour distance between 0 and 0.3".into());
    }
    let filter = if p.has("filter_file") {
        read_filter(Path::new(p.json("filter_file")?.string()?))?
    } else {
        built_in_filter(noise, integer(p.number_or("samples", 60.0)?, "samples", 2)? as i32)?
    };
    let codes = cell_codes(&a)?;
    let (mut kind, mut topology) = (vec![0i32; n], vec![-1i32; n]);
    let mut distinct: BTreeMap<String, i32> = BTreeMap::new();
    let mut topologies: Vec<(String, i64)> = Vec::new();
    for i in 0..n {
        if codes[i].is_empty() {
            continue;
        }
        let k = code_key(&codes[i]);
        if let Some(&t) = filter.types.get(&k) {
            kind[i] = t;
        }
        let entry = *distinct.entry(k.clone()).or_insert_with(|| {
            topologies.push((k.clone(), 0));
            topologies.len() as i32 - 1
        });
        topology[i] = entry;
        topologies[entry as usize].1 += 1;
    }
    let mut result = Json::object();
    let names: Vec<&str> = filter.names.iter().map(|s| s.as_str()).collect();
    add_classification(&mut result, "structure_type", &kind, &names, "");
    let mut order: Vec<usize> = (0..topologies.len()).collect();
    sort_by(&mut order, |&x, &y| topologies[x].1 > topologies[y].1);
    let common: Vec<Json> = order
        .iter()
        .take(10)
        .map(|&k| {
            let (code, count) = &topologies[k];
            let mut item = Json::object();
            item.set("atoms", *count as f64);
            item.set("weinberg_code", format!("({code})").as_str());
            item.set("type", filter.names[*filter.types.get(code).unwrap_or(&0) as usize].as_str());
            item
        })
        .collect();
    result.set("topology_id", ji(&topology));
    result.set("distinct_topologies", topologies.len() as f64);
    result.set("cells_without_topology", codes.iter().filter(|c| c.is_empty()).count() as f64);
    result.set("most_common_topologies", Json::Array(common));
    let filter_name = if p.has("filter_file") { p.json("filter_file")?.string()?.to_string() } else { "built-in (perturbed fcc, bcc, hcp)".to_string() };
    result.set("filter", filter_name.as_str());
    result.set("filter_codes", filter.types.len() as f64);
    if p.boolean("export_filter", false)? {
        let mut text = "# AtomForge VoroTop filter\n".to_string();
        for t in 1..filter.names.len() {
            text += &format!("*\t{}\t{}\n", t, filter.names[t]);
        }
        for (code, t) in &filter.types {
            text += &format!("{t}\t({code})\n");
        }
        let mut files = Json::object();
        files.set("vorotop.filter", text.as_str());
        result.set("files", files);
    }
    Ok(result)
}

pub fn adaptive_cna(p: &Parameters) -> R {
    let a = read_atom_set(p, "positions")?;
    let n = a.positions.len();
    if n == 0 {
        return Err("No atoms".into());
    }
    let fixed = p.number_or("cutoff_A", 0.0)?;
    if fixed < 0.0 {
        return Err("cutoff_A must be nonnegative".into());
    }
    let factor = (1.0 + 2.0f64.sqrt()) / 2.0;
    let neighbors = sorted_neighbor_vectors(&a, if fixed > 0.0 { fixed } else { 1.8 * mean_atomic_spacing(&a)? })?;
    let mut kind = vec![0i32; n];
    for i in 0..n {
        let v = &neighbors[i];
        if fixed > 0.0 {
            kind[i] = cna_type(v, fixed);
            continue;
        }
        if v.len() >= 12 {
            let mut mean = 0.0;
            for k in 0..12 {
                mean += norm(v[k]) / 12.0;
            }
            kind[i] = cna_type(&v[..12], factor * mean);
        }
        if kind[i] == 0 && v.len() >= 14 {
            let (mut first, mut second) = (0.0, 0.0);
            for k in 0..8 {
                first += norm(v[k]) / 8.0;
            }
            for k in 8..14 {
                second += norm(v[k]) / 6.0;
            }
            let cutoff = factor * (2.0 / 3.0f64.sqrt() * first + second) / 2.0;
            if cna_type(&v[..14], cutoff) == 3 {
                kind[i] = 3;
            }
        }
    }
    let mut result = Json::object();
    add_classification(&mut result, "structure_type", &kind, &["other", "fcc", "hcp", "bcc", "icosahedral"], "");
    result.set("mode", if fixed > 0.0 { "fixed cutoff" } else { "adaptive" });
    Ok(result)
}
