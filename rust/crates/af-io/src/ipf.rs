//! Crystal-orientation (IPF) colours: the .atomforge-ipf sidecar written
//! beside saved structures, and the reconstruction from the geometry when
//! there is none (the C++ StructureLoader IPF functions, in single precision
//! and with libstdc++'s sort order so colours equal the C++ ones).

use std::collections::{HashMap, VecDeque};
use std::path::Path;

use af_core::format::format_g;
use af_core::Structure;
use glam::{IVec3, Mat3, Vec3};

pub const SIDECAR_SUFFIX: &str = ".atomforge-ipf";

pub fn sidecar_path(filename: &str) -> String {
    let slash = filename.rfind(['/', '\\']);
    let dot = filename.rfind('.');
    match dot {
        Some(d) if slash.map_or(true, |s| d > s) => format!("{}{SIDECAR_SUFFIX}", &filename[..d]),
        _ => format!("{filename}{SIDECAR_SUFFIX}"),
    }
}

fn legacy_sidecar_path(filename: &str) -> String {
    format!("{filename}{SIDECAR_SUFFIX}")
}

pub fn remove_sidecar(filename: &str) {
    let path = sidecar_path(filename);
    let legacy = legacy_sidecar_path(filename);
    let _ = std::fs::remove_file(&path);
    if legacy != path {
        let _ = std::fs::remove_file(&legacy);
    }
}

/// Writes the colours (and regions) beside a saved file; colours that came
/// from the geometry are not written (they are rebuilt on the next load).
pub fn save_sidecar(structure: &Structure, filename: &str) -> bool {
    let reconstructed = structure.ipf_load_status.starts_with("IPF reconstructed") && structure.grain_region_ids.is_empty();
    if structure.grain_colors.len() != structure.atoms.len() || structure.atoms.is_empty() || reconstructed {
        remove_sidecar(filename);
        return true;
    }
    let regions = structure.grain_region_ids.len() == structure.atoms.len();
    let mut out = String::from(if regions { "ATOMFORGE_IPF_V2\n" } else { "ATOMFORGE_IPF_V1\n" });
    out.push_str(&format!("{}\n", structure.atoms.len()));
    for (i, atom) in structure.atoms.iter().enumerate() {
        let c = structure.grain_colors[i];
        out.push_str(&format!(
            "{} {} {} {} {} {} {}",
            atom.atomic_number,
            format_g(atom.x, 17),
            format_g(atom.y, 17),
            format_g(atom.z, 17),
            format_g(c[0] as f64, 17),
            format_g(c[1] as f64, 17),
            format_g(c[2] as f64, 17)
        ));
        if regions {
            out.push_str(&format!(" {}", structure.grain_region_ids[i]));
        }
        out.push('\n');
    }
    // A C++ text-mode stream: CRLF line ends on Windows.
    let out = if cfg!(windows) { out.replace('\n', "\r\n") } else { out };
    std::fs::write(sidecar_path(filename), out).is_ok()
}

fn quantize(value: f64) -> i64 {
    (value * 500.0).floor() as i64
}

struct Record {
    key: (i32, i64, i64, i64),
    position: [f64; 3],
    region: i32,
    has_region: bool,
    colour: [f32; 3],
}

fn load_records(filename: &str) -> Option<Vec<Record>> {
    let text = std::fs::read_to_string(sidecar_path(filename)).or_else(|_| {
        let legacy = legacy_sidecar_path(filename);
        if legacy == sidecar_path(filename) {
            Err(std::io::Error::other("none"))
        } else {
            std::fs::read_to_string(legacy)
        }
    });
    let text = text.ok()?;
    // Written in text mode by the C++ on Windows: CRLF line ends.
    let mut lines = text.splitn(2, '\n');
    let header = lines.next()?.trim_end_matches('\r').to_string();
    let v1 = header == "ATOMFORGE_IPF_V1";
    let v2 = header == "ATOMFORGE_IPF_V2";
    if !v1 && !v2 {
        return None;
    }
    let mut words = lines.next().unwrap_or("").split_whitespace();
    let count: usize = words.next()?.parse().ok()?;
    if count > 10_000_000 {
        return None;
    }
    let mut records = Vec::with_capacity(count);
    for _ in 0..count {
        let z: i32 = words.next()?.parse().ok()?;
        let mut p = [0.0f64; 3];
        for v in p.iter_mut() {
            *v = words.next()?.parse().ok()?;
        }
        let mut c = [0.0f32; 3];
        for v in c.iter_mut() {
            *v = words.next()?.parse().ok()?;
        }
        let region = if v2 { words.next()?.parse().ok()? } else { -1 };
        if !(1..=118).contains(&z) || p.iter().any(|v| !v.is_finite() || v.abs() > 1e9) || c.iter().any(|v| !v.is_finite()) {
            return None;
        }
        records.push(Record { key: (z, quantize(p[0]), quantize(p[1]), quantize(p[2])), position: p, region, has_region: v2, colour: c });
    }
    Some(records)
}

/// Colours (and regions) from a sidecar written for this file, when every
/// atom matches a record within 0.002 A.
pub fn restore_sidecar(filename: &str, structure: &mut Structure) -> bool {
    if structure.atoms.is_empty() {
        return false;
    }
    let Some(records) = load_records(filename) else { return false };
    if records.len() != structure.atoms.len() {
        return false;
    }
    let mut by_key: HashMap<(i32, i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, a) in structure.atoms.iter().enumerate() {
        by_key.entry((a.atomic_number, quantize(a.x), quantize(a.y), quantize(a.z))).or_default().push(i);
    }
    let n = structure.atoms.len();
    let mut restored = vec![[0.0f32; 3]; n];
    let mut regions = vec![-1i32; n];
    let mut assigned = vec![false; n];
    let mut matched = 0;
    for record in &records {
        let mut chosen = None;
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let key = (record.key.0, record.key.1 + dx, record.key.2 + dy, record.key.3 + dz);
                    for &index in by_key.get(&key).into_iter().flatten() {
                        let a = &structure.atoms[index];
                        let p = record.position;
                        if !assigned[index] && (a.x - p[0]).abs() < 0.002 && (a.y - p[1]).abs() < 0.002 && (a.z - p[2]).abs() < 0.002 {
                            chosen = Some(index);
                            break 'search;
                        }
                    }
                }
            }
        }
        let Some(index) = chosen else { break };
        restored[index] = record.colour;
        regions[index] = record.region;
        assigned[index] = true;
        matched += 1;
    }
    if matched != n {
        return false;
    }
    structure.grain_colors = restored;
    if records.first().is_some_and(|r| r.has_region) {
        structure.grain_region_ids = regions;
    }
    true
}

/// The IPF colour of a crystal direction (cubic standard triangle).
pub fn colour_from_direction(direction: Vec3) -> [f32; 3] {
    let (mut h, mut k, mut l) = (direction.x.abs(), direction.y.abs(), direction.z.abs());
    if h < k {
        std::mem::swap(&mut h, &mut k);
    }
    if h < l {
        std::mem::swap(&mut h, &mut l);
    }
    if k < l {
        std::mem::swap(&mut k, &mut l);
    }
    let len = (h * h + k * k + l * l).sqrt();
    if len < 1e-10 {
        return [0.5, 0.5, 0.5];
    }
    h /= len;
    k /= len;
    l /= len;
    let (mut r, mut g, mut b) = (h - k, k - l, l * 1.7320508);
    let max = r.max(g).max(b).max(1e-6);
    r /= max;
    g /= max;
    b /= max;
    [r, g, b]
}

/// glm::normalize: v * inversesqrt(dot(v, v)).
fn normalize(v: Vec3) -> Vec3 {
    v * (1.0 / v.dot(v).sqrt())
}

/// Orientation colours from the local neighbour frames, smoothed and
/// segmented into regions (for structures of 8 or more atoms).
pub fn recover_from_geometry(structure: &mut Structure) -> bool {
    let n = structure.atoms.len();
    if n < 8 {
        return false;
    }
    const CELL: f32 = 4.0;
    const MAX_DIST: f32 = 4.2;
    const MIN_DIST: f32 = 0.2;
    const TARGET: usize = 16;
    let positions: Vec<Vec3> = structure.atoms.iter().map(|a| Vec3::new(a.x as f32, a.y as f32, a.z as f32)).collect();
    let grid_cell = |p: Vec3| IVec3::new((p.x / CELL).floor() as i32, (p.y / CELL).floor() as i32, (p.z / CELL).floor() as i32);
    let mut grid: HashMap<IVec3, Vec<usize>> = HashMap::with_capacity(n);
    for (i, &p) in positions.iter().enumerate() {
        grid.entry(grid_cell(p)).or_default().push(i);
    }
    let mut dirs = vec![Vec3::Z; n];
    let mut neighbour_ids: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut has_dir = vec![false; n];
    let mut recovered = 0usize;
    for i in 0..n {
        let c = grid_cell(positions[i]);
        let mut neighbours: Vec<(f32, usize)> = Vec::with_capacity(64);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(list) = grid.get(&(c + IVec3::new(dx, dy, dz))) {
                        for &idx in list {
                            if idx == i {
                                continue;
                            }
                            let dist = (positions[idx] - positions[i]).length();
                            if !(MIN_DIST..=MAX_DIST).contains(&dist) {
                                continue;
                            }
                            neighbours.push((dist, idx));
                        }
                    }
                }
            }
        }
        if neighbours.len() < 3 {
            continue;
        }
        af_core::cxx::sort_by(&mut neighbours, |a, b| a.0 < b.0);
        let use_count = neighbours.len().min(TARGET);
        let mut units = Vec::with_capacity(use_count);
        for &(_, ni) in &neighbours[..use_count] {
            let d = positions[ni] - positions[i];
            let dn = d.length();
            if dn < 1e-8 {
                continue;
            }
            neighbour_ids[i].push(ni);
            units.push(d / dn);
        }
        if units.len() < 3 {
            continue;
        }
        let u1 = normalize(units[0]);
        let Some(u2raw) = units[1..].iter().copied().find(|u| u1.dot(*u).abs() < 0.85) else { continue };
        let mut u2 = u2raw - u2raw.dot(u1) * u1;
        let u2n = u2.length();
        if u2n < 1e-6 {
            continue;
        }
        u2 /= u2n;
        let mut u3 = u1.cross(u2);
        let u3n = u3.length();
        if u3n < 1e-6 {
            continue;
        }
        u3 /= u3n;
        let r = Mat3::from_cols(u1, u2, u3);
        let crystal = r.transpose() * Vec3::Z;
        if crystal.length() < 1e-6 {
            continue;
        }
        dirs[i] = normalize(crystal);
        has_dir[i] = true;
        recovered += 1;
    }
    if recovered < n / 4 {
        return false;
    }
    for _ in 0..3 {
        let mut changed = false;
        for i in 0..n {
            if has_dir[i] {
                continue;
            }
            let mut sum = Vec3::ZERO;
            for &ni in &neighbour_ids[i] {
                if !has_dir[ni] {
                    continue;
                }
                let mut a = dirs[ni];
                if a.dot(sum) < 0.0 {
                    a = -a;
                }
                sum += a;
            }
            if sum.length() > 1e-8 {
                dirs[i] = normalize(sum);
                has_dir[i] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    for i in 0..n {
        if !has_dir[i] {
            dirs[i] = Vec3::Z;
        }
    }
    for _ in 0..4 {
        let mut next = dirs.clone();
        for i in 0..n {
            let mut sum = dirs[i] * 2.0;
            for &ni in &neighbour_ids[i] {
                let mut a = dirs[ni];
                if a.dot(sum) < 0.0 {
                    a = -a;
                }
                sum += a;
            }
            if sum.length() > 1e-8 {
                next[i] = normalize(sum);
            }
        }
        dirs = next;
    }
    let cos_threshold = (20.0f32 * 3.1415926535 / 180.0).cos();
    let mut region_id = vec![-1i32; n];
    let mut means: Vec<Vec3> = Vec::new();
    let mut members: Vec<Vec<usize>> = Vec::new();
    for seed in 0..n {
        if region_id[seed] >= 0 {
            continue;
        }
        let rid = means.len();
        means.push(normalize(dirs[seed]));
        members.push(vec![seed]);
        region_id[seed] = rid as i32;
        let mut queue = VecDeque::from([seed]);
        let mut sum = dirs[seed];
        while let Some(u) = queue.pop_front() {
            let reference = if sum.length() > 1e-8 { normalize(sum) } else { normalize(dirs[u]) };
            for &v in &neighbour_ids[u] {
                if region_id[v] >= 0 {
                    continue;
                }
                let dv = normalize(dirs[v]);
                if dv.dot(reference).abs() < cos_threshold {
                    continue;
                }
                region_id[v] = rid as i32;
                members[rid].push(v);
                queue.push_back(v);
                let mut aligned = dv;
                if aligned.dot(sum) < 0.0 {
                    aligned = -aligned;
                }
                sum += aligned;
            }
        }
        if sum.length() > 1e-8 {
            means[rid] = normalize(sum);
        }
    }
    let min_size = 32usize.max(n / 180);
    let mut large: Vec<bool> = members.iter().map(|m| m.len() >= min_size).collect();
    if !large.iter().any(|&l| l) && !members.is_empty() {
        let mut largest = 0;
        for r in 1..members.len() {
            if members[r].len() > members[largest].len() {
                largest = r;
            }
        }
        large[largest] = true;
    }
    let mut remap: Vec<usize> = (0..means.len()).collect();
    // Adjacency counts of small regions to large ones, in first-seen order (the
    // C++ unordered_map order differs only on exact ties of count and angle).
    let mut adjacency: Vec<Vec<(usize, i32)>> = vec![Vec::new(); means.len()];
    let bump = |from: usize, to: usize, adjacency: &mut Vec<Vec<(usize, i32)>>| {
        if let Some(entry) = adjacency[from].iter_mut().find(|e| e.0 == to) {
            entry.1 += 1;
        } else {
            adjacency[from].push((to, 1));
        }
    };
    for i in 0..n {
        let ri = region_id[i];
        if ri < 0 {
            continue;
        }
        let ri = ri as usize;
        for &ni in &neighbour_ids[i] {
            if ni <= i {
                continue;
            }
            let rj = region_id[ni];
            if rj < 0 || rj as usize == ri {
                continue;
            }
            let rj = rj as usize;
            if !large[ri] && large[rj] {
                bump(ri, rj, &mut adjacency);
            }
            if !large[rj] && large[ri] {
                bump(rj, ri, &mut adjacency);
            }
        }
    }
    for r in 0..means.len() {
        if large[r] {
            continue;
        }
        let mut best: i64 = -1;
        let mut best_count = -1;
        let mut best_orient = -1.0f32;
        for &(nr, count) in &adjacency[r] {
            let orient = means[r].dot(means[nr]).abs();
            if count > best_count || (count == best_count && orient > best_orient) {
                best = nr as i64;
                best_count = count;
                best_orient = orient;
            }
        }
        if best < 0 {
            for nr in 0..means.len() {
                if !large[nr] {
                    continue;
                }
                let orient = means[r].dot(means[nr]).abs();
                if orient > best_orient {
                    best_orient = orient;
                    best = nr as i64;
                }
            }
        }
        if best >= 0 {
            remap[r] = best as usize;
        }
    }
    for id in region_id.iter_mut() {
        if *id >= 0 {
            *id = remap[*id as usize] as i32;
        }
    }
    let mut final_means = vec![Vec3::ZERO; means.len()];
    let mut counts = vec![0i32; means.len()];
    for i in 0..n {
        let rid = region_id[i];
        if rid < 0 {
            continue;
        }
        let rid = rid as usize;
        let mut d = normalize(dirs[i]);
        if counts[rid] > 0 && d.dot(final_means[rid]) < 0.0 {
            d = -d;
        }
        final_means[rid] += d;
        counts[rid] += 1;
    }
    for r in 0..final_means.len() {
        final_means[r] = if counts[r] > 0 && final_means[r].length() > 1e-8 { normalize(final_means[r]) } else { Vec3::Z };
    }
    structure.grain_colors = (0..n)
        .map(|i| {
            let rid = if region_id[i] >= 0 { region_id[i] as usize } else { 0 };
            colour_from_direction(final_means.get(rid).copied().unwrap_or(Vec3::Z))
        })
        .collect();
    true
}

pub fn exists(filename: &str) -> bool {
    Path::new(&sidecar_path(filename)).exists() || Path::new(&legacy_sidecar_path(filename)).exists()
}
