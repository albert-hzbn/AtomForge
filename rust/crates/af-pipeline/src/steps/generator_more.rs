//! The rest of the C++ pipeline/GeneratorSteps: adsorbates, near-cubic
//! supercells, interpolation (linear and IDPP) and point defects (antisites,
//! Frenkel pairs). The types are listed in generator.rs.

use std::path::Path;

use af_core::cxx_random::{shuffle, uniform_u64, Mt19937_64};
use af_core::json::Json;
use af_core::Structure;
use af_science::core::{add, cross, determinant, dot, inverse, norm, row_times, scale, sub, Mat3, Vec3};
use af_science::data::read_frames;

use crate::molecule_library::molecule;
use crate::pipeline::{find_modifier_type, PipelineData};
use crate::steps::basic::cxx_u64;
use crate::steps::generator::{element_list, g5, principal_frame, rotated, R, PI};
use crate::structure_edits::{append_atom, rotate, rotation_about, rotation_between};
use crate::support::*;

fn fmin(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

fn fmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

// ---------------------------------------------------------------- adsorbates
pub(crate) fn add_adsorbate(p: &Json, data: &mut PipelineData) -> R {
    let id = "add-adsorbate";
    require_cell(data, "add-adsorbate")?;
    if data.structure.atoms.is_empty() {
        return Err("add-adsorbate needs a slab".into());
    }
    let species_name = text(p, id, "species")?;
    let mut species = molecule(&species_name)?;
    let site = text(p, id, "site")?;
    if site != "top" && site != "bridge" && site != "fcc" && site != "hcp" {
        return Err("site must be top, bridge, fcc or hcp".into());
    }
    let height = number(p, id, "height")?;
    let target = vector(p, id, "position")?;
    let coverage = text(p, id, "coverage")?;
    let side = text(p, id, "surface")?;
    if coverage != "one" && coverage != "all" && coverage != "fraction" {
        return Err("coverage must be one, all or fraction".into());
    }
    if side != "top" && side != "bottom" && side != "both" {
        return Err("surface must be top, bottom or both".into());
    }
    let anchor_index = integer_of(p, id, "anchor", -1)?;
    if anchor_index >= species.len() as i32 {
        return Err(format!("anchor: the molecule has {} atoms", species.len()));
    }
    if flag(p, id, "principal")? && species.len() > 1 {
        let frame = principal_frame(&species);
        species = rotated(species, &frame);
    }
    let tilt = number(p, id, "tilt")? * PI / 180.0;
    let azimuth = number(p, id, "azimuth")? * PI / 180.0;
    species = rotated(species, &rotation_about([0.0, 1.0, 0.0], tilt));
    species = rotated(species, &rotation_about([0.0, 0.0, 1.0], azimuth));
    let s = &data.structure;
    let cell = cell_of(s);
    let mut up = cross(cell[0], cell[1]);
    up = scale(up, 1.0 / norm(up));
    if dot(up, cell[2]) < 0.0 {
        up = scale(up, -1.0);
    }
    let inv = inverse(&cell)?;
    let in_plane = |v: Vec3| {
        let v = sub(v, scale(up, dot(v, up)));
        let mut f = row_times(v, &inv);
        for k in 0..2 {
            f[k] -= f[k].round();
        }
        f[2] = 0.0;
        row_times(f, &cell)
    };
    let layer_tolerance = 0.6;
    let mut placements: Vec<(Vec3, Vec3, f64)> = Vec::new(); // site, normal, surface
    let mut site_count = 0usize;
    let mut extra_notes = Vec::new();
    for sign in [1.0, -1.0] {
        if (sign > 0.0 && side == "bottom") || (sign < 0.0 && side == "top") {
            continue;
        }
        let normal = scale(up, sign);
        let mut top = f64::NEG_INFINITY;
        for a in &s.atoms {
            top = fmax(top, dot(position_of(a), normal));
        }
        let (mut surface, mut second) = (Vec::new(), Vec::new());
        let mut second_height = f64::NEG_INFINITY;
        for a in &s.atoms {
            let h = dot(position_of(a), normal);
            if h < top - layer_tolerance {
                second_height = fmax(second_height, h);
            }
        }
        for (i, a) in s.atoms.iter().enumerate() {
            let h = dot(position_of(a), normal);
            if h >= top - layer_tolerance {
                surface.push(i);
            } else if h >= second_height - layer_tolerance {
                second.push(i);
            }
        }
        let lateral = |x: usize, y: usize, i: i32, j: i32| {
            let mut v = sub(position_of(&s.atoms[y]), position_of(&s.atoms[x]));
            v = add(v, add(scale(cell[0], i as f64), scale(cell[1], j as f64)));
            sub(v, scale(up, dot(v, up)))
        };
        let mut nearest = f64::INFINITY;
        for &x in &surface {
            for &y in &surface {
                for i in -1..=1 {
                    for j in -1..=1 {
                        let v = lateral(x, y, i, j);
                        if norm(v) > 1e-6 {
                            nearest = fmin(nearest, norm(v));
                        }
                    }
                }
            }
        }
        if site != "top" && !(nearest < f64::INFINITY) {
            return Err(format!("A {site} site needs a periodic surface"));
        }
        let mut candidates: Vec<Vec3> = Vec::new();
        let add_candidate = |candidates: &mut Vec<Vec3>, q: Vec3| {
            if !candidates.iter().any(|c| norm(in_plane(sub(*c, q))) < 0.1) {
                candidates.push(q);
            }
        };
        for &x in &surface {
            let px = position_of(&s.atoms[x]);
            if site == "top" {
                add_candidate(&mut candidates, px);
                continue;
            }
            let mut near: Vec<Vec3> = Vec::new();
            for &y in &surface {
                for i in -1..=1 {
                    for j in -1..=1 {
                        let v = lateral(x, y, i, j);
                        let len = norm(v);
                        if len > 1e-6 && len < 1.2 * nearest {
                            near.push(v);
                        }
                    }
                }
            }
            if site == "bridge" {
                for v in &near {
                    add_candidate(&mut candidates, add(px, scale(*v, 0.5)));
                }
                continue;
            }
            for m in 0..near.len() {
                for q in m + 1..near.len() {
                    if norm(sub(near[m], near[q])) > 1.2 * nearest {
                        continue;
                    }
                    let centroid = add(px, scale(add(near[m], near[q]), 1.0 / 3.0));
                    let mut above = false;
                    for &k in &second {
                        above = above || norm(in_plane(sub(position_of(&s.atoms[k]), centroid))) < 0.35 * nearest;
                    }
                    if (site == "hcp") == above {
                        add_candidate(&mut candidates, centroid);
                    }
                }
            }
        }
        if candidates.is_empty() {
            return Err(format!("No {site} site found on the {} surface", if sign > 0.0 { "top" } else { "bottom" }));
        }
        site_count += candidates.len();
        if coverage == "one" {
            let wanted = row_times([target[0], target[1], 0.0], &cell);
            let (mut best, mut best_distance) = (0usize, f64::INFINITY);
            for (k, c) in candidates.iter().enumerate() {
                let dist = norm(in_plane(sub(*c, wanted)));
                if dist < best_distance - 1e-9 {
                    best_distance = dist;
                    best = k;
                }
            }
            candidates = vec![candidates[best]];
        } else if coverage == "fraction" {
            let fraction = number(p, id, "fraction")?;
            let separation = number(p, id, "min_separation")?;
            if !(0.0..=1.0).contains(&fraction) {
                return Err("fraction must lie between 0 and 1".into());
            }
            let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?).wrapping_add(if sign > 0.0 { 0 } else { 7919 }));
            shuffle(&mut candidates, &mut rng);
            let wanted = (fraction * candidates.len() as f64).round() as usize;
            let mut chosen: Vec<Vec3> = Vec::new();
            for c in &candidates {
                if chosen.len() == wanted {
                    break;
                }
                if chosen.iter().all(|o| norm(in_plane(sub(*c, *o))) >= separation) {
                    chosen.push(*c);
                }
            }
            if chosen.len() < wanted {
                extra_notes.push(format!("Only {} of {} sites could be covered at the minimum separation", chosen.len(), wanted));
            }
            candidates = chosen;
        }
        for c in candidates {
            placements.push((c, normal, top));
        }
    }
    data.notes.extend(extra_notes);
    let n = data.structure.atoms.len();
    data.selected = vec![0; n];
    let select_new = flag(p, id, "select")?;
    for &(site_point, normal, surface) in &placements {
        let anchor = add(sub(site_point, scale(normal, dot(site_point, normal))), scale(normal, surface + height));
        let r = rotation_between([0.0, 0.0, 1.0], normal);
        let reference = if anchor_index >= 0 {
            rotate(&r, species[anchor_index as usize].position)
        } else {
            let mut lowest = f64::INFINITY;
            for a in &species {
                lowest = fmin(lowest, dot(rotate(&r, a.position), normal));
            }
            scale(normal, lowest)
        };
        for a in &species {
            append_atom(&mut data.structure, &mut data.selected, &a.symbol, add(anchor, sub(rotate(&r, a.position), reference)), select_new);
        }
    }
    let mut note = format!(
        "{} {} on {} sites ({} available), {} Angstrom above the {}",
        placements.len(),
        species_name,
        site,
        site_count,
        g5(height),
        if side == "both" { "two surfaces".to_string() } else { format!("{side} surface") }
    );
    if coverage != "one" {
        let area = norm(cross(cell[0], cell[1]));
        note += &format!("; coverage {} per nm2", g5(placements.len() as f64 / if side == "both" { 2.0 } else { 1.0 } / area * 100.0));
    }
    data.notes.push(note);
    Ok(())
}

// ---------------------------------------------------------------- supercells
pub(crate) fn near_cubic_supercell(p: &Json, data: &mut PipelineData) -> R {
    let id = "near-cubic-supercell";
    require_cell(data, id)?;
    let s = &data.structure;
    if s.atoms.is_empty() {
        return Err("near-cubic-supercell needs atoms".into());
    }
    let target = integer_of(p, id, "target_atoms", 1)?;
    let range = integer_of(p, id, "range", 1)?;
    if range > 3 {
        return Err("range at most 3 (the search grows as (2r+1)^9)".into());
    }
    let tolerance = number(p, id, "tolerance")? / 100.0;
    if !(0.0..1.0).contains(&tolerance) {
        return Err("tolerance is a percentage from 0 to 100".into());
    }
    let shape = text(p, id, "shape")?;
    let natoms = s.atoms.len() as f64;
    let nearest = 1i64.max((target as f64 / natoms).round() as i64);
    let (mut lowest, mut highest) = (nearest, nearest);
    if tolerance > 0.0 {
        lowest = 1i64.max((target as f64 * (1.0 - tolerance) / natoms).ceil() as i64);
        highest = lowest.max((target as f64 * (1.0 + tolerance) / natoms).floor() as i64);
    }
    let h = cell_of(s);
    let cell_volume = determinant(&h).abs();
    let ideal: [[f64; 3]; 3] = match shape.as_str() {
        "fcc" => [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]],
        "bcc" => [[3.0, -1.0, -1.0], [-1.0, 3.0, -1.0], [-1.0, -1.0, 3.0]],
        "sc" | "image" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        _ => return Err("shape is sc, fcc, bcc or image".into()),
    };
    let ideal_det = ideal[0][0] * (ideal[1][1] * ideal[2][2] - ideal[1][2] * ideal[2][1]) - ideal[0][1] * (ideal[1][0] * ideal[2][2] - ideal[1][2] * ideal[2][0])
        + ideal[0][2] * (ideal[1][0] * ideal[2][1] - ideal[1][1] * ideal[2][0]);
    let shortest = |q: &Mat3| {
        let mut best = f64::INFINITY;
        for i in -2..=2 {
            for j in -2..=2 {
                for k in -2..=2 {
                    if i == 0 && j == 0 && k == 0 {
                        continue;
                    }
                    let v = add(add(scale(q[0], i as f64), scale(q[1], j as f64)), scale(q[2], k as f64));
                    best = fmin(best, norm(v));
                }
            }
        }
        best
    };
    let mut best_score = f64::INFINITY;
    let (mut best_negatives, mut best_trace) = (10, -100);
    let mut best_cells = 0i64;
    let mut best = [[0i32; 3]; 3];
    let mut best_cell = [[0.0; 3]; 3];
    let span = 2 * range + 1;
    let total: i64 = (span as i64).pow(9);
    let mut e = [0i32; 9];
    for code in 0..total {
        let mut c = code;
        for k in 0..9 {
            e[k] = (c % span as i64) as i32 - range;
            c /= span as i64;
        }
        let det = e[0] as i64 * (e[4] * e[8] - e[5] * e[7]) as i64 - e[1] as i64 * (e[3] * e[8] - e[5] * e[6]) as i64 + e[2] as i64 * (e[3] * e[7] - e[4] * e[6]) as i64;
        if det < lowest || det > highest {
            continue;
        }
        let mut q = [[0.0; 3]; 3];
        for r in 0..3 {
            for col in 0..3 {
                q[r][col] = e[3 * r] as f64 * h[0][col] + e[3 * r + 1] as f64 * h[1][col] + e[3 * r + 2] as f64 * h[2][col];
            }
        }
        let volume = cell_volume * det as f64;
        let mut score = 0.0;
        if shape == "image" {
            score = 1.0 - shortest(&q) / (2.0f64.sqrt() * volume).cbrt();
        } else {
            let sc = (volume * volume / ideal_det).cbrt();
            for i in 0..3 {
                for j in 0..3 {
                    let g = dot(q[i], q[j]) - ideal[i][j] * sc;
                    score += g * g;
                }
            }
            score = score.sqrt() / (volume * volume).cbrt();
        }
        score += 1e-6 * (det as f64 - target as f64 / natoms).abs();
        let negatives = e.iter().filter(|&&x| x < 0).count() as i32;
        let trace = e[0] + e[4] + e[8];
        let better = score < best_score - 1e-9 || (score < best_score + 1e-9 && (negatives < best_negatives || (negatives == best_negatives && trace > best_trace)));
        if better {
            best_negatives = negatives;
            best_trace = trace;
            best_score = score;
            best_cells = det;
            best_cell = q;
            for r in 0..3 {
                for col in 0..3 {
                    best[r][col] = e[3 * r + col];
                }
            }
        }
    }
    if !(best_score < f64::INFINITY) {
        return Err(format!(
            "No supercell of {}{} cells within the search range",
            lowest,
            if highest > lowest { format!(" to {highest}") } else { String::new() }
        ));
    }
    let matrix = best.iter().flatten().map(|v| v.to_string()).collect::<Vec<_>>().join(" ");
    let supercell = find_modifier_type("supercell").ok_or("supercell step missing")?;
    let mut q = Json::object();
    q.set("matrix", matrix.as_str());
    (supercell.apply)(&q, data)?;
    let volume = determinant(&best_cell).abs();
    let mut height = f64::INFINITY;
    for k in 0..3 {
        height = fmin(height, volume / norm(cross(best_cell[(k + 1) % 3], best_cell[(k + 2) % 3])));
    }
    let label = match shape.as_str() {
        "image" => "shortest image distance ".to_string(),
        "sc" => "deviation from a cube ".to_string(),
        _ => format!("deviation from an {shape} cell "),
    };
    let note = format!(
        "Supercell matrix {} ({} cells, {} atoms); {}{}; shortest image {} Angstrom, inscribed sphere diameter {} Angstrom",
        matrix,
        best_cells,
        data.structure.atoms.len(),
        label,
        g5(if shape == "image" { shortest(&best_cell) } else { best_score }),
        g5(shortest(&best_cell)),
        g5(height)
    );
    data.notes.insert(0, note);
    Ok(())
}

// ---------------------------------------------------------------- interpolation
pub(crate) fn interpolate(p: &Json, data: &mut PipelineData) -> R {
    let id = "interpolate";
    let file = text(p, id, "file")?;
    if file.is_empty() {
        return Err("Choose the final structure file".into());
    }
    let t = number(p, id, "fraction")?;
    let method = text(p, id, "method")?;
    if method != "linear" && method != "idpp" {
        return Err("method is linear or idpp".into());
    }
    let frames = read_frames(Path::new(&file), None)?;
    let end = &frames.first().ok_or("The final structure file has no frames")?.structure;
    let s = &mut data.structure;
    if end.atoms.len() != s.atoms.len() {
        return Err(format!("The final structure has {} atoms, this one {}", end.atoms.len(), s.atoms.len()));
    }
    for i in 0..s.atoms.len() {
        if end.atoms[i].atomic_number != s.atoms[i].atomic_number {
            return Err(format!("Atom {} is {} here but {} in the final structure", i, s.atoms[i].symbol, end.atoms[i].symbol));
        }
    }
    let start = s.clone();
    for i in 0..s.atoms.len() {
        let mut d = sub(position_of(&end.atoms[i]), position_of(&s.atoms[i]));
        if s.has_unit_cell {
            d = minimum_image(s, d)?;
        }
        let q = add(position_of(&s.atoms[i]), scale(d, t));
        place(&mut s.atoms[i], q);
    }
    if s.has_unit_cell && end.has_unit_cell {
        for r in 0..3 {
            for c in 0..3 {
                s.cell_vectors[r][c] += t * (end.cell_vectors[r][c] - s.cell_vectors[r][c]);
            }
        }
    }
    let name = Path::new(&file).file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default();
    let mut note = format!("{} interpolation {} of the way to {}", if method == "idpp" { "IDPP" } else { "Linear" }, g5(t), name);
    if method == "idpp" {
        let n = s.atoms.len();
        if n > 3000 {
            return Err("IDPP: at most 3000 atoms".into());
        }
        let cutoff = number(p, id, "idpp_cutoff")?;
        let image = |at: &Structure, d: Vec3| -> Result<Vec3, String> {
            if at.has_unit_cell {
                minimum_image(at, d)
            } else {
                Ok(d)
            }
        };
        let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let d0 = norm(image(&start, sub(position_of(&start.atoms[j]), position_of(&start.atoms[i])))?);
                let d1 = norm(image(end, sub(position_of(&end.atoms[j]), position_of(&end.atoms[i])))?);
                if cutoff > 0.0 && fmin(d0, d1) > cutoff {
                    continue;
                }
                pairs.push((i, j, (1.0 - t) * d0 + t * d1));
            }
        }
        let mut x: Vec<Vec3> = s.atoms.iter().map(position_of).collect();
        let current: &Structure = s;
        let evaluate = |q: &[Vec3], force: &mut Vec<Vec3>| -> Result<f64, String> {
            let mut objective = 0.0;
            *force = vec![[0.0; 3]; n];
            for &(i, j, target) in &pairs {
                let v = image(current, sub(q[j], q[i]))?;
                let d = norm(v);
                if d < 1e-8 {
                    continue;
                }
                let w = 1.0 / (d * d * d * d);
                let diff = d - target;
                objective += w * diff * diff;
                let g = 2.0 * w * diff - 4.0 * w * diff * diff / d;
                let f = scale(v, -g / d);
                force[j] = add(force[j], f);
                force[i] = sub(force[i], f);
            }
            Ok(objective)
        };
        let mut velocity = vec![[0.0; 3]; n];
        let mut force = Vec::new();
        let (mut dt, mut alpha) = (0.1, 0.1);
        let before = evaluate(&x, &mut force)?;
        let mut now = before;
        let (mut positive, mut steps) = (0, 0);
        let max_steps = integer_of(p, id, "idpp_steps", 1)?;
        while steps < max_steps {
            let (mut power, mut fnorm, mut vnorm, mut fmx) = (0.0, 0.0, 0.0, 0.0);
            for i in 0..n {
                power += dot(force[i], velocity[i]);
                fnorm += dot(force[i], force[i]);
                vnorm += dot(velocity[i], velocity[i]);
                fmx = fmax(fmx, norm(force[i]));
            }
            if fmx < 1e-5 {
                break;
            }
            fnorm = fnorm.sqrt();
            vnorm = vnorm.sqrt();
            for i in 0..n {
                velocity[i] = add(scale(velocity[i], 1.0 - alpha), scale(force[i], if fnorm > 0.0 { alpha * vnorm / fnorm } else { 0.0 }));
            }
            if power > 0.0 {
                positive += 1;
                if positive > 5 {
                    dt = fmin(dt * 1.1, 1.0);
                    alpha *= 0.99;
                }
            } else {
                positive = 0;
                dt *= 0.5;
                alpha = 0.1;
                velocity.iter_mut().for_each(|v| *v = [0.0; 3]);
            }
            for i in 0..n {
                velocity[i] = add(velocity[i], scale(force[i], dt));
                let mut step = scale(velocity[i], dt);
                let len = norm(step);
                if len > 0.1 {
                    step = scale(step, 0.1 / len);
                }
                x[i] = add(x[i], step);
            }
            now = evaluate(&x, &mut force)?;
            steps += 1;
        }
        for i in 0..n {
            place(&mut s.atoms[i], x[i]);
        }
        note += &format!("; objective {} -> {} in {} steps ({} pairs)", g5(before), g5(now), steps, pairs.len());
    }
    data.notes.push(note);
    Ok(())
}

// ---------------------------------------------------------------- point defects
fn atoms_of(s: &Structure, symbol: &str) -> Vec<usize> {
    (0..s.atoms.len()).filter(|&i| s.atoms[i].symbol == symbol).collect()
}

fn recolour(a: &mut af_core::AtomSite) {
    [a.r, a.g, a.b] = af_core::elements::default_colour(a.atomic_number);
}

pub(crate) fn antisite(p: &Json, data: &mut PipelineData) -> R {
    let id = "create-antisite";
    let pair = element_list(&text(p, id, "elements")?, 2, "create-antisite")?;
    if pair[0] == pair[1] {
        return Err("create-antisite: give two different elements".into());
    }
    let s = &mut data.structure;
    let mut first = atoms_of(s, &pair[0]);
    let second = atoms_of(s, &pair[1]);
    if first.is_empty() || second.is_empty() {
        return Err(format!("create-antisite: the structure needs both {} and {}", pair[0], pair[1]));
    }
    let count = integer_of(p, id, "count", 1)?;
    let separation = number(p, id, "min_separation")?;
    let nearest_partner = text(p, id, "partner")? == "nearest";
    let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
    let distance = |s: &Structure, i: usize, j: usize| -> Result<f64, String> {
        let mut d = sub(position_of(&s.atoms[j]), position_of(&s.atoms[i]));
        if s.has_unit_cell {
            d = minimum_image(s, d)?;
        }
        Ok(norm(d))
    };
    shuffle(&mut first, &mut rng);
    let mut used: Vec<usize> = Vec::new();
    let far_from_used = |s: &Structure, used: &[usize], i: usize| -> Result<bool, String> {
        for &u in used {
            if i == u || distance(s, i, u)? < separation {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let mut swaps = String::new();
    let mut made = 0;
    for &a in &first {
        if made == count {
            break;
        }
        if !far_from_used(s, &used, a)? {
            continue;
        }
        let mut partners = Vec::new();
        for &j in &second {
            if far_from_used(s, &used, j)? {
                partners.push(j);
            }
        }
        if partners.is_empty() {
            continue;
        }
        let mut b = partners[uniform_u64(&mut rng, 0, partners.len() as u64 - 1) as usize];
        if nearest_partner {
            let mut best = f64::INFINITY;
            for &j in &partners {
                if distance(s, a, j)? < best {
                    best = distance(s, a, j)?;
                    b = j;
                }
            }
        }
        let (sa, za) = (s.atoms[a].symbol.clone(), s.atoms[a].atomic_number);
        s.atoms[a].symbol = s.atoms[b].symbol.clone();
        s.atoms[a].atomic_number = s.atoms[b].atomic_number;
        s.atoms[b].symbol = sa;
        s.atoms[b].atomic_number = za;
        recolour(&mut s.atoms[a]);
        recolour(&mut s.atoms[b]);
        used.push(a);
        used.push(b);
        swaps += &format!("{}{}<->{}", if made > 0 { ", " } else { "" }, a, b);
        made += 1;
    }
    if made == 0 {
        return Err("create-antisite: no pair satisfies the minimum separation".into());
    }
    data.selected = vec![0; data.structure.atoms.len()];
    for &u in &used {
        data.selected[u] = 1;
    }
    let mut note = format!("{} antisite pair(s) {}<->{}: atoms {}", made, pair[0], pair[1], swaps);
    if made < count {
        note += &format!(" (only {made} of {count} fit at the minimum separation)");
    }
    data.notes.push(note);
    Ok(())
}

pub(crate) fn frenkel_pair(p: &Json, data: &mut PipelineData) -> R {
    let id = "create-frenkel-pair";
    let symbol = standard_element(&text(p, id, "element")?)?;
    if atoms_of(&data.structure, &symbol).is_empty() {
        return Err(format!("create-frenkel-pair: no {symbol} atom"));
    }
    let kind = text(p, id, "interstitial")?;
    if kind != "void" && kind != "dumbbell" {
        return Err("interstitial is void or dumbbell".into());
    }
    let minimum = number(p, id, "min_distance")?;
    let maximum = number(p, id, "max_distance")?;
    let before = data.structure.atoms.len();
    let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
    let separation = |s: &Structure, i: usize, q: Vec3| -> Result<f64, String> {
        let mut d = sub(position_of(&s.atoms[i]), q);
        if s.has_unit_cell {
            d = minimum_image(s, d)?;
        }
        Ok(norm(d))
    };
    let mut host = before;
    let interstitial;
    if kind == "void" {
        let insert = find_modifier_type("insert-interstitials").ok_or("insert-interstitials step missing")?;
        let mut q = Json::object();
        q.set("element", symbol.as_str());
        q.set("count", 1.0);
        q.set("order", "random");
        q.set("seed", number(p, id, "seed")?);
        (insert.apply)(&q, data)?;
        if data.structure.atoms.len() != before + 1 {
            return Err("create-frenkel-pair: no interstitial void found".into());
        }
        interstitial = position_of(data.structure.atoms.last().unwrap());
    } else {
        let sites = atoms_of(&data.structure, &symbol);
        host = sites[uniform_u64(&mut rng, 0, sites.len() as u64 - 1) as usize];
        interstitial = position_of(&data.structure.atoms[host]);
    }
    let s = &data.structure;
    let mut candidates = Vec::new();
    let (mut farthest, mut farthest_distance) = (before, -1.0);
    for i in 0..before {
        if s.atoms[i].symbol != symbol || i == host {
            continue;
        }
        let dist = separation(s, i, interstitial)?;
        if dist >= minimum && (maximum <= 0.0 || dist <= maximum) {
            candidates.push(i);
        }
        if dist > farthest_distance && (maximum <= 0.0 || dist <= maximum) {
            farthest_distance = dist;
            farthest = i;
        }
    }
    if candidates.is_empty() && farthest == before {
        return Err(format!("create-frenkel-pair: no {symbol} atom within max_distance"));
    }
    let mut pick = Mt19937_64::new(cxx_u64(number(p, id, "seed")?).wrapping_add(7));
    let vacancy = if candidates.is_empty() { farthest } else { candidates[uniform_u64(&mut pick, 0, candidates.len() as u64 - 1) as usize] };
    let distance = separation(s, vacancy, interstitial)?;
    let note;
    if kind == "dumbbell" {
        let mut direction = vector(p, id, "direction")?;
        if flag(p, id, "lattice_direction")? {
            require_cell(data, "A lattice direction")?;
            direction = row_times(direction, &cell_of(&data.structure));
        }
        let length = norm(direction);
        if !(length > 1e-12) {
            return Err("The dumbbell direction must not be zero".into());
        }
        direction = scale(direction, 1.0 / length);
        let mut bond = number(p, id, "dumbbell_length")?;
        if !(bond > 0.0) {
            let mut nearest = f64::INFINITY;
            for i in 0..before {
                if i != host {
                    nearest = fmin(nearest, separation(&data.structure, i, interstitial)?);
                }
            }
            bond = 0.9 * nearest;
        }
        place(&mut data.structure.atoms[host], sub(interstitial, scale(direction, bond / 2.0)));
        if data.selected.len() != data.structure.atoms.len() {
            data.selected = vec![0; data.structure.atoms.len()];
        }
        append_atom(&mut data.structure, &mut data.selected, &symbol, add(interstitial, scale(direction, bond / 2.0)), false);
        note = format!(
            "Frenkel pair: {} dumbbell (length {} Angstrom) at atom {}, vacancy at atom {} {} Angstrom away",
            symbol,
            g5(bond),
            host,
            vacancy,
            g5(distance)
        );
    } else {
        note = format!("Frenkel pair: {} atom {} moved to an interstitial void {} Angstrom away", symbol, vacancy, g5(distance));
    }
    let keep: Vec<usize> = (0..data.structure.atoms.len()).filter(|&i| i != vacancy).collect();
    rebuild(data, &keep);
    data.selected = vec![0; data.structure.atoms.len()];
    *data.selected.last_mut().unwrap() = 1;
    if kind == "dumbbell" {
        data.selected[if host > vacancy { host - 1 } else { host }] = 1;
    }
    data.notes.push(note);
    Ok(())
}
