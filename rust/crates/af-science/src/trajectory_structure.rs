//! Trajectory-averaged structure (the C++ science/TrajectoryStructure):
//! total and partial g(r), coordination numbers and their distribution, and
//! the bond-angle distribution.

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::path::Path;

use af_core::json::Json;

use crate::analysis::jv;
use crate::core::*;
use crate::data::{Parameters, TrajectoryStream};

pub fn trajectory_structure(p: &Parameters) -> Result<Json, String> {
    if !p.has("trajectory_file") {
        return Err("trajectory-structure needs trajectory_file".into());
    }
    let stream = TrajectoryStream::new(Path::new(p.json("trajectory_file")?.string()?), None)?;
    let first = integer(p.number_or("first_frame", 0.0)?, "first_frame", 0)?;
    let every = integer(p.number_or("every", 1.0)?, "every", 1)?;
    let r_max = positive(p.number_or("r_max_A", 8.0)?, "r_max_A")?;
    let bins = integer(p.number_or("bins", 200.0)?, "bins", 10)? as usize;
    let cutoff = positive(p.number_or("cutoff_A", 3.0)?, "cutoff_A")?;
    let angle_bins = integer(p.number_or("angle_bins", 180.0)?, "angle_bins", 10)? as usize;
    if first as usize >= stream.len() {
        return Err("first_frame is beyond the trajectory".into());
    }
    let dr = r_max / bins as f64;

    let mut elements: Vec<String> = Vec::new();
    // pair_counts[a][b][bin]: b atoms in each shell around a atoms, summed over frames.
    let mut pair_counts: BTreeMap<String, BTreeMap<String, Vec<f64>>> = BTreeMap::new();
    // Sum over frames of N_a (N_b - delta_ab) / V.
    let mut pair_normalization: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    let mut coordination_sum: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    let mut centre_count: BTreeMap<String, f64> = BTreeMap::new();
    let mut coordination_histogram: BTreeMap<String, BTreeMap<i64, f64>> = BTreeMap::new();
    let mut angles = vec![0.0; angle_bins];
    let mut angle_total = 0.0;
    let mut used: i64 = 0;
    let mut index = first as usize;
    while index < stream.len() {
        let frame = stream.frame(index)?;
        let s = &frame.structure;
        if !s.has_unit_cell || !(frame.pbc[0] && frame.pbc[1] && frame.pbc[2]) {
            return Err(format!("Trajectory-averaged g(r) needs periodic frames (frame {index})"));
        }
        let cell = s.cell_vectors;
        let volume = cell_volume(&cell);
        let mut positions = Vec::new();
        let mut count: BTreeMap<String, f64> = BTreeMap::new();
        for atom in &s.atoms {
            positions.push([atom.x, atom.y, atom.z]);
            *count.entry(atom.symbol.clone()).or_insert(0.0) += 1.0;
            if !elements.contains(&atom.symbol) {
                elements.push(atom.symbol.clone());
            }
        }
        for (a, &na) in &count {
            for (b, &nb) in &count {
                *pair_normalization.entry(a.clone()).or_default().entry(b.clone()).or_insert(0.0) += na * (nb - if a == b { 1.0 } else { 0.0 }) / volume;
            }
        }
        for (a, &na) in &count {
            *centre_count.entry(a.clone()).or_insert(0.0) += na;
        }
        // Neighbour list out to the larger of r_max and the coordination cutoff.
        let mut bonds: Vec<Vec<Vec3>> = vec![Vec::new(); s.atoms.len()];
        let mut coordination: Vec<BTreeMap<String, i64>> = vec![BTreeMap::new(); s.atoms.len()];
        let reach = if r_max < cutoff { cutoff } else { r_max };
        for n in neighbor_list(&positions, &cell, frame.pbc, reach)? {
            let r = norm(n.vector);
            let a = &s.atoms[n.i].symbol;
            let b = &s.atoms[n.j].symbol;
            if r < r_max {
                let histogram = pair_counts.entry(a.clone()).or_default().entry(b.clone()).or_insert_with(|| vec![0.0; bins]);
                histogram[((r / dr) as usize).min(bins - 1)] += 1.0;
            }
            if r < cutoff {
                bonds[n.i].push(n.vector);
                *coordination[n.i].entry(b.clone()).or_insert(0) += 1;
            }
        }
        for i in 0..s.atoms.len() {
            let a = &s.atoms[i].symbol;
            let mut total: i64 = 0;
            for (b, &c) in &coordination[i] {
                *coordination_sum.entry(a.clone()).or_default().entry(b.clone()).or_insert(0.0) += c as f64;
                total += c;
            }
            *coordination_histogram.entry(a.clone()).or_default().entry(total).or_insert(0.0) += 1.0;
            let v = &bonds[i];
            for j in 0..v.len() {
                for k in j + 1..v.len() {
                    let cosine = (dot(v[j], v[k]) / (norm(v[j]) * norm(v[k]))).clamp(-1.0, 1.0);
                    let degrees = cosine.acos() * 180.0 / PI;
                    angles[((degrees / 180.0 * angle_bins as f64) as usize).min(angle_bins - 1)] += 1.0;
                    angle_total += 1.0;
                }
            }
        }
        used += 1;
        index += every as usize;
    }
    elements.sort();
    let mut r = vec![0.0; bins];
    let mut shell = vec![0.0; bins];
    for k in 0..bins {
        let lo = k as f64 * dr;
        let hi = lo + dr;
        r[k] = lo + 0.5 * dr;
        shell[k] = 4.0 / 3.0 * PI * (hi * hi * hi - lo * lo * lo);
    }
    // g_ab(r) = <counts of b around a> / (N_a rho_b shell volume), averaged
    // over frames; the total g(r) weights pairs by their number.
    let lookup = |m: &BTreeMap<String, BTreeMap<String, f64>>, a: &str, b: &str| m.get(a).and_then(|x| x.get(b)).copied().unwrap_or(0.0);
    let mut partial = Json::object();
    let mut total_counts = vec![0.0; bins];
    let mut total_normalization = 0.0;
    for a in &elements {
        for b in &elements {
            let normalization = lookup(&pair_normalization, a, b);
            let found = pair_counts.get(a).and_then(|x| x.get(b));
            let mut g = vec![0.0; bins];
            for k in 0..bins {
                let counts = found.map_or(0.0, |h| h[k]);
                if normalization > 0.0 {
                    g[k] = counts / (normalization * shell[k]);
                }
                total_counts[k] += counts;
            }
            total_normalization += normalization;
            partial.set(&format!("{a}-{b}"), jv(&g));
        }
    }
    let mut total = vec![0.0; bins];
    for k in 0..bins {
        if total_normalization > 0.0 {
            total[k] = total_counts[k] / (total_normalization * shell[k]);
        }
    }
    let (mut coordination_json, mut pair_coordination, mut distribution) = (Json::object(), Json::object(), Json::object());
    let (mut all_neighbours, mut all_centres) = (0.0, 0.0);
    for a in &elements {
        let centres = centre_count.get(a).copied().unwrap_or(0.0);
        let mut sum = 0.0;
        for b in &elements {
            let value = lookup(&coordination_sum, a, b);
            pair_coordination.set(&format!("{b} around {a}"), value / centres);
            sum += value;
        }
        coordination_json.set(a, sum / centres);
        all_neighbours += sum;
        all_centres += centres;
        let rows: Vec<Json> = coordination_histogram
            .get(a)
            .map(|h| h.iter().map(|(&n, &c)| Json::Array(vec![Json::from(n as f64), Json::from(c / centres)])).collect())
            .unwrap_or_default();
        distribution.set(a, Json::Array(rows));
    }
    let width = 180.0 / angle_bins as f64;
    let angle_grid: Vec<f64> = (0..angle_bins).map(|k| (k as f64 + 0.5) * width).collect();
    let angle_density: Vec<f64> = (0..angle_bins).map(|k| if angle_total > 0.0 { angles[k] / (angle_total * width) } else { 0.0 }).collect();
    let mut result = Json::object();
    result.set("frames_used", used);
    result.set("frame_count", stream.len() as i64);
    result.set("elements", Json::Array(elements.iter().map(|e| Json::from(e.as_str())).collect()));
    result.set("r_A", jv(&r));
    result.set("g_total", jv(&total));
    result.set("g_partial", partial);
    result.set("cutoff_A", cutoff);
    result.set("mean_coordination", if all_centres > 0.0 { all_neighbours / all_centres } else { 0.0 });
    result.set("coordination", coordination_json);
    result.set("pair_coordination", pair_coordination);
    result.set("coordination_distribution", distribution);
    result.set("angle_deg", jv(&angle_grid));
    result.set("angle_density_per_deg", jv(&angle_density));
    result.set("bond_angles_counted", angle_total);
    Ok(result)
}
