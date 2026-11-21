//! Dislocation lines from a reference and a deformed crystal (the C++
//! science/DislocationLines): core atoms whose Ackland-Jones structure type
//! changed, clustered into lines, with Burgers circuits about each line.

use std::collections::{BTreeMap, VecDeque};
use std::f64::consts::PI;

use af_core::cxx::sort_by;
use af_core::eigen::self_adjoint_eigen3;
use af_core::json::Json;
use af_core::structure::Structure;

use crate::analysis::{ackland_jones_type, ji, jv};
use crate::core::*;
use crate::data::{Parameters, StructureInput};
use crate::tools::ToolOutput;

fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

fn positions_of(s: &Structure) -> Vec<Vec3> {
    s.atoms.iter().map(|a| [a.x, a.y, a.z]).collect()
}

fn structure_types(positions: &[Vec3], cell: &Mat3, pbc: Pbc, cutoff: f64) -> Result<Vec<i32>, String> {
    let neighbors = neighbor_list(positions, cell, pbc, cutoff)?;
    let mut bonds: Vec<Vec<Vec3>> = vec![Vec::new(); positions.len()];
    for n in &neighbors {
        bonds[n.i].push(n.vector);
    }
    Ok(bonds.into_iter().map(ackland_jones_type).collect())
}

/// Removes the lattice translation closest to d (greedy descent over the
/// nearest-neighbour vectors, which generate fcc/bcc/hcp-basal lattices).
fn reduce_to_lattice(mut d: Vec3, lattice: &[Vec3]) -> Vec3 {
    for _ in 0..64 {
        let mut best: Option<Vec3> = None;
        let mut best_norm = norm(d) - 1e-9;
        for v in lattice {
            let candidate = norm(sub(d, *v));
            if candidate < best_norm {
                best_norm = candidate;
                best = Some(*v);
            }
        }
        let Some(v) = best else { break };
        d = sub(d, v);
    }
    d
}

fn min(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

struct Line {
    centroid: Vec3,
    direction: Vec3,
    length: f64,
    radius: f64,
    periodic_line: bool,
    source: &'static str,
    points: Vec<Vec3>,
    atoms: usize,
}

pub fn dislocation_lines(reference: &StructureInput, deformed: &StructureInput, p: &Parameters) -> Result<ToolOutput, String> {
    let r = &reference.structure;
    let d = &deformed.structure;
    if r.atoms.len() != d.atoms.len() || r.atoms.is_empty() {
        return Err("Reference and deformed structures need the same atoms in the same order".into());
    }
    let periodic = r.has_unit_cell && reference.pbc[0] && reference.pbc[1] && reference.pbc[2];
    let cell = if r.has_unit_cell { r.cell_vectors } else { identity() };
    let pbc = [periodic; 3];
    let mic = if periodic { MinimumImage::new(&cell, pbc)? } else { MinimumImage::default() };
    let ref_positions = positions_of(r);
    let positions = positions_of(d);
    let n = positions.len();

    // Nearest-neighbour spacing and lattice vectors from the reference crystal.
    let mut nearest = vec![f64::INFINITY; n];
    let probe = neighbor_list(&ref_positions, &cell, pbc, 6.0)?;
    for nb in &probe {
        nearest[nb.i] = min(nearest[nb.i], norm(nb.vector));
    }
    let mut sorted: Vec<f64> = nearest.iter().copied().filter(|d| d.is_finite()).collect();
    if sorted.is_empty() {
        return Err("No neighbours within 6 Angstrom".into());
    }
    sort_by(&mut sorted, |a, b| a < b);
    let spacing = sorted[sorted.len() / 2];
    let mut cutoff = p.number_or("cutoff_A", 0.0)?;
    if cutoff <= 0.0 {
        cutoff = 1.6 * spacing; // beyond the second bcc shell (1.155 d)
    }
    let reference_types = structure_types(&ref_positions, &cell, pbc, cutoff)?;
    let deformed_types = structure_types(&positions, &if d.has_unit_cell { cell } else { identity() }, pbc, cutoff)?;
    let mut lattice: Vec<Vec3> = Vec::new();
    for i in 0..n {
        if !lattice.is_empty() {
            break;
        }
        if reference_types[i] == 0 {
            continue;
        }
        for nb in &probe {
            if nb.i == i && norm(nb.vector) < 1.2 * spacing {
                lattice.push(nb.vector);
            }
        }
    }
    if lattice.is_empty() {
        return Err("The reference has no crystalline (fcc, hcp, bcc) atoms".into());
    }

    // Defect atoms: crystalline in the reference, changed structure type now.
    let defect: Vec<usize> = (0..n).filter(|&i| reference_types[i] != 0 && deformed_types[i] != reference_types[i]).collect();
    if defect.is_empty() {
        return Err("No atoms changed crystal structure: no dislocation cores found".into());
    }
    let cluster_cutoff = if p.has("cluster_cutoff_A") && p.number("cluster_cutoff_A")? > 0.0 { p.number("cluster_cutoff_A")? } else { 1.5 * spacing };
    let minimum_atoms = integer(p.number_or("min_atoms", 4.0)?, "min_atoms", 1)?;
    let defect_positions: Vec<Vec3> = defect.iter().map(|&i| ref_positions[i]).collect();
    let links = neighbor_list(&defect_positions, &cell, pbc, cluster_cutoff)?;
    let mut parent: Vec<usize> = (0..defect.len()).collect();
    for link in &links {
        // C++17 evaluates the right-hand side of the assignment first.
        let rj = find(&mut parent, link.j);
        let ri = find(&mut parent, link.i);
        parent[ri] = rj;
    }
    let mut clusters: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for k in 0..defect.len() {
        let root = find(&mut parent, k);
        clusters.entry(root).or_default().push(k);
    }
    let mut adjacency: Vec<Vec<(usize, Vec3)>> = vec![Vec::new(); defect.len()];
    for link in &links {
        adjacency[link.i].push((link.j, link.vector));
    }
    let mut ordered: Vec<Vec<usize>> = clusters.into_values().filter(|m| m.len() as i64 >= minimum_atoms).collect();
    sort_by(&mut ordered, |a, b| a.len() > b.len());

    let mut found: Vec<Line> = Vec::new();
    let mut line_id = vec![-1i32; n];
    for (c, members) in ordered.iter().enumerate() {
        let mut unwrapped: BTreeMap<usize, Vec3> = BTreeMap::new();
        unwrapped.insert(members[0], defect_positions[members[0]]);
        let mut queue = VecDeque::from([members[0]]);
        let mut period = [0.0; 3];
        let mut period_length = f64::INFINITY;
        while let Some(k) = queue.pop_front() {
            for &(j, vector) in &adjacency[k] {
                let reached = add(unwrapped[&k], vector);
                let Some(&existing) = unwrapped.get(&j) else {
                    unwrapped.insert(j, reached);
                    queue.push_back(j);
                    continue;
                };
                let mismatch = sub(reached, existing);
                if norm(mismatch) > 0.5 * spacing && norm(mismatch) < period_length {
                    period = mismatch;
                    period_length = norm(mismatch);
                }
            }
        }
        let mut mean = [0.0; 3];
        for &k in members {
            mean = add(mean, unwrapped[&k]);
            line_id[defect[k]] = c as i32;
        }
        mean = scale(mean, 1.0 / members.len() as f64);
        let mut covariance = [[0.0; 3]; 3];
        for &k in members {
            let d = sub(unwrapped[&k], mean);
            for i in 0..3 {
                for j in 0..3 {
                    covariance[i][j] += d[i] * d[j];
                }
            }
        }
        let (_, vectors) = self_adjoint_eigen3(&covariance).ok_or_else(|| "Eigen decomposition did not converge".to_string())?;
        let mut t = [vectors[0][2], vectors[1][2], vectors[2][2]];
        let mut source = "shape";
        let periodic_line = period_length.is_finite();
        if periodic_line {
            // Eigen's normalized(): divide by the norm.
            let length = (period[0] * period[0] + period[1] * period[1] + period[2] * period[2]).sqrt();
            t = [period[0] / length, period[1] / length, period[2] / length];
            source = "periodic";
        }
        let mut largest = 0;
        for k in 1..3 {
            if t[k].abs() > t[largest].abs() {
                largest = k;
            }
        }
        if t[largest] < 0.0 {
            t = [-t[0], -t[1], -t[2]];
        }
        let direction = t;
        let (mut low, mut high, mut radius) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
        for &k in members {
            let d = sub(unwrapped[&k], mean);
            let s = dot(d, direction);
            low = min(low, s);
            high = max(high, s);
            radius = max(radius, norm(sub(d, scale(direction, s))));
        }
        let length = if periodic_line { period_length } else { high - low + spacing };
        let bin = max(1.5 * spacing, (high - low) / 20.0);
        let mut bins: BTreeMap<i64, (Vec3, i32)> = BTreeMap::new();
        for &k in members {
            let key = ((dot(sub(unwrapped[&k], mean), direction) - low) / bin).floor() as i64;
            let entry = bins.entry(key).or_insert(([0.0; 3], 0));
            entry.0 = add(entry.0, unwrapped[&k]);
            entry.1 += 1;
        }
        let points = bins.values().map(|(sum, count)| scale(*sum, 1.0 / *count as f64)).collect();
        found.push(Line { centroid: mean, direction, length, radius, periodic_line, source, points, atoms: members.len() });
    }

    // Burgers circuits: right-handed loops about each line, through atoms
    // near a circle in the plane normal to the line, in the reference frame.
    // Displacement steps are reduced modulo lattice vectors, so the slip
    // jump across the cut drops out and the closure failure is b.
    let displacement: Vec<Vec3> = (0..n).map(|i| sub(positions[i], ref_positions[i])).collect();
    let requested_radius = p.number_or("circuit_radius_A", 0.0)?;
    let mut lines = Vec::new();
    let mut other_clusters = 0usize;
    for (c, line) in found.iter().enumerate() {
        let mut radius = if requested_radius > 0.0 { requested_radius } else { line.radius + 2.0 * spacing };
        // Keep the circuit clear of neighbouring lines.
        for (o, other) in found.iter().enumerate() {
            if o == c {
                continue;
            }
            let mut d = mic.apply(sub(other.centroid, line.centroid));
            d = sub(d, scale(line.direction, dot(d, line.direction)));
            radius = min(radius, max(spacing, 0.5 * norm(d)));
        }
        let helper = if line.direction[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let mut e1 = cross(line.direction, helper);
        e1 = scale(e1, 1.0 / norm(e1));
        let e2 = cross(line.direction, e1);
        let mut circuit: Vec<usize> = Vec::new();
        let samples = 16.max((2.0 * PI * radius / (0.7 * spacing)) as i32);
        for k in 0..samples {
            let angle = 2.0 * PI * k as f64 / samples as f64;
            let target = add(line.centroid, add(scale(e1, radius * angle.cos()), scale(e2, radius * angle.sin())));
            let mut best = n;
            let mut best_distance = f64::INFINITY;
            for i in 0..n {
                let d = mic.apply(sub(ref_positions[i], target));
                if dot(d, line.direction).abs() > 0.75 * spacing {
                    continue;
                }
                let distance = norm(d);
                if distance < best_distance {
                    best_distance = distance;
                    best = i;
                }
            }
            if best < n && best_distance < spacing && circuit.last() != Some(&best) {
                circuit.push(best);
            }
        }
        if circuit.len() > 1 && circuit.first() == circuit.last() {
            circuit.pop();
        }
        let mut burgers = [0.0; 3];
        let closed = circuit.len() >= 6;
        if closed {
            for k in 0..circuit.len() {
                let (a, b) = (circuit[k], circuit[(k + 1) % circuit.len()]);
                burgers = add(burgers, reduce_to_lattice(sub(displacement[b], displacement[a]), &lattice));
            }
        }
        let magnitude = norm(burgers);
        // A closed circuit without closure failure encloses no dislocation.
        if closed && magnitude < 0.1 * spacing {
            other_clusters += 1;
            for id in line_id.iter_mut() {
                if *id == c as i32 {
                    *id = -2;
                }
            }
            continue;
        }
        let mut entry = Json::object();
        entry.set("atoms", line.atoms);
        entry.set("direction", jv(&line.direction));
        entry.set("direction_source", line.source);
        entry.set("burgers_vector_A", jv(&burgers));
        entry.set("burgers_magnitude_A", magnitude);
        entry.set(
            "character_angle_deg",
            if magnitude > 1e-6 { Json::from(min(1.0, dot(burgers, line.direction).abs() / magnitude).acos() * 180.0 / PI) } else { Json::Null },
        );
        entry.set("centroid_A", jv(&line.centroid));
        entry.set("length_A", line.length);
        entry.set("periodic_line", line.periodic_line);
        entry.set("core_radius_A", line.radius);
        entry.set("circuit_radius_A", radius);
        entry.set("circuit_atoms", circuit.len());
        entry.set("circuit_closed", closed);
        entry.set("line_points_A", Json::Array(line.points.iter().map(|p| jv(p)).collect()));
        lines.push(entry);
    }
    let mut result = Json::object();
    let count = lines.len();
    result.set("lines", Json::Array(lines));
    result.set("line_count", count);
    result.set("non_dislocation_clusters", other_clusters);
    result.set("defect_atoms", defect.len());
    result.set("nearest_neighbour_A", spacing);
    result.set("cluster_cutoff_A", cluster_cutoff);
    let mut renumber: BTreeMap<i32, i32> = BTreeMap::new();
    let mut kept = 0;
    for c in 0..found.len() as i32 {
        if line_id.contains(&c) {
            renumber.insert(c, kept);
            kept += 1;
        }
    }
    for id in line_id.iter_mut() {
        *id = if *id >= 0 { renumber[id] } else { -1 };
    }
    result.set("structure_type", ji(&deformed_types));
    result.set("line_id", ji(&line_id));
    Ok(ToolOutput { result, ..Default::default() })
}
