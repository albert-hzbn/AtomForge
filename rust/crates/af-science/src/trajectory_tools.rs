//! Trajectory tools (the C++ science/TrajectoryTools): averaged structures,
//! unwrapped (and smoothed) trajectories, per-frame time series of analyses
//! and differences between two frames.

use std::collections::BTreeMap;
use std::path::Path;

use af_core::json::Json;
use af_core::Structure;

use crate::analysis::{jv, R};
use crate::atom_properties::per_atom_properties;
use crate::core::*;
use crate::data::{FrameData, Parameters, TrajectoryStream};
use crate::tools::{per_atom_request, run_tool, ToolOutput};

struct Selection {
    first: usize,
    every: usize,
}

fn frame_selection(p: &Parameters, stream: &TrajectoryStream) -> Result<Selection, String> {
    let first = integer(p.number_or("first_frame", 0.0)?, "first_frame", 0)? as usize;
    let every = integer(p.number_or("every", 1.0)?, "every", 1)? as usize;
    if first >= stream.len() {
        return Err(format!("first_frame is beyond the trajectory ({} frames)", stream.len()));
    }
    Ok(Selection { first, every })
}

fn open_trajectory(p: &Parameters) -> Result<TrajectoryStream, String> {
    if !p.has("trajectory_file") {
        return Err("Choose a trajectory_file".into());
    }
    TrajectoryStream::new(Path::new(p.json("trajectory_file")?.string()?), None)
}

fn positions_of(s: &Structure) -> Vec<Vec3> {
    s.atoms.iter().map(|a| [a.x, a.y, a.z]).collect()
}

fn frame_mic(frame: &FrameData) -> Result<Option<MinimumImage>, String> {
    let periodic = frame.structure.has_unit_cell && frame.pbc.iter().any(|&b| b);
    Ok(if periodic { Some(MinimumImage::new(&frame.structure.cell_vectors, frame.pbc)?) } else { None })
}

/// Continues the unwrapped positions by the minimum-image step of each atom.
fn unwrap_step(unwrapped: &mut [Vec3], previous_wrapped: &[Vec3], current: &[Vec3], frame: &FrameData) -> Result<(), String> {
    let mic = frame_mic(frame)?;
    for i in 0..current.len() {
        let mut step = sub(current[i], previous_wrapped[i]);
        if let Some(m) = &mic {
            step = m.apply(step);
        }
        unwrapped[i] = add(unwrapped[i], step);
    }
    Ok(())
}

fn require_same_atoms(expected: usize, frame: &FrameData, index: usize) -> Result<(), String> {
    if frame.structure.atoms.len() != expected {
        return Err(format!(
            "Frame {} has {} atoms; the first has {} (the atoms must stay the same)",
            index,
            frame.structure.atoms.len(),
            expected
        ));
    }
    Ok(())
}

pub fn trajectory_average(p: &Parameters) -> Result<ToolOutput, String> {
    let stream = open_trajectory(p)?;
    let sel = frame_selection(p, &stream)?;
    let unwrap = p.boolean("unwrap", true)?;
    let first = stream.frame(sel.first)?;
    let n = first.structure.atoms.len();
    let mut wrapped = positions_of(&first.structure);
    let mut unwrapped = wrapped.clone();
    let mut sum = vec![[0.0; 3]; n];
    let mut sum_squares = vec![[0.0; 3]; n];
    let mut cell_sum = [[0.0; 3]; 3];
    let mut used = 0usize;
    let mut index = sel.first;
    while index < stream.len() {
        let frame = if index == sel.first { first.clone() } else { stream.frame(index)? };
        require_same_atoms(n, &frame, index)?;
        let current = positions_of(&frame.structure);
        if index != sel.first {
            if unwrap {
                unwrap_step(&mut unwrapped, &wrapped, &current, &frame)?;
            } else {
                unwrapped = current.clone();
            }
        }
        wrapped = current;
        for i in 0..n {
            for d in 0..3 {
                sum[i][d] += unwrapped[i][d];
                sum_squares[i][d] += unwrapped[i][d] * unwrapped[i][d];
            }
        }
        if frame.structure.has_unit_cell {
            for r in 0..3 {
                cell_sum[r] = add(cell_sum[r], frame.structure.cell_vectors[r]);
            }
        }
        used += 1;
        index += sel.every;
    }
    let mut average = first.structure.clone();
    let mut fluctuation = vec![0.0; n];
    let mut mean_fluctuation = 0.0;
    for i in 0..n {
        let mean = scale(sum[i], 1.0 / used as f64);
        let mut variance = 0.0;
        for d in 0..3 {
            let v = sum_squares[i][d] / used as f64 - mean[d] * mean[d];
            variance += if 0.0 < v { v } else { 0.0 };
        }
        fluctuation[i] = variance.sqrt();
        mean_fluctuation += fluctuation[i] / n as f64;
        average.atoms[i].x = mean[0];
        average.atoms[i].y = mean[1];
        average.atoms[i].z = mean[2];
    }
    if average.has_unit_cell {
        for r in 0..3 {
            average.cell_vectors[r] = scale(cell_sum[r], 1.0 / used as f64);
        }
    }
    let mut result = Json::object();
    result.set("frames_averaged", used as f64);
    result.set("rms_fluctuation_A", jv(&fluctuation));
    result.set("mean_rms_fluctuation_A", mean_fluctuation);
    result.set("unwrapped", unwrap);
    Ok(ToolOutput { result, frames: vec![average], ..Default::default() })
}

pub fn unwrap_trajectory(p: &Parameters) -> Result<ToolOutput, String> {
    let stream = open_trajectory(p)?;
    let sel = frame_selection(p, &stream)?;
    let window = integer(p.number_or("smooth_window", 1.0)?, "smooth_window", 1)?;
    let first = stream.frame(sel.first)?;
    let n = first.structure.atoms.len();
    let mut wrapped = positions_of(&first.structure);
    let mut unwrapped = wrapped.clone();
    let mut frames: Vec<Structure> = Vec::new();
    let mut times = Vec::new();
    let mut index = sel.first;
    while index < stream.len() {
        let mut frame = if index == sel.first { first.clone() } else { stream.frame(index)? };
        require_same_atoms(n, &frame, index)?;
        let current = positions_of(&frame.structure);
        if index != sel.first {
            unwrap_step(&mut unwrapped, &wrapped, &current, &frame)?;
        }
        wrapped = current;
        for i in 0..n {
            frame.structure.atoms[i].x = unwrapped[i][0];
            frame.structure.atoms[i].y = unwrapped[i][1];
            frame.structure.atoms[i].z = unwrapped[i][2];
        }
        times.push(frame.time_fs.unwrap_or(index as f64));
        frames.push(frame.structure);
        index += sel.every;
    }
    if window > 1 {
        let half = (window / 2) as usize;
        let mut smoothed = frames.clone();
        for f in 0..frames.len() {
            let low = if f >= half { f - half } else { 0 };
            let high = (frames.len() - 1).min(f + half);
            for i in 0..n {
                let mut mean = [0.0; 3];
                for g in low..=high {
                    mean = add(mean, [frames[g].atoms[i].x, frames[g].atoms[i].y, frames[g].atoms[i].z]);
                }
                mean = scale(mean, 1.0 / (high - low + 1) as f64);
                smoothed[f].atoms[i].x = mean[0];
                smoothed[f].atoms[i].y = mean[1];
                smoothed[f].atoms[i].z = mean[2];
            }
        }
        frames = smoothed;
    }
    let mut travelled = vec![0.0; n];
    if let (Some(first), Some(last)) = (frames.first(), frames.last()) {
        for i in 0..n {
            travelled[i] = norm(sub([last.atoms[i].x, last.atoms[i].y, last.atoms[i].z], [first.atoms[i].x, first.atoms[i].y, first.atoms[i].z]));
        }
    }
    let mut result = Json::object();
    result.set("frames", frames.len() as f64);
    result.set("smooth_window", window as f64);
    result.set("net_displacement_A", jv(&travelled));
    Ok(ToolOutput { result, frames, times, ..Default::default() })
}

/// Unwrapped paths of some atoms over every stride-th frame (for path plots).
pub fn unwrapped_paths(stream: &TrajectoryStream, atoms: &[usize], stride: usize, frames: &mut Vec<usize>) -> Result<Vec<Vec<Vec3>>, String> {
    let mut paths = vec![Vec::new(); atoms.len()];
    let (mut previous, mut unwrapped): (Vec<Vec3>, Vec<Vec3>) = (Vec::new(), Vec::new());
    frames.clear();
    let mut index = 0;
    while index < stream.len() {
        let frame = stream.frame(index)?;
        let mut current = Vec::new();
        for &atom in atoms {
            let Some(a) = frame.structure.atoms.get(atom) else {
                return Err("The number of atoms changes between frames".into());
            };
            current.push([a.x, a.y, a.z]);
        }
        if frames.is_empty() {
            unwrapped = current.clone();
        } else {
            unwrap_step(&mut unwrapped, &previous, &current, &frame)?;
        }
        previous = current;
        for k in 0..atoms.len() {
            paths[k].push(unwrapped[k]);
        }
        frames.push(index);
        index += stride.max(1);
    }
    Ok(paths)
}

fn optional_text(p: &Parameters, name: &str, fallback: &str) -> Result<String, String> {
    Ok(if p.has(name) { p.json(name)?.string()?.to_string() } else { fallback.to_string() })
}

pub fn time_series(p: &Parameters) -> R {
    let stream = open_trajectory(p)?;
    let sel = frame_selection(p, &stream)?;
    let tool = optional_text(p, "analysis", "cna")?;
    let wanted = optional_text(p, "property", "")?;
    let cutoff = p.number_or("cutoff_A", 3.0)?;
    let (mut frame_index, mut time, mut volume, mut density, mut mean_property) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut counts: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut count_order: Vec<String> = Vec::new();
    let mut property_name = String::new();
    let mut step = 0usize;
    let mut index = sel.first;
    while index < stream.len() {
        let frame = stream.frame(index)?;
        let s = &frame.structure;
        frame_index.push(index as f64);
        time.push(frame.time_fs.unwrap_or(f64::NAN));
        let v = if s.has_unit_cell { cell_volume(&s.cell_vectors) } else { f64::NAN };
        volume.push(v);
        density.push(if s.has_unit_cell { s.atoms.len() as f64 / v } else { f64::NAN });
        if tool != "none" {
            let result = run_tool(&tool, &per_atom_request(&tool, s, &frame.pbc, cutoff)?, Path::new("."), None)?.result;
            if let Some(Json::Object(members)) = result.find("counts") {
                for (name, value) in members {
                    if !counts.contains_key(name) {
                        counts.insert(name.clone(), vec![f64::NAN; step]);
                        count_order.push(name.clone());
                    }
                    counts.get_mut(name).unwrap().push(value.number()?);
                }
            }
            let mut mean = f64::NAN;
            for property in per_atom_properties(&tool, &result) {
                if !wanted.is_empty() && property.name != wanted {
                    continue;
                }
                property_name = property.name.clone();
                let (mut sum, mut valid) = (0.0, 0usize);
                for &x in &property.values {
                    if x.is_finite() {
                        sum += x;
                        valid += 1;
                    }
                }
                mean = if valid > 0 { sum / valid as f64 } else { f64::NAN };
                break;
            }
            mean_property.push(mean);
        }
        index += sel.every;
        step += 1;
    }
    let mut result = Json::object();
    result.set("frame", jv(&frame_index));
    result.set("time_fs", jv(&time));
    result.set("volume_A3", jv(&volume));
    result.set("number_density_per_A3", jv(&density));
    result.set("analysis", tool.as_str());
    if !mean_property.is_empty() {
        result.set("mean_property", jv(&mean_property));
        result.set("property", property_name.as_str());
    }
    let mut count_json = Json::object();
    for name in &count_order {
        count_json.set(name, jv(&counts[name]));
    }
    result.set("counts", count_json);
    Ok(result)
}

fn analyse_frame_property(frame: &FrameData, tool: &str, wanted: &str, cutoff: f64, name: &mut String) -> Result<Vec<f64>, String> {
    let result = run_tool(tool, &per_atom_request(tool, &frame.structure, &frame.pbc, cutoff)?, Path::new("."), None)?.result;
    for property in per_atom_properties(tool, &result) {
        if wanted.is_empty() || property.name == wanted {
            *name = property.name;
            return Ok(property.values);
        }
    }
    Err(format!("{tool}{}", if wanted.is_empty() { " gives no per-atom property".to_string() } else { format!(" has no property {wanted}") }))
}

pub fn frame_difference(p: &Parameters) -> R {
    let stream = open_trajectory(p)?;
    let count = stream.len() as i64;
    let mut reference = p.number_or("reference_frame", 0.0)? as i64;
    let mut current = p.number_or("frame", -1.0)? as i64;
    if reference < 0 {
        reference += count;
    }
    if current < 0 {
        current += count;
    }
    if reference < 0 || reference >= count || current < 0 || current >= count {
        return Err(format!("Frames must lie in the trajectory ({count} frames; negative counts from the end)"));
    }
    let a = stream.frame(reference as usize)?;
    let b = stream.frame(current as usize)?;
    if a.structure.atoms.len() != b.structure.atoms.len() {
        return Err("The two frames have different numbers of atoms".into());
    }
    let tool = optional_text(p, "analysis", "positions")?;
    let n = a.structure.atoms.len();
    let mut difference = vec![0.0; n];
    let (mut before, mut after) = (Vec::new(), Vec::new());
    let mut name = "displacement (A)".to_string();
    if tool == "positions" {
        let mic = frame_mic(&b)?;
        for i in 0..n {
            let (x, y) = (&b.structure.atoms[i], &a.structure.atoms[i]);
            let mut d = sub([x.x, x.y, x.z], [y.x, y.y, y.z]);
            if let Some(m) = &mic {
                d = m.apply(d);
            }
            difference[i] = norm(d);
        }
    } else {
        let wanted = optional_text(p, "property", "")?;
        let cutoff = p.number_or("cutoff_A", 3.0)?;
        before = analyse_frame_property(&a, &tool, &wanted, cutoff, &mut name)?;
        after = analyse_frame_property(&b, &tool, &wanted, cutoff, &mut name)?;
        for i in 0..n {
            difference[i] = after[i] - before[i];
        }
    }
    let (mut sum, mut squares, mut largest) = (0.0, 0.0, 0.0f64);
    let (mut changed, mut valid) = (0usize, 0usize);
    for &d in &difference {
        if !d.is_finite() {
            continue;
        }
        valid += 1;
        sum += d;
        squares += d * d;
        if largest < d.abs() {
            largest = d.abs();
        }
        if d.abs() > 1e-12 {
            changed += 1;
        }
    }
    let mut result = Json::object();
    result.set("difference", jv(&difference));
    if !before.is_empty() {
        result.set("reference_values", jv(&before));
        result.set("frame_values", jv(&after));
    }
    result.set("quantity", name.as_str());
    result.set("reference_frame", reference as f64);
    result.set("frame", current as f64);
    result.set("mean_difference", if valid > 0 { sum / valid as f64 } else { f64::NAN });
    result.set("rms_difference", if valid > 0 { (squares / valid as f64).sqrt() } else { f64::NAN });
    result.set("max_abs_difference", largest);
    result.set("atoms_changed", changed as f64);
    Ok(result)
}
