//! Measurements on the selected atoms (the C++ ui/MeasurementOverlay
//! processMeasurementRequests): distance, angle, dihedral and the Atom Info
//! text, over the scene's instances (periodic images included) in the order
//! the atoms were picked.

use af_core::cell::{BOND_TOLERANCE_FACTOR, MIN_BOND_DISTANCE};
use af_core::elements;
use af_core::glm::{self, length};
use af_core::Structure;
use glam::{DVec3, Mat3, Vec3};

use crate::scene::Scene;

/// A measurement: the value ("2.5456 Å"; empty when it failed), the atoms
/// measured or why it failed, and the instances to draw lines between.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Measurement {
    pub value: String,
    pub message: String,
    pub lines: Vec<usize>,
}

/// A selection made by atom covers every periodic image. When it holds
/// exactly `count` distinct atoms, keep one image of each, as close together
/// as possible (nearestImages); any other selection is returned unchanged.
pub fn nearest_images(selected: &[usize], scene: &Scene, count: usize) -> Vec<usize> {
    if selected.len() <= count {
        return selected.to_vec();
    }
    let positions = &scene.instances.positions;
    let indices = &scene.instances.atom_indices;
    let mut atoms: Vec<i32> = Vec::new();
    let mut images: Vec<Vec<usize>> = Vec::new();
    for &instance in selected {
        if instance >= indices.len() || instance >= positions.len() {
            return selected.to_vec();
        }
        let atom = indices[instance];
        match atoms.iter().position(|&a| a == atom) {
            Some(k) => images[k].push(instance),
            None => {
                atoms.push(atom);
                images.push(vec![instance]);
            }
        }
    }
    if atoms.len() != count {
        return selected.to_vec();
    }
    if count == 2 {
        let mut best = vec![images[0][0], images[1][0]];
        let mut best_distance = f32::MAX;
        for &a in &images[0] {
            for &b in &images[1] {
                let d = length(positions[a] - positions[b]);
                if d < best_distance {
                    best_distance = d;
                    best = vec![a, b];
                }
            }
        }
        return best;
    }
    if count == 4 {
        // A dihedral: a chain, each atom's image nearest the one before it.
        let mut chain = vec![images[0][0]];
        for k in 1..count {
            let mut best = images[k][0];
            let mut best_distance = f32::MAX;
            for &i in &images[k] {
                let d = length(positions[i] - positions[*chain.last().unwrap()]);
                if d < best_distance {
                    best_distance = d;
                    best = i;
                }
            }
            chain.push(best);
        }
        return chain;
    }
    let anchor = if count == 3 { 1 } else { 0 };
    let centre = positions[images[anchor][0]];
    let mut result = Vec::new();
    for (k, list) in images.iter().enumerate() {
        let mut best = list[0];
        let mut best_distance = f32::MAX;
        if k != anchor {
            for &i in list {
                let d = length(positions[i] - centre);
                if d < best_distance {
                    best_distance = d;
                    best = i;
                }
            }
        }
        result.push(best);
    }
    result
}

fn atom_of(scene: &Scene, s: &Structure, instance: usize, fallback: &'static str) -> (i32, String) {
    let base = scene.instances.atom_indices.get(instance).copied().unwrap_or(-1);
    let symbol = if base >= 0 && (base as usize) < s.atoms.len() { s.atoms[base as usize].symbol.clone() } else { fallback.to_string() };
    (base, symbol)
}

pub fn distance(selected: &[usize], scene: &Scene, s: &Structure) -> Measurement {
    let pair = nearest_images(selected, scene, 2);
    let positions = &scene.instances.positions;
    if pair.len() != 2 {
        return Measurement { message: "Select exactly 2 atoms to measure distance.".into(), ..Default::default() };
    }
    if pair.iter().any(|&i| i >= positions.len()) {
        return Measurement { message: "Unable to measure distance for current selection.".into(), ..Default::default() };
    }
    let d = length(positions[pair[0]] - positions[pair[1]]);
    let (a, sa) = atom_of(scene, s, pair[0], "A");
    let (b, sb) = atom_of(scene, s, pair[1], "B");
    Measurement { value: format!("{d:.4} Å"), message: format!("{sa} #{a}  to  {sb} #{b}"), lines: pair }
}

pub fn angle(selected: &[usize], scene: &Scene, s: &Structure) -> Measurement {
    let triple = nearest_images(selected, scene, 3);
    let positions = &scene.instances.positions;
    if triple.len() != 3 {
        return Measurement { message: "Select exactly 3 atoms to measure angle.".into(), ..Default::default() };
    }
    if triple.iter().any(|&i| i >= positions.len()) {
        return Measurement { message: "Unable to measure angle for current selection.".into(), ..Default::default() };
    }
    let (p0, p1, p2) = (positions[triple[0]], positions[triple[1]], positions[triple[2]]);
    let v0 = glm::normalize(p0 - p1);
    let v2 = glm::normalize(p2 - p1);
    let cos = glm::dot(v0, v2).clamp(-1.0, 1.0);
    // glm::degrees in float.
    let degrees = cos.acos() * 57.295_779_513_082_32_f32;
    let (b0, s0) = atom_of(scene, s, triple[0], "A");
    let (b1, s1) = atom_of(scene, s, triple[1], "B");
    let (b2, s2) = atom_of(scene, s, triple[2], "C");
    Measurement { value: format!("{degrees:.2}°"), message: format!("{s0} #{b0}, {s1} #{b1} (vertex), {s2} #{b2}"), lines: triple }
}

pub fn dihedral(selected: &[usize], scene: &Scene, s: &Structure) -> Measurement {
    let four = nearest_images(selected, scene, 4);
    let positions = &scene.instances.positions;
    if four.len() != 4 || four.iter().any(|&i| i >= positions.len()) {
        return Measurement { message: "Select exactly 4 atoms (in order along the chain) to measure a dihedral.".into(), ..Default::default() };
    }
    let p: Vec<DVec3> = four.iter().map(|&i| positions[i].as_dvec3()).collect();
    let (b1, b2, b3) = (p[1] - p[0], p[2] - p[1], p[3] - p[2]);
    let (n1, n2) = (b1.cross(b2), b2.cross(b3));
    if n1.length() < 1e-9 || n2.length() < 1e-9 {
        return Measurement { message: "Three of the atoms are on a line: the dihedral is undefined.".into(), ..Default::default() };
    }
    let m = n1.cross(b2 / b2.length());
    let degrees = m.dot(n2).atan2(n1.dot(n2)).to_degrees();
    let atoms: Vec<String> = four
        .iter()
        .map(|&i| {
            let (base, symbol) = atom_of(scene, s, i, "?");
            format!("{symbol} #{base}")
        })
        .collect();
    Measurement { value: format!("{degrees:.2}°"), message: atoms.join(", "), lines: Vec::new() }
}

/// Coordination and bond lengths of an atom (computeBondStatsForAtom).
fn bond_stats(instance: usize, base: i32, scene: &Scene, s: &Structure, radii: &[f32]) -> (i32, f32, f32, f32) {
    let positions = &scene.instances.positions;
    let scales = &scene.instances.scales;
    if instance >= positions.len() || instance >= scales.len() {
        return (0, 0.0, 0.0, 0.0);
    }
    let use_base = base >= 0 && (base as usize) < s.atoms.len();
    let centre = if use_base {
        let a = &s.atoms[base as usize];
        Vec3::new(a.x as f32, a.y as f32, a.z as f32)
    } else {
        positions[instance]
    };
    let radius = |b: i32, i: i32| -> f32 {
        if b >= 0 && (b as usize) < s.atoms.len() {
            let z = s.atoms[b as usize].atomic_number;
            if z >= 0 && (z as usize) < radii.len() && radii[z as usize] > 0.0 {
                return radii[z as usize];
            }
        }
        if i >= 0 && (i as usize) < scales.len() {
            return scales[i as usize];
        }
        1.0
    };
    let r0 = radius(base, instance as i32);
    let matrices = glm::cell_matrices(s);
    let use_pbc = matrices.is_some();
    let (cell, inv) = matrices.unwrap_or((Mat3::IDENTITY, Mat3::IDENTITY));
    let (mut min, mut max, mut sum, mut count) = (f32::MAX, 0.0f32, 0.0f32, 0);
    let mut consider = |d: f32, r1: f32| {
        if d <= MIN_BOND_DISTANCE || d > (r0 + r1) * BOND_TOLERANCE_FACTOR {
            return;
        }
        min = min.min(d);
        max = max.max(d);
        sum += d;
        count += 1;
    };
    if use_base {
        for (j, a) in s.atoms.iter().enumerate() {
            if j as i32 == base {
                continue;
            }
            let delta = glm::minimum_image_delta(Vec3::new(a.x as f32, a.y as f32, a.z as f32) - centre, use_pbc, &cell, &inv);
            consider(length(delta), radius(j as i32, -1));
        }
    } else {
        for i in 0..positions.len() {
            if i == instance || i >= scales.len() {
                continue;
            }
            let other = scene.instances.atom_indices.get(i).copied().unwrap_or(-1);
            consider(length(positions[i] - centre), radius(other, i as i32));
        }
    }
    if count == 0 {
        return (0, 0.0, 0.0, 0.0);
    }
    (count, sum / count as f32, min, max)
}

/// %.4g / %.5g as printf writes them.
fn g(v: f64, precision: usize) -> String {
    af_core::format::format_g(v, precision)
}

/// The Atom Info text of one selected atom.
pub fn atom_info(selected: &[usize], scene: &Scene, s: &Structure, radii: &[f32]) -> String {
    let single = nearest_images(selected, scene, 1);
    if single.len() != 1 {
        return "Select exactly 1 atom to view info.".into();
    }
    let idx = single[0];
    if idx >= scene.instances.positions.len() || idx >= scene.instances.atom_indices.len() {
        return "Unable to retrieve atom info.".into();
    }
    let base = scene.instances.atom_indices[idx];
    let pos = scene.instances.positions[idx];
    let mut text = String::new();
    if base >= 0 && (base as usize) < s.atoms.len() {
        let a = &s.atoms[base as usize];
        text += &format!("Element:  {} ({})\n", elements::name(a.atomic_number), a.symbol);
        text += &format!("Atomic number:  {}\n", a.atomic_number);
        text += &format!("Index:  #{base}\n");
        if a.fixed != 0 {
            text += &format!(
                "Fixed:  {}{}{} (selective dynamics)\n",
                if a.fixed & 1 != 0 { "x" } else { "" },
                if a.fixed & 2 != 0 { "y" } else { "" },
                if a.fixed & 4 != 0 { "z" } else { "" }
            );
        }
        if a.magmom != 0.0 {
            text += &format!("Magnetic moment:  {} μB\n", g(a.magmom as f64, 4));
        }
        if a.charge != 0.0 {
            text += &format!("Charge:  {} e\n", g(a.charge as f64, 4));
        }
        if a.velocity.iter().any(|&v| v != 0.0) {
            text += &format!("Velocity:  ({}, {}, {}) Å/fs\n", g(a.velocity[0] as f64, 5), g(a.velocity[1] as f64, 5), g(a.velocity[2] as f64, 5));
        }
    }
    let shown = |v: f32| if v.abs() < 5e-6 { 0.0 } else { v };
    text += &format!("Cartesian:  ({:.5}, {:.5}, {:.5}) Å\n", shown(pos.x), shown(pos.y), shown(pos.z));
    if s.has_unit_cell {
        if let Some((_, inv)) = glm::cell_matrices(s) {
            let origin = Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32);
            let f = glm::mul(&inv, pos - origin);
            text += &format!("Direct:  ({:.5}, {:.5}, {:.5})", shown(f.x), shown(f.y), shown(f.z));
        }
    }
    let (count, average, min, max) = bond_stats(idx, base, scene, s, radii);
    text += &format!("\nCoordination number:  {count}\n");
    if count > 0 {
        text += &format!("Avg bond length:  {average:.4} Å\nMin bond length:  {min:.4} Å\nMax bond length:  {max:.4} Å");
    } else {
        text += "Avg bond length:  N/A\nMin bond length:  N/A\nMax bond length:  N/A";
    }
    text
}
