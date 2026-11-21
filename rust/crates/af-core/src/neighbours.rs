//! Neighbour pairs within a cutoff from a cell list (the C++
//! math/NeighbourSearch). The search order is the C++ one, so pair lists are
//! identical, not just equal as sets.

use std::collections::HashMap;

use glam::{DMat3, DVec3};

use crate::structure::Structure;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pair {
    pub i: usize,
    pub j: usize,
    /// From atom i to (the image of) atom j.
    pub d: DVec3,
    pub r: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum NeighbourError {
    NoVolume,
    TooManyImages(i32),
}

impl std::fmt::Display for NeighbourError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NeighbourError::NoVolume => write!(f, "The unit cell has no volume"),
            NeighbourError::TooManyImages(n) => write!(
                f,
                "The cutoff spans more than {n} periodic images of this cell; use a smaller cutoff or a supercell"
            ),
        }
    }
}

impl std::error::Error for NeighbourError {}

/// Calls `visit` for every ordered pair (i, j) with 0 < |d| <= cutoff,
/// grouped by i. Periodic: over all images of the cell (a cell smaller than
/// the cutoff lists every image of a neighbour); positions are wrapped first.
pub fn visit_within(
    structure: &Structure,
    cutoff: f64,
    periodic: bool,
    max_images: i32,
    mut visit: impl FnMut(&Pair),
) -> Result<(), NeighbourError> {
    let n = structure.atoms.len();
    if n == 0 || !(cutoff > 0.0) {
        return Ok(());
    }
    let mut cell = DMat3::IDENTITY;
    let mut inverse = DMat3::IDENTITY;
    let mut range = [0i32; 3];
    if periodic {
        cell = crate::cell::cell_matrix(structure);
        if crate::glm::determinant_d(&cell).abs() < 1e-12 {
            return Err(NeighbourError::NoVolume);
        }
        inverse = crate::glm::inverse_d(&cell);
        // Planes of constant fractional coordinate k lie 1/|row k of the inverse| apart.
        for k in 0..3 {
            let row = inverse.row(k);
            range[k] = (cutoff * row.length()).ceil() as i32;
            if range[k] > max_images {
                return Err(NeighbourError::TooManyImages(max_images));
            }
        }
    }
    let position: Vec<DVec3> = structure
        .atoms
        .iter()
        .map(|a| {
            let p = DVec3::new(a.x, a.y, a.z);
            if periodic {
                let mut f = inverse * p;
                f -= f.floor();
                cell * f
            } else {
                p
            }
        })
        .collect();
    let bin_of = |p: DVec3| -> [i64; 3] { [(p.x / cutoff).floor() as i64, (p.y / cutoff).floor() as i64, (p.z / cutoff).floor() as i64] };
    let mut bins: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (j, &p) in position.iter().enumerate() {
        bins.entry(bin_of(p)).or_default().push(j);
    }
    let cutoff2 = cutoff * cutoff;
    for i in 0..n {
        for a in -range[0]..=range[0] {
            for b in -range[1]..=range[1] {
                for c in -range[2]..=range[2] {
                    let shift = if periodic { cell * DVec3::new(a as f64, b as f64, c as f64) } else { DVec3::ZERO };
                    // Atom j's image at position[j] + shift is near i when position[j] is near position[i] - shift.
                    let centre = bin_of(position[i] - shift);
                    for dx in -1..=1 {
                        for dy in -1..=1 {
                            for dz in -1..=1 {
                                let Some(list) = bins.get(&[centre[0] + dx, centre[1] + dy, centre[2] + dz]) else { continue };
                                for &j in list {
                                    let d = position[j] + shift - position[i];
                                    let d2 = d.dot(d);
                                    if d2 > cutoff2 || d2 < 1e-12 {
                                        continue;
                                    }
                                    visit(&Pair { i, j, d, r: d2.sqrt() });
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// All pairs within the cutoff (at most 20 images of the cell searched).
pub fn within(structure: &Structure, cutoff: f64, periodic: bool) -> Result<Vec<Pair>, NeighbourError> {
    let mut pairs = Vec::new();
    visit_within(structure, cutoff, periodic, 20, |p| pairs.push(*p))?;
    Ok(pairs)
}

/// Calls `visit(i, j)` once for each unordered pair i < j that may lie within
/// the cutoff (the caller applies its own distance test); every pair when the
/// cell is too small to search.
pub fn for_each_candidate_pair(structure: &Structure, cutoff: f64, periodic: bool, mut visit: impl FnMut(usize, usize)) {
    let pairs = match within(structure, cutoff * 1.001 + 1e-3, periodic) {
        Ok(pairs) => pairs,
        Err(_) => {
            let n = structure.atoms.len();
            for i in 0..n {
                for j in i + 1..n {
                    visit(i, j);
                }
            }
            return;
        }
    };
    let mut k = 0;
    let mut partners = Vec::new();
    while k < pairs.len() {
        let i = pairs[k].i;
        partners.clear();
        while k < pairs.len() && pairs[k].i == i {
            if pairs[k].j > i {
                partners.push(pairs[k].j);
            }
            k += 1;
        }
        partners.sort_unstable();
        partners.dedup();
        for &j in &partners {
            visit(i, j);
        }
    }
}
