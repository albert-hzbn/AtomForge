//! Native on-grid QTAIM/Bader partitioning (the C++ electronic/BaderPartition):
//! steepest ascent over the 26 neighbours to each voxel's density maximum
//! (Henkelman, Sanville and Jonsson 2006, on-grid step only).

use std::collections::HashMap;

use af_core::glm::mul_d;
use glam::DVec3;

use crate::volume::{minimum_image, Grid};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BaderBasin {
    pub charge: f64,
    pub volume: f64,
    pub maximum: DVec3,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BaderPartition {
    /// Basin of each voxel, x-fastest order.
    pub basin: Vec<i32>,
    pub basins: Vec<BaderBasin>,
}

fn neighbor_offsets() -> [[i32; 3]; 26] {
    let mut offsets = [[0; 3]; 26];
    let mut n = 0;
    for dz in -1..=1 {
        for dy in -1..=1 {
            for dx in -1..=1 {
                if dx != 0 || dy != 0 || dz != 0 {
                    offsets[n] = [dx, dy, dz];
                    n += 1;
                }
            }
        }
    }
    offsets
}

/// Assigns every voxel to the basin reached by repeatedly stepping to the
/// neighbour of steepest density increase per Cartesian distance.
pub fn bader_on_grid(density: &Grid) -> Result<BaderPartition, String> {
    density.validate()?;
    if density.unit != "e/A^3" {
        return Err("Bader partitioning requires electron density in e/A^3".into());
    }
    if density.values.iter().any(|&v| v < 0.0) {
        return Err("Bader partitioning requires nonnegative density".into());
    }
    let [nx, ny, nz] = density.shape;
    let offsets = neighbor_offsets();
    let steps = density.steps();
    let neighbor_distance: Vec<f64> = offsets.iter().map(|o| mul_d(&steps, DVec3::new(o[0] as f64, o[1] as f64, o[2] as f64)).length()).collect();
    let value_at = |x: i32, y: i32, z: i32| density.values[density.index(x, y, z)];
    let wrap = |v: i32, n: i32| ((v % n) + n) % n;
    let total = nx as usize * ny as usize * nz as usize;
    // next[i]: the steepest-ascent neighbour of voxel i, or -1 at a local maximum.
    let mut next = vec![-1i64; total];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let here = value_at(x, y, z);
                let mut best = 0.0;
                let mut best_neighbor = -1i64;
                for (i, o) in offsets.iter().enumerate() {
                    let (mut px, mut py, mut pz) = (x + o[0], y + o[1], z + o[2]);
                    if density.periodic {
                        px = wrap(px, nx);
                        py = wrap(py, ny);
                        pz = wrap(pz, nz);
                    } else if px < 0 || px >= nx || py < 0 || py >= ny || pz < 0 || pz >= nz {
                        continue;
                    }
                    let slope = (value_at(px, py, pz) - here) / neighbor_distance[i];
                    if slope > best {
                        best = slope;
                        best_neighbor = density.index(px, py, pz) as i64;
                    }
                }
                next[density.index(x, y, z)] = best_neighbor;
            }
        }
    }
    // Follow each ascent chain to its root (a local maximum), memoizing the path.
    let mut root_of = vec![-1i64; total];
    for start in 0..total {
        let mut path = Vec::new();
        let mut current = start;
        while root_of[current] == -1 && next[current] != -1 {
            path.push(current);
            current = next[current] as usize;
        }
        let root = if root_of[current] != -1 { root_of[current] } else { current as i64 };
        root_of[current] = root;
        for voxel in path {
            root_of[voxel] = root;
        }
    }
    let mut root_to_basin: HashMap<i64, usize> = HashMap::new();
    let mut root_list: Vec<i64> = Vec::new();
    let mut result = BaderPartition { basin: vec![0; total], basins: Vec::new() };
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = density.index(x, y, z);
                let id = *root_to_basin.entry(root_of[i]).or_insert_with(|| {
                    root_list.push(root_of[i]);
                    result.basins.push(BaderBasin::default());
                    root_list.len() - 1
                });
                result.basin[i] = id as i32;
                let w = density.weight(x, y, z);
                result.basins[id].charge += w * value_at(x, y, z);
                result.basins[id].volume += w;
            }
        }
    }
    for (i, &r) in root_list.iter().enumerate() {
        let r = r as i32;
        let (x, y, z) = (r % nx, (r / nx) % ny, r / (nx * ny));
        result.basins[i].maximum = density.position(x, y, z);
    }
    Ok(result)
}

/// Each basin's charge added to the site nearest its maximum (minimum image when periodic).
pub fn bader_populations(density: &Grid, sites: &[DVec3]) -> Result<Vec<f64>, String> {
    if sites.is_empty() {
        return Err("Supply at least one site".into());
    }
    let partition = bader_on_grid(density)?;
    let mut populations = vec![0.0; sites.len()];
    for basin in &partition.basins {
        let mut nearest = 0;
        let mut best = f64::INFINITY;
        for (s, site) in sites.iter().enumerate() {
            let mut delta = basin.maximum - *site;
            if density.periodic {
                delta = minimum_image(delta, &density.cell)?;
            }
            let d2 = delta.dot(delta);
            if d2 < best {
                best = d2;
                nearest = s;
            }
        }
        populations[nearest] += basin.charge;
    }
    Ok(populations)
}
