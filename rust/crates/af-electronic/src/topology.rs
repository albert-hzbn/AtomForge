//! Betti curves of a scalar field's superlevel sets (the C++
//! electronic/Topology): components (betti0), loops (betti1, from the Euler
//! characteristic) and enclosed cavities (betti2) of {value >= t} on finite grids.

use std::collections::HashSet;

use glam::DVec4;

use crate::volume::Grid;

struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self { parent: (0..n).collect() }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    fn unite(&mut self, a: usize, b: usize) {
        let a = self.find(a);
        let b = self.find(b);
        if a != b {
            self.parent[a] = b;
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BettiNumbers {
    pub betti0: i32,
    pub betti1: i32,
    pub betti2: i32,
}

pub fn betti_numbers(field: &Grid, threshold: f64) -> Result<BettiNumbers, String> {
    field.validate()?;
    if field.periodic {
        return Err("Betti-number topology currently supports finite (non-periodic) grids only".into());
    }
    if !threshold.is_finite() {
        return Err("threshold must be finite".into());
    }
    let [nx, ny, nz] = field.shape;
    let solid_at = |x: i32, y: i32, z: i32| x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz && field.values[field.index(x, y, z)] >= threshold;
    let voxel = |x: i32, y: i32, z: i32| ((z * ny + y) * nx + x) as usize;
    let total = nx as usize * ny as usize * nz as usize;
    // Betti0: face-connected solid components; betti2: background components
    // that never touch the grid boundary (enclosed cavities).
    let (mut solid_uf, mut bg_uf) = (UnionFind::new(total), UnionFind::new(total));
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if solid_at(x, y, z) {
                    if solid_at(x + 1, y, z) {
                        solid_uf.unite(voxel(x, y, z), voxel(x + 1, y, z));
                    }
                    if solid_at(x, y + 1, z) {
                        solid_uf.unite(voxel(x, y, z), voxel(x, y + 1, z));
                    }
                    if solid_at(x, y, z + 1) {
                        solid_uf.unite(voxel(x, y, z), voxel(x, y, z + 1));
                    }
                } else {
                    if x + 1 < nx && !solid_at(x + 1, y, z) {
                        bg_uf.unite(voxel(x, y, z), voxel(x + 1, y, z));
                    }
                    if y + 1 < ny && !solid_at(x, y + 1, z) {
                        bg_uf.unite(voxel(x, y, z), voxel(x, y + 1, z));
                    }
                    if z + 1 < nz && !solid_at(x, y, z + 1) {
                        bg_uf.unite(voxel(x, y, z), voxel(x, y, z + 1));
                    }
                }
            }
        }
    }
    let (mut solid_roots, mut background, mut unbounded) = (HashSet::new(), HashSet::new(), HashSet::new());
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if solid_at(x, y, z) {
                    solid_roots.insert(solid_uf.find(voxel(x, y, z)));
                    continue;
                }
                let root = bg_uf.find(voxel(x, y, z));
                background.insert(root);
                if x == 0 || x == nx - 1 || y == 0 || y == ny - 1 || z == 0 || z == nz - 1 {
                    unbounded.insert(root);
                }
            }
        }
    }
    let betti2 = background.iter().filter(|r| !unbounded.contains(*r)).count() as i32;
    // Euler characteristic over the closed unit cubes of the solid voxels, in
    // doubled coordinates: odd coordinates count the cell's dimension.
    let stride_x = 2 * nx as i64 + 3;
    let stride_y = 2 * ny as i64 + 3;
    let mut seen = HashSet::new();
    let mut counts = [0i64; 4];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !solid_at(x, y, z) {
                    continue;
                }
                for dz in 0..3 {
                    for dy in 0..3 {
                        for dx in 0..3 {
                            let (cx, cy, cz) = (2 * x + dx, 2 * y + dy, 2 * z + dz);
                            let key = cx as i64 + stride_x * (cy as i64 + stride_y * cz as i64);
                            if seen.insert(key) {
                                counts[((cx & 1) + (cy & 1) + (cz & 1)) as usize] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    let euler = counts[0] - counts[1] + counts[2] - counts[3];
    let betti0 = solid_roots.len() as i32;
    Ok(BettiNumbers { betti0, betti2, betti1: betti0 + betti2 - euler as i32 })
}

/// One (threshold, betti0, betti1, betti2) row per threshold.
pub fn betti_curve(field: &Grid, thresholds: &[f64]) -> Result<Vec<DVec4>, String> {
    if thresholds.is_empty() {
        return Err("Supply at least one threshold".into());
    }
    thresholds
        .iter()
        .map(|&t| {
            let n = betti_numbers(field, t)?;
            Ok(DVec4::new(t, n.betti0 as f64, n.betti1 as f64, n.betti2 as f64))
        })
        .collect()
}
