//! Ambient occlusion per atom (the C++ graphics/AmbientOcclusion): the
//! fraction of the sky an atom's neighbours hide from it.

use std::collections::HashMap;
use std::sync::OnceLock;

use glam::{IVec3, Vec3};
use rayon::prelude::*;

fn cell_of(p: Vec3, size: f32) -> IVec3 {
    IVec3::new((p.x / size).floor() as i32, (p.y / size).floor() as i32, (p.z / size).floor() as i32)
}

/// The 48 directions sampled: a Fibonacci lattice over the sphere.
pub fn occlusion_directions() -> &'static [Vec3] {
    static DIRECTIONS: OnceLock<Vec<Vec3>> = OnceLock::new();
    DIRECTIONS.get_or_init(|| {
        const COUNT: i32 = 48;
        let golden = 3.14159265358979f32 * (3.0 - 5.0f32.sqrt());
        (0..COUNT)
            .map(|i| {
                let y = 1.0 - 2.0 * (i as f32 + 0.5) / COUNT as f32;
                let r = (1.0 - y * y).max(0.0).sqrt();
                let a = golden * i as f32;
                Vec3::new(r * a.cos(), y, r * a.sin())
            })
            .collect()
    })
}

/// The blocked fraction of rays from each atom (0 open, 1 enclosed): rays
/// count as blocked when they meet another sphere within `reach` of the
/// atom's surface. Empty when there are more than `max_atoms` atoms.
pub fn ambient_occlusion(centres: &[Vec3], radii: &[f32], reach: f32, max_atoms: usize) -> Vec<f32> {
    let n = centres.len();
    if n == 0 || n > max_atoms || radii.len() != n || !(reach > 0.0) {
        return Vec::new();
    }
    let largest = radii.iter().fold(0.0f32, |m, &r| m.max(r));
    // Every sphere a ray can meet lies within 2 * largest + reach of the centre.
    let cell_size = (2.0 * largest + reach).max(1e-3);
    let mut grid: HashMap<IVec3, Vec<usize>> = HashMap::with_capacity(n);
    for (i, &c) in centres.iter().enumerate() {
        grid.entry(cell_of(c, cell_size)).or_default().push(i);
    }
    let directions = occlusion_directions();
    let one = |i: usize| -> f32 {
        let c = centres[i];
        let limit = radii[i] + reach;
        let mut near = Vec::new();
        let home = cell_of(c, cell_size);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(list) = grid.get(&(home + IVec3::new(dx, dy, dz))) {
                        for &j in list {
                            if j != i && (centres[j] - c).length() < limit + radii[j] {
                                near.push(j);
                            }
                        }
                    }
                }
            }
        }
        let mut blocked = vec![false; directions.len()];
        for &j in &near {
            let v = centres[j] - c;
            let distance2 = v.dot(v);
            let rj2 = radii[j] * radii[j];
            for (k, d) in directions.iter().enumerate() {
                if blocked[k] {
                    continue;
                }
                let t = v.dot(*d);
                if t <= 0.0 {
                    continue;
                }
                let miss2 = distance2 - t * t;
                if miss2 >= rj2 {
                    continue;
                }
                // The first point of the ray inside sphere j lies within the reach.
                if t - (rj2 - miss2).sqrt() < limit {
                    blocked[k] = true;
                }
            }
        }
        blocked.iter().filter(|&&b| b).count() as f32 / directions.len() as f32
    };
    if n < 2000 {
        (0..n).map(one).collect()
    } else {
        (0..n).into_par_iter().map(one).collect()
    }
}

/// Darkens each colour by strength x its occlusion, rays reaching four times
/// the mean radius (nothing for strength 0 or mismatched sizes).
pub fn darken_occluded(centres: &[Vec3], radii: &[f32], colours: &mut [[f32; 3]], strength: f32) {
    if !(strength > 0.0) || centres.is_empty() || radii.len() != centres.len() || colours.len() != centres.len() {
        return;
    }
    let mean = radii.iter().sum::<f32>() / radii.len() as f32;
    let occlusion = ambient_occlusion(centres, radii, 4.0 * mean, 2_000_000);
    for (colour, o) in colours.iter_mut().zip(occlusion) {
        let f = 1.0 - strength.min(1.0) * o;
        for c in colour.iter_mut() {
            *c *= f;
        }
    }
}
