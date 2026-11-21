//! Grid integration (the C++ electronic/DensityIntegration): compensated
//! sums, spheres, Voronoi basins and local maxima.

use glam::{DVec3, IVec3};

use crate::volume::{minimum_image, Grid};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Basin {
    pub integral: f64,
    pub volume: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Peak {
    pub position: DVec3,
    pub value: f64,
}

/// The integral over the cell (Kahan-compensated, trapezoid ends on open grids).
pub fn integrate(grid: &Grid) -> Result<f64, String> {
    grid.validate()?;
    let (mut sum, mut correction) = (0.0f64, 0.0f64);
    for z in 0..grid.shape[2] {
        for y in 0..grid.shape[1] {
            for x in 0..grid.shape[0] {
                let v = grid.values[grid.index(x, y, z)] * grid.weight(x, y, z) - correction;
                let next = sum + v;
                correction = (next - sum) - v;
                sum = next;
            }
        }
    }
    Ok(sum)
}

pub fn integrate_sphere(grid: &Grid, center: DVec3, radius: f64) -> Result<f64, String> {
    grid.validate()?;
    if !radius.is_finite() || radius <= 0.0 {
        return Err("Sphere radius must be positive".into());
    }
    if !center.is_finite() {
        return Err("Non-finite sphere center".into());
    }
    let mut sum = 0.0;
    for z in 0..grid.shape[2] {
        for y in 0..grid.shape[1] {
            for x in 0..grid.shape[0] {
                let mut d = grid.position(x, y, z) - center;
                if grid.periodic {
                    d = minimum_image(d, &grid.cell)?;
                }
                if d.dot(d) <= radius * radius {
                    sum += grid.values[grid.index(x, y, z)] * grid.weight(x, y, z);
                }
            }
        }
    }
    Ok(sum)
}

/// Charge and volume of each site's Voronoi cell (ties shared equally).
pub fn voronoi_integrate(grid: &Grid, sites: &[DVec3]) -> Result<Vec<Basin>, String> {
    grid.validate()?;
    if sites.is_empty() {
        return Err("Voronoi integration requires sites".into());
    }
    if sites.iter().any(|s| !s.is_finite()) {
        return Err("Non-finite Voronoi site".into());
    }
    if sites.len() as f64 * grid.values.len() as f64 > 1e9 {
        return Err("Voronoi workload too large; reduce grid or site count".into());
    }
    let mut out = vec![Basin::default(); sites.len()];
    for z in 0..grid.shape[2] {
        for y in 0..grid.shape[1] {
            for x in 0..grid.shape[0] {
                let mut best = f64::INFINITY;
                let mut nearest: Vec<usize> = Vec::new();
                let p = grid.position(x, y, z);
                for (i, site) in sites.iter().enumerate() {
                    let mut d = p - *site;
                    if grid.periodic {
                        d = minimum_image(d, &grid.cell)?;
                    }
                    let d2 = d.dot(d);
                    if d2 < best - 1e-10 {
                        best = d2;
                        nearest = vec![i];
                    } else if (d2 - best).abs() <= 1e-10 {
                        nearest.push(i);
                    }
                }
                let w = grid.weight(x, y, z) / nearest.len() as f64;
                for i in nearest {
                    out[i].volume += w;
                    out[i].integral += w * grid.values[grid.index(x, y, z)];
                }
            }
        }
    }
    Ok(out)
}

/// Strict local maxima (ties broken by index), largest first; `limit` 0: all.
pub fn peaks(grid: &Grid, limit: i32) -> Result<Vec<Peak>, String> {
    grid.validate()?;
    if limit < 0 {
        return Err("Peak limit must be nonnegative".into());
    }
    let mut out = Vec::new();
    for z in 0..grid.shape[2] {
        for y in 0..grid.shape[1] {
            for x in 0..grid.shape[0] {
                let index = grid.index(x, y, z);
                let value = grid.values[index];
                let (mut maximum, mut strictly_greater) = (true, false);
                for dz in -1..=1 {
                    for dy in -1..=1 {
                        for dx in -1..=1 {
                            let mut p = IVec3::new(x + dx, y + dy, z + dz);
                            let mut valid = true;
                            for a in 0..3 {
                                if grid.periodic {
                                    p[a] = (p[a] + grid.shape[a]) % grid.shape[a];
                                } else if p[a] < 0 || p[a] >= grid.shape[a] {
                                    valid = false;
                                }
                            }
                            if !valid {
                                continue;
                            }
                            let j = grid.index(p.x, p.y, p.z);
                            if grid.values[j] > value || (grid.values[j] == value && j < index) {
                                maximum = false;
                            }
                            if grid.values[j] < value {
                                strictly_greater = true;
                            }
                        }
                    }
                }
                if maximum && strictly_greater {
                    out.push(Peak { position: grid.position(x, y, z), value });
                }
            }
        }
    }
    out.sort_by(|a, b| b.value.partial_cmp(&a.value).unwrap_or(std::cmp::Ordering::Equal));
    if limit > 0 && out.len() > limit as usize {
        out.truncate(limit as usize);
    }
    Ok(out)
}
