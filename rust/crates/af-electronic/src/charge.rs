//! Charge-density difference analysis (the C++ electronic/ChargeAnalysis):
//! differences, masks, accumulation/depletion and cumulative charge.

use af_core::glm::determinant_d;
use glam::DVec2;

use crate::density::integrate;
use crate::field::planar_average;
use crate::volume::{arithmetic, require_aligned, scale, Grid};

fn require_mask(mask: &Grid) -> Result<(), String> {
    mask.validate()?;
    if mask.unit != "mask" {
        return Err("Create a threshold mask before using Boolean operations".into());
    }
    if mask.values.iter().any(|&v| v != 0.0 && v != 1.0) {
        return Err("Mask values must be zero or one".into());
    }
    Ok(())
}

fn require_density(grid: &Grid) -> Result<(), String> {
    grid.validate()?;
    if grid.unit != "e/A^3" {
        return Err("Charge analysis requires electron density in e/A^3".into());
    }
    Ok(())
}

/// a - weight * b.
pub fn density_difference(a: &Grid, b: &Grid, weight: f64) -> Result<Grid, String> {
    require_density(a)?;
    require_density(b)?;
    let mut out = arithmetic(a, &scale(b, weight)?, "subtract")?;
    out.name = "density_difference".into();
    Ok(out)
}

pub fn threshold_mask(grid: &Grid, low: f64, high: f64) -> Result<Grid, String> {
    grid.validate()?;
    if !low.is_finite() || !high.is_finite() || low > high {
        return Err("Mask bounds must be finite and ordered".into());
    }
    let mut out = grid.clone();
    out.unit = "mask".into();
    out.name = "threshold_mask".into();
    for v in out.values.iter_mut() {
        *v = if *v >= low && *v <= high { 1.0 } else { 0.0 };
    }
    Ok(out)
}

pub fn boolean_mask(a: &Grid, b: &Grid, operation: &str) -> Result<Grid, String> {
    require_aligned(a, b)?;
    require_mask(a)?;
    require_mask(b)?;
    if !["union", "intersection", "difference", "xor"].contains(&operation) {
        return Err("Unknown Boolean mask operation".into());
    }
    let mut out = a.clone();
    out.name = format!("mask_{operation}");
    for i in 0..out.values.len() {
        let (x, y) = (a.values[i] != 0.0, b.values[i] != 0.0);
        let r = match operation {
            "union" => x || y,
            "intersection" => x && y,
            "difference" => x && !y,
            _ => x != y,
        };
        out.values[i] = r as i32 as f64;
    }
    Ok(out)
}

pub fn apply_mask(grid: &Grid, mask: &Grid) -> Result<Grid, String> {
    require_aligned(grid, mask)?;
    require_mask(mask)?;
    let mut out = grid.clone();
    out.name = format!("masked_{}", grid.name);
    for (v, m) in out.values.iter_mut().zip(&mask.values) {
        *v *= m;
    }
    Ok(out)
}

/// Accumulation (positive part) and depletion magnitude (negative part).
pub fn split_density(grid: &Grid) -> Result<[Grid; 2], String> {
    require_density(grid)?;
    let mut out = [grid.clone(), grid.clone()];
    out[0].name = "accumulation".into();
    out[1].name = "depletion_magnitude".into();
    for i in 0..grid.values.len() {
        let v = grid.values[i];
        out[0].values[i] = if 0.0 < v { v } else { 0.0 };
        out[1].values[i] = if 0.0 < -v { -v } else { 0.0 };
    }
    Ok(out)
}

/// Positive accumulation, positive depletion magnitude, net electron change.
pub fn charge_summary(grid: &Grid) -> Result<[f64; 3], String> {
    let parts = split_density(grid)?;
    Ok([integrate(&parts[0])?, integrate(&parts[1])?, integrate(grid)?])
}

/// Normal distance and cumulative electron count from the cell's lower face.
pub fn cumulative_charge(grid: &Grid, axis: i32) -> Result<Vec<DVec2>, String> {
    require_density(grid)?;
    let profile = planar_average(grid, axis)?;
    let axis = axis as usize;
    let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
    let area = grid.cell.col(a).cross(grid.cell.col(b)).length();
    let height = determinant_d(&grid.cell).abs() / area;
    let step = height / (if grid.periodic { grid.shape[axis] } else { grid.shape[axis] - 1 }) as f64;
    let mut out = vec![DVec2::ZERO];
    let mut sum = 0.0;
    // Periodic bins use rectangle quadrature; finite samples use trapezoids,
    // matching integrate() and including the complete cell in the final row.
    if grid.periodic {
        for (i, p) in profile.iter().enumerate() {
            sum += p.y * area * step;
            out.push(DVec2::new((i + 1) as f64 * step, sum));
        }
    } else {
        for i in 1..profile.len() {
            sum += (profile[i - 1].y + profile[i].y) * 0.5 * area * step;
            out.push(DVec2::new(i as f64 * step, sum));
        }
    }
    Ok(out)
}
