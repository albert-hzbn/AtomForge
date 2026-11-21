//! The geometry and colours behind the electronic post-processing previews
//! (the C++ electronic/SlicePlane, electronic/DisplayRange and
//! ElectronicViewport::color): section planes through a grid's cell, their
//! clipped outline and samples, display ranges and the shared palettes.

use af_core::glm::{inverse_d, mul_d, normalize_d};
use af_electronic::Grid;
use glam::{DMat3, DVec2, DVec3, Vec3};

/// The value range of a field and a suggested isovalue (a display
/// heuristic: zero for signed fields, else the 90th percentile of a sample).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayRange {
    pub low: f64,
    pub high: f64,
    pub suggested: f64,
}

impl DisplayRange {
    pub fn contains(&self, value: f64) -> bool {
        value.is_finite() && value >= self.low && value <= self.high
    }
}

pub fn display_range(values: &[f64]) -> Result<DisplayRange, String> {
    if values.is_empty() {
        return Err("Cannot estimate an empty field".into());
    }
    // std::minmax_element: the first minimum and the last maximum (equal values either way).
    let low = values.iter().copied().fold(f64::INFINITY, f64::min);
    let high = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut r = DisplayRange { low, high, suggested: low };
    if low == high {
        return Ok(r);
    }
    if low < 0.0 && high > 0.0 {
        r.suggested = 0.0;
        return Ok(r);
    }
    let stride = 1.max(values.len().div_ceil(262144));
    let mut sample: Vec<f64> = values.iter().step_by(stride).copied().collect();
    let index = (sample.len() - 1) * 9 / 10;
    let (_, nth, _) = sample.select_nth_unstable_by(index, |a, b| a.total_cmp(b));
    r.suggested = *nth;
    if r.suggested <= r.low || r.suggested >= r.high {
        r.suggested = r.low * 0.95 + r.high * 0.05;
    }
    Ok(r)
}

/// ElectronicViewport::color: palette 0 spectrum, 1 blue-white-red, 2 sequential blue.
pub fn palette_colour(t: f32, palette: i32) -> Vec3 {
    let t = t.clamp(0.0, 1.0);
    let mix = |a: Vec3, b: Vec3, t: f32| a * (1.0 - t) + b * t;
    if palette == 1 {
        return if t < 0.5 { mix(Vec3::new(0.10, 0.28, 0.80), Vec3::splat(0.96), 2.0 * t) } else { mix(Vec3::splat(0.96), Vec3::new(0.80, 0.10, 0.12), 2.0 * t - 1.0) };
    }
    if palette == 2 {
        return mix(Vec3::new(0.91, 0.95, 1.0), Vec3::new(0.08, 0.24, 0.70), t);
    }
    let stops = [Vec3::new(0.10, 0.20, 0.85), Vec3::new(0.0, 0.75, 0.95), Vec3::new(0.12, 0.8, 0.3), Vec3::new(1.0, 0.85, 0.08), Vec3::new(0.85, 0.10, 0.08)];
    let i = 3.min((t * 4.0) as usize);
    mix(stops[i], stops[i + 1], t * 4.0 - i as f32)
}

/// A section: its rectangle (origin + spans) and the outline inside the cell.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SlicePlane {
    pub origin: DVec3,
    pub u: DVec3,
    pub v: DVec3,
    pub boundary: Vec<DVec3>,
}

/// A convex polygon clipped to the grid's cell.
pub fn clip_to_cell(grid: &Grid, polygon: Vec<DVec3>) -> Vec<DVec3> {
    let inverse = inverse_d(&grid.cell);
    let mut polygon: Vec<DVec3> = polygon.into_iter().map(|p| mul_d(&inverse, p - grid.origin)).collect();
    for axis in 0..3 {
        for side in 0..2 {
            if polygon.is_empty() {
                return Vec::new();
            }
            let distance = |p: DVec3| {
                let d = if side == 1 { 1.0 - p[axis] } else { p[axis] };
                if d.abs() < 1e-12 {
                    0.0
                } else {
                    d
                }
            };
            let mut clipped = Vec::new();
            let mut previous = *polygon.last().unwrap();
            let mut before = distance(previous);
            for &current in &polygon {
                let after = distance(current);
                if (before >= 0.0) != (after >= 0.0) {
                    clipped.push(previous + (current - previous) * (before / (before - after)));
                }
                if after >= 0.0 {
                    clipped.push(current);
                }
                previous = current;
                before = after;
            }
            polygon.clear();
            for p in clipped {
                if polygon.last().is_none_or(|&q: &DVec3| (p - q).length() > 1e-10) {
                    polygon.push(p);
                }
            }
            if polygon.len() > 1 && (polygon[0] - *polygon.last().unwrap()).length() < 1e-10 {
                polygon.pop();
            }
        }
    }
    polygon.into_iter().map(|p| grid.origin + mul_d(&grid.cell, p)).collect()
}

/// Samples of a grid on a section's rectangle (two layers, so it is a valid
/// grid); finite grids are clamped to their cell for samples outside it.
pub fn sample_slice_plane(grid: &Grid, plane: &SlicePlane, nu: i32, nv: i32) -> Result<Grid, String> {
    grid.validate()?;
    if nu < 2 || nv < 2 || nu as f64 * nv as f64 > 10_000_000.0 {
        return Err("Invalid section preview dimensions".into());
    }
    let mut out = Grid {
        shape: [nu, nv, 2],
        origin: plane.origin,
        cell: DMat3::from_cols(plane.u, plane.v, normalize_d(plane.u.cross(plane.v))),
        unit: grid.unit.clone(),
        name: "section".into(),
        values: vec![0.0; nu as usize * nv as usize * 2],
        periodic: false,
    };
    let inverse = inverse_d(&grid.cell);
    for y in 0..nv {
        for x in 0..nu {
            let mut point = out.position(x, y, 0);
            if !grid.periodic {
                point = grid.origin + mul_d(&grid.cell, mul_d(&inverse, point - grid.origin).clamp(DVec3::ZERO, DVec3::ONE));
            }
            let value = grid.sample(point)?;
            let (a, b) = (out.index(x, y, 0), out.index(x, y, 1));
            out.values[a] = value;
            out.values[b] = value;
        }
    }
    out.validate()?;
    Ok(out)
}

/// A segment clipped to the cell, or None when it lies outside.
pub fn clip_segment_to_cell(grid: &Grid, a: DVec3, b: DVec3) -> Option<(DVec3, DVec3)> {
    let inverse = inverse_d(&grid.cell);
    let start = mul_d(&inverse, a - grid.origin);
    let direction = mul_d(&inverse, b - a);
    let (mut enter, mut leave) = (0.0f64, 1.0f64);
    for axis in 0..3 {
        if direction[axis].abs() < 1e-14 {
            if start[axis] < 0.0 || start[axis] > 1.0 {
                return None;
            }
        } else {
            let x = -start[axis] / direction[axis];
            let y = (1.0 - start[axis]) / direction[axis];
            enter = enter.max(x.min(y));
            leave = leave.min(x.max(y));
        }
    }
    if leave <= enter {
        return None;
    }
    let delta = b - a;
    Some((a + delta * enter, a + delta * leave))
}

/// An arbitrary section: the point is fractional (relative to the cell
/// origin), the normal Cartesian; the rectangle spans the cell's projection.
pub fn slice_plane(grid: &Grid, point: DVec3, normal: DVec3) -> Result<SlicePlane, String> {
    for i in 0..3 {
        if !point[i].is_finite() || !normal[i].is_finite() {
            return Err("Plane point and normal must be finite.".into());
        }
    }
    if normal.length() < 1e-12 {
        return Err("Plane normal must be nonzero.".into());
    }
    let normal = normalize_d(normal);
    let ex = normalize_d(normal.cross(if normal.z.abs() < 0.9 { DVec3::Z } else { DVec3::Y }));
    let ey = normal.cross(ex);
    let center = grid.origin + mul_d(&grid.cell, point);
    let (mut low, mut high) = (DVec2::splat(1e300), DVec2::splat(-1e300));
    for z in 0..2 {
        for y in 0..2 {
            for x in 0..2 {
                let delta = grid.origin + mul_d(&grid.cell, DVec3::new(x as f64, y as f64, z as f64)) - center;
                let p = DVec2::new(delta.dot(ex), delta.dot(ey));
                low = low.min(p);
                high = high.max(p);
            }
        }
    }
    let origin = center + ex * low.x + ey * low.y;
    let (u, v) = (ex * (high.x - low.x), ey * (high.y - low.y));
    let mut plane = SlicePlane { origin, u, v, boundary: clip_to_cell(grid, vec![origin, origin + u, origin + u + v, origin + v]) };
    let mut area = DVec3::ZERO;
    for i in 1..plane.boundary.len().saturating_sub(1) {
        area += (plane.boundary[i] - plane.boundary[0]).cross(plane.boundary[i + 1] - plane.boundary[0]);
    }
    if area.length() < 1e-12 * u.length() * v.length() {
        plane.boundary.clear();
    }
    Ok(plane)
}

/// A lattice plane through the cell at a fraction along one axis (bc, ac or ab).
pub fn axis_plane(grid: &Grid, axis: usize, position: f64) -> SlicePlane {
    let a = if axis == 0 { 1 } else { 0 };
    let b = if axis == 2 { 1 } else { 2 };
    let origin = grid.origin + grid.cell.col(axis) * position;
    let (u, v) = (grid.cell.col(a), grid.cell.col(b));
    SlicePlane { origin, u, v, boundary: vec![origin, origin + u, origin + u + v, origin + v] }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube(n: i32, periodic: bool) -> Grid {
        let mut g = Grid { shape: [n, n, n], cell: DMat3::from_diagonal(DVec3::splat(4.0)), periodic, ..Default::default() };
        g.values = (0..n * n * n).map(|i| (i % n) as f64).collect();
        g
    }

    #[test]
    fn display_range_suggests_zero_for_signed_and_percentile_otherwise() {
        assert_eq!(display_range(&[-1.0, 2.0]).unwrap().suggested, 0.0);
        let values: Vec<f64> = (0..=100).map(|i| i as f64).collect();
        let r = display_range(&values).unwrap();
        assert_eq!((r.low, r.high, r.suggested), (0.0, 100.0, 90.0));
        assert_eq!(display_range(&[3.0, 3.0]).unwrap().suggested, 3.0);
        // The percentile at the top falls back to 5% of the range.
        assert!((display_range(&[0.0, 1.0, 1.0, 1.0]).unwrap().suggested - 0.05).abs() < 1e-12);
        assert!(display_range(&[]).is_err());
    }

    #[test]
    fn diagonal_plane_through_a_cube_is_a_hexagon() {
        let g = cube(4, true);
        let p = slice_plane(&g, DVec3::splat(0.5), DVec3::ONE).unwrap();
        assert_eq!(p.boundary.len(), 6);
        for q in &p.boundary {
            assert!((q.x + q.y + q.z - 6.0).abs() < 1e-9);
        }
        // A plane outside the cell has no outline.
        assert!(slice_plane(&g, DVec3::new(2.0, 2.0, 2.0), DVec3::ONE).unwrap().boundary.is_empty());
        assert!(slice_plane(&g, DVec3::splat(0.5), DVec3::ZERO).is_err());
    }

    #[test]
    fn segments_and_polygons_clip_to_the_cell() {
        let g = cube(4, false);
        let (a, b) = clip_segment_to_cell(&g, DVec3::new(-2.0, 1.0, 1.0), DVec3::new(6.0, 1.0, 1.0)).unwrap();
        assert!((a.x - 0.0).abs() < 1e-12 && (b.x - 4.0).abs() < 1e-12);
        assert!(clip_segment_to_cell(&g, DVec3::new(-2.0, 5.0, 1.0), DVec3::new(6.0, 5.0, 1.0)).is_none());
        let square = vec![DVec3::new(-1.0, -1.0, 2.0), DVec3::new(5.0, -1.0, 2.0), DVec3::new(5.0, 5.0, 2.0), DVec3::new(-1.0, 5.0, 2.0)];
        let clipped = clip_to_cell(&g, square);
        assert_eq!(clipped.len(), 4);
        assert!(clipped.iter().all(|p| p.x >= -1e-12 && p.x <= 4.0 + 1e-12 && p.y >= -1e-12 && p.y <= 4.0 + 1e-12));
    }

    #[test]
    fn section_samples_follow_the_field() {
        // values = x index on a periodic 4-grid with 1 A steps: f(x) = x (periodic wrap at 4).
        let g = cube(4, true);
        let plane = axis_plane(&g, 2, 0.25);
        let s = sample_slice_plane(&g, &plane, 5, 5).unwrap();
        assert_eq!(s.shape, [5, 5, 2]);
        assert!((s.values[s.index(1, 0, 0)] - 1.0).abs() < 1e-12);
        assert!((s.values[s.index(2, 3, 1)] - 2.0).abs() < 1e-12);
        assert!(sample_slice_plane(&g, &plane, 1, 5).is_err());
    }

    #[test]
    fn palettes_match_the_cpp_ramps() {
        assert_eq!(palette_colour(0.0, 0), Vec3::new(0.10, 0.20, 0.85));
        assert_eq!(palette_colour(1.0, 0), Vec3::new(0.85, 0.10, 0.08));
        assert_eq!(palette_colour(0.5, 1), Vec3::splat(0.96));
        assert_eq!(palette_colour(2.0, 2), Vec3::new(0.08, 0.24, 0.70));
    }
}
