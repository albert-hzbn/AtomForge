//! Cell matrices, fractional coordinates and minimum images (the C++
//! math/StructureMath.h, in double precision).

use glam::{DMat3, DVec3};

use crate::structure::Structure;

/// Bond tolerance and shortest bond of the bond search.
pub const BOND_TOLERANCE_FACTOR: f32 = 1.18;
pub const MIN_BOND_DISTANCE: f32 = 0.10;

/// The cell as a matrix whose columns are a, b and c.
pub fn cell_matrix(s: &Structure) -> DMat3 {
    let [a, b, c] = s.cell_vectors;
    DMat3::from_cols(DVec3::from(a), DVec3::from(b), DVec3::from(c))
}

/// Cell and inverse when the structure has a cell with volume.
pub fn cell_matrices(s: &Structure) -> Option<(DMat3, DMat3)> {
    if !s.has_unit_cell {
        return None;
    }
    let cell = cell_matrix(s);
    let det = crate::glm::determinant_d(&cell);
    if !det.is_finite() || det.abs() < 1e-12 {
        return None;
    }
    Some((cell, crate::glm::inverse_d(&cell)))
}

/// Fractional coordinates relative to the cell origin.
pub fn cartesian_to_fractional(s: &Structure, cart: DVec3) -> Option<DVec3> {
    let (_, inverse) = cell_matrices(s)?;
    Some(inverse * (cart - DVec3::from(s.cell_offset)))
}

pub fn fractional_to_cartesian(s: &Structure, frac: DVec3) -> Option<DVec3> {
    let (cell, _) = cell_matrices(s)?;
    Some(DVec3::from(s.cell_offset) + cell * frac)
}

/// The shortest image of a displacement along periodic directions (rounding
/// each fractional component, as the C++ minimumImageDelta).
pub fn minimum_image_delta(delta: DVec3, periodic: [bool; 3], cell: &DMat3, inverse: &DMat3) -> DVec3 {
    if !periodic.iter().any(|&p| p) {
        return delta;
    }
    let mut f = *inverse * delta;
    for k in 0..3 {
        if periodic[k] {
            f[k] -= f[k].round();
        }
    }
    *cell * f
}

pub fn volume(s: &Structure) -> f64 {
    cell_matrix(s).determinant().abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fractional_round_trip() {
        let mut s = Structure::default();
        s.has_unit_cell = true;
        s.cell_vectors = [[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.5, 0.2, 5.0]];
        s.cell_offset = [1.0, -1.0, 0.5];
        let p = DVec3::new(2.3, 0.4, 1.9);
        let f = cartesian_to_fractional(&s, p).unwrap();
        let back = fractional_to_cartesian(&s, f).unwrap();
        assert!((back - p).length() < 1e-12);
        assert!((volume(&s) - 60.0).abs() < 1e-12);
    }
}
