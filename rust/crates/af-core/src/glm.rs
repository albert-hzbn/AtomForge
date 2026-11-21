//! glm's single- and double-precision vector and matrix operations, term by
//! term, for algorithms whose C++ uses glm: glm::determinant, glm::inverse,
//! glm::normalize and matrix products give the same bits as the C++ build.
//! Matrices are glam's (column-major, as glm: m[c][r] is column c, row r).

use glam::{DMat3, DVec3, Mat3, Vec3};

/// glm::determinant(mat3).
pub fn determinant(m: &Mat3) -> f32 {
    let c = |i: usize, j: usize| m.col(i)[j];
    c(0, 0) * (c(1, 1) * c(2, 2) - c(2, 1) * c(1, 2)) - c(1, 0) * (c(0, 1) * c(2, 2) - c(2, 1) * c(0, 2))
        + c(2, 0) * (c(0, 1) * c(1, 2) - c(1, 1) * c(0, 2))
}

/// glm::inverse(mat3): the adjugate times 1/determinant.
pub fn inverse(m: &Mat3) -> Mat3 {
    let c = |i: usize, j: usize| m.col(i)[j];
    let one_over = 1.0 / determinant(m);
    let mut r = [[0.0f32; 3]; 3];
    r[0][0] = (c(1, 1) * c(2, 2) - c(2, 1) * c(1, 2)) * one_over;
    r[1][0] = -(c(1, 0) * c(2, 2) - c(2, 0) * c(1, 2)) * one_over;
    r[2][0] = (c(1, 0) * c(2, 1) - c(2, 0) * c(1, 1)) * one_over;
    r[0][1] = -(c(0, 1) * c(2, 2) - c(2, 1) * c(0, 2)) * one_over;
    r[1][1] = (c(0, 0) * c(2, 2) - c(2, 0) * c(0, 2)) * one_over;
    r[2][1] = -(c(0, 0) * c(2, 1) - c(2, 0) * c(0, 1)) * one_over;
    r[0][2] = (c(0, 1) * c(1, 2) - c(1, 1) * c(0, 2)) * one_over;
    r[1][2] = -(c(0, 0) * c(1, 2) - c(1, 0) * c(0, 2)) * one_over;
    r[2][2] = (c(0, 0) * c(1, 1) - c(1, 0) * c(0, 1)) * one_over;
    Mat3::from_cols_array_2d(&r)
}

/// glm::determinant(dmat3).
pub fn determinant_d(m: &DMat3) -> f64 {
    let c = |i: usize, j: usize| m.col(i)[j];
    c(0, 0) * (c(1, 1) * c(2, 2) - c(2, 1) * c(1, 2)) - c(1, 0) * (c(0, 1) * c(2, 2) - c(2, 1) * c(0, 2))
        + c(2, 0) * (c(0, 1) * c(1, 2) - c(1, 1) * c(0, 2))
}

/// glm::inverse(dmat3).
pub fn inverse_d(m: &DMat3) -> DMat3 {
    let c = |i: usize, j: usize| m.col(i)[j];
    let one_over = 1.0 / determinant_d(m);
    let mut r = [[0.0f64; 3]; 3];
    r[0][0] = (c(1, 1) * c(2, 2) - c(2, 1) * c(1, 2)) * one_over;
    r[1][0] = -(c(1, 0) * c(2, 2) - c(2, 0) * c(1, 2)) * one_over;
    r[2][0] = (c(1, 0) * c(2, 1) - c(2, 0) * c(1, 1)) * one_over;
    r[0][1] = -(c(0, 1) * c(2, 2) - c(2, 1) * c(0, 2)) * one_over;
    r[1][1] = (c(0, 0) * c(2, 2) - c(2, 0) * c(0, 2)) * one_over;
    r[2][1] = -(c(0, 0) * c(2, 1) - c(2, 0) * c(0, 1)) * one_over;
    r[0][2] = (c(0, 1) * c(1, 2) - c(1, 1) * c(0, 2)) * one_over;
    r[1][2] = -(c(0, 0) * c(1, 2) - c(1, 0) * c(0, 2)) * one_over;
    r[2][2] = (c(0, 0) * c(1, 1) - c(1, 0) * c(0, 1)) * one_over;
    DMat3::from_cols_array_2d(&r)
}

/// m * v as glm computes it: m[0]*v.x + m[1]*v.y + m[2]*v.z.
pub fn mul(m: &Mat3, v: Vec3) -> Vec3 {
    m.x_axis * v.x + m.y_axis * v.y + m.z_axis * v.z
}

pub fn mul_d(m: &DMat3, v: DVec3) -> DVec3 {
    m.x_axis * v.x + m.y_axis * v.y + m.z_axis * v.z
}

/// glm::dot: (x*x' + y*y') + z*z'.
pub fn dot(a: Vec3, b: Vec3) -> f32 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// glm::length: sqrt(dot(v, v)).
pub fn length(v: Vec3) -> f32 {
    dot(v, v).sqrt()
}

/// glm::normalize: v * inversesqrt(dot(v, v)), inversesqrt = 1 / sqrt.
pub fn normalize(v: Vec3) -> Vec3 {
    v * (1.0 / dot(v, v).sqrt())
}

pub fn normalize_d(v: DVec3) -> DVec3 {
    v * (1.0 / v.dot(v).sqrt())
}

/// glm::cross.
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y)
}

/// glm::round per component (half away from zero).
pub fn round(v: Vec3) -> Vec3 {
    Vec3::new(v.x.round(), v.y.round(), v.z.round())
}

/// glm::floor per component.
pub fn floor(v: Vec3) -> Vec3 {
    Vec3::new(v.x.floor(), v.y.floor(), v.z.floor())
}

/// glm::min / glm::max per component (std::min: b < a ? b : a).
pub fn min(a: Vec3, b: Vec3) -> Vec3 {
    let m = |x: f32, y: f32| if y < x { y } else { x };
    Vec3::new(m(a.x, b.x), m(a.y, b.y), m(a.z, b.z))
}

pub fn max(a: Vec3, b: Vec3) -> Vec3 {
    let m = |x: f32, y: f32| if x < y { y } else { x };
    Vec3::new(m(a.x, b.x), m(a.y, b.y), m(a.z, b.z))
}

/// The float cell matrix of a structure (columns a, b, c), as makeCellMatrix.
pub fn cell_matrix(s: &crate::Structure) -> Mat3 {
    let v = |r: usize| Vec3::new(s.cell_vectors[r][0] as f32, s.cell_vectors[r][1] as f32, s.cell_vectors[r][2] as f32);
    Mat3::from_cols(v(0), v(1), v(2))
}

/// tryMakeCellMatrices: the float cell and its inverse when the cell has volume.
pub fn cell_matrices(s: &crate::Structure) -> Option<(Mat3, Mat3)> {
    if !s.has_unit_cell {
        return None;
    }
    let cell = cell_matrix(s);
    let det = determinant(&cell);
    if !det.is_finite() || det.abs() <= 1e-8 {
        return None;
    }
    Some((cell, inverse(&cell)))
}

/// minimumImageDelta (float): round every fractional component when `use_pbc`.
pub fn minimum_image_delta(delta: Vec3, use_pbc: bool, cell: &Mat3, inverse: &Mat3) -> Vec3 {
    if !use_pbc {
        return delta;
    }
    let mut frac = mul(inverse, delta);
    frac -= round(frac);
    mul(cell, frac)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_is_an_inverse() {
        let m = Mat3::from_cols(Vec3::new(4.0, 0.1, 0.0), Vec3::new(1.0, 3.0, 0.2), Vec3::new(0.5, 0.2, 5.0));
        let p = m * inverse(&m);
        assert!((p - Mat3::IDENTITY).abs_diff_eq(Mat3::ZERO, 1e-6));
        assert!((determinant(&m) - m.determinant()).abs() < 1e-4);
    }
}
