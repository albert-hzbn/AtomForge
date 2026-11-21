//! Eigen's fixed-size 3x3 formulas, term by term, for code whose C++ uses
//! Eigen::Matrix3d (matrices here are row-major: m[row][col]).

pub type Matrix3 = [[f64; 3]; 3];

fn cofactor(m: &Matrix3, i: usize, j: usize) -> f64 {
    let (i1, i2, j1, j2) = ((i + 1) % 3, (i + 2) % 3, (j + 1) % 3, (j + 2) % 3);
    m[i1][j1] * m[i2][j2] - m[i1][j2] * m[i2][j1]
}

/// Matrix3d::inverse(): the cofactors of column 0 give the determinant.
pub fn inverse3(m: &Matrix3) -> Matrix3 {
    let col0 = [cofactor(m, 0, 0), cofactor(m, 1, 0), cofactor(m, 2, 0)];
    let det = col0[0] * m[0][0] + col0[1] * m[1][0] + col0[2] * m[2][0];
    let invdet = 1.0 / det;
    let mut r = [[0.0; 3]; 3];
    r[1][2] = cofactor(m, 2, 1) * invdet;
    r[2][1] = cofactor(m, 1, 2) * invdet;
    r[2][2] = cofactor(m, 2, 2) * invdet;
    r[1][0] = cofactor(m, 0, 1) * invdet;
    r[1][1] = cofactor(m, 1, 1) * invdet;
    r[2][0] = cofactor(m, 0, 2) * invdet;
    for (j, c) in col0.iter().enumerate() {
        r[0][j] = c * invdet;
    }
    r
}

/// Matrix3d::determinant().
pub fn determinant3(m: &Matrix3) -> f64 {
    let helper = |a: usize, b: usize, c: usize| m[0][a] * (m[1][b] * m[2][c] - m[1][c] * m[2][b]);
    helper(0, 1, 2) - helper(1, 0, 2) + helper(2, 0, 1)
}

/// m * v.
pub fn mul3(m: &Matrix3, v: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2])
}

/// a * b (each entry summed over k = 0, 1, 2 in order).
pub fn matmul3(a: &Matrix3, b: &Matrix3) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]))
}

/// Eigen's fixed 3x3 product (lazy, SSE2): rows 0 and 1 are one packet summed
/// (k0 + k1) + k2, row 2 is a scalar redux summed k0 + (k1 + k2).
pub fn eigen_matmul3(a: &Matrix3, b: &Matrix3) -> Matrix3 {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i < 2 {
                a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]
            } else {
                a[i][0] * b[0][j] + (a[i][1] * b[1][j] + a[i][2] * b[2][j])
            }
        })
    })
}

/// Eigen's Matrix3d * Vector3d: same summation orders as eigen_matmul3.
pub fn eigen_mul3(a: &Matrix3, v: [f64; 3]) -> [f64; 3] {
    [
        a[0][0] * v[0] + a[0][1] * v[1] + a[0][2] * v[2],
        a[1][0] * v[0] + a[1][1] * v[1] + a[1][2] * v[2],
        a[2][0] * v[0] + (a[2][1] * v[1] + a[2][2] * v[2]),
    ]
}

pub fn transpose3(a: &Matrix3) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}

/// u * v^T.
pub fn outer3(u: [f64; 3], v: [f64; 3]) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| u[i] * v[j]))
}

pub fn add3(a: &Matrix3, b: &Matrix3) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}

pub fn sub3(a: &Matrix3, b: &Matrix3) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] - b[i][j]))
}

/// The Frobenius norm.
pub fn norm3(a: &Matrix3) -> f64 {
    a.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()
}

pub const IDENTITY3: Matrix3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// A plane rotation (Eigen's JacobiRotation, real case).
#[derive(Clone, Copy, Debug)]
struct Rotation {
    c: f64,
    s: f64,
}

impl Rotation {
    fn transpose(self) -> Rotation {
        Rotation { c: self.c, s: -self.s }
    }

    fn mul(self, o: Rotation) -> Rotation {
        Rotation { c: self.c * o.c - self.s * o.s, s: self.c * o.s + self.s * o.c }
    }

    /// makeJacobi(x, y, z) of the selfadjoint 2x2 matrix [[x, y], [y, z]].
    fn make_jacobi(x: f64, y: f64, z: f64) -> Rotation {
        let deno = 2.0 * y.abs();
        if deno < f64::MIN_POSITIVE {
            return Rotation { c: 1.0, s: 0.0 };
        }
        let tau = (x - z) / deno;
        let w = (tau * tau + 1.0).sqrt();
        let t = if tau > 0.0 { 1.0 / (tau + w) } else { 1.0 / (tau - w) };
        let sign_t = if t > 0.0 { 1.0 } else { -1.0 };
        let n = 1.0 / (t * t + 1.0).sqrt();
        Rotation { c: n, s: -sign_t * (y / y.abs()) * t.abs() * n }
    }
}

/// apply_rotation_in_the_plane on two vectors given by index accessors.
fn rotate_pair(m: &mut Matrix3, x: [(usize, usize); 3], y: [(usize, usize); 3], j: Rotation) {
    if j.c == 1.0 && j.s == 0.0 {
        return;
    }
    for k in 0..3 {
        let (xi, yi) = (m[x[k].0][x[k].1], m[y[k].0][y[k].1]);
        m[x[k].0][x[k].1] = j.c * xi + j.s * yi;
        m[y[k].0][y[k].1] = -j.s * xi + j.c * yi;
    }
}

/// applyOnTheLeft(p, q, j): rows p and q.
fn apply_left(m: &mut Matrix3, p: usize, q: usize, j: Rotation) {
    rotate_pair(m, [(p, 0), (p, 1), (p, 2)], [(q, 0), (q, 1), (q, 2)], j);
}

/// applyOnTheRight(p, q, j): columns p and q with j transposed.
fn apply_right(m: &mut Matrix3, p: usize, q: usize, j: Rotation) {
    rotate_pair(m, [(0, p), (1, p), (2, p)], [(0, q), (1, q), (2, q)], j.transpose());
}

/// real_2x2_jacobi_svd of the (p, q) block: (left, right) rotations.
fn real_2x2_jacobi_svd(work: &Matrix3, p: usize, q: usize) -> (Rotation, Rotation) {
    let mut m = [[work[p][p], work[p][q]], [work[q][p], work[q][q]]];
    let t = m[0][0] + m[1][1];
    let d = m[1][0] - m[0][1];
    let rot1 = if d.abs() < f64::MIN_POSITIVE {
        Rotation { c: 1.0, s: 0.0 }
    } else {
        let u = t / d;
        let tmp = (1.0 + u * u).sqrt();
        Rotation { c: u / tmp, s: 1.0 / tmp }
    };
    if !(rot1.c == 1.0 && rot1.s == 0.0) {
        for i in 0..2 {
            let (x, y) = (m[0][i], m[1][i]);
            m[0][i] = rot1.c * x + rot1.s * y;
            m[1][i] = -rot1.s * x + rot1.c * y;
        }
    }
    let right = Rotation::make_jacobi(m[0][0], m[0][1], m[1][1]);
    (rot1.mul(right.transpose()), right)
}

/// Eigen's JacobiSVD of a 3x3 matrix (two-sided Jacobi, full U and V):
/// (U, singular values descending, V) with m = U diag(s) V^T.
pub fn jacobi_svd3(matrix: &Matrix3) -> (Matrix3, [f64; 3], Matrix3) {
    let precision = 2.0 * f64::EPSILON;
    let consider_as_zero = f64::MIN_POSITIVE;
    let mut scale = 0.0f64;
    for row in matrix {
        for v in row {
            if scale < v.abs() {
                scale = v.abs();
            }
        }
    }
    if scale == 0.0 {
        scale = 1.0;
    }
    let mut work: Matrix3 = std::array::from_fn(|i| std::array::from_fn(|j| matrix[i][j] / scale));
    let (mut u, mut v) = (IDENTITY3, IDENTITY3);
    let maxi = |a: f64, b: f64| if a < b { b } else { a };
    let mut max_diag = maxi(maxi(work[0][0].abs(), work[1][1].abs()), work[2][2].abs());
    let mut finished = false;
    while !finished {
        finished = true;
        for p in 1..3 {
            for q in 0..p {
                let threshold = maxi(consider_as_zero, precision * max_diag);
                if work[p][q].abs() > threshold || work[q][p].abs() > threshold {
                    finished = false;
                    let (left, right) = real_2x2_jacobi_svd(&work, p, q);
                    apply_left(&mut work, p, q, left);
                    apply_right(&mut u, p, q, left.transpose());
                    apply_right(&mut work, p, q, right);
                    apply_right(&mut v, p, q, right);
                    max_diag = maxi(max_diag, maxi(work[p][p].abs(), work[q][q].abs()));
                }
            }
        }
    }
    let mut singular = [0.0; 3];
    for i in 0..3 {
        let a = work[i][i];
        singular[i] = a.abs();
        if a < 0.0 {
            for r in 0..3 {
                u[r][i] = -u[r][i];
            }
        }
    }
    for s in singular.iter_mut() {
        *s *= scale;
    }
    for i in 0..3 {
        let mut pos = 0;
        for k in 1..3 - i {
            if singular[i + k] > singular[i + pos] {
                pos = k;
            }
        }
        if singular[i + pos] == 0.0 {
            break;
        }
        if pos != 0 {
            let pos = pos + i;
            singular.swap(i, pos);
            for r in 0..3 {
                u[r].swap(i, pos);
                v[r].swap(i, pos);
            }
        }
    }
    (u, singular, v)
}

/// Eigen's JacobiRotation::makeGivens (real): the rotation zeroing q in (p, q).
fn make_givens(p: f64, q: f64) -> Rotation {
    if q == 0.0 {
        Rotation { c: if p < 0.0 { -1.0 } else { 1.0 }, s: 0.0 }
    } else if p == 0.0 {
        Rotation { c: 0.0, s: if q < 0.0 { 1.0 } else { -1.0 } }
    } else if p.abs() > q.abs() {
        let t = q / p;
        let mut u = (1.0 + t * t).sqrt();
        if p < 0.0 {
            u = -u;
        }
        let c = 1.0 / u;
        Rotation { c, s: -t * c }
    } else {
        let t = p / q;
        let mut u = (1.0 + t * t).sqrt();
        if q < 0.0 {
            u = -u;
        }
        let s = -1.0 / u;
        Rotation { c: -t * s, s }
    }
}

/// Eigen's numext::hypot (positive_real_hypot).
fn eigen_hypot(x: f64, y: f64) -> f64 {
    let (x, y) = (x.abs(), y.abs());
    if x.is_infinite() || y.is_infinite() {
        return f64::INFINITY;
    }
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    let p = if x < y { y } else { x };
    if p == 0.0 {
        return 0.0;
    }
    let qp = (if x < y { x } else { y }) / p;
    p * (1.0 + qp * qp).sqrt()
}

/// One implicit symmetric QR step with Wilkinson shift on the tridiagonal
/// (diag, subdiag) block start..=end, accumulating the rotations into q.
fn tridiagonal_qr_step(diag: &mut [f64; 3], subdiag: &mut [f64; 2], start: usize, end: usize, q: &mut Matrix3) {
    let td = (diag[end - 1] - diag[end]) * 0.5;
    let e = subdiag[end - 1];
    let mut mu = diag[end];
    if td == 0.0 {
        mu -= e.abs();
    } else if e != 0.0 {
        let e2 = e * e;
        let h = eigen_hypot(td, e);
        if e2 == 0.0 {
            mu -= e / ((td + if td > 0.0 { h } else { -h }) / e);
        } else {
            mu -= e2 / (td + if td > 0.0 { h } else { -h });
        }
    }
    let mut x = diag[start] - mu;
    let mut z = subdiag[start];
    let mut k = start;
    while k < end && z != 0.0 {
        let rot = make_givens(x, z);
        let (c, s) = (rot.c, rot.s);
        let sdk = s * diag[k] + c * subdiag[k];
        let dkp1 = s * subdiag[k] + c * diag[k + 1];
        diag[k] = c * (c * diag[k] - s * subdiag[k]) - s * (c * subdiag[k] - s * diag[k + 1]);
        diag[k + 1] = s * sdk + c * dkp1;
        subdiag[k] = c * sdk - s * dkp1;
        if k > start {
            subdiag[k - 1] = c * subdiag[k - 1] - s * z;
        }
        x = subdiag[k];
        if k < end - 1 {
            z = -s * subdiag[k + 1];
            subdiag[k + 1] = c * subdiag[k + 1];
        }
        apply_right(q, k, k + 1, rot);
        k += 1;
    }
}

/// Eigen's SelfAdjointEigenSolver<Matrix3d> (compute(), not computeDirect):
/// eigenvalues ascending and the eigenvectors as columns. Only the lower
/// triangle of `matrix` is read. None when the QR iteration does not converge.
pub fn self_adjoint_eigen3(matrix: &Matrix3) -> Option<([f64; 3], Matrix3)> {
    // Map the lower triangle to [-1, 1].
    let mut scale = 0.0f64;
    for i in 0..3 {
        for j in 0..=i {
            if scale < matrix[i][j].abs() {
                scale = matrix[i][j].abs();
            }
        }
    }
    if scale == 0.0 {
        scale = 1.0;
    }
    let m = |i: usize, j: usize| matrix[i][j] / scale;
    // Tridiagonalization (Eigen's 3x3 closed form).
    let mut diag = [m(0, 0), 0.0, 0.0];
    let mut subdiag = [0.0; 2];
    let mut q = IDENTITY3;
    let v1norm2 = m(2, 0) * m(2, 0);
    if v1norm2 <= f64::MIN_POSITIVE {
        diag[1] = m(1, 1);
        diag[2] = m(2, 2);
        subdiag[0] = m(1, 0);
        subdiag[1] = m(2, 1);
    } else {
        let beta = (m(1, 0) * m(1, 0) + v1norm2).sqrt();
        let inv_beta = 1.0 / beta;
        let m01 = m(1, 0) * inv_beta;
        let m02 = m(2, 0) * inv_beta;
        let qq = 2.0 * m01 * m(2, 1) + m02 * (m(2, 2) - m(1, 1));
        diag[1] = m(1, 1) + m02 * qq;
        diag[2] = m(2, 2) - m02 * qq;
        subdiag[0] = beta;
        subdiag[1] = m(2, 1) - m01 * qq;
        q = [[1.0, 0.0, 0.0], [0.0, m01, m02], [0.0, m02, -m01]];
    }
    // computeFromTridiagonal_impl.
    let (n, max_iterations) = (3usize, 30usize);
    let (mut end, mut start, mut iter) = (n - 1, 0usize, 0usize);
    let precision_inv = 1.0 / f64::EPSILON;
    while end > 0 {
        for i in start..end {
            if subdiag[i].abs() < f64::MIN_POSITIVE {
                subdiag[i] = 0.0;
            } else {
                let scaled = precision_inv * subdiag[i];
                if scaled * scaled <= diag[i].abs() + diag[i + 1].abs() {
                    subdiag[i] = 0.0;
                }
            }
        }
        while end > 0 && subdiag[end - 1] == 0.0 {
            end -= 1;
        }
        if end == 0 {
            break;
        }
        iter += 1;
        if iter > max_iterations * n {
            break;
        }
        start = end - 1;
        while start > 0 && subdiag[start - 1] != 0.0 {
            start -= 1;
        }
        tridiagonal_qr_step(&mut diag, &mut subdiag, start, end, &mut q);
    }
    if iter > max_iterations * n {
        return None;
    }
    // Sort ascending (selection, first minimum).
    for i in 0..n - 1 {
        let mut k = 0;
        for j in 1..n - i {
            if diag[i + j] < diag[i + k] {
                k = j;
            }
        }
        if k > 0 {
            diag.swap(i, i + k);
            for row in q.iter_mut() {
                row.swap(i, i + k);
            }
        }
    }
    for d in diag.iter_mut() {
        *d *= scale;
    }
    Some((diag, q))
}

/// The cell matrix with the lattice vectors as columns (cell(j, i) = vector i, component j).
pub fn cell_columns(cell_vectors: &[[f64; 3]; 3]) -> Matrix3 {
    std::array::from_fn(|j| std::array::from_fn(|i| cell_vectors[i][j]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_is_an_inverse() {
        let m = [[4.0, 1.0, 0.5], [0.1, 3.0, 0.2], [0.0, 0.2, 5.0]];
        let r = inverse3(&m);
        for i in 0..3 {
            for j in 0..3 {
                let p: f64 = (0..3).map(|k| m[i][k] * r[k][j]).sum();
                assert!((p - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
        assert!((determinant3(&m) - 59.35).abs() < 1e-9);
    }
}
