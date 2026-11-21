//! Fourth-order elastic tensors (the C++ algorithms/ElasticConstants): from
//! Voigt matrices or cubic/hexagonal constants, rotated to new axes, and the
//! isotropy test.

use glam::DMat3;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ElasticTensor {
    pub c: [[[[f64; 3]; 3]; 3]; 3],
}

/// Voigt index: xx=0, yy=1, zz=2, yz=3, xz=4, xy=5.
fn voigt_index(i: usize, j: usize) -> usize {
    if i == j {
        return i;
    }
    match (i.min(j), i.max(j)) {
        (1, 2) => 3,
        (0, 2) => 4,
        _ => 5,
    }
}

/// A tensor from a Voigt matrix (either triangle may be given).
pub fn make_voigt_elastic_tensor(voigt: &[[f64; 6]; 6]) -> ElasticTensor {
    let mut t = ElasticTensor::default();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let (a, b) = (voigt_index(i, j), voigt_index(k, l));
                    t.c[i][j][k][l] = if voigt[a][b] != 0.0 { voigt[a][b] } else { voigt[b][a] };
                }
            }
        }
    }
    t
}

pub fn make_cubic_elastic_tensor(c11: f64, c12: f64, c44: f64) -> ElasticTensor {
    let mut v = [[0.0; 6]; 6];
    for i in 0..3 {
        v[i][i] = c11;
        v[i + 3][i + 3] = c44;
        for j in 0..3 {
            if i != j {
                v[i][j] = c12;
            }
        }
    }
    make_voigt_elastic_tensor(&v)
}

pub fn make_hexagonal_elastic_tensor(c11: f64, c12: f64, c13: f64, c33: f64, c44: f64) -> ElasticTensor {
    let mut v = [[0.0; 6]; 6];
    v[0][0] = c11;
    v[1][1] = c11;
    v[2][2] = c33;
    v[0][1] = c12;
    v[1][0] = c12;
    for (a, b) in [(0, 2), (2, 0), (1, 2), (2, 1)] {
        v[a][b] = c13;
    }
    v[3][3] = c44;
    v[4][4] = c44;
    v[5][5] = 0.5 * (c11 - c12); // C66 of hexagonal symmetry
    make_voigt_elastic_tensor(&v)
}

/// C'_ijkl = R_pi R_qj R_rk R_sl C_pqrs, the new axes as the columns of `axes`.
pub fn rotate_elastic_tensor(c: &ElasticTensor, axes: &DMat3) -> ElasticTensor {
    let mut r = [[0.0f64; 3]; 3];
    for p in 0..3 {
        for i in 0..3 {
            r[p][i] = axes.col(i)[p];
        }
    }
    let mut out = ElasticTensor::default();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let mut sum = 0.0;
                    for p in 0..3 {
                        for q in 0..3 {
                            for rr in 0..3 {
                                for s in 0..3 {
                                    sum += r[p][i] * r[q][j] * r[rr][k] * r[s][l] * c.c[p][q][rr][s];
                                }
                            }
                        }
                    }
                    out.c[i][j][k][l] = sum;
                }
            }
        }
    }
    out
}

/// Whether the best isotropic fit leaves a relative residual within the tolerance.
pub fn is_nearly_isotropic(c: &ElasticTensor, relative_tolerance: f64) -> bool {
    let mut trace = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            trace += c.c[i][i][j][j];
        }
    }
    let lambda_plus = trace / 9.0;
    let (mut shear_sum, mut shear_count) = (0.0, 0);
    for i in 0..3 {
        for j in 0..3 {
            if i != j {
                shear_sum += c.c[i][j][i][j];
                shear_count += 1;
            }
        }
    }
    let mu = shear_sum / shear_count as f64;
    let lambda = lambda_plus - (2.0 / 3.0) * mu;
    let (mut residual, mut norm) = (0.0, 0.0);
    let d = |a: usize, b: usize| if a == b { 1.0 } else { 0.0 };
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let iso = lambda * d(i, j) * d(k, l) + mu * (d(i, k) * d(j, l) + d(i, l) * d(j, k));
                    let diff = c.c[i][j][k][l] - iso;
                    residual += diff * diff;
                    norm += c.c[i][j][k][l] * c.c[i][j][k][l];
                }
            }
        }
    }
    if norm <= 1e-12 {
        return true;
    }
    (residual / norm).sqrt() <= relative_tolerance
}
