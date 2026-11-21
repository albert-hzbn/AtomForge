//! Straight-dislocation displacement in a fully anisotropic medium (the C++
//! algorithms/AnisotropicDislocation): the Stroh sextic eigenproblem
//! (Stroh 1958/1962; Barnett & Lothe 1973; Hirth & Lothe ch. 13), with the
//! coefficients D fixed by the Burgers vector and single-valued stresses.
//!
//! The C++ uses Eigen's EigenSolver and FullPivLU; here nalgebra's Schur
//! eigenvalues with inverse-iteration eigenvectors and full-pivot LU. The
//! displacement is independent of root order and eigenvector scaling, so it
//! agrees with the C++ to rounding (tested to 1e-9 relative).

use af_core::cxx_random::{uniform_real, Mt19937};
use af_core::eigen;
use glam::DVec3;
use nalgebra::{Matrix6, SMatrix, Vector6};
use num_complex::Complex64;

use crate::elastic::{is_nearly_isotropic, ElasticTensor};

const PI: f64 = 3.14159265358979323846;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StrohSolution {
    pub valid: bool,
    pub error: String,
    /// Roots with Im > 0 and their displacement (A) and traction (L) vectors.
    pub p: [Complex64; 3],
    pub a: [[Complex64; 3]; 3],
    pub l: [[Complex64; 3]; 3],
    pub d: [Complex64; 3],
}

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

/// Perturbs the 21 independent Voigt entries (BABEL's CVoigt_noise).
fn add_voigt_noise(c: &ElasticTensor, seed: u32, amplitude: f64) -> ElasticTensor {
    if amplitude <= 0.0 {
        return *c;
    }
    let mut norm = 0.0;
    for x in c.c.iter().flatten().flatten().flatten() {
        norm += x * x;
    }
    let scale = amplitude * norm.sqrt();
    let mut rng = Mt19937::new(seed);
    let mut delta = [[0.0f64; 6]; 6];
    for a in 0..6 {
        for b in a..6 {
            let v = uniform_real(&mut rng, -1.0, 1.0) * scale;
            delta[a][b] = v;
            delta[b][a] = v;
        }
    }
    let mut out = *c;
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    out.c[i][j][k][l] += delta[voigt_index(i, j)][voigt_index(k, l)];
                }
            }
        }
    }
    out
}

fn matmul(a: &eigen::Matrix3, b: &eigen::Matrix3) -> eigen::Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}

fn transpose(a: &eigen::Matrix3) -> eigen::Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}

/// An eigenvector of `n` for the eigenvalue `lambda` by inverse iteration.
fn eigenvector(n: &Matrix6<f64>, lambda: Complex64) -> Vector6<Complex64> {
    let shift = lambda + Complex64::new(1e-10, 1e-10) * (1.0 + lambda.norm());
    let m: SMatrix<Complex64, 6, 6> = n.map(|x| Complex64::new(x, 0.0)) - SMatrix::<Complex64, 6, 6>::identity() * shift;
    let lu = m.lu();
    let mut x = Vector6::from_fn(|i, _| Complex64::new(1.0 + 0.1 * i as f64, 0.05 * i as f64));
    for _ in 0..4 {
        if let Some(y) = lu.solve(&x) {
            let norm = y.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
            if norm > 0.0 && norm.is_finite() {
                x = y / Complex64::new(norm, 0.0);
            }
        }
    }
    x
}

/// Solves the sextic problem for constants and Burgers vector in the
/// dislocation's frame (x = e1, y = e2, z = line). Nearly isotropic input is
/// refused unless `noise_amplitude` perturbs it.
pub fn solve_stroh(elastic_local_in: &ElasticTensor, burgers_local: DVec3, noise_seed: u32, noise_amplitude: f64) -> StrohSolution {
    let mut sol = StrohSolution::default();
    if noise_amplitude <= 0.0 && is_nearly_isotropic(elastic_local_in, 1e-6) {
        sol.error = "Elastic constants are (nearly) isotropic; the Stroh sextic formalism is singular there (repeated roots). Use the isotropic model instead, or pass a small noiseAmplitude (e.g. 1e-4) to perturb the constants, matching BABEL's own documented CVoigt_noise workaround for this case.".into();
        return sol;
    }
    let c = add_voigt_noise(elastic_local_in, noise_seed, noise_amplitude);
    let (mut q, mut r, mut t) = ([[0.0; 3]; 3], [[0.0; 3]; 3], [[0.0; 3]; 3]);
    for i in 0..3 {
        for k in 0..3 {
            q[i][k] = c.c[i][0][k][0];
            r[i][k] = c.c[i][0][k][1];
            t[i][k] = c.c[i][1][k][1];
        }
    }
    let max_abs = t.iter().flatten().fold(0.0f64, |m, &x| if m < x.abs() { x.abs() } else { m });
    if eigen::determinant3(&t).abs() < 1e-9 * max_abs.powi(3) {
        sol.error = "Elastic tensor is singular in the plane perpendicular to the line direction.".into();
        return sol;
    }
    let t_inv = eigen::inverse3(&t);
    let rt = transpose(&r);
    let n1 = matmul(&t_inv, &rt).map(|row| row.map(|x| -x));
    let n2 = t_inv;
    let mut n3 = matmul(&matmul(&r, &t_inv), &rt);
    for i in 0..3 {
        for j in 0..3 {
            n3[i][j] -= q[i][j];
        }
    }
    let n = Matrix6::from_fn(|i, j| match (i < 3, j < 3) {
        (true, true) => n1[i][j],
        (true, false) => n2[i][j - 3],
        (false, true) => n3[i - 3][j],
        (false, false) => n1[j - 3][i - 3],
    });
    let values = n.complex_eigenvalues();
    let mut found = 0;
    for value in values.iter() {
        if found == 3 {
            break;
        }
        if value.im > 1e-8 {
            let v = eigenvector(&n, *value);
            sol.p[found] = *value;
            for comp in 0..3 {
                sol.a[found][comp] = v[comp];
                sol.l[found][comp] = v[comp + 3];
            }
            found += 1;
        }
    }
    if found != 3 {
        sol.error = "Sextic roots are degenerate (elastically isotropic or a high-symmetry orientation); retry with a small noiseAmplitude to break the degeneracy.".into();
        return sol;
    }
    // b = 2 Re[A D] and 0 = 2 Re[L D] as one real system for Re(D), Im(D).
    let system = Matrix6::from_fn(|i, j| {
        let (row, alpha) = (i % 3, j % 3);
        let entry = if i < 3 { sol.a[alpha][row] } else { sol.l[alpha][row] };
        if j < 3 {
            entry.re
        } else {
            -entry.im
        }
    });
    let rhs = Vector6::new(0.5 * burgers_local.x, 0.5 * burgers_local.y, 0.5 * burgers_local.z, 0.0, 0.0, 0.0);
    let Some(solution) = system.full_piv_lu().solve(&rhs) else {
        sol.error = "Burgers-vector linear system was singular or ill-conditioned.".into();
        return sol;
    };
    let residual = (system * solution - rhs).norm();
    if !(residual < 1e-6 * rhs.norm().max(1.0)) {
        sol.error = "Burgers-vector linear system was singular or ill-conditioned.".into();
        return sol;
    }
    for a in 0..3 {
        sol.d[a] = Complex64::new(solution[a], solution[a + 3]);
    }
    sol.valid = true;
    sol
}

/// Displacement u(x, y) in the solution's frame (z component zero for a
/// pure Volterra dislocation up to rounding).
pub fn anisotropic_displacement(solution: &StrohSolution, x: f64, y: f64) -> DVec3 {
    if !solution.valid {
        return DVec3::ZERO;
    }
    // Floor |z| at a small fraction of the in-plane scale so near-degenerate
    // roots stay finite away from the core.
    let floor = 1e-6 * x.abs().max(y.abs()).max(1.0);
    let mut u = DVec3::ZERO;
    for a in 0..3 {
        let mut z = Complex64::new(x, 0.0) + solution.p[a] * y;
        if z.norm() < floor {
            z = z / z.norm().max(1e-300) * floor;
        }
        let log = z.ln();
        u.x += (solution.a[a][0] * solution.d[a] * log).im;
        u.y += (solution.a[a][1] * solution.d[a] * log).im;
        u.z += (solution.a[a][2] * solution.d[a] * log).im;
    }
    u / PI
}
