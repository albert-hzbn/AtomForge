//! Wannier90 tight-binding Hamiltonians (the C++ electronic/Wannier):
//! seedname_hr.dat reader, Bloch Hamiltonians, bands, Kubo Berry curvature
//! and Fukui-Hatsugai Chern numbers.
//!
//! Diagonalisation uses nalgebra's Hermitian solver (the C++ uses Eigen's
//! SelfAdjointEigenSolver), so bands agree to rounding; curvature and Chern
//! numbers are gauge invariant (docs/rust-port-plan.md 6a).

use glam::{DVec3, IVec3};
use nalgebra::{Complex, DMatrix, DVector};

use crate::text::{lines_of, stod_prefix, stoi_prefix, words};
use crate::volume::PI;

type C = Complex<f64>;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WannierHamiltonian {
    pub num_wann: i32,
    /// R, file order.
    pub vectors: Vec<IVec3>,
    /// H(R), row-major num_wann x num_wann, aligned with vectors.
    pub hamiltonians: Vec<Vec<C>>,
    /// Wigner-Seitz ndegen(R), aligned with vectors.
    pub degeneracies: Vec<i32>,
}

/// GCC's complex product: (ac - bd, ad + bc).
fn cmul(a: C, b: C) -> C {
    C::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

pub fn read_wannier_hr(path: &str) -> Result<WannierHamiltonian, String> {
    let lines = lines_of(path).ok_or_else(|| "Cannot open Wannier90 Hamiltonian file".to_string())?;
    let mut next = lines.iter();
    let mut line = || next.next().ok_or_else(|| "Truncated Wannier90 Hamiltonian file".to_string());
    line()?; // comment/date line, unused
    let header1 = words(line()?);
    let header2 = words(line()?);
    if header1.len() != 1 || header2.len() != 1 {
        return Err("Invalid Wannier90 Hamiltonian header".into());
    }
    let num_wann = stoi_prefix(header1[0])?;
    let num_r = stoi_prefix(header2[0])?;
    if !(1..=2000).contains(&num_wann) || !(1..=1_000_000).contains(&num_r) {
        return Err("Invalid Wannier90 Hamiltonian dimensions".into());
    }
    let mut degeneracy_values = Vec::new();
    while (degeneracy_values.len() as i32) < num_r {
        for token in words(line()?) {
            degeneracy_values.push(stoi_prefix(token)?);
        }
    }
    if degeneracy_values.len() as i32 != num_r || degeneracy_values.iter().any(|&d| d < 1) {
        return Err("Invalid Wannier90 degeneracy block".into());
    }
    let n = num_wann as usize;
    let mut result = WannierHamiltonian { num_wann, ..Default::default() };
    for degeneracy in degeneracy_values {
        let mut matrix = vec![C::new(0.0, 0.0); n * n];
        let mut r_vector = IVec3::ZERO;
        let mut first = true;
        for _ in 0..n * n {
            let tokens = words(line()?);
            if tokens.len() != 7 {
                return Err("Invalid Wannier90 Hamiltonian row".into());
            }
            let r = IVec3::new(stoi_prefix(tokens[0])?, stoi_prefix(tokens[1])?, stoi_prefix(tokens[2])?);
            let m = stoi_prefix(tokens[3])? - 1;
            let k = stoi_prefix(tokens[4])? - 1;
            if first {
                r_vector = r;
                first = false;
            } else if r != r_vector {
                return Err("Hamiltonian block mixes R vectors".into());
            }
            if m < 0 || m >= num_wann || k < 0 || k >= num_wann {
                return Err("Wannier index outside num_wann".into());
            }
            matrix[m as usize * n + k as usize] = C::new(stod_prefix(tokens[5])?, stod_prefix(tokens[6])?);
        }
        if result.vectors.contains(&r_vector) {
            return Err("Repeated R vector in Wannier90 Hamiltonian".into());
        }
        result.vectors.push(r_vector);
        result.hamiltonians.push(matrix);
        result.degeneracies.push(degeneracy);
    }
    Ok(result)
}

/// H(k), or dH(k)/dk along 0, 1 or 2 (`derivative` -1 for H(k) itself).
fn bloch_matrix(hr: &WannierHamiltonian, k: DVec3, derivative: i32) -> Result<DMatrix<C>, String> {
    if !k.is_finite() {
        return Err("kpoint must have three finite fractional components".into());
    }
    let n = hr.num_wann as usize;
    let mut result = DMatrix::<C>::zeros(n, n);
    for (i, r) in hr.vectors.iter().enumerate() {
        let dot = k.x * r.x as f64 + k.y * r.y as f64 + k.z * r.z as f64;
        let angle = 2.0 * PI * dot;
        let degeneracy = hr.degeneracies[i] as f64;
        let mut weight = C::new(angle.cos() / degeneracy, angle.sin() / degeneracy);
        if derivative >= 0 {
            let component = [r.x, r.y, r.z][derivative as usize] as f64;
            weight = cmul(weight, C::new(0.0, 2.0 * PI * component));
        }
        let h = &hr.hamiltonians[i];
        for row in 0..n {
            for col in 0..n {
                result[(row, col)] += cmul(weight, h[row * n + col]);
            }
        }
    }
    Ok(result)
}

/// Eigenvalues ascending with their eigenvectors (columns).
fn diagonalize(h: &DMatrix<C>) -> Result<(Vec<f64>, DMatrix<C>), String> {
    let mut largest = 0.0f64;
    let difference = h - h.adjoint();
    for v in difference.iter() {
        let magnitude = v.re.hypot(v.im);
        if largest < magnitude {
            largest = magnitude;
        }
    }
    if largest > 1e-8 {
        return Err("Interpolated Hamiltonian is not Hermitian; check R/-R hopping consistency".into());
    }
    let solver = h.clone().symmetric_eigen();
    let mut order: Vec<usize> = (0..solver.eigenvalues.len()).collect();
    order.sort_by(|&a, &b| solver.eigenvalues[a].partial_cmp(&solver.eigenvalues[b]).unwrap_or(std::cmp::Ordering::Equal));
    let values = order.iter().map(|&i| solver.eigenvalues[i]).collect();
    let vectors = DMatrix::from_columns(&order.iter().map(|&i| solver.eigenvectors.column(i).into_owned()).collect::<Vec<DVector<C>>>());
    Ok((values, vectors))
}

fn check_plane(a: i32, b: i32) -> Result<(), String> {
    if a == b || !(0..=2).contains(&a) || !(0..=2).contains(&b) {
        return Err("plane must select two distinct reciprocal directions".into());
    }
    Ok(())
}

/// H(k), row-major.
pub fn bloch_hamiltonian(hr: &WannierHamiltonian, k: DVec3) -> Result<Vec<C>, String> {
    let m = bloch_matrix(hr, k, -1)?;
    let n = hr.num_wann as usize;
    Ok((0..n * n).map(|i| m[(i / n, i % n)]).collect())
}

/// num_wann energies per k point (fractional), ascending.
pub fn wannier_bands(hr: &WannierHamiltonian, kpoints: &[DVec3]) -> Result<Vec<f64>, String> {
    let mut energies = Vec::with_capacity(kpoints.len() * hr.num_wann as usize);
    for &k in kpoints {
        energies.extend(diagonalize(&bloch_matrix(hr, k, -1)?)?.0);
    }
    Ok(energies)
}

/// Kubo Berry curvature of each band in the (a, b) plane at k.
pub fn berry_curvature(hr: &WannierHamiltonian, k: DVec3, plane_a: i32, plane_b: i32) -> Result<Vec<f64>, String> {
    check_plane(plane_a, plane_b)?;
    let (values, states) = diagonalize(&bloch_matrix(hr, k, -1)?)?;
    let va = states.adjoint() * bloch_matrix(hr, k, plane_a)? * &states;
    let vb = states.adjoint() * bloch_matrix(hr, k, plane_b)? * &states;
    let n = hr.num_wann as usize;
    let mut curvature = vec![0.0; n];
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            if i == j {
                continue;
            }
            let gap = values[i] - values[j];
            let contribution = 2.0 * cmul(va[(i, j)], vb[(j, i)]).im / (gap * gap);
            if !contribution.is_finite() {
                return Err(format!("Band {i} is degenerate at this k point; Kubo curvature diverges"));
            }
            sum += contribution;
        }
        curvature[i] = sum;
    }
    Ok(curvature)
}

/// Chern number of a band over a grid x grid mesh of the (a, b) plane at `fixed`.
pub fn chern_number(hr: &WannierHamiltonian, band: i32, plane_a: i32, plane_b: i32, grid: i32, fixed: f64) -> Result<f64, String> {
    check_plane(plane_a, plane_b)?;
    if grid < 2 {
        return Err("grid must be an integer >= 2".into());
    }
    if band < 0 || band >= hr.num_wann {
        return Err("band index is outside num_wann".into());
    }
    let g = grid as usize;
    let mut states = Vec::with_capacity(g * g);
    for i in 0..grid {
        for j in 0..grid {
            let mut k = [fixed; 3];
            k[plane_a as usize] = i as f64 / grid as f64;
            k[plane_b as usize] = j as f64 / grid as f64;
            let (_, vectors) = diagonalize(&bloch_matrix(hr, DVec3::from_array(k), -1)?)?;
            states.push(vectors.column(band as usize).into_owned());
        }
    }
    let mut flux = 0.0;
    for i in 0..g {
        for j in 0..g {
            let (i1, j1) = ((i + 1) % g, (j + 1) % g);
            let (s00, s10, s11, s01) = (&states[i * g + j], &states[i1 * g + j], &states[i1 * g + j1], &states[i * g + j1]);
            let mut product = C::new(1.0, 0.0);
            for link in [s00.dotc(s10), s10.dotc(s11), s11.dotc(s01), s01.dotc(s00)] {
                let magnitude = link.re.hypot(link.im);
                if magnitude < 1e-10 {
                    return Err(format!("Band {band} touches a degeneracy on this mesh; refine grid or shift fixed"));
                }
                product = cmul(product, C::new(link.re / magnitude, link.im / magnitude));
            }
            flux += product.im.atan2(product.re);
        }
    }
    Ok(flux / (2.0 * PI))
}
