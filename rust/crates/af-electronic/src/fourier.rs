//! Fourier analysis of periodic grids (the C++ electronic/FourierAnalysis and
//! Electrostatics): structure factors, Fourier synthesis, Patterson maps,
//! atomic structure factors and Ewald site potentials. Complex arithmetic
//! follows GCC's inline expansion of std::complex<double>.

use af_core::glm::{determinant_d, inverse_d, mul_d};
use glam::{DMat3, DVec3, IVec3};

use crate::volume::{Grid, BOHR, HARTREE, PI};

extern "C" {
    #[link_name = "erfc"]
    fn c_erfc(x: f64) -> f64;
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub const ZERO: Complex = Complex { re: 0.0, im: 0.0 };
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    /// std::polar(1.0, theta).
    pub fn polar(theta: f64) -> Self {
        Self { re: theta.cos(), im: theta.sin() }
    }
    pub fn mul(self, o: Complex) -> Complex {
        Complex { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
    pub fn scale(self, s: f64) -> Complex {
        Complex { re: self.re * s, im: self.im * s }
    }
    pub fn add(self, o: Complex) -> Complex {
        Complex { re: self.re + o.re, im: self.im + o.im }
    }
    pub fn sub(self, o: Complex) -> Complex {
        Complex { re: self.re - o.re, im: self.im - o.im }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Reflection {
    pub hkl: IVec3,
    pub value: Complex,
}

fn fft(a: &mut [Complex], sign: f64) {
    let n = a.len();
    if n & (n - 1) != 0 {
        // Non-power-of-two axes use an exact DFT; the workload is checked by the caller.
        let original = a.to_vec();
        for k in 0..n {
            let mut sum = Complex::ZERO;
            for (j, value) in original.iter().enumerate() {
                sum = sum.add(value.mul(Complex::polar(sign * 2.0 * PI * k as f64 * j as f64 / n as f64)));
            }
            a[k] = sum;
        }
        return;
    }
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut length = 2;
    while length <= n {
        let root = Complex::polar(sign * 2.0 * PI / length as f64);
        let mut start = 0;
        while start < n {
            let mut w = Complex::new(1.0, 0.0);
            for j in 0..length / 2 {
                let even = a[start + j];
                let odd = a[start + j + length / 2].mul(w);
                a[start + j] = even.add(odd);
                a[start + j + length / 2] = even.sub(odd);
                w = w.mul(root);
            }
            start += length;
        }
        length <<= 1;
    }
}

fn transform(grid: &Grid, data: &mut [Complex], sign: f64) -> Result<(), String> {
    let mut work = 0.0;
    for &n in &grid.shape {
        work += (if n & (n - 1) != 0 { n as f64 } else { (n as f64).log2() }) * data.len() as f64;
    }
    if work > 2e8 {
        return Err("Fourier workload too large; resample to power-of-two dimensions".into());
    }
    for axis in 0..3 {
        let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
        let mut row = vec![Complex::ZERO; grid.shape[axis] as usize];
        for j in 0..grid.shape[b] {
            for k in 0..grid.shape[c] {
                let mut p = IVec3::ZERO;
                p[b] = j;
                p[c] = k;
                for i in 0..grid.shape[axis] {
                    p[axis] = i;
                    row[i as usize] = data[grid.index(p.x, p.y, p.z)];
                }
                fft(&mut row, sign);
                for i in 0..grid.shape[axis] {
                    p[axis] = i;
                    data[grid.index(p.x, p.y, p.z)] = row[i as usize];
                }
            }
        }
    }
    Ok(())
}

fn require_periodic(grid: &Grid) -> Result<(), String> {
    grid.validate()?;
    if !grid.periodic {
        return Err("Fourier tools require periodic endpoint-excluded grids; convert explicitly".into());
    }
    Ok(())
}

fn as_dvec(h: IVec3) -> DVec3 {
    DVec3::new(h.x as f64, h.y as f64, h.z as f64)
}

/// The grid's position in fractional cell coordinates.
fn fractional_origin(grid: &Grid) -> DVec3 {
    mul_d(&inverse_d(&grid.cell), grid.origin)
}

fn unique_index(grid: &Grid, h: IVec3) -> Result<IVec3, String> {
    let mut q = h;
    for a in 0..3 {
        if h[a] < -grid.shape[a] / 2 || h[a] > (grid.shape[a] - 1) / 2 {
            return Err("Reflection outside unique sampled Nyquist range".into());
        }
        q[a] = (h[a] + grid.shape[a]) % grid.shape[a];
    }
    Ok(q)
}

pub fn structure_factors(grid: &Grid, indices: &[IVec3]) -> Result<Vec<Reflection>, String> {
    require_periodic(grid)?;
    let mut data: Vec<Complex> = grid.values.iter().map(|&v| Complex::new(v, 0.0)).collect();
    transform(grid, &mut data, 1.0)?;
    let voxel = determinant_d(&grid.cell).abs() / grid.values.len() as f64;
    let origin = fractional_origin(grid);
    let mut out = Vec::new();
    for &h in indices {
        let q = unique_index(grid, h)?;
        let value = data[grid.index(q.x, q.y, q.z)].scale(voxel).mul(Complex::polar(2.0 * PI * as_dvec(h).dot(origin)));
        out.push(Reflection { hkl: h, value });
    }
    Ok(out)
}

pub fn fourier_synthesis(geometry: &Grid, reflections: &[Reflection]) -> Result<Grid, String> {
    require_periodic(geometry)?;
    let mut data = vec![Complex::ZERO; geometry.values.len()];
    let mut occupied = vec![false; data.len()];
    let origin = fractional_origin(geometry);
    for r in reflections {
        let q = unique_index(geometry, r.hkl)?;
        let i = geometry.index(q.x, q.y, q.z);
        if occupied[i] {
            return Err("Duplicate reflection".into());
        }
        if !r.value.re.is_finite() || !r.value.im.is_finite() {
            return Err("Non-finite structure factor".into());
        }
        occupied[i] = true;
        data[i] = r.value.mul(Complex::polar(-2.0 * PI * as_dvec(r.hkl).dot(origin)));
    }
    transform(geometry, &mut data, -1.0)?;
    let mut out = geometry.clone();
    out.name = "Fourier synthesis".into();
    let volume = determinant_d(&out.cell).abs();
    let (mut real_max, mut imaginary_max) = (0.0f64, 0.0f64);
    for v in &data {
        if real_max < v.re.abs() {
            real_max = v.re.abs();
        }
        if imaginary_max < v.im.abs() {
            imaginary_max = v.im.abs();
        }
    }
    if imaginary_max > 1e-8 * if 1.0 < real_max { real_max } else { 1.0 } {
        return Err("Structure factors do not produce a real density; supply conjugate/Friedel pairs explicitly".into());
    }
    for (o, v) in out.values.iter_mut().zip(&data) {
        *o = v.re / volume;
    }
    out.validate()?;
    Ok(out)
}

pub fn patterson(grid: &Grid) -> Result<Grid, String> {
    require_periodic(grid)?;
    let mut data: Vec<Complex> = grid.values.iter().map(|&v| Complex::new(v, 0.0)).collect();
    transform(grid, &mut data, 1.0)?;
    for v in data.iter_mut() {
        *v = Complex::new(v.re * v.re + v.im * v.im, 0.0);
    }
    transform(grid, &mut data, -1.0)?;
    let mut out = grid.clone();
    out.origin = DVec3::ZERO;
    out.name = "Patterson autocorrelation".into();
    out.unit = format!("({})^2*A^3", grid.unit);
    let factor = determinant_d(&grid.cell).abs() / (data.len() as f64 * data.len() as f64);
    for (o, v) in out.values.iter_mut().zip(&data) {
        *o = v.re * factor;
    }
    out.validate()?;
    Ok(out)
}

/// Structure factors from explicit atomic scattering factors (one per
/// reflection and site), occupancies and isotropic B (A^2).
pub fn atomic_structure_factors(
    geometry: &Grid,
    positions: &[DVec3],
    occupancy: &[f64],
    b_factors: &[f64],
    indices: &[IVec3],
    scattering: &[Complex],
) -> Result<Vec<Reflection>, String> {
    require_periodic(geometry)?;
    let count = positions.len();
    if occupancy.len() != count || b_factors.len() != count || scattering.len() != indices.len() * count {
        return Err("Atomic factors require one occupancy/B per site and one scattering factor per reflection/site".into());
    }
    let inverse = inverse_d(&geometry.cell);
    let transposed = inverse.transpose();
    let mut out = Vec::new();
    for (r, &hkl) in indices.iter().enumerate() {
        let h = as_dvec(hkl);
        let reciprocal = mul_d(&transposed, h);
        let s2 = reciprocal.dot(reciprocal) / 4.0;
        let mut sum = Complex::ZERO;
        for j in 0..count {
            if !occupancy[j].is_finite() || occupancy[j] < 0.0 || occupancy[j] > 1.0 || !b_factors[j].is_finite() || b_factors[j] < 0.0 {
                return Err("Occupancy must be 0..1 and B must be nonnegative".into());
            }
            let factor = scattering[r * count + j];
            let p = mul_d(&inverse, positions[j]);
            if !p.is_finite() {
                return Err("Non-finite atomic position".into());
            }
            if !factor.re.is_finite() || !factor.im.is_finite() {
                return Err("Non-finite scattering factor".into());
            }
            let term = factor.scale(occupancy[j]).scale((-b_factors[j] * s2).exp()).mul(Complex::polar(2.0 * PI * h.dot(p)));
            sum = sum.add(term);
        }
        out.push(Reflection { hkl, value: sum });
    }
    Ok(out)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct EwaldResult {
    /// Site potentials (V).
    pub potentials: Vec<f64>,
    /// eV.
    pub energy: f64,
}

/// Ewald site potentials of point charges in a neutral periodic cell
/// (conducting boundary), cell vectors as columns.
pub fn ewald(cell: &DMat3, sites: &[DVec3], charges: &[f64], alpha: f64, real_cutoff: f64, reciprocal_cutoff: f64) -> Result<EwaldResult, String> {
    let check = Grid { cell: *cell, values: vec![0.0; 8], ..Default::default() };
    check.validate()?;
    if sites.is_empty() || sites.len() != charges.len() {
        return Err("Ewald requires one explicit charge per site".into());
    }
    if !alpha.is_finite() || !real_cutoff.is_finite() || !reciprocal_cutoff.is_finite() || alpha <= 0.0 || real_cutoff <= 0.0 || reciprocal_cutoff <= 0.0 {
        return Err("Ewald alpha and both cutoffs must be positive and finite".into());
    }
    let (mut total, mut magnitude) = (0.0f64, 0.0f64);
    for i in 0..sites.len() {
        if !charges[i].is_finite() {
            return Err("Non-finite charge".into());
        }
        if !sites[i].is_finite() {
            return Err("Non-finite site position".into());
        }
        total += charges[i];
        magnitude += charges[i].abs();
    }
    if total.abs() > 1e-10 * if 1.0 < magnitude { magnitude } else { 1.0 } {
        return Err("Ewald requires a neutral cell; no implicit compensating background is used".into());
    }
    let inverse = inverse_d(cell);
    let reciprocal = inverse.transpose() * (2.0 * PI);
    let volume = determinant_d(cell).abs();
    let (mut real_bounds, mut reciprocal_bounds) = (IVec3::ZERO, IVec3::ZERO);
    let mut real_work = sites.len() as f64 * sites.len() as f64;
    let mut reciprocal_work = sites.len() as f64;
    for a in 0..3 {
        let rb = (real_cutoff * DVec3::new(inverse.col(0)[a], inverse.col(1)[a], inverse.col(2)[a]).length()).ceil() + 1.0;
        let kb = (reciprocal_cutoff * cell.col(a).length() / (2.0 * PI)).ceil();
        if rb > 1000.0 || kb > 1000.0 {
            return Err("Excessive Ewald image bounds".into());
        }
        real_bounds[a] = rb as i32;
        reciprocal_bounds[a] = kb as i32;
        real_work *= 2.0 * rb + 1.0;
        reciprocal_work *= 2.0 * kb + 1.0;
    }
    if real_work + reciprocal_work > 2e8 {
        return Err("Ewald workload too large; reduce cutoffs or use a dedicated solver".into());
    }
    let positions: Vec<DVec3> = sites
        .iter()
        .map(|&site| {
            let f = mul_d(&inverse, site);
            mul_d(cell, f - DVec3::new(f.x.floor(), f.y.floor(), f.z.floor()))
        })
        .collect();
    let mut out = EwaldResult { potentials: vec![0.0; sites.len()], energy: 0.0 };
    for i in 0..sites.len() {
        for j in 0..sites.len() {
            for z in -real_bounds.z..=real_bounds.z {
                for y in -real_bounds.y..=real_bounds.y {
                    for x in -real_bounds.x..=real_bounds.x {
                        if i == j && x == 0 && y == 0 && z == 0 {
                            continue;
                        }
                        let r = (positions[i] - positions[j] + mul_d(cell, DVec3::new(x as f64, y as f64, z as f64))).length();
                        if r < 1e-10 {
                            return Err("Coincident periodic Ewald sites".into());
                        }
                        if r <= real_cutoff {
                            out.potentials[i] += charges[j] * unsafe { c_erfc(alpha * r) } / r;
                        }
                    }
                }
            }
        }
        out.potentials[i] -= 2.0 * alpha * charges[i] / PI.sqrt();
    }
    for z in -reciprocal_bounds.z..=reciprocal_bounds.z {
        for y in -reciprocal_bounds.y..=reciprocal_bounds.y {
            for x in -reciprocal_bounds.x..=reciprocal_bounds.x {
                if x == 0 && y == 0 && z == 0 {
                    continue;
                }
                let k = mul_d(&reciprocal, DVec3::new(x as f64, y as f64, z as f64));
                let k2 = k.dot(k);
                if k2 > reciprocal_cutoff * reciprocal_cutoff {
                    continue;
                }
                let mut rho = Complex::ZERO;
                for j in 0..sites.len() {
                    rho = rho.add(Complex::polar(k.dot(positions[j])).scale(charges[j]));
                }
                let factor = 4.0 * PI / volume * (-k2 / (4.0 * alpha * alpha)).exp() / k2;
                for i in 0..sites.len() {
                    out.potentials[i] += factor * rho.mul(Complex::polar(-k.dot(positions[i]))).re;
                }
            }
        }
    }
    // Conducting (tin-foil) boundary: no dipole surface term. Potentials are V.
    let coulomb = HARTREE * BOHR;
    for i in 0..sites.len() {
        out.potentials[i] *= coulomb;
        out.energy += 0.5 * charges[i] * out.potentials[i];
    }
    Ok(out)
}
