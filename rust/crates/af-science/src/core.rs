//! Numerical basics of the science tools (the C++ science/ScienceCore):
//! row-vector cell algebra (cell[i] is the i-th lattice vector), the general
//! minimum image, the periodic neighbour list, FFTs and numeric arrays from
//! JSON. Formulas follow the C++ term by term so results are identical.

use af_core::json::Json;
use af_core::task::checkpoint;
use num_complex::Complex64;

pub type Vec3 = [f64; 3];
pub type Mat3 = [[f64; 3]; 3];
pub type Pbc = [bool; 3];

pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

pub fn determinant(m: &Mat3) -> f64 {
    dot(m[0], cross(m[1], m[2]))
}

pub fn inverse(m: &Mat3) -> Result<Mat3, String> {
    let det = determinant(m);
    if det == 0.0 {
        return Err("Singular 3x3 matrix".into());
    }
    Ok([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det,
        ],
    ])
}

pub fn multiply(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

pub fn transpose(m: &Mat3) -> Mat3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = m[j][i];
        }
    }
    r
}

/// Row vector times matrix: v M.
pub fn row_times(v: Vec3, m: &Mat3) -> Vec3 {
    [
        v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
        v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
        v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
    ]
}

pub fn fractional(cartesian: Vec3, cell: &Mat3) -> Result<Vec3, String> {
    Ok(row_times(cartesian, &inverse(cell)?))
}

pub fn cell_volume(cell: &Mat3) -> f64 {
    determinant(cell).abs()
}

pub fn identity() -> Mat3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// The shortest image along periodic axes (Gauss-reduced basis, one shell).
#[derive(Clone, Debug, Default)]
pub struct MinimumImage {
    active: bool,
    basis: Mat3,
    inverse: Mat3,
    pbc: Pbc,
}

impl MinimumImage {
    pub fn new(cell: &Mat3, pbc: Pbc) -> Result<Self, String> {
        let active = pbc.iter().any(|&p| p);
        if !active {
            return Ok(Self { active, pbc, ..Default::default() });
        }
        let mut basis = *cell;
        for _ in 0..100 {
            let mut changed = false;
            for a in 0..3 {
                for b in 0..3 {
                    if a == b || !pbc[a] || !pbc[b] {
                        continue;
                    }
                    let length = dot(basis[a], basis[a]);
                    let factor = (dot(basis[a], basis[b]) / length).round();
                    let reduced = sub(basis[b], scale(basis[a], factor));
                    if factor != 0.0 && dot(basis[b], basis[b]) > dot(reduced, reduced) + 1e-12 {
                        basis[b] = reduced;
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        let inverse = inverse(&basis)?;
        Ok(Self { active, basis, inverse, pbc })
    }

    pub fn active(&self) -> bool {
        self.active
    }

    pub fn apply(&self, vector: Vec3) -> Vec3 {
        if !self.active {
            return vector;
        }
        let mut f = row_times(vector, &self.inverse);
        for k in 0..3 {
            if self.pbc[k] {
                f[k] -= f[k].round();
            }
        }
        let base = row_times(f, &self.basis);
        let mut best = base;
        let mut best_norm = dot(base, base);
        let range = self.pbc.map(|p| if p { 1 } else { 0 });
        for i in -range[0]..=range[0] {
            for j in -range[1]..=range[1] {
                for k in -range[2]..=range[2] {
                    if i == 0 && j == 0 && k == 0 {
                        continue;
                    }
                    let candidate = add(base, row_times([i as f64, j as f64, k as f64], &self.basis));
                    let n = dot(candidate, candidate);
                    if n < best_norm - 1e-12 {
                        best = candidate;
                        best_norm = n;
                    }
                }
            }
        }
        best
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Neighbor {
    pub i: usize,
    pub j: usize,
    pub shift: [i32; 3],
    pub vector: Vec3,
}

/// All pairs i -> j within the cutoff, including periodic images (no self
/// pair at zero shift); vector = r_j + shift*cell - r_i of the given positions.
pub fn neighbor_list(positions: &[Vec3], cell: &Mat3, pbc: Pbc, cutoff: f64) -> Result<Vec<Neighbor>, String> {
    let mut result = Vec::new();
    if positions.is_empty() {
        return Ok(result);
    }
    let periodic = pbc.iter().any(|&p| p);
    let inv = if periodic { inverse(cell)? } else { identity() };
    let mut wrapped = positions.to_vec();
    let mut offsets = vec![[0i32; 3]; positions.len()];
    if periodic {
        for i in 0..positions.len() {
            let f = row_times(positions[i], &inv);
            let mut shift = [0.0; 3];
            for k in 0..3 {
                if pbc[k] {
                    offsets[i][k] = f[k].floor() as i32;
                    shift[k] = offsets[i][k] as f64;
                }
            }
            wrapped[i] = sub(positions[i], row_times(shift, cell));
        }
    }
    let mut range = [0i32; 3];
    for k in 0..3 {
        if pbc[k] {
            let reciprocal = [inv[0][k], inv[1][k], inv[2][k]];
            range[k] = (cutoff * norm(reciprocal)).ceil() as i32;
            if range[k] > 200 {
                return Err("Neighbor cutoff is too large for the periodic cell".into());
            }
        }
    }
    let mut low = wrapped[0];
    let mut high = wrapped[0];
    for p in &wrapped {
        for k in 0..3 {
            low[k] = low[k].min(p[k]);
            high[k] = high[k].max(p[k]);
        }
    }
    struct Image {
        atom: usize,
        shift: [i32; 3],
        position: Vec3,
    }
    let mut images = Vec::new();
    for a in -range[0]..=range[0] {
        for b in -range[1]..=range[1] {
            for c in -range[2]..=range[2] {
                let translation = row_times([a as f64, b as f64, c as f64], cell);
                for (j, w) in wrapped.iter().enumerate() {
                    let p = add(*w, translation);
                    if (0..3).all(|k| p[k] >= low[k] - cutoff && p[k] <= high[k] + cutoff) {
                        images.push(Image { atom: j, shift: [a, b, c], position: p });
                    }
                }
                checkpoint()?;
            }
        }
    }
    let key = |p: Vec3, dx: i64, dy: i64, dz: i64| -> i64 {
        let x = ((p[0] - low[0] + cutoff) / cutoff).floor() as i64 + dx;
        let y = ((p[1] - low[1] + cutoff) / cutoff).floor() as i64 + dy;
        let z = ((p[2] - low[2] + cutoff) / cutoff).floor() as i64 + dz;
        (x.wrapping_mul(2000003).wrapping_add(y)).wrapping_mul(2000003).wrapping_add(z)
    };
    let mut bins: std::collections::HashMap<i64, Vec<usize>> = std::collections::HashMap::new();
    for (n, image) in images.iter().enumerate() {
        bins.entry(key(image.position, 0, 0, 0)).or_default().push(n);
    }
    let cutoff2 = cutoff * cutoff;
    for i in 0..wrapped.len() {
        if i % 256 == 0 {
            checkpoint()?;
        }
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(found) = bins.get(&key(wrapped[i], dx, dy, dz)) else { continue };
                    for &index in found {
                        let image = &images[index];
                        if image.atom == i && image.shift == [0, 0, 0] {
                            continue;
                        }
                        let d = sub(image.position, wrapped[i]);
                        if dot(d, d) >= cutoff2 {
                            continue;
                        }
                        let mut shift = [0i32; 3];
                        for k in 0..3 {
                            shift[k] = image.shift[k] + offsets[image.atom][k] - offsets[i][k];
                        }
                        result.push(Neighbor { i, j: image.atom, shift, vector: d });
                    }
                }
            }
        }
    }
    Ok(result)
}

fn radix2(data: &mut [Complex64], inverse: bool) {
    let n = data.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            data.swap(i, j);
        }
    }
    let pi = (-1.0f64).acos();
    let mut length = 2;
    while length <= n {
        let angle = 2.0 * pi / length as f64 * if inverse { 1.0 } else { -1.0 };
        for start in (0..n).step_by(length) {
            for k in 0..length / 2 {
                let w = Complex64::from_polar(1.0, angle * k as f64);
                let u = data[start + k];
                let v = data[start + k + length / 2] * w;
                data[start + k] = u + v;
                data[start + k + length / 2] = u - v;
            }
        }
        length <<= 1;
    }
}

/// Discrete Fourier transform of any length (Bluestein beyond powers of two).
pub fn fft(mut data: Vec<Complex64>, inverse: bool) -> Vec<Complex64> {
    let n = data.len();
    if n <= 1 {
        return data;
    }
    if n & (n - 1) == 0 {
        radix2(&mut data, inverse);
        if inverse {
            for v in data.iter_mut() {
                *v /= n as f64;
            }
        }
        return data;
    }
    let mut m = 1;
    while m < 2 * n - 1 {
        m <<= 1;
    }
    let pi = (-1.0f64).acos();
    let chirp: Vec<Complex64> = (0..n)
        .map(|k| {
            let phase = pi * ((k * k) % (2 * n)) as f64 / n as f64 * if inverse { 1.0 } else { -1.0 };
            Complex64::from_polar(1.0, phase)
        })
        .collect();
    let mut a = vec![Complex64::new(0.0, 0.0); m];
    let mut b = vec![Complex64::new(0.0, 0.0); m];
    for k in 0..n {
        a[k] = data[k] * chirp[k];
    }
    b[0] = chirp[0].conj();
    for k in 1..n {
        b[k] = chirp[k].conj();
        b[m - k] = chirp[k].conj();
    }
    radix2(&mut a, false);
    radix2(&mut b, false);
    for k in 0..m {
        a[k] *= b[k];
    }
    radix2(&mut a, true);
    for k in 0..n {
        data[k] = a[k] / m as f64 * chirp[k];
        if inverse {
            data[k] /= n as f64;
        }
    }
    data
}

/// Dense row-major array with a shape.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NdArray {
    pub shape: Vec<usize>,
    pub values: Vec<f64>,
}

impl NdArray {
    pub fn new(shape: Vec<usize>, fill: f64) -> Self {
        let count = shape.iter().product();
        Self { shape, values: vec![fill; count] }
    }
    pub fn ndim(&self) -> usize {
        self.shape.len()
    }
    pub fn size(&self) -> usize {
        self.values.len()
    }
    pub fn at2(&self, i: usize, j: usize) -> f64 {
        self.values[i * self.shape[1] + j]
    }
    pub fn at3(&self, i: usize, j: usize, k: usize) -> f64 {
        self.values[(i * self.shape[1] + j) * self.shape[2] + k]
    }
}

fn flatten(value: &Json, depth: usize, shape: &[usize], out: &mut Vec<f64>, name: &str) -> Result<(), String> {
    if depth < shape.len() {
        if !value.is_array() || value.size() != shape[depth] {
            return Err(format!("{name} must be a rectangular numeric array"));
        }
        for item in value.items()? {
            flatten(item, depth + 1, shape, out, name)?;
        }
        return Ok(());
    }
    match value {
        Json::Number(n) => out.push(*n),
        Json::Bool(b) => out.push(if *b { 1.0 } else { 0.0 }),
        Json::Array(_) => return Err(format!("{name} must be a rectangular numeric array")),
        _ => return Err(format!("{name} must contain only numbers")),
    }
    Ok(())
}

/// A rectangular numeric JSON array (or a scalar, shape []).
pub fn to_array(value: &Json, name: &str) -> Result<NdArray, String> {
    let mut result = NdArray::default();
    let mut current = value;
    while let Json::Array(items) = current {
        result.shape.push(items.len());
        if items.is_empty() {
            break;
        }
        current = &items[0];
    }
    flatten(value, 0, &result.shape.clone(), &mut result.values, name)?;
    let count: usize = result.shape.iter().product();
    if count != result.values.len() {
        return Err(format!("{name} must be a rectangular numeric array"));
    }
    Ok(result)
}

pub fn array_json(array: &NdArray) -> Json {
    fn part(a: &NdArray, axis: usize, offset: usize) -> Json {
        if axis == a.shape.len() {
            return Json::Number(a.values[offset]);
        }
        let stride: usize = a.shape[axis + 1..].iter().product();
        Json::Array((0..a.shape[axis]).map(|i| part(a, axis + 1, offset + i * stride)).collect())
    }
    if array.shape.is_empty() {
        return array.values.first().map_or(Json::Null, |&v| Json::Number(v));
    }
    part(array, 0, 0)
}

pub fn vec3_json(v: Vec3) -> Json {
    Json::Array(v.iter().map(|&x| Json::Number(x)).collect())
}

pub fn mat3_json(m: &Mat3) -> Json {
    Json::Array(m.iter().map(|r| vec3_json(*r)).collect())
}

pub fn finite_array(array: NdArray, name: &str, ndim: i32) -> Result<NdArray, String> {
    let finite = array.size() > 0 && array.values.iter().all(|v| v.is_finite());
    if !finite || (ndim >= 0 && array.ndim() != ndim as usize) {
        return Err(format!("{name} must be a nonempty finite array of the required dimension"));
    }
    Ok(array)
}

pub fn positive(value: f64, name: &str) -> Result<f64, String> {
    if !value.is_finite() || value <= 0.0 {
        return Err(format!("{name} must be finite and positive"));
    }
    Ok(value)
}

pub fn integer(value: f64, name: &str, minimum: i64) -> Result<i64, String> {
    if !value.is_finite() || value != value.floor() || value < minimum as f64 || value.abs() > 9e15 {
        return Err(format!("{name} must be an integer >= {minimum}"));
    }
    Ok(value as i64)
}

/// A 3x3 cell of independent, well-conditioned rows.
pub fn cell_matrix(cell: NdArray) -> Result<Mat3, String> {
    let m = finite_array(cell, "cell", 2)?;
    let bad = || "cell must contain three independent, well-conditioned row vectors".to_string();
    if m.shape[0] != 3 || m.shape[1] != 3 {
        return Err(bad());
    }
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = m.at2(i, j);
        }
    }
    let matrix = nalgebra::Matrix3::from_fn(|i, j| result[i][j]);
    let singular = matrix.singular_values();
    let (largest, smallest) = (singular.max(), singular.min());
    if matrix.determinant().abs() < 1e-12 || smallest <= 0.0 || largest / smallest > 1e10 {
        return Err(bad());
    }
    Ok(result)
}

pub fn points(array: NdArray, name: &str) -> Result<Vec<Vec3>, String> {
    let checked = finite_array(array, name, 2)?;
    if checked.shape[1] != 3 {
        return Err(format!("{name} must have shape (atoms, 3)"));
    }
    Ok((0..checked.shape[0]).map(|i| [checked.at2(i, 0), checked.at2(i, 1), checked.at2(i, 2)]).collect())
}

pub fn parse_pbc(value: &Json) -> Result<Pbc, String> {
    let bad = || "pbc must contain three booleans".to_string();
    let items = value.items().map_err(|_| bad())?;
    if items.len() != 3 {
        return Err(bad());
    }
    let mut result = [false; 3];
    for k in 0..3 {
        result[k] = items[k].boolean().map_err(|_| bad())?;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fft_inverts() {
        for n in [1usize, 2, 7, 8, 12] {
            let data: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, (i * i) as f64 * 0.5)).collect();
            let back = fft(fft(data.clone(), false), true);
            for (a, b) in data.iter().zip(&back) {
                assert!((a - b).norm() < 1e-9);
            }
        }
    }

    #[test]
    fn minimum_image_of_skewed_cell() {
        let cell = [[3.0, 0.0, 0.0], [2.9, 0.5, 0.0], [0.0, 0.0, 4.0]];
        let image = MinimumImage::new(&cell, [true, true, true]).unwrap();
        let d = image.apply([2.9, 0.5, 0.0]);
        assert!(norm(d) < 0.6, "{d:?}");
    }
}
