//! A sampled scalar field on a (possibly periodic) grid (the C++
//! electronic/Volume Grid), with glm's double-precision arithmetic.

use af_core::glm::{determinant_d, inverse_d, mul_d};
use glam::{DMat3, DVec3};

pub const BOHR: f64 = 0.529177210903;
pub const HARTREE: f64 = 27.211386245988;
pub const PI: f64 = 3.14159265358979323846;

#[derive(Clone, Debug, PartialEq)]
pub struct Grid {
    pub shape: [i32; 3],
    /// Columns are the cell vectors (glm).
    pub cell: DMat3,
    pub origin: DVec3,
    pub periodic: bool,
    pub name: String,
    pub unit: String,
    /// x fastest: index = (z * ny + y) * nx + x.
    pub values: Vec<f64>,
}

impl Default for Grid {
    fn default() -> Self {
        Self { shape: [2, 2, 2], cell: DMat3::IDENTITY, origin: DVec3::ZERO, periodic: false, name: "field".into(), unit: "raw".into(), values: Vec::new() }
    }
}

/// Triangles (three vertices each) with a scalar per vertex.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    pub vertices: Vec<DVec3>,
    pub colors: Vec<f64>,
}

impl Grid {
    pub fn validate(&self) -> Result<(), String> {
        let mut count: usize = 1;
        for &n in &self.shape {
            if !(2..=100000).contains(&n) || count > 100_000_000 / n as usize {
                return Err("Grid dimensions must be >= 2 and total samples <= 100 million".into());
            }
            count *= n as usize;
        }
        if self.values.len() != count {
            return Err("Grid value count does not match dimensions".into());
        }
        for i in 0..3 {
            if !self.origin[i].is_finite() {
                return Err("Non-finite grid origin".into());
            }
            for j in 0..3 {
                if !self.cell.col(i)[j].is_finite() {
                    return Err("Non-finite lattice".into());
                }
            }
        }
        let determinant = determinant_d(&self.cell);
        if !determinant.is_finite() || determinant.abs() < 1e-12 {
            return Err("Singular grid lattice".into());
        }
        if self.values.iter().any(|v| !v.is_finite()) {
            return Err("Non-finite grid value".into());
        }
        Ok(())
    }

    pub fn index(&self, x: i32, y: i32, z: i32) -> usize {
        (z as usize * self.shape[1] as usize + y as usize) * self.shape[0] as usize + x as usize
    }

    /// The cell vectors divided into grid steps.
    pub fn steps(&self) -> DMat3 {
        let d = |i: usize| (self.shape[i] - if self.periodic { 0 } else { 1 }) as f64;
        DMat3::from_cols(self.cell.x_axis / d(0), self.cell.y_axis / d(1), self.cell.z_axis / d(2))
    }

    pub fn position(&self, x: i32, y: i32, z: i32) -> DVec3 {
        self.origin + mul_d(&self.steps(), DVec3::new(x as f64, y as f64, z as f64))
    }

    /// Integration weight of a sample (trapezoid ends on open grids).
    pub fn weight(&self, x: i32, y: i32, z: i32) -> f64 {
        let mut w = determinant_d(&self.steps()).abs();
        if !self.periodic {
            if x == 0 || x == self.shape[0] - 1 {
                w *= 0.5;
            }
            if y == 0 || y == self.shape[1] - 1 {
                w *= 0.5;
            }
            if z == 0 || z == self.shape[2] - 1 {
                w *= 0.5;
            }
        }
        w
    }

    /// Trilinear interpolation at a point.
    pub fn sample(&self, point: DVec3) -> Result<f64, String> {
        let mut q = mul_d(&inverse_d(&self.steps()), point - self.origin);
        let (mut low, mut high, mut t) = ([0i32; 3], [0i32; 3], [0.0f64; 3]);
        for i in 0..3 {
            if !q[i].is_finite() {
                return Err("Non-finite sampling point".into());
            }
            let n = self.shape[i];
            if self.periodic {
                q[i] -= (q[i] / n as f64).floor() * n as f64;
            } else {
                if q[i] < -1e-8 || q[i] > (n - 1) as f64 + 1e-8 {
                    return Err("Sampling point lies outside finite grid".into());
                }
                q[i] = q[i].clamp(0.0, (n - 1) as f64);
            }
            low[i] = (q[i].floor() as i32).min(n - 1);
            high[i] = if self.periodic { (low[i] + 1) % n } else { (low[i] + 1).min(n - 1) };
            t[i] = q[i] - low[i] as f64;
        }
        let mut value = 0.0;
        for z in 0..2 {
            for y in 0..2 {
                for x in 0..2 {
                    let v = self.values[self.index(if x == 1 { high[0] } else { low[0] }, if y == 1 { high[1] } else { low[1] }, if z == 1 { high[2] } else { low[2] })];
                    value += v * (if x == 1 { t[0] } else { 1.0 - t[0] }) * (if y == 1 { t[1] } else { 1.0 - t[1] }) * (if z == 1 { t[2] } else { 1.0 - t[2] });
                }
            }
        }
        Ok(value)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Site {
    pub number: i32,
    pub position: DVec3,
}

/// Atoms and fields of a volumetric file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Volume {
    pub sites: Vec<Site>,
    pub fields: Vec<Grid>,
}

pub fn require_aligned(a: &Grid, b: &Grid) -> Result<(), String> {
    a.validate()?;
    b.validate()?;
    if a.shape != b.shape || a.periodic != b.periodic || (a.origin - b.origin).length() > 1e-8 {
        return Err("Grids are not aligned; resample explicitly first".into());
    }
    for i in 0..3 {
        if (a.cell.col(i) - b.cell.col(i)).length() > 1e-8 {
            return Err("Grid lattices differ; resample explicitly first".into());
        }
    }
    Ok(())
}

pub fn arithmetic(a: &Grid, b: &Grid, operation: &str) -> Result<Grid, String> {
    require_aligned(a, b)?;
    if (operation == "add" || operation == "subtract") && a.unit != b.unit {
        return Err("Addition and subtraction require matching units".into());
    }
    let mut out = a.clone();
    out.name = operation.to_string();
    match operation {
        "multiply" => out.unit = format!("({})*({})", a.unit, b.unit),
        "divide" => out.unit = if a.unit == b.unit { "1".into() } else { format!("({})/({})", a.unit, b.unit) },
        "add" | "subtract" => {}
        _ => return Err("Unknown arithmetic operation".into()),
    }
    for (v, w) in out.values.iter_mut().zip(&b.values) {
        match operation {
            "add" => *v += w,
            "subtract" => *v -= w,
            "multiply" => *v *= w,
            _ => {
                if *w == 0.0 {
                    return Err("Division by zero in grid arithmetic".into());
                }
                *v /= w;
            }
        }
    }
    out.validate()?;
    Ok(out)
}

pub fn scale(a: &Grid, factor: f64) -> Result<Grid, String> {
    a.validate()?;
    if !factor.is_finite() {
        return Err("Scale must be finite".into());
    }
    let mut out = a.clone();
    for v in out.values.iter_mut() {
        *v *= factor;
    }
    out.validate()?;
    Ok(out)
}

/// `source` sampled at the points of `target` (trilinear).
pub fn resample(source: &Grid, target: &Grid) -> Result<Grid, String> {
    source.validate()?;
    target.validate()?;
    let mut out = target.clone();
    out.name = source.name.clone();
    out.unit = source.unit.clone();
    for z in 0..out.shape[2] {
        for y in 0..out.shape[1] {
            for x in 0..out.shape[0] {
                let index = out.index(x, y, z);
                out.values[index] = source.sample(out.position(x, y, z))?;
            }
        }
    }
    Ok(out)
}

/// The shortest periodic image of `delta` (cell columns), for any lattice.
pub fn minimum_image(delta: DVec3, cell: &DMat3) -> Result<DVec3, String> {
    // Enumerate all translations that can improve the initial candidate; unlike
    // component rounding alone, this also works for non-reduced triclinic cells.
    let inverse = inverse_d(cell);
    let f = mul_d(&inverse, delta);
    let rounded = DVec3::new(f.x.round(), f.y.round(), f.z.round());
    let mut best = mul_d(cell, f - rounded);
    let radius = best.length();
    let (mut low, mut high) = ([0i32; 3], [0i32; 3]);
    let mut count: usize = 1;
    for i in 0..3 {
        let row = DVec3::new(inverse.col(0)[i], inverse.col(1)[i], inverse.col(2)[i]);
        let bound = row.length() * radius + 1e-10;
        if !f[i].is_finite() || bound > 10000.0 || f[i].abs() > 1e8 {
            return Err("Ill-conditioned lattice or invalid coordinates".into());
        }
        low[i] = (f[i] - bound).ceil() as i32;
        high[i] = (f[i] + bound).floor() as i32;
        count *= 1.max(high[i] - low[i] + 1) as usize;
        if count > 1_000_000 {
            return Err("Lattice requires excessive periodic images; reduce the cell".into());
        }
    }
    for z in low[2]..=high[2] {
        for y in low[1]..=high[1] {
            for x in low[0]..=high[0] {
                let candidate = delta - mul_d(cell, DVec3::new(x as f64, y as f64, z as f64));
                if candidate.dot(candidate) < best.dot(best) {
                    best = candidate;
                }
            }
        }
    }
    Ok(best)
}

/// A finite grid whose opposite end planes repeat, as a periodic grid without them.
pub fn periodic_endpoints(source: &Grid, tolerance: f64) -> Result<Grid, String> {
    source.validate()?;
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err("Endpoint tolerance must be nonnegative and finite".into());
    }
    if source.periodic {
        return Ok(source.clone());
    }
    let mut out = source.clone();
    for n in out.shape.iter_mut() {
        if *n < 3 {
            return Err("Need >=3 samples to remove periodic endpoint planes".into());
        }
        *n -= 1;
    }
    for z in 0..source.shape[2] {
        for y in 0..source.shape[1] {
            for x in 0..source.shape[0] {
                let a = source.values[source.index(x, y, z)];
                let b = source.values[source.index(x % out.shape[0], y % out.shape[1], z % out.shape[2])];
                let mut largest = 1.0f64;
                for c in [a.abs(), b.abs()] {
                    if largest < c {
                        largest = c;
                    }
                }
                if (a - b).abs() > tolerance * largest {
                    return Err("Opposite endpoint planes do not match".into());
                }
            }
        }
    }
    out.periodic = true;
    out.values = vec![0.0; out.shape[0] as usize * out.shape[1] as usize * out.shape[2] as usize];
    for z in 0..out.shape[2] {
        for y in 0..out.shape[1] {
            for x in 0..out.shape[0] {
                let index = out.index(x, y, z);
                out.values[index] = source.values[source.index(x, y, z)];
            }
        }
    }
    Ok(out)
}
