//! Molecular point groups (the C++ chemistry/PointGroup): Schoenflies symbols
//! within a distance tolerance, symmetrisation (averaging every atom over the
//! group so the symmetry is exact) and the standard orientation (principal
//! axis along z).

use af_core::crmath;
use af_core::eigen::self_adjoint_eigen3;
use af_core::elements::atomic_mass;

use crate::molecule::{fmax, fmin};

pub type Mat3x3 = [[f64; 3]; 3];
type V = [f64; 3];

const PI: f64 = 3.14159265358979323846;

#[derive(Clone, Debug, PartialEq)]
pub struct PointGroupResult {
    /// e.g. C2v, D6h, Td, Cinfv.
    pub name: String,
    /// Number of operations (0 for linear groups).
    pub order: i32,
    pub rotational_symmetry_number: i32,
    /// About the centre of mass.
    pub operations: Vec<Mat3x3>,
    pub centre: V,
    pub principal_axis: V,
    /// A C2 or mirror-plane direction perpendicular to the principal axis.
    pub secondary_axis: V,
    /// Largest distance of an atom image from its partner (Angstrom).
    pub max_deviation: f64,
}

impl Default for PointGroupResult {
    fn default() -> Self {
        Self {
            name: String::new(),
            order: 1,
            rotational_symmetry_number: 1,
            operations: Vec::new(),
            centre: [0.0; 3],
            principal_axis: [0.0, 0.0, 1.0],
            secondary_axis: [1.0, 0.0, 0.0],
            max_deviation: 0.0,
        }
    }
}

fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: V, b: V) -> V {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: V, s: f64) -> V {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V, b: V) -> V {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: V) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: V) -> V {
    let n = norm(a);
    if n > 1e-12 {
        scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 1.0]
    }
}
fn transform(m: &Mat3x3, v: V) -> V {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
fn multiply(a: &Mat3x3, b: &Mat3x3) -> Mat3x3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
fn transposed(a: &Mat3x3) -> Mat3x3 {
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}
const IDENTITY: Mat3x3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
fn rotation(axis: V, angle: f64) -> Mat3x3 {
    let u = unit(axis);
    let (c, s) = (angle.cos(), angle.sin());
    let t = 1.0 - c;
    [
        [t * u[0] * u[0] + c, t * u[0] * u[1] - s * u[2], t * u[0] * u[2] + s * u[1]],
        [t * u[0] * u[1] + s * u[2], t * u[1] * u[1] + c, t * u[1] * u[2] - s * u[0]],
        [t * u[0] * u[2] - s * u[1], t * u[1] * u[2] + s * u[0], t * u[2] * u[2] + c],
    ]
}
fn reflection(normal: V) -> Mat3x3 {
    let n = unit(normal);
    let mut m = IDENTITY;
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] -= 2.0 * n[i] * n[j];
        }
    }
    m
}
fn distance_of(a: &Mat3x3, b: &Mat3x3) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            s += (a[i][j] - b[i][j]) * (a[i][j] - b[i][j]);
        }
    }
    s.sqrt()
}
fn determinant(m: &Mat3x3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

struct Frame<'a> {
    z: &'a [i32],
    /// Relative to the centre of mass.
    x: Vec<V>,
    tolerance: f64,
}

/// The largest displacement between an atom's image and its nearest partner
/// of the same element, or a negative value when some atom has none.
fn mismatch(f: &Frame, op: &Mat3x3) -> f64 {
    let mut worst = 0.0;
    for i in 0..f.x.len() {
        let image = transform(op, f.x[i]);
        let mut best = 1e300;
        for j in 0..f.x.len() {
            if f.z[j] != f.z[i] {
                continue;
            }
            best = fmin(best, norm(sub(image, f.x[j])));
        }
        if best > f.tolerance {
            return -1.0;
        }
        worst = fmax(worst, best);
    }
    worst
}

fn holds(f: &Frame, op: &Mat3x3) -> bool {
    mismatch(f, op) >= 0.0
}

fn add_direction(list: &mut Vec<V>, v: V) {
    if norm(v) < 1e-6 {
        return;
    }
    let u = unit(v);
    if list.iter().any(|&w| (dot(u, w).abs() - 1.0).abs() < 1e-4) {
        return;
    }
    list.push(u);
}

#[derive(Clone, Copy)]
struct Axis {
    direction: V,
    order: i32,
    error: f64,
}

struct Classified {
    name: String,
    principal: V,
    secondary: V,
}

/// Classification (Schoenflies) from rotation axes, mirror planes, the
/// inversion and the presence of improper rotations.
fn classify(principal: &[V; 3], ax: &[Axis], mir: &[V], inv: bool, improper: bool) -> Classified {
    let mut c = Classified { name: "C1".into(), principal: principal[2], secondary: principal[0] };
    let mut highest = 1;
    for a in ax {
        highest = highest.max(a.order);
    }
    // Cubic and icosahedral groups: two C3 axes at the tetrahedral angle
    // (|cos| = 1/3) or, for I and Ih, at 41.8 degrees (cos 0.745).
    let mut cubic = false;
    'outer: for x in 0..ax.len() {
        for y in x + 1..ax.len() {
            if ax[x].order % 3 != 0 || ax[y].order % 3 != 0 {
                continue;
            }
            let cs = dot(ax[x].direction, ax[y].direction).abs();
            if (cs - 1.0 / 3.0).abs() < 0.06 || (cs - 0.7454).abs() < 0.06 {
                cubic = true;
                break 'outer;
            }
        }
    }
    if cubic {
        let five = ax.iter().any(|a| a.order == 5);
        let four = ax.iter().any(|a| a.order == 4);
        c.name = if five {
            if inv { "Ih" } else { "I" }
        } else if four {
            if inv { "Oh" } else { "O" }
        } else if inv {
            "Th"
        } else if mir.is_empty() {
            "T"
        } else {
            "Td"
        }
        .into();
        let wanted = if five { 5 } else if four { 4 } else { 2 };
        if let Some(a) = ax.iter().find(|a| a.order == wanted || (wanted == 2 && a.order % 2 == 0)) {
            c.principal = a.direction;
        }
        if let Some(a) = ax.iter().find(|a| a.order % 2 == 0 && dot(a.direction, c.principal).abs() < 0.05) {
            c.secondary = a.direction;
        }
        return c;
    }
    if highest == 1 {
        if !mir.is_empty() {
            c.name = "Cs".into();
            c.principal = mir[0];
        } else if inv {
            c.name = "Ci".into();
        }
        return c;
    }
    if let Some(a) = ax.iter().find(|a| a.order == highest) {
        c.principal = a.direction;
    }
    let mut perpendicular_c2 = 0;
    for a in ax {
        if a.order % 2 == 0 && dot(a.direction, c.principal).abs() < 0.05 {
            if perpendicular_c2 == 0 {
                c.secondary = a.direction;
            }
            perpendicular_c2 += 1;
        }
    }
    let mut horizontal = false;
    let mut vertical = 0;
    for &m in mir {
        if (dot(m, c.principal).abs() - 1.0).abs() < 0.05 {
            horizontal = true;
        } else if dot(m, c.principal).abs() < 0.05 {
            if perpendicular_c2 == 0 && vertical == 0 {
                c.secondary = unit(cross(c.principal, m));
            }
            vertical += 1;
        }
    }
    c.name = if perpendicular_c2 >= highest {
        if horizontal {
            format!("D{highest}h")
        } else if vertical >= highest {
            format!("D{highest}d")
        } else {
            format!("D{highest}")
        }
    } else if horizontal {
        format!("C{highest}h")
    } else if vertical >= highest {
        format!("C{highest}v")
    } else if improper {
        format!("S{}", 2 * highest)
    } else {
        format!("C{highest}")
    };
    c
}

/// Point group of atoms (atomic numbers z, positions) within `tolerance` Angstrom.
pub fn point_group(z: &[i32], positions: &[V], tolerance: f64) -> PointGroupResult {
    let mut result = PointGroupResult::default();
    let n = positions.len();
    if n == 0 {
        result.name = "C1".into();
        return result;
    }
    let mut mass = 0.0;
    let mut com = [0.0; 3];
    for i in 0..n {
        let w = atomic_mass(z[i]);
        mass += w;
        com = add(com, scale(positions[i], w));
    }
    com = scale(com, 1.0 / mass);
    result.centre = com;
    let f = Frame { z, x: positions.iter().map(|&p| sub(p, com)).collect(), tolerance };
    if n == 1 {
        result.name = "Kh".into();
        result.order = 0;
        return result;
    }
    // Principal axes of inertia.
    let mut inertia = [[0.0; 3]; 3];
    for i in 0..n {
        let w = atomic_mass(z[i]);
        let r = f.x[i];
        let r2 = dot(r, r);
        for a in 0..3 {
            for b in 0..3 {
                inertia[a][b] += w * ((if a == b { r2 } else { 0.0 }) - r[a] * r[b]);
            }
        }
    }
    let vectors = self_adjoint_eigen3(&inertia).map(|(_, v)| v).unwrap_or(IDENTITY);
    let principal: [V; 3] = std::array::from_fn(|k| [vectors[0][k], vectors[1][k], vectors[2][k]]);
    let inversion: Mat3x3 = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
    let has_inversion = holds(&f, &inversion);
    // Linear molecules: all atoms on a line through the centre.
    {
        let line = principal[0];
        let linear = f.x.iter().all(|&r| !(norm(cross(r, line)) > tolerance));
        if linear {
            result.name = if has_inversion { "Dinfh" } else { "Cinfv" }.into();
            result.order = 0;
            result.rotational_symmetry_number = if has_inversion { 2 } else { 1 };
            result.principal_axis = line;
            result.operations = vec![IDENTITY];
            if has_inversion {
                result.operations.push(inversion);
            }
            return result;
        }
    }
    // Candidate axes and plane normals.
    let mut candidates: Vec<V> = Vec::new();
    for &p in &principal {
        add_direction(&mut candidates, p);
    }
    let mut i = 0;
    while i < n && candidates.len() < 400 {
        add_direction(&mut candidates, f.x[i]);
        i += 1;
    }
    let mut i = 0;
    while i < n && candidates.len() < 1200 {
        let mut j = i + 1;
        while j < n && candidates.len() < 1200 {
            if z[i] == z[j] && !((norm(f.x[i]) - norm(f.x[j])).abs() > tolerance) {
                add_direction(&mut candidates, add(f.x[i], f.x[j])); // midpoints: C2 axes
                add_direction(&mut candidates, sub(f.x[i], f.x[j])); // mirror normals
                add_direction(&mut candidates, cross(f.x[i], f.x[j])); // normals of planes through two atoms
            }
            j += 1;
        }
        i += 1;
    }
    // Axes and mirror normals, merging directions within a few degrees (the
    // candidates of a slightly distorted molecule scatter about each element).
    let same = crmath::cos(4.0 * PI / 180.0);
    let mut axes: Vec<Axis> = Vec::new();
    for &c in &candidates {
        for order in (2..=8).rev() {
            let e = mismatch(&f, &rotation(c, 2.0 * PI / f64::from(order)));
            if e < 0.0 {
                continue;
            }
            let mut merged = false;
            for a in axes.iter_mut() {
                if dot(a.direction, c).abs() > same {
                    merged = true;
                    if order > a.order || (order == a.order && e < a.error) {
                        *a = Axis { direction: c, order, error: e };
                    }
                }
            }
            if !merged {
                axes.push(Axis { direction: c, order, error: e });
            }
            break;
        }
    }
    let mut mirrors: Vec<V> = Vec::new();
    let mut mirror_error: Vec<f64> = Vec::new();
    for &c in &candidates {
        let e = mismatch(&f, &reflection(c));
        if e < 0.0 {
            continue;
        }
        let mut merged = false;
        for k in 0..mirrors.len() {
            if dot(mirrors[k], c).abs() > same {
                merged = true;
                if e < mirror_error[k] {
                    mirrors[k] = c;
                    mirror_error[k] = e;
                }
            }
        }
        if !merged {
            mirrors.push(c);
            mirror_error.push(e);
        }
    }
    let mut generators: Vec<Mat3x3> = Vec::new();
    for a in &axes {
        generators.push(rotation(a.direction, 2.0 * PI / f64::from(a.order)));
    }
    for &m in &mirrors {
        generators.push(reflection(m));
    }
    if has_inversion {
        generators.push(inversion);
    }
    {
        // S2n along the highest axis: rotation by pi/n with the horizontal reflection.
        let mut highest = 1;
        let mut along = principal[2];
        for a in &axes {
            if a.order > highest {
                highest = a.order;
                along = a.direction;
            }
        }
        if highest > 1 {
            let s_op = multiply(&reflection(along), &rotation(along, PI / f64::from(highest)));
            if holds(&f, &s_op) {
                generators.push(s_op);
            }
        }
    }
    // The full group from the generators (closure), then the operations that hold.
    let mut group: Vec<Mat3x3> = vec![IDENTITY];
    for g in &generators {
        if !group.iter().any(|m| distance_of(m, g) < 0.15) {
            group.push(*g);
        }
    }
    let mut a = 0;
    while a < group.len() && group.len() < 240 {
        let mut b = 0;
        while b < group.len() && group.len() < 240 {
            let pr = multiply(&group[a], &group[b]);
            if !group.iter().any(|m| distance_of(m, &pr) < 0.15) {
                group.push(pr);
            }
            b += 1;
        }
        a += 1;
    }
    let mut worst = 0.0;
    let mut valid: Vec<Mat3x3> = Vec::new();
    for m in &group {
        let d = mismatch(&f, m);
        if d >= 0.0 {
            valid.push(*m);
            worst = fmax(worst, d);
        }
    }
    // The elements of the closed group (axes such as the C3 of an octahedron
    // need not be among the candidate directions).
    let mut group_axes: Vec<Axis> = Vec::new();
    let mut group_mirrors: Vec<V> = Vec::new();
    let (mut group_inversion, mut group_improper) = (false, false);
    let largest_column = |m: &Mat3x3| {
        let mut best = [0.0, 0.0, 1.0];
        let mut best_norm = -1.0;
        for col in 0..3 {
            let v = [m[0][col], m[1][col], m[2][col]];
            if norm(v) > best_norm {
                best_norm = norm(v);
                best = v;
            }
        }
        unit(best)
    };
    for m in &valid {
        let det = determinant(m);
        let trace = m[0][0] + m[1][1] + m[2][2];
        if det > 0.0 {
            let angle = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0).acos();
            if angle < 1e-2 {
                continue;
            }
            let order = (2.0 * PI / angle).round() as i32;
            if !(2..=8).contains(&order) {
                continue;
            }
            let mut axis = [m[2][1] - m[1][2], m[0][2] - m[2][0], m[1][0] - m[0][1]];
            if norm(axis) < 1e-3 {
                let mut plus = *m;
                for d in 0..3 {
                    plus[d][d] += 1.0;
                }
                axis = largest_column(&plus);
            }
            axis = unit(axis);
            let mut merged = false;
            for a in group_axes.iter_mut() {
                if dot(a.direction, axis).abs() > same {
                    merged = true;
                    a.order = a.order.max(order);
                }
            }
            if !merged {
                group_axes.push(Axis { direction: axis, order, error: 0.0 });
            }
        } else if (trace - 1.0).abs() < 0.1 {
            let mut minus = IDENTITY;
            for r in 0..3 {
                for c2 in 0..3 {
                    minus[r][c2] -= m[r][c2];
                }
            }
            let normal = largest_column(&minus);
            if !group_mirrors.iter().any(|&v| dot(v, normal).abs() > same) {
                group_mirrors.push(normal);
            }
        } else if (trace + 3.0).abs() < 0.1 {
            group_inversion = true;
        } else {
            group_improper = true;
        }
    }
    let c = classify(&principal, &group_axes, &group_mirrors, group_inversion, group_improper);
    result.name = c.name;
    result.order = valid.len() as i32;
    result.rotational_symmetry_number = valid.iter().filter(|m| determinant(m) > 0.0).count() as i32;
    result.operations = valid;
    result.principal_axis = unit(c.principal);
    let p = result.principal_axis;
    let mut sec = sub(c.secondary, scale(p, dot(c.secondary, p)));
    if norm(sec) < 1e-6 {
        sec = if p[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    }
    sec = sub(sec, scale(p, dot(sec, p)));
    result.secondary_axis = unit(sec);
    result.max_deviation = worst;
    result
}

/// Positions averaged over the group so the symmetry holds exactly (repeated
/// detection with a shrinking tolerance). Returns the group found.
pub fn symmetrize(z: &[i32], positions: &mut [V], tolerance: f64) -> PointGroupResult {
    let mut group = point_group(z, positions, tolerance);
    for _ in 0..4 {
        if group.operations.len() <= 1 {
            break;
        }
        let n = positions.len();
        let x: Vec<V> = positions.iter().map(|&p| sub(p, group.centre)).collect();
        let mut sum = vec![[0.0; 3]; n];
        for op in &group.operations {
            let inverse = transposed(op);
            for i in 0..n {
                // Partner j of i under op: op x_i ~ x_j; then op^-1 x_j estimates x_i.
                let image = transform(op, x[i]);
                let mut best = i;
                let mut best_distance = 1e300;
                for j in 0..n {
                    if z[j] != z[i] {
                        continue;
                    }
                    let d = norm(sub(image, x[j]));
                    if d < best_distance {
                        best_distance = d;
                        best = j;
                    }
                }
                sum[i] = add(sum[i], transform(&inverse, x[best]));
            }
        }
        let count = group.operations.len() as f64;
        for i in 0..n {
            positions[i] = add(scale(sum[i], 1.0 / count), group.centre);
        }
        let name = group.name.clone();
        group = point_group(z, positions, fmax(1e-3, tolerance * 0.1));
        if group.name != name {
            group = point_group(z, positions, tolerance);
        }
        if group.max_deviation < 1e-6 {
            break;
        }
    }
    group
}

/// Rotates (and translates) positions so the centre of mass is at the origin,
/// the principal symmetry axis is along z and a secondary element along x;
/// without symmetry axes, the principal axes of inertia.
pub fn standard_orientation(positions: &mut [V], group: &PointGroupResult) {
    let e3 = group.principal_axis;
    let e1 = unit(sub(group.secondary_axis, scale(e3, dot(group.secondary_axis, e3))));
    let e2 = cross(e3, e1);
    for p in positions.iter_mut() {
        let r = sub(*p, group.centre);
        *p = [dot(r, e1), dot(r, e2), dot(r, e3)];
    }
}
