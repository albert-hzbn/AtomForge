//! Shapes for carving crystals, modelled as in CAD tools (the C++
//! algorithms/ShapeModel): primitives, booleans (also smooth), transforms,
//! modifiers, patterns and element regions, written as a small OpenSCAD-like
//! script and evaluated as a signed distance function (negative inside, A).

use af_core::glm::{inverse_d, mul_d, normalize_d};
use glam::{DMat3, DVec2, DVec3};

const PI: f64 = 3.14159265358979323846;
const DEGREE: f64 = PI / 180.0;
const HUGE: f64 = 1e30;

/// A script error with its 1-based line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShapeError {
    pub line: i32,
    pub message: String,
}

impl std::fmt::Display for ShapeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for ShapeError {}

fn err<T>(line: i32, message: impl Into<String>) -> Result<T, ShapeError> {
    Err(ShapeError { line, message: message.into() })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Bounds {
    pub low: DVec3,
    pub high: DVec3,
    pub finite: bool,
    /// No volume (regions).
    pub empty: bool,
}

// std::min / std::max / std::clamp on doubles.
fn min_d(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}
fn max_d(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}
fn clamp_d(v: f64, lo: f64, hi: f64) -> f64 {
    if v < lo {
        lo
    } else if hi < v {
        hi
    } else {
        v
    }
}
fn vmin(a: DVec3, b: DVec3) -> DVec3 {
    DVec3::new(min_d(a.x, b.x), min_d(a.y, b.y), min_d(a.z, b.z))
}
fn vmax(a: DVec3, b: DVec3) -> DVec3 {
    DVec3::new(max_d(a.x, b.x), max_d(a.y, b.y), max_d(a.z, b.z))
}
fn len3(v: DVec3) -> f64 {
    v.dot(v).sqrt()
}
fn len2(v: DVec2) -> f64 {
    v.dot(v).sqrt()
}

// ---------------------------------------------------------------- values
#[derive(Clone, Debug, PartialEq)]
enum Value {
    Number(f64),
    Vector(Vec<Value>),
    Text(String),
    Boolean(bool),
}

fn describe(v: &Value) -> &'static str {
    match v {
        Value::Number(_) => "a number",
        Value::Vector(_) => "a vector",
        Value::Text(_) => "text",
        Value::Boolean(_) => "true/false",
    }
}

// ---------------------------------------------------------------- lexer
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TokenKind {
    Identifier,
    Number,
    Text,
    Symbol,
    End,
}

#[derive(Clone, Debug)]
struct Token {
    kind: TokenKind,
    text: String,
    number: f64,
    line: i32,
}

fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn tokenize(source: &str) -> Result<Vec<Token>, ShapeError> {
    let s = source.as_bytes();
    let mut out = Vec::new();
    let mut line = 1;
    let mut i = 0;
    while i < s.len() {
        let c = s[i];
        if c == b'\n' {
            line += 1;
            i += 1;
            continue;
        }
        if is_space(c) {
            i += 1;
            continue;
        }
        if c == b'#' || (c == b'/' && s.get(i + 1) == Some(&b'/')) {
            while i < s.len() && s[i] != b'\n' {
                i += 1;
            }
            continue;
        }
        if c == b'/' && s.get(i + 1) == Some(&b'*') {
            i += 2;
            while i + 1 < s.len() && !(s[i] == b'*' && s[i + 1] == b'/') {
                if s[i] == b'\n' {
                    line += 1;
                }
                i += 1;
            }
            i += 2;
            continue;
        }
        let mut t = Token { kind: TokenKind::End, text: String::new(), number: 0.0, line };
        if c.is_ascii_alphabetic() || c == b'_' {
            let mut j = i;
            while j < s.len() && (s[j].is_ascii_alphanumeric() || s[j] == b'_') {
                j += 1;
            }
            t.kind = TokenKind::Identifier;
            t.text = String::from_utf8_lossy(&s[i..j]).into_owned();
            i = j;
        } else if c.is_ascii_digit() || (c == b'.' && s.get(i + 1).is_some_and(|d| d.is_ascii_digit())) {
            let rest = &source[i..];
            let (value, used) = af_core::format::strtod(rest);
            t.kind = TokenKind::Number;
            t.number = value;
            t.text = rest[..used].to_string();
            i += used;
        } else if c == b'"' || c == b'\'' {
            let mut j = i + 1;
            while j < s.len() && s[j] != c && s[j] != b'\n' {
                j += 1;
            }
            if j >= s.len() || s[j] != c {
                return err(line, "unterminated text");
            }
            t.kind = TokenKind::Text;
            t.text = String::from_utf8_lossy(&s[i + 1..j]).into_owned();
            i = j + 1;
        } else if b"(){}[],;=+-*/^".contains(&c) {
            t.kind = TokenKind::Symbol;
            t.text = (c as char).to_string();
            i += 1;
        } else {
            let shown = String::from_utf8_lossy(&[c]).into_owned();
            return err(line, format!("unexpected character '{shown}'"));
        }
        out.push(t);
    }
    out.push(Token { kind: TokenKind::End, text: String::new(), number: 0.0, line });
    Ok(out)
}

// ---------------------------------------------------------------- shape tree
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Kind {
    Sphere,
    Box,
    Cylinder,
    Cone,
    Torus,
    Capsule,
    Ellipsoid,
    Polygon,
    HalfSpace,
    Slab,
    Facets,
    Union,
    Intersection,
    Difference,
    SmoothUnion,
    SmoothIntersection,
    SmoothDifference,
    Transform,
    Shell,
    Offset,
    LinearArray,
    PolarArray,
    GridArray,
    Region,
}

#[derive(Clone, Debug)]
struct Node {
    kind: Kind,
    line: i32,
    size: DVec3,
    r: f64,
    r2: f64,
    h: f64,
    k: f64,
    polygon: Vec<DVec2>,
    normals: Vec<DVec3>,
    offsets: Vec<f64>,
    forward: DMat3,
    inverse: DMat3,
    translation: DVec3,
    distance_scale: f64,
    step: DVec3,
    counts: [i32; 3],
    axis: DVec3,
    angle: f64,
    element: String,
    children: Vec<Node>,
}

impl Default for Node {
    fn default() -> Self {
        Self {
            kind: Kind::Union,
            line: 0,
            size: DVec3::ZERO,
            r: 0.0,
            r2: 0.0,
            h: 0.0,
            k: 0.0,
            polygon: Vec::new(),
            normals: Vec::new(),
            offsets: Vec::new(),
            forward: DMat3::IDENTITY,
            inverse: DMat3::IDENTITY,
            translation: DVec3::ZERO,
            distance_scale: 1.0,
            step: DVec3::ZERO,
            counts: [1, 1, 1],
            axis: DVec3::new(0.0, 0.0, 1.0),
            angle: 360.0,
            element: String::new(),
            children: Vec::new(),
        }
    }
}

fn max_component(v: DVec3) -> f64 {
    max_d(v.x, max_d(v.y, v.z))
}

/// The right-handed turn by `angle` (radians) about a unit axis.
fn rotation(axis: DVec3, angle: f64) -> DMat3 {
    let (c, s) = (angle.cos(), angle.sin());
    let mut r = [[0.0f64; 3]; 3]; // r[column][row]
    for (i, col) in r.iter_mut().enumerate() {
        col[i] = c;
    }
    for i in 0..3 {
        for j in 0..3 {
            r[j][i] += (1.0 - c) * axis[i] * axis[j];
        }
    }
    r[1][0] -= s * axis.z;
    r[0][1] += s * axis.z;
    r[2][0] += s * axis.y;
    r[0][2] -= s * axis.y;
    r[2][1] -= s * axis.x;
    r[1][2] += s * axis.x;
    DMat3::from_cols_array_2d(&r)
}

fn polygon_distance(v: &[DVec2], p: DVec2) -> f64 {
    let mut d = (p - v[0]).dot(p - v[0]);
    let mut s = 1.0;
    let mut j = v.len() - 1;
    for i in 0..v.len() {
        let e = v[j] - v[i];
        let w = p - v[i];
        let b = w - e * clamp_d(w.dot(e) / e.dot(e), 0.0, 1.0);
        d = min_d(d, b.dot(b));
        let c1 = p.y >= v[i].y;
        let c2 = p.y < v[j].y;
        let c3 = e.x * w.y > e.y * w.x;
        if (c1 && c2 && c3) || (!c1 && !c2 && !c3) {
            s = -s;
        }
        j = i;
    }
    s * d.sqrt()
}

fn smooth_min(a: f64, b: f64, k: f64) -> f64 {
    if k <= 0.0 {
        return min_d(a, b);
    }
    let h = clamp_d(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    b + (a - b) * h - k * h * (1.0 - h)
}

/// Rodrigues: p turned by `turn` (radians) about the unit axis.
fn turned(p: DVec3, axis: DVec3, turn: f64) -> DVec3 {
    let (c, s) = (turn.cos(), turn.sin());
    p * c + axis.cross(p) * s + axis * axis.dot(p) * (1.0 - c)
}

/// The copies of an array (grid arrays: only those near p).
fn for_each_copy(n: &Node, p: DVec3, mut visit: impl FnMut(DVec3)) {
    match n.kind {
        Kind::LinearArray => {
            for i in 0..n.counts[0] {
                visit(p - n.step * i as f64);
            }
        }
        Kind::PolarArray => {
            let copies = n.counts[0];
            let full = (n.angle.abs() - 360.0).abs() < 1e-9;
            for i in 0..copies {
                let turn = if copies <= 1 { 0.0 } else { n.angle * DEGREE * i as f64 / (if full { copies } else { copies - 1 }) as f64 };
                visit(turned(p, n.axis, -turn));
            }
        }
        _ => {
            let (mut low, mut high) = ([0i32; 3], [0i32; 3]);
            for a in 0..3 {
                if n.counts[a] <= 1 || n.step[a].abs() < 1e-12 {
                    low[a] = 0;
                    high[a] = n.counts[a] - 1;
                    continue;
                }
                let nearest = (p[a] / n.step[a]).round() as i32;
                low[a] = (nearest - 1).clamp(0, n.counts[a] - 1);
                high[a] = (nearest + 1).clamp(0, n.counts[a] - 1);
            }
            for i in low[0]..=high[0] {
                for j in low[1]..=high[1] {
                    for k in low[2]..=high[2] {
                        visit(p - DVec3::new(i as f64 * n.step.x, j as f64 * n.step.y, k as f64 * n.step.z));
                    }
                }
            }
        }
    }
}

fn union_of(children: &[Node], p: DVec3) -> f64 {
    let mut d = HUGE;
    for c in children {
        d = min_d(d, distance_of(c, p));
    }
    d
}

fn distance_of(n: &Node, p: DVec3) -> f64 {
    match n.kind {
        Kind::Sphere => len3(p) - n.r,
        Kind::Box => {
            let q = p.abs() - n.size * 0.5;
            len3(vmax(q, DVec3::ZERO)) + min_d(max_component(q), 0.0)
        }
        Kind::Cylinder => {
            let d = DVec2::new(len2(DVec2::new(p.x, p.y)) - n.r, p.z.abs() - n.h / 2.0);
            min_d(max_d(d.x, d.y), 0.0) + len2(DVec2::new(max_d(d.x, 0.0), max_d(d.y, 0.0)))
        }
        Kind::Cone => {
            // Radius r at z = -h/2, r2 at z = +h/2.
            let (hh, r1, r2) = (n.h / 2.0, n.r, n.r2);
            let q = DVec2::new(len2(DVec2::new(p.x, p.y)), p.z);
            let k1 = DVec2::new(r2, hh);
            let k2 = DVec2::new(r2 - r1, 2.0 * hh);
            let ca = DVec2::new(q.x - min_d(q.x, if q.y < 0.0 { r1 } else { r2 }), q.y.abs() - hh);
            let cb = q - k1 + k2 * clamp_d((k1 - q).dot(k2) / k2.dot(k2), 0.0, 1.0);
            let s = if cb.x < 0.0 && ca.y < 0.0 { -1.0 } else { 1.0 };
            s * min_d(ca.dot(ca), cb.dot(cb)).sqrt()
        }
        Kind::Torus => len2(DVec2::new(len2(DVec2::new(p.x, p.y)) - n.r, p.z)) - n.r2,
        Kind::Capsule => {
            let mut q = p;
            q.z -= clamp_d(q.z, -n.h / 2.0, n.h / 2.0);
            len3(q) - n.r
        }
        Kind::Ellipsoid => {
            let k0 = len3(p / n.size);
            if k0 < 1e-12 {
                return -min_d(n.size.x, min_d(n.size.y, n.size.z));
            }
            let k1 = len3(p / (n.size * n.size));
            k0 * (k0 - 1.0) / k1
        }
        Kind::Polygon => {
            let w = DVec2::new(polygon_distance(&n.polygon, DVec2::new(p.x, p.y)), p.z.abs() - n.h / 2.0);
            min_d(max_d(w.x, w.y), 0.0) + len2(DVec2::new(max_d(w.x, 0.0), max_d(w.y, 0.0)))
        }
        Kind::HalfSpace => p.dot(n.normals[0]) - n.offsets[0],
        Kind::Slab => (p.dot(n.normals[0]) - n.offsets[0]).abs() - n.h / 2.0,
        Kind::Facets => {
            let mut d = -HUGE;
            for (normal, offset) in n.normals.iter().zip(&n.offsets) {
                d = max_d(d, p.dot(*normal) - offset);
            }
            d
        }
        Kind::Union => union_of(&n.children, p),
        // Regions label atoms of the shape; they add no volume.
        Kind::Region => HUGE,
        Kind::Intersection => {
            let mut d = -HUGE;
            for c in &n.children {
                d = max_d(d, distance_of(c, p));
            }
            if n.children.is_empty() {
                HUGE
            } else {
                d
            }
        }
        Kind::Difference => {
            let Some(first) = n.children.first() else { return HUGE };
            let mut d = distance_of(first, p);
            for c in &n.children[1..] {
                d = max_d(d, -distance_of(c, p));
            }
            d
        }
        Kind::SmoothUnion => {
            let mut d = HUGE;
            for c in &n.children {
                d = if d >= HUGE { distance_of(c, p) } else { smooth_min(d, distance_of(c, p), n.k) };
            }
            d
        }
        Kind::SmoothIntersection => {
            let mut d = HUGE;
            for c in &n.children {
                d = if d >= HUGE { -distance_of(c, p) } else { smooth_min(d, -distance_of(c, p), n.k) };
            }
            if n.children.is_empty() {
                HUGE
            } else {
                -d
            }
        }
        Kind::SmoothDifference => {
            let Some(first) = n.children.first() else { return HUGE };
            let mut d = -distance_of(first, p);
            for c in &n.children[1..] {
                d = smooth_min(d, distance_of(c, p), n.k);
            }
            -d
        }
        Kind::Transform => union_of(&n.children, mul_d(&n.inverse, p - n.translation)) * n.distance_scale,
        Kind::Shell => {
            let d = union_of(&n.children, p);
            max_d(d, -(d + n.h))
        }
        Kind::Offset => union_of(&n.children, p) - n.r,
        Kind::LinearArray | Kind::PolarArray | Kind::GridArray => {
            let mut d = HUGE;
            for_each_copy(n, p, |q| d = min_d(d, union_of(&n.children, q)));
            d
        }
    }
}

fn material_of<'a>(n: &'a Node, p: DVec3, out: &mut Option<&'a str>) {
    match n.kind {
        Kind::Region => {
            if union_of(&n.children, p) <= 0.0 {
                *out = Some(&n.element);
            }
            for c in &n.children {
                material_of(c, p, out);
            }
        }
        Kind::Transform => {
            let q = mul_d(&n.inverse, p - n.translation);
            for c in &n.children {
                material_of(c, q, out);
            }
        }
        Kind::LinearArray | Kind::PolarArray | Kind::GridArray => {
            let mut copies = Vec::new();
            for_each_copy(n, p, |q| copies.push(q));
            for q in copies {
                for c in &n.children {
                    material_of(c, q, out);
                }
            }
        }
        _ => {
            for c in &n.children {
                material_of(c, p, out);
            }
        }
    }
}

// ---------------------------------------------------------------- bounds
fn merged(mut a: Bounds, b: Bounds) -> Bounds {
    if a.empty {
        return b;
    }
    if b.empty {
        return a;
    }
    if !a.finite || !b.finite {
        return Bounds::default();
    }
    a.low = vmin(a.low, b.low);
    a.high = vmax(a.high, b.high);
    a
}

fn bounds_box(half: DVec3) -> Bounds {
    Bounds { low: -half, high: half, finite: true, empty: false }
}

fn mapped(b: Bounds, f: impl Fn(DVec3) -> DVec3) -> Bounds {
    if !b.finite || b.empty {
        return b;
    }
    let mut out = Bounds { low: DVec3::splat(HUGE), high: DVec3::splat(-HUGE), finite: true, empty: false };
    for c in 0..8 {
        let corner = DVec3::new(
            if c & 1 != 0 { b.high.x } else { b.low.x },
            if c & 2 != 0 { b.high.y } else { b.low.y },
            if c & 4 != 0 { b.high.z } else { b.low.z },
        );
        let q = f(corner);
        out.low = vmin(out.low, q);
        out.high = vmax(out.high, q);
    }
    out
}

fn facet_bounds(n: &Node) -> Bounds {
    // Bounded only if the normals surround the origin.
    let m = n.normals.len();
    for t in 0..2000 {
        let z = 1.0 - 2.0 * (t as f64 + 0.5) / 2000.0;
        let rr = max_d(0.0, 1.0 - z * z).sqrt();
        let phi = t as f64 * 2.399963229728653;
        let u = DVec3::new(rr * phi.cos(), rr * phi.sin(), z);
        if n.normals.iter().all(|normal| normal.dot(u) <= 1e-9) {
            return Bounds::default();
        }
    }
    let mut b = Bounds { low: DVec3::splat(HUGE), high: DVec3::splat(-HUGE), finite: false, empty: false };
    for i in 0..m {
        for j in i + 1..m {
            for k in j + 1..m {
                let a = DMat3::from_cols(n.normals[i], n.normals[j], n.normals[k]).transpose();
                if af_core::glm::determinant_d(&a).abs() < 1e-10 {
                    continue;
                }
                let x = mul_d(&inverse_d(&a), DVec3::new(n.offsets[i], n.offsets[j], n.offsets[k]));
                if !(0..m).all(|q| n.normals[q].dot(x) <= n.offsets[q] + 1e-7) {
                    continue;
                }
                b.low = vmin(b.low, x);
                b.high = vmax(b.high, x);
                b.finite = true;
            }
        }
    }
    b
}

fn bounds_of(n: &Node) -> Bounds {
    match n.kind {
        Kind::Sphere => bounds_box(DVec3::splat(n.r)),
        Kind::Box => bounds_box(n.size * 0.5),
        Kind::Cylinder => bounds_box(DVec3::new(n.r, n.r, n.h / 2.0)),
        Kind::Cone => {
            let r = max_d(n.r, n.r2);
            bounds_box(DVec3::new(r, r, n.h / 2.0))
        }
        Kind::Torus => bounds_box(DVec3::new(n.r + n.r2, n.r + n.r2, n.r2)),
        Kind::Capsule => bounds_box(DVec3::new(n.r, n.r, n.h / 2.0 + n.r)),
        Kind::Ellipsoid => bounds_box(n.size),
        Kind::Polygon => {
            let mut b = bounds_box(DVec3::new(0.0, 0.0, n.h / 2.0));
            b.low.x = HUGE;
            b.low.y = HUGE;
            b.high.x = -HUGE;
            b.high.y = -HUGE;
            for v in &n.polygon {
                b.low.x = min_d(b.low.x, v.x);
                b.low.y = min_d(b.low.y, v.y);
                b.high.x = max_d(b.high.x, v.x);
                b.high.y = max_d(b.high.y, v.y);
            }
            b
        }
        Kind::HalfSpace | Kind::Slab => Bounds::default(),
        Kind::Facets => facet_bounds(n),
        Kind::Union => bounds_of_union(&n.children),
        Kind::Region => Bounds { finite: true, empty: true, ..Default::default() },
        Kind::Intersection | Kind::SmoothIntersection => {
            let mut b = Bounds::default();
            for c in &n.children {
                let cb = bounds_of(c);
                if !cb.finite || cb.empty {
                    continue;
                }
                if !b.finite {
                    b = cb;
                } else {
                    b.low = vmax(b.low, cb.low);
                    b.high = vmin(b.high, cb.high);
                }
            }
            if b.finite {
                b.high = vmax(b.high, b.low);
            }
            b
        }
        Kind::Difference | Kind::SmoothDifference => n.children.first().map_or(Bounds::default(), bounds_of),
        Kind::SmoothUnion => {
            let mut b = bounds_of_union(&n.children);
            if b.finite {
                b.low -= DVec3::splat(n.k);
                b.high += DVec3::splat(n.k);
            }
            b
        }
        Kind::Transform => mapped(bounds_of_union(&n.children), |q| mul_d(&n.forward, q) + n.translation),
        Kind::Shell => bounds_of_union(&n.children),
        Kind::Offset => {
            let mut b = bounds_of_union(&n.children);
            if b.finite && n.r > 0.0 {
                b.low -= DVec3::splat(n.r);
                b.high += DVec3::splat(n.r);
            }
            b
        }
        Kind::LinearArray => {
            let b = bounds_of_union(&n.children);
            merged(b, mapped(b, |q| q + n.step * (n.counts[0] - 1) as f64))
        }
        Kind::GridArray => {
            let b = bounds_of_union(&n.children);
            let last = DVec3::new(
                n.step.x * (n.counts[0] - 1) as f64,
                n.step.y * (n.counts[1] - 1) as f64,
                n.step.z * (n.counts[2] - 1) as f64,
            );
            merged(b, mapped(b, |q| q + last))
        }
        Kind::PolarArray => {
            let child = bounds_of_union(&n.children);
            if !child.finite {
                return child;
            }
            let mut b = child;
            let full = (n.angle.abs() - 360.0).abs() < 1e-9;
            for i in 1..n.counts[0] {
                let turn = n.angle * DEGREE * i as f64 / (if full { n.counts[0] } else { n.counts[0] - 1 }) as f64;
                b = merged(b, mapped(child, |p| turned(p, n.axis, turn)));
            }
            b
        }
    }
}

fn bounds_of_union(children: &[Node]) -> Bounds {
    if children.is_empty() {
        return Bounds::default();
    }
    let mut b = Bounds { finite: true, empty: true, ..Default::default() };
    for c in children {
        b = merged(b, bounds_of(c));
    }
    b
}

/// Whether a node adds volume (anything but regions).
fn has_volume(n: &Node) -> bool {
    match n.kind {
        Kind::Region => false,
        Kind::Union | Kind::Transform | Kind::LinearArray | Kind::PolarArray | Kind::GridArray => n.children.iter().any(has_volume),
        _ => true,
    }
}

fn collect_regions(n: &Node, out: &mut Vec<String>) {
    if n.kind == Kind::Region && !out.contains(&n.element) {
        out.push(n.element.clone());
    }
    for c in &n.children {
        collect_regions(c, out);
    }
}

// ---------------------------------------------------------------- parser
struct Arguments<'a> {
    command: &'a Token,
    values: Vec<(String, Value)>,
    /// Names of positional arguments, in order.
    positional: Vec<&'static str>,
    used: Vec<String>,
}

impl Arguments<'_> {
    fn find(&mut self, name: &str) -> Option<Value> {
        let mut p = 0;
        for (key, value) in &self.values {
            let resolved = if key.is_empty() { self.positional.get(p).map_or("", |s| *s) } else { key.as_str() };
            if key.is_empty() {
                p += 1;
            }
            if resolved == name {
                if !self.used.iter().any(|u| u == name) {
                    self.used.push(name.to_string());
                }
                return Some(value.clone());
            }
        }
        None
    }
    fn has(&mut self, name: &str) -> bool {
        self.find(name).is_some()
    }
    fn number(&mut self, name: &str, fallback: f64) -> Result<f64, ShapeError> {
        match self.find(name) {
            None => Ok(fallback),
            Some(Value::Number(x)) => Ok(x),
            Some(Value::Boolean(b)) => Ok(if b { 1.0 } else { 0.0 }),
            Some(v) => err(self.command.line, format!("{}(): {name} must be a number, not {}", self.command.text, describe(&v))),
        }
    }
    fn positive_number(&mut self, name: &str, fallback: f64) -> Result<f64, ShapeError> {
        let x = self.number(name, fallback)?;
        if !(x > 0.0) {
            return err(self.command.line, format!("{}(): {name} must be positive", self.command.text));
        }
        Ok(x)
    }
    fn vector(&mut self, name: &str, fallback: DVec3, scalar_fills: bool) -> Result<DVec3, ShapeError> {
        let Some(v) = self.find(name) else { return Ok(fallback) };
        if let (Value::Number(x), true) = (&v, scalar_fills) {
            return Ok(DVec3::splat(*x));
        }
        match v {
            Value::Vector(items) if (2..=3).contains(&items.len()) => {
                let mut out = DVec3::ZERO;
                for (i, item) in items.iter().enumerate() {
                    out[i] = number(item, self.command.line)?;
                }
                Ok(out)
            }
            _ => err(self.command.line, format!("{}(): {name} must be a vector [x, y, z]", self.command.text)),
        }
    }
    fn text(&mut self, name: &str, fallback: &str) -> Result<String, ShapeError> {
        match self.find(name) {
            None => Ok(fallback.to_string()),
            Some(Value::Text(t)) => Ok(t),
            Some(_) => err(self.command.line, format!("{}(): {name} must be text in quotes", self.command.text)),
        }
    }
    fn finish(&self) -> Result<(), ShapeError> {
        let mut p = 0;
        for (key, _) in &self.values {
            if key.is_empty() {
                if p >= self.positional.len() {
                    return err(
                        self.command.line,
                        format!("{}() takes at most {} unnamed argument(s)", self.command.text, self.positional.len()),
                    );
                }
                p += 1;
                continue;
            }
            if !self.used.contains(key) {
                return err(self.command.line, format!("{}() has no argument '{key}'", self.command.text));
            }
        }
        Ok(())
    }
}

fn number(v: &Value, line: i32) -> Result<f64, ShapeError> {
    match v {
        Value::Number(x) => Ok(*x),
        Value::Boolean(b) => Ok(if *b { 1.0 } else { 0.0 }),
        _ => err(line, format!("expected a number, found {}", describe(v))),
    }
}

fn arithmetic(a: &Value, b: &Value, op: u8, line: i32) -> Result<Value, ShapeError> {
    let (av, bv) = (matches!(a, Value::Vector(_)), matches!(b, Value::Vector(_)));
    if !av && !bv {
        let (x, y) = (number(a, line)?, number(b, line)?);
        return Ok(Value::Number(match op {
            b'+' => x + y,
            b'-' => x - y,
            b'*' => x * y,
            _ => {
                if y == 0.0 {
                    return err(line, "division by zero");
                }
                x / y
            }
        }));
    }
    if let (Value::Vector(ai), Value::Vector(bi)) = (a, b) {
        if op != b'+' && op != b'-' {
            return err(line, "vectors can only be added or subtracted");
        }
        if ai.len() != bi.len() {
            return err(line, "vectors of different lengths");
        }
        return Ok(Value::Vector(ai.iter().zip(bi).map(|(x, y)| arithmetic(x, y, op, line)).collect::<Result<_, _>>()?));
    }
    let (vector, scalar) = if av { (a, b) } else { (b, a) };
    if op == b'+' || op == b'-' || (op == b'/' && !av) {
        return err(line, "a vector and a number can only be multiplied or divided");
    }
    let Value::Vector(items) = vector else { unreachable!() };
    let items = items
        .iter()
        .map(|item| if av { arithmetic(item, scalar, op, line) } else { arithmetic(scalar, item, op, line) })
        .collect::<Result<_, _>>()?;
    Ok(Value::Vector(items))
}

fn function(t: &Token, a: &[f64]) -> Result<f64, ShapeError> {
    let need = |n: usize| -> Result<(), ShapeError> {
        if a.len() != n {
            return err(t.line, format!("{}() takes {n} argument{}", t.text, if n == 1 { "" } else { "s" }));
        }
        Ok(())
    };
    let f = t.text.as_str();
    Ok(match f {
        "sqrt" => {
            need(1)?;
            a[0].sqrt()
        }
        "sin" => {
            need(1)?;
            (a[0] * DEGREE).sin()
        }
        "cos" => {
            need(1)?;
            (a[0] * DEGREE).cos()
        }
        "tan" => {
            need(1)?;
            (a[0] * DEGREE).tan()
        }
        "asin" => {
            need(1)?;
            a[0].asin() / DEGREE
        }
        "acos" => {
            need(1)?;
            a[0].acos() / DEGREE
        }
        "atan" => {
            need(1)?;
            a[0].atan() / DEGREE
        }
        "atan2" => {
            need(2)?;
            a[0].atan2(a[1]) / DEGREE
        }
        "abs" => {
            need(1)?;
            a[0].abs()
        }
        "floor" => {
            need(1)?;
            a[0].floor()
        }
        "ceil" => {
            need(1)?;
            a[0].ceil()
        }
        "round" => {
            need(1)?;
            a[0].round()
        }
        "exp" => {
            need(1)?;
            a[0].exp()
        }
        "ln" | "log" => {
            need(1)?;
            a[0].ln()
        }
        "pow" => {
            need(2)?;
            a[0].powf(a[1])
        }
        "min" | "max" => {
            if a.is_empty() {
                return err(t.line, format!("{f}() needs values"));
            }
            // min_element / max_element: the first extreme.
            let mut best = a[0];
            for &x in &a[1..] {
                if (f == "min" && x < best) || (f == "max" && best < x) {
                    best = x;
                }
            }
            best
        }
        _ => return err(t.line, format!("unknown function '{f}'")),
    })
}

struct Parser {
    tokens: Vec<Token>,
    at: usize,
    cell: Option<DMat3>,
    variables: std::collections::BTreeMap<String, Value>,
    primitives: i32,
}

impl Parser {
    fn peek(&self, ahead: usize) -> &Token {
        &self.tokens[(self.at + ahead).min(self.tokens.len() - 1)]
    }
    fn next(&mut self) -> Token {
        let t = self.tokens[self.at.min(self.tokens.len() - 1)].clone();
        self.at += 1;
        t
    }
    fn symbol(&self, s: &str) -> bool {
        let t = self.peek(0);
        t.kind == TokenKind::Symbol && t.text == s
    }
    fn found(&self) -> String {
        let t = self.peek(0);
        if t.kind == TokenKind::End {
            " at the end of the script".into()
        } else {
            format!(" before '{}'", t.text)
        }
    }
    fn expect(&mut self, s: &str) -> Result<(), ShapeError> {
        if !self.symbol(s) {
            return err(self.peek(0).line, format!("expected '{s}'{}", self.found()));
        }
        self.at += 1;
        Ok(())
    }

    /// statement: name = expr ; | command ( args ) ( { statements } | statement | ; )
    fn statement(&mut self, out: &mut Vec<Node>) -> Result<(), ShapeError> {
        if self.symbol(";") {
            self.at += 1;
            return Ok(());
        }
        let name = self.next();
        if name.kind != TokenKind::Identifier {
            return err(name.line, format!("expected a command or a variable, found '{}'", name.text));
        }
        if self.symbol("=") {
            self.at += 1;
            let v = self.expression()?;
            self.variables.insert(name.text.clone(), v);
            if self.symbol(";") {
                self.at += 1;
            }
            return Ok(());
        }
        self.expect("(")?;
        let mut arguments = Vec::new();
        while !self.symbol(")") {
            if self.peek(0).kind == TokenKind::End {
                return err(name.line, format!("'{}(' is not closed with ')'", name.text));
            }
            let mut key = String::new();
            if self.peek(0).kind == TokenKind::Identifier && self.peek(1).kind == TokenKind::Symbol && self.peek(1).text == "=" {
                key = self.next().text;
                self.at += 1;
            }
            let v = self.expression()?;
            arguments.push((key, v));
            if self.peek(0).kind == TokenKind::End {
                return err(name.line, format!("'{}(' is not closed with ')'", name.text));
            }
            if !self.symbol(")") && !self.symbol(",") {
                return err(self.peek(0).line, format!("expected ',' or ')' in {}(){}", name.text, self.found()));
            }
            if self.symbol(",") {
                self.at += 1;
            }
        }
        self.expect(")")?;
        let mut node = self.make(&name, arguments)?;
        if node.kind >= Kind::Union {
            if self.symbol("{") {
                self.at += 1;
                while !self.symbol("}") {
                    if self.peek(0).kind == TokenKind::End {
                        return err(name.line, format!("'{}' block is not closed with '}}'", name.text));
                    }
                    self.statement(&mut node.children)?;
                }
                self.at += 1;
            } else if !self.symbol(";") && self.peek(0).kind != TokenKind::End {
                self.statement(&mut node.children)?; // one child: translate([..]) sphere(1);
            } else if self.symbol(";") {
                self.at += 1;
            }
            if node.children.is_empty() {
                return err(name.line, format!("'{}' needs shapes to act on, in {{ }} or after it", name.text));
            }
            if node.kind == Kind::Transform && node.distance_scale <= 0.0 {
                return err(name.line, "invalid transform");
            }
        } else if self.symbol(";") {
            self.at += 1;
        }
        out.push(node);
        Ok(())
    }

    fn expression(&mut self) -> Result<Value, ShapeError> {
        let mut v = self.term()?;
        while self.symbol("+") || self.symbol("-") {
            let op = self.next();
            let rhs = self.term()?;
            v = arithmetic(&v, &rhs, op.text.as_bytes()[0], op.line)?;
        }
        Ok(v)
    }

    fn term(&mut self) -> Result<Value, ShapeError> {
        let mut v = self.unary()?;
        while self.symbol("*") || self.symbol("/") {
            let op = self.next();
            let rhs = self.unary()?;
            v = arithmetic(&v, &rhs, op.text.as_bytes()[0], op.line)?;
        }
        Ok(v)
    }

    fn unary(&mut self) -> Result<Value, ShapeError> {
        if self.symbol("-") {
            let line = self.next().line;
            let v = self.unary()?;
            return arithmetic(&Value::Number(-1.0), &v, b'*', line);
        }
        if self.symbol("+") {
            self.at += 1;
            return self.unary();
        }
        let base = self.primary()?;
        if self.symbol("^") {
            let line = self.next().line;
            let exponent = self.unary()?;
            return Ok(Value::Number(number(&base, line)?.powf(number(&exponent, line)?)));
        }
        Ok(base)
    }

    fn primary(&mut self) -> Result<Value, ShapeError> {
        let t = self.next();
        match t.kind {
            TokenKind::Number => return Ok(Value::Number(t.number)),
            TokenKind::Text => return Ok(Value::Text(t.text)),
            TokenKind::Symbol if t.text == "(" => {
                let v = self.expression()?;
                self.expect(")")?;
                return Ok(v);
            }
            TokenKind::Symbol if t.text == "[" => {
                let mut items = Vec::new();
                while !self.symbol("]") {
                    items.push(self.expression()?);
                    if !self.symbol("]") {
                        self.expect(",")?;
                    }
                }
                self.at += 1;
                return Ok(Value::Vector(items));
            }
            TokenKind::Identifier => {
                if self.symbol("(") {
                    self.at += 1;
                    let mut args = Vec::new();
                    while !self.symbol(")") {
                        let v = self.expression()?;
                        args.push(number(&v, t.line)?);
                        if !self.symbol(")") {
                            self.expect(",")?;
                        }
                    }
                    self.at += 1;
                    return Ok(Value::Number(function(&t, &args)?));
                }
                return match self.variables.get(&t.text) {
                    Some(v) => Ok(v.clone()),
                    None => err(t.line, format!("unknown variable '{}' (set it first: {} = 10;)", t.text, t.text)),
                };
            }
            _ => {}
        }
        err(t.line, format!("expected a value, found '{}'", t.text))
    }

    /// The normal of {hkl} (reciprocal lattice) or a Cartesian direction.
    fn plane_normal(&self, v: DVec3, miller: bool) -> DVec3 {
        match (&self.cell, miller) {
            (Some(cell), true) => mul_d(&inverse_d(cell).transpose(), v),
            _ => v,
        }
    }

    fn make(&mut self, command: &Token, values: Vec<(String, Value)>) -> Result<Node, ShapeError> {
        let mut n = Node { line: command.line, ..Default::default() };
        let mut a = Arguments { command, values, positional: Vec::new(), used: Vec::new() };
        let c = command.text.as_str();
        let line = command.line;
        let radius = |a: &mut Arguments, name: &str, diameter: &str, fallback: f64| -> Result<f64, ShapeError> {
            if a.has(diameter) {
                return Ok(a.positive_number(diameter, 2.0 * fallback)? / 2.0);
            }
            a.positive_number(name, fallback)
        };
        match c {
            "sphere" => {
                a.positional = vec!["r"];
                n.kind = Kind::Sphere;
                n.r = radius(&mut a, "r", "d", 10.0)?;
            }
            "box" | "cube" => {
                a.positional = vec!["size", "center"];
                n.kind = Kind::Box;
                n.size = a.vector("size", DVec3::splat(10.0), true)?;
                if !(n.size.x > 0.0 && n.size.y > 0.0 && n.size.z > 0.0) {
                    return err(line, format!("{c}(): size must be positive"));
                }
                if a.number("center", 1.0)? == 0.0 {
                    let shift = n.size * 0.5;
                    return self.corner(n, shift, &a);
                }
            }
            "cylinder" => {
                a.positional = vec!["h", "r", "r2", "center"];
                n.h = a.positive_number("h", 10.0)?;
                if a.has("r1") || a.has("r2") || a.has("d1") || a.has("d2") {
                    n.kind = Kind::Cone;
                    n.r = if a.has("d1") {
                        a.number("d1", 0.0)? / 2.0
                    } else {
                        let r = a.number("r", 5.0)?;
                        a.number("r1", r)?
                    };
                    n.r2 = if a.has("d2") {
                        a.number("d2", 0.0)? / 2.0
                    } else {
                        let r = a.number("r", 5.0)?;
                        a.number("r2", r)?
                    };
                    if n.r < 0.0 || n.r2 < 0.0 || (n.r == 0.0 && n.r2 == 0.0) {
                        return err(line, "cylinder(): radii must not be negative or both zero");
                    }
                } else {
                    n.kind = Kind::Cylinder;
                    n.r = radius(&mut a, "r", "d", 5.0)?;
                }
                if a.number("center", 1.0)? == 0.0 {
                    let shift = DVec3::new(0.0, 0.0, n.h / 2.0);
                    return self.corner(n, shift, &a);
                }
            }
            "cone" => {
                a.positional = vec!["r1", "r2", "h"];
                n.kind = Kind::Cone;
                n.r = a.number("r1", 5.0)?;
                n.r2 = a.number("r2", 0.0)?;
                n.h = a.positive_number("h", 10.0)?;
                if n.r < 0.0 || n.r2 < 0.0 || (n.r == 0.0 && n.r2 == 0.0) {
                    return err(line, "cone(): radii must not be negative or both zero");
                }
            }
            "torus" => {
                a.positional = vec!["R", "r"];
                n.kind = Kind::Torus;
                n.r = a.positive_number("R", 10.0)?;
                n.r2 = a.positive_number("r", 3.0)?;
            }
            "capsule" => {
                a.positional = vec!["r", "h"];
                n.kind = Kind::Capsule;
                n.r = a.positive_number("r", 5.0)?;
                n.h = a.number("h", 10.0)?;
                if n.h < 0.0 {
                    return err(line, "capsule(): h must not be negative");
                }
            }
            "ellipsoid" => {
                a.positional = vec!["r"];
                n.kind = Kind::Ellipsoid;
                n.size = a.vector("r", DVec3::new(10.0, 6.0, 4.0), true)?;
                if !(n.size.x > 0.0 && n.size.y > 0.0 && n.size.z > 0.0) {
                    return err(line, "ellipsoid(): radii must be positive");
                }
            }
            "prism" => {
                // A regular polygon of n sides (circumradius r, or apothem) extruded along z.
                a.positional = vec!["n", "r", "h"];
                n.kind = Kind::Polygon;
                let sides = a.number("n", 6.0)? as i32;
                if !(3..=1000).contains(&sides) {
                    return err(line, "prism(): n must be 3 or more sides");
                }
                let r = if a.has("apothem") { a.positive_number("apothem", 1.0)? / (PI / sides as f64).cos() } else { a.positive_number("r", 10.0)? };
                n.h = a.positive_number("h", 20.0)?;
                let start = a.number("rotation", 0.0)? * DEGREE;
                for i in 0..sides {
                    let angle = start + 2.0 * PI * i as f64 / sides as f64;
                    n.polygon.push(r * DVec2::new(angle.cos(), angle.sin()));
                }
            }
            "polygon" | "extrude" => {
                a.positional = vec!["points", "h"];
                n.kind = Kind::Polygon;
                let points = match a.find("points") {
                    Some(Value::Vector(items)) if items.len() >= 3 => items,
                    _ => return err(line, format!("{c}(): points must list at least three [x, y] corners")),
                };
                for p in &points {
                    match p {
                        Value::Vector(xy) if xy.len() == 2 => n.polygon.push(DVec2::new(number(&xy[0], line)?, number(&xy[1], line)?)),
                        _ => return err(line, format!("{c}(): each point is [x, y]")),
                    }
                }
                n.h = a.positive_number("h", 10.0)?;
            }
            "halfspace" | "plane" | "slab" => {
                let slab = c == "slab";
                a.positional = if slab { vec!["normal", "thickness", "d"] } else { vec!["normal", "d"] };
                n.kind = if slab { Kind::Slab } else { Kind::HalfSpace };
                let miller = a.has("hkl");
                let v = if miller { a.vector("hkl", DVec3::new(0.0, 0.0, 1.0), false)? } else { a.vector("normal", DVec3::new(0.0, 0.0, 1.0), false)? };
                let direction = self.plane_normal(v, miller);
                if len3(direction) < 1e-12 {
                    return err(line, if slab { "slab(): the normal must not be zero".to_string() } else { format!("{c}(): the normal must not be zero") });
                }
                n.normals = vec![normalize_d(direction)];
                if slab {
                    n.h = a.positive_number("thickness", 10.0)?;
                }
                n.offsets = vec![a.number("d", 0.0)?];
            }
            "facets" | "polyhedron" => {
                // Planes {hkl} at distances d (A); with family, every permutation and sign.
                a.positional = vec!["hkl", "d"];
                n.kind = Kind::Facets;
                let hkl = match a.find("hkl") {
                    Some(Value::Vector(items)) if !items.is_empty() => items,
                    _ => return err(line, format!("{c}(): hkl must be [h, k, l] or a list of them")),
                };
                let mut planes = Vec::new();
                if matches!(hkl[0], Value::Vector(_)) {
                    for p in &hkl {
                        match p {
                            Value::Vector(v) if v.len() == 3 => planes.push(DVec3::new(number(&v[0], line)?, number(&v[1], line)?, number(&v[2], line)?)),
                            _ => return err(line, format!("{c}(): each plane is [h, k, l]")),
                        }
                    }
                } else {
                    if hkl.len() != 3 {
                        return err(line, format!("{c}(): hkl is [h, k, l]"));
                    }
                    planes.push(DVec3::new(number(&hkl[0], line)?, number(&hkl[1], line)?, number(&hkl[2], line)?));
                }
                let distances: Vec<f64> = match a.find("d") {
                    None => return err(line, format!("{c}(): give the facet distances d (Angstrom)")),
                    Some(Value::Vector(items)) => items.iter().map(|x| number(x, line)).collect::<Result<_, _>>()?,
                    Some(v) => vec![number(&v, line)?; planes.len()],
                };
                if distances.len() != planes.len() {
                    return err(line, format!("{c}(): one distance per plane"));
                }
                let family = a.number("family", 1.0)? != 0.0;
                const ORDERS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
                for (plane, &distance) in planes.iter().zip(&distances) {
                    if !(distance > 0.0) {
                        return err(line, format!("{c}(): distances must be positive"));
                    }
                    let mut members = vec![*plane];
                    if family {
                        members.clear();
                        for order in ORDERS {
                            for signs in 0..8 {
                                let mut m = DVec3::new(plane[order[0]], plane[order[1]], plane[order[2]]);
                                for s in 0..3 {
                                    if signs & (1 << s) != 0 {
                                        m[s] = -m[s];
                                    }
                                }
                                if !members.iter().any(|q| len3(*q - m) < 1e-9) {
                                    members.push(m);
                                }
                            }
                        }
                    }
                    for m in members {
                        let direction = self.plane_normal(m, true);
                        if len3(direction) < 1e-12 {
                            return err(line, format!("{c}(): (0 0 0) is not a plane"));
                        }
                        n.normals.push(normalize_d(direction));
                        n.offsets.push(distance);
                    }
                }
            }
            "octahedron" => {
                a.positional = vec!["d"];
                n.kind = Kind::Facets;
                let d = a.positive_number("d", 10.0)?;
                for s in 0..8 {
                    let sign = |bit: i32| if s & bit != 0 { -1.0 } else { 1.0 };
                    n.normals.push(normalize_d(DVec3::new(sign(1), sign(2), sign(4))));
                    n.offsets.push(d);
                }
            }
            "union" | "group" => n.kind = Kind::Union,
            "intersection" | "intersect" => n.kind = Kind::Intersection,
            "difference" | "subtract" => n.kind = Kind::Difference,
            "smooth_union" | "blend" | "smooth_intersection" | "smooth_difference" => {
                a.positional = vec!["k"];
                n.kind = match c {
                    "smooth_intersection" => Kind::SmoothIntersection,
                    "smooth_difference" => Kind::SmoothDifference,
                    _ => Kind::SmoothUnion,
                };
                n.k = a.positive_number("k", 2.0)?;
            }
            "translate" | "move" => {
                a.positional = vec!["v"];
                n.kind = Kind::Transform;
                n.translation = a.vector("v", DVec3::ZERO, false)?;
            }
            "rotate" => {
                // rotate([ax, ay, az]): about x, then y, then z (degrees); rotate(a, v): about an axis.
                a.positional = vec!["a", "v"];
                n.kind = Kind::Transform;
                let Some(angle) = a.find("a") else {
                    return err(line, "rotate(): give the angles [ax, ay, az] or an angle a and an axis v");
                };
                if matches!(angle, Value::Vector(_)) {
                    let e = a.vector("a", DVec3::ZERO, false)? * DEGREE;
                    n.forward = rotation(DVec3::new(0.0, 0.0, 1.0), e.z) * rotation(DVec3::new(0.0, 1.0, 0.0), e.y) * rotation(DVec3::new(1.0, 0.0, 0.0), e.x);
                } else {
                    let t = number(&angle, line)? * DEGREE;
                    let axis = a.vector("v", DVec3::new(0.0, 0.0, 1.0), false)?;
                    if len3(axis) < 1e-12 {
                        return err(line, "rotate(): the axis must not be zero");
                    }
                    n.forward = rotation(normalize_d(axis), t);
                }
                n.inverse = n.forward.transpose();
            }
            "scale" => {
                a.positional = vec!["v"];
                n.kind = Kind::Transform;
                let s = a.vector("v", DVec3::ONE, true)?;
                if s.x == 0.0 || s.y == 0.0 || s.z == 0.0 {
                    return err(line, "scale(): factors must not be zero");
                }
                n.forward = DMat3::from_diagonal(s);
                n.inverse = inverse_d(&n.forward);
                n.distance_scale = min_d(s.x.abs(), min_d(s.y.abs(), s.z.abs()));
            }
            "mirror" => {
                a.positional = vec!["v"];
                n.kind = Kind::Transform;
                let normal = a.vector("v", DVec3::new(1.0, 0.0, 0.0), false)?;
                if len3(normal) < 1e-12 {
                    return err(line, "mirror(): the normal must not be zero");
                }
                let u = normalize_d(normal);
                // I - 2 u u^T (outerProduct(u, u): column i = u * u[i]).
                let outer = DMat3::from_cols(u * u.x, u * u.y, u * u.z);
                n.forward = DMat3::IDENTITY - outer * 2.0;
                n.inverse = n.forward;
            }
            "shell" | "hollow" => {
                a.positional = vec!["t"];
                n.kind = Kind::Shell;
                n.h = a.positive_number("t", 2.0)?;
            }
            "offset" | "grow" => {
                a.positional = vec!["r"];
                n.kind = Kind::Offset;
                n.r = a.number("r", 1.0)?;
            }
            "linear_array" => {
                a.positional = vec!["n", "step"];
                n.kind = Kind::LinearArray;
                n.counts[0] = count(&mut a, "n", 2.0)?;
                n.step = a.vector("step", DVec3::new(10.0, 0.0, 0.0), false)?;
            }
            "polar_array" => {
                a.positional = vec!["n", "angle", "axis"];
                n.kind = Kind::PolarArray;
                n.counts[0] = count(&mut a, "n", 6.0)?;
                n.angle = a.number("angle", 360.0)?;
                n.axis = a.vector("axis", DVec3::new(0.0, 0.0, 1.0), false)?;
                if len3(n.axis) < 1e-12 {
                    return err(line, "polar_array(): the axis must not be zero");
                }
                n.axis = normalize_d(n.axis);
            }
            "grid_array" => {
                a.positional = vec!["n", "step"];
                n.kind = Kind::GridArray;
                let counts = a.vector("n", DVec3::new(2.0, 2.0, 1.0), false)?;
                for i in 0..3 {
                    if counts[i] < 1.0 || counts[i] != counts[i].floor() || counts[i] > 10000.0 {
                        return err(line, "grid_array(): n must be whole numbers from 1");
                    }
                    n.counts[i] = counts[i] as i32;
                }
                n.step = a.vector("step", DVec3::splat(10.0), false)?;
            }
            "region" | "material" => {
                a.positional = vec!["element"];
                n.kind = Kind::Region;
                n.element = a.text("element", "")?;
                if n.element.is_empty() {
                    return err(line, format!("{c}(): give the element, e.g. region(\"Au\") {{ ... }}"));
                }
            }
            _ => return err(line, format!("unknown command '{c}'")),
        }
        a.finish()?;
        if n.kind < Kind::Union {
            self.primitives += 1;
        }
        Ok(n)
    }

    /// center=false: the primitive's corner (or base) at the origin.
    fn corner(&mut self, n: Node, shift: DVec3, a: &Arguments) -> Result<Node, ShapeError> {
        a.finish()?;
        self.primitives += 1;
        Ok(Node { kind: Kind::Transform, line: n.line, translation: shift, children: vec![n], ..Default::default() })
    }
}

fn count(a: &mut Arguments, name: &str, fallback: f64) -> Result<i32, ShapeError> {
    let v = a.number(name, fallback)?;
    if v < 1.0 || v != v.floor() || v > 100000.0 {
        return err(a.command.line, format!("{}(): {name} must be a whole number from 1", a.command.text));
    }
    Ok(v as i32)
}

/// A parsed shape.
#[derive(Clone, Debug, Default)]
pub struct Model {
    root: Option<Box<Node>>,
    primitives: i32,
}

impl Model {
    /// `cell` (columns a, b, c) turns Miller indices of facets into plane
    /// normals; without it they are Cartesian directions.
    pub fn parse(script: &str, cell: Option<&DMat3>) -> Result<Model, ShapeError> {
        let mut parser = Parser { tokens: tokenize(script)?, at: 0, cell: cell.copied(), variables: Default::default(), primitives: 0 };
        parser.variables.insert("pi".into(), Value::Number(PI));
        parser.variables.insert("PI".into(), Value::Number(PI));
        parser.variables.insert("true".into(), Value::Boolean(true));
        parser.variables.insert("false".into(), Value::Boolean(false));
        let mut root = Node { kind: Kind::Union, ..Default::default() };
        while parser.peek(0).kind != TokenKind::End {
            parser.statement(&mut root.children)?;
        }
        if root.children.is_empty() {
            return err(1, "the script has no shape (try: sphere(r = 15);)");
        }
        if !has_volume(&root) {
            return err(1, "the script has only regions: give the shape itself too (regions only assign elements inside it)");
        }
        Ok(Model { root: Some(Box::new(root)), primitives: parser.primitives })
    }

    /// Signed distance (A): negative inside.
    pub fn distance(&self, p: DVec3) -> f64 {
        self.root.as_ref().map_or(HUGE, |r| distance_of(r, p))
    }

    pub fn contains(&self, p: DVec3) -> bool {
        self.distance(p) <= 0.0
    }

    /// The element of the last region containing p ("none": remove the atom).
    pub fn material(&self, p: DVec3) -> Option<&str> {
        let mut out = None;
        if let Some(r) = &self.root {
            material_of(r, p, &mut out);
        }
        out
    }

    /// The shape's surface: marching tetrahedra on the signed distance,
    /// sampled at `spacing` (Angstrom) over the bounds plus two samples.
    pub fn mesh(&self, spacing: f64) -> Result<af_electronic::volume::Mesh, String> {
        let b = self.bounds();
        if !b.finite {
            return Err("The shape is unbounded: intersect it with a box to mesh it".into());
        }
        if !(spacing > 0.0) {
            return Err("The mesh spacing must be positive".into());
        }
        let low = b.low - DVec3::splat(2.0 * spacing);
        let extent = (b.high - b.low) + DVec3::splat(4.0 * spacing);
        let mut grid = af_electronic::volume::Grid::default();
        for a in 0..3 {
            grid.shape[a] = 2.max((extent[a] / spacing).ceil() as i32 + 1);
        }
        let samples = grid.shape[0] as f64 * grid.shape[1] as f64 * grid.shape[2] as f64;
        if samples > 6e7 {
            return Err("The mesh grid is too fine for this shape; use a larger spacing".into());
        }
        grid.cell = glam::DMat3::from_cols(DVec3::new(extent.x, 0.0, 0.0), DVec3::new(0.0, extent.y, 0.0), DVec3::new(0.0, 0.0, extent.z));
        grid.origin = low;
        grid.periodic = false;
        grid.values = vec![0.0; samples as usize];
        for z in 0..grid.shape[2] {
            for y in 0..grid.shape[1] {
                for x in 0..grid.shape[0] {
                    let index = grid.index(x, y, z);
                    grid.values[index] = self.distance(grid.position(x, y, z));
                }
            }
        }
        af_electronic::isosurface::isosurface(&grid, 0.0, None)
    }

    pub fn bounds(&self) -> Bounds {
        let b = self.root.as_ref().map_or(Bounds::default(), |r| bounds_of(r));
        if b.empty {
            Bounds::default()
        } else {
            b
        }
    }

    /// Region elements in script order, without repeats.
    pub fn region_elements(&self) -> Vec<String> {
        let mut out = Vec::new();
        if let Some(r) = &self.root {
            collect_regions(r, &mut out);
        }
        out
    }

    pub fn primitives(&self) -> i32 {
        self.primitives
    }

    pub fn is_empty(&self) -> bool {
        self.root.is_none()
    }

    /// Volume (A^3) counted on a grid of the given spacing.
    pub fn volume(&self, spacing: f64) -> f64 {
        let b = self.bounds();
        if !b.finite || !(spacing > 0.0) {
            return 0.0;
        }
        let mut inside = 0i64;
        let mut x = b.low.x + spacing / 2.0;
        while x < b.high.x {
            let mut y = b.low.y + spacing / 2.0;
            while y < b.high.y {
                let mut z = b.low.z + spacing / 2.0;
                while z < b.high.z {
                    if self.contains(DVec3::new(x, y, z)) {
                        inside += 1;
                    }
                    z += spacing;
                }
                y += spacing;
            }
            x += spacing;
        }
        inside as f64 * spacing * spacing * spacing
    }
}

/// A command for inserting into a script.
pub struct Command {
    pub group: &'static str,
    pub name: &'static str,
    pub snippet: &'static str,
    pub help: &'static str,
}

/// Ready-made shapes: (title, script).
pub fn templates() -> &'static [(&'static str, &'static str)] {
    &[
        ("Sphere", "sphere(r = 15);\n"),
        ("Core-shell particle", "// A particle of the crystal with a core of another element.\nR = 20;\nsphere(r = R);\nregion(\"Au\") sphere(r = 0.6 * R);\n"),
        ("Hollow sphere (nanoshell)", "shell(t = 5) sphere(r = 20);\n"),
        ("Hollow cylinder (tube)", "difference() {\n    cylinder(r = 15, h = 60);\n    cylinder(r = 10, h = 62);\n}\n"),
        ("Truncated octahedron {111} + {100}", "// Distances of the facets from the centre (Angstrom); family: all equivalent planes (cubic).\nfacets(hkl = [[1, 1, 1], [1, 0, 0]], d = [18, 20], family = true);\n"),
        ("Faceted rod {110} x {001}", "intersection() {\n    facets(hkl = [[1, 1, 0], [1, -1, 0]], d = [15, 15], family = false) ;\n    facets(hkl = [[-1, -1, 0], [-1, 1, 0]], d = [15, 15], family = false);\n    slab(normal = [0, 0, 1], thickness = 60);\n}\n"),
        ("Hexagonal nanowire segment", "prism(n = 6, r = 15, h = 60);\n"),
        ("Cube with a spherical void", "difference() {\n    box(size = [40, 40, 40]);\n    sphere(r = 10);\n}\n"),
        ("Nanoporous block (pore lattice)", "difference() {\n    box(size = [60, 60, 30]);\n    translate([-22.5, -22.5, 0])\n        grid_array(n = [4, 4, 1], step = [15, 15, 0]) cylinder(r = 4, h = 32);\n}\n"),
        ("Dumbbell (sintered particles)", "smooth_union(k = 4) {\n    translate([-12, 0, 0]) sphere(r = 10);\n    translate([12, 0, 0]) sphere(r = 10);\n}\n"),
        ("Nanoring (torus)", "torus(R = 25, r = 6);\n"),
        ("Janus particle", "R = 20;\nsphere(r = R);\n// The upper half (z >= 0) becomes silver.\nregion(\"Ag\") halfspace(normal = [0, 0, -1], d = 0);\n"),
        ("Probe tip (cone on a shank)", "union() {\n    translate([0, 0, 15]) cone(r1 = 15, r2 = 2, h = 30);\n    translate([0, 0, -10]) cylinder(r = 15, h = 20);\n}\n"),
        ("Matrix with an ellipsoidal precipitate", "box(size = [60, 60, 60]);\nregion(\"Ni\") rotate([0, 0, 30]) ellipsoid(r = [15, 8, 8]);\n"),
        ("Notched beam (fracture specimen)", "difference() {\n    box(size = [80, 20, 30]);\n    // A V-notch from the top surface, along y.\n    translate([0, 0, 15]) rotate([90, 0, 0]) prism(n = 4, r = 8, h = 22);\n}\n"),
        ("Gear-like disc (polar array)", "difference() {\n    cylinder(r = 25, h = 10);\n    polar_array(n = 8) translate([25, 0, 0]) cylinder(r = 5, h = 12);\n}\n"),
        ("Embedded nanoparticle array", "box(size = [80, 80, 30]);\nregion(\"Au\") translate([-25, -25, 0]) grid_array(n = [3, 3, 1], step = [25, 25, 0]) sphere(r = 6);\n"),
    ]
}

pub fn commands() -> &'static [Command] {
    const fn c(group: &'static str, name: &'static str, snippet: &'static str, help: &'static str) -> Command {
        Command { group, name, snippet, help }
    }
    const LIST: &[Command] = &[
        c("Primitives", "sphere", "sphere(r = 10);", "Sphere of radius r (or diameter d) about the origin"),
        c("Primitives", "box", "box(size = [20, 20, 20]);", "Box of the given edge lengths, centred (center = false: corner at the origin)"),
        c("Primitives", "cylinder", "cylinder(r = 8, h = 30);", "Cylinder along z, centred; r1 and r2 make a frustum"),
        c("Primitives", "cone", "cone(r1 = 10, r2 = 0, h = 20);", "Cone or frustum along z: radius r1 at the bottom, r2 at the top"),
        c("Primitives", "torus", "torus(R = 20, r = 5);", "Torus about z: ring radius R, tube radius r"),
        c("Primitives", "capsule", "capsule(r = 6, h = 20);", "Cylinder of length h with hemispherical ends"),
        c("Primitives", "ellipsoid", "ellipsoid(r = [15, 10, 6]);", "Ellipsoid with semi-axes along x, y, z"),
        c("Primitives", "prism", "prism(n = 6, r = 12, h = 30);", "Regular n-sided prism along z (circumradius r, or apothem)"),
        c("Primitives", "polygon", "polygon(points = [[0, 0], [20, 0], [10, 15]], h = 10);", "Any polygon in the xy plane extruded along z"),
        c("Primitives", "facets", "facets(hkl = [[1, 1, 1], [1, 0, 0]], d = [15, 17]);", "Crystal facets {hkl} at distances d (Angstrom); family: all equivalent planes"),
        c("Primitives", "octahedron", "octahedron(d = 12);", "Regular octahedron, faces at distance d from the centre"),
        c("Primitives", "halfspace", "halfspace(normal = [0, 0, 1], d = 0);", "Everything below a plane (n . p <= d); also hkl = [h, k, l]"),
        c("Primitives", "slab", "slab(normal = [0, 0, 1], thickness = 10);", "A slab of a thickness about a plane (unbounded: intersect it)"),
        c("Booleans", "union", "union() {\n    \n}", "All of the shapes inside"),
        c("Booleans", "intersection", "intersection() {\n    \n}", "Where all of the shapes overlap"),
        c("Booleans", "difference", "difference() {\n    \n}", "The first shape minus the others"),
        c("Booleans", "smooth_union", "smooth_union(k = 3) {\n    \n}", "Union blended over about k Angstrom (necks, fillets)"),
        c("Booleans", "smooth_difference", "smooth_difference(k = 3) {\n    \n}", "Difference with rounded edges"),
        c("Transforms", "translate", "translate([0, 0, 10])", "Move the shapes that follow"),
        c("Transforms", "rotate", "rotate([0, 0, 45])", "Turn about x, then y, then z (degrees); rotate(a = 30, v = [1, 1, 1]) about an axis"),
        c("Transforms", "scale", "scale([1, 1, 2])", "Stretch along x, y, z"),
        c("Transforms", "mirror", "mirror([1, 0, 0])", "Reflect through the plane with this normal"),
        c("Modifiers", "shell", "shell(t = 3)", "Hollow: keep a skin of thickness t inside the surface"),
        c("Modifiers", "offset", "offset(r = 2)", "Grow (r > 0) or shrink (r < 0) the shape"),
        c("Patterns", "linear_array", "linear_array(n = 3, step = [15, 0, 0])", "n copies, each moved by step"),
        c("Patterns", "polar_array", "polar_array(n = 6, angle = 360, axis = [0, 0, 1])", "n copies turned about an axis through the origin"),
        c("Patterns", "grid_array", "grid_array(n = [3, 3, 1], step = [15, 15, 0])", "A 3D grid of copies"),
        c("Regions", "region", "region(\"Au\")", "Atoms of the shape inside get this element (\"none\": removed); adds no volume; later regions win"),
    ];
    LIST
}

/// Writes a mesh as Wavefront OBJ, triangle corners at the same place (to
/// 1e-6 Angstrom) shared as one vertex.
pub fn save_obj(mesh: &af_electronic::volume::Mesh, path: &str) -> Result<(), String> {
    use std::fmt::Write as _;
    let mut out = format!("# AtomForge shape: {} triangles (Angstrom)\n", mesh.vertices.len() / 3);
    let mut index: std::collections::HashMap<[i64; 3], usize> = std::collections::HashMap::new();
    let mut faces = Vec::with_capacity(mesh.vertices.len());
    let mut vertices: Vec<DVec3> = Vec::new();
    for v in &mesh.vertices {
        let key = [(v.x * 1e6).round() as i64, (v.y * 1e6).round() as i64, (v.z * 1e6).round() as i64];
        let id = *index.entry(key).or_insert_with(|| {
            vertices.push(*v);
            vertices.len() - 1
        });
        faces.push(id);
    }
    let g = |x: f64| af_core::format::format_g(x, 10);
    for v in &vertices {
        let _ = writeln!(out, "v {} {} {}", g(v.x), g(v.y), g(v.z));
    }
    for f in faces.chunks_exact(3) {
        let _ = writeln!(out, "f {} {} {}", f[0] + 1, f[1] + 1, f[2] + 1);
    }
    // The C++ writes through a text-mode stream: CRLF line ends on Windows.
    if cfg!(windows) {
        out = out.replace('\n', "\r\n");
    }
    std::fs::write(path, out).map_err(|_| format!("Cannot write {path}"))
}
