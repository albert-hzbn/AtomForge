//! Builders that make a structure from parameters (the C++
//! pipeline/GeneratorSteps): nanotubes, 2D layers, ribbons, twisted bilayers,
//! clusters, molecules and packings. Adsorbates, supercells, interpolation
//! and point defects are in generator_more.rs.

use std::collections::{BTreeSet, HashMap};

use af_core::cxx::{sort_by, stable_sort_by};
use af_core::cxx_random::{shuffle, uniform_real, Mt19937_64};
use af_core::format::format_g;
use af_core::json::Json;
use af_core::Structure;
use af_science::core::{add, dot, identity, inverse, norm, row_times, scale, sub, Mat3, Vec3};
use af_science::data::atomic_number;

use crate::molecule_library::{molecule, molecule_mass, principal_axes, MoleculeAtom};
use crate::pipeline::{params, ModifierType, PipelineData};
use crate::steps::basic::cxx_u64;
use crate::steps::edit::strtol;
use crate::steps::generator_more::{add_adsorbate, antisite, frenkel_pair, interpolate, near_cubic_supercell};
use crate::structure_edits::{append_atom, bond_neighbours, covalent_radius, rotate, saturate, SaturationOptions};
use crate::support::*;

pub(crate) type R = Result<(), String>;
pub(crate) const PI: f64 = 3.14159265358979323846;

pub(crate) fn g5(x: f64) -> String {
    format_g(x, 5)
}

fn g6(x: f64) -> String {
    format_g(x, 6)
}

// ---------------------------------------------------------------- helpers
pub(crate) fn words(text: &str) -> Vec<String> {
    text.split_ascii_whitespace().map(|w| w.to_string()).collect()
}

pub(crate) fn element_list(text: &str, count: usize, what: &str) -> Result<Vec<String>, String> {
    let mut out = words(text).iter().map(|w| standard_element(w)).collect::<Result<Vec<_>, _>>()?;
    if out.len() == 1 && count > 1 {
        out.resize(count, out[0].clone());
    }
    if out.len() != count {
        return Err(format!("{what}: give {count} element symbol{}", if count > 1 { "s" } else { "" }));
    }
    Ok(out)
}

fn produce(data: &mut PipelineData, s: Structure) {
    data.structure = s;
    data.clear_selection();
}

fn add_atom(s: &mut Structure, symbol: &str, position: Vec3) {
    let mut ignored = vec![0u8; s.atoms.len()];
    append_atom(s, &mut ignored, symbol, position, false);
}

fn set_cell(s: &mut Structure, a: Vec3, b: Vec3, c: Vec3, pbc: [bool; 3]) {
    s.has_unit_cell = true;
    s.cell_vectors = [a, b, c];
    s.cell_offset = [0.0; 3];
    s.pbc = pbc;
}

fn shift(s: &mut Structure, d: Vec3) {
    for a in s.atoms.iter_mut() {
        let q = add(position_of(a), d);
        place(a, q);
    }
}

fn rotate_z(v: Vec3, angle: f64) -> Vec3 {
    let (c, s) = (angle.cos(), angle.sin());
    [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]]
}

fn fmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

fn fmin(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

// ---------------------------------------------------------------- 2D layers
#[derive(Clone, Debug)]
struct Site {
    symbol: String,
    u: f64,
    v: f64,
    z: f64,
}

#[derive(Clone, Debug)]
struct Layer {
    a: f64,
    a1: Vec3,
    a2: Vec3,
    hexagonal: bool,
    sites: Vec<Site>,
    thickness: f64,
    interlayer: f64,
    kind: String,
}

impl Default for Layer {
    fn default() -> Self {
        Self { a: 0.0, a1: [0.0; 3], a2: [0.0; 3], hexagonal: true, sites: Vec::new(), thickness: 0.0, interlayer: 3.35, kind: String::new() }
    }
}

fn site(symbol: &str, u: f64, v: f64, z: f64) -> Site {
    Site { symbol: symbol.into(), u, v, z }
}

fn hexagonal_layer(a: f64) -> Layer {
    Layer { a, a1: [a, 0.0, 0.0], a2: [-a / 2.0, a * 3.0f64.sqrt() / 2.0, 0.0], ..Default::default() }
}

fn xene(element: &str, a: f64, buckling: f64) -> Layer {
    let mut l = hexagonal_layer(a);
    l.sites = vec![site(element, 1.0 / 3.0, 2.0 / 3.0, buckling / 2.0), site(element, 2.0 / 3.0, 1.0 / 3.0, -buckling / 2.0)];
    l.thickness = buckling;
    l
}

fn monolayer(kind: &str, elements_text: &str, a: f64, thickness: f64) -> Result<Layer, String> {
    let or = |default: &str| if elements_text.is_empty() { default.to_string() } else { elements_text.to_string() };
    let mut layer;
    match kind {
        "graphene" => {
            let e = element_list(&or("C"), 1, "graphene")?;
            layer = hexagonal_layer(if a > 0.0 { a } else { 2.4612 });
            layer.sites = vec![site(&e[0], 1.0 / 3.0, 2.0 / 3.0, 0.0), site(&e[0], 2.0 / 3.0, 1.0 / 3.0, 0.0)];
            layer.interlayer = 3.35;
        }
        "hbn" => {
            let e = element_list(&or("B N"), 2, "h-BN")?;
            layer = hexagonal_layer(if a > 0.0 { a } else { 2.504 });
            layer.sites = vec![site(&e[0], 1.0 / 3.0, 2.0 / 3.0, 0.0), site(&e[1], 2.0 / 3.0, 1.0 / 3.0, 0.0)];
            layer.interlayer = 3.33;
        }
        "silicene" | "germanene" | "stanene" => {
            let (element, a0, b0) = match kind {
                "silicene" => ("Si", 3.866, 0.44),
                "germanene" => ("Ge", 4.02, 0.69),
                _ => ("Sn", 4.67, 0.85),
            };
            let e = element_list(&or(element), 1, kind)?;
            layer = xene(&e[0], if a > 0.0 { a } else { a0 }, if thickness > 0.0 { thickness } else { b0 });
            layer.interlayer = layer.thickness + 3.0;
        }
        "mx2-2h" | "mx2-1t" => {
            let e = element_list(&or("Mo S"), 2, "MX2")?;
            layer = hexagonal_layer(if a > 0.0 { a } else { 3.18 });
            layer.thickness = if thickness > 0.0 { thickness } else { 3.19 };
            let h = layer.thickness / 2.0;
            layer.sites = vec![site(&e[0], 0.0, 0.0, 0.0), site(&e[1], 1.0 / 3.0, 2.0 / 3.0, h)];
            if kind == "mx2-2h" {
                layer.sites.push(site(&e[1], 1.0 / 3.0, 2.0 / 3.0, -h));
            } else {
                layer.sites.push(site(&e[1], 2.0 / 3.0, 1.0 / 3.0, -h));
            }
            layer.interlayer = 6.15;
        }
        "phosphorene" => {
            let e = element_list(&or("P"), 1, "phosphorene")?;
            let za = if a > 0.0 { a } else { 3.3136 };
            let zb = 4.3763 * if a > 0.0 { a / 3.3136 } else { 1.0 };
            layer = Layer { a: za, a1: [za, 0.0, 0.0], a2: [0.0, zb, 0.0], hexagonal: false, ..Default::default() };
            let h = if thickness > 0.0 { thickness / 2.0 } else { 0.10168 * 10.478 };
            layer.sites = vec![site(&e[0], 0.0, 0.08056, h), site(&e[0], 0.0, 0.91944, -h), site(&e[0], 0.5, 0.58056, -h), site(&e[0], 0.5, 0.41944, h)];
            layer.thickness = 2.0 * h;
            layer.interlayer = 5.239;
        }
        _ => return Err(format!("Unknown layer '{kind}': graphene, hbn, silicene, germanene, stanene, mx2-2h, mx2-1t or phosphorene")),
    }
    layer.kind = kind.to_string();
    Ok(layer)
}

fn strain_layer(l: &mut Layer, strain: f64) {
    l.a *= 1.0 + strain;
    l.a1 = scale(l.a1, 1.0 + strain);
    l.a2 = scale(l.a2, 1.0 + strain);
}

#[derive(Clone, Copy, Debug, Default)]
struct Registry {
    du: f64,
    dv: f64,
    inverted: bool,
}

fn registry_of(stacking: &str, k: i32, layer: &Layer) -> Result<Registry, String> {
    if stacking == "AA" || k == 0 {
        return Ok(Registry::default());
    }
    if !layer.hexagonal {
        return Err(format!("Only AA stacking is available for {}", layer.kind));
    }
    match stacking {
        "AB" => Ok(if k % 2 != 0 { Registry { du: 2.0 / 3.0, dv: 1.0 / 3.0, inverted: false } } else { Registry::default() }),
        "ABC" => Ok(Registry { du: (k as f64 * 2.0 / 3.0) % 1.0, dv: (k as f64 * 1.0 / 3.0) % 1.0, inverted: false }),
        "2H" => Ok(if k % 2 != 0 { Registry { du: 1.0 / 3.0, dv: 2.0 / 3.0, inverted: true } } else { Registry::default() }),
        _ => Err("stacking must be AA, AB, ABC or 2H".into()),
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_layer(s: &mut Structure, layer: &Layer, angle: f64, l1: Vec3, l2: Vec3, z: f64, count: &mut usize, registry: Registry) -> R {
    let a1 = rotate_z(layer.a1, angle);
    let a2 = rotate_z(layer.a2, angle);
    let inv = inverse(&[l1, l2, [0.0, 0.0, 1.0]])?;
    let shortest = fmin(norm(layer.a1), norm(layer.a2));
    let reach = (norm(l1) + norm(l2)) / shortest + 3.0;
    let range = reach.ceil() as i32 + 2;
    let eps = 1e-7;
    for i in -2 * range..=2 * range {
        for j in -2 * range..=2 * range {
            for st in &layer.sites {
                let u = (if registry.inverted { -st.u } else { st.u }) + registry.du;
                let v = (if registry.inverted { -st.v } else { st.v }) + registry.dv;
                let q = add(scale(a1, i as f64 + u), scale(a2, j as f64 + v));
                let f = row_times([q[0], q[1], 0.0], &inv);
                if f[0] < -eps || f[0] >= 1.0 - eps || f[1] < -eps || f[1] >= 1.0 - eps {
                    continue;
                }
                add_atom(s, &st.symbol, [q[0], q[1], z + if registry.inverted { -st.z } else { st.z }]);
                *count += 1;
            }
        }
    }
    Ok(())
}

fn build_monolayer(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-monolayer";
    let thickness = number(p, id, "thickness")?;
    let a = number(p, id, "a")?;
    let elements_text = text(p, id, "elements")?;
    let mut layer = monolayer(&text(p, id, "kind")?, &elements_text, a, thickness)?;
    strain_layer(&mut layer, number(p, id, "strain")? / 100.0);
    let na = integer_of(p, id, "na", 1)?;
    let nb = integer_of(p, id, "nb", 1)?;
    let layers = integer_of(p, id, "layers", 1)?;
    let vacuum = number(p, id, "vacuum")?;
    if vacuum < 0.0 {
        return Err("vacuum must be nonnegative".into());
    }
    let interlayer = if number(p, id, "interlayer")? > 0.0 { number(p, id, "interlayer")? } else { layer.interlayer };
    let stacking = text(p, id, "stacking")?;
    let (mut l1, mut l2) = (layer.a1, layer.a2);
    let rectangular = flag(p, id, "rectangular")? && layer.hexagonal;
    if rectangular {
        l2 = add(layer.a1, scale(layer.a2, 2.0));
    }
    l1 = scale(l1, na as f64);
    l2 = scale(l2, nb as f64);
    let stack = (layers - 1) as f64 * interlayer;
    let height = stack + layer.thickness + vacuum;
    let mut s = Structure::default();
    let mut count = 0;
    for k in 0..layers {
        let registry = registry_of(&stacking, k, &layer)?;
        fill_layer(&mut s, &layer, 0.0, l1, l2, (height - stack) / 2.0 + k as f64 * interlayer, &mut count, registry)?;
    }
    set_cell(&mut s, l1, l2, [0.0, 0.0, fmax(height, 1.0)], [true, true, flag(p, id, "periodic_c")?]);
    produce(data, s);
    data.notes.push(format!(
        "{} atoms, {} layer{}, a = {} Angstrom{}",
        data.structure.atoms.len(),
        layers,
        if layers > 1 { format!("s {stacking}") } else { String::new() },
        g6(layer.a),
        if rectangular { ", rectangular cell" } else { "" }
    ));
    Ok(())
}

// ---------------------------------------------------------------- nanotubes
fn gcd(a: i32, b: i32) -> i32 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

#[allow(clippy::too_many_arguments)]
fn add_tube_wall(s: &mut Structure, n: i32, m: i32, bond: f64, species: &[String], cells: i32, radius: &mut f64, period: &mut f64) -> R {
    let a = bond * 3.0f64.sqrt();
    let a1 = [a * 3.0f64.sqrt() / 2.0, a / 2.0, 0.0];
    let a2 = [a * 3.0f64.sqrt() / 2.0, -a / 2.0, 0.0];
    let c = add(scale(a1, n as f64), scale(a2, m as f64));
    let d_r = gcd(2 * m + n, 2 * n + m);
    let t1 = (2 * m + n) / d_r;
    let t2 = -(2 * n + m) / d_r;
    let t = add(scale(a1, t1 as f64), scale(a2, t2 as f64));
    let (lc, lt) = (norm(c), norm(t));
    *radius = lc / (2.0 * PI);
    *period = lt;
    let basis = [[0.0; 3], scale(add(a1, a2), 1.0 / 3.0)];
    let range = n.abs() + m.abs() + t1.abs() + t2.abs() + 2;
    let eps = 1e-9;
    let mut found = 0usize;
    for i in -range..=range {
        for j in -range..=range {
            for b in 0..2 {
                let q = add(add(scale(a1, i as f64), scale(a2, j as f64)), basis[b]);
                let u = dot(q, c) / (lc * lc);
                let v = dot(q, t) / (lt * lt);
                if u < -eps || u >= 1.0 - eps || v < -eps || v >= 1.0 - eps {
                    continue;
                }
                found += 1;
                for k in 0..cells {
                    let theta = 2.0 * PI * u;
                    add_atom(s, &species[b], [*radius * theta.cos(), *radius * theta.sin(), (v + k as f64) * lt]);
                }
            }
        }
    }
    let expected = (4 * (n * n + n * m + m * m) / d_r) as usize;
    if found != expected {
        return Err(format!("Nanotube unit cell has {found} atoms, expected {expected}"));
    }
    Ok(())
}

/// Removes atoms with fewer than two bonds, repeatedly.
fn remove_dangling(s: &mut Structure) -> Result<usize, String> {
    let mut removed = 0;
    for _ in 0..20 {
        let nb = bond_neighbours(s, 1.2)?;
        let kept: Vec<_> = s.atoms.iter().enumerate().filter(|(i, _)| nb.atoms[*i].len() >= 2).map(|(_, a)| a.clone()).collect();
        if kept.len() == s.atoms.len() {
            break;
        }
        removed += s.atoms.len() - kept.len();
        s.atoms = kept;
    }
    Ok(removed)
}

fn axis_to(s: &mut Structure, axis: &str) {
    if axis == "z" {
        return;
    }
    let permute = |v: [f64; 3]| if axis == "x" { [v[2], v[0], v[1]] } else { [v[1], v[2], v[0]] };
    for a in s.atoms.iter_mut() {
        let q = permute([a.x, a.y, a.z]);
        place(a, q);
    }
    if s.has_unit_cell {
        let (old, old_pbc) = (s.cell_vectors, s.pbc);
        for k in 0..3 {
            let from = if axis == "x" { (k + 2) % 3 } else { (k + 1) % 3 };
            s.cell_vectors[k] = permute(old[from]);
            s.pbc[k] = old_pbc[from];
        }
    }
}

fn build_nanotube(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-nanotube";
    let n = integer_of(p, id, "n", 0)?;
    let m = integer_of(p, id, "m", 0)?;
    if n == 0 && m == 0 {
        return Err("The chiral indices (n, m) must not both be zero".into());
    }
    if m > n {
        return Err("Use n >= m (the tube (m, n) is the mirror image of (n, m))".into());
    }
    let mut cells = integer_of(p, id, "length", 1)?;
    let walls = integer_of(p, id, "walls", 1)?;
    let species = element_list(&text(p, id, "elements")?, 2, "nanotube")?;
    let mut bond = number(p, id, "bond")?;
    if !(bond > 0.0) {
        bond = if species[0] == "C" && species[1] == "C" {
            1.421
        } else if (species[0] == "B" && species[1] == "N") || (species[0] == "N" && species[1] == "B") {
            1.446
        } else {
            covalent_radius(atomic_number(&species[0])) + covalent_radius(atomic_number(&species[1]))
        };
    }
    let spacing = number(p, id, "wall_spacing")?;
    let vacuum = number(p, id, "vacuum")?;
    if !(spacing > 0.0) || vacuum < 0.0 {
        return Err("wall spacing must be positive, vacuum nonnegative".into());
    }
    if walls > 1 && !(n == m || m == 0) {
        return Err("Multi-wall tubes need armchair (n, n) or zigzag (n, 0) indices, which share the axial period".into());
    }
    {
        let (mut radius, mut period) = (0.0, 0.0);
        let mut probe = Structure::default();
        add_tube_wall(&mut probe, n, m, bond, &species, 1, &mut radius, &mut period)?;
        let wanted = number(p, id, "length_A")?;
        if wanted > 0.0 {
            cells = 1.max((wanted / period).round() as i32);
        }
    }
    let mut s = Structure::default();
    let (mut radius, mut period, mut outer) = (0.0, 0.0, 0.0);
    let mut notes = String::new();
    for w in 0..walls {
        let a = bond * 3.0f64.sqrt();
        let step = if n == m { (spacing * 2.0 * PI / (a * 3.0f64.sqrt())).round() as i32 } else { (spacing * 2.0 * PI / a).round() as i32 };
        let nw = n + w * step;
        let mw = if m == 0 { 0 } else { m + w * step };
        add_tube_wall(&mut s, nw, mw, bond, &species, cells, &mut radius, &mut period)?;
        outer = radius;
        notes += &format!("{}({},{}) d = {}", if w > 0 { ", " } else { "" }, nw, mw, g5(2.0 * radius));
    }
    let finite = flag(p, id, "finite")?;
    let bx = 2.0 * outer + vacuum + 2.0 * bond;
    shift(&mut s, [bx / 2.0, bx / 2.0, if finite { vacuum / 2.0 } else { 0.0 }]);
    let (mut capped, mut dropped) = (0usize, 0usize);
    if finite {
        set_cell(&mut s, [bx, 0.0, 0.0], [0.0, bx, 0.0], [0.0, 0.0, period * cells as f64 + vacuum], [false; 3]);
        dropped = remove_dangling(&mut s)?;
        let cap = text(p, id, "end_cap")?;
        if cap != "none" {
            let mut selected = vec![0u8; s.atoms.len()];
            let options = SaturationOptions { terminator: standard_element(&cap)?, ..Default::default() };
            let set: BTreeSet<String> = [species[0].clone(), species[1].clone()].into_iter().collect();
            capped = saturate(&mut s, &mut selected, &set, false, 1.2, false, &options)?;
        }
    } else {
        set_cell(&mut s, [bx, 0.0, 0.0], [0.0, bx, 0.0], [0.0, 0.0, period * cells as f64], [false, false, true]);
    }
    axis_to(&mut s, &text(p, id, "axis")?);
    produce(data, s);
    let chiral = (3.0f64.sqrt() * m as f64).atan2(2.0 * n as f64 + m as f64) * 180.0 / PI;
    let mut note = format!(
        "{} atoms; walls {} Angstrom; chiral angle {} degrees; period {} Angstrom x {}",
        data.structure.atoms.len(),
        notes,
        g5(chiral),
        g5(period),
        cells
    );
    if species[0] == "C" && species[1] == "C" {
        note += &format!("; {} (zone folding)", if (n - m) % 3 == 0 { "metallic" } else { "semiconducting" });
    }
    if finite {
        note += &format!(
            "; finite, {} end atoms {}{}",
            capped,
            text(p, id, "end_cap")?,
            if dropped > 0 { format!(", {dropped} dangling atoms removed") } else { String::new() }
        );
    }
    data.notes.push(note);
    Ok(())
}

// ---------------------------------------------------------------- nanoribbons
fn build_nanoribbon(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-nanoribbon";
    let edge = text(p, id, "edge")?;
    let width = integer_of(p, id, "width", 1)?;
    let length = integer_of(p, id, "length", 1)?;
    let d = number(p, id, "bond")?;
    let vacuum = number(p, id, "vacuum")?;
    if !(d > 0.0) || vacuum < 0.0 {
        return Err("bond must be positive, vacuum nonnegative".into());
    }
    let species = element_list(&text(p, id, "elements")?, 2, "nanoribbon")?;
    let (lx, ly) = (3.0f64.sqrt() * d, 3.0 * d);
    let cell_sites: [Vec3; 4] = [[0.0, 0.0, 0.0], [lx / 2.0, d / 2.0, 0.0], [lx / 2.0, 1.5 * d, 0.0], [0.0, 2.0 * d, 0.0]];
    let sublattice = [0usize, 1, 0, 1];
    let zigzag = edge == "zigzag";
    if !zigzag && edge != "armchair" {
        return Err("edge must be armchair or zigzag".into());
    }
    // (along, across, sublattice)
    let mut sites: Vec<(f64, f64, usize)> = Vec::new();
    for r in 0..width + 2 {
        for k in 0..4 {
            let q = add(cell_sites[k], if zigzag { [0.0, r as f64 * ly, 0.0] } else { [r as f64 * lx, 0.0, 0.0] });
            sites.push(if zigzag { (q[0], q[1], sublattice[k]) } else { (q[1], q[0], sublattice[k]) });
        }
    }
    let mut levels: Vec<f64> = Vec::new();
    for st in &sites {
        if !levels.iter().any(|l| (l - st.1).abs() < 1e-6) {
            levels.push(st.1);
        }
    }
    sort_by(&mut levels, |a, b| a < b);
    let rows = if zigzag { 2 * width as usize } else { width as usize };
    if levels.len() < rows {
        return Err("ribbon too wide".into());
    }
    let top = levels[rows - 1];
    let period = if zigzag { lx } else { ly };
    let ribbon_width = top - levels[0];
    let mut s = Structure::default();
    for k in 0..length {
        for st in &sites {
            if st.1 <= top + 1e-6 {
                add_atom(&mut s, &species[st.2], [st.0 + k as f64 * period, st.1, 0.0]);
            }
        }
    }
    let across = ribbon_width + vacuum + 2.2;
    let height = fmax(vacuum, 1.0);
    let finite = flag(p, id, "finite")?;
    shift(&mut s, [if finite { vacuum / 2.0 } else { 0.0 }, (across - ribbon_width) / 2.0, height / 2.0]);
    if finite {
        set_cell(&mut s, [period * length as f64 + vacuum, 0.0, 0.0], [0.0, across, 0.0], [0.0, 0.0, height], [false; 3]);
    } else {
        set_cell(&mut s, [period * length as f64, 0.0, 0.0], [0.0, across, 0.0], [0.0, 0.0, height], [true, false, false]);
    }
    let dropped = if finite { remove_dangling(&mut s)? } else { 0 };
    let mut selected = vec![0u8; s.atoms.len()];
    let mut terminated = 0;
    let termination = text(p, id, "termination")?;
    if flag(p, id, "passivate")? && termination != "none" {
        let options = SaturationOptions { terminator: standard_element(&termination)?, bond_length: number(p, id, "termination_bond")?, ..Default::default() };
        let set: BTreeSet<String> = [species[0].clone(), species[1].clone()].into_iter().collect();
        terminated = saturate(&mut s, &mut selected, &set, false, 1.2, false, &options)?;
    }
    produce(data, s);
    let mut note = format!(
        "{}{} edges, width {} Angstrom; {} atoms ({} {})",
        width,
        if zigzag { " chains, zigzag" } else { " dimer lines, armchair" },
        g5(ribbon_width),
        data.structure.atoms.len(),
        terminated,
        termination
    );
    if finite {
        note += &format!("; finite flake{}", if dropped > 0 { format!(", {dropped} dangling atoms removed") } else { String::new() });
    }
    if !zigzag {
        note += &format!("; family {} (3p, 3p+1, 3p+2 differ in band gap)", width % 3);
    }
    data.notes.push(note);
    Ok(())
}

// ---------------------------------------------------------------- twisted bilayers
fn build_twisted_bilayer(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-twisted-bilayer";
    let thickness = number(p, id, "thickness")?;
    let a = number(p, id, "a")?;
    let elements_text = text(p, id, "elements")?;
    let layer = monolayer(&text(p, id, "kind")?, &elements_text, a, thickness)?;
    if !layer.hexagonal {
        return Err("Twisted bilayers need hexagonal layers".into());
    }
    let top_kind = text(p, id, "top_kind")?;
    let mut top = if top_kind == "same" { layer.clone() } else { monolayer(&top_kind, &text(p, id, "top_elements")?, 0.0, 0.0)? };
    if !top.hexagonal {
        return Err("Twisted bilayers need hexagonal layers".into());
    }
    let natural_top = top.a;
    let mismatch = layer.a / natural_top - 1.0;
    strain_layer(&mut top, mismatch);
    let angle_of = |i: i32| {
        let i = i as f64;
        let c = (3.0 * i * i + 3.0 * i + 0.5) / (3.0 * i * i + 3.0 * i + 1.0);
        fmin(1.0, c).acos()
    };
    let mut index = integer_of(p, id, "index", 0)?;
    let wanted = number(p, id, "angle")?;
    if wanted > 0.0 {
        index = 1;
        for i in 1..=300 {
            if (angle_of(i) * 180.0 / PI - wanted).abs() < (angle_of(index) * 180.0 / PI - wanted).abs() {
                index = i;
            }
        }
    }
    if index < 1 {
        return Err("Give a twist angle or a commensurate index >= 1".into());
    }
    if index > 300 {
        return Err("index at most 300 (the cell grows as 3 i^2)".into());
    }
    let theta = angle_of(index);
    let b1 = layer.a1;
    let b2 = add(layer.a1, layer.a2);
    let l1 = add(scale(b1, index as f64), scale(b2, (index + 1) as f64));
    let l2 = rotate_z(l1, PI / 3.0);
    let interlayer = if number(p, id, "interlayer")? > 0.0 { number(p, id, "interlayer")? } else { (layer.interlayer + top.interlayer) / 2.0 };
    let vacuum = number(p, id, "vacuum")?;
    let height = interlayer + (layer.thickness + top.thickness) / 2.0 + vacuum;
    let mut s = Structure::default();
    let (mut bottom_count, mut top_count) = (0usize, 0usize);
    let z0 = (height - interlayer) / 2.0;
    fill_layer(&mut s, &layer, 0.0, l1, l2, z0, &mut bottom_count, Registry::default())?;
    let v2 = add(scale(b1, (index + 1) as f64), scale(b2, index as f64));
    let turn = l1[1].atan2(l1[0]) - v2[1].atan2(v2[0]);
    let stacking = text(p, id, "stacking")?;
    let registry = if stacking == "AB" { Registry { du: 2.0 / 3.0, dv: 1.0 / 3.0, inverted: false } } else { Registry::default() };
    if stacking != "AA" && stacking != "AB" {
        return Err("stacking must be AA or AB".into());
    }
    fill_layer(&mut s, &top, turn, l1, l2, z0 + interlayer, &mut top_count, registry)?;
    let cells = (3 * index * index + 3 * index + 1) as usize;
    let (expected_bottom, expected_top) = (layer.sites.len() * cells, top.sites.len() * cells);
    if bottom_count != expected_bottom || top_count != expected_top {
        return Err(format!("Twisted bilayer layers have {bottom_count} and {top_count} atoms, expected {expected_bottom} and {expected_top}"));
    }
    set_cell(&mut s, l1, l2, [0.0, 0.0, height], [true; 3]);
    produce(data, s);
    let mut note = format!(
        "Twist {} degrees (index {}, {} before the twist), {} atoms, |L| = {} Angstrom, moire period {} Angstrom",
        g6(theta * 180.0 / PI),
        index,
        stacking,
        data.structure.atoms.len(),
        g6(norm(l1)),
        g6(layer.a / (2.0 * (theta / 2.0).sin()))
    );
    if top_kind != "same" {
        note += &format!("; top {} strained by {}% from a = {}", top_kind, g6(100.0 * mismatch), g6(natural_top));
    }
    data.notes.push(note);
    Ok(())
}

// ---------------------------------------------------------------- clusters
fn build_cluster(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-cluster";
    let kind = text(p, id, "kind")?;
    let symbol = standard_element(&text(p, id, "element")?)?;
    let shells = integer_of(p, id, "shells", 1)?;
    let d = number(p, id, "distance")?;
    let vacuum = number(p, id, "vacuum")?;
    if !(d > 0.0) || vacuum < 0.0 {
        return Err("distance must be positive, vacuum nonnegative".into());
    }
    let mut points: Vec<Vec3> = vec![[0.0; 3]];
    let mut seen: BTreeSet<[i64; 3]> = BTreeSet::from([[0, 0, 0]]);
    let add_unique = |points: &mut Vec<Vec3>, seen: &mut BTreeSet<[i64; 3]>, q: Vec3| {
        let key = [(q[0] / (1e-4 * d)).round() as i64, (q[1] / (1e-4 * d)).round() as i64, (q[2] / (1e-4 * d)).round() as i64];
        if seen.insert(key) {
            points.push(q);
        }
    };
    match kind.as_str() {
        "icosahedron" => {
            let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
            let mut v: Vec<Vec3> = Vec::new();
            for s1 in [-1.0, 1.0] {
                for s2 in [-1.0, 1.0] {
                    v.push([0.0, s1, phi * s2]);
                    v.push([s1, phi * s2, 0.0]);
                    v.push([phi * s2, 0.0, s1]);
                }
            }
            let edge = |x: usize, y: usize| (norm(sub(v[x], v[y])) - 2.0).abs() < 1e-6;
            let mut faces: Vec<[usize; 3]> = Vec::new();
            for i in 0..12 {
                for j in i + 1..12 {
                    for k in j + 1..12 {
                        if edge(i, j) && edge(j, k) && edge(i, k) {
                            faces.push([i, j, k]);
                        }
                    }
                }
            }
            for shell in 1..shells {
                for f in &faces {
                    let to_vertex = shell as f64 * d / norm(v[0]);
                    let (a, b, c) = (scale(v[f[0]], to_vertex), scale(v[f[1]], to_vertex), scale(v[f[2]], to_vertex));
                    for i in 0..=shell {
                        for j in 0..=shell - i {
                            add_unique(
                                &mut points,
                                &mut seen,
                                add(a, add(scale(sub(b, a), i as f64 / shell as f64), scale(sub(c, a), j as f64 / shell as f64))),
                            );
                        }
                    }
                }
            }
        }
        "decahedron" => {
            if shells < 2 {
                return Err("A decahedron needs at least 2 atoms along an edge (shells >= 2)".into());
            }
            let n = shells - 1;
            let edge = n as f64 * d;
            let r = edge / (2.0 * (PI / 5.0).sin());
            let h = fmax(0.0, edge * edge - r * r).sqrt();
            let (top, bottom) = ([0.0, 0.0, h], [0.0, 0.0, -h]);
            points.clear();
            seen.clear();
            for k in 0..5 {
                let ang1 = 2.0 * PI * k as f64 / 5.0;
                let ang2 = 2.0 * PI * (k + 1) as f64 / 5.0;
                let p1 = [r * ang1.cos(), r * ang1.sin(), 0.0];
                let p2 = [r * ang2.cos(), r * ang2.sin(), 0.0];
                for i in 0..=n {
                    for j in 0..=n - i {
                        for l in 0..=n - i - j {
                            let m = n - i - j - l;
                            let q = scale(add(add(scale(top, i as f64), scale(bottom, j as f64)), add(scale(p1, l as f64), scale(p2, m as f64))), 1.0 / n as f64);
                            add_unique(&mut points, &mut seen, q);
                        }
                    }
                }
            }
        }
        "cuboctahedron" | "octahedron" => {
            let unit = d / 2.0f64.sqrt();
            let k = shells - 1;
            points.clear();
            for x in -k..=k {
                for y in -k..=k {
                    for z in -k..=k {
                        let keep = if kind == "cuboctahedron" {
                            ((x + y + z) % 2 + 2) % 2 == 0 && x.abs() + y.abs() + z.abs() <= 2 * k
                        } else {
                            ((x + y + z - k) % 2 + 2) % 2 == 0 && x.abs() + y.abs() + z.abs() <= k
                        };
                        if keep {
                            points.push([x as f64 * unit, y as f64 * unit, z as f64 * unit]);
                        }
                    }
                }
            }
        }
        _ => return Err(format!("Unknown cluster '{kind}': icosahedron, decahedron, cuboctahedron or octahedron")),
    }
    let mut extent = 0.0;
    for q in &points {
        extent = fmax(extent, norm(*q));
    }
    let bx = 2.0 * extent + vacuum;
    let mut symbols = vec![symbol.clone(); points.len()];
    let second = text(p, id, "second_element")?;
    let mut substituted = 0usize;
    if !second.is_empty() {
        let other = standard_element(&second)?;
        let fraction = number(p, id, "fraction")?;
        if !(0.0..=1.0).contains(&fraction) {
            return Err("fraction must lie between 0 and 1".into());
        }
        let distribution = text(p, id, "distribution")?;
        let mut order: Vec<usize> = (0..points.len()).collect();
        if distribution == "random" {
            let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
            shuffle(&mut order, &mut rng);
        } else {
            let key: Vec<f64> = points
                .iter()
                .map(|q| match distribution.as_str() {
                    "core" => norm(*q),
                    "shell" => -norm(*q),
                    "janus" => -q[2],
                    _ => 0.0,
                })
                .collect();
            if distribution != "core" && distribution != "shell" && distribution != "janus" {
                return Err("distribution must be random, core, shell or janus".into());
            }
            stable_sort_by(&mut order, |&a, &b| key[a] < key[b] - 1e-9);
        }
        substituted = (fraction * points.len() as f64).round() as usize;
        for k in 0..substituted {
            symbols[order[k]] = other.clone();
        }
    }
    let mut s = Structure::default();
    for (i, q) in points.iter().enumerate() {
        add_atom(&mut s, &symbols[i], add(*q, [bx / 2.0, bx / 2.0, bx / 2.0]));
    }
    set_cell(&mut s, [bx, 0.0, 0.0], [0.0, bx, 0.0], [0.0, 0.0, bx], [false; 3]);
    produce(data, s);
    let mut note = format!("{} of {} shells: {} atoms, diameter {} Angstrom", kind, shells, points.len(), g5(2.0 * extent));
    if substituted > 0 {
        note += &format!("; {} {} ({})", substituted, standard_element(&second)?, text(p, id, "distribution")?);
    }
    data.notes.push(note);
    Ok(())
}

// ---------------------------------------------------------------- molecules
pub(crate) fn rotated(mut atoms: Vec<MoleculeAtom>, r: &Mat3) -> Vec<MoleculeAtom> {
    for a in atoms.iter_mut() {
        a.position = rotate(r, a.position);
    }
    atoms
}

/// Rows: the principal axes onto z, x, y (long molecules along z).
pub(crate) fn principal_frame(atoms: &[MoleculeAtom]) -> Mat3 {
    let axes = principal_axes(atoms);
    [axes[1], axes[2], axes[0]]
}

/// A uniformly random rotation from a random unit quaternion.
fn random_rotation(rng: &mut Mt19937_64) -> Mat3 {
    let u1 = uniform_real(rng, 0.0, 1.0);
    let u2 = uniform_real(rng, 0.0, 1.0) * 2.0 * PI;
    let u3 = uniform_real(rng, 0.0, 1.0) * 2.0 * PI;
    let w = (1.0 - u1).sqrt() * u2.sin();
    let x = (1.0 - u1).sqrt() * u2.cos();
    let y = u1.sqrt() * u3.sin();
    let z = u1.sqrt() * u3.cos();
    [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
    ]
}

fn build_molecule(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-molecule";
    let mut atoms = molecule(&text(p, id, "name")?)?;
    let orientation = text(p, id, "orientation")?;
    if orientation == "principal" && atoms.len() > 1 {
        let frame = principal_frame(&atoms);
        atoms = rotated(atoms, &frame);
    } else if orientation == "random" {
        let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
        let r = random_rotation(&mut rng);
        atoms = rotated(atoms, &r);
    } else if orientation != "library" && orientation != "principal" {
        return Err("orientation must be library, principal or random".into());
    }
    let vacuum = number(p, id, "vacuum")?;
    let mut s = Structure::default();
    let mut extent = 0.0;
    for a in &atoms {
        extent = fmax(extent, norm(a.position));
    }
    let bx = 2.0 * extent + vacuum;
    for a in &atoms {
        add_atom(&mut s, &a.symbol, if vacuum > 0.0 { add(a.position, [bx / 2.0, bx / 2.0, bx / 2.0]) } else { a.position });
    }
    if vacuum > 0.0 {
        let periodic = flag(p, id, "periodic")?;
        set_cell(&mut s, [bx, 0.0, 0.0], [0.0, bx, 0.0], [0.0, 0.0, bx], [periodic; 3]);
    }
    produce(data, s);
    data.notes.push(format!("{} atoms, {} u, extent {} Angstrom", atoms.len(), g6(molecule_mass(&atoms)), g6(2.0 * extent)));
    Ok(())
}

// ---------------------------------------------------------------- packings
struct Region {
    kind: String,
    low: Vec3,
    high: Vec3,
    centre: Vec3,
    radius: f64,
    axis: usize,
}

impl Region {
    fn contains(&self, q: Vec3) -> bool {
        match self.kind.as_str() {
            "all" => true,
            "box" => q[0] >= self.low[0] && q[0] <= self.high[0] && q[1] >= self.low[1] && q[1] <= self.high[1] && q[2] >= self.low[2] && q[2] <= self.high[2],
            "sphere" => norm(sub(q, self.centre)) <= self.radius,
            "outside-sphere" => norm(sub(q, self.centre)) >= self.radius,
            "slab" => q[self.axis] >= self.low[self.axis] && q[self.axis] <= self.high[self.axis],
            "cylinder" => {
                let mut d = sub(q, self.centre);
                d[self.axis] = 0.0;
                norm(d) <= self.radius
            }
            _ => true,
        }
    }

    fn volume(&self, bx: Vec3) -> f64 {
        match self.kind.as_str() {
            "box" => fmax(0.0, self.high[0] - self.low[0]) * fmax(0.0, self.high[1] - self.low[1]) * fmax(0.0, self.high[2] - self.low[2]),
            "sphere" => 4.0 / 3.0 * PI * self.radius * self.radius * self.radius,
            "outside-sphere" => bx[0] * bx[1] * bx[2] - 4.0 / 3.0 * PI * self.radius * self.radius * self.radius,
            "slab" => (self.high[self.axis] - self.low[self.axis]) * bx[0] * bx[1] * bx[2] / bx[self.axis],
            "cylinder" => PI * self.radius * self.radius * bx[self.axis],
            _ => bx[0] * bx[1] * bx[2],
        }
    }
}

#[derive(Clone)]
struct Kind {
    atoms: Vec<MoleculeAtom>,
    count: i32,
    name: String,
    automatic: bool,
}

fn build_packing(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-packing";
    let list = words(&text(p, id, "molecules")?);
    if list.is_empty() || list.len() % 2 != 0 {
        return Err("molecules: give pairs of a molecule and a count, e.g. H2O 100".into());
    }
    let mut kinds: Vec<Kind> = Vec::new();
    for k in (0..list.len()).step_by(2) {
        let (count, whole) = strtol(&list[k + 1]);
        if !whole || count < 0 {
            return Err(format!("molecules: '{}' is not a count", list[k + 1]));
        }
        kinds.push(Kind { atoms: molecule(&list[k])?, count: count as i32, name: list[k].clone(), automatic: count == 0 });
    }
    let solvate = text(p, id, "mode")? == "solvate";
    let density = number(p, id, "density")?;
    let mut bx = vector(p, id, "box")?;
    let mut s = Structure::default();
    if solvate {
        require_cell(data, "Solvation")?;
        s = data.structure.clone();
        let cell = cell_of(&s);
        for r in 0..3 {
            for c in 0..3 {
                if r != c && cell[r][c].abs() > 1e-6 {
                    return Err("Solvation needs an orthogonal cell (use shrink-wrap or set-cell)".into());
                }
            }
        }
        bx = [cell[0][0], cell[1][1], cell[2][2]];
    }
    let axis_name = {
        let region_kind = text(p, id, "region")?;
        let region = Region {
            kind: region_kind,
            low: vector(p, id, "region_low")?,
            high: vector(p, id, "region_high")?,
            centre: vector(p, id, "region_center")?,
            radius: number(p, id, "region_radius")?,
            axis: 2,
        };
        (region, text(p, id, "region_axis")?)
    };
    let (mut region, axis_name) = axis_name;
    region.axis = match axis_name.as_str() {
        "x" => 0,
        "y" => 1,
        _ => 2,
    };
    let mut fixed_mass = 0.0;
    let mut automatic = false;
    for k in &kinds {
        fixed_mass += molecule_mass(&k.atoms) * k.count as f64;
        automatic = automatic || k.count == 0;
    }
    if bx[0] <= 0.0 || bx[1] <= 0.0 || bx[2] <= 0.0 {
        if !(density > 0.0) || automatic {
            return Err("Give the box lengths, or a density with every count".into());
        }
        let side = (fixed_mass * 1.66053906660 / density).cbrt();
        bx = [side; 3];
    }
    if automatic {
        if !(density > 0.0) {
            return Err("A count of 0 needs a density to fill to".into());
        }
        let n: [i32; 3] = std::array::from_fn(|k| ((bx[k] / 0.4).ceil() as i32).clamp(1, 400));
        let mut taken = vec![false; n[0] as usize * n[1] as usize * n[2] as usize];
        for a in &s.atoms {
            let r = covalent_radius(a.atomic_number) + number(p, id, "solute_distance")? / 2.0;
            let q = [a.x, a.y, a.z];
            let lo: [i32; 3] = std::array::from_fn(|k| ((q[k] - r) / bx[k] * n[k] as f64).floor() as i32);
            let hi: [i32; 3] = std::array::from_fn(|k| ((q[k] + r) / bx[k] * n[k] as f64).ceil() as i32);
            for x in lo[0]..=hi[0] {
                for y in lo[1]..=hi[1] {
                    for z in lo[2]..=hi[2] {
                        let c = [(x as f64 + 0.5) * bx[0] / n[0] as f64, (y as f64 + 0.5) * bx[1] / n[1] as f64, (z as f64 + 0.5) * bx[2] / n[2] as f64];
                        if norm(sub(c, q)) > r {
                            continue;
                        }
                        let (wx, wy, wz) = (((x % n[0]) + n[0]) % n[0], ((y % n[1]) + n[1]) % n[1], ((z % n[2]) + n[2]) % n[2]);
                        taken[(wx as usize * n[1] as usize + wy as usize) * n[2] as usize + wz as usize] = true;
                    }
                }
            }
        }
        let (mut open, mut inside) = (0usize, 0usize);
        for x in 0..n[0] {
            for y in 0..n[1] {
                for z in 0..n[2] {
                    let c = [(x as f64 + 0.5) * bx[0] / n[0] as f64, (y as f64 + 0.5) * bx[1] / n[1] as f64, (z as f64 + 0.5) * bx[2] / n[2] as f64];
                    if region.contains(c) {
                        inside += 1;
                        if !taken[(x as usize * n[1] as usize + y as usize) * n[2] as usize + z as usize] {
                            open += 1;
                        }
                    }
                }
            }
        }
        let free = if inside > 0 { region.volume(bx) * open as f64 / inside as f64 } else { 0.0 };
        let mass_to_place = fmax(0.0, density * free / 1.66053906660 - fixed_mass);
        let mut per_set = 0.0;
        for k in &kinds {
            if k.count == 0 {
                per_set += molecule_mass(&k.atoms);
            }
        }
        let sets = (mass_to_place / per_set).floor() as i32;
        for k in kinds.iter_mut() {
            if k.count == 0 {
                k.count = sets;
            }
        }
    }
    let minimum = number(p, id, "min_distance")?;
    let solute_minimum = number(p, id, "solute_distance")?;
    let attempts = integer_of(p, id, "attempts", 1)?;
    let fixed_orientation = flag(p, id, "fixed_orientation")?;
    let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
    let mut reach = minimum;
    for v in [solute_minimum, 0.5] {
        if reach < v {
            reach = v;
        }
    }
    let g: [i32; 3] = std::array::from_fn(|k| 1.max((bx[k] / reach) as i32));
    let mut grid: HashMap<i64, Vec<(Vec3, bool)>> = HashMap::new();
    let key = |x: i32, y: i32, z: i32| {
        let (x, y, z) = (((x % g[0]) + g[0]) % g[0], ((y % g[1]) + g[1]) % g[1], ((z % g[2]) + g[2]) % g[2]);
        (x as i64 * g[1] as i64 + y as i64) * g[2] as i64 + z as i64
    };
    let wrapped = |mut q: Vec3| {
        for k in 0..3 {
            q[k] -= bx[k] * (q[k] / bx[k]).floor();
        }
        q
    };
    let cell_of_point = |q: Vec3, k: usize| (g[k] - 1).min((q[k] / bx[k] * g[k] as f64) as i32);
    let insert = |grid: &mut HashMap<i64, Vec<(Vec3, bool)>>, q: Vec3, solute: bool| {
        let w = wrapped(q);
        grid.entry(key(cell_of_point(w, 0), cell_of_point(w, 1), cell_of_point(w, 2))).or_default().push((w, solute));
    };
    let clear = |grid: &HashMap<i64, Vec<(Vec3, bool)>>, q: Vec3| {
        let (cx, cy, cz) = (cell_of_point(q, 0), cell_of_point(q, 1), cell_of_point(q, 2));
        let span = (fmax(minimum, solute_minimum) / (bx[0] / g[0] as f64)).ceil() as i32 + 1;
        for dx in -span..=span {
            for dy in -span..=span {
                for dz in -span..=span {
                    let Some(list) = grid.get(&key(cx + dx, cy + dy, cz + dz)) else { continue };
                    for &(r, solute) in list {
                        let mut v = sub(q, r);
                        for k in 0..3 {
                            v[k] -= bx[k] * (v[k] / bx[k]).round();
                        }
                        if norm(v) < if solute { solute_minimum } else { minimum } {
                            return false;
                        }
                    }
                }
            }
        }
        true
    };
    for a in &s.atoms {
        insert(&mut grid, [a.x, a.y, a.z], true);
    }
    stable_sort_by(&mut kinds, |a, b| a.atoms.len() > b.atoms.len());
    let solute_atoms = s.atoms.len();
    let mut summary = String::new();
    let mut placed_mass = 0.0;
    for kind in kinds.iter_mut() {
        let mut c = 0;
        while c < kind.count {
            let mut placed = false;
            let mut attempt = 0;
            while attempt < attempts && !placed {
                attempt += 1;
                let r = if fixed_orientation { identity() } else { random_rotation(&mut rng) };
                let centre = [uniform_real(&mut rng, 0.0, 1.0) * bx[0], uniform_real(&mut rng, 0.0, 1.0) * bx[1], uniform_real(&mut rng, 0.0, 1.0) * bx[2]];
                let mut trial = Vec::new();
                let mut ok = true;
                for a in &kind.atoms {
                    let raw = add(centre, rotate(&r, a.position));
                    if !region.contains(raw) || !clear(&grid, wrapped(raw)) {
                        ok = false;
                        break;
                    }
                    trial.push(raw);
                }
                if !ok {
                    continue;
                }
                for (k, q) in trial.iter().enumerate() {
                    add_atom(&mut s, &kind.atoms[k].symbol, *q);
                    insert(&mut grid, *q, false);
                }
                placed = true;
            }
            if !placed && kind.automatic {
                data.notes.push(format!(
                    "Placed {} of the {} {} molecules the density asks for (random insertion jammed; lower min_distance or relax the result)",
                    c, kind.count, kind.name
                ));
                kind.count = c;
                break;
            }
            if !placed {
                return Err(format!(
                    "Placed {} of {} {} molecules: lower the density or min_distance, enlarge the region, or raise attempts",
                    c, kind.count, kind.name
                ));
            }
            c += 1;
        }
        placed_mass += molecule_mass(&kind.atoms) * kind.count as f64;
        summary += &format!("{}{} {}", if summary.is_empty() { "" } else { ", " }, kind.count, kind.name);
    }
    if !solvate {
        set_cell(&mut s, [bx[0], 0.0, 0.0], [0.0, bx[1], 0.0], [0.0, 0.0, bx[2]], [true; 3]);
    }
    produce(data, s);
    if solvate {
        for i in solute_atoms..data.selected.len() {
            data.selected[i] = 1;
        }
    }
    let mut note = format!(
        "{}{} in {} x {} x {} Angstrom",
        summary,
        if solvate { format!(" around the {solute_atoms}-atom solute") } else { String::new() },
        g5(bx[0]),
        g5(bx[1]),
        g5(bx[2])
    );
    if region.kind != "all" {
        note += &format!(", inside the {} region", region.kind);
    }
    note += &format!(" ({} g/cm3 of molecules in the region)", g5(placed_mass * 1.66053906660 / region.volume(bx)));
    data.notes.push(note);
    Ok(())
}

fn t(id: &str, title: &str, category: &str, help: &str, parameters: Vec<crate::pipeline::ModifierParameter>, apply: fn(&Json, &mut PipelineData) -> R) -> ModifierType {
    ModifierType::new(id, title, category, help, parameters, apply)
}

pub fn types() -> Vec<ModifierType> {
    let layers = "graphene|hbn|silicene|germanene|stanene|mx2-2h|mx2-1t|phosphorene";
    let top_choices = "same|graphene|hbn|silicene|germanene|stanene|mx2-2h|mx2-1t";
    vec![
        t("build-nanotube", "Nanotube", "Build",
          "Rolls a honeycomb sheet into a tube with chiral indices (n, m): armchair (n, n), zigzag (n, 0) or chiral. Two elements make a BN-type tube (alternating sublattices). Multi-wall tubes (armchair or zigzag) add walls about wall_spacing apart. Periodic along the tube axis, or finite with the ends capped (H, F, ...). The notes give diameters, the chiral angle and, for carbon, metallic or semiconducting.",
          params(&[("n", "Chiral index n", "int", "6", ""), ("m", "Chiral index m", "int", "6", ""), ("length", "Unit cells along the axis", "int", "1", ""),
                   ("length_A", "Length in Angstrom (0 = use the unit cells)", "float", "0", ""),
                   ("elements", "Elements (one, or two for BN)", "elements", "C", ""), ("bond", "Bond length (0 = C-C 1.421, B-N 1.446)", "float", "0", ""),
                   ("walls", "Walls", "int", "1", ""), ("wall_spacing", "Wall spacing (Angstrom)", "float", "3.4", ""),
                   ("finite", "Finite tube (open ends)", "bool", "false", ""), ("end_cap", "End atoms of a finite tube", "choice", "H", "H|F|Cl|none"),
                   ("axis", "Tube axis", "choice", "z", "x|y|z"), ("vacuum", "Vacuum around the tube (Angstrom)", "float", "10", "")]),
          build_nanotube),
        t("build-monolayer", "2D layer", "Build",
          "A 2D layer or few-layer stack: graphene, hexagonal BN, buckled silicene, germanene and stanene, 2H or 1T MX2 dichalcogenides, or phosphorene; na x nb primitive cells (or the rectangular a x a sqrt3 cell of hexagonal layers), biaxial strain in %, and stacking of several layers (AA, AB Bernal, ABC rhombohedral, or 2H for MX2). a = 0 and thickness = 0 take the defaults.",
          params(&[("kind", "Layer", "choice", "graphene", layers), ("elements", "Elements (MX2: M X)", "elements", "", ""),
                   ("a", "Lattice constant a (0 = default)", "float", "0", ""),
                   ("thickness", "MX2 X-X thickness or X-ene buckling (0 = default)", "float", "0", ""),
                   ("strain", "Biaxial strain (%)", "float", "0", ""), ("na", "Cells along a", "int", "1", ""), ("nb", "Cells along b", "int", "1", ""),
                   ("rectangular", "Rectangular cell (hexagonal layers)", "bool", "false", ""), ("layers", "Layers", "int", "1", ""),
                   ("stacking", "Stacking", "choice", "AB", "AA|AB|ABC|2H"), ("interlayer", "Interlayer distance (0 = default)", "float", "0", ""),
                   ("periodic_c", "Periodic along c (bulk-like stack)", "bool", "true", ""), ("vacuum", "Vacuum (Angstrom)", "float", "15", "")]),
          build_monolayer),
        t("build-nanoribbon", "Nanoribbon", "Build",
          "A graphene (or BN) ribbon with armchair edges (width in dimer lines, N-AGNR) or zigzag edges (width in zigzag chains, N-ZGNR), periodic along x or a finite flake, with H, F or Cl on the edge atoms.",
          params(&[("edge", "Edge", "choice", "armchair", "armchair|zigzag"), ("width", "Width (dimer lines or chains)", "int", "7", ""),
                   ("length", "Periods along the ribbon", "int", "1", ""), ("elements", "Elements (one, or two for BN)", "elements", "C", ""),
                   ("bond", "Bond length (Angstrom)", "float", "1.421", ""), ("passivate", "Terminate the edges", "bool", "true", ""),
                   ("termination", "Edge atoms", "choice", "H", "H|F|Cl|none"), ("termination_bond", "Edge bond length (0 = default)", "float", "0", ""),
                   ("finite", "Finite flake (open ends)", "bool", "false", ""), ("vacuum", "Vacuum (Angstrom)", "float", "10", "")]),
          build_nanoribbon),
        t("build-twisted-bilayer", "Twisted bilayer", "Build",
          "Two hexagonal layers turned by a commensurate angle (cos t = (3i^2+3i+1/2)/(3i^2+3i+1): 21.79, 13.17, 9.43, ... degrees), in their moire supercell. Give the angle (the nearest commensurate one is used) or the index i. The top layer may be another layer, strained to the bottom lattice (a heterobilayer, the strain reported), and start in AA or AB registry.",
          params(&[("kind", "Layer", "choice", "graphene", layers), ("elements", "Elements (MX2: M X)", "elements", "", ""),
                   ("angle", "Twist angle (degrees; 0 = use the index)", "float", "0", ""), ("index", "Commensurate index i", "int", "3", ""),
                   ("a", "Lattice constant a (0 = default)", "float", "0", ""), ("thickness", "MX2 X-X thickness (0 = 3.19)", "float", "0", ""),
                   ("top_kind", "Top layer", "choice", "same", top_choices),
                   ("top_elements", "Top layer elements", "elements", "", ""), ("stacking", "Registry before the twist", "choice", "AA", "AA|AB"),
                   ("interlayer", "Interlayer distance (0 = default)", "float", "0", ""), ("vacuum", "Vacuum (Angstrom)", "float", "15", "")]),
          build_twisted_bilayer),
        t("build-cluster", "Magic-number cluster", "Build",
          "Closed-shell clusters: Mackay icosahedra and fcc cuboctahedra (1, 13, 55, 147, ... atoms for 1, 2, 3, 4 shells), regular decahedra (pentagonal bipyramids of five tetrahedra: 7, 23, 54, ... atoms for 2, 3, 4 atoms per edge) and fcc octahedra (shells = atoms per edge: 6, 19, 44, 85 atoms), with the nearest-neighbour distance along the edges. Nanoalloys: a fraction of the atoms becomes a second element, at random, in the core, in the shell (outermost atoms first) or on one side (Janus, along z).",
          params(&[("kind", "Shape", "choice", "icosahedron", "icosahedron|decahedron|cuboctahedron|octahedron"), ("element", "Element", "element", "Cu", ""),
                   ("shells", "Shells (1 = one atom)", "int", "3", ""), ("distance", "Nearest-neighbour distance (Angstrom)", "float", "2.556", ""),
                   ("second_element", "Alloying element (empty = none)", "string", "", ""), ("fraction", "Alloying fraction", "float", "0.5", ""),
                   ("distribution", "Alloying element in", "choice", "random", "random|core|shell|janus"), ("seed", "Seed (random)", "int", "1", ""),
                   ("vacuum", "Vacuum (Angstrom)", "float", "10", "")]),
          build_cluster),
        t("build-molecule", "Molecule", "Build",
          "A molecule from the library (experimental geometries): H2, N2, O2, F2, Cl2, CO, NO, HF, HCl, H2O, H2S, CO2, SO2, NO2, O3, NH3, PH3, CH4, SiH4, CF4, C2H2, C2H4, C2H6, CH3OH, HCN, H2CO, C6H6, C60 (or common names: water, benzene, ...), or any structure file (XYZ, POSCAR, CIF, ...). Orientation: as in the library or file, along the principal axes of inertia (long axis on z, a planar molecule in the xz plane), or random. With vacuum, in a box (periodic or open).",
          params(&[("name", "Molecule (name or file)", "string", "H2O", ""), ("orientation", "Orientation", "choice", "library", "library|principal|random"),
                   ("seed", "Seed (random orientation)", "int", "1", ""), ("vacuum", "Vacuum (Angstrom; 0 = no cell)", "float", "10", ""),
                   ("periodic", "Periodic box", "bool", "false", "")]),
          build_molecule),
        t("build-packing", "Molecular packing", "Build",
          "Packs molecules (library names, single atoms or structure files) at random positions and orientations into a periodic box, keeping atoms of different molecules at least min_distance apart (as Packmol). The box follows from the density unless given. Solvate mode fills the input structure's orthogonal cell around it, keeping solute_distance from its atoms. A count of 0 fills the region to the density. Regions (all atoms of a molecule inside): box (low, high corners), sphere or outside-sphere (centre, radius), slab (low, high along the axis) or cylinder (centre, radius, along the axis).",
          params(&[("molecules", "Molecules and counts (e.g. H2O 100 CO2 10; 0 = to the density)", "list", "H2O 64", ""),
                   ("mode", "Mode", "choice", "new", "new|solvate"),
                   ("density", "Density (g/cm3)", "float", "1.0", ""), ("box", "Box lengths (Angstrom; 0 0 0 = from the density)", "vector", "0 0 0", ""),
                   ("min_distance", "Minimum distance between molecules (Angstrom)", "float", "2.0", ""),
                   ("solute_distance", "Minimum distance to the solute (Angstrom)", "float", "2.5", ""),
                   ("region", "Region", "choice", "all", "all|box|sphere|outside-sphere|slab|cylinder"),
                   ("region_low", "Region low corner (Angstrom)", "vector", "0 0 0", ""), ("region_high", "Region high corner (Angstrom)", "vector", "0 0 0", ""),
                   ("region_center", "Region centre (Angstrom)", "vector", "0 0 0", ""), ("region_radius", "Region radius (Angstrom)", "float", "0", ""),
                   ("region_axis", "Region axis (slab, cylinder)", "choice", "z", "x|y|z"),
                   ("fixed_orientation", "Keep the library orientation", "bool", "false", ""), ("seed", "Seed", "int", "1", ""),
                   ("attempts", "Attempts per molecule", "int", "20000", "")]),
          build_packing),
        t("add-adsorbate", "Add adsorbate", "Modification",
          "Places a molecule (library name, element or structure file) on a slab surface (normal a x b; top towards c, bottom, or both) at top, bridge, fcc-hollow or hcp-hollow sites, its anchor atom (or lowest atom) `height` above the outer layer. Coverage: one adsorbate at the site nearest the fractional position (u, v), all symmetry-distinct sites of the kind, or a random fraction of them kept min_separation apart. The molecule's z axis points away from the surface, after an optional principal-axes alignment, a tilt about y and an azimuthal turn about z.",
          params(&[("species", "Molecule, element or file", "string", "CO", ""), ("site", "Site", "choice", "top", "top|bridge|fcc|hcp"),
                   ("height", "Height above the surface (Angstrom)", "float", "1.8", ""), ("position", "Near fractional position (u v, w ignored)", "vector", "0.5 0.5 0", ""),
                   ("surface", "Surface", "choice", "top", "top|bottom|both"), ("coverage", "Coverage", "choice", "one", "one|all|fraction"),
                   ("fraction", "Fraction of the sites", "float", "0.5", ""), ("min_separation", "Minimum adsorbate separation (Angstrom)", "float", "0", ""),
                   ("seed", "Seed", "int", "1", ""), ("anchor", "Anchor atom of the molecule (0-based; -1 = lowest)", "int", "-1", ""),
                   ("principal", "Align the principal axes first", "bool", "false", ""), ("tilt", "Tilt from the normal (degrees)", "float", "0", ""),
                   ("azimuth", "Azimuth (degrees)", "float", "0", ""),
                   ("select", "Select the adsorbate", "bool", "true", "")]),
          add_adsorbate),
        t("near-cubic-supercell", "Near-cubic supercell", "Cell",
          "The supercell (integer matrix with entries up to +-range) of about target_atoms atoms (or any count within the tolerance) whose cell is closest to a cube, an fcc or bcc primitive cell (as ASE's find_optimal_cell_shape), or that has the largest shortest image distance, for point-defect calculations. The notes give the shortest image distance and the inscribed sphere.",
          params(&[("target_atoms", "Target number of atoms", "int", "64", ""), ("shape", "Target shape", "choice", "sc", "sc|fcc|bcc|image"),
                   ("tolerance", "Atom-count tolerance (%)", "float", "0", ""), ("range", "Matrix entries up to +-", "int", "2", "")]),
          near_cubic_supercell),
        t("interpolate", "Interpolate to a structure", "Modification",
          "Moves every atom the fraction t of the way to its position in another structure with the same atoms in the same order (minimum image across periodic boundaries); the cell is interpolated too. Images for a nudged elastic band are t = k/(n+1). IDPP (image-dependent pair potential, Smidstrup et al. 2014) then relaxes the linear guess towards interpolated interatomic distances, avoiding atoms passing too close (better NEB starting paths).",
          params(&[("file", "Final structure file", "string", "", ""), ("fraction", "Fraction t (0 = this, 1 = final)", "float", "0.5", ""),
                   ("method", "Method", "choice", "linear", "linear|idpp"),
                   ("idpp_cutoff", "IDPP pair cutoff (Angstrom; 0 = all pairs)", "float", "0", ""), ("idpp_steps", "IDPP steps", "int", "2000", "")]),
          interpolate),
        t("create-antisite", "Antisite defect", "Modification",
          "Swaps atoms of the first element (chosen at random) with the nearest (or a random) atom of the second, count times, every swapped atom at least min_separation from the others (isolated defects); the swapped atoms are selected.",
          params(&[("elements", "Elements A B", "elements", "Ni Al", ""), ("partner", "Partner", "choice", "nearest", "nearest|random"),
                   ("count", "Antisite pairs", "int", "1", ""), ("min_separation", "Minimum separation between defects (Angstrom)", "float", "0", ""),
                   ("seed", "Seed", "int", "1", "")]),
          antisite),
        t("create-frenkel-pair", "Frenkel pair", "Modification",
          "Moves an atom of the element into an interstitial void (as Insert interstitial atoms) or into a split (dumbbell) interstitial on a lattice site, e.g. <100> in fcc or <110> in bcc metals, leaving a vacancy between min_distance and max_distance away; the interstitial atoms are selected.",
          params(&[("element", "Element", "element", "Cu", ""), ("interstitial", "Interstitial", "choice", "void", "void|dumbbell"),
                   ("direction", "Dumbbell direction", "vector", "1 0 0", ""), ("lattice_direction", "Lattice direction [u v w]", "bool", "true", ""),
                   ("dumbbell_length", "Dumbbell length (Angstrom; 0 = 0.9 x nearest neighbour)", "float", "0", ""),
                   ("min_distance", "Minimum vacancy-interstitial distance (Angstrom)", "float", "5", ""),
                   ("max_distance", "Maximum vacancy-interstitial distance (Angstrom; 0 = any)", "float", "0", ""),
                   ("seed", "Seed", "int", "1", "")]),
          frenkel_pair),
    ]
}

