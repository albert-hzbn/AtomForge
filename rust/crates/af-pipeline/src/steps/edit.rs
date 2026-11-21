//! Geometric and chemical edits (the C++ pipeline/EditSteps): mirror,
//! inversion, rotation and alignment, sorting, hydrogen saturation, setting
//! and measuring bond geometry, duplicating the selection, cell reduction and
//! scaling, swapping elements, making molecules whole, selecting fragments,
//! per-atom simulation data and stacking layers.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::Path;

use af_core::cxx::{sort_by, stable_sort_by};
use af_core::cxx_random::{shuffle, Mt19937_64, Normal};
use af_core::elements;
use af_core::format::format_g;
use af_core::json::Json;
use af_science::core::{add, cross, determinant, dot, inverse, multiply, norm, row_times, scale, sub, Mat3, Vec3};
use af_science::data::read_frames;

use crate::pipeline::{params, ModifierParameter, ModifierType, PipelineData};
use crate::steps::basic::cxx_u64;
use crate::structure_edits::{append_atom, bond_neighbours, rotate, rotation_about, rotation_between, saturate, Neighbours, SaturationOptions};
use crate::support::*;

type R = Result<(), String>;

const PI: f64 = 3.14159265358979323846;
const DEGREE: f64 = PI / 180.0;

/// ostream << double at precision 6.
fn g(x: f64) -> String {
    format_g(x, 6)
}

extern "C" {
    #[link_name = "remainder"]
    fn c_remainder(x: f64, y: f64) -> f64;
}

/// std::remainder (IEEE remainder).
fn remainder(x: f64, y: f64) -> f64 {
    unsafe { c_remainder(x, y) }
}

fn unit(v: Vec3) -> Result<Vec3, String> {
    let n = norm(v);
    if !(n > 1e-12) {
        return Err("A direction must not be zero".into());
    }
    Ok(scale(v, 1.0 / n))
}

fn wrap_all(s: &mut af_core::Structure) -> R {
    if !s.has_unit_cell {
        return Ok(());
    }
    let cell = cell_of(s);
    let inv = inverse(&cell)?;
    let (offset, pbc) = (s.cell_offset, s.pbc);
    for a in s.atoms.iter_mut() {
        let mut f = row_times(sub(position_of(a), offset), &inv);
        for k in 0..3 {
            if pbc[k] {
                f[k] -= f[k].floor();
            }
        }
        place(a, add(row_times(f, &cell), offset));
    }
    Ok(())
}

fn centroid(s: &af_core::Structure, atoms: &[usize]) -> Vec3 {
    let mut c = [0.0; 3];
    for &i in atoms {
        c = add(c, position_of(&s.atoms[i]));
    }
    if atoms.is_empty() {
        c
    } else {
        scale(c, 1.0 / atoms.len() as f64)
    }
}

fn cell_centre(s: &af_core::Structure) -> Vec3 {
    let cell = cell_of(s);
    add(s.cell_offset, scale(add(add(cell[0], cell[1]), cell[2]), 0.5))
}

fn pivot_of(p: &Json, id: &str, data: &PipelineData, atoms: &[usize]) -> Result<Vec3, String> {
    let centre = text(p, id, "center")?;
    match centre.as_str() {
        "origin" => Ok([0.0; 3]),
        "point" => vector(p, id, "point"),
        "cell" => {
            require_cell(data, "The cell centre")?;
            Ok(cell_centre(&data.structure))
        }
        _ => Ok(centroid(&data.structure, atoms)),
    }
}

/// strtol(w, &end, 10) on Windows (32-bit long): (value, whole word parsed).
pub(crate) fn strtol(w: &str) -> (i64, bool) {
    let b = w.as_bytes();
    let mut i = 0;
    let negative = match b.first() {
        Some(b'-') => {
            i = 1;
            true
        }
        Some(b'+') => {
            i = 1;
            false
        }
        _ => false,
    };
    let start = i;
    let mut value: i64 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        value = (value * 10 + (b[i] - b'0') as i64).min(1 << 40);
        i += 1;
    }
    if i == start {
        return (0, w.is_empty());
    }
    let v = if negative { -value } else { value };
    (v.clamp(i32::MIN as i64, i32::MAX as i64), i == b.len())
}

fn parse_indices(list: &str, atoms: usize, id: &str) -> Result<Vec<usize>, String> {
    let mut out = Vec::new();
    for w in list.split_ascii_whitespace() {
        let (value, whole) = strtol(w);
        if !whole || value < 0 || value as usize >= atoms {
            return Err(format!("{id}: '{w}' is not an atom index (0 to {})", if atoms > 0 { atoms - 1 } else { 0 }));
        }
        out.push(value as usize);
    }
    Ok(out)
}

// ---------------------------------------------------------------- symmetry edits
fn mirror(p: &Json, data: &mut PipelineData) -> R {
    let id = "mirror";
    let selected_only = text(p, id, "target")? == "selected";
    let atoms = targets(data, selected_only, id)?;
    let n = unit(vector(p, id, "normal")?)?;
    let point = pivot_of(p, id, data, &atoms)?;
    let reflect = |x: Vec3| sub(x, scale(n, 2.0 * dot(sub(x, point), n)));
    for &i in &atoms {
        let q = reflect(position_of(&data.structure.atoms[i]));
        place(&mut data.structure.atoms[i], q);
    }
    let s = &mut data.structure;
    if !selected_only && s.has_unit_cell {
        for row in s.cell_vectors.iter_mut() {
            *row = sub(*row, scale(n, 2.0 * dot(*row, n)));
        }
        s.cell_vectors[2] = scale(s.cell_vectors[2], -1.0);
        s.cell_offset = reflect(s.cell_offset);
        s.cell_offset = sub(s.cell_offset, s.cell_vectors[2]);
        wrap_all(s)?;
    } else if flag(p, id, "wrap")? {
        wrap_all(s)?;
    }
    Ok(())
}

fn invert(p: &Json, data: &mut PipelineData) -> R {
    let id = "invert";
    let atoms = targets(data, text(p, id, "target")? == "selected", id)?;
    let point = pivot_of(p, id, data, &atoms)?;
    for &i in &atoms {
        let q = sub(scale(point, 2.0), position_of(&data.structure.atoms[i]));
        place(&mut data.structure.atoms[i], q);
    }
    if flag(p, id, "wrap")? {
        wrap_all(&mut data.structure)?;
    }
    Ok(())
}

fn euler_rotation(angles: Vec3, convention: &str) -> Result<Mat3, String> {
    let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let (a, b, c) = (angles[0] * DEGREE, angles[1] * DEGREE, angles[2] * DEGREE);
    match convention {
        "zyz" => Ok(multiply(&rotation_about(z, a), &multiply(&rotation_about(y, b), &rotation_about(z, c)))),
        "zxz" => Ok(multiply(&rotation_about(z, a), &multiply(&rotation_about(x, b), &rotation_about(z, c)))),
        "xyz" => Ok(multiply(&rotation_about(z, c), &multiply(&rotation_about(y, b), &rotation_about(x, a)))),
        _ => Err("Euler convention must be zyz, zxz or xyz".into()),
    }
}

fn rotate_step(p: &Json, data: &mut PipelineData) -> R {
    let id = "rotate";
    let selected_only = text(p, id, "target")? == "selected";
    let atoms = targets(data, selected_only, id)?;
    let mode = text(p, id, "mode")?;
    if mode != "axis" && mode != "euler" {
        return Err("rotate: mode is axis or euler".into());
    }
    let r = if mode == "euler" {
        let convention = text(p, id, "convention")?;
        euler_rotation(vector(p, id, "euler")?, &convention)?
    } else {
        let angle = number(p, id, "angle")? * DEGREE;
        rotation_about(vector(p, id, "axis")?, angle)
    };
    let with_cell = !selected_only && flag(p, id, "rotate_cell")? && data.structure.has_unit_cell;
    let point = if with_cell { data.structure.cell_offset } else { pivot_of(p, id, data, &atoms)? };
    for &i in &atoms {
        let q = add(point, rotate(&r, sub(position_of(&data.structure.atoms[i]), point)));
        place(&mut data.structure.atoms[i], q);
    }
    if with_cell {
        for row in data.structure.cell_vectors.iter_mut() {
            *row = rotate(&r, *row);
        }
    }
    if mode == "euler" {
        let trace = r[0][0] + r[1][1] + r[2][2];
        let angle = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0).acos() / DEGREE;
        let mut axis = [r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]];
        if norm(axis) > 1e-12 {
            axis = scale(axis, 1.0 / norm(axis));
        }
        data.notes.push(format!("Equivalent to a turn of {} degrees about ({}, {}, {})", g(angle), g(axis[0]), g(axis[1]), g(axis[2])));
    }
    Ok(())
}

fn align_direction(p: &Json, data: &mut PipelineData) -> R {
    let id = "align-direction";
    let lattice = flag(p, id, "lattice")?;
    let cartesian = |data: &PipelineData, v: Vec3| -> Result<Vec3, String> {
        if !lattice {
            return Ok(v);
        }
        require_cell(data, "A lattice direction")?;
        Ok(row_times(v, &cell_of(&data.structure)))
    };
    let axis_of = |name: &str| -> Vec3 {
        match name {
            "x" => [1.0, 0.0, 0.0],
            "y" => [0.0, 1.0, 0.0],
            _ => [0.0, 0.0, 1.0],
        }
    };
    let d = unit(cartesian(data, vector(p, id, "direction")?)?)?;
    let e = axis_of(&text(p, id, "axis")?);
    let mut r = rotation_between(d, e);
    let second_axis = text(p, id, "second_axis")?;
    if second_axis != "none" {
        let f = axis_of(&second_axis);
        if norm(cross(e, f)) < 0.5 {
            return Err("The second axis must differ from the first".into());
        }
        let mut t = rotate(&r, unit(cartesian(data, vector(p, id, "second_direction")?)?)?);
        t = sub(t, scale(e, dot(t, e)));
        if norm(t) < 1e-6 {
            return Err("The second direction is parallel to the first".into());
        }
        let angle = dot(cross(unit(t)?, f), e).atan2(dot(unit(t)?, f));
        r = multiply(&rotation_about(e, angle), &r);
    }
    let s = &mut data.structure;
    for a in s.atoms.iter_mut() {
        let q = rotate(&r, position_of(a));
        place(a, q);
    }
    if s.has_unit_cell {
        for row in s.cell_vectors.iter_mut() {
            *row = rotate(&r, *row);
        }
        s.cell_offset = rotate(&r, s.cell_offset);
    }
    Ok(())
}

// ---------------------------------------------------------------- ordering
fn sort_atoms(p: &Json, data: &mut PipelineData) -> R {
    let id = "sort-atoms";
    let c = centroid(&data.structure, &targets(data, false, id)?);
    let mut rank: BTreeMap<i32, i32> = BTreeMap::new();
    for (k, w) in text(p, id, "element_order")?.replace(',', " ").split_ascii_whitespace().enumerate() {
        rank.entry(af_science::data::atomic_number(&standard_element(w)?)).or_insert(k as i32);
    }
    let s = &data.structure;
    let values = |key: &str| -> Result<Vec<f64>, String> {
        let mut value = vec![0.0; s.atoms.len()];
        for (i, a) in s.atoms.iter().enumerate() {
            value[i] = match key {
                "none" => 0.0,
                "element" => match rank.get(&a.atomic_number) {
                    Some(&r) => r as f64 - 1000.0,
                    None => a.atomic_number as f64,
                },
                "x" => a.x,
                "y" => a.y,
                "z" => a.z,
                "distance" => norm(sub(position_of(a), c)),
                "magmom" => a.magmom,
                "charge" => a.charge,
                "mass" => elements::atomic_mass(a.atomic_number),
                "selected" => {
                    if i < data.selected.len() && data.selected[i] != 0 {
                        0.0
                    } else {
                        1.0
                    }
                }
                "property" => {
                    if s.atom_property.len() != s.atoms.len() {
                        return Err("sort-atoms by property needs a per-atom property".into());
                    }
                    if s.atom_property[i].is_nan() {
                        f64::INFINITY
                    } else {
                        s.atom_property[i]
                    }
                }
                _ => return Err(format!("Unknown sort key {key}")),
            };
        }
        Ok(value)
    };
    // GCC evaluates `first = values(key), second = values(then)` in declaration order.
    let first = values(&text(p, id, "key")?)?;
    let second = values(&text(p, id, "then")?)?;
    let tolerance = number(p, id, "tolerance")?;
    let descending = flag(p, id, "descending")?;
    let mut order: Vec<usize> = (0..s.atoms.len()).collect();
    stable_sort_by(&mut order, |&x, &y| {
        if (first[x] - first[y]).abs() > tolerance {
            return if descending { first[x] > first[y] } else { first[x] < first[y] };
        }
        second[x] < second[y]
    });
    rebuild(data, &order);
    Ok(())
}

// ---------------------------------------------------------------- chemistry
fn saturate_hydrogen(p: &Json, data: &mut PipelineData) -> R {
    let id = "saturate-hydrogen";
    let mut elements_set = BTreeSet::new();
    for w in text(p, id, "elements")?.split_ascii_whitespace() {
        elements_set.insert(standard_element(w)?);
    }
    if data.selected.len() != data.structure.atoms.len() {
        data.selected = vec![0; data.structure.atoms.len()];
    }
    let mask = data.selected.clone();
    let selected_target = text(p, id, "target")? == "selected";
    if selected_target && !mask.iter().any(|&c| c != 0) {
        return Err("saturate-hydrogen: nothing is selected".into());
    }
    let options = SaturationOptions {
        terminator: standard_element(&text(p, id, "terminator")?)?,
        bond_length: number(p, id, "bond_length")?,
        max_per_atom: integer_of(p, id, "max_per_atom", 1)?,
        mask: if selected_target { Some(&mask) } else { None },
    };
    let select_new = flag(p, id, "select")?;
    if select_new {
        data.selected.iter_mut().for_each(|c| *c = 0);
    }
    let bond_scale = number(p, id, "bond_scale")?;
    let sp3 = text(p, id, "hybridization")? == "sp3";
    let added = saturate(&mut data.structure, &mut data.selected, &elements_set, sp3, bond_scale, select_new, &options)?;
    let name = if options.terminator == "H" { "hydrogen".to_string() } else { options.terminator.clone() };
    data.notes.push(format!("{added} {name} atoms added"));
    Ok(())
}

fn swap_elements(p: &Json, data: &mut PipelineData) -> R {
    let id = "swap-elements";
    let a = standard_element(&text(p, id, "first")?)?;
    let b = standard_element(&text(p, id, "second")?)?;
    let (mut na, mut nb) = (0usize, 0usize);
    for atom in data.structure.atoms.iter_mut() {
        if atom.symbol == a {
            atom.symbol = b.clone();
            na += 1;
        } else if atom.symbol == b {
            atom.symbol = a.clone();
            nb += 1;
        } else {
            continue;
        }
        atom.atomic_number = af_science::data::atomic_number(&atom.symbol);
        [atom.r, atom.g, atom.b] = elements::default_colour(atom.atomic_number);
    }
    data.notes.push(format!("{na} {a} -> {b}, {nb} {b} -> {a}"));
    Ok(())
}

// ---------------------------------------------------------------- bond geometry
fn fragment(nb: &Neighbours, start: usize, barrier: usize) -> Result<Vec<(usize, Vec3)>, String> {
    let mut out = vec![(start, [0.0; 3])];
    let mut seen = vec![false; nb.atoms.len()];
    seen[start] = true;
    let mut queue = VecDeque::from([(start, [0.0; 3])]);
    while let Some((atom, offset)) = queue.pop_front() {
        for k in 0..nb.atoms[atom].len() {
            let next = nb.atoms[atom][k];
            if atom == start && next == barrier {
                continue;
            }
            if next == barrier {
                return Err("The atoms are in a ring: the side to move is bonded to the fixed side (use move = atom)".into());
            }
            if seen[next] {
                continue;
            }
            seen[next] = true;
            let to = add(offset, nb.vectors[atom][k]);
            out.push((next, to));
            queue.push_back((next, to));
        }
    }
    Ok(out)
}

/// Unwrapped positions of the listed atoms (minimum image from the first).
fn geometry_of(s: &af_core::Structure, atoms: &[usize]) -> Result<Vec<Vec3>, String> {
    let mut previous = position_of(&s.atoms[atoms[0]]);
    let mut points = vec![previous];
    for k in 1..atoms.len() {
        let mut d = sub(position_of(&s.atoms[atoms[k]]), position_of(&s.atoms[atoms[k - 1]]));
        if s.has_unit_cell {
            d = minimum_image(s, d)?;
        }
        previous = add(previous, d);
        points.push(previous);
    }
    Ok(points)
}

fn angle_of(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let (u, v) = (sub(a, b), sub(c, b));
    (dot(u, v) / (norm(u) * norm(v))).clamp(-1.0, 1.0).acos() / DEGREE
}

fn dihedral_of(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f64 {
    let (b1, b2, b3) = (sub(b, a), sub(c, b), sub(d, c));
    let (n1, n2) = (cross(b1, b2), cross(b2, b3));
    let m = cross(n1, scale(b2, 1.0 / norm(b2)));
    dot(m, n2).atan2(dot(n1, n2)) / DEGREE
}

fn measure(q: &[Vec3]) -> f64 {
    match q.len() {
        2 => norm(sub(q[1], q[0])),
        3 => angle_of(q[0], q[1], q[2]),
        _ => dihedral_of(q[0], q[1], q[2], q[3]),
    }
}

fn quantity(atoms: usize) -> &'static str {
    match atoms {
        2 => "distance",
        3 => "angle",
        _ => "dihedral",
    }
}

fn measure_geometry(p: &Json, data: &mut PipelineData) -> R {
    let id = "measure-geometry";
    let list = text(p, id, "atoms")?;
    let atoms = parse_indices(&list, data.structure.atoms.len(), id)?;
    if atoms.len() < 2 || atoms.len() > 4 {
        return Err("measure-geometry: give 2, 3 or 4 atom indices".into());
    }
    let value = measure(&geometry_of(&data.structure, &atoms)?);
    let unit_name = if atoms.len() == 2 { " Angstrom" } else { " degrees" };
    data.notes.push(format!("{} {} = {}{}", quantity(atoms.len()), list, g(value), unit_name));
    Ok(())
}

fn set_geometry(p: &Json, data: &mut PipelineData) -> R {
    let id = "set-geometry";
    let atoms = parse_indices(&text(p, id, "atoms")?, data.structure.atoms.len(), id)?;
    if atoms.len() < 2 || atoms.len() > 4 {
        return Err("set-geometry: give 2 (distance), 3 (angle) or 4 (dihedral) atom indices".into());
    }
    for x in 0..atoms.len() {
        for y in x + 1..atoms.len() {
            if atoms[x] == atoms[y] {
                return Err("set-geometry: the atoms must differ".into());
            }
        }
    }
    let target = number(p, id, "value")?;
    let s = &data.structure;
    let gpoints = geometry_of(s, &atoms)?;
    let current = measure(&gpoints);
    let four = atoms.len() == 4;
    let mover = if four { atoms[2] } else { *atoms.last().unwrap() };
    let fixed = if four { atoms[1] } else { atoms[atoms.len() - 2] };
    let mv = text(p, id, "move")?;
    if mv != "fragment" && mv != "atom" && mv != "both" {
        return Err("set-geometry: move is fragment, atom or both".into());
    }
    let nb = if mv == "atom" { Neighbours::default() } else { bond_neighbours(s, number(p, id, "bond_scale")?)? };
    let mut moving = if mv != "atom" { fragment(&nb, mover, fixed)? } else { vec![(mover, [0.0; 3])] };
    if four && mv == "atom" {
        moving = vec![(atoms[3], sub(gpoints[3], gpoints[2]))];
    }
    let mover_point = if four { gpoints[2] } else { *gpoints.last().unwrap() };
    let other_index = if four { 1 } else { 0 };
    let other = if mv == "both" { fragment(&nb, atoms[other_index], if four { atoms[2] } else { atoms[1] })? } else { Vec::new() };
    let other_point = gpoints[other_index];
    for &(atom, _) in &moving {
        for k in 0..atoms.len() {
            if atoms[k] == atom && atom != mover && !(four && k == 3) {
                return Err(format!("set-geometry: the moving side contains atom {atom} of the fixed side"));
            }
        }
        if other.iter().any(|o| o.0 == atom) {
            return Err("set-geometry: the two sides are bonded to each other (a ring); use move = atom".into());
        }
    }
    let share = if mv == "both" { 0.5 } else { 1.0 };
    let result: Vec<(usize, Vec3)> = if atoms.len() == 2 {
        if !(target > 0.0) {
            return Err("set-geometry: a distance must be positive".into());
        }
        let direction = unit(sub(gpoints[1], gpoints[0]))?;
        let mut r = Vec::new();
        for &(atom, offset) in &moving {
            r.push((atom, add(add(mover_point, offset), scale(direction, share * (target - current)))));
        }
        for &(atom, offset) in &other {
            r.push((atom, add(add(other_point, offset), scale(direction, -(1.0 - share) * (target - current)))));
        }
        r
    } else {
        let (pivot, axis) = if atoms.len() == 3 {
            let mut axis = cross(sub(gpoints[0], gpoints[1]), sub(gpoints[2], gpoints[1]));
            if norm(axis) < 1e-8 {
                let d01 = sub(gpoints[0], gpoints[1]);
                let trial = if (gpoints[0][0] - gpoints[1][0]).abs() < 0.9 * norm(d01) { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                axis = cross(d01, trial);
            }
            (gpoints[1], axis)
        } else {
            (gpoints[2], sub(gpoints[2], gpoints[1]))
        };
        let mut delta = target - current;
        if four {
            delta = remainder(delta, 360.0);
        }
        let turned = |degrees: f64| {
            let forward = rotation_about(axis, share * degrees * DEGREE);
            let back = rotation_about(axis, -(1.0 - share) * degrees * DEGREE);
            let mut out = Vec::new();
            for &(atom, offset) in &moving {
                out.push((atom, add(pivot, rotate(&forward, sub(add(mover_point, offset), pivot)))));
            }
            for &(atom, offset) in &other {
                out.push((atom, add(pivot, rotate(&back, sub(add(other_point, offset), pivot)))));
            }
            out
        };
        let check = |moved: &[(usize, Vec3)]| {
            let mut h = gpoints.clone();
            for &(atom, position) in moved {
                for k in 0..atoms.len() {
                    if atoms[k] == atom {
                        h[k] = position;
                    }
                }
            }
            let value = measure(&h);
            if four {
                remainder(value - target, 360.0).abs()
            } else {
                (value - target).abs()
            }
        };
        let (forward, backward) = (turned(delta), turned(-delta));
        let chosen = if check(&forward) <= check(&backward) { forward } else { backward };
        if atoms.len() == 3 && (target < 0.0 || target > 180.0) {
            return Err("set-geometry: an angle lies between 0 and 180 degrees".into());
        }
        chosen
    };
    for &(atom, position) in &result {
        place(&mut data.structure.atoms[atom], position);
    }
    let after = measure(&geometry_of(&data.structure, &atoms)?);
    data.notes.push(format!("{} {} -> {}; moved {} atom(s)", quantity(atoms.len()), g(current), g(after), result.len()));
    Ok(())
}

// ---------------------------------------------------------------- duplication
fn duplicate_selected(p: &Json, data: &mut PipelineData) -> R {
    let id = "duplicate-selected";
    let atoms = targets(data, true, id)?;
    let offset = vector(p, id, "offset")?;
    let copies = number(p, id, "copies")? as i32;
    if copies < 1 {
        return Err("duplicate-selected: copies must be at least 1".into());
    }
    let select_new = flag(p, id, "select")?;
    if select_new {
        data.selected.iter_mut().for_each(|c| *c = 0);
    }
    for c in 1..=copies {
        for &i in &atoms {
            let original = data.structure.atoms[i].clone();
            append_atom(&mut data.structure, &mut data.selected, &original.symbol, add(position_of(&original), scale(offset, c as f64)), select_new);
            let copy = data.structure.atoms.last_mut().unwrap();
            copy.r = original.r;
            copy.g = original.g;
            copy.b = original.b;
            copy.fixed = original.fixed;
        }
    }
    data.notes.push(format!("{} atoms added", atoms.len() * copies as usize));
    Ok(())
}

// ---------------------------------------------------------------- cell
/// Niggli reduction (Krivy & Gruber 1976, tolerances of Grosse-Kunstleve,
/// Sauter & Adams 2004) of the cell vectors.
fn niggli(mut v: Mat3, eps: f64) -> Result<Mat3, String> {
    let dot_rows = |v: &Mat3, i: usize, j: usize| dot(v[i], v[j]);
    let flip_for = |v: &mut Mat3, positive: bool| {
        let (xi, eta, zeta) = (2.0 * dot_rows(v, 1, 2), 2.0 * dot_rows(v, 0, 2), 2.0 * dot_rows(v, 0, 1));
        for pass in 0..2 {
            for mask in 0..8 {
                let sa = if mask & 1 != 0 { -1.0 } else { 1.0 };
                let sb = if mask & 2 != 0 { -1.0 } else { 1.0 };
                let sc = if mask & 4 != 0 { -1.0 } else { 1.0 };
                if pass == 0 && sa * sb * sc < 0.0 {
                    continue;
                }
                let (x, y, z) = (sb * sc * xi, sa * sc * eta, sa * sb * zeta);
                let ok = if positive { x > eps && y > eps && z > eps } else { x <= eps && y <= eps && z <= eps };
                if !ok {
                    continue;
                }
                v[0] = scale(v[0], sa);
                v[1] = scale(v[1], sb);
                v[2] = scale(v[2], sc);
                if sa * sb * sc < 0.0 {
                    for row in v.iter_mut() {
                        *row = scale(*row, -1.0);
                    }
                }
                return;
            }
        }
    };
    let sign = |x: f64| if x > 0.0 { 1.0 } else { -1.0 };
    let sign_of = |x: f64| if x > eps { 1 } else if x < -eps { -1 } else { 0 };
    for _ in 0..1000 {
        let (a, b) = (dot_rows(&v, 0, 0), dot_rows(&v, 1, 1));
        let c = dot_rows(&v, 2, 2);
        let (mut xi, mut eta, mut zeta) = (2.0 * dot_rows(&v, 1, 2), 2.0 * dot_rows(&v, 0, 2), 2.0 * dot_rows(&v, 0, 1));
        if a > b + eps || ((a - b).abs() < eps && xi.abs() > eta.abs() + eps) {
            v.swap(0, 1);
            v[2] = scale(v[2], -1.0);
            continue;
        }
        if b > c + eps || ((b - c).abs() < eps && eta.abs() > zeta.abs() + eps) {
            v.swap(1, 2);
            v[0] = scale(v[0], -1.0);
            continue;
        }
        flip_for(&mut v, sign_of(xi) * sign_of(eta) * sign_of(zeta) == 1);
        xi = 2.0 * dot_rows(&v, 1, 2);
        eta = 2.0 * dot_rows(&v, 0, 2);
        zeta = 2.0 * dot_rows(&v, 0, 1);
        if xi.abs() > b + eps || ((xi - b).abs() < eps && 2.0 * eta < zeta - eps) || ((xi + b).abs() < eps && zeta < -eps) {
            v[2] = sub(v[2], scale(v[1], sign(xi)));
            continue;
        }
        if eta.abs() > a + eps || ((eta - a).abs() < eps && 2.0 * xi < zeta - eps) || ((eta + a).abs() < eps && zeta < -eps) {
            v[2] = sub(v[2], scale(v[0], sign(eta)));
            continue;
        }
        if zeta.abs() > a + eps || ((zeta - a).abs() < eps && 2.0 * xi < eta - eps) || ((zeta + a).abs() < eps && eta < -eps) {
            v[1] = sub(v[1], scale(v[0], sign(zeta)));
            continue;
        }
        if xi + eta + zeta + a + b < -eps || ((xi + eta + zeta + a + b).abs() < eps && 2.0 * (a + eta) + zeta > eps) {
            v[2] = add(v[2], add(v[0], v[1]));
            continue;
        }
        return Ok(v);
    }
    Err("Niggli reduction did not converge".into())
}

fn niggli_reduce(_: &Json, data: &mut PipelineData) -> R {
    require_cell(data, "niggli-reduce")?;
    let s = &mut data.structure;
    if !(s.pbc[0] && s.pbc[1] && s.pbc[2]) {
        return Err("niggli-reduce needs a cell periodic along a, b and c".into());
    }
    let cell = cell_of(s);
    let volume = determinant(&cell).abs();
    let reduced = niggli(cell, 1e-5 * volume.powf(2.0 / 3.0))?;
    s.cell_vectors = reduced;
    wrap_all(s)?;
    let len = |k: usize| norm(reduced[k]);
    let ang = |i: usize, j: usize| (dot(reduced[i], reduced[j]) / (len(i) * len(j))).acos() / DEGREE;
    data.notes.push(format!(
        "a b c = {} {} {} Angstrom; alpha beta gamma = {} {} {}",
        g(len(0)),
        g(len(1)),
        g(len(2)),
        g(ang(1, 2)),
        g(ang(0, 2)),
        g(ang(0, 1))
    ));
    Ok(())
}

fn scale_cell(p: &Json, data: &mut PipelineData) -> R {
    let id = "scale-cell";
    require_cell(data, id)?;
    let mode = text(p, id, "mode")?;
    let value = number(p, id, "value")?;
    if !(value > 0.0) {
        return Err("scale-cell: the value must be positive".into());
    }
    let s = &mut data.structure;
    let volume = determinant(&cell_of(s)).abs();
    let mut mass = 0.0;
    for a in &s.atoms {
        mass += elements::atomic_mass(a.atomic_number);
    }
    if mode == "axes" {
        let f = vector(p, id, "factors")?;
        if !(f[0] > 0.0 && f[1] > 0.0 && f[2] > 0.0) {
            return Err("scale-cell: the factors must be positive".into());
        }
        let inv = inverse(&cell_of(s))?;
        if flag(p, id, "scale_atoms")? {
            let (offset, cell) = (s.cell_offset, cell_of(s));
            for a in s.atoms.iter_mut() {
                let mut fr = row_times(sub(position_of(a), offset), &inv);
                for k in 0..3 {
                    fr[k] *= f[k];
                }
                place(a, add(offset, row_times(fr, &cell)));
            }
        }
        for k in 0..3 {
            s.cell_vectors[k] = scale(s.cell_vectors[k], f[k]);
        }
        let v = determinant(&cell_of(s)).abs();
        let mut note = format!(
            "a, b, c x {}, {}, {} (strains {}%, {}%, {}%); volume {} Angstrom^3",
            g(f[0]),
            g(f[1]),
            g(f[2]),
            g(100.0 * (f[0] - 1.0)),
            g(100.0 * (f[1] - 1.0)),
            g(100.0 * (f[2] - 1.0)),
            g(v)
        );
        if !s.atoms.is_empty() {
            note += &format!(", density {} g/cm3", g(mass * 1.66053906660 / v));
        }
        data.notes.push(note);
        return Ok(());
    }
    let factor = match mode.as_str() {
        "volume" => (value / volume).cbrt(),
        "density" => {
            if s.atoms.is_empty() {
                return Err("scale-cell to a density needs atoms".into());
            }
            (mass * 1.66053906660 / value / volume).cbrt()
        }
        "factor" => value,
        _ => return Err("scale-cell: mode is factor, volume, density or axes".into()),
    };
    if flag(p, id, "scale_atoms")? {
        let offset = s.cell_offset;
        for a in s.atoms.iter_mut() {
            let q = add(offset, scale(sub(position_of(a), offset), factor));
            place(a, q);
        }
    }
    for row in s.cell_vectors.iter_mut() {
        *row = scale(*row, factor);
    }
    let mut note = format!("Lengths x {}; volume {} Angstrom^3", g(factor), g(volume * factor * factor * factor));
    if !s.atoms.is_empty() {
        note += &format!(", density {} g/cm3", g(mass * 1.66053906660 / (volume * factor * factor * factor)));
    }
    data.notes.push(note);
    Ok(())
}

// ---------------------------------------------------------------- molecules
fn components(nb: &Neighbours) -> Vec<i32> {
    let mut label = vec![-1; nb.atoms.len()];
    let mut next = 0;
    for start in 0..nb.atoms.len() {
        if label[start] >= 0 {
            continue;
        }
        let mut queue = VecDeque::from([start]);
        label[start] = next;
        while let Some(atom) = queue.pop_front() {
            for &other in &nb.atoms[atom] {
                if label[other] < 0 {
                    label[other] = next;
                    queue.push_back(other);
                }
            }
        }
        next += 1;
    }
    label
}

fn unwrap_molecules(p: &Json, data: &mut PipelineData) -> R {
    let id = "unwrap-molecules";
    let nb = bond_neighbours(&data.structure, number(p, id, "bond_scale")?)?;
    let s = &mut data.structure;
    let mut done = vec![false; s.atoms.len()];
    let (mut molecules, mut moved) = (0, 0usize);
    for start in 0..s.atoms.len() {
        if done[start] {
            continue;
        }
        molecules += 1;
        done[start] = true;
        let mut queue = VecDeque::from([start]);
        while let Some(atom) = queue.pop_front() {
            for k in 0..nb.atoms[atom].len() {
                let other = nb.atoms[atom][k];
                if done[other] {
                    continue;
                }
                done[other] = true;
                let target = add(position_of(&s.atoms[atom]), nb.vectors[atom][k]);
                if norm(sub(target, position_of(&s.atoms[other]))) > 1e-6 {
                    moved += 1;
                }
                place(&mut s.atoms[other], target);
                queue.push_back(other);
            }
        }
    }
    data.notes.push(format!("{molecules} molecules or fragments; {moved} atoms moved next to their bonded neighbours"));
    Ok(())
}

fn select_fragment(p: &Json, data: &mut PipelineData) -> R {
    let id = "select-fragment";
    let nb = bond_neighbours(&data.structure, number(p, id, "bond_scale")?)?;
    let label = components(&nb);
    let count = label.iter().max().map_or(0, |&m| m + 1) as usize;
    let mut size = vec![0usize; count];
    for &l in &label {
        size[l as usize] += 1;
    }
    let mut chosen = vec![0u8; count];
    if text(p, id, "rule")? == "touching" {
        for i in 0..label.len() {
            if data.selected[i] != 0 {
                chosen[label[i] as usize] = 1;
            }
        }
    } else {
        let (low, high) = (number(p, id, "min_atoms")?, number(p, id, "max_atoms")?);
        for l in 0..count {
            let n = size[l] as f64;
            chosen[l] = (n >= low && (high <= 0.0 || n <= high)) as u8;
        }
    }
    let fragments: i32 = chosen.iter().map(|&c| c as i32).sum();
    let mask: Vec<u8> = label.iter().map(|&l| chosen[l as usize]).collect();
    select(data, &mask, &text(p, id, "mode")?);
    data.notes.push(format!("{fragments} of {count} fragments"));
    Ok(())
}

// ---------------------------------------------------------------- per-atom simulation data
fn magmom_line(s: &af_core::Structure) -> String {
    let mut out = String::new();
    let mut i = 0;
    while i < s.atoms.len() {
        let mut j = i;
        while j < s.atoms.len() && s.atoms[j].magmom == s.atoms[i].magmom {
            j += 1;
        }
        out += &format!("{}{}*{}", if i > 0 { " " } else { "" }, j - i, g(s.atoms[i].magmom));
        i = j;
    }
    out
}

fn element_values(list: &str, what: &str) -> Result<BTreeMap<i32, f64>, String> {
    let list = list.replace(',', " ");
    let mut words = list.split_ascii_whitespace();
    let mut out = BTreeMap::new();
    while let Some(symbol) = words.next() {
        let Some(value) = words.next() else {
            return Err(format!("{what}: give pairs of an element and a value, e.g. Fe 2.5 O 0"));
        };
        let (v, used) = af_core::format::strtod(value);
        if used != value.len() {
            return Err(format!("{what}: '{value}' is not a number"));
        }
        out.insert(af_science::data::atomic_number(&standard_element(symbol)?), v);
    }
    if out.is_empty() {
        return Err(format!("{what}: no element values given"));
    }
    Ok(out)
}

fn layers_along(s: &af_core::Structure, atoms: &[usize], axis: usize, tolerance: f64) -> Vec<i32> {
    let coordinate = |i: usize| position_of(&s.atoms[i])[axis];
    let mut order = atoms.to_vec();
    sort_by(&mut order, |&a, &b| coordinate(a) < coordinate(b));
    let mut layer: BTreeMap<usize, i32> = BTreeMap::new();
    let (mut current, mut start) = (-1, f64::NEG_INFINITY);
    for &i in &order {
        if coordinate(i) > start + tolerance {
            current += 1;
            start = coordinate(i);
        }
        layer.insert(i, current);
    }
    atoms.iter().map(|i| layer[i]).collect()
}

fn axis_index(name: &str) -> usize {
    match name {
        "x" => 0,
        "y" => 1,
        _ => 2,
    }
}

fn set_magnetic_moments(p: &Json, data: &mut PipelineData) -> R {
    let id = "set-magnetic-moments";
    let atoms = targets(data, text(p, id, "target")? == "selected", id)?;
    let mode = text(p, id, "mode")?;
    let value = number(p, id, "value")?;
    let what = match mode.as_str() {
        "value" => {
            for &i in &atoms {
                data.structure.atoms[i].magmom = value;
            }
            format!("set to {} Bohr magnetons", g(value))
        }
        "elements" => {
            let values = element_values(&text(p, id, "moments")?, id)?;
            for &i in &atoms {
                if let Some(&v) = values.get(&data.structure.atoms[i].atomic_number) {
                    data.structure.atoms[i].magmom = v;
                }
            }
            "set per element".to_string()
        }
        "layered" => {
            let axis = text(p, id, "axis")?;
            let layer = layers_along(&data.structure, &atoms, axis_index(&axis), number(p, id, "layer_tolerance")?);
            let period = integer_of(p, id, "period", 1)?;
            for k in 0..atoms.len() {
                data.structure.atoms[atoms[k]].magmom = if (layer[k] / period) % 2 != 0 { -value } else { value };
            }
            let layers = if atoms.is_empty() { 0 } else { layer.iter().max().unwrap() + 1 };
            format!("+-{} alternating every {} layer(s) along {} ({} layers)", g(value), period, axis, layers)
        }
        "checkerboard" => {
            let tolerance = number(p, id, "layer_tolerance")?;
            let s = &data.structure;
            let (lx, ly, lz) = (layers_along(s, &atoms, 0, tolerance), layers_along(s, &atoms, 1, tolerance), layers_along(s, &atoms, 2, tolerance));
            for k in 0..atoms.len() {
                data.structure.atoms[atoms[k]].magmom = if (lx[k] + ly[k] + lz[k]) % 2 != 0 { -value } else { value };
            }
            format!("+-{} in a checkerboard (G-type)", g(value))
        }
        "random" => {
            let mut order = atoms.clone();
            let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
            shuffle(&mut order, &mut rng);
            let half = order.len() / 2;
            for (k, &i) in order.iter().enumerate() {
                data.structure.atoms[i].magmom = if k < half { -value } else { value };
            }
            format!("+-{} at random (half up, half down)", g(value))
        }
        _ => return Err("set-magnetic-moments: mode is value, elements, layered, checkerboard or random".into()),
    };
    let (mut total, mut absolute) = (0.0, 0.0);
    for a in &data.structure.atoms {
        total += a.magmom;
        absolute += a.magmom.abs();
    }
    let line = magmom_line(&data.structure);
    data.notes.push(format!("{} atoms {}; total {}, absolute {}; MAGMOM = {} (in this atom order)", atoms.len(), what, g(total), g(absolute), line));
    Ok(())
}

fn set_charges(p: &Json, data: &mut PipelineData) -> R {
    let id = "set-charges";
    let atoms = targets(data, text(p, id, "target")? == "selected", id)?;
    let mode = text(p, id, "mode")?;
    let s = &mut data.structure;
    match mode.as_str() {
        "value" => {
            let value = number(p, id, "value")?;
            for &i in &atoms {
                s.atoms[i].charge = value;
            }
        }
        "elements" => {
            let values = element_values(&text(p, id, "charges")?, id)?;
            for &i in &atoms {
                if let Some(&v) = values.get(&s.atoms[i].atomic_number) {
                    s.atoms[i].charge = v;
                }
            }
        }
        _ => return Err("set-charges: mode is value or elements".into()),
    }
    let mut total = 0.0;
    for a in &s.atoms {
        total += a.charge;
    }
    let mut note = String::new();
    let neutralize = text(p, id, "neutralize")?;
    if neutralize == "shift" && !atoms.is_empty() && total.abs() > 1e-12 {
        let shift = total / atoms.len() as f64;
        for &i in &atoms {
            s.atoms[i].charge -= shift;
        }
        note += &format!("shifted by {} e to neutralise; ", g(-shift));
    } else if neutralize == "scale" && total.abs() > 1e-12 {
        let (mut positive, mut negative, mut other) = (0.0, 0.0, total);
        for &i in &atoms {
            let c = s.atoms[i].charge;
            if c > 0.0 {
                positive += c;
            } else {
                negative += c;
            }
            other -= c;
        }
        if total > 0.0 && positive > 0.0 {
            let f = (-negative - other) / positive;
            for &i in &atoms {
                if s.atoms[i].charge > 0.0 {
                    s.atoms[i].charge *= f;
                }
            }
            note += &format!("positive charges x {} to neutralise; ", g(f));
        } else if total < 0.0 && negative < 0.0 {
            let f = (-positive - other) / negative;
            for &i in &atoms {
                if s.atoms[i].charge < 0.0 {
                    s.atoms[i].charge *= f;
                }
            }
            note += &format!("negative charges x {} to neutralise; ", g(f));
        }
    } else if neutralize != "none" && neutralize != "shift" && neutralize != "scale" {
        return Err("set-charges: neutralize is none, shift or scale".into());
    }
    total = 0.0;
    for a in &s.atoms {
        total += a.charge;
    }
    note += &format!("{} atoms set; total charge {} e{}", atoms.len(), g(total), if total.abs() > 1e-9 { " (not neutral)" } else { " (neutral)" });
    data.notes.push(note);
    Ok(())
}

fn set_velocities(p: &Json, data: &mut PipelineData) -> R {
    let id = "set-velocities";
    let atoms = targets(data, text(p, id, "target")? == "selected", id)?;
    let temperature = number(p, id, "temperature")?;
    if temperature < 0.0 {
        return Err("set-velocities: the temperature must not be negative".into());
    }
    const KB: f64 = 8.617333262e-5;
    const EV_PER_AMU_A2_FS2: f64 = 103.6427;
    let mut rng = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
    let mut gauss = Normal::new(0.0, 1.0);
    let s = &data.structure;
    let mass: Vec<f64> = atoms.iter().map(|&i| elements::atomic_mass(s.atoms[i].atomic_number)).collect();
    let mut v: Vec<Vec3> = Vec::with_capacity(atoms.len());
    for k in 0..atoms.len() {
        let sigma = (KB * temperature / (mass[k] * EV_PER_AMU_A2_FS2)).sqrt();
        let x = gauss.sample(&mut rng) * sigma;
        let y = gauss.sample(&mut rng) * sigma;
        let z = gauss.sample(&mut rng) * sigma;
        v.push([x, y, z]);
    }
    let zero_momentum = flag(p, id, "zero_momentum")? && atoms.len() > 1;
    if zero_momentum {
        let (mut momentum, mut total) = ([0.0; 3], 0.0);
        for k in 0..atoms.len() {
            momentum = add(momentum, scale(v[k], mass[k]));
            total += mass[k];
        }
        let drift = scale(momentum, 1.0 / total);
        for vk in v.iter_mut() {
            *vk = sub(*vk, drift);
        }
    }
    let zero_rotation = flag(p, id, "zero_rotation")? && atoms.len() > 2;
    if zero_rotation {
        let (mut com, mut total) = ([0.0; 3], 0.0);
        for k in 0..atoms.len() {
            com = add(com, scale(position_of(&s.atoms[atoms[k]]), mass[k]));
            total += mass[k];
        }
        com = scale(com, 1.0 / total);
        let mut angular = [0.0; 3];
        let mut inertia = [[0.0; 3]; 3];
        for k in 0..atoms.len() {
            let r = sub(position_of(&s.atoms[atoms[k]]), com);
            angular = add(angular, scale(cross(r, v[k]), mass[k]));
            for i in 0..3 {
                for j in 0..3 {
                    inertia[i][j] += mass[k] * ((if i == j { dot(r, r) } else { 0.0 }) - r[i] * r[j]);
                }
            }
        }
        if determinant(&inertia).abs() > 1e-12 {
            let omega = rotate(&inverse(&inertia)?, angular);
            for k in 0..atoms.len() {
                v[k] = sub(v[k], cross(omega, sub(position_of(&s.atoms[atoms[k]]), com)));
            }
        }
    }
    let mut kinetic = 0.0;
    for k in 0..atoms.len() {
        kinetic += 0.5 * mass[k] * dot(v[k], v[k]) * EV_PER_AMU_A2_FS2;
    }
    let dof = 3.0 * atoms.len() as f64 - if zero_momentum { 3.0 } else { 0.0 } - if zero_rotation { 3.0 } else { 0.0 };
    if temperature > 0.0 && kinetic > 0.0 && dof > 0.0 {
        let factor = (temperature / (2.0 * kinetic / (dof * KB))).sqrt();
        for vk in v.iter_mut() {
            *vk = scale(*vk, factor);
        }
    }
    for k in 0..atoms.len() {
        data.structure.atoms[atoms[k]].velocity = v[k];
    }
    data.notes.push(format!(
        "{} velocities at {} K (Angstrom/fs{}{})",
        atoms.len(),
        g(temperature),
        if zero_momentum { ", no net momentum" } else { "" },
        if zero_rotation { ", no rotation" } else { "" }
    ));
    Ok(())
}

// ---------------------------------------------------------------- stacking
fn stack_layers(p: &Json, data: &mut PipelineData) -> R {
    let id = "stack-layers";
    require_cell(data, id)?;
    let file = text(p, id, "file")?;
    if file.is_empty() {
        return Err("Choose the structure file to stack".into());
    }
    let frames = read_frames(Path::new(&file), None)?;
    let top = &frames.first().ok_or("The structure to stack has no atoms")?.structure;
    if top.atoms.is_empty() {
        return Err("The structure to stack has no atoms".into());
    }
    let gap = number(p, id, "gap")?;
    let vacuum = number(p, id, "vacuum")?;
    let strain = flag(p, id, "strain")? && top.has_unit_cell;
    let matching = text(p, id, "match")?;
    if matching != "layer" && matching != "average" && matching != "substrate" {
        return Err("stack-layers: match is layer, average or substrate".into());
    }
    let copies = integer_of(p, id, "copies", 1)?;
    let shift = vector(p, id, "shift")?;
    let twist = number(p, id, "twist")?;
    let select_new = flag(p, id, "select")?;
    let s = &mut data.structure;
    let mut cell = cell_of(s);
    let mut normal = cross(cell[0], cell[1]);
    normal = scale(normal, 1.0 / norm(normal));
    if dot(normal, cell[2]) < 0.0 {
        normal = scale(normal, -1.0);
    }
    let mut placed: Vec<Vec3> = Vec::new();
    if strain {
        let other_cell = top.cell_vectors;
        let la = norm(other_cell[0]) / norm(cell[0]);
        let lb = norm(other_cell[1]) / norm(cell[1]);
        let fa = match matching.as_str() {
            "layer" => 1.0,
            "substrate" => la,
            _ => 0.5 * (1.0 + la),
        };
        let fb = match matching.as_str() {
            "layer" => 1.0,
            "substrate" => lb,
            _ => 0.5 * (1.0 + lb),
        };
        if fa != 1.0 || fb != 1.0 {
            let inv = inverse(&cell)?;
            let offset = s.cell_offset;
            for a in s.atoms.iter_mut() {
                let f = row_times(sub(position_of(a), offset), &inv);
                place(a, add(offset, add(add(scale(cell[0], f[0] * fa), scale(cell[1], f[1] * fb)), scale(cell[2], f[2]))));
            }
            s.cell_vectors[0] = scale(cell[0], fa);
            s.cell_vectors[1] = scale(cell[1], fb);
            cell = cell_of(s);
        }
        let mut other_normal = cross(other_cell[0], other_cell[1]);
        other_normal = scale(other_normal, 1.0 / norm(other_normal));
        let inv = inverse(&other_cell)?;
        for a in &top.atoms {
            let f = row_times(position_of(a), &inv);
            let in_plane = add(scale(cell[0], f[0]), scale(cell[1], f[1]));
            placed.push(add(in_plane, scale(normal, dot(position_of(a), other_normal))));
        }
        let g4 = |x: f64| format_g(x, 4);
        data.notes.push(format!(
            "Strains: stacked layer {}% along a and {}% along b; this structure {}% and {}%",
            g4(100.0 * (fa / la - 1.0)),
            g4(100.0 * (fb / lb - 1.0)),
            g4(100.0 * (fa - 1.0)),
            g4(100.0 * (fb - 1.0))
        ));
    } else {
        placed = top.atoms.iter().map(position_of).collect();
    }
    if twist != 0.0 {
        let mut c = [0.0; 3];
        for q in &placed {
            c = add(c, *q);
        }
        c = scale(c, 1.0 / placed.len() as f64);
        let r = rotation_about(normal, twist * DEGREE);
        for q in placed.iter_mut() {
            *q = add(c, rotate(&r, sub(*q, c)));
        }
        data.notes.push("The stacked layer is turned about the normal; unless the angle is commensurate with the cell, the stack is not periodic".into());
    }
    let s = &mut data.structure;
    let before = s.atoms.len();
    if data.selected.len() != before {
        data.selected = vec![0; before];
    }
    data.selected.iter_mut().for_each(|c| *c = 0);
    let mut low = f64::INFINITY;
    for q in &placed {
        let d = dot(*q, normal);
        if d < low {
            low = d;
        }
    }
    for copy in 1..=copies {
        let mut high = f64::NEG_INFINITY;
        for a in &s.atoms {
            let d = dot(position_of(a), normal);
            if high < d {
                high = d;
            }
        }
        if s.atoms.is_empty() {
            high = 0.0;
        }
        let lateral = add(scale(cell[0], shift[0] * copy as f64), scale(cell[1], shift[1] * copy as f64));
        let lift = add(scale(normal, high + gap - low), lateral);
        for (k, t) in top.atoms.iter().enumerate() {
            append_atom(s, &mut data.selected, &t.symbol, add(placed[k], lift), select_new);
            let last = s.atoms.last_mut().unwrap();
            last.magmom = t.magmom;
            last.charge = t.charge;
            last.fixed = t.fixed;
        }
    }
    if shift[0] != 0.0 || shift[1] != 0.0 {
        let inv = inverse(&cell)?;
        let offset = s.cell_offset;
        for a in s.atoms[before..].iter_mut() {
            let mut f = row_times(sub(position_of(a), offset), &inv);
            f[0] -= f[0].floor();
            f[1] -= f[1].floor();
            place(a, add(row_times(f, &cell), offset));
        }
    }
    let (mut new_high, mut new_low) = (f64::NEG_INFINITY, f64::INFINITY);
    for a in &s.atoms {
        let d = dot(position_of(a), normal);
        if new_high < d {
            new_high = d;
        }
        if d < new_low {
            new_low = d;
        }
    }
    let needed = new_high - new_low + vacuum;
    let along = dot(cell[2], normal);
    if needed > along {
        s.cell_vectors[2] = scale(cell[2], needed / along);
    }
    let c_len = norm(s.cell_vectors[2]);
    data.notes.push(format!("{} x {} atoms stacked {} Angstrom apart; c = {} Angstrom", copies, top.atoms.len(), g(gap), g(c_len)));
    Ok(())
}

fn t(id: &str, title: &str, category: &str, help: &str, parameters: Vec<ModifierParameter>, apply: fn(&Json, &mut PipelineData) -> R) -> ModifierType {
    ModifierType::new(id, title, category, help, parameters, apply)
}

fn with(mut list: Vec<ModifierParameter>, extra: &[&ModifierParameter]) -> Vec<ModifierParameter> {
    list.extend(extra.iter().map(|&p| p.clone()));
    list
}

pub fn types() -> Vec<ModifierType> {
    let target = ModifierParameter::new("target", "Atoms", "choice", "all", "all|selected");
    let centre = ModifierParameter::new("center", "About", "choice", "centroid", "centroid|cell|origin|point");
    let point = ModifierParameter::new("point", "Point (Angstrom; for center = point)", "vector", "0 0 0", "");
    let wrap = ModifierParameter::new("wrap", "Wrap into the cell", "bool", "false", "");
    vec![
        t("mirror", "Mirror", "Modification",
          "Reflects atoms through the plane with the given normal through the atoms' centroid, the cell centre, the origin or a point. Reflecting all atoms reflects the cell too (made right-handed again), so a crystal stays periodic.",
          with(params(&[("normal", "Plane normal", "vector", "0 0 1", "")]), &[&centre, &point, &target, &wrap]), mirror),
        t("invert", "Invert through a point", "Modification", "Inversion x -> 2p - x through the centroid, the cell centre, the origin or a point.",
          with(Vec::new(), &[&centre, &point, &target, &wrap]), invert),
        t("rotate", "Rotate", "Modification",
          "Turns atoms by an angle (degrees, right-handed) about an axis, or by Euler angles (intrinsic z-y-z or z-x-z, or extrinsic x-y-z), about the centroid, the cell centre, the origin or a point. All atoms turn with the cell by default, about the cell origin, so a crystal keeps its periodicity and fractional coordinates.",
          with(
              with(
                  params(&[("axis", "Axis", "vector", "0 0 1", ""), ("angle", "Angle (degrees)", "float", "90", ""), ("mode", "Rotation", "choice", "axis", "axis|euler"),
                           ("euler", "Euler angles (degrees)", "vector", "0 0 0", ""), ("convention", "Euler convention", "choice", "zyz", "zyz|zxz|xyz")]),
                  &[&centre, &point, &target],
              ),
              &[&ModifierParameter::new("rotate_cell", "Rotate the cell too (all atoms)", "bool", "true", "")],
          ),
          rotate_step),
        t("align-direction", "Align a direction with an axis", "Cell",
          "Turns the whole structure and its cell so that a direction (lattice [u v w] or Cartesian) lies along x, y or z, and optionally a second direction in the plane of that axis and a second one (as Atomsk -orient).",
          params(&[("direction", "Direction", "vector", "1 1 1", ""), ("lattice", "Lattice direction [u v w]", "bool", "true", ""),
                   ("axis", "Along", "choice", "z", "x|y|z"), ("second_direction", "Second direction", "vector", "1 -1 0", ""),
                   ("second_axis", "Second direction towards", "choice", "none", "none|x|y|z")]),
          align_direction),
        t("sort-atoms", "Sort atoms", "Modification",
          "Reorders the atoms by element (atomic number, or a given species order such as O Ti for POSCAR blocks), a coordinate, distance from the centroid, magnetic moment, charge, mass, selection or the per-atom property; ties (within the tolerance, e.g. one layer) are broken by a second key.",
          params(&[("key", "Sort by", "choice", "element", "element|x|y|z|distance|magmom|charge|mass|selected|property"),
                   ("then", "Then by", "choice", "none", "none|element|x|y|z|distance|magmom|charge|mass|selected|property"),
                   ("element_order", "Species order (elements first; empty = atomic number)", "string", "", ""),
                   ("tolerance", "Tie tolerance", "float", "0", ""), ("descending", "Descending", "bool", "false", "")]),
          sort_atoms),
        t("saturate-hydrogen", "Saturate with hydrogen", "Modification",
          "Adds hydrogen (or another terminator: F, Cl, OH-like O, ...) to under-coordinated atoms (bonds from covalent radii): tetrahedral (sp3) directions, or in-plane for planar carbon, boron and nitrogen with two bonds near 120 degrees; X-H lengths C 1.09, Si 1.48, N 1.01, O 0.97 Angstrom, other terminators at the sum of covalent radii unless a length is given. Only selected atoms, or at most a number per atom, on request.",
          params(&[("elements", "Elements to saturate", "elements", "C Si Ge N O B P S", ""), ("hybridization", "Geometry", "choice", "auto", "auto|sp3"),
                   ("terminator", "Terminating element", "element", "H", ""),
                   ("bond_length", "Bond length (Angstrom; 0 = table or covalent radii)", "float", "0", ""),
                   ("max_per_atom", "At most per atom", "int", "4", ""), ("target", "Atoms to saturate", "choice", "all", "all|selected"),
                   ("bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""), ("select", "Select the added atoms", "bool", "true", "")]),
          saturate_hydrogen),
        t("swap-elements", "Swap elements", "Modification", "Exchanges two elements everywhere (A becomes B and B becomes A).",
          params(&[("first", "Element A", "element", "Ni", ""), ("second", "Element B", "element", "Al", "")]), swap_elements),
        t("set-geometry", "Set distance, angle or dihedral", "Modification",
          "Sets the distance (2 atoms), angle (3) or dihedral (4; atom indices from 0) by moving the last atom, or the third and everything beyond the middle bond of a dihedral, with its bonded fragment (as in Avogadro); or both sides by half each.",
          params(&[("atoms", "Atom indices", "list", "0 1", ""), ("value", "Value (Angstrom or degrees)", "float", "1.5", ""),
                   ("move", "Move", "choice", "fragment", "fragment|atom|both"), ("bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", "")]),
          set_geometry),
        t("measure-geometry", "Measure distance, angle or dihedral", "Analysis",
          "Reports the distance (2 atoms), angle (3) or dihedral angle (4) between atoms given by index (from 0), with minimum images in periodic cells.",
          params(&[("atoms", "Atom indices", "list", "0 1 2 3", "")]), measure_geometry),
        t("duplicate-selected", "Duplicate selection", "Modification", "Adds copies of the selected atoms shifted by k x offset (k = 1 ... copies).",
          params(&[("offset", "Offset (Angstrom)", "vector", "0 0 3", ""), ("copies", "Copies", "int", "1", ""), ("select", "Select the copies", "bool", "true", "")]),
          duplicate_selected),
        t("niggli-reduce", "Niggli-reduced cell", "Cell",
          "Replaces the cell by its Niggli-reduced form (shortest, most orthogonal basis of the same lattice; Krivy-Gruber) and wraps the atoms into it.",
          Vec::new(), niggli_reduce),
        t("scale-cell", "Scale cell", "Cell",
          "Scales the cell and atom positions uniformly by a factor, or to a volume (Angstrom^3) or a density (g/cm3); or a, b and c by separate factors (uniaxial and biaxial strain). Without scaling the atoms, their Cartesian positions stay.",
          params(&[("mode", "Scale to", "choice", "factor", "factor|volume|density|axes"), ("value", "Value", "float", "1.01", ""),
                   ("factors", "Factors for a, b, c (axes mode)", "vector", "1 1 1", ""), ("scale_atoms", "Scale the atom positions", "bool", "true", "")]),
          scale_cell),
        t("unwrap-molecules", "Make molecules whole", "Modification",
          "Moves atoms across periodic boundaries so that bonded atoms sit next to each other (molecules split by the cell become whole).",
          params(&[("bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", "")]), unwrap_molecules),
        t("set-magnetic-moments", "Set magnetic moments", "Modification",
          "Sets the initial magnetic moment (Bohr magnetons) of the selected (or all) atoms, for spin-polarised calculations: written to extended XYZ (initial_magmoms), CASTEP (SPIN=) and Quantum ESPRESSO (species per moment, nspin = 2); the notes give VASP's MAGMOM line. Modes: one value; a value per element (Fe 2.5 O 0); layered antiferromagnets (sign alternating every `period` layers along an axis: A-type); checkerboard (G-type); random signs with zero net moment. Other orders: select a sublattice and set it.",
          params(&[("value", "Moment (Bohr magnetons)", "float", "2", ""), ("mode", "Moments", "choice", "value", "value|elements|layered|checkerboard|random"),
                   ("moments", "Per element (e.g. Fe 2.5 O 0)", "string", "", ""), ("axis", "Layer axis", "choice", "z", "x|y|z"),
                   ("period", "Layers per sign", "int", "1", ""), ("layer_tolerance", "Layer tolerance (Angstrom)", "float", "0.5", ""),
                   ("seed", "Seed (random)", "int", "1", ""), ("target", "Atoms", "choice", "all", "all|selected")]),
          set_magnetic_moments),
        t("set-charges", "Set charges", "Modification",
          "Sets the charge (e) of the selected (or all) atoms, e.g. formal ionic charges for a LAMMPS charge or full model; written to LAMMPS data (atom style charge) and extended XYZ (initial_charges); the notes give the total charge. One value, or a value per element; neutralising shifts the target atoms' charges equally, or scales the side in excess.",
          params(&[("value", "Charge (e)", "float", "1", ""), ("mode", "Charges", "choice", "value", "value|elements"),
                   ("charges", "Per element (e.g. Ti 2.4 O -1.2)", "string", "", ""),
                   ("neutralize", "Neutralise", "choice", "none", "none|shift|scale"), ("target", "Atoms", "choice", "all", "all|selected")]),
          set_charges),
        t("set-velocities", "Set velocities (temperature)", "Modification",
          "Gives the selected (or all) atoms Maxwell-Boltzmann velocities at a temperature (Angstrom/fs), without net momentum and scaled to that temperature exactly, reproducibly for a seed; written to LAMMPS data (Velocities, Angstrom/ps) and extended XYZ (velo). A temperature of 0 sets them to zero.",
          params(&[("temperature", "Temperature (K)", "float", "300", ""), ("seed", "Seed", "int", "1", ""),
                   ("zero_momentum", "Remove net momentum", "bool", "true", ""), ("zero_rotation", "Remove rotation (clusters, molecules)", "bool", "false", ""),
                   ("target", "Atoms", "choice", "all", "all|selected")]),
          set_velocities),
        t("stack-layers", "Stack a layer", "Modification",
          "Places the atoms of another structure file on top of this one along c, `gap` above the top atom (van der Waals heterostructures, multilayers), copies times. With strain, the two share one in-plane lattice: this cell's (the layer strains), the layer's (this structure strains) or the average (both strain; the strains are reported). Each copy can be shifted laterally (fractional u v, for AB stacking) and turned about the normal; c grows to hold the stack and the vacuum.",
          params(&[("file", "Structure file to stack", "string", "", ""), ("gap", "Gap (Angstrom)", "float", "3.35", ""),
                   ("copies", "Copies", "int", "1", ""), ("shift", "Lateral shift per copy (fractional u v)", "vector", "0 0 0", ""),
                   ("twist", "Turn about the normal (degrees)", "float", "0", ""),
                   ("strain", "Strain to a common lattice", "bool", "true", ""), ("match", "Common lattice", "choice", "layer", "layer|average|substrate"),
                   ("vacuum", "Vacuum (Angstrom)", "float", "15", ""),
                   ("select", "Select the stacked atoms", "bool", "true", "")]),
          stack_layers),
        t("select-fragment", "Select fragments", "Selection",
          "Selects whole bonded fragments (molecules, clusters): those touching the current selection, or those with a number of atoms in a range.",
          with(
              params(&[("rule", "Fragments", "choice", "touching", "touching|size"), ("min_atoms", "Minimum atoms (size rule)", "int", "1", ""),
                       ("max_atoms", "Maximum atoms (size rule; 0 = any)", "int", "0", ""), ("bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", "")]),
              &[&mode_parameter()],
          ),
          select_fragment),
    ]
}
