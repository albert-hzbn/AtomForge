//! The basic modifier library (the C++ pipeline/Modifiers.cpp): selection,
//! modification and cell steps. Each entry describes its parameters (shown by
//! the desktop panel and accepted by the text syntax) and applies itself.

use af_core::cxx_random::{shuffle, Mt19937_64, Normal};
use af_core::elements;
use af_core::json::Json;
use af_core::AtomSite;
use af_science::core::{add, determinant, dot, identity, inverse, multiply, neighbor_list, norm, row_times, scale, sub, Mat3, Vec3};
use af_science::data::{periodic_axes, read_frames};

use crate::expression::Expression;
use crate::pipeline::{params, ModifierParameter, ModifierType, PipelineData};
use crate::support::*;

type R = Result<(), String>;

/// static_cast<std::uint64_t>(double) as x86-64 GCC compiles it (negative
/// values wrap through the signed conversion).
pub fn cxx_u64(x: f64) -> u64 {
    if x < 9_223_372_036_854_775_808.0 {
        x as i64 as u64
    } else {
        ((x - 9_223_372_036_854_775_808.0) as i64 as u64) ^ (1 << 63)
    }
}

/// std::to_string(double): "%f".
pub fn to_string_f(x: f64) -> String {
    format!("{x:.6}")
}

fn colour_atom(symbol: &str, position: Vec3) -> AtomSite {
    let z = elements::atomic_number(symbol);
    let mut atom = AtomSite { symbol: symbol.into(), atomic_number: z, ..Default::default() };
    place(&mut atom, position);
    [atom.r, atom.g, atom.b] = elements::default_colour(z);
    atom
}

/// New atoms extend the per-atom metadata (no grain colour, region -1, no property).
fn push_atom(data: &mut PipelineData, atom: AtomSite, selected: bool) {
    let s = &mut data.structure;
    let n = s.atoms.len();
    if s.grain_colors.len() == n {
        s.grain_colors.push([0.6, 0.6, 0.6]);
    }
    if s.grain_region_ids.len() == n {
        s.grain_region_ids.push(-1);
    }
    if s.atom_property.len() == n && n > 0 {
        s.atom_property.push(f64::NAN);
    }
    if s.frozen_property.len() == n && n > 0 {
        s.frozen_property.push(f64::NAN);
    }
    s.atoms.push(atom);
    data.selected.push(selected as u8);
}

// ---------------------------------------------------------------- selection
fn select_element(p: &Json, data: &mut PipelineData) -> R {
    let list = text(p, "select-element", "elements")?.replace(',', " ");
    let elements = list.split_ascii_whitespace().map(standard_element).collect::<Result<Vec<_>, _>>()?;
    if elements.is_empty() {
        return Err("Name at least one element".into());
    }
    let mask: Vec<u8> = data.structure.atoms.iter().map(|a| elements.contains(&a.symbol) as u8).collect();
    select(data, &mask, &text(p, "select-element", "mode")?);
    Ok(())
}

fn select_expression(p: &Json, data: &mut PipelineData) -> R {
    let expression = Expression::new(&text(p, "select-expression", "expression")?)?;
    let values = evaluate_for_atoms(&expression, data)?;
    let mask: Vec<u8> = values.iter().map(|&v| (v != 0.0) as u8).collect();
    select(data, &mask, &text(p, "select-expression", "mode")?);
    Ok(())
}

fn select_slab(p: &Json, data: &mut PipelineData) -> R {
    let axis = text(p, "select-slab", "axis")?;
    let low = number(p, "select-slab", "minimum")?;
    let high = number(p, "select-slab", "maximum")?;
    if !(low <= high) {
        return Err("The slab minimum must not exceed its maximum".into());
    }
    let fractional = matches!(axis.as_str(), "a" | "b" | "c");
    if fractional {
        require_cell(data, "A fractional slab")?;
    }
    let k = match axis.as_str() {
        "x" | "a" => 0,
        "y" | "b" => 1,
        _ => 2,
    };
    let inv = if fractional { inverse(&cell_of(&data.structure))? } else { [[0.0; 3]; 3] };
    let mask: Vec<u8> = data
        .structure
        .atoms
        .iter()
        .map(|a| {
            let x = position_of(a);
            let mut value = if fractional { row_times(x, &inv)[k] } else { x[k] };
            if fractional {
                value -= value.floor();
            }
            (value >= low && value <= high) as u8
        })
        .collect();
    select(data, &mask, &text(p, "select-slab", "mode")?);
    Ok(())
}

fn select_sphere(p: &Json, data: &mut PipelineData) -> R {
    let centre = vector(p, "select-sphere", "center")?;
    let radius = number(p, "select-sphere", "radius")?;
    if !(radius > 0.0) {
        return Err("The sphere radius must be positive".into());
    }
    let mut mask = Vec::with_capacity(data.structure.atoms.len());
    for a in &data.structure.atoms {
        mask.push((norm(minimum_image(&data.structure, sub(position_of(a), centre))?) <= radius) as u8);
    }
    select(data, &mask, &text(p, "select-sphere", "mode")?);
    Ok(())
}

fn select_property(p: &Json, data: &mut PipelineData) -> R {
    let s = &data.structure;
    if s.atom_property.len() != s.atoms.len() || s.atoms.is_empty() {
        return Err("No per-atom property: add compute-property above (or colour atoms by a result first)".into());
    }
    let low = number(p, "select-property", "minimum")?;
    let high = number(p, "select-property", "maximum")?;
    let mask: Vec<u8> = s.atom_property.iter().map(|&v| (v.is_finite() && v >= low && v <= high) as u8).collect();
    select(data, &mask, &text(p, "select-property", "mode")?);
    Ok(())
}

fn select_random(p: &Json, data: &mut PipelineData) -> R {
    let fraction = number(p, "select-random", "fraction")?;
    if !(0.0..=1.0).contains(&fraction) {
        return Err("The fraction must be between 0 and 1".into());
    }
    // Exactly round(fraction * N) atoms, reproducible for a seed.
    let n = data.structure.atoms.len();
    let mut order: Vec<usize> = (0..n).collect();
    let mut random = Mt19937_64::new(cxx_u64(number(p, "select-random", "seed")?));
    shuffle(&mut order, &mut random);
    let mut mask = vec![0u8; n];
    let count = (fraction * n as f64).round() as usize;
    for &i in order.iter().take(count) {
        mask[i] = 1;
    }
    select(data, &mask, &text(p, "select-random", "mode")?);
    Ok(())
}

fn expand_selection(p: &Json, data: &mut PipelineData) -> R {
    let cutoff = number(p, "expand-selection", "cutoff")?;
    let steps = number(p, "expand-selection", "steps")? as i64;
    let nearest = number(p, "expand-selection", "nearest")? as i64;
    if !(cutoff > 0.0) || steps < 1 || nearest < 0 {
        return Err("expand-selection needs a positive cutoff and at least one step".into());
    }
    let s = &data.structure;
    let positions: Vec<Vec3> = s.atoms.iter().map(position_of).collect();
    let pbc = periodic_axes(s);
    let mut neighbours = neighbor_list(&positions, &if s.has_unit_cell { cell_of(s) } else { identity() }, pbc, cutoff)?;
    if nearest > 0 {
        // Only each atom's N nearest neighbours within the cutoff (stable sort).
        neighbours.sort_by(|a, b| {
            if a.i != b.i {
                a.i.cmp(&b.i)
            } else {
                dot(a.vector, a.vector).partial_cmp(&dot(b.vector, b.vector)).unwrap_or(std::cmp::Ordering::Equal)
            }
        });
        let mut kept = Vec::new();
        let mut run = 0usize;
        for k in 0..neighbours.len() {
            run = if k > 0 && neighbours[k].i == neighbours[k - 1].i { run + 1 } else { 0 };
            if (run as i64) < nearest {
                kept.push(neighbours[k]);
            }
        }
        neighbours = kept;
    }
    for _ in 0..steps {
        let mut grown = data.selected.clone();
        for n in &neighbours {
            if data.selected[n.i] != 0 {
                grown[n.j] = 1;
            }
        }
        data.selected = grown;
    }
    let count = data.selected_count();
    data.notes.push(format!("{count} atoms selected after expansion"));
    Ok(())
}

fn invert_selection(_: &Json, data: &mut PipelineData) -> R {
    for s in data.selected.iter_mut() {
        *s = (*s == 0) as u8;
    }
    let count = data.selected_count();
    data.notes.push(format!("{count} atoms selected"));
    Ok(())
}

fn clear_selection(_: &Json, data: &mut PipelineData) -> R {
    data.clear_selection();
    Ok(())
}

// ---------------------------------------------------------------- modification
fn delete_selected(_: &Json, data: &mut PipelineData) -> R {
    let keep: Vec<usize> = (0..data.structure.atoms.len()).filter(|&i| data.selected[i] == 0).collect();
    let removed = data.structure.atoms.len() - keep.len();
    rebuild(data, &keep);
    data.notes.push(format!("{removed} atoms deleted"));
    Ok(())
}

/// Selective dynamics: hold coordinates of the target atoms fixed.
fn fix_atoms(p: &Json, data: &mut PipelineData) -> R {
    let axes = text(p, "fix-atoms", "axes")?;
    let mask = axes.bytes().fold(0u8, |m, c| {
        m | match c {
            b'x' => 1,
            b'y' => 2,
            b'z' => 4,
            _ => 0,
        }
    });
    let indices = targets(data, text(p, "fix-atoms", "target")? == "selected", "fix-atoms")?;
    for &i in &indices {
        data.structure.atoms[i].fixed |= mask;
    }
    data.notes.push(format!("{} atoms fixed along {axes}", indices.len()));
    Ok(())
}

fn free_atoms(p: &Json, data: &mut PipelineData) -> R {
    let indices = targets(data, text(p, "free-atoms", "target")? == "selected", "free-atoms")?;
    for &i in &indices {
        data.structure.atoms[i].fixed = 0;
    }
    data.notes.push(format!("{} atoms free to move", indices.len()));
    Ok(())
}

fn assign_element(p: &Json, data: &mut PipelineData) -> R {
    let symbol = standard_element(&text(p, "assign-element", "element")?)?;
    let z = elements::atomic_number(&symbol);
    let [r, g, b] = elements::default_colour(z);
    let indices = targets(data, text(p, "assign-element", "target")? == "selected", "assign-element")?;
    for &i in &indices {
        let a = &mut data.structure.atoms[i];
        a.symbol = symbol.clone();
        a.atomic_number = z;
        a.r = r;
        a.g = g;
        a.b = b;
    }
    data.notes.push(format!("{} atoms changed to {symbol}", indices.len()));
    Ok(())
}

fn wrap_atoms(s: &mut af_core::Structure) -> R {
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

fn displace(p: &Json, data: &mut PipelineData) -> R {
    let mut shift = vector(p, "displace", "vector")?;
    if flag(p, "displace", "fractional")? {
        require_cell(data, "A fractional displacement")?;
        shift = row_times(shift, &cell_of(&data.structure));
    }
    for i in targets(data, text(p, "displace", "target")? == "selected", "displace")? {
        let a = &mut data.structure.atoms[i];
        place(a, add(position_of(a), shift));
    }
    if flag(p, "displace", "wrap")? && data.structure.has_unit_cell {
        wrap_atoms(&mut data.structure)?;
    }
    Ok(())
}

fn random_displacement(p: &Json, data: &mut PipelineData) -> R {
    let amplitude = number(p, "random-displacement", "amplitude")?;
    if !(amplitude >= 0.0) {
        return Err("The amplitude must not be negative".into());
    }
    let mut random = Mt19937_64::new(cxx_u64(number(p, "random-displacement", "seed")?));
    let mut gaussian = Normal::new(0.0, amplitude);
    for i in targets(data, text(p, "random-displacement", "target")? == "selected", "random-displacement")? {
        let d = [gaussian.sample(&mut random), gaussian.sample(&mut random), gaussian.sample(&mut random)];
        let a = &mut data.structure.atoms[i];
        place(a, add(position_of(a), d));
    }
    Ok(())
}

fn slice(p: &Json, data: &mut PipelineData) -> R {
    let mut normal = vector(p, "slice", "normal")?;
    let length = norm(normal);
    if !(length > 0.0) {
        return Err("The slice normal must not be zero".into());
    }
    normal = scale(normal, 1.0 / length);
    let distance = number(p, "slice", "distance")?;
    let width = number(p, "slice", "width")?;
    let invert = flag(p, "slice", "invert")?;
    // Without a width the plane cuts away the half-space above it; with a width a slab is cut.
    let inside: Vec<u8> = data
        .structure
        .atoms
        .iter()
        .map(|a| {
            let d = dot(position_of(a), normal) - distance;
            let hit = if width > 0.0 { d.abs() <= 0.5 * width } else { d > 0.0 };
            (hit != invert) as u8
        })
        .collect();
    if text(p, "slice", "action")? == "select" {
        select(data, &inside, "replace");
        return Ok(());
    }
    let keep: Vec<usize> = (0..inside.len()).filter(|&i| inside[i] == 0).collect();
    let removed = inside.len() - keep.len();
    rebuild(data, &keep);
    data.notes.push(format!("{removed} atoms removed by the slice"));
    Ok(())
}

fn add_atom(p: &Json, data: &mut PipelineData) -> R {
    let symbol = standard_element(&text(p, "add-atom", "element")?)?;
    let mut position = vector(p, "add-atom", "position")?;
    if flag(p, "add-atom", "fractional")? {
        require_cell(data, "A fractional position")?;
        position = add(row_times(position, &cell_of(&data.structure)), data.structure.cell_offset);
    }
    let atom = colour_atom(&symbol, position);
    let selected = flag(p, "add-atom", "select")?;
    push_atom(data, atom, selected);
    Ok(())
}

fn insert_interstitials(p: &Json, data: &mut PipelineData) -> R {
    use af_algorithms::interstitial_voids::{detect_interstitial_voids, VoidKind, VoidParams};
    let id = "insert-interstitials";
    let symbol = standard_element(&text(p, id, "element")?)?;
    let params = VoidParams {
        grid_resolution: number(p, id, "resolution")? as i32,
        min_clearance: number(p, id, "clearance")? as f32,
        min_separation: number(p, id, "separation")? as f32,
        ..Default::default()
    };
    let detection = detect_interstitial_voids(&data.structure, &params);
    if !detection.success {
        return Err(if detection.message.is_empty() { "No interstitial sites were found".into() } else { detection.message });
    }
    // Candidate sites of the requested kind (the voids of Edit > Add Interstitial Atoms).
    let kind = text(p, id, "kind")?;
    let mut sites: Vec<&af_algorithms::interstitial_voids::VoidRegion> = detection
        .regions
        .iter()
        .filter(|r| match kind.as_str() {
            "any" => true,
            "tetrahedral" => r.kind == VoidKind::Tetrahedral,
            "octahedral" => r.kind == VoidKind::Octahedral,
            "irregular" => r.kind == VoidKind::Irregular,
            _ => false,
        })
        .collect();
    if sites.is_empty() {
        return Err(format!("No {kind} interstitial sites in this structure"));
    }
    let count = number(p, id, "count")? as i64;
    if count < 0 {
        return Err("The count must not be negative".into());
    }
    if text(p, id, "order")? == "random" {
        let mut random = Mt19937_64::new(cxx_u64(number(p, id, "seed")?));
        shuffle(&mut sites, &mut random);
    } else {
        // std::stable_sort, largest clearance first.
        sites.sort_by(|a, b| b.clearance.partial_cmp(&a.clearance).unwrap_or(std::cmp::Ordering::Equal));
    }
    let n = if count == 0 { sites.len() } else { sites.len().min(count as usize) };
    for site in sites.iter().take(n) {
        let position = [site.position.x as f64, site.position.y as f64, site.position.z as f64];
        push_atom(data, colour_atom(&symbol, position), true);
    }
    data.notes.push(format!("{n} of {} {kind} sites filled with {symbol}", sites.len()));
    Ok(())
}

fn merge(p: &Json, data: &mut PipelineData) -> R {
    let file = text(p, "merge", "file")?;
    if file.is_empty() {
        return Err("Choose the structure file to merge".into());
    }
    let path = std::path::Path::new(&file);
    let frames = read_frames(path, None)?;
    let offset = vector(p, "merge", "offset")?;
    let select_new = flag(p, "merge", "select")?;
    let use_cell = flag(p, "merge", "use_cell")?;
    let s = &mut data.structure;
    let n = s.atoms.len();
    let colours = s.grain_colors.len() == n;
    let regions = s.grain_region_ids.len() == n;
    let property = s.atom_property.len() == n && n > 0;
    let frozen = s.frozen_property.len() == n && n > 0;
    let other = &frames[0].structure;
    for atom in &other.atoms {
        let mut atom = atom.clone();
        let moved = add(position_of(&atom), offset);
        place(&mut atom, moved);
        s.atoms.push(atom);
        if colours {
            s.grain_colors.push([0.6, 0.6, 0.6]);
        }
        if regions {
            s.grain_region_ids.push(-1);
        }
        if property {
            s.atom_property.push(f64::NAN);
        }
        if frozen {
            s.frozen_property.push(f64::NAN);
        }
        data.selected.push(select_new as u8);
    }
    if !s.has_unit_cell && other.has_unit_cell && use_cell {
        s.has_unit_cell = true;
        s.cell_vectors = other.cell_vectors;
    }
    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    data.notes.push(format!("{} atoms merged from {name}", other.atoms.len()));
    Ok(())
}

// ---------------------------------------------------------------- cell
fn set_cell(p: &Json, data: &mut PipelineData) -> R {
    let cell: Mat3 = [vector(p, "set-cell", "a")?, vector(p, "set-cell", "b")?, vector(p, "set-cell", "c")?];
    if determinant(&cell).abs() < 1e-9 {
        return Err("The cell vectors are coplanar".into());
    }
    if flag(p, "set-cell", "scale_atoms")? {
        require_cell(data, "Scaling atoms with the cell")?;
        let s = &mut data.structure;
        let inv = inverse(&cell_of(s))?;
        let offset = s.cell_offset;
        for a in s.atoms.iter_mut() {
            place(a, row_times(row_times(sub(position_of(a), offset), &inv), &cell));
        }
        s.cell_offset = [0.0; 3];
    }
    data.structure.cell_vectors = cell;
    data.structure.has_unit_cell = true;
    Ok(())
}

/// Copies of the atoms with their aligned metadata, by (atom, translation).
fn fill(data: &mut PipelineData, copies: &[(usize, Vec3)], cell: Mat3) {
    let s = &data.structure;
    let count = s.atoms.len();
    let colours = s.grain_colors.len() == count;
    let regions = s.grain_region_ids.len() == count;
    let property = s.atom_property.len() == count;
    let frozen = s.frozen_property.len() == count;
    let mut out = s.clone();
    out.atoms.clear();
    out.grain_colors.clear();
    out.grain_region_ids.clear();
    out.atom_property.clear();
    out.frozen_property.clear();
    let mut selected = Vec::with_capacity(copies.len());
    for &(a, position) in copies {
        let mut atom = s.atoms[a].clone();
        place(&mut atom, position);
        out.atoms.push(atom);
        if colours {
            out.grain_colors.push(s.grain_colors[a]);
        }
        if regions {
            out.grain_region_ids.push(s.grain_region_ids[a]);
        }
        if property {
            out.atom_property.push(s.atom_property[a]);
        }
        if frozen {
            out.frozen_property.push(s.frozen_property[a]);
        }
        selected.push(data.selected[a]);
    }
    out.cell_vectors = cell;
    out.dislocation_loop_points.clear();
    out.dislocation_detection_done = false;
    data.structure = out;
    data.selected = selected;
}

fn replicate(p: &Json, data: &mut PipelineData) -> R {
    require_cell(data, "replicate")?;
    let counts = vector(p, "replicate", "counts")?;
    let mut n = [0i64; 3];
    for k in 0..3 {
        n[k] = counts[k].round() as i64;
        if n[k] < 1 || n[k] > 100 {
            return Err("Replication counts must be between 1 and 100".into());
        }
    }
    if (n[0] * n[1] * n[2]) as f64 * data.structure.atoms.len() as f64 > 50e6 {
        return Err("The replicated structure would exceed 50 million atoms".into());
    }
    let s = &data.structure;
    let v = s.cell_vectors;
    let mut copies = Vec::new();
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                let shift = add(add(scale(v[0], i as f64), scale(v[1], j as f64)), scale(v[2], k as f64));
                for (a, atom) in s.atoms.iter().enumerate() {
                    copies.push((a, add(position_of(atom), shift)));
                }
            }
        }
    }
    let cell = [scale(v[0], n[0] as f64), scale(v[1], n[1] as f64), scale(v[2], n[2] as f64)];
    fill(data, &copies, cell);
    Ok(())
}

/// Integer supercell: the new cell vectors are the rows of M as combinations
/// of a, b and c; the atoms fill the new cell.
fn supercell(p: &Json, data: &mut PipelineData) -> R {
    require_cell(data, "supercell")?;
    let matrix = text(p, "supercell", "matrix")?.replace(',', " ");
    let bad = || "The supercell matrix needs nine integers, row by row".to_string();
    // istream >> double: whitespace-separated words, each read as a number prefix.
    let mut words = matrix.split_ascii_whitespace();
    let mut m: Mat3 = [[0.0; 3]; 3];
    let mut integers = [[0i64; 3]; 3];
    let mut leftover = String::new();
    for r in 0..3 {
        for k in 0..3 {
            let word = if leftover.is_empty() { words.next().map(str::to_string) } else { Some(std::mem::take(&mut leftover)) };
            let Some(word) = word else { return Err(bad()) };
            let Some((value, used)) = af_core::format::parse_prefix_f64(&word) else { return Err(bad()) };
            leftover = word[used..].to_string();
            if (value - value.round()).abs() > 1e-9 {
                return Err(bad());
            }
            integers[r][k] = value.round() as i64;
            m[r][k] = integers[r][k] as f64;
        }
    }
    let det = determinant(&m);
    if det.abs() < 0.5 {
        return Err("The supercell matrix is singular".into());
    }
    let s = &data.structure;
    let count = s.atoms.len();
    if det.abs() * count as f64 > 50e6 {
        return Err("The supercell would exceed 50 million atoms".into());
    }
    let cell = cell_of(s);
    let inverse_cell = inverse(&cell)?;
    let inverse_m = inverse(&m)?;
    let origin = s.cell_offset;
    // The old-cell translations that reach the new cell: its corners.
    let (mut low, mut high) = ([0i64; 3], [0i64; 3]);
    for k in 0..3 {
        for r in 0..3 {
            if integers[r][k] < 0 {
                low[k] += integers[r][k];
            } else {
                high[k] += integers[r][k];
            }
        }
    }
    const EPS: f64 = 1e-7;
    let mut copies = Vec::new();
    for (a, atom) in s.atoms.iter().enumerate() {
        let mut f = row_times(sub(position_of(atom), origin), &inverse_cell);
        for x in f.iter_mut() {
            *x -= x.floor();
        }
        for i in low[0]..=high[0] {
            for j in low[1]..=high[1] {
                for k in low[2]..=high[2] {
                    let shifted = [f[0] + i as f64, f[1] + j as f64, f[2] + k as f64];
                    let g = row_times(shifted, &inverse_m);
                    if g.iter().any(|&x| x < -EPS || x >= 1.0 - EPS) {
                        continue;
                    }
                    copies.push((a, add(origin, row_times(shifted, &cell))));
                }
            }
        }
    }
    let new_cell = multiply(&m, &cell);
    let cells = det.abs().round() as i64;
    fill(data, &copies, new_cell);
    data.notes.push(format!("Supercell: {} atoms ({cells} cells)", data.structure.atoms.len()));
    Ok(())
}

fn transform(p: &Json, data: &mut PipelineData) -> R {
    let m: Mat3 = [vector(p, "transform", "row1")?, vector(p, "transform", "row2")?, vector(p, "transform", "row3")?];
    if determinant(&m).abs() < 1e-12 {
        return Err("The transformation matrix is singular".into());
    }
    let selected_only = text(p, "transform", "target")? == "selected";
    // Positions and cell vectors are row vectors: x' = x M.
    for i in targets(data, selected_only, "transform")? {
        let a = &mut data.structure.atoms[i];
        place(a, row_times(position_of(a), &m));
    }
    if flag(p, "transform", "transform_cell")? && !selected_only && data.structure.has_unit_cell {
        for row in data.structure.cell_vectors.iter_mut() {
            *row = row_times(*row, &m);
        }
    }
    Ok(())
}

fn strain(p: &Json, data: &mut PipelineData) -> R {
    let normal = vector(p, "strain", "normal")?;
    let shear = vector(p, "strain", "shear")?;
    // Engineering-Voigt shears (yz, xz, xy) are twice the tensor components.
    let mut f = identity();
    f[0][0] += normal[0];
    f[1][1] += normal[1];
    f[2][2] += normal[2];
    f[1][2] = 0.5 * shear[0];
    f[2][1] = f[1][2];
    f[0][2] = 0.5 * shear[1];
    f[2][0] = f[0][2];
    f[0][1] = 0.5 * shear[2];
    f[1][0] = f[0][1];
    for a in data.structure.atoms.iter_mut() {
        place(a, row_times(position_of(a), &f));
    }
    if data.structure.has_unit_cell {
        for row in data.structure.cell_vectors.iter_mut() {
            *row = row_times(*row, &f);
        }
    }
    Ok(())
}

fn set_pbc(p: &Json, data: &mut PipelineData) -> R {
    require_cell(data, "set-pbc")?;
    let pbc = [flag(p, "set-pbc", "a")?, flag(p, "set-pbc", "b")?, flag(p, "set-pbc", "c")?];
    data.structure.pbc = pbc;
    let axes: Vec<&str> = ["a", "b", "c"].iter().zip(pbc).filter(|(_, on)| *on).map(|(n, _)| *n).collect();
    data.notes.push(if axes.is_empty() { "No periodic cell vector (open boundaries)".into() } else { format!("Periodic along {}", axes.join(", ")) });
    Ok(())
}

fn wrap(_: &Json, data: &mut PipelineData) -> R {
    require_cell(data, "wrap")?;
    wrap_atoms(&mut data.structure)
}

fn center(_: &Json, data: &mut PipelineData) -> R {
    require_cell(data, "center")?;
    if data.structure.atoms.is_empty() {
        return Ok(());
    }
    let cell = cell_of(&data.structure);
    let inv = inverse(&cell)?;
    let (mut low, mut high) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for a in &data.structure.atoms {
        let f = row_times(position_of(a), &inv);
        for k in 0..3 {
            // std::min(low, f) keeps low unless f < low.
            if f[k] < low[k] {
                low[k] = f[k];
            }
            if high[k] < f[k] {
                high[k] = f[k];
            }
        }
    }
    let shift = [0.5 - 0.5 * (low[0] + high[0]), 0.5 - 0.5 * (low[1] + high[1]), 0.5 - 0.5 * (low[2] + high[2])];
    let cartesian = add(row_times(shift, &cell), data.structure.cell_offset);
    for a in data.structure.atoms.iter_mut() {
        place(a, add(position_of(a), cartesian));
    }
    Ok(())
}

fn add_vacuum(p: &Json, data: &mut PipelineData) -> R {
    require_cell(data, "add-vacuum")?;
    let axis = text(p, "add-vacuum", "axis")?;
    let k = match axis.as_str() {
        "a" => 0,
        "b" => 1,
        _ => 2,
    };
    let thickness = number(p, "add-vacuum", "thickness")?;
    if thickness < 0.0 {
        return Err("The vacuum thickness must not be negative".into());
    }
    let centre = flag(p, "add-vacuum", "center")?;
    let row = data.structure.cell_vectors[k];
    let length = norm(row);
    let direction = scale(row, 1.0 / length);
    data.structure.cell_vectors[k] = scale(row, (length + thickness) / length);
    if centre {
        for a in data.structure.atoms.iter_mut() {
            place(a, add(position_of(a), scale(direction, 0.5 * thickness)));
        }
    }
    Ok(())
}

/// An orthogonal cell just enclosing the atoms (plus padding).
fn shrink_wrap(p: &Json, data: &mut PipelineData) -> R {
    if data.structure.atoms.is_empty() {
        return Err("shrink-wrap needs atoms".into());
    }
    let padding = number(p, "shrink-wrap", "padding")?;
    if padding < 0.0 {
        return Err("padding must be nonnegative".into());
    }
    let s = &mut data.structure;
    let mut low = position_of(&s.atoms[0]);
    let mut high = low;
    for a in &s.atoms {
        let x = position_of(a);
        for d in 0..3 {
            if x[d] < low[d] {
                low[d] = x[d];
            }
            if high[d] < x[d] {
                high[d] = x[d];
            }
        }
    }
    for a in s.atoms.iter_mut() {
        let x = [a.x - low[0] + padding, a.y - low[1] + padding, a.z - low[2] + padding];
        place(a, x);
    }
    s.has_unit_cell = true;
    for r in 0..3 {
        for c in 0..3 {
            let extent = high[r] - low[r] + 2.0 * padding;
            // std::max(extent, 1e-3)
            s.cell_vectors[r][c] = if r == c { if extent < 1e-3 { 1e-3 } else { extent } } else { 0.0 };
        }
    }
    s.cell_offset = [0.0; 3];
    let note = format!(
        "Cell {} x {} x {} Angstrom",
        to_string_f(s.cell_vectors[0][0]),
        to_string_f(s.cell_vectors[1][1]),
        to_string_f(s.cell_vectors[2][2])
    );
    data.notes.push(note);
    Ok(())
}

fn t(
    id: &str,
    title: &str,
    category: &str,
    help: &str,
    parameters: Vec<ModifierParameter>,
    apply: fn(&Json, &mut PipelineData) -> R,
) -> ModifierType {
    ModifierType::new(id, title, category, help, parameters, apply)
}

fn with_mode(mut list: Vec<ModifierParameter>) -> Vec<ModifierParameter> {
    list.push(mode_parameter());
    list
}

/// The basic built-in types (before the other families and the category sort).
pub fn types() -> Vec<ModifierType> {
    vec![
        // Selection
        t("select-element", "Select by element", "Selection", "Selects atoms of the listed elements.",
          with_mode(params(&[("elements", "Elements (e.g. Cu O)", "elements", "Cu", "")])), select_element),
        t("select-expression", "Select by expression", "Selection",
          "Selects atoms for which a condition holds, e.g. fz > 0.5 && element == O.\nVariables: x y z (Angstrom), fx fy fz (fractional), index, count, element,\nZ, property (per-atom property), frozen, selected, fixed, magmom, charge, vx vy vz.\nOperators + - * / % (floored remainder) ^, < <= > >= == !=,\n&& || ! (and or not), functions abs sqrt exp log sin cos tan floor ceil min max.",
          with_mode(params(&[("expression", "Condition", "expression", "x > 0", "")])), select_expression),
        t("select-slab", "Select a slab", "Selection", "Selects atoms between two coordinates along x, y or z (Angstrom) or along a, b or c (fractional, wrapped into [0, 1)).",
          with_mode(params(&[("axis", "Axis", "choice", "z", "x|y|z|a|b|c"), ("minimum", "From", "float", "0", ""), ("maximum", "To", "float", "5", "")])), select_slab),
        t("select-sphere", "Select a sphere", "Selection", "Selects atoms within a radius of a point (minimum-image distance in periodic cells).",
          with_mode(params(&[("center", "Centre (Angstrom)", "vector", "0 0 0", ""), ("radius", "Radius (Angstrom)", "float", "3", "")])), select_sphere),
        t("select-property", "Select by property range", "Selection", "Selects atoms whose per-atom property lies in a range (needs compute-property or a coloured result).",
          with_mode(params(&[("minimum", "Minimum", "float", "0", ""), ("maximum", "Maximum", "float", "1", "")])), select_property),
        t("select-random", "Select a random fraction", "Selection", "Selects round(fraction x N) atoms at random, reproducibly for a seed.",
          with_mode(params(&[("fraction", "Fraction", "float", "0.1", ""), ("seed", "Seed", "int", "1", "")])), select_random),
        t("expand-selection", "Expand selection to neighbours", "Selection", "Adds the neighbours (within the cutoff, or the N nearest) of selected atoms, repeated for the number of steps.",
          params(&[("cutoff", "Neighbour cutoff (Angstrom)", "float", "3", ""), ("steps", "Steps", "int", "1", ""),
                   ("nearest", "Only the N nearest neighbours (0 = all within the cutoff)", "int", "0", "")]), expand_selection),
        t("invert-selection", "Invert selection", "Selection", "Selects the unselected atoms and deselects the rest.", vec![], invert_selection),
        t("clear-selection", "Clear selection", "Selection", "Deselects all atoms.", vec![], clear_selection),
        // Modification
        t("delete-selected", "Delete selected atoms", "Modification", "Removes the selected atoms.", vec![], delete_selected),
        t("assign-element", "Assign element", "Modification", "Changes the element of the selected (or all) atoms, e.g. to create substitutions.",
          params(&[("element", "Element", "element", "Cu", ""), ("target", "Atoms", "choice", "selected", "selected|all")]), assign_element),
        t("fix-atoms", "Fix atoms", "Modification",
          "Holds coordinates of the selected (or all) atoms fixed in a relaxation; POSCAR files keep them as selective dynamics (F).",
          params(&[("axes", "Fixed coordinates", "choice", "xyz", "xyz|xy|xz|yz|x|y|z"), ("target", "Atoms", "choice", "selected", "selected|all")]), fix_atoms),
        t("free-atoms", "Free atoms", "Modification", "Lets the selected (or all) atoms move along every coordinate again.",
          params(&[("target", "Atoms", "choice", "selected", "selected|all")]), free_atoms),
        t("displace", "Displace", "Modification",
          "Moves atoms by a fixed vector, in Angstrom or in fractions of the cell vectors (shifting the origin); optionally wraps them into the cell.",
          params(&[("vector", "Displacement (Angstrom, or fractional)", "vector", "0 0 0", ""), ("target", "Atoms", "choice", "all", "all|selected"),
                   ("fractional", "Fractional displacement", "bool", "false", ""), ("wrap", "Wrap into the cell", "bool", "false", "")]), displace),
        t("random-displacement", "Random displacement", "Modification", "Adds Gaussian displacements (standard deviation per component), reproducibly for a seed.",
          params(&[("amplitude", "Standard deviation (Angstrom)", "float", "0.05", ""), ("seed", "Seed", "int", "1", ""), ("target", "Atoms", "choice", "all", "all|selected")]), random_displacement),
        t("add-atom", "Add atom", "Modification", "Adds one atom at a Cartesian (Angstrom) or fractional position.",
          params(&[("element", "Element", "element", "Cu", ""), ("position", "Position", "vector", "0 0 0", ""), ("fractional", "Fractional position", "bool", "false", ""),
                   ("select", "Select the new atom", "bool", "true", "")]), add_atom),
        t("insert-interstitials", "Insert interstitial atoms", "Modification", "Fills interstitial voids (found as in Edit > Add Interstitial Atoms) with atoms: the largest voids first or at random, all of them (count 0) or a number; the new atoms are selected.",
          params(&[("element", "Element", "element", "H", ""), ("kind", "Void kind", "choice", "any", "any|tetrahedral|octahedral|irregular"),
                   ("count", "Atoms to insert (0 = every site)", "int", "1", ""), ("order", "Choose sites", "choice", "largest", "largest|random"),
                   ("seed", "Seed (random order)", "int", "1", ""), ("resolution", "Grid resolution", "int", "14", ""),
                   ("clearance", "Minimum clearance (Angstrom)", "float", "0.6", ""), ("separation", "Minimum separation (Angstrom)", "float", "0.5", "")]), insert_interstitials),
        t("merge", "Merge a structure file", "Modification", "Appends the atoms of another structure file (POSCAR, extended XYZ, JSON or LAMMPS dump), shifted by an offset.",
          params(&[("file", "Structure file", "string", "", ""), ("offset", "Offset (Angstrom)", "vector", "0 0 0", ""), ("select", "Select the merged atoms", "bool", "true", ""),
                   ("use_cell", "Take its cell when this structure has none", "bool", "true", "")]), merge),
        t("slice", "Slice", "Modification", "Cuts with a plane n.r = distance: removes the half-space above it, or with a width the slab |n.r - distance| <= width/2. Invert keeps that region instead; action select only selects it.",
          params(&[("normal", "Plane normal", "vector", "0 0 1", ""), ("distance", "Distance from the origin (Angstrom)", "float", "0", ""),
                   ("width", "Slab width (0 = half-space)", "float", "0", ""), ("invert", "Invert", "bool", "false", ""), ("action", "Action", "choice", "delete", "delete|select")]), slice),
        // Cell
        t("replicate", "Replicate", "Cell", "Repeats the cell along a, b and c.",
          params(&[("counts", "Copies along a, b, c", "ivec3", "2 2 2", "")]), replicate),
        t("supercell", "Transform structure (supercell)", "Cell", "Builds the supercell whose cell vectors are integer combinations of a, b and c (rows of the matrix, as in Edit > Transform Structure); the atoms fill the new cell.",
          params(&[("matrix", "Matrix rows (nine integers)", "string", "2 0 0 0 2 0 0 0 2", "")]), supercell),
        t("transform", "Affine transformation", "Cell", "Applies x' = x M to atoms (row vectors) and, for all atoms, optionally to the cell.",
          params(&[("row1", "Matrix row 1", "vector", "1 0 0", ""), ("row2", "Matrix row 2", "vector", "0 1 0", ""), ("row3", "Matrix row 3", "vector", "0 0 1", ""),
                   ("transform_cell", "Transform the cell too", "bool", "true", ""), ("target", "Atoms", "choice", "all", "all|selected")]), transform),
        t("strain", "Strain", "Cell", "Homogeneous strain of atoms and cell: normal components xx yy zz and engineering shears yz xz xy.",
          params(&[("normal", "Normal strain xx, yy, zz", "vector", "0 0 0", ""), ("shear", "Engineering shear yz, xz, xy", "vector", "0 0 0", "")]), strain),
        t("set-cell", "Set cell vectors", "Cell", "Replaces the lattice vectors; with scale atoms, atoms keep their fractional coordinates.",
          params(&[("a", "Vector a (Angstrom)", "vector", "1 0 0", ""), ("b", "Vector b", "vector", "0 1 0", ""), ("c", "Vector c", "vector", "0 0 1", ""),
                   ("scale_atoms", "Scale atoms with the cell", "bool", "true", "")]), set_cell),
        t("wrap", "Wrap into cell", "Cell", "Maps atoms outside the cell back into it (along the periodic cell vectors).", vec![], wrap),
        t("set-pbc", "Periodic boundaries", "Cell",
          "Sets which cell vectors are periodic (a slab: a and b; a wire: c). Analyses, neighbour searches, wrapping and extended XYZ files follow it.",
          params(&[("a", "Periodic along a", "bool", "true", ""), ("b", "Periodic along b", "bool", "true", ""), ("c", "Periodic along c", "bool", "true", "")]), set_pbc),
        t("shrink-wrap", "Shrink-wrap cell", "Cell", "Replaces the cell by the orthogonal box that just encloses the atoms, plus padding on every side.",
          params(&[("padding", "Padding (Angstrom)", "float", "0", "")]), shrink_wrap),
        t("center", "Center in cell", "Cell", "Shifts the atoms so their extent (in fractional coordinates) is centred in the cell.", vec![], center),
        t("add-vacuum", "Add vacuum", "Cell", "Lengthens one cell vector by a thickness, e.g. to turn a bulk cell into a slab; optionally centres the atoms.",
          params(&[("axis", "Cell vector", "choice", "c", "a|b|c"), ("thickness", "Vacuum (Angstrom)", "float", "10", ""), ("center", "Centre the atoms", "bool", "true", "")]), add_vacuum),
    ]
}
