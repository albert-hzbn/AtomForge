//! The Shape Modeller step (the C++ pipeline/ShapeSteps): carves a crystal
//! into a shape written as a CAD-like script (af_algorithms::shape), or keeps,
//! deletes or selects the atoms of a structure inside such a shape. Regions
//! give the atoms inside them another element.

use std::collections::BTreeMap;

use af_algorithms::shape::Model;
use af_core::format::format_g;
use af_core::glm;
use af_core::json::Json;
use af_core::{AtomSite, Structure};
use af_science::core::{add, scale, Vec3};
use glam::{DMat3, DVec3};

use crate::pipeline::{params, ModifierParameter, ModifierType, PipelineData};
use crate::support::*;

type R = Result<(), String>;

fn to_glm(v: Vec3) -> DVec3 {
    DVec3::new(v[0], v[1], v[2])
}

fn set_element(a: &mut AtomSite, symbol: &str) {
    if a.symbol == symbol {
        return;
    }
    a.symbol = symbol.to_string();
    a.atomic_number = af_science::data::atomic_number(symbol);
    [a.r, a.g, a.b] = af_core::elements::default_colour(a.atomic_number);
}

fn region_symbols(model: &Model) -> Result<BTreeMap<String, String>, String> {
    let mut out = BTreeMap::new();
    for name in model.region_elements() {
        let symbol = if name == "none" || name == "vacancy" { "none".to_string() } else { standard_element(&name)? };
        out.insert(name, symbol);
    }
    Ok(out)
}

fn vmin(a: DVec3, b: DVec3) -> DVec3 {
    let m = |x: f64, y: f64| if y < x { y } else { x };
    DVec3::new(m(a.x, b.x), m(a.y, b.y), m(a.z, b.z))
}

fn vmax(a: DVec3, b: DVec3) -> DVec3 {
    let m = |x: f64, y: f64| if x < y { y } else { x };
    DVec3::new(m(a.x, b.x), m(a.y, b.y), m(a.z, b.z))
}

/// The shape's surface, placed by `shift`, written as OBJ when mesh_file is set.
fn write_mesh(p: &Json, model: &Model, shift: DVec3, data: &mut PipelineData) -> R {
    let id = "build-shape";
    let file = text(p, id, "mesh_file")?;
    if file.is_empty() {
        return Ok(());
    }
    let mut mesh = model.mesh(number(p, id, "mesh_spacing")?)?;
    for v in mesh.vertices.iter_mut() {
        *v += shift;
    }
    af_algorithms::shape::save_obj(&mesh, &file)?;
    data.notes.push(format!("Shape surface written to {file} ({} triangles)", mesh.vertices.len() / 3));
    Ok(())
}

fn build_shape(p: &Json, data: &mut PipelineData) -> R {
    let id = "build-shape";
    let mode = text(p, id, "mode")?;
    if mode != "carve" && mode != "keep" && mode != "delete" && mode != "select" {
        return Err("build-shape: mode is carve, keep, delete or select".into());
    }
    let has_cell = data.structure.has_unit_cell;
    let mut lattice = DMat3::IDENTITY;
    if has_cell {
        let cell = cell_of(&data.structure);
        lattice = DMat3::from_cols(to_glm(cell[0]), to_glm(cell[1]), to_glm(cell[2]));
    }
    let model = Model::parse(&text(p, id, "script")?, if has_cell { Some(&lattice) } else { None }).map_err(|e| format!("Shape script, {e}"))?;
    let bounds = model.bounds();
    let symbols = region_symbols(&model)?;
    let mut anchor = text(p, id, "anchor")?;
    if anchor == "auto" {
        anchor = (if mode == "carve" {
            "origin"
        } else if has_cell {
            "cell-centre"
        } else {
            "centroid"
        })
        .into();
    }
    let input = &data.structure;
    let mut origin: Vec3 = [0.0; 3];
    if anchor == "cell-centre" {
        require_cell(data, "The cell centre")?;
        let cell = cell_of(input);
        origin = add(input.cell_offset, scale(add(add(cell[0], cell[1]), cell[2]), 0.5));
    } else if anchor == "centroid" {
        for a in &input.atoms {
            origin = add(origin, position_of(a));
        }
        if !input.atoms.is_empty() {
            origin = scale(origin, 1.0 / input.atoms.len() as f64);
        }
    } else if anchor != "origin" {
        return Err("build-shape: anchor is auto, origin, cell-centre or centroid".into());
    }
    origin = add(origin, vector(p, id, "position")?);
    let shape_origin = to_glm(origin);
    let mut region_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut removed_by_region = 0usize;
    let mut note;

    if mode == "carve" {
        require_cell(data, "Carving a shape (open a crystal first, e.g. Build > Bulk Crystal)")?;
        let input = &data.structure;
        if input.atoms.is_empty() {
            return Err("build-shape: the crystal has no atoms".into());
        }
        if !bounds.finite {
            return Err("build-shape: the shape is unbounded; intersect it with a box or other finite shape".into());
        }
        let inverse = glm::inverse_d(&lattice);
        let cell_origin = to_glm(input.cell_offset);
        let (mut low, mut high) = (DVec3::splat(1e300), DVec3::splat(-1e300));
        for c in 0..8 {
            let corner = DVec3::new(
                if c & 1 != 0 { bounds.high.x } else { bounds.low.x },
                if c & 2 != 0 { bounds.high.y } else { bounds.low.y },
                if c & 4 != 0 { bounds.high.z } else { bounds.low.z },
            );
            let f = glm::mul_d(&inverse, corner + shape_origin - cell_origin);
            low = vmin(low, f);
            high = vmax(high, f);
        }
        let mut range = [[0i64; 2]; 3];
        let mut images = 1.0;
        for k in 0..3 {
            range[k][0] = low[k].floor() as i64 - 1;
            range[k][1] = high[k].ceil() as i64 + 1;
            images *= (range[k][1] - range[k][0] + 1) as f64;
        }
        let max_atoms = number(p, id, "max_atoms")?;
        let cell_volume = glm::determinant_d(&lattice).abs();
        let span = ((bounds.high.x - bounds.low.x) * (bounds.high.y - bounds.low.y) * (bounds.high.z - bounds.low.z)).cbrt() / 80.0;
        let spacing = if 0.5 < span { span } else { 0.5 };
        let expected = model.volume(spacing) / cell_volume * input.atoms.len() as f64;
        if expected > max_atoms {
            return Err(format!("build-shape: about {} atoms, more than max_atoms ({})", expected as i64, max_atoms as i64));
        }
        if images * input.atoms.len() as f64 > 4e9 {
            return Err("build-shape: the shape is far too large for this cell".into());
        }
        let mut out = Structure::default();
        for i in range[0][0]..=range[0][1] {
            for j in range[1][0]..=range[1][1] {
                for k in range[2][0]..=range[2][1] {
                    let shift = glm::mul_d(&lattice, DVec3::new(i as f64, j as f64, k as f64));
                    for a in &input.atoms {
                        let position = DVec3::new(a.x, a.y, a.z) + shift;
                        let q = position - shape_origin;
                        if !model.contains(q) {
                            continue;
                        }
                        let mut atom = a.clone();
                        place(&mut atom, position.to_array());
                        if let Some(region) = model.material(q) {
                            let symbol = &symbols[region];
                            if symbol == "none" {
                                removed_by_region += 1;
                                continue;
                            }
                            set_element(&mut atom, symbol);
                            *region_counts.entry(symbol.clone()).or_insert(0) += 1;
                        }
                        out.atoms.push(atom);
                    }
                }
            }
        }
        if out.atoms.is_empty() {
            return Err("build-shape: no atom of the crystal lies inside the shape".into());
        }
        let vacuum = number(p, id, "vacuum")?;
        let (mut lo, mut hi) = (DVec3::splat(1e300), DVec3::splat(-1e300));
        for a in &out.atoms {
            lo = vmin(lo, DVec3::new(a.x, a.y, a.z));
            hi = vmax(hi, DVec3::new(a.x, a.y, a.z));
        }
        let mv = DVec3::splat(vacuum / 2.0) - lo;
        for a in out.atoms.iter_mut() {
            a.x += mv.x;
            a.y += mv.y;
            a.z += mv.z;
        }
        let size = hi - lo + DVec3::splat(vacuum);
        out.has_unit_cell = true;
        out.cell_vectors = [[size.x, 0.0, 0.0], [0.0, size.y, 0.0], [0.0, 0.0, size.z]];
        out.cell_offset = [0.0; 3];
        out.pbc = [false; 3];
        data.structure = out;
        data.clear_selection();
        note = format!("{} atoms carved from the crystal", data.structure.atoms.len());
        write_mesh(p, &model, shape_origin + mv, data)?;
    } else {
        let input = &data.structure;
        let mut inside = vec![0u8; input.atoms.len()];
        let mut count = 0usize;
        for (i, a) in input.atoms.iter().enumerate() {
            if model.contains(to_glm(position_of(a)) - shape_origin) {
                inside[i] = 1;
                count += 1;
            }
        }
        if mode == "select" {
            select(data, &inside, &text(p, id, "combine")?);
            note = format!("{count} atoms inside the shape selected");
        } else {
            let keeping = mode == "keep";
            let mut keep = Vec::new();
            for i in 0..data.structure.atoms.len() {
                if keeping != (inside[i] != 0) {
                    continue;
                }
                if keeping {
                    let q = to_glm(position_of(&data.structure.atoms[i])) - shape_origin;
                    if let Some(region) = model.material(q) {
                        let symbol = &symbols[region];
                        if symbol == "none" {
                            removed_by_region += 1;
                            continue;
                        }
                        set_element(&mut data.structure.atoms[i], symbol);
                        *region_counts.entry(symbol.clone()).or_insert(0) += 1;
                    }
                }
                keep.push(i);
            }
            rebuild(data, &keep);
            note = if keeping { format!("{} atoms inside the shape kept", keep.len()) } else { format!("{count} atoms inside the shape deleted") };
        }
        write_mesh(p, &model, shape_origin, data)?;
    }
    let primitives = model.primitives();
    note += &format!(" ({} primitive{}", primitives, if primitives == 1 { "" } else { "s" });
    if bounds.finite {
        let g5 = |x: f64| format_g(x, 5);
        note += &format!(
            ", {} x {} x {} Angstrom",
            g5(bounds.high.x - bounds.low.x),
            g5(bounds.high.y - bounds.low.y),
            g5(bounds.high.z - bounds.low.z)
        );
    }
    note += ")";
    for (symbol, n) in &region_counts {
        note += &format!("; {n} {symbol} in regions");
    }
    if removed_by_region > 0 {
        note += &format!("; {removed_by_region} removed by \"none\" regions");
    }
    data.notes.insert(0, note);
    Ok(())
}

pub fn types() -> Vec<ModifierType> {
    let mut list: Vec<ModifierParameter> = params(&[
        ("script", "Shape script", "string", "sphere(r = 15);", ""),
        ("mode", "Mode", "choice", "carve", "carve|keep|delete|select"),
        ("anchor", "Shape origin at", "choice", "auto", "auto|origin|cell-centre|centroid"),
        ("position", "Shape offset (Angstrom)", "vector", "0 0 0", ""),
        ("vacuum", "Vacuum around a carved shape (Angstrom)", "float", "10", ""),
        ("max_atoms", "At most this many atoms", "float", "20000000", ""),
        ("combine", "Combine with the selection (select)", "choice", "replace", MODES),
    ]);
    list.extend(params(&[
        ("mesh_file", "Write the shape surface to (OBJ; empty = no)", "string", "", ""),
        ("mesh_spacing", "Surface mesh spacing (Angstrom)", "float", "0.5", ""),
    ]));
    vec![ModifierType::new(
        "build-shape",
        "Shape Modeller",
        "Build",
        "Carves a crystal into any shape modelled as in CAD: primitives (sphere, box, cylinder, cone, torus, capsule, ellipsoid, prism, extruded polygon, half-space, slab and crystallographic {hkl} facets), booleans (union, intersection, difference, smooth union), transforms (translate, rotate, scale, mirror), shells, offsets and linear, polar and grid arrays, written as an OpenSCAD-like script with variables. Regions give the atoms inside them another element (core-shell particles, inclusions, Janus particles; \"none\" removes them). Carve fills the shape with the crystal (this structure's cell, repeated) in an open box; keep, delete and select use the structure's own atoms (voids, notches, pores). The shape's surface can be written as OBJ.",
        list,
        build_shape,
    )]
}
