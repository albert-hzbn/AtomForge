//! Helpers shared by the modifier implementations (the C++
//! pipeline/ModifierSupport): parameter access, atom bookkeeping and selections.

use af_core::elements;
use af_core::json::Json;
use af_core::{AtomSite, Structure};
use af_science::core::{fractional, inverse, row_times, Mat3, Vec3};

use crate::expression::{AtomVariables, Expression};
use crate::pipeline::{find_modifier_type, parameter, ModifierParameter, PipelineData};

/// A parameter of the step `id` (its value, or its default).
pub fn value(p: &Json, id: &str, name: &str) -> Result<Json, String> {
    let t = find_modifier_type(id).ok_or_else(|| format!("Unknown modifier '{id}'"))?;
    parameter(p, &t, name)
}

pub fn number(p: &Json, id: &str, name: &str) -> Result<f64, String> {
    value(p, id, name)?.number()
}

pub fn flag(p: &Json, id: &str, name: &str) -> Result<bool, String> {
    value(p, id, name)?.boolean()
}

pub fn text(p: &Json, id: &str, name: &str) -> Result<String, String> {
    Ok(value(p, id, name)?.string()?.to_string())
}

/// An integer parameter of at least `minimum`.
pub fn integer_of(p: &Json, id: &str, name: &str, minimum: i32) -> Result<i32, String> {
    let v = number(p, id, name)?;
    if v != v.floor() || v < minimum as f64 {
        return Err(format!("{id}: {name} must be an integer >= {minimum}"));
    }
    Ok(v as i32)
}

pub fn vector(p: &Json, id: &str, name: &str) -> Result<Vec3, String> {
    let v = value(p, id, name)?;
    let items = v.items()?;
    let at = |k: usize| -> Result<f64, String> { items.get(k).unwrap_or(&Json::Null).number() };
    Ok([at(0)?, at(1)?, at(2)?])
}

/// How a selection step combines with the current selection.
pub const MODES: &str = "replace|add|subtract|intersect";

pub fn mode_parameter() -> ModifierParameter {
    ModifierParameter::new("mode", "Combine with the current selection", "choice", "replace", MODES)
}

pub fn cell_of(s: &Structure) -> Mat3 {
    s.cell_vectors
}

pub fn position_of(a: &AtomSite) -> Vec3 {
    [a.x, a.y, a.z]
}

pub fn place(a: &mut AtomSite, x: Vec3) {
    a.x = x[0];
    a.y = x[1];
    a.z = x[2];
}

pub fn require_cell(data: &PipelineData, what: &str) -> Result<(), String> {
    if !data.structure.has_unit_cell {
        return Err(format!("{what} needs a structure with a unit cell"));
    }
    Ok(())
}

/// Keeps atoms by index (in the given order), carrying aligned per-atom metadata.
pub fn rebuild(data: &mut PipelineData, keep: &[usize]) {
    let s = &data.structure;
    let n = s.atoms.len();
    let mut out = s.clone();
    out.atoms.clear();
    out.grain_colors.clear();
    out.grain_region_ids.clear();
    out.atom_property.clear();
    out.frozen_property.clear();
    let mut selected = Vec::with_capacity(keep.len());
    for &i in keep {
        out.atoms.push(s.atoms[i].clone());
        if s.grain_colors.len() == n {
            out.grain_colors.push(s.grain_colors[i]);
        }
        if s.grain_region_ids.len() == n {
            out.grain_region_ids.push(s.grain_region_ids[i]);
        }
        if s.atom_property.len() == n {
            out.atom_property.push(s.atom_property[i]);
        }
        if s.frozen_property.len() == n {
            out.frozen_property.push(s.frozen_property[i]);
        }
        selected.push(data.selected.get(i).copied().unwrap_or(0));
    }
    out.dislocation_loop_points.clear();
    out.dislocation_detection_done = false;
    data.structure = out;
    data.selected = selected;
}

/// Indices the operation applies to: the selection, or every atom.
pub fn targets(data: &PipelineData, selected_only: bool, what: &str) -> Result<Vec<usize>, String> {
    let out: Vec<usize> = (0..data.structure.atoms.len()).filter(|&i| !selected_only || data.selected[i] != 0).collect();
    if selected_only && out.is_empty() {
        return Err(format!("{what}: no atoms are selected; add a selection modifier above"));
    }
    Ok(out)
}

/// Combines a new selection with the current one by `mode`.
pub fn select(data: &mut PipelineData, mask: &[u8], mode: &str) {
    for (i, &m) in mask.iter().enumerate() {
        let s = &mut data.selected[i];
        *s = match mode {
            "replace" => m,
            "add" => (*s != 0 || m != 0) as u8,
            "subtract" => (*s != 0 && m == 0) as u8,
            _ => (*s != 0 && m != 0) as u8,
        };
    }
    let count = data.selected_count();
    data.notes.push(format!("{count} of {} atoms selected", mask.len()));
}

/// The shortest image of a vector along the periodic cell vectors.
pub fn minimum_image(s: &Structure, d: Vec3) -> Result<Vec3, String> {
    if !s.has_unit_cell {
        return Ok(d);
    }
    let cell = cell_of(s);
    let mut f = fractional(d, &cell)?;
    for k in 0..3 {
        if s.pbc[k] {
            f[k] -= f[k].round();
        }
    }
    Ok(row_times(f, &cell))
}

/// An element symbol in any case ("cu", "CU" -> "Cu").
pub fn standard_element(name: &str) -> Result<String, String> {
    let lower = name.to_ascii_lowercase();
    for z in 1..=118 {
        if elements::symbol(z).to_ascii_lowercase() == lower {
            return Ok(elements::symbol(z).to_string());
        }
    }
    Err(format!("Unknown element {name}"))
}

/// The value of an expression for every atom (variables as in select-expression).
pub fn evaluate_for_atoms(expression: &Expression, data: &PipelineData) -> Result<Vec<f64>, String> {
    let s = &data.structure;
    let inv = if s.has_unit_cell { inverse(&cell_of(s))? } else { [[0.0; 3]; 3] };
    let n = s.atoms.len();
    let mut values = Vec::with_capacity(n);
    let mut v = AtomVariables { count: n as f64, ..Default::default() };
    for (i, a) in s.atoms.iter().enumerate() {
        v.x = a.x;
        v.y = a.y;
        v.z = a.z;
        if s.has_unit_cell {
            let f = row_times(position_of(a), &inv);
            v.fx = f[0];
            v.fy = f[1];
            v.fz = f[2];
        }
        v.index = i as f64;
        v.property = if s.atom_property.len() == n { s.atom_property[i] } else { f64::NAN };
        v.frozen = if s.frozen_property.len() == n { s.frozen_property[i] } else { f64::NAN };
        v.selected = if data.selected[i] != 0 { 1.0 } else { 0.0 };
        v.atomic_number = a.atomic_number as f64;
        v.fixed = ((a.fixed & 1) + ((a.fixed >> 1) & 1) + ((a.fixed >> 2) & 1)) as f64;
        v.element = a.symbol.clone();
        v.magmom = a.magmom;
        v.charge = a.charge;
        v.vx = a.velocity[0];
        v.vy = a.velocity[1];
        v.vz = a.velocity[2];
        values.push(expression.evaluate(&v)?);
    }
    Ok(values)
}
