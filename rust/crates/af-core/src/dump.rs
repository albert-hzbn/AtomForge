//! A complete JSON form of a structure (every field the algorithms touch),
//! as the golden C++ dumpers write it: lossless numbers ("nan" for NaN
//! per-atom values), so test inputs and outputs compare exactly.

use crate::json::Json;
use crate::{AtomSite, Structure};

fn vec3(v: [f64; 3]) -> Json {
    Json::Array(v.iter().map(|&x| Json::Number(x)).collect())
}

fn maybe_nan(v: f64) -> Json {
    if v.is_nan() {
        Json::String("nan".into())
    } else {
        Json::Number(v)
    }
}

pub fn structure_dump(s: &Structure) -> Json {
    let mut o = Json::object();
    let atoms = s
        .atoms
        .iter()
        .map(|a| {
            let mut atom = Json::object();
            atom.set("symbol", a.symbol.as_str());
            atom.set("z", a.atomic_number);
            atom.set("position", vec3([a.x, a.y, a.z]));
            atom.set("colour", vec3([a.r as f64, a.g as f64, a.b as f64]));
            atom.set("fixed", a.fixed as i32);
            atom.set("magmom", a.magmom);
            atom.set("charge", a.charge);
            atom.set("velocity", vec3(a.velocity));
            atom
        })
        .collect();
    o.set("atoms", Json::Array(atoms));
    o.set("has_cell", s.has_unit_cell);
    o.set("cell", Json::Array(s.cell_vectors.iter().map(|r| vec3(*r)).collect()));
    o.set("pbc", Json::Array(s.pbc.iter().map(|&b| Json::Bool(b)).collect()));
    o.set("offset", vec3(s.cell_offset));
    o.set("grain_colors", Json::Array(s.grain_colors.iter().map(|c| vec3([c[0] as f64, c[1] as f64, c[2] as f64])).collect()));
    o.set("grain_region_ids", Json::Array(s.grain_region_ids.iter().map(|&r| Json::Number(r as f64)).collect()));
    o.set("atom_property", Json::Array(s.atom_property.iter().map(|&v| maybe_nan(v)).collect()));
    o.set("frozen_property", Json::Array(s.frozen_property.iter().map(|&v| maybe_nan(v)).collect()));
    o.set("loop_points", s.dislocation_loop_points.len());
    o.set("detection_done", s.dislocation_detection_done);
    o
}

fn numbers(j: &Json) -> Result<Vec<f64>, String> {
    j.items()?.iter().map(|v| if let Json::String(_) = v { Ok(f64::NAN) } else { v.number() }).collect()
}

fn v3(j: &Json) -> Result<[f64; 3], String> {
    let n = numbers(j)?;
    if n.len() != 3 {
        return Err("Expected three numbers".into());
    }
    Ok([n[0], n[1], n[2]])
}

pub fn structure_from_dump(o: &Json) -> Result<Structure, String> {
    let mut s = Structure::default();
    for a in o.at("atoms")?.items()? {
        let colour = v3(a.at("colour")?)?;
        let position = v3(a.at("position")?)?;
        s.atoms.push(AtomSite {
            symbol: a.at("symbol")?.string()?.into(),
            atomic_number: a.at("z")?.number()? as i32,
            x: position[0],
            y: position[1],
            z: position[2],
            r: colour[0] as f32,
            g: colour[1] as f32,
            b: colour[2] as f32,
            fixed: a.at("fixed")?.number()? as u8,
            magmom: a.at("magmom")?.number()?,
            charge: a.at("charge")?.number()?,
            velocity: v3(a.at("velocity")?)?,
        });
    }
    s.has_unit_cell = o.at("has_cell")?.boolean()?;
    let cell = o.at("cell")?.items()?;
    s.cell_vectors = [v3(&cell[0])?, v3(&cell[1])?, v3(&cell[2])?];
    let pbc = o.at("pbc")?.items()?;
    s.pbc = [pbc[0].boolean()?, pbc[1].boolean()?, pbc[2].boolean()?];
    s.cell_offset = v3(o.at("offset")?)?;
    s.grain_colors = o.at("grain_colors")?.items()?.iter().map(|c| v3(c).map(|v| v.map(|x| x as f32))).collect::<Result<_, _>>()?;
    s.grain_region_ids = numbers(o.at("grain_region_ids")?)?.iter().map(|&r| r as i32).collect();
    s.atom_property = numbers(o.at("atom_property")?)?;
    s.frozen_property = numbers(o.at("frozen_property")?)?;
    Ok(s)
}

/// Whether two golden numbers agree: exactly on Windows, where the C++
/// references were made; elsewhere (another maths library, another spglib
/// build) to rounding, 2e-6 relative (single precision).
pub fn golden_number_matches(x: f64, y: f64) -> bool {
    if x == y || (x.is_nan() && y.is_nan()) {
        return true;
    }
    !cfg!(windows) && (x - y).abs() <= 2e-6 * x.abs().max(y.abs()).max(1.0)
}

/// Differences between two JSON values, as "path: a vs b" (numbers compared
/// by `golden_number_matches`).
pub fn json_diff(path: &str, a: &Json, b: &Json, out: &mut Vec<String>) {
    if out.len() > 40 {
        return;
    }
    match (a, b) {
        (Json::Object(x), Json::Object(y)) => {
            let kx: Vec<_> = x.iter().map(|(k, _)| k).collect();
            let ky: Vec<_> = y.iter().map(|(k, _)| k).collect();
            if kx != ky {
                out.push(format!("{path}: keys {kx:?} vs {ky:?}"));
                return;
            }
            for ((k, va), (_, vb)) in x.iter().zip(y) {
                json_diff(&format!("{path}.{k}"), va, vb, out);
            }
        }
        (Json::Array(x), Json::Array(y)) => {
            if x.len() != y.len() {
                out.push(format!("{path}: length {} vs {}", x.len(), y.len()));
                return;
            }
            for (i, (va, vb)) in x.iter().zip(y).enumerate() {
                json_diff(&format!("{path}[{i}]"), va, vb, out);
            }
        }
        (Json::Number(x), Json::Number(y)) if golden_number_matches(*x, *y) => {}
        _ if a == b => {}
        _ => out.push(format!("{path}: {} vs {}", a.dump(-1), b.dump(-1))),
    }
}
