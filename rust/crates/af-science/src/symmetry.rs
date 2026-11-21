//! Space-group analysis (the C++ science/Symmetry): spglib dataset, Wyckoff
//! orbits, symmetry operations, symmetrized positions and the conventional,
//! primitive and symmetrized cells.

use std::collections::BTreeMap;

use af_core::elements;
use af_core::json::Json;
use af_core::structure::{AtomSite, Structure};

use crate::analysis::ji;
use crate::core::*;
use crate::data::{structure_json, Parameters, StructureInput};
use crate::tools::ToolOutput;

fn crystal_system(number: i32) -> &'static str {
    match number {
        ..=2 => "triclinic",
        ..=15 => "monoclinic",
        ..=74 => "orthorhombic",
        ..=142 => "tetragonal",
        ..=167 => "trigonal",
        ..=194 => "hexagonal",
        _ => "cubic",
    }
}

fn wrap(x: f64) -> f64 {
    x - x.round()
}

/// A periodic structure from a cell (rows), fractional positions and species.
pub fn make_structure(cell: &Mat3, fractions: &[Vec3], numbers: &[i32]) -> Structure {
    let mut structure = Structure { has_unit_cell: true, cell_vectors: *cell, ..Default::default() };
    for (f, &z) in fractions.iter().zip(numbers) {
        let position = row_times(*f, cell);
        let [r, g, b] = elements::default_colour(z);
        structure.atoms.push(AtomSite {
            atomic_number: z,
            symbol: elements::symbol(z).to_string(),
            x: position[0],
            y: position[1],
            z: position[2],
            r,
            g,
            b,
            ..Default::default()
        });
    }
    structure
}

fn lattice_parameters(cell: &Mat3) -> Json {
    let (a, b, c) = (norm(cell[0]), norm(cell[1]), norm(cell[2]));
    let degrees = 180.0 / (-1.0f64).acos();
    let mut result = Json::object();
    result.set("a_A", a);
    result.set("b_A", b);
    result.set("c_A", c);
    result.set("alpha_deg", (dot(cell[1], cell[2]) / (b * c)).clamp(-1.0, 1.0).acos() * degrees);
    result.set("beta_deg", (dot(cell[0], cell[2]) / (a * c)).clamp(-1.0, 1.0).acos() * degrees);
    result.set("gamma_deg", (dot(cell[0], cell[1]) / (a * b)).clamp(-1.0, 1.0).acos() * degrees);
    result.set("volume_A3", cell_volume(cell));
    result
}

/// spglib standardization (conventional or primitive, idealized), positions
/// wrapped into [0, 1).
pub fn standardized(cell: &Mat3, fractions: &[Vec3], numbers: &[i32], primitive: bool, symprec: f64) -> Result<Structure, String> {
    let mut c = af_spglib::standardize_cell(cell, fractions, numbers, primitive, false, symprec)
        .ok_or_else(|| "spglib could not standardize the cell".to_string())?;
    for f in c.positions.iter_mut() {
        for x in f.iter_mut() {
            *x -= x.floor();
        }
    }
    Ok(make_structure(&c.lattice, &c.positions, &c.types))
}

/// The transform taking Miller indices of `cell` to spglib's conventional
/// setting (h_s = h P^-1), or None.
pub fn conventional_index_transform(cell: &Mat3, fractions: &[Vec3], numbers: &[i32], symprec: f64) -> Option<Mat3> {
    if fractions.is_empty() || fractions.len() != numbers.len() {
        return None;
    }
    let dataset = af_spglib::get_dataset(cell, fractions, numbers, symprec)?;
    let p = dataset.transformation_matrix;
    if determinant(&p).abs() < 1e-9 {
        return None;
    }
    inverse(&p).ok()
}

pub fn symmetry_analysis(input: &StructureInput, p: &Parameters) -> Result<ToolOutput, String> {
    let structure = &input.structure;
    if !structure.has_unit_cell || !(input.pbc[0] && input.pbc[1] && input.pbc[2]) || structure.atoms.is_empty() {
        return Err("Symmetry analysis needs a nonempty periodic structure".into());
    }
    let symprec = positive(p.number_or("symprec_A", 1e-3)?, "symprec_A")?;
    let angle_tolerance = p.number_or("angle_tolerance_deg", -1.0)?;
    let output_cell = if p.has("output_cell") { p.json("output_cell")?.string()?.to_string() } else { "conventional".to_string() };
    if output_cell != "conventional" && output_cell != "primitive" && output_cell != "symmetrized" {
        return Err("output_cell must be conventional, primitive or symmetrized".into());
    }
    let cell = structure.cell_vectors;
    let inv = inverse(&cell)?;
    let fractions: Vec<Vec3> = structure.atoms.iter().map(|a| row_times([a.x, a.y, a.z], &inv)).collect();
    let numbers: Vec<i32> = structure.atoms.iter().map(|a| a.atomic_number).collect();
    let d = af_spglib::get_dataset_with_angle(&cell, &fractions, &numbers, symprec, angle_tolerance)
        .ok_or_else(|| "spglib could not detect the symmetry; try a larger symprec_A".to_string())?;

    let mut result = Json::object();
    let number = d.spacegroup_number;
    let n_operations = d.rotations.len();
    result.set("spacegroup_number", number);
    result.set("international_symbol", d.international_symbol.clone());
    result.set("hall_symbol", d.hall_symbol.clone());
    result.set("hall_number", d.hall_number);
    result.set("point_group", d.pointgroup_symbol.clone());
    result.set("crystal_system", crystal_system(number));
    result.set("operation_count", n_operations);
    result.set("symprec_A", symprec);

    let letter = |w: i32| ((b'a' as i32 + w) as u8 as char).to_string();
    let (mut wyckoff, mut site, mut orbit) = (Vec::new(), Vec::new(), Vec::new());
    let mut orbits = vec![0i32; numbers.len()];
    let mut letters = vec![0i32; numbers.len()];
    let mut orbit_index: BTreeMap<i32, i32> = BTreeMap::new();
    for i in 0..numbers.len() {
        letters[i] = d.wyckoffs[i];
        wyckoff.push(Json::from(letter(d.wyckoffs[i])));
        site.push(Json::from(d.site_symmetry_symbols[i].clone()));
        let representative = d.equivalent_atoms[i];
        let next = orbit_index.len() as i32;
        let index = *orbit_index.entry(representative).or_insert(next);
        orbits[i] = index;
        orbit.push(Json::from(index));
    }
    result.set("wyckoff_letters", Json::Array(wyckoff));
    result.set("site_symmetry", Json::Array(site));
    result.set("orbit_id", Json::Array(orbit));
    result.set("wyckoff_index", ji(&letters));
    result.set("orbit_count", orbit_index.len());
    let mut sites = Vec::new();
    for (&representative, &index) in &orbit_index {
        let multiplicity = orbits.iter().filter(|&&v| v == index).count();
        let r = representative as usize;
        let mut row = Json::object();
        row.set("element", elements::symbol(numbers[r]));
        row.set("wyckoff", letter(d.wyckoffs[r]));
        row.set("site_symmetry", d.site_symmetry_symbols[r].clone());
        row.set("multiplicity", multiplicity);
        row.set("representative_atom", representative);
        sites.push(row);
    }
    result.set("sites", Json::Array(sites));
    let mut operations = Vec::new();
    for k in 0..n_operations {
        let rotation = Json::Array(d.rotations[k].iter().map(|r| Json::Array(r.iter().map(|&v| Json::from(v)).collect())).collect());
        let mut op = Json::object();
        op.set("rotation", rotation);
        op.set("translation", Json::Array(d.translations[k].iter().map(|&v| Json::from(v)).collect()));
        operations.push(op);
    }
    result.set("operations", Json::Array(operations));

    // Symmetrized positions in the input cell: average each atom's images under all operations.
    let mut sum = vec![[0.0f64; 3]; numbers.len()];
    for k in 0..n_operations {
        let (rot, t) = (&d.rotations[k], &d.translations[k]);
        for i in 0..numbers.len() {
            let f = fractions[i];
            let image: Vec3 = std::array::from_fn(|r| rot[r][0] as f64 * f[0] + rot[r][1] as f64 * f[1] + rot[r][2] as f64 * f[2] + t[r]);
            let mut best = 0;
            let mut best_distance = f64::INFINITY;
            for j in 0..numbers.len() {
                if numbers[j] != numbers[i] {
                    continue;
                }
                let delta = [wrap(image[0] - fractions[j][0]), wrap(image[1] - fractions[j][1]), wrap(image[2] - fractions[j][2])];
                let distance = norm(row_times(delta, &cell));
                if distance < best_distance {
                    best_distance = distance;
                    best = j;
                }
            }
            let delta = [wrap(image[0] - fractions[best][0]), wrap(image[1] - fractions[best][1]), wrap(image[2] - fractions[best][2])];
            sum[best] = add(sum[best], add(fractions[best], delta));
        }
    }
    let mut largest_shift = 0.0f64;
    let symmetrized: Vec<Vec3> = sum.iter().map(|s| scale(*s, 1.0 / n_operations as f64)).collect();
    for i in 0..numbers.len() {
        let shift = norm(row_times(sub(symmetrized[i], fractions[i]), &cell));
        if largest_shift < shift {
            largest_shift = shift;
        }
    }
    result.set("max_symmetrization_shift_A", largest_shift);

    let conventional = standardized(&cell, &fractions, &numbers, false, symprec)?;
    let primitive = standardized(&cell, &fractions, &numbers, true, symprec)?;
    let symmetrized_structure = make_structure(&cell, &symmetrized, &numbers);
    result.set("input_cell", lattice_parameters(&cell));
    result.set("conventional_cell", lattice_parameters(&conventional.cell_vectors));
    result.set("primitive_cell", lattice_parameters(&primitive.cell_vectors));
    result.set("conventional_atom_count", conventional.atoms.len());
    result.set("primitive_atom_count", primitive.atoms.len());
    result.set("conventional_structure", structure_json(&conventional));
    result.set("primitive_structure", structure_json(&primitive));
    result.set("output_cell", output_cell.clone());
    // The requested cell last (std::stable_partition), the others in order.
    let cells = [("conventional", conventional), ("primitive", primitive), ("symmetrized", symmetrized_structure)];
    let (chosen, others): (Vec<_>, Vec<_>) = cells.into_iter().partition(|(name, _)| *name == output_cell);
    let mut output = ToolOutput::default();
    let mut order = Vec::new();
    for (name, frame) in others.into_iter().chain(chosen) {
        order.push(Json::from(name));
        output.frames.push(frame);
    }
    result.set("frame_order", Json::Array(order));
    output.result = result;
    Ok(output)
}
