//! Potential-driven workflows (the C++ science/Workflows): elastic constants
//! from strain-stress samples, equation-of-state scans, defect formation
//! energies and planar defect energies.

use std::collections::BTreeMap;
use std::path::Path;

use af_core::json::Json;

use crate::analysis::{elastic_tensor, equation_of_state, jv};
use crate::core::*;
use crate::data::{Parameters, StructureInput};
use crate::potentials::{configuration_from, structure_from, Configuration, Potential, PotentialResult};
use crate::simulation::{relax_configuration, RelaxOptions};
use crate::tools::ToolOutput;

/// The C++ Workflows' own conversion (CODATA 2014).
const EV_PER_A3_TO_GPA: f64 = 160.21766208;
const EV_PER_A2_TO_J_PER_M2: f64 = 16.02176634;

fn periodic_configuration(input: &StructureInput, tool: &str) -> Result<Configuration, String> {
    if !input.structure.has_unit_cell || !(input.pbc[0] && input.pbc[1] && input.pbc[2]) || input.structure.atoms.is_empty() {
        return Err(format!("{tool} needs a nonempty fully periodic structure"));
    }
    configuration_from(&input.structure, &input.pbc)
}

/// Homogeneous deformation: positions and cell rows times (I + strain).
fn deformed(c: &Configuration, deformation: &Mat3) -> Configuration {
    let mut d = c.clone();
    for x in d.positions.iter_mut() {
        *x = row_times(*x, deformation);
    }
    for r in 0..3 {
        d.cell[r] = row_times(c.cell[r], deformation);
    }
    d
}

/// Ions relaxed at a fixed cell (when requested): the final energy and stress, and convergence.
fn evaluate(c: &mut Configuration, potential: &dyn Potential, relax_ions: bool, fmax: f64, steps: i64) -> Result<(PotentialResult, bool), String> {
    if !relax_ions {
        return Ok((potential.compute(c, true)?, true));
    }
    let options = RelaxOptions { fmax, steps, ..Default::default() };
    let relaxed = relax_configuration(c, potential, &options)?;
    *c = relaxed.configuration;
    let forces = if relaxed.forces.has_stress { relaxed.forces } else { potential.compute(c, true)? };
    Ok((forces, relaxed.converged))
}

pub fn elastic_constants(input: &StructureInput, potential: &dyn Potential, p: &Parameters) -> Result<ToolOutput, String> {
    let reference = periodic_configuration(input, "Elastic constants")?;
    let max_strain = positive(p.number_or("max_strain", 0.01)?, "max_strain")?;
    if max_strain > 0.05 {
        return Err("max_strain should be at most 0.05 for the linear (small-strain) regime".into());
    }
    let amplitudes = integer(p.number_or("strain_steps", 2.0)?, "strain_steps", 1)?;
    let relax_ions = p.boolean("relax_ions", true)?;
    let fmax = positive(p.number_or("fmax", 1e-4)?, "fmax")?;
    let steps = integer(p.number_or("steps", 1000.0)?, "steps", 1)?;

    // Engineering-Voigt strains xx, yy, zz, yz, xz, xy; shear components are 2 eps_ij.
    let mut strains = Vec::new();
    let mut stresses = Vec::new();
    let mut unconverged: i64 = 0;
    let mut add_sample = |voigt: [f64; 6]| -> Result<(), String> {
        let mut deformation = identity();
        let pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let e = if k < 3 { voigt[k] } else { 0.5 * voigt[k] };
            deformation[a][b] += e;
            if a != b {
                deformation[b][a] += e;
            }
        }
        let mut c = deformed(&reference, &deformation);
        let (result, converged) = evaluate(&mut c, potential, relax_ions, fmax, steps)?;
        unconverged += (!converged) as i64;
        let s = &result.stress;
        strains.push(jv(&voigt));
        stresses.push(jv(&[
            s[0][0] * EV_PER_A3_TO_GPA,
            s[1][1] * EV_PER_A3_TO_GPA,
            s[2][2] * EV_PER_A3_TO_GPA,
            s[1][2] * EV_PER_A3_TO_GPA,
            s[0][2] * EV_PER_A3_TO_GPA,
            s[0][1] * EV_PER_A3_TO_GPA,
        ]));
        Ok(())
    };
    add_sample([0.0; 6])?;
    for component in 0..6 {
        for n in 1..=amplitudes {
            for sign in [-1.0, 1.0] {
                let mut voigt = [0.0; 6];
                voigt[component] = sign * max_strain * n as f64 / amplitudes as f64;
                add_sample(voigt)?;
            }
        }
    }
    let strains = Json::Array(strains);
    let stresses = Json::Array(stresses);
    let mut request = Json::object();
    request.set("strains", strains.clone());
    request.set("stresses_GPa", stresses.clone());
    let mut result = elastic_tensor(&Parameters::new("elastic-tensor", &request, Path::new("."), None)?)?;
    result.set("strains", strains);
    result.set("stresses_GPa", stresses);
    result.set("max_strain", max_strain);
    result.set("relaxed_ions", relax_ions);
    result.set("unconverged_relaxations", unconverged);
    result.set("potential", potential.description());
    Ok(ToolOutput { result, ..Default::default() })
}

pub fn equation_of_state_scan(input: &StructureInput, potential: &dyn Potential, p: &Parameters) -> Result<ToolOutput, String> {
    let reference = periodic_configuration(input, "The equation-of-state scan")?;
    let range = positive(p.number_or("volume_range", 0.1)?, "volume_range")?;
    if range >= 0.5 {
        return Err("volume_range must be below 0.5".into());
    }
    let points = integer(p.number_or("points", 11.0)?, "points", 5)?;
    let relax_ions = p.boolean("relax_ions", true)?;
    let fmax = positive(p.number_or("fmax", 1e-3)?, "fmax")?;
    let steps = integer(p.number_or("steps", 500.0)?, "steps", 1)?;
    let volume0 = cell_volume(&reference.cell);
    let (mut volumes, mut energies, mut pressures) = (Vec::new(), Vec::new(), Vec::new());
    let mut unconverged: i64 = 0;
    for n in 0..points {
        let fraction = 1.0 - range + 2.0 * range * n as f64 / (points - 1) as f64;
        let mut m = identity();
        for k in 0..3 {
            m[k][k] = fraction.cbrt();
        }
        let mut c = deformed(&reference, &m);
        let (result, converged) = evaluate(&mut c, potential, relax_ions, fmax, steps)?;
        unconverged += (!converged) as i64;
        volumes.push(volume0 * fraction);
        energies.push(result.energy);
        pressures.push(-(result.stress[0][0] + result.stress[1][1] + result.stress[2][2]) / 3.0 * EV_PER_A3_TO_GPA);
    }
    let mut request = Json::object();
    request.set("volumes_A3", jv(&volumes));
    request.set("energies_eV", jv(&energies));
    let mut result = equation_of_state(&Parameters::new("equation-of-state", &request, Path::new("."), None)?)?;
    result.set("volumes_A3", jv(&volumes));
    result.set("energies_eV", jv(&energies));
    result.set("pressures_GPa", jv(&pressures));
    result.set("atom_count", reference.len() as i64);
    result.set("relaxed_ions", relax_ions);
    result.set("unconverged_relaxations", unconverged);
    result.set("potential", potential.description());
    Ok(ToolOutput { result, ..Default::default() })
}

struct Relaxed {
    configuration: Configuration,
    energy: f64,
    converged: bool,
}

fn relaxed(start: &Configuration, potential: &dyn Potential, p: &Parameters, relax_cell: bool) -> Result<Relaxed, String> {
    if !p.boolean("relax", true)? && !relax_cell {
        return Ok(Relaxed { configuration: start.clone(), energy: potential.compute(start, false)?.energy, converged: true });
    }
    let options = RelaxOptions {
        fmax: positive(p.number_or("fmax", 1e-3)?, "fmax")?,
        steps: integer(p.number_or("steps", 2000.0)?, "steps", 1)?,
        relax_cell,
        ..Default::default()
    };
    let result = relax_configuration(start, potential, &options)?;
    Ok(Relaxed { energy: result.forces.energy, converged: result.converged, configuration: result.configuration })
}

fn composition(c: &Configuration) -> BTreeMap<String, i64> {
    let mut counts = BTreeMap::new();
    for symbol in &c.symbols {
        *counts.entry(symbol.clone()).or_insert(0) += 1;
    }
    counts
}

pub fn formation_energy(defect_input: &StructureInput, bulk_input: &StructureInput, potential: &dyn Potential, p: &Parameters) -> Result<ToolOutput, String> {
    let defect = periodic_configuration(defect_input, "The defect cell")?;
    let bulk = periodic_configuration(bulk_input, "The bulk reference")?;
    let bulk_relaxed = relaxed(&bulk, potential, p, false)?;
    let defect_relaxed = relaxed(&defect, potential, p, p.boolean("relax_cell", false)?)?;
    let (n_bulk, n_defect) = (composition(&bulk), composition(&defect));
    // Chemical potentials: given, else the bulk energy per atom of a one-element reference.
    let mut mu: BTreeMap<String, f64> = BTreeMap::new();
    if p.has("chemical_potentials_eV") {
        for (element, value) in p.json("chemical_potentials_eV")?.members()? {
            mu.insert(element.clone(), value.number()?);
        }
    }
    let mut elements = n_bulk.clone();
    for element in n_defect.keys() {
        elements.entry(element.clone()).or_insert(0);
    }
    let mut formation = defect_relaxed.energy - bulk_relaxed.energy;
    let (mut change, mut used) = (Json::object(), Json::object());
    for element in elements.keys() {
        let dn = n_defect.get(element).copied().unwrap_or(0) - n_bulk.get(element).copied().unwrap_or(0);
        change.set(element, dn);
        if dn == 0 {
            continue;
        }
        if !mu.contains_key(element) {
            if n_bulk.len() == 1 && n_bulk.contains_key(element) {
                mu.insert(element.clone(), bulk_relaxed.energy / bulk.len() as f64);
            } else {
                return Err(format!("Give chemical_potentials_eV for {element} (the bulk reference does not fix it)"));
            }
        }
        used.set(element, mu[element]);
        formation -= dn as f64 * mu[element];
    }
    let mut result = Json::object();
    result.set("formation_energy_eV", formation);
    result.set("defect_energy_eV", defect_relaxed.energy);
    result.set("bulk_energy_eV", bulk_relaxed.energy);
    result.set("bulk_energy_per_atom_eV", bulk_relaxed.energy / bulk.len() as f64);
    result.set("atom_change", change);
    result.set("chemical_potentials_eV", used);
    result.set("defect_converged", defect_relaxed.converged);
    result.set("bulk_converged", bulk_relaxed.converged);
    result.set("relaxation_volume_A3", cell_volume(&defect_relaxed.configuration.cell) - cell_volume(&defect.cell));
    result.set("potential", potential.description());
    Ok(ToolOutput { result, frames: vec![structure_from(&defect_relaxed.configuration)], ..Default::default() })
}

pub fn planar_defect_energy(cell_input: &StructureInput, bulk_input: &StructureInput, potential: &dyn Potential, p: &Parameters) -> Result<ToolOutput, String> {
    if !cell_input.structure.has_unit_cell || cell_input.structure.atoms.is_empty() {
        return Err("The defect cell needs atoms and lattice vectors".into());
    }
    let cell = configuration_from(&cell_input.structure, &cell_input.pbc)?;
    let bulk = periodic_configuration(bulk_input, "The bulk reference")?;
    let interfaces = integer(p.number_or("interfaces", 2.0)?, "interfaces", 1)?;
    let axis = integer(p.number_or("normal_axis", 2.0)?, "normal_axis", 0)? as i32 as usize;
    if axis > 2 {
        return Err("normal_axis must be 0 (a), 1 (b) or 2 (c)".into());
    }
    // The cell must hold a whole multiple of the bulk composition.
    let (n_cell, n_bulk) = (composition(&cell), composition(&bulk));
    let ratio = cell.len() as f64 / bulk.len() as f64;
    let mut stoichiometric = n_cell.len() == n_bulk.len();
    for (element, &count) in &n_bulk {
        stoichiometric = stoichiometric && n_cell.get(element).is_some_and(|&found| (found as f64 - ratio * count as f64).abs() < 1e-6);
    }
    if !stoichiometric {
        return Err("The defect cell is not stoichiometric with the bulk reference; planar energies need the same composition".into());
    }
    let bulk_relaxed = relaxed(&bulk, potential, p, false)?;
    let cell_relaxed = relaxed(&cell, potential, p, false)?;
    let (u, v) = (cell.cell[(axis + 1) % 3], cell.cell[(axis + 2) % 3]);
    let area = norm(cross(u, v));
    let excess = cell_relaxed.energy - ratio * bulk_relaxed.energy;
    let gamma = excess / (interfaces as f64 * area);
    let mut result = Json::object();
    result.set("energy_eV_per_A2", gamma);
    result.set("energy_J_per_m2", gamma * EV_PER_A2_TO_J_PER_M2);
    result.set("energy_mJ_per_m2", gamma * EV_PER_A2_TO_J_PER_M2 * 1000.0);
    result.set("excess_energy_eV", excess);
    result.set("area_A2", area);
    result.set("interfaces", interfaces);
    result.set("cell_energy_eV", cell_relaxed.energy);
    result.set("bulk_energy_per_atom_eV", bulk_relaxed.energy / bulk.len() as f64);
    result.set("cell_converged", cell_relaxed.converged);
    result.set("bulk_converged", bulk_relaxed.converged);
    result.set("potential", potential.description());
    Ok(ToolOutput { result, frames: vec![structure_from(&cell_relaxed.configuration)], ..Default::default() })
}
