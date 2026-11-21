//! The scientific tool dispatcher (the C++ science/ScienceTools): a tool id
//! and a JSON request are validated against the catalog and run; results
//! are described for the report pane. Tools are added as their modules are
//! ported; the rest say so.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::sync::OnceLock;

use af_core::json::Json;
use af_core::suggest;
use af_core::Structure;

use crate::analysis;
use crate::atom_properties::{per_atom_properties, AtomProperty};
use crate::bonds_regions;
use crate::clusters;
use crate::diffraction;
use crate::dislocation_lines;
use crate::lammps_export;
use crate::local_environment;
use crate::property_analysis;
use crate::reciprocal_path;
use crate::surface_mesh;
use crate::symmetry;
use crate::template_matching;
use crate::trajectory_structure;
use crate::trajectory_tools;
use crate::vasp_electronic;
use crate::voronoi;
use crate::workflows;
use crate::{phonons, potentials::potential_factory, simulation};
use crate::catalog::{catalog, find_tool, ToolDef};
use crate::core::{Pbc, Vec3};
use crate::data::{Parameters, StructureReader};

/// A tool's numerical result and, for some tools, structures.
#[derive(Clone, Debug, Default)]
pub struct ToolOutput {
    pub result: Json,
    /// Images, trajectory frames or a primitive cell.
    pub frames: Vec<Structure>,
    /// Angstrom/fs per frame (dynamics only).
    pub velocities: Vec<Vec<Vec3>>,
    /// fs per frame (dynamics only).
    pub times: Vec<f64>,
}

type Runner = fn(&Parameters, &Path) -> Result<ToolOutput, String>;

fn numeric(analyse: fn(&Parameters) -> Result<Json, String>, p: &Parameters) -> Result<ToolOutput, String> {
    Ok(ToolOutput { result: analyse(p)?, ..Default::default() })
}

macro_rules! numeric_tool {
    ($f:path) => {
        (|p: &Parameters, _: &Path| numeric($f, p)) as Runner
    };
}

fn runners() -> &'static BTreeMap<&'static str, Runner> {
    static RUNNERS: OnceLock<BTreeMap<&'static str, Runner>> = OnceLock::new();
    RUNNERS.get_or_init(|| {
        BTreeMap::from([
            ("msd", numeric_tool!(analysis::mean_square_displacement)),
            ("diffusion", numeric_tool!(analysis::diffusion_coefficient)),
            ("vacf", numeric_tool!(analysis::velocity_autocorrelation)),
            ("vibrational-spectrum", numeric_tool!(analysis::vibrational_spectrum)),
            ("local-strain", numeric_tool!(analysis::local_strain)),
            ("centrosymmetry", numeric_tool!(analysis::centrosymmetry)),
            ("structure-type", numeric_tool!(analysis::structure_type)),
            ("diamond-structure", numeric_tool!(analysis::diamond_structure)),
            ("bond-order", numeric_tool!(analysis::bond_order)),
            ("wigner-seitz", numeric_tool!(analysis::wigner_seitz)),
            ("structure-factor", numeric_tool!(analysis::static_structure_factor)),
            ("band-gap", numeric_tool!(analysis::band_gap)),
            ("effective-mass", numeric_tool!(analysis::effective_mass)),
            ("work-function", numeric_tool!(analysis::work_function)),
            ("equation-of-state", numeric_tool!(analysis::equation_of_state)),
            ("elastic-tensor", numeric_tool!(analysis::elastic_tensor)),
            ("phonon-dos", numeric_tool!(analysis::phonon_dos)),
            ("harmonic-thermodynamics", numeric_tool!(analysis::harmonic_thermodynamics)),
            ("trajectory-average", (|p: &Parameters, _: &Path| trajectory_tools::trajectory_average(p)) as Runner),
            ("unwrap-trajectory", (|p: &Parameters, _: &Path| trajectory_tools::unwrap_trajectory(p)) as Runner),
            ("frame-difference", numeric_tool!(trajectory_tools::frame_difference)),
            ("time-series", numeric_tool!(trajectory_tools::time_series)),
            ("spatial-binning", numeric_tool!(property_analysis::spatial_binning)),
            ("property-histogram", numeric_tool!(property_analysis::property_histogram)),
            ("scatter-plot", numeric_tool!(property_analysis::scatter_plot)),
            ("ambient-occlusion", numeric_tool!(property_analysis::ambient_occlusion)),
            ("spatial-correlation", numeric_tool!(property_analysis::spatial_correlation)),
            ("chill-plus", numeric_tool!(local_environment::chill_plus)),
            ("local-entropy", numeric_tool!(local_environment::local_entropy)),
            ("find-rings", numeric_tool!(local_environment::find_rings)),
            ("displacements", numeric_tool!(local_environment::displacement_vectors)),
            ("bond-distributions", numeric_tool!(local_environment::bond_distributions)),
            ("cluster-analysis", numeric_tool!(clusters::cluster_analysis)),
            ("void-analysis", numeric_tool!(clusters::void_analysis)),
            ("bond-orders", numeric_tool!(bonds_regions::bond_orders)),
            ("lammps-regions", numeric_tool!(bonds_regions::lammps_regions)),
            ("surface-mesh", numeric_tool!(surface_mesh::surface_mesh)),
            ("voronoi", numeric_tool!(voronoi::voronoi_analysis)),
            ("vorotop", numeric_tool!(voronoi::voro_top)),
            ("cna", numeric_tool!(voronoi::adaptive_cna)),
            ("ptm", numeric_tool!(template_matching::polyhedral_template_matching)),
            ("grain-segmentation", numeric_tool!(template_matching::grain_segmentation)),
            ("elastic-strain", numeric_tool!(template_matching::elastic_strain)),
            ("fcc-planar-faults", numeric_tool!(template_matching::fcc_planar_faults)),
            ("powder-xrd", (|p: &Parameters, _: &Path| diffraction::powder_xrd(p.structure("structure")?, p)) as Runner),
            ("electron-diffraction", (|p: &Parameters, _: &Path| diffraction::electron_diffraction(p.structure("structure")?, p)) as Runner),
            ("dft-inputs", (|p: &Parameters, _: &Path| {
                // GCC evaluates the C++ call's arguments right to left.
                let time_reversal = p.boolean("time_reversal", true)?;
                let symprec = p.number_or("symprec_A", 1e-5)?;
                let points = crate::core::integer(p.number_or("points_per_segment", 40.0)?, "points_per_segment", 2)? as i32 as i64;
                reciprocal_path::dft_inputs(p.structure("structure")?, points, symprec, time_reversal)
            }) as Runner),
            ("reciprocal-path", (|p: &Parameters, _: &Path| {
                let time_reversal = p.boolean("time_reversal", true)?;
                let symprec = p.number_or("symprec_A", 1e-5)?;
                let spacing = p.number_or("spacing_inv_A", 0.025)?;
                reciprocal_path::reciprocal_path(p.structure("structure")?, spacing, symprec, time_reversal)
            }) as Runner),
            ("dislocation-lines", (|p: &Parameters, _: &Path| {
                // GCC evaluates the C++ call's arguments right to left.
                let deformed = p.structure("structure")?;
                dislocation_lines::dislocation_lines(p.structure("reference")?, deformed, p)
            }) as Runner),
            ("trajectory-structure", numeric_tool!(trajectory_structure::trajectory_structure)),
            ("lammps-export", (|p: &Parameters, _: &Path| lammps_export::lammps_export(p)) as Runner),
            ("vasp-electronic", (|p: &Parameters, _: &Path| vasp_electronic::vasp_electronic(p)) as Runner),
            // GCC evaluates the C++ calls' arguments right to left: the calculator before the structures.
            ("elastic-constants", (|p: &Parameters, base: &Path| {
                let potential = potential_factory(p.json("calculator")?, base)?;
                workflows::elastic_constants(p.structure("structure")?, potential.as_ref(), p)
            }) as Runner),
            ("eos-scan", (|p: &Parameters, base: &Path| {
                let potential = potential_factory(p.json("calculator")?, base)?;
                workflows::equation_of_state_scan(p.structure("structure")?, potential.as_ref(), p)
            }) as Runner),
            ("formation-energy", (|p: &Parameters, base: &Path| {
                let potential = potential_factory(p.json("calculator")?, base)?;
                let bulk = p.structure("bulk")?;
                workflows::formation_energy(p.structure("structure")?, bulk, potential.as_ref(), p)
            }) as Runner),
            ("planar-defect-energy", (|p: &Parameters, base: &Path| {
                let potential = potential_factory(p.json("calculator")?, base)?;
                let bulk = p.structure("bulk")?;
                workflows::planar_defect_energy(p.structure("structure")?, bulk, potential.as_ref(), p)
            }) as Runner),
            ("neb", simulation::run_neb as Runner),
            ("relax", simulation::run_relax as Runner),
            ("phonons", (|p: &Parameters, base: &Path| {
                let potential = potential_factory(p.json("calculator")?, base)?;
                phonons::phonon_calculation(p.structure("structure")?, potential.as_ref(), p)
            }) as Runner),
            ("nvt", (|p: &Parameters, base: &Path| simulation::run_dynamics(p, base, false)) as Runner),
            ("npt", (|p: &Parameters, base: &Path| simulation::run_dynamics(p, base, true)) as Runner),
            ("symmetry", (|p: &Parameters, _: &Path| symmetry::symmetry_analysis(p.structure("structure")?, p)) as Runner),
        ])
    })
}

fn definition(tool: &str) -> Result<&'static ToolDef, String> {
    find_tool(tool).ok_or_else(|| {
        let ids: Vec<&str> = catalog().iter().map(|t| t.id).collect();
        format!("Unknown scientific tool: {tool}{}; AtomForge --science --help lists them", suggest::did_you_mean(tool, ids))
    })
}

fn validate_keys(tool: &ToolDef, request: &Json) -> Result<(), String> {
    let Json::Object(members) = request else {
        return Err("parameters must be an object".into());
    };
    let known: BTreeSet<&str> = tool.parameters.iter().map(|p| p.name).collect();
    for (name, _) in members {
        if !known.contains(name.as_str()) {
            return Err(format!(
                "{} has no parameter '{}'{}; it takes {}",
                tool.id,
                name,
                suggest::did_you_mean(name, known.iter().copied()),
                suggest::list(known.iter().copied())
            ));
        }
    }
    for parameter in tool.parameters {
        if parameter.required && !members.iter().any(|(k, _)| k == parameter.name) {
            let example = if parameter.value.is_empty() { String::new() } else { format!(", for example \"{}\": {}", parameter.name, parameter.value) };
            return Err(format!("{} requires parameter '{}'{}", tool.id, parameter.name, example));
        }
    }
    Ok(())
}

/// Runs a catalog tool on a JSON request (files relative to `base`).
pub fn run_tool(tool: &str, request: &Json, base: &Path, reader: StructureReader) -> Result<ToolOutput, String> {
    let def = definition(tool)?;
    validate_keys(def, request)?;
    let p = Parameters::new(tool, request, base, reader)?;
    let Some(runner) = runners().get(tool) else {
        return Err(format!("{tool} is not ported to the Rust build yet"));
    };
    runner(&p, base)
}

/// The tools that run (sorted).
pub fn runnable_tools() -> Vec<String> {
    runners().keys().map(|k| k.to_string()).collect()
}

pub fn periodic_axes(structure: &Structure) -> Pbc {
    let c = structure.has_unit_cell;
    [c && structure.pbc[0], c && structure.pbc[1], c && structure.pbc[2]]
}

/// {"positions", "cell", "pbc"} of a structure.
pub fn atoms_request(structure: &Structure, pbc: &Pbc) -> Json {
    let mut request = Json::object();
    request.set("positions", Json::Array(structure.atoms.iter().map(|a| Json::Array(vec![a.x.into(), a.y.into(), a.z.into()])).collect()));
    if structure.has_unit_cell {
        request.set("cell", Json::Array(structure.cell_vectors.iter().map(|r| Json::Array(vec![r[0].into(), r[1].into(), r[2].into()])).collect()));
        request.set("pbc", Json::Array(vec![pbc[0].into(), pbc[1].into(), pbc[2].into()]));
    }
    request
}

/// The request of a per-atom analysis of a structure (cutoff for tools that
/// need one, species for tools that take them).
pub fn per_atom_request(tool: &str, structure: &Structure, pbc: &Pbc, cutoff: f64) -> Result<Json, String> {
    let def = definition(tool)?;
    let parameter = |name: &str| def.parameters.iter().find(|p| p.name == name);
    let mut request = atoms_request(structure, pbc);
    if let Some(c) = parameter("cutoff_A") {
        let automatic_needs_cell = tool == "structure-type" || tool == "diamond-structure";
        if c.required || (automatic_needs_cell && !structure.has_unit_cell) {
            request.set("cutoff_A", cutoff);
        }
    }
    if parameter("types").is_some() {
        request.set("types", Json::Array(structure.atoms.iter().map(|a| Json::Number(a.atomic_number as f64)).collect()));
    }
    if parameter("local_density").is_some() && !(pbc[0] && pbc[1] && pbc[2]) {
        request.set("local_density", true);
    }
    Ok(request)
}

/// One per-atom property of a structure (the first, or the named one).
pub fn analyse_per_atom(structure: &Structure, tool: &str, cutoff: f64, property: &str) -> Result<AtomProperty, String> {
    let mut request = per_atom_request(tool, structure, &periodic_axes(structure), cutoff)?;
    if tool == "bond-order" {
        request.set("degrees", Json::Array(vec![Json::Number(6.0)]));
    }
    let properties = per_atom_properties(tool, &run_tool(tool, &request, Path::new("."), None)?.result);
    for candidate in properties {
        if property.is_empty() || candidate.name == property {
            return Ok(candidate);
        }
    }
    Err(format!("{tool}{}", if property.is_empty() { " returned no per-atom values".to_string() } else { format!(" has no per-atom property {property}") }))
}

// ---------------------------------------------------------------- reports
fn is_structure(value: &Json) -> bool {
    matches!(value, Json::Object(_)) && value.find("symbols").is_some() && value.find("positions").is_some()
}

fn is_numeric(value: &Json) -> bool {
    match value {
        Json::Array(items) => items.iter().all(is_numeric),
        Json::Number(_) | Json::Null | Json::Bool(_) => true,
        _ => false,
    }
}

fn collect(value: &Json, out: &mut Vec<f64>) {
    match value {
        Json::Array(items) => items.iter().for_each(|i| collect(i, out)),
        Json::Number(n) => out.push(*n),
        Json::Bool(b) => out.push(if *b { 1.0 } else { 0.0 }),
        _ => out.push(f64::NAN),
    }
}

fn number(value: f64) -> String {
    if !value.is_finite() {
        return "nan".into();
    }
    af_core::format::format_g(value, 7)
}

fn compact(value: &Json) -> String {
    match value {
        Json::Array(items) => format!("[{}]", items.iter().map(compact).collect::<Vec<_>>().join(" ")),
        Json::Bool(b) => (if *b { "True" } else { "False" }).into(),
        Json::Number(n) => number(*n),
        _ => "nan".into(),
    }
}

fn shape(value: &Json) -> String {
    let mut dims = Vec::new();
    let mut current = value;
    while let Json::Array(items) = current {
        dims.push(items.len().to_string());
        if items.is_empty() {
            break;
        }
        current = &items[0];
    }
    match dims.len() {
        0 => "()".into(),
        1 => format!("({},)", dims[0]),
        _ => format!("({})", dims.join(", ")),
    }
}

fn describe(key: &str, value: &Json, lines: &mut Vec<String>) {
    if is_structure(value) {
        lines.push(format!("{key}: {} atoms", value.find("symbols").map_or(0, |s| s.size())));
        return;
    }
    match value {
        Json::Object(members) => {
            for (child, item) in members {
                describe(&if key.is_empty() { child.clone() } else { format!("{key}.{child}") }, item, lines);
            }
        }
        Json::Array(items) if !items.is_empty() && is_structure(&items[0]) => lines.push(format!("{key}: {} structural frames", items.len())),
        Json::Array(items) if !items.is_empty() && is_numeric(value) => {
            let mut data = Vec::new();
            collect(value, &mut data);
            lines.push(format!("{key} (shape {}):", shape(value)));
            if data.len() <= 36 {
                lines.push(compact(value));
            } else {
                let (mut low, mut high, mut sum, mut valid) = (f64::INFINITY, f64::NEG_INFINITY, 0.0, 0usize);
                for &v in &data {
                    if v.is_finite() {
                        if v < low {
                            low = v;
                        }
                        if high < v {
                            high = v;
                        }
                        sum += v;
                        valid += 1;
                    }
                }
                if valid > 0 {
                    lines.push(format!(
                        "  min={}, max={}, mean={}; {} valid / {} values",
                        number(low),
                        number(high),
                        number(sum / valid as f64),
                        valid,
                        data.len()
                    ));
                } else {
                    lines.push("  No valid values".into());
                }
            }
            let _ = items;
        }
        Json::Array(_) => {
            let mut text = value.dump(-1);
            text.truncate(1000);
            lines.push(format!("{key}: {text}"));
        }
        Json::String(text) => {
            let shown = if text.contains('\n') {
                format!("({}-line file; use Save generated files)", text.matches('\n').count())
            } else {
                text.clone()
            };
            lines.push(format!("{key}: {shown}"));
        }
        Json::Bool(b) => lines.push(format!("{key}: {}", if *b { "True" } else { "False" })),
        Json::Number(n) => lines.push(format!("{key}: {}", af_core::format::format_g(*n, 12))),
        Json::Null => lines.push(format!("{key}: None")),
    }
}

/// The report text of a result.
pub fn result_report(tool: &str, result: &Json) -> Result<String, String> {
    let mut lines = vec![definition(tool)?.title.to_string(), String::new(), "Full numerical data are available through Save results.".into(), String::new()];
    describe("", result, &mut lines);
    Ok(lines.iter().map(|l| format!("{l}\n")).collect())
}

/// {"tool", "parameters", "result"}.
pub fn result_document(tool: &str, request: &Json, output: &ToolOutput) -> Json {
    let mut document = Json::object();
    document.set("tool", tool);
    document.set("parameters", request.clone());
    document.set("result", output.result.clone());
    document
}
