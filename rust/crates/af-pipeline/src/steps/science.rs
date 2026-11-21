//! Scientific tools that produce a structure, as pipeline steps (the C++
//! scienceStepTypes of pipeline/Modifiers): the request is the step's
//! parameters (JSON-valued ones as text) plus the step's structure.

use std::path::Path;

use af_core::json::Json;
use af_science::catalog::find_tool;
use af_science::data::structure_json;
use af_science::tools::run_tool;

use crate::pipeline::{find_modifier_type, parameter, ModifierParameter, ModifierType, PipelineData};
use crate::support::{place, position_of};

struct ScienceStep {
    /// Modifier id.
    id: &'static str,
    /// Catalog tool.
    tool: &'static str,
    title: &'static str,
    /// Which output frame: "last", or a symmetry cell chosen by the "cell" parameter.
    frame: &'static str,
}

static STEPS: [ScienceStep; 4] = [
    ScienceStep { id: "relax", tool: "relax", title: "Relax structure (interatomic potential)", frame: "last" },
    ScienceStep { id: "nvt-dynamics", tool: "nvt", title: "NVT molecular dynamics (final frame)", frame: "last" },
    ScienceStep { id: "npt-dynamics", tool: "npt", title: "NPT molecular dynamics (final frame)", frame: "last" },
    ScienceStep { id: "standardize-cell", tool: "symmetry", title: "Standardize cell (symmetry)", frame: "cell" },
];

fn scalar_kind(kind: &str) -> bool {
    kind == "float" || kind == "int" || kind == "bool"
}

fn run_science_step(step: &ScienceStep, p: &Json, data: &mut PipelineData) -> Result<(), String> {
    let tool = find_tool(step.tool).ok_or_else(|| format!("Unknown scientific tool: {}", step.tool))?;
    let t = find_modifier_type(step.id).ok_or_else(|| format!("Unknown step {}", step.id))?;
    let mut request = Json::object();
    for definition in tool.parameters {
        let kind = definition.kind;
        if kind == "structure" {
            request.set(definition.name, structure_json(&data.structure));
            continue;
        }
        let value = parameter(p, &t, definition.name)?;
        if scalar_kind(kind) {
            request.set(definition.name, value);
            continue;
        }
        let text = if value.is_string() { value.string()?.to_string() } else { value.dump(-1) };
        if text.chars().all(|c| c == ' ' || c == '\t') {
            continue; // optional and empty
        }
        let parsed = Json::parse(&text).map_err(|error| format!("{}: {error}", definition.label))?;
        request.set(definition.name, parsed);
    }
    let output = run_tool(step.tool, &request, Path::new("."), None)?;
    if output.frames.is_empty() {
        return Err(format!("{} produced no structure", tool.title));
    }
    let mut index = output.frames.len() - 1;
    if step.frame == "cell" {
        let wanted = parameter(p, &t, "cell")?.string()?.to_string();
        if let Some(order) = output.result.find("frame_order") {
            for (i, name) in order.items()?.iter().enumerate() {
                if name.string()? == wanted {
                    index = i;
                }
            }
        }
    }
    let mut frames = output.frames;
    let result = frames.swap_remove(index);
    // Same atoms (relaxation, dynamics): keep the selection and per-atom metadata.
    if result.atoms.len() == data.structure.atoms.len() {
        let s = &mut data.structure;
        for (atom, new) in s.atoms.iter_mut().zip(&result.atoms) {
            place(atom, position_of(new));
        }
        s.cell_vectors = result.cell_vectors;
        s.has_unit_cell = result.has_unit_cell;
    } else {
        data.structure = result;
        data.clear_selection();
    }
    data.notes.push(format!("{}: {} atoms", tool.title, data.structure.atoms.len()));
    Ok(())
}

/// The steps' types, with the tool's own parameters (structure inputs come from the pipeline).
pub fn types() -> Vec<ModifierType> {
    let mut types = Vec::new();
    for step in &STEPS {
        let Some(tool) = find_tool(step.tool) else { continue };
        let mut parameters: Vec<ModifierParameter> = Vec::new();
        for definition in tool.parameters {
            if definition.kind == "structure" {
                continue;
            }
            // JSON-valued inputs (potentials, lists) are edited as their JSON text.
            let kind = if scalar_kind(definition.kind) { definition.kind } else { "string" };
            parameters.push(ModifierParameter::new(definition.name, definition.label, kind, definition.value, ""));
        }
        if step.frame == "cell" {
            parameters.push(ModifierParameter::new("cell", "Cell to keep", "choice", "primitive", "primitive|conventional|symmetrized"));
        }
        let first_line = tool.help.split('\n').next().unwrap_or("");
        let help = format!(
            "Runs {} (Analyze > {}) on the structure and keeps its {}. {first_line}",
            tool.title,
            tool.category,
            if step.frame == "cell" { "chosen standardized cell" } else { "final structure" }
        );
        types.push(ModifierType::new(step.id, step.title, "Simulation", &help, parameters, move |p: &Json, data: &mut PipelineData| run_science_step(step, p, data)));
    }
    types
}
