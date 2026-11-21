//! `--pipe` and `--pipeline` (the C++ cli/PipelineCLI): non-destructive
//! structure-editing pipelines from text or JSON, on files or extended XYZ
//! streams.

use std::io::Read;

use af_core::json::Json;
use af_core::structure::Structure;
use af_io::loader::{load_structure_from_file, save_structure};
use af_pipeline::pipeline::{make_modifier, modifier_types, Pipeline, PipelineData};
use af_science::data::{read_extxyz_frames, LineReader};

use crate::args::{detect_format, err, format_g12, out, text_mode, unknown_option, Args};
use crate::help;
use crate::menu::menu_path_of;

fn list_modifiers() {
    let mut text = String::new();
    let mut category = String::new();
    for t in modifier_types() {
        if category != t.category {
            category = t.category.clone();
            text += &format!("\n{category}\n");
        }
        text += &format!("  {:<20}{}\n", t.id, t.title);
        for p in &t.parameters {
            text += &format!("      {:<16}{:<12}default {}", p.name, p.kind, p.value);
            if !p.options.is_empty() {
                text += &format!("  ({})", p.options);
            }
            text += "\n";
        }
    }
    text += "\n";
    out(&text);
}

fn list_json() -> Result<String, String> {
    let mut types = Vec::new();
    for t in modifier_types() {
        let mut entry = Json::object();
        entry.set("id", t.id.as_str());
        entry.set("title", t.title.as_str());
        entry.set("category", t.category.as_str());
        entry.set("help", t.help.as_str());
        entry.set("menu", menu_path_of(&t.id));
        let defaults = make_modifier(&t.id)?.parameters;
        let mut parameters = Vec::new();
        for p in &t.parameters {
            let mut item = Json::object();
            item.set("name", p.name.as_str());
            item.set("label", p.label.as_str());
            item.set("kind", p.kind.as_str());
            item.set("default", defaults.at(&p.name)?.clone());
            item.set("options", p.options.as_str());
            parameters.push(item);
        }
        entry.set("parameters", Json::Array(parameters));
        types.push(entry);
    }
    Ok(format!("{}\n", Json::Array(types).dump(-1)))
}

fn read_input(path: &str) -> Result<Structure, String> {
    if path == "-" {
        let mut bytes = Vec::new();
        std::io::stdin().read_to_end(&mut bytes).map_err(|e| e.to_string())?;
        let text = String::from_utf8_lossy(&bytes).replace("\r\n", "\n");
        let mut frames = read_extxyz_frames(&mut LineReader::new(&text), 1)?;
        if frames.is_empty() {
            return Err("No structure on the standard input".into());
        }
        return Ok(frames.swap_remove(0).structure);
    }
    load_structure_from_file(path).map_err(|e| format!("Cannot read {path}: {e}"))
}

/// Extended XYZ with the pipeline's selection and per-atom property.
fn extxyz(data: &PipelineData) -> String {
    let s = &data.structure;
    let property = s.atom_property.len() == s.atoms.len() && !s.atoms.is_empty();
    // Rounding residue such as -2.5e-32 is written as 0.
    let clean = |v: f64| format_g12(if v.abs() < 5e-11 { 0.0 } else { v });
    let mut text = format!("{}\n", s.atoms.len());
    if s.has_unit_cell {
        let values: Vec<String> = s.cell_vectors.iter().flatten().map(|&v| clean(v)).collect();
        text += &format!("Lattice=\"{}\" pbc=\"T T T\" ", values.join(" "));
    }
    text += &format!("Properties=species:S:1:pos:R:3:selection:I:1{}\n", if property { ":property:R:1" } else { "" });
    for (i, a) in s.atoms.iter().enumerate() {
        text += &format!("{} {} {} {} {}", a.symbol, clean(a.x), clean(a.y), clean(a.z), data.selected[i]);
        if property {
            let v = s.atom_property[i];
            text += &format!(" {}", format_g12(if v.is_finite() { v } else { 0.0 }));
        }
        text += "\n";
    }
    text
}

pub fn run(args: &Args) -> i32 {
    if args.has("--help") || args.has("-h") {
        out(help::PIPE);
        return 0;
    }
    let problem = unknown_option(&args.0, help::PIPE, "--pipe", "AtomForge --pipe --help");
    if !problem.is_empty() {
        err(&format!("Error: {problem}\n"));
        return 1;
    }
    if args.has("--list") {
        list_modifiers();
        return 0;
    }
    if args.has("--list-json") {
        return match list_json() {
            Ok(text) => {
                out(&text);
                0
            }
            Err(e) => {
                err(&format!("Error: {e}\n"));
                1
            }
        };
    }
    match pipe(args) {
        Ok(code) => code,
        Err(e) => {
            err(&format!("Error: {e}\n"));
            1
        }
    }
}

fn pipe(args: &Args) -> Result<i32, String> {
    let pipeline = match (args.find("--pipe"), args.find("--pipeline")) {
        (Some(text), _) if text != "--list" => Pipeline::parse(text)?,
        (_, Some(file)) => {
            let bytes = std::fs::read(file).map_err(|_| format!("Cannot read {file}"))?;
            let body = String::from_utf8_lossy(&bytes).replace("\r\n", "\n");
            if body.trim_start_matches([' ', '\t', '\r', '\n']).starts_with('{') {
                Pipeline::from_json(&Json::parse(&body)?)?
            } else {
                Pipeline::parse(&body)?
            }
        }
        _ => return Err("--pipe \"STEPS\" or --pipeline FILE is required (see --pipe --help)".into()),
    };
    let (Some(input), Some(output)) = (args.find("--input"), args.find("--output")) else {
        return Err("--input FILE|- and --output FILE|- are required".into());
    };
    let quiet = args.has("--quiet");
    if let Some(save) = args.find("--save-pipeline") {
        std::fs::write(save, text_mode(&format!("{}\n", pipeline.to_json().dump(2)))).map_err(|_| format!("Cannot write {save}"))?;
    }
    let source = read_input(input)?;
    let (data, stages) = pipeline.evaluate(&source)?;
    if let Some(report) = args.find("--report") {
        // Per-step results as JSON (atoms, selection, time, error, notes).
        let mut list = Vec::new();
        for (i, stage) in stages.iter().enumerate() {
            let mut item = Json::object();
            item.set("type", pipeline.modifiers[i].type_id.as_str());
            item.set("atoms", stage.atoms);
            item.set("selected", stage.selected);
            item.set("milliseconds", stage.milliseconds);
            item.set("skipped", stage.skipped);
            item.set("error", stage.error.as_str());
            item.set("notes", Json::Array(stage.notes.iter().map(|n| Json::from(n.as_str())).collect()));
            list.push(item);
        }
        let mut document = Json::object();
        document.set("stages", Json::Array(list));
        document.set("property_name", data.structure.atom_property_name.as_str());
        let _ = std::fs::write(report, text_mode(&format!("{}\n", document.dump(2))));
    }
    // The report goes to stderr so that stdout can carry the structure.
    for (i, stage) in stages.iter().enumerate() {
        let id = &pipeline.modifiers[i].type_id;
        if !stage.error.is_empty() {
            err(&format!("Error: step {} ({id}): {}\n", i + 1, stage.error));
            return Ok(1);
        }
        if quiet {
            continue;
        }
        let mut line = format!("  {} {id:<20}", i + 1);
        if stage.skipped {
            line += "skipped (disabled)\n";
        } else {
            line += &format!("{} atoms, {} selected ({:.1} ms)\n", stage.atoms, stage.selected, stage.milliseconds);
        }
        for note in &stage.notes {
            line += &format!("      {note}\n");
        }
        err(&line);
    }
    if output == "-" {
        out(&extxyz(&data));
    } else {
        if data.structure.atoms.is_empty() {
            return Err("The pipeline removed every atom; nothing to write".into());
        }
        save_structure(&data.structure, output, &detect_format(output)).map_err(|_| format!("Cannot write {output}"))?;
        if !quiet {
            err(&format!("Saved {} atoms to {output}\n", data.structure.atoms.len()));
        }
    }
    Ok(0)
}
