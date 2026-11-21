//! The python-script step (the C++ pipeline/PythonScriptStep): runs a user's
//! Python script on the structure, as OVITO's Python script modifier. The
//! script is called as  python SCRIPT INPUT.json OUTPUT.json  and writes the
//! modified data: {"structure": {"symbols", "positions", "cell"}, "pbc",
//! "selected", "property", "property_name"}; atomforge.pipeline.script_step
//! does the reading and writing.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Stdio};

use af_core::elements::default_colour;
use af_core::json::Json;
use af_science::data::{periodic_axes, structure_from_json, structure_json};

use crate::pipeline::{ModifierParameter, ModifierType, PipelineData};
use crate::support::text;

const ID: &str = "python-script";

fn scratch_directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    std::env::temp_dir().join(format!("AtomForge-python-{}-{stamp}", std::process::id()))
}

fn python_script(p: &Json, data: &mut PipelineData) -> Result<(), String> {
    let script = text(p, ID, "script")?;
    if script.is_empty() || !std::path::Path::new(&script).is_file() {
        return Err("Choose an existing Python script file".into());
    }
    let directory = scratch_directory();
    fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
    let (input, output, log) = (directory.join("input.json"), directory.join("output.json"), directory.join("log.txt"));
    let s = &data.structure;
    let mut document = Json::object();
    document.set("structure", structure_json(s));
    let pbc = periodic_axes(s);
    document.set("pbc", Json::Array(pbc.iter().map(|&b| b.into()).collect()));
    document.set("selected", Json::Array(data.selected.iter().map(|&c| Json::from(if c != 0 { 1.0 } else { 0.0 })).collect()));
    if s.atom_property.len() == s.atoms.len() && !s.atoms.is_empty() {
        document.set("property", Json::Array(s.atom_property.iter().map(|&v| Json::from(v)).collect()));
        document.set("property_name", s.atom_property_name.as_str());
    }
    document.set("arguments", text(p, ID, "arguments")?);
    fs::write(&input, document.dump(-1)).map_err(|e| e.to_string())?;
    let mut interpreter = text(p, ID, "python")?;
    if interpreter.is_empty() {
        interpreter = match std::env::var("ATOMFORGE_PYTHON") {
            Ok(v) if !v.is_empty() => v,
            _ => if cfg!(windows) { "python" } else { "python3" }.to_string(),
        };
    }
    let status = (|| -> std::io::Result<bool> {
        let file = fs::File::create(&log)?;
        let status = Command::new(&interpreter)
            .arg(&script)
            .arg(&input)
            .arg(&output)
            .stdin(Stdio::null())
            .stdout(file.try_clone()?)
            .stderr(file)
            .status()?;
        Ok(status.success())
    })();
    // Read as the C++ text-mode stream does (CRLF to LF).
    let mut messages = fs::read_to_string(&log).unwrap_or_default().replace("\r\n", "\n");
    let ok = match status {
        Ok(ok) => ok,
        Err(e) => {
            messages.push_str(&e.to_string());
            false
        }
    };
    if !ok || !output.exists() {
        let _ = fs::remove_dir_all(&directory);
        let tail = if messages.len() > 1500 {
            let mut start = messages.len() - 1500;
            while !messages.is_char_boundary(start) {
                start += 1;
            }
            &messages[start..]
        } else {
            &messages[..]
        };
        return Err(format!("The Python script failed ({interpreter}):\n{tail}"));
    }
    let result_text = fs::read_to_string(&output).map_err(|e| e.to_string());
    let _ = fs::remove_dir_all(&directory);
    let result = Json::parse(&result_text?)?;
    if let Some(structure) = result.find("structure") {
        let s = &data.structure;
        let mut changed = structure_from_json(structure)?;
        changed.pbc = s.pbc;
        // Atoms keep their colours where they still exist.
        for (i, atom) in changed.atoms.iter_mut().enumerate() {
            if i < s.atoms.len() && s.atoms[i].symbol == atom.symbol {
                (atom.r, atom.g, atom.b) = (s.atoms[i].r, s.atoms[i].g, s.atoms[i].b);
                atom.fixed = s.atoms[i].fixed;
            } else {
                let c = default_colour(atom.atomic_number);
                (atom.r, atom.g, atom.b) = (c[0], c[1], c[2]);
            }
        }
        if changed.atoms.len() == s.atoms.len() {
            changed.grain_colors = s.grain_colors.clone();
            changed.grain_region_ids = s.grain_region_ids.clone();
            changed.atom_property = s.atom_property.clone();
            changed.atom_property_name = s.atom_property_name.clone();
            changed.frozen_property = s.frozen_property.clone();
            changed.frozen_property_name = s.frozen_property_name.clone();
        } else {
            data.selected = vec![0; changed.atoms.len()];
        }
        changed.cell_offset = s.cell_offset;
        data.structure = changed;
    }
    if let Some(pbc) = result.find("pbc") {
        if let Ok(items) = pbc.items() {
            if items.len() == 3 {
                for (k, v) in items.iter().enumerate() {
                    data.structure.pbc[k] = v.boolean()?;
                }
            }
        }
    }
    let n = data.structure.atoms.len();
    if let Some(mask) = result.find("selected") {
        if let Ok(items) = mask.items() {
            if items.len() != n {
                return Err(format!("The script's selection has {} entries for {n} atoms", items.len()));
            }
            data.selected = items.iter().map(|v| u8::from(if v.is_bool() { v.boolean().unwrap_or(false) } else { v.number().unwrap_or(0.0) != 0.0 })).collect();
        }
    }
    if let Some(values) = result.find("property") {
        if let Ok(items) = values.items() {
            if items.len() != n {
                return Err(format!("The script's property has {} values for {n} atoms", items.len()));
            }
            data.structure.atom_property = items.iter().map(|v| if v.is_number() { v.number().unwrap_or(f64::NAN) } else { f64::NAN }).collect();
            if let Some(name) = result.find("property_name") {
                if name.is_string() {
                    data.structure.atom_property_name = name.string()?.to_string();
                }
            }
        }
    }
    for line in messages.split('\n') {
        if !line.is_empty() {
            data.notes.push(line.to_string());
        }
    }
    Ok(())
}

pub fn types() -> Vec<ModifierType> {
    vec![ModifierType::new(
        ID,
        "Python script",
        "Modification",
        "Runs a Python script on the structure: python SCRIPT INPUT.json OUTPUT.json. With atomforge.pipeline.script_step the script is a function of the\nstructure (and its selection and property) returning the changed data; what it prints becomes the step's notes.",
        vec![
            ModifierParameter::new("script", "Python script file", "string", "", ""),
            ModifierParameter::new("arguments", "Arguments for the script (text)", "string", "", ""),
            ModifierParameter::new("python", "Python interpreter (empty: ATOMFORGE_PYTHON or python)", "string", "", ""),
        ],
        python_script,
    )]
}
