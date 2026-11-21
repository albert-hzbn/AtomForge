//! Build and Edit operations as structure-pipeline steps (the C++
//! cli/BuildModifiers): each `--build` mode (and the cell sculptor) runs on
//! the pipeline's structure with its command-line options.

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};

use af_core::json::Json;
use af_pipeline::pipeline::{add_family, ModifierParameter, ModifierType, PipelineData};

use crate::args::{captured, load, save, unknown_option, Args};
use crate::help;

#[derive(Clone, Copy, PartialEq)]
enum Input {
    None,
    Required,
    Optional,
    LayerA,
}

struct Mode {
    /// --build mode (or "sculpt" for --analyze sculpt).
    id: &'static str,
    title: &'static str,
    /// "Build" (makes a structure) or "Build & edit" (changes it).
    category: &'static str,
    input: Input,
    /// Default option flags.
    defaults: &'static str,
    help: &'static str,
}

const MODES: [Mode; 18] = [
    Mode { id: "bulk", title: "Bulk crystal", category: "Build", input: Input::None, defaults: "--system cubic --spacegroup 225 --a 3.615 --atom \"Cu 0 0 0\"", help: help::BULK },
    Mode { id: "custom", title: "Custom structure", category: "Build", input: Input::Optional, defaults: "", help: help::CUSTOM },
    Mode { id: "gb", title: "Grain boundary", category: "Build", input: Input::Optional, defaults: "", help: help::GB },
    Mode { id: "poly", title: "Polycrystal", category: "Build", input: Input::Optional, defaults: "", help: help::POLY },
    Mode { id: "nano", title: "Nanocrystal", category: "Build", input: Input::Optional, defaults: "", help: help::NANO },
    Mode { id: "amorphous", title: "Amorphous structure", category: "Build", input: Input::None, defaults: "--element \"Si 64\" --density 2.33", help: help::AMORPHOUS },
    Mode { id: "interface", title: "Heterogeneous interface", category: "Build", input: Input::LayerA, defaults: "", help: help::INTERFACE },
    Mode { id: "nanowire", title: "Nanowire", category: "Build", input: Input::Required, defaults: "--axis 2 --radius 10", help: help::NANOWIRE },
    Mode {
        id: "core-shell",
        title: "Core-shell particle",
        category: "Build & edit",
        input: Input::Required,
        defaults: "--core-radius 8 --core-element Cu --shell-element Ag",
        help: help::CORESHELL,
    },
    Mode { id: "surface", title: "Surface slab", category: "Build & edit", input: Input::Required, defaults: "--h 1 --k 1 --l 1 --layers 4 --vacuum 15", help: help::SURFACE },
    Mode { id: "vacancy", title: "Vacancies", category: "Build & edit", input: Input::Required, defaults: "--count 1 --seed 1", help: help::VACANCY },
    Mode {
        id: "sss",
        title: "Substitutional solid solution",
        category: "Build & edit",
        input: Input::Required,
        defaults: "--frac \"Cu=0.7,Zn=0.3\" --seed 42",
        help: help::SSS,
    },
    Mode {
        id: "sqs",
        title: "Special quasirandom structure",
        category: "Build & edit",
        input: Input::Required,
        defaults: "--element \"Cu 0.5\" --element \"Au 0.5\"",
        help: help::SQS,
    },
    Mode { id: "strain", title: "Strain the cell", category: "Build & edit", input: Input::Required, defaults: "--exx 0.01", help: help::STRAIN },
    Mode { id: "primitive", title: "Primitive cell", category: "Build & edit", input: Input::Required, defaults: "", help: help::PRIMITIVE },
    Mode { id: "dislocation", title: "Insert dislocation", category: "Build & edit", input: Input::Optional, defaults: "", help: help::DISLOCATION },
    Mode { id: "stacking-fault", title: "Stacking fault", category: "Build & edit", input: Input::Required, defaults: "--plane 1 --layers 9", help: help::STACKINGFAULT },
    Mode {
        id: "sculpt",
        title: "Cell sculptor",
        category: "Build & edit",
        input: Input::Required,
        defaults: "--slabs \"1 0 0 0 1 1;0 1 0 0 1 1;0 0 1 0 1 1\"",
        help: help::SCULPT,
    },
];

/// Shell-style tokens (quotes group words and are removed).
pub fn tokens(text: &str) -> Result<Vec<String>, String> {
    let mut out = Vec::new();
    let mut current = String::new();
    let (mut in_token, mut quote) = (false, None::<char>);
    for c in text.chars() {
        if let Some(q) = quote {
            if c == q {
                quote = None;
            } else {
                current.push(c);
            }
            continue;
        }
        if c == '"' || c == '\'' {
            quote = Some(c);
            in_token = true;
            continue;
        }
        if c.is_ascii() && crate::args::is_space(c as u8) {
            if in_token {
                out.push(std::mem::take(&mut current));
                in_token = false;
            }
            continue;
        }
        current.push(c);
        in_token = true;
    }
    if quote.is_some() {
        return Err("Unterminated quote in the options".into());
    }
    if in_token {
        out.push(current);
    }
    Ok(out)
}

fn run_mode(mode: &Mode, parameters: &Json, data: &mut PipelineData) -> Result<(), String> {
    static COUNTER: AtomicU32 = AtomicU32::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    let folder: PathBuf = std::env::temp_dir().join("AtomForge-pipeline").join(format!("{stamp}-{}", COUNTER.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir_all(&folder).map_err(|e| e.to_string())?;
    let result = run_in(mode, parameters, data, &folder);
    let _ = std::fs::remove_dir_all(&folder);
    result
}

fn run_in(mode: &Mode, parameters: &Json, data: &mut PipelineData, folder: &std::path::Path) -> Result<(), String> {
    let sculpt = mode.id == "sculpt";
    let mut args: Vec<String> = vec!["AtomForge".into()];
    if sculpt {
        args.extend(["--analyze".into(), "sculpt".into()]);
    } else {
        args.extend(["--build".into(), mode.id.into()]);
    }
    let options = match parameters.find("options") {
        Some(v) => v.string()?.to_string(),
        None => mode.defaults.to_string(),
    };
    args.extend(tokens(&options)?);
    // The pipeline's structure is the input (extended XYZ keeps the cell, if
    // any, the atom order and full precision).
    let have_atoms = !data.structure.atoms.is_empty();
    let use_input = match parameters.find("use_input") {
        Some(v) => v.boolean()?,
        None => true,
    };
    let wants = mode.input == Input::Required || mode.input == Input::LayerA || (mode.input == Input::Optional && use_input);
    if wants && have_atoms {
        let flag = if mode.input == Input::LayerA { "--layerA" } else { "--input" };
        if !args.iter().any(|a| a == flag) {
            let input = folder.join("input.xyz");
            let path = input.to_string_lossy().into_owned();
            af_io::loader::save_structure(&data.structure, &path, "xyz").map_err(|e| format!("Cannot pass the structure to {}: {e}", mode.id))?;
            args.extend([flag.to_string(), path]);
        }
    } else if mode.input == Input::Required {
        return Err(format!("{} needs a structure from earlier steps", mode.title));
    }
    // A build may make a structure without a cell (core-shell, a finite nanocrystal).
    let output = folder.join(if sculpt { "output.vasp" } else { "output.xyz" });
    let output_path = output.to_string_lossy().into_owned();
    args.extend(["--output".into(), output_path.clone()]);
    if !sculpt {
        let problem = unknown_option(&args, mode.help, mode.id, "");
        if !problem.is_empty() {
            return Err(problem);
        }
    }
    // The mode reports on stdout/stderr; keep that out of the pipeline's own output.
    let argv = Args(args[1..].to_vec());
    let (status, out, err) = captured(|| -> Result<i32, String> {
        if sculpt {
            Ok(crate::analysis::run(&argv))
        } else {
            let entry = crate::BUILD_MODES.iter().find(|m| m.name == mode.id).expect("registered mode");
            (entry.run)(&argv)
        }
    });
    let status = status?;
    if status != 0 || !output.exists() {
        let mut message = if err.is_empty() { out } else { err };
        while message.ends_with(|c: char| c.is_ascii_whitespace()) {
            message.pop();
        }
        if let Some(rest) = message.strip_prefix("Error: ") {
            message = rest.to_string();
        }
        return Err(if message.is_empty() { format!("{} failed", mode.title) } else { message });
    }
    let result = load(&output_path).map_err(|e| format!("Cannot read the result of {}: {e}", mode.id))?;
    data.structure = result;
    data.clear_selection();
    data.notes.push(format!("{}: {} atoms", mode.title, data.structure.atoms.len()));
    let _ = save; // (the input is written by af_io directly above)
    Ok(())
}

fn types() -> Vec<ModifierType> {
    MODES
        .iter()
        .map(|mode| {
            let lead = match mode.input {
                Input::None => "Makes a new structure (the pipeline's structure is replaced).",
                Input::Optional => "Uses the pipeline's structure as its input when 'use input' is on.",
                _ => "Takes the pipeline's structure as its input.",
            };
            let help = format!("{lead} Options are the command-line flags of AtomForge --build {}:\n\n{}", mode.id, mode.help);
            let mut parameters = vec![ModifierParameter::new("options", "Options (command-line flags)", "options", mode.defaults, "")];
            if mode.input == Input::Optional {
                parameters.push(ModifierParameter::new("use_input", "Use the pipeline's structure as input", "bool", "true", ""));
            }
            ModifierType::new(&format!("build-{}", mode.id), mode.title, mode.category, &help, parameters, move |p: &Json, data: &mut PipelineData| {
                run_mode(mode, p, data)
            })
        })
        .collect()
}

/// Registers the build-* steps (before the pipeline registry is first used).
pub fn register() {
    add_family(types);
}
