//! `--science TOOL` and `--science-batch` (the C++ cli/ScienceCLI): catalog
//! tools on JSON requests, batches of runs, reports and generated files.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use af_core::json::Json;
use af_core::structure::Structure;
use af_io::loader::load_structure_from_file;
use af_science::batch::{batch_csv, lexically_normal, run_batch};
use af_science::catalog::catalog;
use af_science::data::write_extxyz;
use af_science::tools::{result_document, result_report, run_tool};

use crate::args::{err, g, out, Args};

fn print_help() {
    let mut text = String::from(concat!(
        "AtomForge scientific tools (native, no Python required)\n",
        "\n",
        "Usage:\n",
        "  AtomForge --science <tool> --input REQUEST.json --output RESULT.json\n",
        "            [--report REPORT.txt] [--structures FRAMES.extxyz] [--files DIR] [--overwrite]\n",
        "  AtomForge --science --catalog      Print tools, parameters and defaults as JSON\n",
        "  AtomForge --science-batch BATCH.json --output RESULTS.json [--csv TABLE.csv]\n",
        "            BATCH.json: {\"tool\": ..., \"base\": {...}, \"sweep\": {\"param\": [values]},\n",
        "                         \"files\": {\"param\": \"data/*.xyz\"}, \"collect\": [\"result.key\"]}\n",
        "            Runs the Cartesian product of sweeps and file matches; exit 2 if a run failed.\n",
        "\n",
        "REQUEST.json is an object of tool parameters. Data values are inline JSON\n",
        "arrays or {\"file\": PATH, \"field\": NAME, \"column\": N} references to JSON,\n",
        "NPY, CSV/TXT/DAT, EIGENVAL, XYZ/extXYZ, POSCAR/XDATCAR, LAMMPS dump or other\n",
        "structure files; paths are relative to REQUEST.json.\n",
        "Simulation tools take an interatomic potential such as {\"potential\": \"EMT\"}\n",
        "or {\"potential\": \"LennardJones\", \"epsilon\": 0.0104, \"sigma\": 3.4, \"cutoff\": 8.5}\n",
        "or {\"potential\": \"EAM\", \"file\": \"Cu.eam.alloy\"} (setfl, eam/fs or funcfl tables).\n",
        "\n",
        "Tools:\n",
    ));
    for tool in catalog() {
        text += &format!("  {}{}{}\n", tool.id, " ".repeat(26 - tool.id.len().min(25)), tool.title);
    }
    text += "\n";
    out(&text);
}

fn catalog_json() -> Result<Json, String> {
    let mut tools = Vec::new();
    for (priority, tool) in catalog().iter().enumerate() {
        let mut entry = Json::object();
        entry.set("priority", priority + 1);
        entry.set("id", tool.id);
        entry.set("name", tool.title);
        entry.set("category", tool.category);
        entry.set("description", tool.help);
        let mut parameters = Vec::new();
        for p in tool.parameters {
            let mut item = Json::object();
            item.set("name", p.name);
            item.set("label", p.label);
            item.set("kind", p.kind);
            item.set("required", p.required);
            item.set("default", if p.value.is_empty() { Json::Null } else { Json::parse(p.value)? });
            parameters.push(item);
        }
        entry.set("parameters", Json::Array(parameters));
        tools.push(entry);
    }
    Ok(Json::Array(tools))
}

/// Writes beside the destination first so a failure never truncates an existing file (binary).
fn write_text(path: &Path, text: &str) -> Result<(), String> {
    let parent = path.parent().filter(|p| !p.as_os_str().is_empty());
    if let Some(parent) = parent {
        std::fs::create_dir_all(parent).map_err(|_| format!("Cannot write {}", path.display()))?;
    }
    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let temporary = parent.map(Path::to_path_buf).unwrap_or_default().join(format!(".atomforge-{name}.tmp"));
    std::fs::write(&temporary, text).map_err(|_| format!("Cannot write {}", path.display()))?;
    let _ = std::fs::remove_file(path);
    std::fs::rename(&temporary, path).map_err(|e| e.to_string())
}

fn absolute(path: &Path) -> PathBuf {
    lexically_normal(&std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf()))
}

fn reader(path: &Path) -> Result<Structure, String> {
    load_structure_from_file(&path.to_string_lossy())
}

pub fn run(argv: &Args) -> i32 {
    let mut args: BTreeMap<String, String> = BTreeMap::new();
    let (mut tool, mut batch_file) = (String::new(), String::new());
    let a = &argv.0;
    let mut i = 0;
    while i < a.len() {
        let key = a[i].as_str();
        match key {
            "--science" => {
                if i + 1 < a.len() && !a[i + 1].starts_with('-') {
                    i += 1;
                    tool = a[i].clone();
                }
            }
            "--science-batch" => {
                if i + 1 >= a.len() {
                    err("Error: --science-batch requires a batch file\n");
                    return 1;
                }
                i += 1;
                batch_file = a[i].clone();
            }
            "--csv" => {
                if i + 1 >= a.len() {
                    err("Error: --csv requires a value\n");
                    return 1;
                }
                i += 1;
                args.insert(key.into(), a[i].clone());
            }
            "--catalog" | "--overwrite" | "--help" | "-h" => {
                args.insert(key.into(), String::new());
            }
            "--input" | "--output" | "--report" | "--structures" | "--files" => {
                if i + 1 >= a.len() {
                    err(&format!("Error: {key} requires a value\n"));
                    return 1;
                }
                i += 1;
                args.insert(key.into(), a[i].clone());
            }
            _ => {
                err(&format!("Error: unknown science option {key}\n"));
                return 1;
            }
        }
        i += 1;
    }
    match science(&args, &tool, &batch_file) {
        Ok(code) => code,
        Err(e) => {
            err(&format!("AtomForge science: {e}\n"));
            1
        }
    }
}

fn science(args: &BTreeMap<String, String>, tool: &str, batch_file: &str) -> Result<i32, String> {
    if args.contains_key("--catalog") {
        out(&format!("{}\n", catalog_json()?.dump(2)));
        return Ok(0);
    }
    let structure_reader: &dyn Fn(&Path) -> Result<Structure, String> = &reader;
    if !batch_file.is_empty() {
        let Some(output) = args.get("--output") else { return Err("--output is required for a batch".into()) };
        let input = PathBuf::from(batch_file);
        let text = std::fs::read(&input).map_err(|_| format!("Cannot open {}", input.display()))?;
        let output = PathBuf::from(output);
        if output.exists() && !args.contains_key("--overwrite") {
            return Err(format!("{} exists; use --overwrite", output.display()));
        }
        let base = std::path::absolute(&input).unwrap_or(input.clone()).parent().map(Path::to_path_buf).unwrap_or_default();
        let result = run_batch(&Json::parse(&String::from_utf8_lossy(&text))?, &base, Some(structure_reader))?;
        write_text(&output, &format!("{}\n", result.dump(2)))?;
        if let Some(csv) = args.get("--csv") {
            write_text(Path::new(csv), &batch_csv(&result)?)?;
        }
        let (runs, failures) = (result.at("run_count")?.number()?, result.at("failures")?.number()?);
        out(&format!("Batch: {} runs, {} failed; saved {}\n", g(runs), g(failures), output.display()));
        return Ok(if failures > 0.0 { 2 } else { 0 });
    }
    let help = args.contains_key("--help") || args.contains_key("-h");
    if help || tool.is_empty() {
        print_help();
        return Ok(if tool.is_empty() && !help { 1 } else { 0 });
    }
    let (Some(input), Some(output)) = (args.get("--input"), args.get("--output")) else {
        return Err("--input and --output are required".into());
    };
    let (input, output) = (PathBuf::from(input), PathBuf::from(output));
    let overwrite = args.contains_key("--overwrite");
    let mut destinations = vec![output.clone()];
    if let Some(report) = args.get("--report") {
        destinations.push(PathBuf::from(report));
    }
    if let Some(structures) = args.get("--structures") {
        destinations.push(PathBuf::from(structures));
    }
    for i in 0..destinations.len() {
        if destinations[i].exists() && !overwrite {
            return Err(format!("{} exists; choose another file or use --overwrite", destinations[i].display()));
        }
        for j in 0..i {
            if absolute(&destinations[i]) == absolute(&destinations[j]) {
                return Err("Result, report and structure files must have distinct paths".into());
            }
        }
        if absolute(&destinations[i]) == absolute(&input) {
            return Err("Input and output files must differ".into());
        }
    }
    let text = std::fs::read(&input).map_err(|_| format!("Cannot open {}", input.display()))?;
    let request = Json::parse(&String::from_utf8_lossy(&text))?;
    let base = std::path::absolute(&input).unwrap_or(input.clone()).parent().map(Path::to_path_buf).unwrap_or_default();
    let result = run_tool(tool, &request, &base, Some(structure_reader))?;
    if let Some(structures) = args.get("--structures") {
        if !result.frames.is_empty() {
            write_extxyz(Path::new(structures), &result.frames, &result.velocities, &result.times)?;
        }
    }
    write_text(&output, &format!("{}\n", result_document(tool, &request, &result).dump(2)))?;
    if let Some(report) = args.get("--report") {
        write_text(Path::new(report), &result_report(tool, &result.result)?)?;
    }
    // Generated input files (e.g. KPOINTS, POSCAR) go to the --files directory.
    if let Some(directory) = args.get("--files") {
        if let Some(files) = result.result.find("files") {
            if files.is_object() {
                for (name, content) in files.members()? {
                    let target = Path::new(directory).join(name);
                    if target.exists() && !overwrite {
                        return Err(format!("{} exists; use --overwrite", target.display()));
                    }
                    write_text(&target, content.string()?)?;
                }
            }
        }
    }
    out(&format!("Saved {}\n", output.display()));
    Ok(0)
}
