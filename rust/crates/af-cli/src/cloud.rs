//! `--cloud build|generate|info` (the C++ cli/CloudCLI): out-of-core atom
//! clouds of up to billions of atoms.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicI32, Ordering};
use std::time::Instant;

use af_cloud::{build_from_lammps_dump, build_from_xyz, generate_crystal, BuildOptions, CloudBuilder, CloudFile, CloudInfo, SourceAtom};
use af_core::format::stod;

use crate::args::{err, g, load, out, Args};
use crate::help;

fn bytes(value: f64) -> String {
    if value >= 1e9 {
        format!("{:.2} GB", value / 1e9)
    } else if value >= 1e6 {
        format!("{:.1} MB", value / 1e6)
    } else {
        format!("{:.1} kB", value / 1e3)
    }
}

fn print_info(info: &CloudInfo, path: &Path) {
    let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(u64::MAX);
    let mut text = format!("{}\n  atoms    {}\n  chunks   {}\n  species ", path.display(), info.atoms, info.chunks.len());
    for s in &info.species {
        text += &format!(" {}", s.name);
    }
    let (l, u) = (info.lower, info.upper);
    text += &format!(
        "\n  box      [{}, {}, {}] - [{}, {}, {}] Angstrom\n  size     {}\n",
        g(l[0]),
        g(l[1]),
        g(l[2]),
        g(u[0]),
        g(u[1]),
        g(u[2]),
        bytes(size as f64)
    );
    if !info.source.is_empty() {
        text += &format!("  source   {}\n", info.source);
    }
    out(&text);
}

/// std::stoull / std::stoul: a leading unsigned integer, else the function's name.
fn stoull(text: &str, name: &str) -> Result<u64, String> {
    let t = text.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let t = t.strip_prefix('+').unwrap_or(t);
    let digits: String = t.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().map_err(|_| name.to_string())
}

fn stod_prefix(text: &str) -> Result<f64, String> {
    stod(text).map(|(v, _)| v).map_err(|_| "stod".to_string())
}

pub fn run(args: &Args) -> i32 {
    let command = args.find("--cloud");
    if command.is_none() || args.has("--help") || args.has("-h") {
        out(help::CLOUD);
        return if command.is_some() { 0 } else { 1 };
    }
    match cloud(args, command.unwrap_or_default()) {
        Ok(()) => 0,
        Err(e) => {
            err(&format!("\nError: {e}\n"));
            1
        }
    }
}

fn cloud(args: &Args, mode: &str) -> Result<(), String> {
    let mut options = BuildOptions::default();
    if let Some(chunk) = args.find("--chunk") {
        options.chunk_atoms = stoull(chunk, "stoull")?;
    }
    if let Some(temp) = args.find("--temp") {
        options.temporary = Some(PathBuf::from(temp));
    }
    if let Some(threads) = args.find("--threads") {
        options.threads = stoull(threads, "stoul")? as u32;
    }
    let quiet = args.has("--quiet");
    let start = Instant::now();
    let last_percent = AtomicI32::new(-1);
    let report = |fraction: f64, message: &str| -> bool {
        if quiet || fraction < 0.0 {
            return true;
        }
        let percent = (fraction * 100.0) as i32;
        if last_percent.swap(percent, Ordering::Relaxed) != percent {
            let mut stderr = std::io::stderr().lock();
            let _ = write!(stderr, "\r  {percent:>3}%  {message}          ");
            let _ = stderr.flush();
        }
        true
    };
    let progress: af_cloud::Progress = Some(&report);
    let finished = |info: CloudInfo, output: &Path| {
        let seconds = start.elapsed().as_secs_f64();
        if !quiet {
            err(&(format!(
                "\r  100%  done in {seconds:.1} s ({:.1} M atoms/s)          \n",
                info.atoms as f64 / seconds.max(1e-9) / 1e6
            )));
        }
        print_info(&info, output);
    };
    if mode == "info" {
        // The file is the argument after "info", or --input.
        let mut file = args.find("--input").map(str::to_string);
        if file.is_none() {
            let a = &args.0;
            file = (0..a.len().saturating_sub(1)).find(|&i| a[i] == "info").map(|i| a[i + 1].clone());
        }
        let Some(file) = file else { return Err("--cloud info FILE.afcloud".into()) };
        let path = PathBuf::from(file);
        let cloud = CloudFile::open(&path)?;
        print_info(cloud.info(), &path);
        return Ok(());
    }
    let Some(output) = args.find("--output") else { return Err("--output OUT.afcloud is required".into()) };
    let output = PathBuf::from(output);
    match mode {
        "generate" => {
            let (Some(lattice), Some(a), Some(element)) = (args.find("--lattice"), args.find("--a"), args.find("--element")) else {
                return Err("--lattice, --a and --element are required".into());
            };
            let mut cells = [1u64; 3];
            let list = &args.0;
            for i in 0..list.len() {
                if i + 3 < list.len() && list[i] == "--cells" {
                    for k in 0..3 {
                        cells[k] = stoull(&list[i + 1 + k], "stoull")?;
                    }
                }
            }
            let c = match args.find("--c") {
                Some(c) => stod_prefix(c)?,
                None => 0.0,
            };
            let a = stod_prefix(a)?;
            let info = generate_crystal(&output, lattice, a, c, cells, element, args.find("--second").unwrap_or(""), progress, options)?;
            finished(info, &output);
            Ok(())
        }
        "build" => {
            let Some(input) = args.find("--input") else { return Err("--input FILE is required".into()) };
            let input = PathBuf::from(input);
            let mut types = Vec::new();
            if let Some(names) = args.find("--types") {
                let mut pieces: Vec<&str> = names.split(',').collect();
                if pieces.last() == Some(&"") {
                    pieces.pop();
                }
                types = pieces.into_iter().map(str::to_string).collect();
            }
            let frame = match args.find("--frame") {
                Some(f) => stoull(f, "stoull")? as usize,
                None => 0,
            };
            // Streamed formats by content, others through the structure loaders.
            let first = std::fs::read(&input).map(|b| String::from_utf8_lossy(&b[..b.len().min(64)]).lines().next().unwrap_or("").to_string()).unwrap_or_default();
            let extension = input.extension().map(|e| format!(".{}", e.to_string_lossy().to_ascii_lowercase())).unwrap_or_default();
            let info = if first.starts_with("ITEM: TIMESTEP") {
                build_from_lammps_dump(&input, &output, &types, frame, progress, options)?
            } else if extension == ".xyz" || extension == ".extxyz" {
                build_from_xyz(&input, &output, frame, progress, options)?
            } else {
                let name = input.to_string_lossy().into_owned();
                let structure = load(&name).map_err(|e| format!("Cannot read {name}: {e}"))?;
                let (mut lower, mut upper) = ([1e300f64; 3], [-1e300f64; 3]);
                for atom in &structure.atoms {
                    let p = [atom.x, atom.y, atom.z];
                    for k in 0..3 {
                        lower[k] = lower[k].min(p[k]);
                        upper[k] = upper[k].max(p[k]);
                    }
                }
                if structure.atoms.is_empty() {
                    return Err(format!("{name} has no atoms"));
                }
                let mut builder = CloudBuilder::new(&output, lower, upper, structure.atoms.len() as u64, options)?;
                builder.set_source(&input.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default());
                if structure.has_unit_cell {
                    builder.set_cell(structure.cell_vectors, structure.cell_offset);
                }
                for atom in &structure.atoms {
                    let species = builder.species(&atom.symbol)?;
                    builder.add(&SourceAtom { x: atom.x, y: atom.y, z: atom.z, species })?;
                }
                builder.finish(progress)?
            };
            finished(info, &output);
            Ok(())
        }
        _ => Err(format!("Unknown --cloud command {mode} (build, generate or info)")),
    }
}
