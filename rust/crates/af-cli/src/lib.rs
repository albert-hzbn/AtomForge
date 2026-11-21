//! The AtomForge command line (the C++ main.cpp CLI path and src/cli) as a
//! library (the desktop registers the same Build steps):
//! builders, analyses and conversion, scientific tools, structure pipelines,
//! atom clouds and headless rendering. Messages, help pages, exit codes and
//! output files match the C++ program.

mod analysis;
pub mod args;
mod build_crystal;
mod build_defect;
pub mod build_modifiers;
mod build_nano;
mod cloud;
mod help;
pub mod menu;
mod pipe;
mod render;
mod science;

use args::{err, out, unknown_option, Args};

const VERSION: &str = "0.4.0";

type RunFn = fn(&Args) -> Result<i32, String>;

/// A `--build` mode: its runner and help page.
pub struct BuildMode {
    pub name: &'static str,
    pub run: RunFn,
    pub help: &'static str,
}

pub const BUILD_MODES: [BuildMode; 17] = [
    BuildMode { name: "bulk", run: build_crystal::run_bulk, help: help::BULK },
    BuildMode { name: "gb", run: build_defect::run_gb, help: help::GB },
    BuildMode { name: "poly", run: build_nano::run_poly, help: help::POLY },
    BuildMode { name: "nano", run: build_nano::run_nano, help: help::NANO },
    BuildMode { name: "amorphous", run: build_nano::run_amorphous, help: help::AMORPHOUS },
    BuildMode { name: "sss", run: build_crystal::run_sss, help: help::SSS },
    BuildMode { name: "dislocation", run: build_defect::run_dislocation, help: help::DISLOCATION },
    BuildMode { name: "custom", run: build_crystal::run_custom, help: help::CUSTOM },
    BuildMode { name: "interface", run: build_defect::run_interface, help: help::INTERFACE },
    BuildMode { name: "stacking-fault", run: build_defect::run_stacking_fault, help: help::STACKINGFAULT },
    BuildMode { name: "vacancy", run: build_crystal::run_vacancy, help: help::VACANCY },
    BuildMode { name: "strain", run: build_crystal::run_strain, help: help::STRAIN },
    BuildMode { name: "primitive", run: build_crystal::run_primitive, help: help::PRIMITIVE },
    BuildMode { name: "surface", run: build_crystal::run_surface, help: help::SURFACE },
    BuildMode { name: "sqs", run: build_crystal::run_sqs, help: help::SQS },
    BuildMode { name: "nanowire", run: build_nano::run_nanowire, help: help::NANOWIRE },
    BuildMode { name: "core-shell", run: build_nano::run_core_shell, help: help::CORESHELL },
];

fn find_build_mode(name: &str) -> Option<&'static BuildMode> {
    BUILD_MODES.iter().find(|m| m.name == name)
}

fn build_mode_names() -> String {
    BUILD_MODES.iter().map(|m| m.name).collect::<Vec<_>>().join(" | ")
}

/// Whether the arguments ask for the command line (else the desktop opens).
pub fn is_cli_mode(args: &Args) -> bool {
    args.0.iter().any(|a| {
        matches!(
            a.as_str(),
            "--build" | "--analyze" | "--convert" | "--render" | "--science" | "--science-batch" | "--pipe" | "--pipeline" | "--cloud" | "--help" | "-h" | "--version" | "-v"
        )
    })
}

/// Runs a build mode (unknown options refused first); errors print as "Error: ...".
pub fn run_build_mode(mode: &BuildMode, args: &Args) -> i32 {
    let problem = unknown_option(&args.0, mode.help, mode.name, "");
    if !problem.is_empty() {
        err(&format!("Error: {problem}\n"));
        return 1;
    }
    match (mode.run)(args) {
        Ok(code) => code,
        Err(e) => {
            err(&format!("Error: {e}\n"));
            1
        }
    }
}

pub fn run_cli(args: &Args) -> i32 {
    if args.has("--analyze") || args.has("--convert") {
        return analysis::run(args);
    }
    if args.has("--render") {
        return render::run(args);
    }
    if args.has("--science") || args.has("--science-batch") {
        return science::run(args);
    }
    if args.has("--pipe") || args.has("--pipeline") {
        return pipe::run(args);
    }
    if args.has("--cloud") {
        return cloud::run(args);
    }
    if args.has("--version") || args.has("-v") {
        out(&format!("AtomForge {VERSION}\n"));
        return 0;
    }
    if args.has("--help") || args.has("-h") {
        let topic = args.find("--help").or_else(|| args.find("-h"));
        if let Some(topic) = topic {
            if let Some(mode) = find_build_mode(topic) {
                out(mode.help);
                return 0;
            }
            err(&format!("Unknown help topic '{topic}'.  Valid topics: {}\n", build_mode_names()));
            return 1;
        }
        out(help::MAIN);
        return 0;
    }
    let Some(mode) = args.find("--build") else {
        err("Error: --build <mode> is required.  Use --help for usage.\n");
        return 1;
    };
    match find_build_mode(mode) {
        Some(entry) => run_build_mode(entry, args),
        None => {
            err(&format!("Error: unknown build mode '{mode}'.  Valid modes: {}\n", build_mode_names()));
            1
        }
    }
}

/// Runs the command line on the process arguments; returns the exit code.
pub fn main_entry() -> i32 {
    // Build and Edit operations are also structure-pipeline steps.
    build_modifiers::register();
    let args = Args(std::env::args().skip(1).collect());
    if !is_cli_mode(&args) {
        err("This is the AtomForge command line; give --help for its modes. The desktop is the AtomForge program.
");
        return 1;
    }
    run_cli(&args)
}
