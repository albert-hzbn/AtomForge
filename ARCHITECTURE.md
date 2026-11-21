# AtomForge architecture

AtomForge is a Cargo workspace in `rust/`. The crates depend on each other in
layers; nothing below `af-render` knows about the GPU, and nothing below
`af-app` knows about the interface.

| Crate | Responsibility |
| --- | --- |
| `af-core` | `Structure` and per-atom data, elements (names, radii, masses, colours), cells and supercells, neighbour lists, picking, JSON, number formatting, `BackgroundTask`, display settings and exactness helpers (`cxx` sorting, `crmath`, `glm`, `eigen`) |
| `af-io` | Structure files: CIF, POSCAR, XYZ/extXYZ, PDB, SDF/MOL/MOL2, LAMMPS data and dumps, the inputs of Quantum ESPRESSO, CASTEP, DFTB+, XCrySDen and CP2K; trajectories with random frame access |
| `af-spglib` | Bindings to spglib 2.7 (symmetry, standardisation, primitive cells) |
| `af-chemistry` | Molecules: SMILES, the compound library, UFF force field, embedding, point groups |
| `af-algorithms` | Builders and structural analyses: bulk crystals (230 space groups), solid solutions, SQS, CSL grain boundaries, interfaces, stacking faults, dislocations (isotropic and Stroh), nanocrystals and Wulff shapes, polycrystals, amorphous packing, cell sculptor, Shape Modeller, CNA, RDF, ADF, SRO, Voronoi, interstitial voids, Nye tensors, elastic fields |
| `af-electronic` | Volumetric data (CHGCAR, LOCPOT, cube, XSF), field calculus, integration, Fourier analysis, Bader, isosurfaces, LOBSTER, Wannier90, and AtomForge projects (`project.rs`, formats 1–7) |
| `af-science` | The scientific-tool registry: the catalog, one runner per tool, result plots, per-atom properties, batches |
| `af-pipeline` | The non-destructive structure pipeline: modifiers, expressions, the editor with its stage cache, and option text (`options.rs`) |
| `af-cloud` | Out-of-core atom clouds (`.afcloud`): chunked, randomly ordered 8-byte atoms and the converters from LAMMPS dumps and XYZ |
| `af-render` | wgpu rendering: the structure scene (sphere impostors, bonds, cell), colouring, measurements, selection regions, grab, polyhedra and the cloud renderer |
| `af-cli` | The command line (`--build`, `--analyze`, `--convert`, `--render`, `--science`, `--pipe`, `--cloud`) and the Build steps of the pipeline |
| `af-app` | The desktop (egui/eframe on wgpu); its binary, `AtomForge`, also runs the command line |
| `af-ffi-science`, `af-ffi-electronic` | C ABIs of the science tools and electronic post-processing for the Python package |

## Desktop

`af-app/src/main.rs` registers the Build steps, runs the command line when the
arguments name a mode (`af_cli::is_cli_mode`) and otherwise opens the window
(`--renderer` picks the wgpu backend). `app.rs` holds the tabs
(`document.rs`: a structure, its undo history, camera, selection, scene and
pipeline), the menus, toolbar, status bar, viewport and every window. The
viewport is drawn by `af-render` on the interface's own wgpu device into a
texture that egui shows, so there is no copy through the CPU.

Menus come from one table: `af_cli::menu::MENU_LOCATIONS` gives the menu,
submenu and label of every pipeline step, and the Analysis menu lists the
science catalog by category. A step without a dedicated dialog opens the
generic operation window (`operation.rs`, its form from the step's parameter
kinds in `params_form.rs`).

Dedicated dialogs implement `structure_dialog::StructureDialog`: the step ids
they serve, `show` (returning an `Outcome`: a built structure, an edit, a
selection, steps to add to the pipeline or a pipeline step update), dropped
files, the electronic workspace for projects, and `edit_step` for
**Edit in Dialog** on a pipeline step (`StepEdit`: the step, its parameters and
the structure entering it; Update Step writes the parameters back to the tab it
came from). `structure_dialog::all()` lists them; adding a dialog is a new
module and one row there. Bulk Crystal and Insert Dislocation are app fields
because they predate the trait.

Long calculations run on threads with captured snapshots (`af_core::task`,
or a channel polled each frame); script commands run synchronously.

`ATOMFORGE_UI_SCRIPT=FILE` drives the desktop one command per frame
(`rust/tests/ui/*.afscript`, run by `rust/tools/run_ui_tests.py`): opening
files, menu steps, dialogs (`dialog STEP COMMAND ARGS`), the pipeline, picking,
camera, projects, screenshots, and `expect-...` checks. A failed check exits
with code 3. Autosave and animations are off in scripted runs.

## Structure pipeline

`af-pipeline` holds the non-destructive pipeline. `Pipeline` is an ordered list
of `Modifier` values (type id, JSON parameters, enabled flag, label);
evaluation runs them on an input `Structure` and returns `PipelineData`
(structure plus selection), with a `StageResult` per step. The editor
(`editor.rs`) caches each stage, so an edit to step *i* re-runs only steps *i*
and later, and gives the list operations the panel uses. The text syntax
(`replicate 2 2 2 | select-expression fx < 0.5 | delete-selected`) and JSON
round-trip.

Modifier families live in `steps/`: selection, editing and cell steps
(`basic.rs`, `edit.rs`), generators (`generator*.rs`), analysis and property
steps (`analysis.rs`), the Shape Modeller (`shape.rs`), the Python-script step
(`python.rs`) and the steps generated from the science catalog (`science.rs`).
`af_cli::build_modifiers::register` adds every `--build` mode as a Build step at
start-up, running the same code as the menus. Build steps keep their
command-line flags in an `options` parameter; `af_pipeline::options` reads and
writes them for the dialogs.

To add a modifier, add an entry to its family: id, title, category, help,
parameters (name, label, kind, default, choices) and an apply function on
`PipelineData`. The panel, text syntax, files and command line accept it with
no other change; add a menu row in `af-cli/src/menu.rs` (tests fail without one).

## Scientific tool registry

Each tool is described once in `af-science/src/catalog_table.rs`: id, title,
description and method notes, category, parameters (request keys with labels,
defaults and kinds) and its window view. Behaviour is attached by id: a runner
(required), result plots (`plots.rs`) and per-atom properties
(`atom_properties.rs`). The catalog drives request validation, the CLI, batches,
the Python module and the desktop's tool window.

## Large datasets

Structures above five million atoms, LAMMPS dumps and `.afcloud` files open as
atom clouds instead of editable structures: `af-cloud` stores atoms by space in
chunks of 8-byte atoms in random order (a chunk prefix is a uniform
subsample); `af_render::cloud` draws them as sphere impostors with culling,
prefix level of detail, point and GPU-memory budgets and a streaming thread.

## Exactness and tests

The port from C++ was checked against the C++ program it replaced: golden
dumpers (`rust/tools/golden`, built against the C++ sources while they existed)
wrote reference outputs into `rust/tests/golden`, and the Rust tests compare
with them, most of them bit for bit. Helpers reproduce what the C++ results
depended on: libstdc++ sorting (`af_core::cxx`), correctly rounded
trigonometry (`af_core::crmath`), glm arithmetic (`af_core::glm`) and the C++
random engines. `cargo test --workspace --release` runs them; the desktop is
covered by the UI scripts and the Python package by `rust/tools/run_python_suite.sh`.

## Python package

The public API remains `Atom`, `Structure`, `load`, `save`, and `view`.

| Module | Responsibility |
| --- | --- |
| `_structure.py` | Domain data, editing, and transforms |
| `_elements.py` | Colors, masses, and symbol normalization |
| `_cell.py` | Lattice parameters and coordinate conversion |
| `_formats/xyz.py` | XYZ/extXYZ codec |
| `_formats/vasp.py` | POSCAR/CONTCAR codec |
| `_formats/pdb.py` | PDB codec |
| `_formats/cif.py` | P1/pre-expanded CIF codec |
| `_formats/lammps.py` | LAMMPS data codec |
| `_io.py` | Filename dispatch and compatibility imports |
| `_viewer.py` | External GUI process lifecycle |
| `_notebook.py` | Inline notebook renderer |
| `native.py` | The science tools through `atomforge_science_native` |
| `electronic/` | Electronic post-processing through `atomforge_electronic` |

To add a format, implement its reader/writer under `_formats/`, register the
extensions in `_io._FORMAT_MAP`, and add round-trip and malformed-input tests.
The package itself is pure Python with a Python 3.8 minimum; the native
libraries are found beside it, beside the `AtomForge` program, or through
`ATOMFORGE_SCIENCE_LIBRARY` / `ATOMFORGE_ELECTRONIC_LIBRARY`.
