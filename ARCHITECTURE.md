# AtomForge architecture

## C++ dependency boundaries

- `src/model/Structure.h` owns atom and structure data, including per-atom metadata
  and deletion rules. Include it when only domain types are needed.
- `src/util/ElementData.*` owns element names, radii, masses, and default colors.
- `src/algorithms/` owns structure generation and scientific calculations.
  CNA and RDF expose typed parameter/result objects in `atomforge::analysis`.
  Their implementations have no ImGui dependencies.
- `src/ui/` adapts controls and presentation to the algorithms. CNA, RDF, and ADF
  use `atomforge::BackgroundTask<Result>` to run work on captured snapshots.
- `src/io/StructureLoader.*` handles file conversion and metadata persistence.
  Its header still includes the model and element declarations for compatibility;
  new domain-only code should include those headers directly.
- `src/app/` coordinates tabs, editing, rendering, and file workflows.
- `src/graphics/` manages graphics resources and scene rendering. Ray picking is
  a geometry operation and is also available through the core target.
- `src/cli/CLIMode.cpp` registers build modes in `kBuildModes`, pairing each mode
  with its execution and help functions. Dispatch and valid-mode diagnostics use
  that registry. The modes themselves live in family files declared by
  `src/cli/BuildModes.h`: `BuildCrystalCLI.cpp` (bulk, custom, solid solution,
  primitive, surface, SQS, vacancy, strain), `BuildDefectCLI.cpp` (dislocation,
  grain boundary, interface, stacking fault) and `BuildNanoCLI.cpp` (polycrystal,
  nanocrystal, amorphous, nanowire, core-shell). `src/cli/CliArgs.*` holds the
  argument parsing shared by every CLI mode; numeric flags reject malformed or
  non-finite values with an error naming the flag. To add a mode, implement
  `printHelpX`/`runX` in the matching family file, declare them in
  `BuildModes.h` and add a `kBuildModes` row. The manual's CLI reference is
  generated from these help pages.
- `src/science/` is the native scientific-tools library (`atomforge_science`):
  analyses, simulations, file readers and writers, with no GUI dependency. It is
  driven by JSON requests from the desktop, `AtomForge --science` and batches.

### Reusable core

`AtomForge::Core` is the CMake alias for the `atomforge_core` static library. It
contains element data, picking, amorphous/solid-solution builders, CNA, RDF, ADF,
SRO, void analysis, Voronoi computation, mesh loading, and cell sculpting.
It requires C++17, GLM, and thread support. OpenMP is optional.

The core has no OpenGL, GLFW, ImGui, Open Babel, or spglib dependency. The desktop
target owns builders that still depend on its symmetry or scene-building code.
This is an incremental boundary: move additional implementations into the core
when their dependencies and independent tests support it.

```sh
cmake -S . -B build-core -DATOMFORGE_BUILD_APP=OFF
cmake --build build-core --parallel 1
ctest --test-dir build-core --output-on-failure
```

`cmake/AtomForgeCore.cmake` defines the reusable target;
`cmake/AtomForgeTests.cmake` defines tests for both build modes. The normal desktop
build remains the default. `BUILD_TESTING=OFF` omits tests. If Python is not found,
set `Python3_EXECUTABLE` to an interpreter path to enable the Python/CLI tests.

### Scientific tool registry

Each native tool is described once in `src/science/ScienceCatalog.cpp`: id,
title, one-line description and method notes, Analysis-menu category,
parameters (JSON request keys with labels, defaults and input kinds) and the
desktop `ScienceToolView`. Behaviour is attached through id-keyed tables, one
named function per tool:

| Table | File | Purpose |
|---|---|---|
| `toolRunners()` | `ScienceTools.cpp` | runs a validated request (required) |
| `plotBuilders()` | `ResultPlots.cpp` | result plots (optional) |
| `propertyBuilders()` | `AtomProperties.cpp` | per-atom properties for colouring (optional) |

`runTool`, `resultPlots` and `perAtomProperties` look tools up in these tables;
the catalog drives request validation, the CLI catalog, the batch runner and
the desktop window. The registry test in `tests/science_tools_regressions.cpp`
fails if a catalog tool has no runner or a table names an unknown tool.

Parameters of kind `structure` in the catalog are parsed as structures
(`Parameters::structure(name)`), so a tool may name several structure inputs
(for example `structure` and `bulk`). `TrajectoryStream` gives random access to
trajectory frames without loading them all; tools that read whole trajectories
(such as `trajectory-structure`) and Trajectory playback use it.

`atomforge_science_native` is a shared library around the same registry
(`src/science/NativeAPI.cpp`: `afs_catalog`, `afs_run`, `afs_free`) for the
Python module `atomforge.native`, so every catalog tool is callable from Python
without further bindings.

To add a tool: add its catalog entry, implement the calculation (taking
`Parameters` and returning `ToolOutput`), add a runner row and, if useful,
plot and property builders; then describe its method, units and limits in the
manual (`docs/manual/chapters/14-workflows.tex`). The desktop picks the window
layout from the entry's view: `ui/ScientificToolsDialog.cpp` holds dialog
state, persistence and running, and `ui/ScientificToolsLayouts.cpp` draws the
parameter widgets, result views and one layout per view.

### Adding an analysis

1. Define parameter and result types under `src/algorithms/`. Prefer named fields
   over long lists of positional options. Return a useful invalid-result message
   for invalid user input. Validate inputs at the calculation boundary, since
   callers may not use the GUI's constrained controls.
2. Keep calculation functions independent of ImGui and file dialogs. Add suitable
   independent implementations to `AtomForge::Core`.
3. Add headless tests with known geometries and expected observables.
4. In a dialog, capture the structure and parameters by value when starting a
   `BackgroundTask`. Call `poll()` on the UI thread and render its typed result
   or error. Do not capture dialog state by reference in a worker.

`BackgroundTask` is single-owner state, not a general concurrent container.
`running()` remains true until `poll()` consumes completion. A running task
rejects replacement; worker exceptions become owner-readable error strings.
Destruction waits for unfinished work. Cancellation is not provided, so closing
an owner can wait for a long calculation. Tests cover these lifecycle contracts.

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

To add a format, implement its reader/writer under `_formats/`, register the
extensions in `_io._FORMAT_MAP`, and add round-trip and malformed-input tests.
Both reading and writing use `_resolve_format`, including extensionless POSCAR
and CONTCAR names. Existing private codec/math imports from `_io` are re-exported
for compatibility. Shared numerical helpers must not depend on Structure or I/O.

The package remains dependency-free and retains its Python 3.8 minimum. Its
existing `find_packages()` setup discovers the new `_formats` subpackage.

## Refactor validation and scope

The Windows desktop build with OpenMP/spglib and a separate core-only build with
OpenMP disabled both compile. Six desktop CTest suites and five core-only suites
cover core regressions, extracted analyses, worker lifecycle/error handling,
existing Python APIs, notebook rendering, file-format regressions, and CLI smoke
tests where the application is built.

Fifteen exports (five formats across cluster, cubic, and triclinic structures)
match the pre-refactor files byte for byte. All package modules parse using
Python 3.8 syntax rules; tests executed on the available Python 3.14 interpreter.
`git diff --check` passes.

Numerical implementations and supported file formats were preserved during
extraction. Invalid non-finite inputs are additionally rejected at the new CNA
and RDF entry points. CNA/RDF settings now belong to each dialog instance;
worker exceptions in CNA/RDF/ADF are reported to the UI instead of escaping a
raw thread.

The refactor does not replace every legacy UI or algorithm implementation.
Vendored dependencies are unchanged. Interactive GUI behavior, Linux/static
packaging, and wheel installation have not been validated in this pass; the
available Python interpreter does not include setuptools or IPython.

## Electronic post-processing extension

`src/electronic` contains the independent scalar-grid model, format adapters,
field calculus, density integration, Fourier analysis, electrostatics and mesh
extraction. It is part of `AtomForge::Core`. `PythonAPI.cpp` exposes a versioned,
exception-safe C boundary in `atomforge_electronic`; the lazy ctypes facade in
`python/atomforge/electronic` shares the native calculations with the desktop.
`ElectronicPostProcessingDialog` owns fields/results and uses `BackgroundTask`
for loading and computation.
