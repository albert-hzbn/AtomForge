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
  that registry.

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
cmake --build build-core --parallel
ctest --test-dir build-core --output-on-failure
```

`cmake/AtomForgeCore.cmake` defines the reusable target;
`cmake/AtomForgeTests.cmake` defines tests for both build modes. The normal desktop
build remains the default. `BUILD_TESTING=OFF` omits tests. If Python is not found,
set `Python3_EXECUTABLE` to an interpreter path to enable the Python/CLI tests.

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
