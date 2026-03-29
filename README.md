# AtomForge

AtomForge is an atomic structure builder for metallurgical simulation and atomistic modeling. It provides interactive tools for constructing, editing, and analyzing periodic crystals, grain boundaries, and nanostructures—all essential for preparing input structures for molecular dynamics and first-principles computational materials science.

## Core capabilities

- **Build structures from scratch** using crystal-system-aware builders for bulk crystals, cubic grain boundaries, and nanocrystals.
- **Edit and refine** atoms and lattice parameters interactively, with real-time visualization.
- **Generate polycrystalline microstructures** with Voronoi-based grain construction from a reference single crystal.
- **Create metallurgical defects** such as grain boundaries with explicit Sigma selection and in-plane replication.
- **Carve finite nanostructures** from any loaded reference structure using geometric shapes.
- **Analyze local order** with common-neighbour analysis and radial distribution functions.
- **Export for simulation** in VASP, LAMMPS, and other popular MD and DFT formats.
- **Inspect structure quality** through coordination statistics, symmetry detection, and lattice metrics.

## Platform support

- Linux
- Windows via MSYS2 UCRT64 / MinGW-w64

## Quick start

### Linux

```bash
sudo apt update
sudo apt install build-essential cmake pkg-config libglfw3-dev libglew-dev libglm-dev \
                 libopenbabel-dev libopenbabel3 libsymspg-dev

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/AtomForge
./build/AtomForge structure.cif
```

### Windows (MSYS2 UCRT64)

Install MSYS2, open the UCRT64 shell, then run:

```bash
pacman -Syu
# close and reopen the UCRT64 shell, then run
pacman -Su

pacman -S --needed mingw-w64-ucrt-x86_64-toolchain \
                  mingw-w64-ucrt-x86_64-cmake \
                  mingw-w64-ucrt-x86_64-glfw \
                  mingw-w64-ucrt-x86_64-glew \
                  mingw-w64-ucrt-x86_64-glm \
                  mingw-w64-ucrt-x86_64-openbabel \
                  mingw-w64-ucrt-x86_64-pkgconf

# optional: enables symmetry features through spglib
pacman -S --needed mingw-w64-ucrt-x86_64-spglib

cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/AtomForge.exe
./build/AtomForge.exe structure.cif
```

Note: run the Windows build from the MSYS2 UCRT64 shell. PowerShell will often miss the required toolchain and DLL paths.

## Build notes

- The project uses CMake and targets C++11.
- OpenGL 3.x, GLFW, GLEW, GLM, and Open Babel are required.
- spglib or symspg is optional; if found, symmetry-related features are enabled automatically.

To clean the current build tree:

```bash
cmake --build build --target clean
```

To remove all generated artifacts:

```bash
rm -rf build
```

## Portable Builds

Portable builds now bundle required runtime libraries and Open Babel plugins into the package folder automatically.

### Windows Portable

Build a portable ZIP package:

```bash
# In MSYS2 UCRT64 shell:
cmake -S . -B build-mingw -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_PORTABLE=ON
cmake --build build-mingw -- -j1

# Create package artifact in build-mingw/package:
cpack --config build-mingw/CPackConfig.cmake -B build-mingw/package
```

This generates a portable archive containing:
- `AtomForge.exe`
- required runtime DLLs
- Open Babel plugins under `openbabel/<version>/`
- `README.md`

Archive output path:
- `build-mingw/package/AtomForge-1.0.0-win64.zip`

Extract and run `AtomForge.exe` directly.

Note: do not reuse a build folder configured with another generator (for example NMake). Use a dedicated folder like `build-mingw` for MinGW.

### Linux Portable (Tarball)

Build a portable tarball:

```bash
cmake -S . -B build-portable -DCMAKE_BUILD_TYPE=Release -DBUILD_PORTABLE=ON
cmake --build build-portable -j
cpack --config build-portable/CPackConfig.cmake -B build-portable/package
```

This generates a portable archive containing:
- `AtomForge/bin/AtomForge`
- required `.so` runtime libraries copied next to the executable
- Open Babel plugins under `AtomForge/bin/openbabel/`
- `AtomForge/README.md`

To run the portable Linux package:
```bash
tar xzf build-portable/package/AtomForge-*.tar.gz
cd AtomForge/bin
./AtomForge
```

### Alternative: Build for System Installation

For standard system-wide installation (requires dependencies on target machine):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j install
```

On Linux, this installs to `/usr/local/bin` (may require `sudo`).

## Supported formats

**Structure files** (open and save): `.xyz`, `.cif`, `.pdb`, `.sdf`, `.mol`, `.vasp`, `.mol2`, `.pwi`, `.gjf`

**Rendered images**: `.png`, `.jpg`, `.svg`

You can also open a structure directly at launch by passing the file path as the first argument, for example `AtomForge structure.cif`.

## Core controls

### Scene navigation

| Action | Input |
| --- | --- |
| Rotate scene | Left drag |
| Zoom | Scroll wheel |
| Reset fitted default view | View -> Reset Default View |

### Selection

| Action | Input |
| --- | --- |
| Select one atom | Left click |
| Add or remove from selection | Ctrl + left click |
| Select all | Ctrl+A |
| Clear selection | Ctrl+D or Escape |
| Delete selection | Delete |

Box selection is available through Edit -> Box Select Mode. When enabled, right-drag draws a selection rectangle. Hold Ctrl to add to the current selection.

### File shortcuts

| Action | Shortcut |
| --- | --- |
| Open structure | Ctrl+O |
| Save structure as | Ctrl+S |
| Export rendered image | Ctrl+Shift+S |
| Undo | Ctrl+Z |
| Redo | Ctrl+Y or Ctrl+Shift+Z |

## Main workflows

### Open and inspect a structure

1. Use File -> Open.
2. Rotate and zoom the scene.
3. Toggle View -> Show Bonds and View -> Show Element as needed.
4. Open View -> Structure Info for composition, lattice, positions, and symmetry.

### Edit atoms and structure data

- Right-click a selection to substitute atoms, insert an atom at a midpoint, measure, or delete.
- Use Edit -> Edit Structure to modify lattice vectors and atom positions.
- Use Edit -> Atomic Sizes and Edit -> Element Colors to tune radii, colors, and shininess.
- Use Edit -> Transform Structure to apply a 3x3 matrix to periodic structures.

### Build and construct structures

- **Bulk Crystal**: Generate periodic cells from crystal system, space group, lattice parameters, and asymmetric-unit atoms. Ideal for preparing supercells and commensurate structures.
- **CSL Grain Boundary**: Construct cubic bicrystals with explicit Sigma selection from sc, bcc, fcc, or diamond source lattices. Control GB plane, in-plane replication, rigid translation, and overlap removal.
- **Nanocrystal**: Carve finite particles from any loaded reference structure using sphere, ellipsoid, box, cylinder, octahedron, truncated octahedron, or cuboctahedron geometry. Apply auto-centering and optional unit-cell padding.
- **Polycrystal**: Build a polycrystalline box from a loaded reference single crystal using Voronoi grain partitioning. Control box size, grain count, random seed, and grain orientations.

## Structure builders and analysis

AtomForge's builder system is designed to prepare structures for atomistic simulation—from perfect crystals to defected metallurgical systems.

### Bulk crystal builder

- Organizes space groups by crystal system (triclinic through cubic).
- Applies system-specific lattice constraints and trigonal settings.
- Expands asymmetric-unit atoms using symmetry operations to build the full unit cell.
- Perfect for creating supercells and preparing input for vasp, LAMMPS, and quantum espresso.

### CSL grain boundary builder

- Targeted at metallurgical simulation of grain boundaries in cubic systems.
- Generates sigma-list candidates from a rotation axis; select explicitly by Sigma, m, n, and angle.
- Controls GB plane, supercell dimensions, in-plane replication, overlap removal, and rigid translation.
- Outputs bicrystal structures ready for MD or DFT study of interface properties.

### Nanocrystal builder

- Carves finite particles from any loaded reference structure.
- Supports geometric shapes: sphere, ellipsoid, box, cylinder, octahedron, truncated octahedron, cuboctahedron.
- Auto-replicates periodic inputs for efficient particle generation and can pad with vacuum.
- Ideal for preparing nanoparticle input for cluster simulations and nanoscale property calculations.

### Polycrystal builder

- Uses a loaded single-crystal reference structure with unit-cell information.
- Builds grains by Voronoi tessellation inside a user-defined simulation box.
- Supports fully random grain orientations, fully specified Euler-angle orientations, or mixed specified/random orientation sets.
- Shows a live 3D preview of the reference crystal inside the builder dialog.
- Writes per-atom IPF-Z orientation colors for visualization and export-side metadata.

### Crystal orientation coloring

- View -> Color Structure By -> Crystal Orientation switches atom colors from element colors to cubic IPF-Z colors.
- An IPF triangle legend is shown in the main view when crystal-orientation coloring is active.
- When a structure is saved, AtomForge writes a companion `basename.atomforge-ipf` file when IPF data is available.
- On load, AtomForge first restores IPF from that saved metadata and only falls back to geometry-based reconstruction if no sidecar is found.

### Structure analysis

- **Common Neighbour Analysis**: Identifies local crystalline environments (FCC, HCP, BCC, ICO) and reports Honeycutt-Andersen-style pair signatures.
- **Radial Distribution Function**: Plots RDF with species filters, normalization, smoothing, and first-peak analysis for structure validation.
- **Coordination and bonding**: Automatic bond detection, PBC-aware coordination counting, and bond-length statistics for quality assurance.

## Display and measurements

- Bonds are inferred from covalent radii and rendered as split-color cylinders.
- Periodic image atoms are drawn at cell boundaries for periodic visualization.
- Distance and angle tools draw helper overlays in the scene.
- Atom Info reports element data, Cartesian and direct coordinates, and bond statistics.

## Menus at a glance

### File

- Open structures.
- Save the current structure.
- Export a rendered image.
- Close the current structure.

### Edit

- Undo and redo.
- Box selection mode.
- Structure editing.
- Atomic-size and element-color overrides.
- Supercell transform.

### View

- Element labels and bond visibility.
- Isometric and orthographic camera modes.
- Color by element type or crystal orientation.
- Structure information and measurement dialogs.
- Reset default view.

### Build

- Bulk Crystal.
- CSL Grain Boundary.
- Nanocrystal.
- Polycrystal.
- Nanocrystal.

### Analysis

- Common Neighbour Analysis.
- Radial Distribution Function.

### Help

- Manual.
- About.

## Windows troubleshooting

### Open Babel plugin path is missing

If AtomForge reports that Open Babel plugins cannot be found, set the plugin path before launch:

```bash
export BABEL_LIBDIR=/ucrt64/lib/openbabel/3.1.0/
```

If you moved/cloned the project and loading any structure fails, this is usually the reason.
Use the actual installed plugin directory for your Open Babel version, for example:

```bash
# Linux examples
export BABEL_LIBDIR=/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/
export BABEL_LIBDIR=/usr/local/lib/openbabel/3.1.1/

# MSYS2 UCRT64 examples
export BABEL_LIBDIR=/ucrt64/lib/openbabel/3.1.1/
export BABEL_LIBDIR=/ucrt64/lib/openbabel/3.1.0/
```

Quick check:

```bash
echo $BABEL_LIBDIR
ls "$BABEL_LIBDIR"
```

### Missing DLLs at startup

Run from the UCRT64 shell, or ensure `C:/msys64/ucrt64/bin` is present in `PATH`.

### CMake cannot find dependencies

Verify that `pkgconf` and the required MSYS2 packages are installed in the same UCRT64 environment used for the build.

## Project layout

```text
CMakeLists.txt              cross-platform build configuration
src/
  main.cpp                  application entry point and frame loop
  ElementData.cpp/.h        element symbols, colors, radii, metadata
  app/                      editor coordination, file handling, interaction glue
  camera/                   orbit camera and callbacks
  graphics/                 renderer, scene buffers, picking, meshes, shaders
  io/                       structure loading through Open Babel
  math/                     structure and geometry helpers
  ui/                       Dear ImGui windows, dialogs, and builders
imgui/                      bundled Dear ImGui sources and backend bindings
```
