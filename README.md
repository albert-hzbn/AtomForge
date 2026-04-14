# AtomForge

AtomForge is an interactive atomic structure builder for metallurgical simulation and atomistic modeling. It is designed to help researchers and engineers create, edit, inspect, and export structures used in molecular dynamics (MD) and first-principles workflows.

## Highlights

- Build structures from scratch using crystal-aware workflows for bulk crystals, cubic grain boundaries, and nanocrystals.
- Generate custom finite structures by filling imported 3D mesh volumes with atoms from a reference crystal.
- Combine multiple loaded structures in a dedicated 3D Merge Structures workflow with per-structure transform gizmos.
- Create polycrystalline microstructures via Voronoi grain construction.
- Analyze local order with common neighbour analysis (CNA) and radial distribution function (RDF) tools.
- Export structures for simulation (VASP and other supported formats).
- Export publication-ready images (PNG, JPEG, SVG), including high-resolution scaling.
- Use dark/light themes with HiDPI-aware UI scaling for readable overlays and dialogs.

## Platform support

- Linux
- Windows (MSYS2 UCRT64 / MinGW-w64)

## Quick start

### Linux

```bash
sudo apt update
sudo apt install build-essential cmake pkg-config libglfw3-dev libglew-dev libglm-dev \
                 libopenbabel-dev libsymspg-dev

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

Note: Build and run from the MSYS2 UCRT64 shell. PowerShell often misses required toolchain and DLL paths.

## Build notes

- The project uses CMake and targets C++17.
- OpenGL 3.x, GLFW, GLEW, GLM, and Open Babel are required.
- spglib/symspg is optional; when available, symmetry-related features are enabled automatically.

Clean current build artifacts:

```bash
cmake --build build --target clean
```

Remove the full build directory:

```bash
rm -rf build
```

## Portable builds

Portable packaging bundles runtime libraries and Open Babel plugins into the output package.

### Windows portable ZIP

```bash
# In MSYS2 UCRT64 shell:
cmake -S . -B build-mingw -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_PORTABLE=ON
cmake --build build-mingw -- -j1

# Create package artifact in build-mingw/package:
cpack --config build-mingw/CPackConfig.cmake -B build-mingw/package
```

Archive output:

- `build-mingw/package/AtomForge-1.0.0-win64.zip`

Contains:

- `AtomForge.exe`
- required runtime DLLs
- Open Babel plugins under `openbabel/<version>/`
- `README.md`

Extract and run `AtomForge.exe` directly.

Note: Do not reuse a build folder configured with a different generator (for example, NMake). Use a dedicated folder such as `build-mingw`.

### Linux portable tarball

```bash
cmake -S . -B build-portable -DCMAKE_BUILD_TYPE=Release -DBUILD_PORTABLE=ON
cmake --build build-portable -j
cpack --config build-portable/CPackConfig.cmake -B build-portable/package
```

Contains:

- `AtomForge/bin/AtomForge`
- required `.so` libraries copied next to the executable
- Open Babel plugins under `AtomForge/bin/openbabel/`
- `AtomForge/README.md`

Run:

```bash
tar xzf build-portable/package/AtomForge-*.tar.gz
cd AtomForge/bin
./AtomForge
```

### Alternative: system installation

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j install
```

On Linux, this installs to `/usr/local/bin` (may require `sudo`).

## Supported formats

Structure files (open/save): `.xyz`, `.cif`, `.pdb`, `.sdf`, `.mol`, `.vasp`, `.mol2`, `.pwi`, `.gjf`

Rendered image export: `.png`, `.jpg`, `.svg`

You can also open a structure at launch by passing a file path, for example:

```bash
AtomForge structure.cif
```

## Core controls

### Scene navigation

| Action | Input |
| --- | --- |
| Rotate scene | Left drag |
| Zoom | Scroll wheel |
| Reset fitted default view | `R` or View -> Reset Default View |

### Selection

| Action | Input |
| --- | --- |
| Select one atom | Left click |
| Add/remove from selection | Ctrl + left click |
| Select all | Ctrl + A |
| Clear selection | Ctrl + D or Escape |
| Delete selection | Delete |

Box selection is available from Edit -> Box Select Mode. When enabled, right-drag draws a selection rectangle. Hold Ctrl to add to current selection.

### File shortcuts

| Action | Shortcut |
| --- | --- |
| Open structure | Ctrl + O |
| Save structure as | Ctrl + S |
| Export rendered image | Ctrl + Shift + S |
| Undo | Ctrl + Z |
| Redo | Ctrl + Y or Ctrl + Shift + Z |

## Main workflows

### Open and inspect

1. Use File -> Open.
2. Navigate with rotate/zoom controls.
3. Toggle View -> Show Bonds and View -> Show Element as needed.
4. Open View -> Structure Info for composition, lattice, positions, and symmetry.

### Edit structure data

- Right-click a selection to substitute atoms, insert midpoint atoms, measure, or delete.
- Use Edit -> Edit Structure to modify lattice vectors and atom positions.
- Use Edit -> Atomic Sizes and Edit -> Element Colors to adjust visual properties.
- Use Edit -> Transform Structure to apply a 3x3 matrix to periodic structures.
- Use Edit -> Merge Structures to load, arrange, and merge multiple structures in an interactive 3D preview.

### Build structures

- **Bulk Crystal**: Create periodic cells from crystal system, space group, lattice parameters, and asymmetric-unit atoms.
- **CSL Grain Boundary**: Build cubic bicrystals with explicit Sigma selection and control over GB plane, in-plane replication, translation, and overlap handling.
- **Nanocrystal**: Carve finite particles from loaded reference structures (sphere, ellipsoid, box, cylinder, octahedron, truncated octahedron, cuboctahedron).
- **Custom Structure**: Fill imported mesh volumes (OBJ/STL) with atoms from a reference crystal, with side-by-side 3D previews.
- **Interface Builder**: Build interfaces from two input structures with candidate matching and preview.
- **Polycrystal**: Generate Voronoi-based polycrystalline structures from a reference crystal.

### Merge structures

- Open **Edit -> Merge Structures** to combine multiple structures into one.
- Drag-and-drop supported structure files while the dialog is open.
- Select individual structures in the list or preview, then use the 3D gizmo for translate/rotate.
- Orbit and zoom the preview, optionally show/hide the merged bounding box, then apply **Merge Structures**.

## Builder and analysis details

### Bulk crystal builder

- Organizes space groups by crystal system (triclinic to cubic).
- Applies system-specific lattice constraints and trigonal settings.
- Expands asymmetric-unit atoms using symmetry operations to produce full unit cells.

### CSL grain boundary builder

- Designed for cubic grain boundary studies in metallurgy.
- Generates Sigma candidates from rotation axis input.
- Controls GB plane, dimensions, replication, overlap removal, and rigid translation.

### Nanocrystal builder

- Carves finite structures from loaded references.
- Supports sphere, ellipsoid, box, cylinder, octahedron, truncated octahedron, and cuboctahedron geometries.
- Can auto-replicate periodic inputs and apply vacuum padding.

### Polycrystal builder

- Uses unit-cell-aware reference structures.
- Builds grains via Voronoi tessellation in a user-defined box.
- Supports random, explicit, and mixed orientation assignment.

### Custom structure builder

- Accepts drag-and-drop crystal and mesh inputs.
- Provides live 3D previews for both reference crystal and mesh model.
- Produces finite atomistic structures constrained by the model volume.

### Crystal orientation coloring

- View -> Color Structure By -> Crystal Orientation switches from element colors to cubic IPF-Z colors.
- Displays an IPF triangle legend in the main view when active.
- Saves companion `basename.atomforge-ipf` metadata when IPF data is available.
- Restores from sidecar metadata on load when present, with geometry-based fallback otherwise.

### Analysis tools

- **Common Neighbour Analysis (CNA)**: Classifies local motifs (FCC, HCP, BCC, ICO) with Honeycutt-Andersen-style signatures.
- **Radial Distribution Function (RDF)**: Provides species filtering, normalization, smoothing, and first-peak analysis.
- **Coordination and bonding**: Includes automatic bond detection, PBC-aware coordination counting, and bond-length statistics.

## Display and measurement

- Bonds are inferred from covalent radii and rendered as split-color cylinders.
- Periodic image atoms are shown at cell boundaries for periodic context.
- Distance and angle tools draw overlays directly in the scene.
- Atom Info reports element, Cartesian/direct coordinates, and bond statistics.
- View -> Select Theme supports dark and light themes with matching overlay colors.
- UI scales automatically for high-resolution and HiDPI displays.

## Menus at a glance

### File

- Open structures
- Save structures
- Export rendered images
- Close current structure

### Edit

- Undo/redo
- Box selection mode
- Structure editing
- Atomic-size and element-color overrides
- Supercell transform

### View

- Element labels and bonds
- Isometric/orthographic camera modes
- Element/crystal-orientation color modes
- Structure and measurement dialogs
- Reset default view

### Build

- Bulk Crystal
- CSL Grain Boundary
- Nanocrystal
- Custom Structure
- Interface Builder
- Polycrystal

### Analysis

- Common Neighbour Analysis
- Radial Distribution Function

### Help

- Manual
- About

## Windows troubleshooting

### Open Babel plugin path missing

If AtomForge reports missing Open Babel plugins, set `BABEL_LIBDIR` before launch:

```bash
export BABEL_LIBDIR=/ucrt64/lib/openbabel/3.1.0/
```

If loading fails after moving/cloning the project, this is often the cause. Use the installed Open Babel plugin directory for your version, for example:

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

Verify `pkgconf` and the required MSYS2 packages are installed in the same UCRT64 environment used for build and run.

## Project layout

```text
CMakeLists.txt              cross-platform build configuration
src/
  main.cpp                  application entry point and frame loop
  app/                      editor coordination, workflow control, undo/redo, interactions
  algorithms/               structure generation and geometry/model processing
  camera/                   orbit camera and callbacks
  graphics/                 renderer, scene buffers, picking, meshes, shaders
  io/                       structure loading and saving through Open Babel
  math/                     structure and geometry helpers
  ui/                       Dear ImGui windows, dialogs, and builders
  util/                     shared utilities and element metadata
imgui/                      bundled Dear ImGui sources and backend bindings
```
