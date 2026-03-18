# AtomForge

An OpenGL-based molecular structure viewer and editor using Dear ImGui and Open Babel.

## Platform support

- Linux
- Windows (MSYS2 UCRT64 / MinGW-w64)

## Prerequisites

### Linux (Debian/Ubuntu or derivatives)

```bash
sudo apt update
sudo apt install build-essential cmake libglfw3-dev libglew-dev libglm-dev \
                 libopenbabel-dev libopenbabel3 libsymspg-dev
```

### Windows (MSYS2 UCRT64)

Install MSYS2, open the **UCRT64** shell, then install dependencies:

```bash
pacman -Syu
# close and reopen UCRT64 shell, then run
pacman -Su

pacman -S --needed mingw-w64-ucrt-x86_64-toolchain \
                  mingw-w64-ucrt-x86_64-glfw \
                  mingw-w64-ucrt-x86_64-glew \
                  mingw-w64-ucrt-x86_64-glm \
                  mingw-w64-ucrt-x86_64-openbabel \
                  mingw-w64-ucrt-x86_64-pkgconf
```

Optional (enables symmetry features via spglib):

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-spglib
```

> **Note:** This project uses OpenGL 3.0+ (GLSL 130) and GLFW.

## Build

### Linux

From the `src/` directory:

```bash
make
```

Then run:

```bash
./AtomForge
```

### Windows (MSYS2 UCRT64)

From the `src/` directory in the UCRT64 shell:

```bash
mingw32-make
```

Then run:

```bash
./AtomForge.exe
```

## Windows troubleshooting

- **Open Babel plugin error** (`Unable to find OpenBabel plugins`):
  set the plugin path before launching.

```bash
export BABEL_LIBDIR=/ucrt64/lib/openbabel/3.1.0/
```

- **Application does not start due to missing DLLs**:
  run from the UCRT64 shell, or add `C:/msys64/ucrt64/bin` to `PATH`.

- **`make: command not found` in UCRT64**:
  use `mingw32-make` instead of `make`.

## Usage

### Opening files

Use **File → Open…** (`Ctrl+O`) to browse and load a structure file.  
Supported formats: `.xyz`, `.cif`, `.pdb`, `.sdf`, `.mol`, `.vasp`, `.mol2`, `.pwi`, `.gjf`.

Use **File → Save As…** (`Ctrl+S`) to export the current structure.

### Camera

| Action | Input |
|---|---|
| Rotate | Left drag |
| Zoom | Scroll wheel |

Default camera behavior:
- New structures are framed automatically so the full structure is visible.
- Default orientation is isometric.
- Use **View → Reset Default View** to re-apply isometric framing at any time.

### Atom selection

| Action | Input |
|---|---|
| Select atom | Left click |
| Add / remove from selection | Ctrl + click |
| Select all | Ctrl+A |
| Deselect all | Ctrl+D or Escape |
| Delete selected | Delete |

Box selection:
- Enable **Edit → Box Select Mode**, then right-drag to select atoms in a screen rectangle.
- Hold `Ctrl` while box-selecting to add to the current selection.

### Right-click context menu

- **Substitute Atom…** — replace all selected atoms with a chosen element (periodic table picker)
- **Insert Atom at Midpoint…** — insert a new atom at the centroid of ≥ 2 selected atoms
- **Measure Distance** — available when exactly 2 atoms are selected
- **Measure Angle** — available when exactly 3 atoms are selected (angle at the 2nd selected atom)
- **Atom Info** — available when exactly 1 atom is selected (element, Cartesian, direct coordinates)
- **Delete** / **Deselect**

### View menu

- **Show Element** — toggle on-screen element symbol labels
- **Show Bonds** — display bond cylinders split by atom colours
- **Structure Info** — show composition, lattice metrics, positions, and symmetry summary
- **Measure Distance (2 selected)**
- **Measure Angle (3 selected)**
- **Atom Info (1 selected)**
- **Reset Default View** — restore isometric view and fit structure to window

### Edit menu

- **Undo / Redo** — full edit-history snapshots (`Ctrl+Z`, `Ctrl+Y`, `Ctrl+Shift+Z`)
- **Box Select Mode** — switch right mouse drag to rectangle selection mode
- **Edit Structure...** — modify lattice vectors and atom list (add/edit/delete)
- **Atomic Sizes…** — adjust per-element covalent radii (literature defaults: Cordero et al., *Dalton Trans.* 2008)
- **Element Colors…** — override CPK colours per element and tune material shininess
- **Transform Structure…** — apply a 3×3 matrix transformation to all atom positions (only available when the structure has a unit cell)

### Help menu

- **Manual** — complete in-app reference for controls and menu actions
- **About** — project overview, libraries, and references

### Analysis menu

- **Common Neighbour Analysis…** — run CNA from a dedicated dialog.
- Pair-signature distribution is reported in Honeycutt-Andersen style `1-j-k-l`.
- Per-atom environment summary is reported (FCC-like/HCP-like/BCC-like/ICO-like/Unknown).
- Per-atom CNA details are listed (coordination, dominant signature, environment).
- PBC-aware analysis is used when a unit cell is present (toggle in dialog).
- **Radial Distribution Function…** — compute and plot RDF from a dedicated dialog.
- RDF controls include species filters, PBC toggle, normalization, radius range, bin count, smoothing, and plot overlays.
- RDF results include the plotted curve, per-bin table, density/volume summary, and first-peak / first-minimum estimates.

### Build menu

- **Bulk Crystal…** — build a full periodic unit cell from a selected crystal system, space group, lattice parameters, and asymmetric-unit atoms.
- Space-group selection is grouped by crystal system (triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic).
- Lattice inputs are constrained by the selected crystal system; trigonal currently uses the hexagonal setting.
- Asymmetric-unit atoms can be added, edited, deleted, and substituted directly in the dialog before symmetry expansion.
- Generation applies the selected space-group symmetry operations to fill the full unit cell and replaces the current structure.
- **CSL Grain Boundary…** — create a bicrystal grain-boundary structure from cubic lattice templates only (`sc`, `bcc`, `fcc`, `diamond`).
- The builder does not use the loaded structure as the source lattice; it generates from selected basis, lattice parameter, and element.
- Sigma can be chosen explicitly by generating a Sigma candidate list from the selected axis and choosing one entry (`Sigma`, `m`, `n`, angle).
- In-plane replication controls (`r2 r3`) are available to expand the generated GB cell without creating multiple GB images along the normal.
- Builder controls include rotation axis `[u v w]`, `(m, n)` misorientation parameters (Sigma/angle), GB plane `(h k l)`, dimensions, overlap removal, optional rigid translation, and box padding.
- Generation replaces the current structure and reports atoms in/out, overlap removals, Sigma, and misorientation angle.

### Measurements and overlays

- Distance and angle tools open result dialogs and draw dashed helper lines.
- **Atom Info** shows element name, atomic number, Cartesian/direct coordinates, and coordination number with average/min/max bond length (bond stats use PBC when a unit cell exists).
- Dashed helper lines are cleared when selection changes or when the dialog is confirmed.
- **Show Element** labels are rendered for periodic image atoms as well.

### Periodic and bonding display

- Boundary atoms are duplicated across unit-cell faces/edges/vertices for periodic visualization.
- Bonds are inferred from covalent radii and rendered as thicker cylinders with half-colour per bonded atom.

### Structure information

- **View → Structure Info** includes:
  - Total atoms, element counts, and formula
  - Lattice lengths, angles, and volume
  - Atomic positions in Cartesian and fractional coordinates
  - Space group and point group using spglib when available

### Save/export formats

- **File → Save As** supports: `.xyz`, `.cif`, `.vasp`, `.pdb`, `.sdf`, `.mol2`, `.pwi`, `.gjf`
- When transform mode is active, save uses the expanded supercell representation.

## Project layout

```
src/
  main.cpp                   — application entry point, render loop
  ElementData.cpp/.h         — element symbols, radii, colours
  camera/                    — camera orbit controller
  graphics/                  — OpenGL renderer, scene buffers, picking, shaders
  math/                      — math utilities
  ui/                        — Dear ImGui panels and dialogs
  io/                        — structure loading (Open Babel)
imgui/                       — bundled Dear ImGui sources
backends/                    — ImGui GLFW / OpenGL backends
```
