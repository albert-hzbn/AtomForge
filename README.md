# Atoms Editor

An OpenGL-based molecular structure viewer and editor using Dear ImGui and Open Babel.

## Prerequisites

On Debian/Ubuntu (or derivatives):

```bash
sudo apt update
sudo apt install build-essential cmake libglfw3-dev libglew-dev libglm-dev \
                 libopenbabel-dev libopenbabel3 libsymspg-dev
```

> **Note:** This project uses OpenGL 3.0+ (GLSL 130) and GLFW.

## Build

From the `src/` directory:

```bash
make
```

Then run:

```bash
./atoms_editor
```

## Usage

### Opening files

Use **File → Open…** (`Ctrl+O`) to browse and load a structure file.  
Supported formats: `.cif`, `.mol`, `.pdb`, `.xyz`, `.sdf`.

### Camera

| Action | Input |
|---|---|
| Rotate | Left drag |
| Pan | Right drag |
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
- **Edit Structure...** — modify lattice vectors and atom list (add/edit/delete)
- **Atomic Sizes…** — adjust per-element covalent radii (literature defaults: Cordero et al., *Dalton Trans.* 2008)
- **Element Colors…** — override CPK colours per element and tune material shininess
- **Transform Structure…** — apply a 3×3 matrix transformation to all atom positions (only available when the structure has a unit cell)

### Measurements and overlays

- Distance and angle tools open result dialogs and draw dashed helper lines.
- **Atom Info** shows element name, atomic number, Cartesian coordinates, and direct coordinates (when a unit cell exists).
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
