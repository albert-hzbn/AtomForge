<p align="center">
  <img src="assets/icon/atomforge-icon.svg" alt="AtomForge logo" width="80" />
</p>

<h1 align="center">AtomForge</h1>

<p align="center">Build atomic structures. Explore materials. Analyse electronic fields.</p>

<p align="center">
  <a href="https://github.com/albert-hzbn/AtomForge/releases/latest">Download AtomForge</a> ·
  <a href="docs/manual/AtomForge-manual.pdf">User manual</a> ·
  <a href="https://github.com/albert-hzbn/AtomForge/releases/tag/v0.3.0">What's new</a> ·
  <a href="https://pypi.org/project/atomforge-py/">Python package</a>
</p>

AtomForge is a desktop application for creating, editing and exploring atomic structures for materials research. It brings structure building, interactive 3D inspection and electronic post-processing into one workspace, helping researchers prepare inputs and examine results for molecular dynamics and first-principles studies.

## Get started

Download the package for your platform from the [releases page](https://github.com/albert-hzbn/AtomForge/releases/latest), extract the entire archive, and open the application. Keep the extracted folders and libraries together.

| Platform | What you need |
| --- | --- |
| **Windows x64** | Open AtomForge.exe in the extracted AtomForge folder. Runtime libraries are included; no compiler or MSYS2 installation is needed. |
| **Linux x86_64** | Open AtomForge in the extracted AtomForge/bin folder. The current build targets Ubuntu 24.04 or compatible systems with glibc 2.39 or later, X11/XWayland and working OpenGL drivers. |

Each package includes the illustrated manual. Open it through **Help > Manual** or **Help > About > Manual**, using your default PDF viewer. You can also [read it here](docs/manual/AtomForge-manual.pdf).

For a first session, open an existing structure through **File > Open**, or choose a builder from **Build**. Rotate and zoom the view to inspect the result, then save your structure or export an image. The manual includes 62 short tutorials with examples to try.

## Build the structure you need

| Tool | What you can create |
| --- | --- |
| **Bulk crystals** | Periodic crystals from crystal systems, space groups, lattice parameters and atomic sites. |
| **Solid solutions** | Substitutional alloys with a chosen composition on an existing host lattice. |
| **Grain boundaries** | Coincidence-site lattice (CSL) bicrystals with controls for orientation, replication, translation and overlap removal. |
| **Nanocrystals** | Particles cut from geometric shapes, or Wulff constructions based on crystal facets and relative surface energies. |
| **Custom shapes** | Atomic structures filling an imported OBJ or STL mesh. |
| **Polycrystals** | Multiple grains arranged through Voronoi-based construction. |
| **Amorphous structures** | Random atomic packings with specified composition, density and box dimensions. |

Dislocation and interstitial tools support defect studies. Interactive merging lets you arrange and combine structures, while **Cell Sculptor** removes atoms within a chosen region to shape an existing model.

## Edit and inspect interactively

- Select individual atoms or use box and free-form lasso selection.
- Add, remove or substitute atoms, edit coordinates and lattice vectors, and transform periodic structures.
- Position and rotate structures in a 3D preview before merging them.
- Measure distances and angles directly in the scene.
- Inspect composition, unit-cell information, atomic positions and bonding information.
- Undo and redo edits as you refine a model.

## Understand structure and local environments

Explore radial distribution functions, short-range chemical order, interstitial sites and voids. Display bonds, element labels, periodic boundaries and crystal-orientation colouring to help interpret the structure.

Voronoi and coordination-polyhedron overlays reveal local environments. Lattice planes, Miller directions and dislocation outlines provide additional geometric context. Export the view as an image for reports, presentations or further editing.

## Explore electronic calculations

The **Electronic Post-processing** workspace loads charge densities and other scalar fields from VASP, Gaussian Cube and XSF files. File browsers and drag-and-drop loading make it easy to bring in primary and reference data.

**See the field from different perspectives**

- Independent 3D and 2D views with orbit, pan, zoom and reset controls.
- Density-volume and isosurface displays, adjustable transparency and density-based colours.
- Axis-aligned or arbitrary plane sections, with an optional translucent plane guide in 3D.
- Separate display levels, initial estimates from the file and warnings for out-of-range values.

**Analyse and compare fields**

- Field arithmetic, scaling, smoothing, gradients and Laplacians.
- Line profiles, planar and macroscopic averages, contours and peak searches.
- Whole-field, spherical-region and Voronoi-region integration.
- Grid structure factors, Fourier synthesis, Patterson maps and point-charge electrostatics.
- Explicit resampling, reusable result fields and exports for further analysis.

## Study charge redistribution

Compare a combined system with reference densities using weighted density differences. Separate accumulation and depletion, inspect charge summaries, and follow cumulative charge along a cell direction.

Threshold masks and Boolean operations let you define regions, combine or subtract them, apply them to densities and integrate the selected values. These tools support charge-transfer and fragment-comparison studies when the inputs use consistent geometry and grid conventions. Geometric region integrals are not Bader charges; the manual explains how to interpret the results.

## Work comfortably and share results

Light and dark themes, adjustable interface size and responsive dialogs support different screen resolutions. Atom sizes, colours and display settings can be tailored to the task.

| Data | Supported formats |
| --- | --- |
| **Desktop structures** | XYZ, CIF, PDB, SDF, MOL, VASP, MOL2, Quantum ESPRESSO PWI and Gaussian GJF |
| **Custom meshes** | OBJ and STL |
| **Electronic input grids** | VASP charge, potential and ELF files; Gaussian Cube; XSF |
| **Rendered images** | PNG, JPEG and SVG |
| **Analysis results** | Volumetric files, CSV tables and OBJ/PLY surfaces, depending on the operation |

## Python and repeatable workflows

The [atomforge-py package](https://pypi.org/project/atomforge-py/) brings structure loading, editing and notebook viewing to Python. It also connects to AtomForge's native builders and electronic calculation library for repeatable workflows. The desktop application supports automated structure generation for batch work.

The current Python release is **0.2.0**. Basic structure operations do not require third-party Python dependencies; native builders and electronic calculations need the corresponding AtomForge executable or library. Setup and worked examples are available in the [Python guide](python/README.md) and [user manual](docs/manual/AtomForge-manual.pdf).

## Help, availability and contributing

- **Learn a tool:** the [user manual](docs/manual/AtomForge-manual.pdf) explains controls, examples, interpretation and troubleshooting.
- **See release changes:** the [release notes](https://github.com/albert-hzbn/AtomForge/releases) describe new features, fixes and platform requirements.
- **Report a problem or suggest a feature:** use [GitHub Issues](https://github.com/albert-hzbn/AtomForge/issues). Include the application version, operating system and steps to reproduce the issue; attach a shareable example when possible.
- **Build or extend AtomForge:** see the [installation guide](INSTALL.md) and [architecture guide](ARCHITECTURE.md).

The standard desktop release does not expose common-neighbour analysis, angular-distribution analysis or the interface builder. The stacking-fault builder is optional and disabled in the distributed builds. These availability limits are documented in the manual.

## Citation and licence

If AtomForge contributes to your research, cite the software and record the version you used. The [Zenodo record for v0.2.0](https://doi.org/10.5281/zenodo.20054535) provides citation details and downloadable citation formats.

See [LICENSE.txt](LICENSE.txt) for the software licence.
