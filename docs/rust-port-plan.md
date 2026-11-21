# AtomForge: port to Rust

Status: phase 0 done (toolchain, renderer, UI toolkit chosen by benchmark); phase 1 in progress. This document is the contract for the
port: what moves where, in which order, and what has to be true before each
piece replaces its C++ original.

## 1. What is being ported

| C++ module | Lines | Contents | Rust crate |
|---|---:|---|---|
| `src/model`, `src/math`, `src/util` | 2,450 | `Structure`, `AtomSite`, elements, colours, neighbour search, workspace, display settings | `af-core` |
| `src/io` | 4,400 | XYZ/extXYZ, POSCAR, CIF, PDB, LAMMPS data/dump, XSF, GEN, CELL, PWI, CP2K, trajectories | `af-io` |
| `src/algorithms` | 15,900 | bulk, SSS, SQS, stacking faults, CSL, interfaces, nano/poly/amorphous, Shape Modeller, dislocations, CNA/RDF/ADF/SRO/Voronoi | `af-algorithms` |
| `src/science` | 15,300 | analyses, potentials, MD/relaxation, phonons, diffraction, symmetry, batch, workflows, surface meshes | `af-science` |
| `src/electronic` | 3,200 | volumetric data, Bader, Wannier, LOBSTER, topology, Fourier analysis | `af-electronic` |
| `src/chemistry` | 6,100 | molecular graphs, SMILES, UFF, embedding, point groups, organic builders and edits | `af-chemistry` |
| `src/pipeline` | 6,400 | modifiers, expressions, text syntax, all build/edit/analysis steps | `af-pipeline` |
| `src/cloud` | 1,000 | out-of-core atom clouds | `af-cloud` |
| `src/cli` | 4,300 | the `AtomForge` command line (build, analyse, pipeline, render, cloud) | `af-cli` (binary) |
| `src/graphics`, `src/camera` | 11,000 | OpenGL renderer: impostor spheres, cylinders, AO, shadows, picking, overlays, offscreen capture | `af-render` (wgpu) |
| `src/app`, `src/ui` | 44,900 | GLFW window, Dear ImGui menus, 50+ dialogs, undo, file service, interaction | `af-app` (egui/eframe, see section 3a) |
| `python/atomforge` | 17,700 | Python package over `atomforge_science_native` and `atomforge_electronic` | unchanged; loads the Rust `cdylib` |

About 115,000 lines of C++. The Python package stays Python: it talks to the
native libraries through a small C ABI (JSON in, JSON out), which the Rust
`af-ffi` crate re-exports under the same names, so the 220-test Python suite
becomes the first full regression test of the Rust core.

## 2. Workspace layout

```
rust/
  Cargo.toml                 workspace
  crates/
    af-core/                 Structure, elements, math, neighbours, JSON model
    af-io/                   file formats
    af-algorithms/           builders and structure analyses
    af-science/              science tools, potentials, simulations, symmetry
    af-electronic/           volumetric and electronic-structure tools
    af-chemistry/            organic and coordination chemistry
    af-pipeline/             modifiers and steps
    af-cloud/                atom clouds
    af-render/               wgpu renderer (also headless, for the CLI)
    af-ffi/                  cdylib: atomforge_science_native, atomforge_electronic
    af-cli/                  binary: atomforge
    af-app/                  binary: AtomForge (egui/eframe on wgpu)
  tests/
    golden/                  reference outputs produced by the C++ build
    fixtures/                input structures and volumes
```

## 3. Dependency mapping

| C++ | Rust |
|---|---|
| glm, Eigen | `glam` (f64 vectors, 3x3) and `nalgebra` (dense eigen/SVD) |
| OpenMP | `rayon` |
| JSON (own) | `serde`, `serde_json` |
| zlib, bzip2, xz | `flate2`, `bzip2`, `xz2` |
| spglib | `af-spglib`: bindings to spglib 2.7.0 itself (the version the C++ links), the MSVC build vendored from the spglib wheel (BSD-3), `symspg.dll` shipped next to the executable; moyo was tried first but its primitive standardisation differs from spglib (origin, atom order, wrapping) |
| Open Babel (file formats in the desktop app) | native MOL/SDF/MOL2/GJF readers and writers in `af-io` |
| GLFW, GLEW, OpenGL | eframe window and `wgpu` (the viewport renders on the UI's own device) |
| Dear ImGui | egui (immediate mode, like Dear ImGui) |
| stb (images), GIF writer | `image`, `gif` |
| Python script step | same approach: spawn the Python interpreter |

## 3a. UI toolkit: measured, GPUI rejected

GPUI was to be used only if it performed better than the alternatives. The
viewport is the demanding part of this application, so phase 0 measured the
same af-render scene (103,823 impostor spheres, rotating camera, NVIDIA RTX
2000 Ada, vsync off) in GPUI 0.2.2 and in egui/eframe 0.37 (both on wgpu 30):

| Window | GPUI 0.2.2 | egui/eframe 0.37 |
|---|---:|---:|
| 640 x 372 | 716 fps | 1,801 fps |
| 1280 x 772 | 236 fps | 734 fps |
| 1900 x 1122 | 127 fps | 395 fps |

GPUI has no way to show a wgpu texture: every frame must be read back to the
CPU and uploaded again as an image (on Windows into its DirectX atlas), so it
is 2.5 to 3.1 times slower and the gap grows with the window. eframe renders
the viewport on its own wgpu device and shows the texture directly. With
vsync both run at the display rate (32 Hz over Remote Desktop here). The UI is
therefore egui/eframe; its immediate-mode model also matches the existing
Dear ImGui code, which keeps the 45,000-line UI port a close translation.
(GPUI 0.2.2 also paces frames from its own vsync thread with no option to
turn it off; the GPUI figures used a benchmark-only patch of that wait.)

## 4. How each piece is verified before it replaces the C++ one

1. **Golden outputs.** A small C++ program (`tests/golden_dump.cpp`), built
   against the existing libraries, writes JSON reference results for every
   public function on the fixture set: builders' atoms and cells, analyses'
   per-atom values, file round trips, pipeline step outputs and notes. The Rust
   tests read the same fixtures and compare within stated tolerances (1e-9
   relative for deterministic geometry, statistical bounds for random
   builders, exact for integer classifications).
2. **Unit tests** ported from `tests/*.cpp` (30 regression files) to Rust
   `#[test]`s, one file per C++ test file.
3. **Python suite** (`python/tests`, 220 tests): run against the Rust `cdylib`
   once `af-ffi` exists; it must pass unchanged.
4. **CLI smoke tests** (`tests/cli_smoke.py`, `render_cli_smoke.py`): run
   against the Rust binary with the same arguments.
5. **Rendering:** headless renders of reference scenes compared with the
   OpenGL renders by image difference (mean absolute error threshold).
6. **UI:** egui `egui_kittest` tests for each dialog's state logic, plus the
   UI tour (scripted screenshots of every menu and dialog) re-implemented for
   the Rust app.

A crate is "ported" only when its golden comparisons, ported unit tests and
the relevant Python/CLI tests pass on Windows and Linux. Until the final
switch-over the C++ build stays as the reference and keeps working.

## 5. Phases

| Phase | Work | Exit criterion |
|---|---|---|
| 0 | Toolchain (rustup, MSVC target), workspace skeleton; viewport spike of 10^5 spheres in GPUI and in egui, benchmarked; UI toolkit chosen | done: egui/eframe (section 3a) |
| 1 | `af-core`: Structure, elements, colours, cell maths, neighbour search, workspace/display-settings JSON | golden + ported `core_regressions`, `display_settings`, `workspace` tests |
| 2 | `af-io`: every reader/writer, trajectories | round trips and golden parses of all fixtures; `simulation_formats` tests |
| 3 | `af-pipeline` framework (expressions, text syntax, parameters, registry) with the selection and basic edit steps | `pipeline_regressions` subset |
| 4 | `af-algorithms`: builders, then structure analyses | golden per builder/analysis; `builder_tools`, `analysis`, `pattern_match`, `dislocation*`, `nye_tensor`, `lattice_family` tests |
| 5 | `af-science`, `af-electronic`, `af-cloud` | golden; `science_tools`, `electronic`, `bader`, `lobster`, `wannier`, `topology`, `cloud` tests |
| 6 | `af-chemistry` and the remaining pipeline steps | golden; chemistry tests (CIP, point groups, UFF gradients, Open Babel comparisons) |
| 7 | `af-ffi`: the C ABI of both native libraries | full Python suite passes against the Rust libraries |
| 8 | `af-cli` | CLI smoke tests and golden CLI outputs |
| 9 | `af-render`: scene buffers, impostors, cylinders, cell, AO, shadows, depth cueing, outlines, picking, headless capture | image comparisons; `render_cli_smoke` |
| 10 | `af-app`: window, menus, viewport, interaction, undo/redo, file service | UI tests; manual check list |
| 11 | dialogs, one at a time (builders, edits, analyses, pipeline editor, settings) | per-dialog UI tests and tour screenshots |
| 12 | switch-over: packaging, manual and docs updated, C++ sources removed; the Rust autosave back to `recovery.afproject` | all of the above green on Windows and Linux |

## 6. Risks

- **egui** changes its API between minor versions (0.37 replaced `App::update`
  with `App::ui` and the panel types); the workspace pins one version.
- **wgpu dependency resolution:** `gpu-allocator` may resolve against an older
  `windows` crate than wgpu's DX12 backend; the lockfile pins it to the one
  wgpu uses.
- **Numerical parity:** floating-point order differs (rayon reductions,
  different eigen solvers). Tolerances are stated per test; random builders
  are compared statistically, not atom by atom.
- **Open Babel replacement:** formats it read (MOL/SDF/MOL2/GJF) get native
  readers; exotic formats it covered are listed in the switch-over notes.
- **Scale:** about 115,000 lines. Porting proceeds crate by crate, with the C++
  build kept working until phase 12.

## 6a. Parity policy

Everything is compared bit for bit with the C++ except where an Eigen algorithm is not reproduced: dense (complex) eigen-solvers and least-squares fits use nalgebra and are compared to 1e-9 relative plus 1e-10 absolute (the physical results, such as a dislocation's displacement field, are independent of eigenvector scaling and root order). Phonons are compared to 1e-5 relative plus 1e-6 absolute: their Gamma acoustic modes are square roots of rounding noise (about 1e-7 THz) and the lowest of them also places the DOS energy grid. Eigen's closed-form 3x3 inverse and determinant, its 3x3 product summation order (af_core::eigen::eigen_matmul3, eigen_mul3), JacobiSVD (jacobi_svd3) and the 3x3 SelfAdjointEigenSolver (self_adjoint_eigen3) are reproduced exactly; each has its own golden test against Eigen 3.5.

Distance-geometry embedding (af-chemistry `embed`, SMILES to 3D) diagonalises an n x n metric matrix with Eigen's dynamic SelfAdjointEigenSolver (blocked Householder tridiagonalisation with alignment-dependent SIMD reductions), which is not reproduced; Rust uses nalgebra, so the random starting coordinates differ and the refined geometry is another (equally valid) embedding. Embeddings, and the builders that embed fragments (polymers, complexes, functional-group substitution), are therefore compared by everything that does not depend on the coordinates: the molecular graph, formula, charges, stereo descriptors, the count of unspecified stereocentres, the stereochemistry of the result and the success count, exactly; and the UFF energy of the kept embedding after minimisation, which may be at most 2 kcal/mol (or 10 %) higher than the C++ one, or up to 5 kcal/mol (25 %) lower (another conformer of a flexible molecule). In the golden set every optimised embedding but one reaches the same minimum as the C++ to about 1e-4 kcal/mol; the flexible hexitol settles 1.9 kcal/mol lower.

## 7. Progress

| Phase | State |
|---|---|
| 0 | done: Rust 1.100 (MSVC), workspace, `af-render` sphere impostors (headless test), eframe shell at 724 fps uncapped |
| 1 | done: `af-core` (structure, elements, cell maths, neighbour search with C++ pair order, display settings and looks, ambient occlusion, picking, C-compatible number formatting, libstdc++ `std::sort` order, task cancellation); 15 tests incl. the ported display_settings and core_regressions checks |
| 2 | mostly done: `af-io` readers and writers for every structure format (Open Babel formats now native), CIF space groups from spglib's exported table, partial occupancy, IPF sidecars and geometry reconstruction, simulation-code inputs, XYZ and LAMMPS trajectories; golden comparison over 92 file loads, 105 saves and 9 trajectory cases, the ported simulation_formats test, Open Babel reading the port's MOL/SDF/MOL2/PDB; project files wait for `af-electronic` (volumes) |
| 3 | done: `af-science` foundations (row-vector cell algebra, general minimum image, periodic neighbour list with the C++ pair order, FFT, numeric arrays; frame readers for extXYZ/POSCAR/XDATCAR/LAMMPS/JSON, streamed trajectories, NPY, text tables, EIGENVAL with gzip/bzip2/xz, tool `Parameters`, the 62-tool catalog generated from the C++ source); `af-pipeline` (expression language, text syntax, JSON form, registry with category order, cached evaluation, the 28 selection, modification and cell steps of Modifiers.cpp). Golden: science.json (algebra, minimum images and neighbour lists on 12 random cells, FFTs, array checks, 40 data fixtures through `readFrames`, `TrajectoryStream` and `Parameters`), pipeline.json (40 texts, 15 JSON forms, 55 expressions x 4 inputs, 39 pipelines x 4 inputs including notes and stage results), all identical; 9 ported pipeline_regressions tests. insert-interstitials waits for `af-algorithms`, the relax/dynamics/standardize steps for the science tools. Fixed on the way: number parsing in af-core was quadratic in the text length (15 s -> 10 ms on a 1.5 MB JSON file). |
| 4 | done (the Python script step waits for phase 7): `af-algorithms`: interstitial voids, bulk crystal builder (230 space groups, both origins), CNA, RDF, ADF, Warren-Cowley/Rao-Curtin SRO, Voronoi cells, substitutional solid solutions, amorphous packing, cell sculptor, Shape Modeller scripts (meshing waits for af-electronic's isosurface), OBJ/STL meshes, elastic tensors, strain, vacancies, SQS annealing, drag interpolation, nanowires, core-shell, pattern matching, anisotropic (Stroh) dislocation fields, stacking-fault sequences, dislocation insertion (edge/screw/mixed, five shapes, dipoles, anisotropic cubic/hexagonal/general), Vitek maps, Nye tensors and Nye fits; af-core supercells; symmetry through spglib 2.7.0 itself (af-spglib, af-science::spg); CSL tables and bicrystals; primitive reduction; surface slabs; interface lattice matching and assembly; Voronoi polycrystals; nanocrystals (7 shapes, mesh fill, Euler/Miller orientation, Wulff polyhedra). Golden over 10 inputs, all identical: 30 void detections (3022 sites), 268 bulk builds, 30 CNA, 40 RDF, 30 ADF, 20 SRO, 20 Voronoi, 24 solid solutions, 10 amorphous packings, 42 sculptor cuts, 136 shape scripts, 8 mesh files, 58 small-tool cases, 15 stacking-fault sequences (FCC, BCC (110)/(112), HCP basal/prismatic/pyramidal, diamond, CNA fallback), 20 dislocations (the 4 anisotropic ones and 9 Nye/fit runs to 1e-9 relative, being Eigen-based), 9 Vitek maps, 8 CSL tables, 10 CSL bicrystals, 10 primitive reductions, 15 surface slabs, 9 interface searches (slabs, every match, 3 assemblies each) and 7 helper cases, 11 polycrystals, 20 nanocrystals with Wulff previews. Pipeline families ported: StructureEdits (bonds, rotations, saturation), EditSteps (18 steps), ShapeSteps (build-shape), GeneratorSteps (12 builders) and the molecule library (28 molecules bit-identical); golden runs on four inputs: 79 edit, 17 shape and 89 generator pipelines, all identical. New exactness helpers: af_core::cxx::stable_sort_by (libstdc++ merge sort, also for non-strict tolerance comparators; 208 golden cases) and af_core::crmath (correctly rounded sin/cos, as GCC constant-folds them with MPFR; 3000 golden values). AnalysisSteps (14 steps: overlap, surface and region selections, colouring, ambient occlusion, per-atom properties from 18 analyses, expressions, elastic strain, displacements, binning, statistics, grain segmentation) and the Simulation steps (relax, NVT, NPT, standardize-cell through the catalog tools): 57 more golden pipelines on four inputs, identical; build-shape writes its OBJ surface (af_algorithms::shape Model::mesh/save_obj on the ported isosurface), 3 files byte-identical. **C++ fix found by the comparison:** StackingFaultBuilder, NanoCrystalBuilder (Wulff point group) and StructureInfoDialog passed spglib the cell vectors as rows, but spglib takes columns: wrong for every non-symmetric cell matrix (HCP was reported FCC-like, the Wulff holohedry and the displayed space group of hexagonal cells were wrong). Fixed in the C++ too. |
| 5 | done: `af-science::tools` dispatcher (catalog validation with did-you-mean, reports, per-atom requests and properties) and all 62 catalog tools: the 18 of Analysis; trajectory average, unwrap, frame difference, time series and trajectory structure (g(r), coordination, bond angles); 5 property analyses; 5 local-environment tools; clusters and voids; bond orders and LAMMPS regions; surface meshes; Voronoi/VoroTop/CNA; template matching (PTM, grains, elastic strain, planar faults); dislocation lines (Burgers circuits); symmetry (spglib, Wyckoff orbits, symmetrized cells); powder XRD and electron diffraction; HPKOT reciprocal paths and DFT inputs (KPOINTS, POSCAR, QE cards); LAMMPS export; VASP electronic (EIGENVAL, DOSCAR, KPOINTS, PROCAR); the native potentials (EMT, Lennard-Jones, EAM funcfl/setfl/fs, Tersoff, Stillinger-Weber, Buckingham with Ewald) and relax, NVT, NPT, NEB, phonons, elastic constants, EOS scans, formation and planar defect energies. `af-electronic` done (the C ABI of its *API.cpp waits for phase 7): grids, Cube/VASP/XSF reading and writing, grid algebra and resampling, integration (spheres, Voronoi, peaks), gradient/Laplacian/energy densities and RDG/sign(lambda2)rho/DORI, smoothing, profiles, sections and contours, FFT structure factors, Fourier synthesis, Patterson maps, atomic structure factors, Ewald potentials, charge-difference masks and cumulative charge, on-grid Bader, Betti curves, LOBSTER COHPCAR/ICOHPLIST, Wannier90 bands, Berry curvature and Chern numbers, marching-tetrahedra isosurfaces. Golden: 146 cases over 30 synthetic fixtures (tools/golden/make_electronic_fixtures.py, golden_electronic.cpp) and 6 saved volume files, all bit-identical except the Wannier90 diagonalisations (nalgebra, 1e-9); GCC-folded pow constants are stored as their correctly rounded bits. Golden: 335 requests (results with reports, frames and files, or error messages; synthetic EAM/Tersoff/SW/VASP fixtures), generated by tools/golden/make_tool_requests.py and golden_tools.cpp (built by tools/golden/build.sh); all bit-identical except the 6a classes (11 Eigen-fit tools to 1e-9, phonons to 1e-5). Relaxations, MD trajectories and NEB chains match the C++ bit for bit. `af-cloud`: .afcloud reading, the bucketed out-of-core builder (libstdc++ nth_element, mt19937_64 and shuffle reproduced), LAMMPS-dump/XYZ converters and the crystal generator; golden: 37 cases and 11 cloud files byte-identical with one build thread, plus a test that multi-threaded builds hold the same chunks and atoms. |
| 6 | done: `af-chemistry`: molecular graphs (perception from coordinates with periodic bonds, bond orders and formal charges, rings, aromaticity, Kekule structures, hybridisation, hydrogen placement, canonical ranks, CIP R/S and E/Z), SMILES reading and writing, the UFF force field (127 atom types, bond/angle/torsion/inversion/vdW/Coulomb terms with periodic images, restraints, L-BFGS with fixed coordinates), distance-geometry embedding, point groups (Schoenflies classification, symmetrisation, standard orientation), editing (functional-group substitution, bond orders with hydrogens, (de)protonation, ionisation at pH 7, stereocentre and double-bond inversion, dihedrals, conformer searches, geometry scans, Gasteiger charges, Kabsch RMSD), builders (217 named compounds, peptides from sequences, polymers with tacticity and periodic chains, coordination complexes, Z-matrix reading and writing, ORCA/NWChem/Psi4/Q-Chem/MOPAC/Turbomole/Gaussian/GAMESS inputs). Golden: chemistry.json (91 cases: 80 SMILES with properties and explicit hydrogens, 10 perceptions, element tables) and chemistry_ops.json (307 cases: UFF energies, gradients, terms, restraints and minimisations on 11 molecules incl. a periodic chain; 18 point-group cases; 67 editing cases incl. conformer searches and scans; 12 peptides, 12 polymers, 19 complexes, 16 Z-matrices read and 6 written, 80 quantum inputs; 50 embeddings), all bit-identical except the embedding class of 6a. The occasional crash of the C++ golden dumpers was the shell's PATH (Git's libstdc++ before MSYS2's); build.sh now puts the UCRT64 runtime first. |
| 7 | done: `af-ffi-science` (cdylib `atomforge_science_native`: afs_version/free/catalog/run, afs_pipeline_modifiers/run/parse/text) and `af-ffi-electronic` (cdylib `atomforge_electronic`: the volume handles and 40 operations of PythonAPI, Bader, LOBSTER COHPCAR/ICOHPLIST and Wannier90 APIs, with thread-local error strings; panics never cross the boundary); the `python-script` pipeline step (af-pipeline steps/python.rs, registered between shape and science steps as in C++). The Python suite (`rust/tools/run_python_suite.sh`) against the Rust libraries: 220 tests, 190 pass, 29 skipped (optional packages, external fixtures, or the AtomForge executable), 1 error from `test_cli_rejects_unknown_options`, which needs `ATOMFORGE_PATH` (phase 8) and does not skip without it. spglib ships as `symspg.dll` next to the libraries. |
| 8 | done: `af-cli` (binary `atomforge`): the dispatcher, all 17 `--build` modes, `--analyze` (13 analyses, sculptor, drag) and `--convert`, `--science` and `--science-batch` (af-science::batch ported from science/Batch.cpp), `--pipe`/`--pipeline` (lists, JSON, stdin/stdout extended XYZ, reports), `--cloud build|generate|info`, the `build-*` pipeline steps of BuildModifiers (modes run in process with their output captured), help pages generated from the C++ sources (tools/gen_cli_help.py) and the menu-location table (tools/gen_cli_menu.py). Comparison with the C++ executable (tools/compare_cli.py: 42 cases, about 140 commands, exit codes, stdout, stderr and every written file byte for byte, timings masked): identical, except the Nye-tensor and dislocation-fit CSVs (Eigen least squares, compared to 1e-9 as in 6a) and PDB/MOL2 written by Open Babel in the C++ (native writers in Rust; listed for the switch-over notes). Found on the way: the orientation sidecar must be written in text mode (CRLF on Windows) as the C++ does. `tests/cli_smoke.py` and `tests/render_cli_smoke.py` pass with the Rust binary; the Python suite with the Rust executable and both Rust libraries: 220 tests, 204 pass, 16 skipped (optional packages, external fixtures), none fail. |
| 9 | done for the scene and offscreen rendering (the desktop's interactive extras come with phases 10-11): `af-render::scene` (instances with periodic boundary images and integer supercell transforms, cell box, bonds by spatial hash in the C++ order, ambient-occlusion darkening, the default framing, frame view matrices with glm's formulas) bit-identical to the C++ (render.json: 6 structures x instances, 2 supercell transforms, AO, 3 bond styles and 9 views; tools/golden/golden_render.cpp); `af-render::structure` (wgpu): the GLSL ported to WGSL (40x40 sphere meshes, 32-sided two-coloured bond cylinders, 2048x2048 PCF shadow map with front-face culling and polygon offset, six-light shading, saturation and contrast, depth cueing, outlines, GL-style wide cell lines, depth clamping); `af-science::gif` (GifWriter and reader) byte-identical to the C++ (tools/golden/golden_gif.cpp); `atomforge --render` (PNG with optional pHYs DPI, turntables, trajectory and turntable GIFs, display-settings files, colour and radius overrides). Images compared with the C++ executable (tools/compare_render.py, 9 scenes incl. display settings with depth cueing, outlines and AO): mean difference 0-0.9 levels, at most 0.4 % of pixels off by more than 24 levels (cell lines and bond edges rasterise slightly differently), three scenes pixel-identical. Not yet ported (desktop only): low-poly and billboard modes for huge structures, GPU culling, the atom-cloud renderer, the electronic-structure viewport, selection halos, gizmo and SVG export. |
| 10 | done: `af-app` (binary `atomforge-rs`, egui/eframe on wgpu): structures in tabs, menus generated from the shared menu table (every pipeline step and catalog tool in its C++ location), the af-render viewport drawn on the UI device and shown as a texture, the C++ camera controls (orbit, pan, zoom with decay, look along axes and lattice vectors), picking, undo/redo with the ported UndoRedo, the generic operation dialog (All / view selection / condition targets), the scientific-tool window, file browser, recent files, unsaved-changes prompt, image export; View menu (atom styles, colour by element / crystal orientation / grain boundaries / atom property, labels, element legend, orientation axes, projection), Display Settings and Atomic Sizes (looks, lighting, sizes, bond rules, element radii and colour schemes, saved defaults), Structure Info (spglib in the background), Atom Info and distance/angle/dihedral measurements with result cards, Find a Command (C++ ranking), Keyboard Shortcuts, interface size and theme; projects (Open Project, Save Project As, autosave every minute, Recover Last Autosave) through `af_electronic::project`, which reads and writes the C++ formats 1-7 byte for byte (`golden_project.cpp`; the C++ loader reads the Rust app's projects; the Rust autosave uses `rust-recovery.afproject` until phase 12 so it never overwrites the C++ desktop's recovery file); `af_render::colouring` and `af_render::measure` compared with the C++ (`golden_colouring.cpp`, `golden_measure.cpp`: grain-boundary flags and 102 measurement texts identical); the view toolbar (look along axes and lattice vectors, Rotate as the C++ camera orbit, bit-exact over 432 cases), status bar, notices, hover tooltip, box and lasso selection (`af_render::region`), the selected atoms' menu (substitute, insert at midpoint, measure, delete) with S and I, grab mode G with axis locks and periodic wrapping (`af_render::grab`), Help > Manual; UI scripts (`ATOMFORGE_UI_SCRIPT`, `rust/tests/ui/*.afscript`: editing, view, measure, palette, project, select, grab) drive and check it |
| 11 | done: Edit > Structure Pipeline (`af_pipeline::editor::PipelineEditor`, the C++ PipelineEditor with its regression tests; the panel with per-step results, bypass, drag/arrow/Alt reordering, replace and insert menus, the step form shared with the operation dialog, Edit as Text, Copy as Command, pipeline files, colour by selection, Bake/Discard, kept in projects; UI test `pipeline`); result plots (`af_science::plots`, the C++ science/ResultPlots: every plot, CSV and overlay of the 335 golden tool results identical, `golden_plots.cpp`) drawn by the ported SciencePlot widget in the tool window, with plot tabs, Copy CSV, comparison with earlier runs, per-atom values coloured onto the atoms, and the open structure supplying positions, cell, periodic axes and species (UI test `tools`); Build > Bulk Crystal (system, space group with both origins, free lattice parameters, asymmetric-unit table with an element picker, Copy as Pipeline Step that the CLI runs; UI test `builders`); Analysis > Structure and Defects dialogs: Common Neighbour Analysis (signatures, environments, per-atom table), Radial Distribution Function (species pair, PBC, g(r) or counts, sampling, overlays, first-shell distortion metrics), Angular Distribution Function (centre/triplet selection, smoothing, reference geometries, peaks, coordination, guide) and Short-Range Order (Warren-Cowley and Rao-Curtin per shell with the alpha colour scale, shells, composition, guide), all run in the background (UI test `analysis`; SRO values equal the CLI's); Trajectory Playback (background indexing and frames read on demand, transport and slider with the ported TrajectoryPlayback, LAMMPS type names, per-frame analysis colouring, unwrapped trajectory lines over the view, GIF export rendered offscreen with the current camera; UI test `trajectory`); View > Lattice Planes and Miller Directions overlays with their settings windows (cell-clipped (hkl) polygons, [uvw] arrows; UI test `overlays`); Voronoi volume cells and dislocation loop lines (View > Show); polyhedral atom style and Polyhedral Viewer with Polyhedral Settings (`af_render::polyhedra`); Edit > Insert Dislocation (lattice detection, character, region shape, manual directions, line point, isotropic or Stroh elasticity, dipoles, background generation, result summary, Replace Main Scene with the loop shown; UI test `dislocation`); a shared `StructureDialog` interface (`af-app/src/structure_dialog.rs`) through which the remaining dedicated dialogs plug into menus, palette and `dialog STEP ...` script commands; the system font (Segoe UI / DejaVu Sans, as the C++) for unit symbols. Every other operation and tool is reachable through the generic operation dialog and the tool window; dedicated dialogs follow one by one; every remaining dialog through the `StructureDialog` interface: Build > Substitutional Solid Solution, Stacking Fault, CSL Grain Boundary, Interface Builder, Nanocrystal, Custom Structure, Polycrystal, Amorphous Structure and Shape Modeller; Edit > Edit Structure, Transform, Merge, Cell Sculptor and Add Interstitial Atoms; Analysis > Electronic Post-Processing, Wannier, LOBSTER and Dislocation Post-Processing (UI tests `grain_builders`, `nano_poly_amorphous`, `composition_builders`, `edit_dialogs`, `interstitial_sculptor`, `electronic`, `bonding`, `dislocation_analysis`); Edit in Dialog on every pipeline step with each dialog's editStep/stepOptions ported (`af_pipeline::options`, the C++ Options.h; option strings identical to the C++ dialogs in golden comparisons; the `*_steps` UI tests); dropped files routed to open dialogs, projects keeping the electronic workspace and the tool window; the scientific-tool window at C++ parity (layouts, inputs, potentials, comparison, saving, surface meshes over the view); the large-data cloud viewer (`af_render::cloud`, streaming, level of detail, budgets, occlusion queries; UI test `cloud`); Export Image as PNG, JPEG and SVG with scale, background and gizmo (UI test `export`). The interface matches the C++ one: the C++ UI tour (`ATOMFORGE_UI_TOUR`) and its Rust counterpart `rust/tools/ui_tour.afscript` give the same 76 screenshots, compared side by side by `rust/tools/compare_tour.py`: the C++ light and dark palettes, metrics and fonts (`theme.rs`, the Codicons icon font), toolbar, tabs, menus with check marks, welcome panel, modal dialogs with dimming and minimise, notifications, selection halos, legends, the Atom Property Colouring panel and every dialog's layout (`dialog_layout.rs`, `dialog_widgets.rs`) |
| 12 | done: the C++ sources removed (src, imgui, cmake, CMakeLists.txt, the C++ tests; CI builds and tests the Rust workspace on Linux and Windows); one program, `AtomForge`, for the desktop and the command line (`--renderer`, wildcards, Settings > Renderer and `ATOMFORGE_RENDERER`, the icon, manifest and version resource); `rust/tools/package.py` builds the portable packages (tested from a fresh extraction: desktop, command line and the Python package); INSTALL.md, ARCHITECTURE.md and the manual updated for the Rust build; the formats Open Babel used to read are read natively (CIF with Cartesian coordinates, Quantum ESPRESSO input with any `ibrav`); the autosave is `recovery.afproject` and the recent files are shared with the C++ desktop |
