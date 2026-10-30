# Structure builder and editor: feature audit

Audited 2026-10-30 against the building and editing tools of ASE (`ase.build`),
Atomsk, pymatgen transformations, VESTA, Avogadro, OVITO and Materials Studio.
AtomForge already has: bulk crystals from space groups, CSL grain boundaries,
polycrystals, nanocrystals (sphere, ellipsoid, box, cylinder, octahedron,
cuboctahedron, Wulff), amorphous packing of atoms, solid solutions, SQS,
dislocations, interfaces, stacking faults, mesh-filled shapes, vacancies,
strain, primitive cells, surface slabs, nanowires, core-shell particles,
interstitial insertion, the cell sculptor, replicate/supercell/affine
transforms, slicing, selections, undo, and CIF, POSCAR, LAMMPS data, XYZ,
PDB, Gaussian, mol2 output (more through Open Babel).

Status: **done** (implemented and tested), **open**.

Every item is a pipeline step (desktop menus, `AtomForge --pipe`, Python
`Pipeline`, project files) or a native file format, and is tested in
`python/tests/test_structure_tools.py` (against ASE: nanotubes, monolayers,
nanoribbons, icosahedra, cuboctahedra, adsorption sites, Niggli reduction,
near-cubic supercells; and against exact invariants) and
`tests/simulation_formats.cpp` (round trips and hand-written inputs).

## Builders

| # | Feature | Reference | Status |
|---|---------|-----------|--------|
| A1 | Nanotubes with chiral indices (n, m), single- and multi-wall, any honeycomb 2D crystal (C, BN) | `ase.build.nanotube` | done |
| A2 | 2D monolayers: graphene, h-BN, 2H/1T transition-metal dichalcogenides, with vacuum | `ase.build.graphene`, `mx2` | done |
| A3 | Graphene/h-BN nanoribbons, armchair and zigzag, hydrogen-saturated edges | `graphene_nanoribbon` | done |
| A4 | Twisted bilayers: commensurate moiré supercells for a twist angle | moiré builders | done |
| A5 | Magic-number clusters: Mackay icosahedra, regular decahedra, cuboctahedra by shell count | `ase.cluster` | done |
| A6 | Molecule library (water, CO2, methane, benzene, C60, ...) | `ase.build.molecule`, Avogadro | done |
| A7 | Molecular packing: N copies of a molecule in a box at a density, random orientations, minimum distance | Packmol | done |
| A8 | Adsorbates on a slab at top, bridge and hollow sites at a height | `ase.build.add_adsorbate` | done |
| A9 | Near-cubic supercell of a target size (for point-defect calculations) | `find_optimal_cell_shape` | done |
| A10 | Interpolated images between two structures (NEB start, morphs) | `ase.neb` interpolate | done |
| A11 | Antisite pairs and Frenkel pairs | pymatgen defects | done |

## Editing

| # | Feature | Reference | Status |
|---|---------|-----------|--------|
| B1 | Mirror through a plane, invert through a point | VESTA, Atomsk `-mirror` | done |
| B2 | Rotate about an axis by an angle; orient a crystal direction along x, y or z | Atomsk `-rotate`, `-orient` | done |
| B3 | Shift by a fractional vector (move the origin) | Atomsk `-shift` | done |
| B4 | Sort atoms by element, coordinate or property | Atomsk `-sort` | done |
| B5 | Saturate dangling bonds with hydrogen | Materials Studio, Avogadro | done |
| B6 | Set a distance, angle or dihedral between atoms, moving one side | Avogadro | done |
| B7 | Duplicate the selection with an offset | VESTA, Avogadro | done |
| B8 | Niggli-reduced cell | spglib, pymatgen | done |
| B9 | Scale the cell to a volume, density or lattice factor | ASE, Atomsk `-deform` | done |
| B10 | Swap two elements | Atomsk `-substitute` | done |
| B11 | Make molecules whole across periodic boundaries | OVITO unwrap, ASE | done |
| B12 | Select connected fragments (molecules, clusters) | Avogadro, OVITO | done |
| B13 | Measure dihedral angles | VESTA, Avogadro | done (measure-geometry step; the interactive Measure tools cover distance and angle) |

## File formats (write, and read where simple)

| # | Format | Status |
|---|--------|--------|
| D1 | Quantum ESPRESSO pw.x input (`.pwi`) | done |
| D2 | CASTEP `.cell` | done |
| D3 | DFTB+ `.gen` | done |
| D4 | XCrySDen `.xsf` | done |
| D5 | CP2K `&SUBSYS` (`.cp2k`) | done |

## Per-atom simulation data, stacking, measurement and databases

A second pass over the same tools found these further gaps.

| # | Feature | Reference | Status |
|---|---------|-----------|--------|
| E1 | Initial magnetic moments per atom: set by selection (antiferromagnetic orders), kept in projects and edits; extended XYZ `initial_magmoms`, CASTEP `SPIN=`, Quantum ESPRESSO species per moment with `nspin = 2`, VASP `MAGMOM` line in the notes | ASE, VESTA, pymatgen | done |
| E2 | Charges per atom: LAMMPS data in the charge style (read and written), extended XYZ `initial_charges`, total charge reported | ASE, Atomsk | done |
| E3 | Maxwell-Boltzmann velocities at a temperature (no net momentum, exact kinetic temperature): LAMMPS `Velocities` (read and written), extended XYZ `velo` | ASE `MaxwellBoltzmannDistribution`, LAMMPS `velocity create` | done |
| E4 | Stacking layers into heterostructures and multilayers (gap, in-plane strain to the substrate reported) | Atomsk `-merge`, heterostructure builders | done |
| E5 | Interactive dihedral measurement (View > Measure Dihedral, the atom context menu, the command palette) | VESTA, Avogadro | done |
| E6 | Downloads from the Crystallography Open Database and the Materials Project (with its magnetic moments), as `atomforge.databases` and `python -m atomforge.databases` | ASE, pymatgen `MPRester` | done (Python; the desktop opens the downloaded file) |
| E7 | Expressions: a remainder operator (`index % 2` for sublattices) and the variables `magmom`, `charge`, `vx`, `vy`, `vz` | OVITO expressions | done |
| E8 | Fixes found on the way: Python's XYZ writer always wrote `pbc="T T T"`; the Python bridge dropped fixed (selective dynamics) flags; Atom Info could write past its text buffer | | done |

Tests: `python/tests/test_structure_data.py` (antiferromagnetic orders, data through
edits, neutrality, the Maxwell-Boltzmann distribution and temperature, Python and
ASE extended XYZ, C++ LAMMPS and XYZ round trips, CASTEP and QE output, a graphene
on MoS2 heterostructure, database downloads with the network mocked),
`tests/simulation_formats.cpp`, `tests/workspace_regressions.cpp` (project format 7)
and `tests/pipeline_regressions.cpp` (the new operator and variables).

## Not planned

Desktop-side downloads would need an HTTP client in the C++ application; the
Python module covers them and the desktop opens the result. Polymer and
protein builders, isotopes and per-atom masses are outside the scope of a
crystal and nanostructure editor.
