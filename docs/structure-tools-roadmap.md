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

## Beyond this round

Online structure databases (COD, Materials Project) need network access and
accounts; per-atom magnetic moments and charges for input files need new
per-atom data in projects. They are recorded here and left for a later round.
