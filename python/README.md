# atomforge-py 0.2.0

Python tools for atomic structures, AtomForge builders, notebook viewing, and electronic post-processing.

```sh
python -m pip install --upgrade atomforge-py==0.2.0
```

Import the package as `atomforge`, not `atomforge-py`. Python 3.8 or newer is required. Basic structure I/O, editing, and notebook HTML generation have no required third-party Python dependencies.

## Native application and library

The platform-independent wheel contains Python code, not the desktop executable or native calculation library. Install/build the matching current AtomForge application for native features. Configure its location before using builders or electronic calculations:

```python
import os
os.environ['ATOMFORGE_PATH'] = r'C:\path\to\AtomForge.exe'
# Optional explicit override, including for a headless library-only installation:
# os.environ['ATOMFORGE_ELECTRONIC_LIBRARY'] = r'C:\path\to\atomforge_electronic.dll'
```

On Linux/macOS, use the corresponding executable and `.so`/`.dylib`. The electronic loader also searches alongside the executable found on PATH and in standard repository/release layouts. Importing the package does not load native binaries. Missing libraries/executables produce an actionable error when a native feature is called. Wulff scripting requires the updated executable with `--shape wulff` support.

## Build and edit structures

```python
import atomforge as af

result = af.build('bulk', ['--system', 'cubic', '--spacegroup', '225',
                           '--a', '3.61', '--atom', 'Cu 0 0 0'])
cu = result.structure
assert len(cu) == 4
host = cu.repeat(5, 5, 5)
alloy = af.build('sss', ['--frac', 'Cu=0.7,Ni=0.3', '--seed', '42'],
                 source=host, output='cu_ni.cif')
print(alloy.stdout, alloy.stderr)  # retain native warnings and provenance
alloy.structure.view()

particle = af.wulff(cu, [(1, 0, 0, 1.0), (1, 1, 1, 1.0)],
                    radius=20, vacuum=5, output='wulff.cif')
print(len(particle.structure))
```

`build(mode, options=(), *, source=None, output=None, executable=None, timeout=300)` accepts all eight native modes: `bulk`, `gb`, `poly`, `nano`, `amorphous`, `sss`, `dislocation`, and `custom`. Use `af.builder_help('custom')` for the exact options supported by the installed executable. Multi-component values occupy one argument (`'Cu 0 0 0'`); repeated `--atom`, `--element`, `--euler`, and `--facet` flags are supported. The call runs without a command shell or a GUI window.

`source` accepts a `Structure` or a path. Without `output`, temporary files are cleaned after the result is loaded. An explicit output is written/replaced and retains native sidecars. `BuildResult` contains `structure`, `stdout`, `stderr`, and `output` (an absolute Path or None). Python structures contain atoms and cell vectors; grain-region/IPF sidecar metadata is retained on disk but not imported into the Python object. Native failures raise `subprocess.CalledProcessError` with diagnostic streams; time limits raise `subprocess.TimeoutExpired`.

```python
sphere = af.build('nano', ['--shape', 'sphere', '--radius', '10'], source=cu)
gb = af.build('gb', ['--axis', '0 0 1', '--sigma', '5', '--uca', '3',
                      '--ucb', '3', '--overlap', '1.5'], source=cu)
poly = af.build('poly', ['--sizex', '25', '--sizey', '25', '--sizez', '25',
                         '--grains', '4', '--seed', '7'], source=cu)
glass = af.build('amorphous', ['--element', 'Si 40', '--element', 'O 80',
                               '--density', '1.5', '--seed', '7'])
# Supply your closed OBJ/STL mesh:
# filled = af.build('custom', ['--mesh', 'octahedron.obj', '--scale', '12'], source=cu)
# Defects require physically appropriate host geometry and parameters:
# defect = af.build('dislocation', ['--character', 'edge', '--shape', 'cylinder',
#                                   '--cyl-radius', '5', '--core', '1.2', '--cutoff', '8'], source=host)

# Anisotropic elasticity (Stroh sextic formalism) instead of the default
# isotropic model -- elastic constants in GPa, same axes as the input cell:
cu_aniso = af.build('dislocation', ['--character', 'edge', '--shape', 'cylinder', '--cyl-radius', '15',
                                    '--anisotropic', '--elastic-c11', '168.4',
                                    '--elastic-c12', '121.4', '--elastic-c44', '75.4'], source=host)
```

Wulff `radius` is the maximum facet-plane distance in angstrom, not the farthest vertex radius. Facets are `(h,k,l,positive_relative_energy)` rows, expanded using the native symmetry engine. Geometric builders do not relax structures.

`--anisotropic` replaces the isotropic (`--nu`) dislocation displacement field with the anisotropic Stroh sextic formalism (see `AnisotropicDislocation.h` for the underlying literature and `AtomForge --help dislocation` for all elastic-constant flags, including `--elastic-symmetry hexagonal`). Use real single-crystal elastic constants for accuracy; elastically isotropic or high-symmetry orientations need the automatic `--elastic-noise` perturbation (on by default) since the sextic formalism is mathematically singular exactly there. `--dipole --dipole-offset "dx dy"` adds a second, opposite-Burgers-vector dislocation, giving a periodicity-compatible structure (zero net Burgers vector).

```python
# Nye (dislocation density) tensor: compare a dislocated structure against
# its undeformed reference (same atom count/order in both).
rows = af.nye_tensor(dislo.structure, reference, cutoff=3.0, no_pbc=True)
core_atoms = sorted(rows, key=lambda r: -r["norm"])[:10]
```
`nye_tensor` returns one row per atom (`index`, `symbol`, the nine `alpha_xx..alpha_zz` tensor components in 1/Angstrom, and their Frobenius `norm`), following the same Hartley & Mishin lattice-correspondence algorithm BABEL's own `nyeTensor.f90` implements. Atoms far from any lattice defect have a norm near zero; the tensor is largest right at a dislocation core.

A few more BABEL-equivalent post-processing tools, each a native reimplementation of the corresponding BABEL program (see the header comment in each source file for the exact literature/algorithm followed):

```python
# Differential-displacement (Vitek) map, for plotting a dislocation core's
# characteristic pattern (e.g. the BCC screw core's 3-fold arrangement).
pairs = af.vitek_map(dislo.structure, reference, line=(1, 1, 1), burgers=2.48)

# Pattern-based defect detection: flag atoms whose local neighbor
# environment no longer matches a perfect reference (near a core, fault, or
# surface).
pattern = af.build_pattern(reference, cutoff=3.2)
rows = af.detect_pattern(dislo.structure, pattern, angle_threshold=10.0)

# Prepare a constrained-minimization ("drag") migration-barrier calculation
# by interpolating between two configurations (AtomForge does not itself
# run the constrained minimization, matching BABEL's own prepareDrag).
midpoint, constraint_directions = af.prepare_drag(initial, final, zeta=0.5)

# Recover a dislocation's position/Burgers vector from a measured field
# (the reverse of af.build inserting a known one).
fit = af.fit_dislocation(dislo.structure, reference, line=(0, 0, 1))
```

`af.load` / `af.save` support XYZ/extXYZ, VASP POSCAR/CONTCAR, PDB, explicit-site CIF, and LAMMPS data. Python CIF reading does not expand asymmetric-unit symmetry; use the native bulk builder for that. `Structure.copy`, `repeat`, and `filter_species` return separate objects; `translate` and `scale` mutate in place. `view()` opens the desktop; `view_notebook()` displays the existing interactive notebook viewer.

## Render snapshots without the GUI

`af.render` produces a PNG snapshot of a structure headlessly (an offscreen OpenGL context is still required — a GPU/display driver, not a visible window):

```python
cu.set_element_color('Cu', 0.9, 0.5, 0.2)
af.render(cu, 'cu.png', width=1600, height=1200, yaw=30, pitch=20,
          radii={'Cu': 1.4}, background=(1, 1, 1))

# Turntable: writes cu_turn-000.png .. cu_turn-011.png
af.render(cu, 'cu_turn.png', frames=12, yaw_step=30)
```

`render(structure, output, *, width=1600, height=1200, yaw=0.0, pitch=0.0, roll=0.0, distance=None, orthographic=False, background=(1,1,1), show_bonds=True, show_box=True, colors=None, radii=None, radius_scale=1.0, dpi=None, frames=None, yaw_step=None, executable=None, timeout=120)` returns the output path. `colors`/`radii` override one element's appearance at a time (`{'Fe': (0.8, 0.4, 0.1)}`, `{'Fe': 1.4}` in angstrom); any color already set via `Structure.set_element_color` or direct `atom.r/g/b` edits is used unless overridden, since the structure file written for the native renderer doesn't itself carry per-atom color. `distance` defaults to an auto-fit view. `dpi` embeds a physical resolution (dots per inch) in the saved PNG's metadata for print/publication use — it does not change the pixel dimensions (`width`/`height` do that); omitted, no resolution metadata is written. Color/size customization is per-element, matching the desktop GUI's own Edit Structure dialog. The same options are available from the CLI: `AtomForge --render --input FILE --output FILE.png [options]` (see `AtomForge --render --help`).

## Charge transfer and electronic analysis

The existing native-backed electronic APIs are included in this release. Read VASP CHGCAR/LOCPOT-style grids, Gaussian Cube, and XSF:

```python
from atomforge.electronic import load_volume

rho = load_volume('CHGCAR-AB').fields[0]
a = load_volume('CHGCAR-A').fields[0]
b = load_volume('CHGCAR-B').fields[0]
delta = rho.density_difference(a, b)
print(delta.charge_summary())
delta.save('difference.xsf')
accumulation, depletion = delta.split_density()
mask = delta.threshold_mask(0.01, max(delta.values))
print(delta.apply_mask(mask).integrate())
```

Use fragment calculations on the same cell, geometry, grid, and units. `resample(target)` explicitly interpolates compatible fields; it does not register displaced atoms. Boolean mask operations include union, intersection, difference, xor, and complement. Charge summaries integrate electron redistribution and do not assign per-atom transfer.

Available tools also include arithmetic, smoothing, gradients/Laplacians, approximate energy densities, profiles, planar/macroscopic averages, cumulative charge, arbitrary-plane sections, contours, isosurfaces with scalar colouring, sphere/Voronoi integrals, maxima, Fourier/structure-factor/Patterson tools, model densities, scattering factors, and Ewald calculations. Surface export supports OBJ and PLY; density export supports VASP, Cube, and XSF.

```python
# Cartesian origin and span vectors in angstrom; use your own cell geometry.
section = rho.section(rho.origin, rho.cell[0], rho.cell[1], shape=(100, 100))
values_2d = section.section_values
surface = rho.isosurface(level=0.2)
surface.save('density.ply')
profile = delta.cumulative_charge(axis=2)
```

VASP density values are converted to e/angstrom³. Cube/XSF `quantity='auto'` retains raw scalar units; specify `quantity='density'` only when appropriate. Cube coordinates default to bohr. Check `shape`, `cell`, `origin`, `periodic`, and `unit` before combining fields. Desktop camera, transparency, lighting, file-dialog, and interface-scale controls are GUI presentation settings, not numerical calculation parameters.

## Documentation and release checks

The current source also adds structural analyses (`cna`, `rdf`, `adf`, `sro`,
`interstitial_sites`, `sculpt`), interface/stacking-fault builders, grain metadata,
and electronic recipes/batch processing. These require the matching native build.

The optional `atomforge.science` package adds bands/DOS/PDOS, projected-band plots,
WAVECAR orbital reconstruction, reciprocal-space BXSF surfaces, diffraction,
trajectories/GIFs, Hirshfeld integration, and ASE relaxation/dynamics/phonons.
Install the `science` extra using Python 3.12 for the tested dependency set.
Bader and DDEC use external Henkelman and Chargemol solvers; DFT requires a
caller-configured ASE calculator. The updated manual explains conventions and
examples. Existing PyPI wheels do not acquire these source additions automatically.

The detailed [AtomForge PDF manual](https://github.com/albert-hzbn/AtomForge/blob/main/docs/manual/AtomForge-manual.pdf) includes complete signatures, physical interpretation, tutorials, and native CLI options. Help > Manual and About > Manual in the updated desktop open the packaged PDF.

For maintainers, run CTest against the native build, then build and validate both wheel and source distribution:

```sh
python -m build python
python -m twine check python/dist/atomforge_py-0.2.0*
```

Version 0.2.0 adds native builder scripting, Wulff scripting, shared executable discovery for electronic calculations, and updated package documentation. Existing structure and electronic API names remain available.

### Condensed-matter tools in the current source

The desktop's **Analysis** menu sections, Python API and scientific
CLI share twenty tools: MSD, diffusion fitting, velocity autocorrelation,
vibrational spectra, local strain/D2min, centrosymmetry, Steinhardt order,
Wigner-Seitz defects, static structure factors, band gaps, effective masses,
work functions, equation-of-state fitting, elastic tensors, phonon DOS, harmonic
thermodynamics, NEB, NVT, NPT and symmetry-based reciprocal paths.

Install the matching source with its `science` extra into Python 3.12. The
desktop lets you browse for that interpreter. These tools require NumPy >=2,
SciPy >=1.15, ASE >=3.26 and SeeK-path >=2.1 in addition to the existing science
dependencies. Simulation tools require an explicitly chosen ASE calculator.

```sh
python -m atomforge.science --catalog
python -m atomforge.science band-gap --input request.json --output gap.json
```

Here `request.json` can contain
`{"energies_eV":{"file":"EIGENVAL"},"fermi_eV":5.4}`. Relative data paths are
resolved beside the request. Numeric CSV, NPY and JSON arrays are also accepted.
JSON references can select nested fields, such as `result.msd_A2`, and numeric
tables can select a zero-based `column`. Result files preserve the full arrays
and parameters; invalid local environments have a validity mask and JSON nulls.
Existing result files require `--overwrite`.

Use `--structures frames.extxyz` to export simulation images or a returned
primitive cell. The desktop provides the same export and can open the final
structure in a new tab. The manual lists the priority order, input shapes,
units, examples and physical limits for every tool. Numerical regression tests
are in `tests/test_condensed_matter.py`; `tests/test_science_workflows.py` runs
all twenty through the entry point shared with the desktop.
