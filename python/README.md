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
```

Wulff `radius` is the maximum facet-plane distance in angstrom, not the farthest vertex radius. Facets are `(h,k,l,positive_relative_energy)` rows, expanded using the native symmetry engine. Geometric builders do not relax structures.

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

`render(structure, output, *, width=1600, height=1200, yaw=0.0, pitch=0.0, roll=0.0, distance=None, orthographic=False, background=(1,1,1), show_bonds=True, show_box=True, colors=None, radii=None, radius_scale=1.0, frames=None, yaw_step=None, executable=None, timeout=120)` returns the output path. `colors`/`radii` override one element's appearance at a time (`{'Fe': (0.8, 0.4, 0.1)}`, `{'Fe': 1.4}` in angstrom); any color already set via `Structure.set_element_color` or direct `atom.r/g/b` edits is used unless overridden, since the structure file written for the native renderer doesn't itself carry per-atom color. `distance` defaults to an auto-fit view. Color/size customization is per-element, matching the desktop GUI's own Edit Structure dialog. The same options are available from the CLI: `AtomForge --render --input FILE --output FILE.png [options]` (see `AtomForge --render --help`).

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
