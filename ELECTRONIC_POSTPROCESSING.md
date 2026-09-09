# Electronic post-processing

AtomForge provides a native scalar-grid toolkit, a desktop panel under
**Analysis → Electronic Post-processing**, and Python bindings in
`atomforge.electronic`. Desktop and Python operations share the C++ engine.
Existing Python structure APIs remain independent of native libraries.

## Coverage against VESTA

The reference is the [VESTA utilities manual](https://jp-minerals.org/vesta/en/doc/VESTAch14.html)
and [2D display documentation](https://jp-minerals.org/vesta/en/doc/VESTAch15.html).
This is an original implementation of the electronic/scalar-grid workflows,
not a complete VESTA replacement or a bundled copy of its external programs.

| Workflow | AtomForge support | Interface |
| --- | --- | --- |
| VASP, Cube and XSF scalar datasets | Multiple channels, explicit quantity and coordinate conventions | Desktop, Python |
| Field arithmetic / difference densities | Add, subtract, multiply, divide, scale; alignment and units checked | Desktop, Python |
| Interpolation / smoothing | Trilinear resampling; Cartesian-distance Gaussian kernel | Desktop, Python |
| Isosurfaces | Marching tetrahedra, scalar coloring from another field, OBJ/PLY export | Desktop, Python |
| Multiple surface levels / transparency | Append surfaces and adjust preview opacity | Desktop; repeated Python calls |
| 2D sections and contours | Arbitrary Cartesian planes, heatmap preview, triangulated contours | Desktop, Python |
| Miller-plane sections | Deterministic plane frame and sampling | Python helper |
| Line profiles | Cartesian endpoints, distance/value table | Desktop, Python |
| Planar / macroscopic averages | Triclinic plane-normal distances; centered box average | Desktop, Python |
| Density integration | Total, sampled sphere, nearest-site Voronoi partitions | Desktop, Python |
| Peak search | 26-neighbor sample maxima; positions and values | Desktop, Python |
| Peak-centered integrations | Pass peak positions to Voronoi integration | Python composition |
| Density-derived fields | Cartesian gradient, Laplacian, approximate kinetic/potential/total energy density | Desktop, Python |
| Structure factors from grids | Positive-phase discrete Fourier transform, cell-volume normalization | Desktop, Python |
| Fourier synthesis / difference maps | Explicit complex reflection list; grid arithmetic for differences | Desktop, Python |
| Atomic structure factors | Supplied scattering factors, occupancies and isotropic B | C++, Python |
| Model electron / nuclear densities | Fourier synthesis using supplied Gaussian form-factor coefficients or neutron lengths | Python composition of native calculations |
| Patterson functions | Periodic autocorrelation; also applicable to synthesized model densities | Desktop, Python |
| Site potentials / Madelung energy | Neutral-cell point-charge Ewald sum, explicit charges and cutoffs | Desktop, Python |

VESTA's RIETAN-FP powder-diffraction workflow, MADEL executable, ORFFE geometry
utilities, STRUCTURE TIDY, space-group reflection expansion, bundled scattering
tables, anomalous-dispersion databases, SHELX/MEM reflection-file readers,
anisotropic displacement parameters, and exact VESTA rendering are not included.
Supply full reflection lists and the scattering data appropriate to the experiment.
Band structures, DOS, Bader partitioning and XSF BANDGRID/Fermi surfaces are also
outside this real-space scalar-grid implementation.

## Coordinates, ordering and units

The internal grid stores three lattice/span vectors as columns. Python accepts
them as three vector rows, consistently with `Structure.cell`. Values are
x-fastest: `values[(z * ny + y) * nx + x]`. Coordinates are in Angstrom.
Periodic grids exclude endpoint planes and use steps `cell[i]/shape[i]`.
Finite grids include both endpoints and use `cell[i]/(shape[i]-1)`.
All three dimensions must be at least two. Invalid numbers, singular cells,
truncated data and excessive dimensions raise errors.

### VASP

POSCAR headers support direct/Cartesian coordinates, selective dynamics,
negative-volume scaling and three component scales. One, two or four grid
channels are retained. PAW augmentation records are consumed and skipped;
they are not scalar-grid samples. VASP 4 headers retain unknown atomic numbers.
Trajectories, velocities between coordinates and grids, and more than four
channels are rejected.

For CHGCAR/CHG density, each file value is divided by the cell volume. Thus
the integrated electron count is the sum of file values divided by the number
of grid samples. LOCPOT values stay in eV; ELF is dimensionless. The documented
[CHGCAR ordering and electron-count check](https://vasp.at/wiki/CHGCAR),
[VASP normalization clarification](https://vasp.at/forum/viewtopic.php?t=3865)
and [LOCPOT convention](https://vasp.at/wiki/LOCPOT) inform these choices.

`quantity="auto"` recognizes uppercase names starting with CHGCAR, LOCPOT or
ELFCAR, and the exact name CHG. Renamed files need an explicit quantity. Spin
labels are `total, magnetization` or `total, mx, my, mz`; potential-file channel
semantics depend on the VASP calculation. Do not apply electron-density energy
conversion to signed spin or deformation density.

Export writes scalar grids without PAW augmentation. These are post-processing
files, not complete VASP restart files. Density values are multiplied by the
cell volume on output. Output requires a periodic grid, zero origin and atoms.

### Gaussian Cube

The reader handles z-fastest ordering, NVAL channels, negative atom counts
with orbital IDs, and Fortran D exponents. Coordinates default to Bohr;
`cube_coordinates="angstrom"` explicitly supports nonstandard writers. A
negative axis count alone never silently changes units because conventions vary.
The [Cube community specification](https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html)
documents the common layouts and their ambiguities.

Scalar values remain raw unless a quantity is specified. `density` converts
e/Bohr³ to e/Å³; `potential` converts Hartree to eV. The coordinate-unit option
does not change the assumed scalar units. Orbital amplitudes remain raw; they
are not silently squared or interpreted as electron densities.

### XSF / Quantum ESPRESSO exports

The reader supports multiple real-space DATAGRID_3D blocks, PRIMCOORD and
molecular ATOMS sections, names or atomic numbers, and x-fastest ordering.
As described in the [XCrySDen XSF specification](https://web.mit.edu/xcrysden_v1.5.60/www/XCRYSDEN/doc/XSF.html),
DATAGRID is a general grid: span endpoints are included. Animated XSF and
reciprocal-space BANDGRID data are rejected. This imports exported grids,
not Quantum ESPRESSO binary save directories.

XSF does not standardize the physical unit of every scalar quantity. Auto
leaves values raw. Explicit `density` means the producer supplied e/Å³ and
`potential` means eV. If a producer writes atomic/Rydberg units, load raw,
apply the documented conversion, and construct a grid with the resulting unit.

For Cube/XSF spanning a complete periodic cell, call `as_periodic()` or use
the desktop endpoint tool. It verifies opposite planes before removing them.
Export from periodic grids to Cube/XSF adds the endpoint planes. Format
interchange does not preserve arbitrary unit labels or Cube orbital IDs;
pass the correct quantity when reloading. Cube density/potential export uses
atomic units; XSF export retains internal scalar units.

## Numerical methods and limitations

- Finite-grid integrals use tensor-product trapezoidal weights. Periodic
  integrals use uniform voxel volumes. Sphere and Voronoi boundaries are
  sampled, not analytically clipped: refine the grid to establish convergence.
  Equidistant Voronoi samples split weight among sites, conserving the total.
  Minimum-image distances enumerate all potentially closer translations in
  triclinic cells; highly ill-conditioned cells are rejected.
- Gradients use central differences internally and second-order one-sided
  differences at finite boundaries. Laplacians include mixed derivatives and
  the full lattice metric. Finite derivatives need at least three samples per
  axis; finite Laplacians need four. Smoothing uses a finite Cartesian Gaussian
  kernel. Finite boundaries renormalize the available kernel and need not
  conserve the global integral.
- Energy conversion follows the gradient expansion and local-virial relations
  in the VESTA reference (lambda=1/72, k=1/6), evaluated in atomic units and
  returned in eV/Å³. It is an approximation, not orbital-derived kinetic energy.
  Negative density is rejected; density at or below the specified positive
  floor is masked to zero. Inspect sensitivity to floor and grid spacing.
- FFTs use radix two for power-of-two axes and direct axis transforms otherwise.
  A workload limit prevents very expensive non-power-of-two transforms.
  Reflections are restricted to the unique Nyquist range. Synthesis requires
  explicit conjugate symmetry; it rejects material imaginary residuals.
  No space-group expansion is inferred. Atomic model synthesis omits even-grid
  Nyquist planes and has finite-resolution termination effects.
- Ewald uses real and reciprocal spherical cutoffs, self-potential subtraction,
  and a conducting boundary with no dipole surface term. This boundary is also
  described by [LAMMPS](https://docs.lammps.org/kspace_style.html). The cell must
  be neutral; there is no implicit background. Potentials are volts and the
  per-cell energy is eV. Results require cutoff/alpha convergence checks and
  are point-charge model values, not a reconstruction of DFT electrostatics.
- Contours use a fixed triangulation to resolve saddle ambiguity. Isosurfaces
  use six tetrahedra per voxel and return unwelded triangle vertices. PLY stores
  the scalar as a vertex property; OBJ stores it as comments. The desktop uses
  the full mesh in an OpenGL viewport with a light background, smooth normals,
  diffuse/specular lighting and depth testing. Transparent triangles are sorted
  back to front; intersecting transparent surfaces can still have sorting artifacts.
  Full meshes are exported.
- Sections return a two-layer grid with identical values; the first layer is
  the actual 2D section. Use `section_values` and `contours`, not the artificial
  slab's volume integral. Peak search operates on samples with a deterministic
  adjacent tie break; it does not refine critical points or identify Bader basins.

## Desktop controls

The electronic workspace uses the application's standard button padding and
label/value table spacing beside the viewport. Its title bar cannot collapse.
**Open volume**,
**Reference** and **Save as** reuse the application's directory navigation and
file list, including location/drive shortcuts and overwrite confirmation.
Click folders and select a file to load it; no path entry is required. The optional
address bar supports direct folder navigation. **Save as** starts in the source
folder, supplies a filename, and remembers the chosen destination across format
changes. Source and reference filenames are displayed separately; **Reload volume**
always reloads the source. Import options and lighting controls expand when needed.

Drag in the viewport to orbit, right-drag to pan, and use the wheel to zoom.
**Fit view** or a double-click resets the camera. The viewport keeps the computed
surface visible while other analyses run. Before a surface is calculated, it
shows a density slice.

The **Spectrum** palette maps low values through blue, cyan, green and yellow to
red. Diverging blue-white-red and sequential blue palettes are also available.
The legend shows the scalar range and units; disable **Automatic range** to set
consistent bounds across calculations. A single isosurface has a constant density
and therefore one scalar color. Append multiple levels or enable **Color from
reference** to compare levels or map another field, such as electrostatic potential,
onto the surface. Opacity, specular strength and shininess are adjustable.

## Build and Python setup

The shared library is built with the desktop app, or independently:

```sh
cmake -S . -B build -DATOMFORGE_BUILD_APP=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python -m pip install ./python
```

GLM and a C++17 compiler are required for the core. Set
`ATOMFORGE_ELECTRONIC_LIBRARY` to the absolute path of
`atomforge_electronic.dll`, `.so`, or `.dylib` when it is not next to the Python
package, AtomForge executable, or a supported `build`/`lib` directory.
`ATOMFORGE_PATH` also helps locate a sibling library. The library and Python
must have the same architecture; runtime dependencies must be installed or
bundled. Windows portable packages include these dependencies. A pip install
of the source Python package does not compile or download the native library.
No NumPy or SciPy dependency is required.

```python
from atomforge.electronic import load_volume

rho = load_volume("CHGCAR").fields[0]
reference = load_volume("reference/CHGCAR").fields[0]
difference = rho - reference
difference.save("difference.xsf")
print("Electrons:", rho.integrate())
print("Atomic Voronoi populations:", rho.voronoi_integrate())

potential = load_volume("LOCPOT").fields[0]
profile = potential.planar_average(axis=2)
surface = rho.isosurface(0.1, color=potential)
surface.save("density-potential.ply")
kinetic, potential_energy, total_energy = rho.energy_density()
total_energy.save("total-energy.xsf")
```

For fragment/reference grids with different sampling, explicitly resample onto
the target first. Ensure the physical origins and structures are aligned;
interpolation does not register unrelated calculations.

```python
from atomforge.electronic import miller_section, model_density

section = miller_section(rho, (1, 1, 0), width=5, height=5)
segments = section.contours(0.1)
peak_sites = [row[:3] for row in rho.peaks()]
peak_populations = rho.voronoi_integrate(peak_sites) if peak_sites else []

# Provide the actual published coefficients for every species in your model.
# coefficients = {atomic_number: (a_coefficients, b_coefficients, constant)}
# free_atom_density = model_density(rho, coefficients=coefficients)
# nuclear = model_density(rho, scattering_lengths=scattering_lengths)
```

The command-line interface is available without the desktop application:

```sh
python -m atomforge.electronic CHGCAR --operation integrate
python -m atomforge.electronic LOCPOT --operation planar --axis 2 --output potential.csv
python -m atomforge.electronic density.cube --quantity density --operation isosurface --level 0.1 --output density.ply
```

## Validation

CTest includes native analytical regressions and the Python electronic suite.
Fixtures check VASP normalization, scaling and spin/augmentation handling;
Cube channel ordering and units; XSF endpoints and molecular atoms; all three
format round trips; finite/periodic integration; triclinic polynomial
derivatives; Fourier round trips with shifted origins and odd/even sizes;
Patterson normalization; plane contours/meshes; density-energy identities;
NaCl Madelung energy and cutoff convergence; error recovery and CLI behavior.
CI runs the suites on Linux and Windows. Platform CI results are separate from
local Windows validation.

Local Windows/MinGW validation includes all 9 desktop CTest suites, including
15 Python electronic test cases and an actual OpenGL rendering test covering
background, lighting, density colors, orbit and graphics-state restoration.
The renderer test skips explicitly when no graphics context is available.
The earlier toolkit also passed all 7 headless suites and loaded its native
library from a built/installed wheel using CPython 3.12 and the portable runtime.
Desktop checks cover loading and saving through the shared browser, the lit
isosurface and wheel zoom. Reloading the exported XSF preserves its density integral.
Linux execution remains for platform CI.
