# Code audit — 2026-09-09

The audit covered the first-party C++ application, command-line paths, Python
package, build configuration, and existing tests. Review focused on memory
lifetime, asynchronous results, structure editing, numerical inputs, and file
round trips. Vendored libraries were compiled but were not independently audited.
Existing untracked build, package, and temporary files were preserved.

## Changes

- Close Structure no longer uses a destroyed tab's state or leaks its scene
  buffers. Closing tabs maintains the correct active index and restores the
  replacement tab's camera.
- CNA and RDF only read worker result pointers after observing completion through
  the atomic computing flag; completion is published before the worker becomes idle.
- Picking compares ray/sphere surface intersections, handles non-unit directions
  and rays inside atoms, and rejects zero directions.
- Atom deletion preserves the corresponding grain colors and region identifiers
  and invalidates the cached dislocation geometry.
- Atom fields and absent unit-cell vectors have defined initial values.
- Solid-solution generation rejects non-finite fractions and uses double precision
  for count allocation. CLI fractions reject trailing garbage.
- Amorphous generation rejects invalid elements, overflowing count sums, and
  non-finite box dimensions; colors have a defined fallback.
- Cell matrix inversion rejects non-finite determinants. Python fractional
  conversion rejects singular cells instead of returning Cartesian values.
- Python XYZ input rejects truncated records and negative counts. POSCAR input
  rejects mismatched species/count lists; saving POSCAR/CONTCAR without a filename
  extension now selects VASP output.
- CIF export uses Cartesian coordinate tags for structures without a cell.
  Import handles blank lines and wrapped atom-site rows, rejects incomplete rows,
  and validates lattice parameters consistently with Structure.set_cell().
- LAMMPS input handles declared atomic/full/charge/molecular styles, optional image
  flags and inline comments without shifting coordinate columns. It normalizes
  the box origin to the Python structure's implicit origin zero.
- LAMMPS output preserves triclinic geometry and rotates coordinates with the
  cell. Cluster boxes enclose their atoms. Atomic masses cover all 118 elements
  using the existing C++ table; unknown symbols raise an error.
- Viewer temporary files are removed when serialization fails.
- README installation/format guidance reflects the dependency-free Python API.
  CTest now includes C++ regressions and, when Python is found, the Python suites
  and headless CLI integration tests.

## Validation

- Full Windows MinGW Release build succeeded with OpenMP and spglib enabled;
  the final build log contained no compiler warnings or errors.
- All five CTest suites passed: core regressions, existing Python API tests,
  notebook tests, twelve new Python regression tests, and CLI smoke tests.
- CLI coverage includes version output, eight help modes, generation through
  bulk/GB/poly/nano/SSS builders, and five invalid-input cases.
- Running the new Python regression suite against HEAD's original Python files
  failed; the modified package passes all twelve tests.
- `git diff --check` passed.

Reproduce with a configured build:

```text
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

If CMake does not discover Python, configure with
`-DPython3_EXECUTABLE=<absolute path to python>` before running CTest.

## Limits

This is a broad audit and targeted repair pass, not proof that the entire
codebase is bug-free. Interactive GUI workflows were reviewed and compiled but
not exercised manually. ThreadSanitizer, Linux builds, packaging installations,
large-structure stress tests, and exhaustive scientific validation of all
builders and analysis algorithms were not run. Notebook tests skipped IPython
integration because IPython is unavailable in the test interpreter.

LAMMPS import retains wrapped coordinates and ignores image counters; bonds,
charges, and molecule IDs are outside the current Python Structure model.
Left-handed or degenerate cells are rejected by the LAMMPS writer.
