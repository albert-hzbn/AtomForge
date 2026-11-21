#!/bin/sh
# Runs the Python package's test suite (python/tests) against the Rust native
# libraries (af-ffi-science -> atomforge_science_native.dll, af-ffi-electronic
# -> atomforge_electronic.dll). Run from the repository root. PYTHON selects
# an interpreter with numpy/scipy/ase (default: python); ATOMFORGE_PATH, when
# set, enables the tests that need the executable.
set -e
PYTHON=${PYTHON:-python}
cd rust
cargo build --release -p af-ffi-science -p af-ffi-electronic
# spglib's DLL must sit next to the libraries (Python does not search PATH).
cp crates/af-spglib/vendor/lib/symspg.dll target/release/
LIB=$(cd target/release && pwd -W 2>/dev/null || pwd)
cd ../python
ATOMFORGE_SCIENCE_LIBRARY="$LIB/atomforge_science_native.dll" \
ATOMFORGE_ELECTRONIC_LIBRARY="$LIB/atomforge_electronic.dll" \
    "$PYTHON" -m unittest discover -s tests "$@"
