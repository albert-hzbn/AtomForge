#!/bin/sh
# Rebuilds and runs only the algorithms dumper (see build.sh).
set -e
BUILD=${1:-build-global-safety}
APP_ALGORITHMS="BulkCrystalBuilder ElasticConstants StrainTool VacancyBuilder SQSBuilder DragPrep NanostructureTools PatternMatch StackingFaultBuilder DislocationBuilder AnisotropicDislocation VitekMap CSLComputation SurfaceBuilder InterfaceBuilder InterfaceMatching PolyCrystalBuilder NanoCrystalBuilder"
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    -DATOMS_ENABLE_SPGLIB -Isrc/util -Isrc/algorithms rust/tools/golden/golden_algorithms.cpp $(for f in $APP_ALGORITHMS; do printf "src/algorithms/%s.cpp " $f; done) src/graphics/StructureInstanceBuilder.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_algorithms
rust/target/golden/golden_algorithms rust/tests/golden
