#!/bin/sh
# Builds the golden-output dumper against the C++ libraries and writes the
# reference JSON files into rust/tests/golden. Run from the repository root
# after building the C++ project (default build directory: build-global-safety).
set -e
# The dumpers need the MSYS2 UCRT64 runtime DLLs (libstdc++, libgomp) first on
# PATH; other copies on PATH (e.g. Git for Windows) make them crash at start.
if [ -d /c/msys64/ucrt64/bin ]; then export PATH=/c/msys64/ucrt64/bin:$PATH; fi
BUILD=${1:-build-global-safety}
OUT=rust/tests/golden
mkdir -p "$OUT" rust/target/golden
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_dump.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -o rust/target/golden/golden_dump
rust/target/golden/golden_dump "$OUT"
# The structure loader as in the desktop build (Open Babel, spglib).
g++ -O2 -std=c++17 -fopenmp -DATOMS_ENABLE_SPGLIB -Isrc -Isrc/io -I/c/msys64/ucrt64/include/openbabel3 -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_io.cpp src/io/StructureLoader.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lopenbabel -lsymspg -o rust/target/golden/golden_io
BABEL_LIBDIR=$(ls -d /c/msys64/ucrt64/lib/openbabel/*/ 2>/dev/null | head -1) \
BABEL_DATADIR=$(ls -d /c/msys64/ucrt64/share/openbabel/*/ 2>/dev/null | head -1) \
    rust/target/golden/golden_io rust/tests/fixtures/io docs/manual/examples "$OUT"
# Science core and data readers (fixtures from tools/make_science_fixtures.py).
python rust/tools/make_science_fixtures.py rust/tests/fixtures/science
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_science.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lz -lbz2 -llzma -o rust/target/golden/golden_science
rust/target/golden/golden_science rust/tests/fixtures/science "$OUT"
# Pipeline: text and JSON forms, expressions, step outputs.
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_pipeline.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_pipeline
rust/target/golden/golden_pipeline "$OUT"
# Algorithms: builders and analyses (the app-only sources compiled in).
APP_ALGORITHMS="BulkCrystalBuilder ElasticConstants StrainTool VacancyBuilder SQSBuilder DragPrep NanostructureTools PatternMatch StackingFaultBuilder DislocationBuilder AnisotropicDislocation VitekMap CSLComputation SurfaceBuilder InterfaceBuilder InterfaceMatching PolyCrystalBuilder NanoCrystalBuilder"
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    -DATOMS_ENABLE_SPGLIB -Isrc/util -Isrc/algorithms rust/tools/golden/golden_algorithms.cpp $(for f in $APP_ALGORITHMS; do printf "src/algorithms/%s.cpp " $f; done) src/graphics/StructureInstanceBuilder.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_algorithms
rust/target/golden/golden_algorithms "$OUT"
# Science tools: every request of make_tool_requests.py through runTool.
python rust/tools/golden/make_tool_requests.py
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_tools.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_tools
rust/target/golden/golden_tools rust/tests/fixtures/tools/requests.json "$OUT"
# Eigen's 3x3 SelfAdjointEigenSolver (af_core::eigen::self_adjoint_eigen3).
g++ -O3 -DNDEBUG -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 rust/tools/golden/golden_selfadjoint.cpp -o rust/target/golden/golden_selfadjoint
rust/target/golden/golden_selfadjoint > "$OUT/self_adjoint_eigen.txt"
# (golden_pipeline above also writes the build-shape OBJ surfaces rust/tests/golden/shape_*.obj.)
# Electronic structure: volumes, field analysis, Fourier, charge, Bader, topology, LOBSTER, Wannier.
python rust/tools/golden/make_electronic_fixtures.py
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_electronic.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_electronic
rust/target/golden/golden_electronic rust/tests/fixtures/electronic "$OUT"
# Atom clouds (.afcloud), built with one thread so files are byte-reproducible.
python rust/tools/golden/make_cloud_fixtures.py
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_cloud.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_cloud
rust/target/golden/golden_cloud rust/tests/fixtures/cloud "$OUT"
# Chemistry: SMILES, molecular properties, perception, CIP labels, hydrogens.
python rust/tools/golden/make_chemistry_fixtures.py
g++ -O2 -std=c++17 -fopenmp -Isrc -Isrc/util -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_chemistry.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_chemistry
rust/target/golden/golden_chemistry rust/tests/fixtures/chemistry "$OUT"
# Chemistry operations: UFF, point groups, editing, builders, Z-matrices,
# quantum-chemistry inputs, embedding (graph, stereo and energy).
g++ -O2 -std=c++17 -fopenmp -Isrc -Isrc/util -I/c/msys64/ucrt64/include/eigen3 -I/usr/include/eigen3 \
    rust/tools/golden/golden_chemistry_ops.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_chemistry_ops
rust/target/golden/golden_chemistry_ops rust/tests/fixtures/chemistry "$OUT"
# Renderer scene data and the grain-boundary colouring (app sources compiled in).
APP_RENDER="src/graphics/SceneBuffers.cpp src/graphics/StructureInstanceBuilder.cpp src/app/SceneView.cpp src/graphics/GraphicsBackend.cpp src/camera/Camera.cpp imgui/imgui.cpp imgui/imgui_draw.cpp imgui/imgui_tables.cpp imgui/imgui_widgets.cpp"
APP_FLAGS="-O2 -std=c++17 -fopenmp -Iimgui -Isrc -Isrc/util -Isrc/camera -Isrc/graphics -Isrc/app -I/c/msys64/ucrt64/include/eigen3"
APP_LIBS="$BUILD/libatomforge_science.a $BUILD/libatomforge_core.a -lglfw3 -lglew32 -lopengl32 -lsymspg -lz -lbz2 -llzma"
python rust/tools/golden/make_render_fixtures.py
g++ $APP_FLAGS rust/tools/golden/golden_render.cpp src/app/EditorOps.cpp $APP_RENDER $APP_LIBS -o rust/target/golden/golden_render
rust/target/golden/golden_render rust/tests/fixtures/render "$OUT"
python rust/tools/golden/make_colouring_fixtures.py
g++ $APP_FLAGS rust/tools/golden/golden_colouring.cpp $APP_RENDER $APP_LIBS -o rust/target/golden/golden_colouring
rust/target/golden/golden_colouring rust/tests/fixtures/colouring "$OUT"
# Measurements and Atom Info (ui/MeasurementOverlay) over the render fixtures.
g++ $APP_FLAGS rust/tools/golden/golden_measure.cpp src/ui/MeasurementOverlay.cpp src/ui/ResponsiveLayout.cpp src/app/EditorOps.cpp $APP_RENDER $APP_LIBS -o rust/target/golden/golden_measure
rust/target/golden/golden_measure rust/tests/fixtures/render "$OUT"
# Projects (model/Workspace) in formats 1-7.
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 rust/tools/golden/golden_project.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_project
rust/target/golden/golden_project "$OUT"
# Plots of the tool results above (science/ResultPlots).
g++ -O2 -std=c++17 -fopenmp -Isrc -I/c/msys64/ucrt64/include/eigen3 rust/tools/golden/golden_plots.cpp "$BUILD/libatomforge_science.a" "$BUILD/libatomforge_core.a" \
    -lsymspg -lz -lbz2 -llzma -o rust/target/golden/golden_plots
rust/target/golden/golden_plots "$OUT/tools.json" "$OUT"
# Step options of the solid solution and custom structure dialogs (header-only Options.h).
g++ -O2 -std=c++17 -Isrc rust/tools/golden/golden_step_options_composition.cpp -o rust/target/golden/golden_step_options_composition
rust/target/golden/golden_step_options_composition "$OUT"
# editStep + stepOptions of the Nanocrystal, Polycrystal and Amorphous builder
# dialogs, linked against the desktop app's objects in $BUILD (all but main).
(cd "$BUILD" && sed 's#CMakeFiles/AtomForge.dir/src/main.cpp.obj##' CMakeFiles/AtomForge.dir/objects1.rsp > ../rust/target/golden/step_builders_objects.rsp && \
    g++ -O2 -std=c++17 -fopenmp -DATOMFORGE_ENABLE_SFE_BUILDER=1 -DATOMFORGE_ENABLE_SSS_BUILDER=1 -DATOMFORGE_VULKAN -DATOMS_ENABLE_SPGLIB \
    -DIMGUI_IMPL_VULKAN_NO_PROTOTYPES -DVK_NO_PROTOTYPES @CMakeFiles/AtomForge.dir/includes_CXX.rsp ../rust/tools/golden/golden_step_builders.cpp \
    @../rust/target/golden/step_builders_objects.rsp @CMakeFiles/AtomForge.dir/linkLibs.rsp -o ../rust/target/golden/golden_step_builders)
rust/target/golden/golden_step_builders rust/tools/golden/step_builders_cases.tsv "$OUT/step_builders.tsv"
# editStep + stepOptions of the CSL Grain Boundary, Interface Builder and
# Stacking Fault dialogs (the same objects; cases in rust/tests/fixtures).
(cd "$BUILD" && sed 's#CMakeFiles/AtomForge.dir/src/main.cpp.obj##' CMakeFiles/AtomForge.dir/objects1.rsp > ../rust/target/golden/step_grain_objects.rsp &&     g++ -O2 -std=c++17 -fopenmp -DATOMFORGE_ENABLE_SFE_BUILDER=1 -DATOMFORGE_ENABLE_SSS_BUILDER=1 -DATOMFORGE_VULKAN -DATOMS_ENABLE_SPGLIB     -DIMGUI_IMPL_VULKAN_NO_PROTOTYPES -DVK_NO_PROTOTYPES @CMakeFiles/AtomForge.dir/includes_CXX.rsp ../rust/tools/golden/golden_step_grain.cpp     @../rust/target/golden/step_grain_objects.rsp @CMakeFiles/AtomForge.dir/linkLibs.rsp -o ../rust/target/golden/golden_step_grain)
rust/target/golden/golden_step_grain rust/tests/fixtures/step_grain_cases.tsv "$OUT/step_grain.tsv"
