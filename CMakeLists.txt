cmake_minimum_required(VERSION 3.16)

project(AtomForge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ATOMFORGE_ENABLE_SPGLIB "Enable spglib/symspg support when available" ON)
option(BUILD_PORTABLE "Build portable package with bundled dependencies" OFF)

option(ATOMFORGE_ENABLE_SSS_BUILDER "Enable Substitutional Solid Solution builder UI" ON)
option(ATOMFORGE_ENABLE_SFE_BUILDER "Enable Stacking Fault (SFE) builder UI" OFF)

if(UNIX AND NOT APPLE)
    if(BUILD_PORTABLE)
        # Portable packaging copies shared libraries next to the executable,
        # so enforce shared-link mode here even if a previous cache had ON.
        set(ATOMFORGE_LINUX_STATIC_LINK OFF CACHE BOOL "Statically link Linux build dependencies" FORCE)
    else()
        option(ATOMFORGE_LINUX_STATIC_LINK "Statically link Linux build dependencies" OFF)
    endif()
else()
    option(ATOMFORGE_LINUX_STATIC_LINK "Statically link Linux build dependencies" OFF)
endif()

if(UNIX AND NOT APPLE AND ATOMFORGE_LINUX_STATIC_LINK)
    set(PKG_CONFIG_ARGN --static)
    # Prefer static archives when available for find_package() results.
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)

    find_library(ATOMFORGE_STATIC_GLFW_LIB NAMES glfw libglfw)
    find_library(ATOMFORGE_STATIC_OPENBABEL_LIB NAMES openbabel libopenbabel)
    if(NOT ATOMFORGE_STATIC_GLFW_LIB OR NOT ATOMFORGE_STATIC_OPENBABEL_LIB)
        message(FATAL_ERROR
            "ATOMFORGE_LINUX_STATIC_LINK=ON requires static archives for GLFW and Open Babel. "
            "Install static .a variants (for example libglfw.a and libopenbabel.a), or configure with "
            "-DATOMFORGE_LINUX_STATIC_LINK=OFF to build with shared system libraries.")
    endif()
endif()

find_package(PkgConfig REQUIRED)

pkg_check_modules(GLFW3 REQUIRED IMPORTED_TARGET glfw3)
pkg_check_modules(GLEW REQUIRED IMPORTED_TARGET glew)
pkg_check_modules(OPENBABEL REQUIRED IMPORTED_TARGET openbabel-3)
pkg_get_variable(OPENBABEL_PLUGINDIR openbabel-3 plugindir)
if(NOT OPENBABEL_PLUGINDIR)
    pkg_get_variable(OPENBABEL_PLUGINDIR openbabel-3 pkglibdir)
endif()
if(NOT OPENBABEL_PLUGINDIR)
    pkg_get_variable(OPENBABEL_LIBDIR openbabel-3 libdir)
    if(OPENBABEL_LIBDIR)
        set(OPENBABEL_PLUGINDIR "${OPENBABEL_LIBDIR}/openbabel")
    endif()
endif()

set(ATOMFORGE_SOURCES
    src/main.cpp
    src/app/EditorApplication.cpp
    src/app/FileDropHandler.cpp
    src/app/ImageExport.cpp
    src/app/EditorOps.cpp
    src/app/InteractionHandlers.cpp
    src/app/SceneView.cpp
    src/app/StructureFileService.cpp
    src/app/WindowSetup.cpp
    src/util/ElementData.cpp
    src/camera/Camera.cpp
    src/graphics/BillboardMesh.cpp
    src/graphics/CylinderMesh.cpp
    src/graphics/CustomStructureShaders.cpp
    src/graphics/CellSculptorShaders.cpp
    src/graphics/LowPolyMesh.cpp
    src/graphics/Renderer.cpp
    src/graphics/Picking.cpp
    src/graphics/SceneBuffers.cpp
    src/graphics/Shader.cpp
    src/graphics/ShadowMap.cpp
    src/graphics/SphereMesh.cpp
    src/graphics/StructureInstanceBuilder.cpp
    src/ui/AtomContextMenu.cpp
    src/ui/BulkCrystalBuilderDialog.cpp
    src/ui/CSLGrainBoundaryDialog.cpp
    src/ui/NanoCrystalBuilderDialog.cpp
    src/ui/CustomStructureDialog.cpp
    src/ui/MergeStructuresDialog.cpp
    src/ui/PolyCrystalBuilderDialog.cpp
    src/ui/InterfaceBuilderDialog.cpp
    src/ui/CommonNeighbourAnalysis.cpp
    src/ui/EditMenuDialogs.cpp
    src/ui/FileBrowser.cpp
    src/ui/ImGuiSetup.cpp
    src/ui/LatticePlaneOverlay.cpp
    src/ui/MillerDirectionOverlay.cpp
    src/ui/PolyhedralOverlay.cpp
    src/ui/VoronoiOverlay.cpp
    src/ui/MeasurementOverlay.cpp
    src/ui/PeriodicTableDialog.cpp
    src/ui/RadialDistributionAnalysis.cpp
    src/ui/AngularDistributionAnalysis.cpp
    src/ui/ShortRangeOrderDialog.cpp
    src/ui/StackingFaultBuilderDialog.cpp
    src/ui/StructureInfoDialog.cpp
    src/ui/TransformAtomsDialog.cpp
    src/ui/SubstitutionalSolidSolutionDialog.cpp
    src/ui/AmorphousBuilderDialog.cpp
    src/ui/InterstitialAtomsDialog.cpp
    src/ui/CellSculptorDialog.cpp
    src/util/PathUtils.cpp
    src/algorithms/CSLComputation.cpp
    src/algorithms/MeshLoader.cpp
    src/algorithms/NanoCrystalBuilder.cpp
    src/algorithms/PolyCrystalBuilder.cpp
    src/algorithms/BulkCrystalBuilder.cpp
    src/algorithms/InterfaceBuilder.cpp
    src/algorithms/StackingFaultBuilder.cpp
    src/algorithms/VoronoiComputation.cpp
    src/algorithms/SubstitutionalSolidSolutionBuilder.cpp
    src/algorithms/AmorphousBuilder.cpp
    src/algorithms/CellSculptorAlgo.cpp
    src/algorithms/ShortRangeOrderAnalysis.cpp
    src/algorithms/AngularDistributionAnalysis.cpp
    src/algorithms/InterstitialVoidAnalysis.cpp
    src/cli/CLIMode.cpp
    src/io/StructureLoader.cpp
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/imgui_impl_glfw.cpp
    imgui/imgui_impl_opengl3.cpp
)

add_executable(AtomForge ${ATOMFORGE_SOURCES})

target_compile_definitions(AtomForge PRIVATE
    ATOMFORGE_ENABLE_SSS_BUILDER=$<BOOL:${ATOMFORGE_ENABLE_SSS_BUILDER}>
    ATOMFORGE_ENABLE_SFE_BUILDER=$<BOOL:${ATOMFORGE_ENABLE_SFE_BUILDER}>
)

# On Windows, embed the application manifest via a resource script so that
# the executable declares itself as non-elevated and properly identified.
if(WIN32)
    target_sources(AtomForge PRIVATE src/app/AtomForge.rc)
endif()

target_include_directories(AtomForge PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/src/camera
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util
    ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithms
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cli
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui
)

# Vendored ImGui sources are third-party code; suppress all warnings for them
# so that our strict project warnings are not buried in third-party noise.
set(ATOMFORGE_IMGUI_SOURCES
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/imgui_impl_glfw.cpp
    imgui/imgui_impl_opengl3.cpp
)
set_source_files_properties(${ATOMFORGE_IMGUI_SOURCES} PROPERTIES COMPILE_FLAGS "-w")

target_compile_options(AtomForge PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wformat>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wshadow>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wnon-virtual-dtor>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wpedantic>
)

target_link_libraries(AtomForge PRIVATE
    PkgConfig::GLFW3
    PkgConfig::GLEW
    PkgConfig::OPENBABEL
)

if(OPENBABEL_PLUGINDIR)
    file(TO_CMAKE_PATH "${OPENBABEL_PLUGINDIR}" OPENBABEL_PLUGINDIR_NORMALIZED)
    target_compile_definitions(AtomForge PRIVATE
        ATOMFORGE_OPENBABEL_PLUGIN_DIR="${OPENBABEL_PLUGINDIR_NORMALIZED}"
    )
    message(STATUS "Open Babel plugin hint: ${OPENBABEL_PLUGINDIR_NORMALIZED}")
else()
    message(STATUS "Open Babel plugin hint not found from pkg-config; runtime auto-discovery will be used")
endif()

if(WIN32)
    target_link_libraries(AtomForge PRIVATE
        opengl32
        gdi32
        user32
        kernel32
        shell32
        comdlg32
        ole32
        uuid
    )
elseif(UNIX)
    find_package(OpenGL REQUIRED)
    find_package(Threads REQUIRED)
    target_link_libraries(AtomForge PRIVATE OpenGL::GL Threads::Threads dl)
    if(ATOMFORGE_LINUX_STATIC_LINK)
        target_link_options(AtomForge PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
endif()

if(ATOMFORGE_ENABLE_SPGLIB)
    pkg_check_modules(SPGLIB QUIET IMPORTED_TARGET spglib)
    if(NOT SPGLIB_FOUND)
        pkg_check_modules(SPGLIB QUIET IMPORTED_TARGET symspg)
    endif()

    if(SPGLIB_FOUND)
        target_compile_definitions(AtomForge PRIVATE ATOMS_ENABLE_SPGLIB)
        target_link_libraries(AtomForge PRIVATE PkgConfig::SPGLIB)
        message(STATUS "spglib/symspg found: enabling ATOMS_ENABLE_SPGLIB")
    else()
        message(STATUS "spglib/symspg not found: building without symmetry extensions")
    endif()
endif()

# OpenMP (optional) – enables multi-threaded atom generation in structure builders
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(AtomForge PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found: parallel atom generation enabled")
else()
    message(STATUS "OpenMP not found: structure builders will run single-threaded")
endif()

set(ATOMFORGE_PORTABLE_ROOT "AtomForge")
if(WIN32)
    set(ATOMFORGE_PORTABLE_BIN_DIR "${ATOMFORGE_PORTABLE_ROOT}")
else()
    set(ATOMFORGE_PORTABLE_BIN_DIR "${ATOMFORGE_PORTABLE_ROOT}/bin")
endif()

if(BUILD_PORTABLE)
    message(STATUS "Building portable package: runtime dependencies and Open Babel plugins will be bundled")

    if(UNIX AND NOT APPLE)
        # Linux portable build loads copied .so files from executable folder.
        set_target_properties(AtomForge PROPERTIES
            BUILD_RPATH "\$ORIGIN"
            INSTALL_RPATH "\$ORIGIN"
        )
    endif()
endif()

# ============================================================================
# Installation Rules
# ============================================================================
if(BUILD_PORTABLE)
    install(TARGETS AtomForge DESTINATION "${ATOMFORGE_PORTABLE_BIN_DIR}")
    install(FILES README.md DESTINATION "${ATOMFORGE_PORTABLE_ROOT}")

    # Bundle Open Babel plugin files (.obf), including versioned subdirectories.
    if(OPENBABEL_PLUGINDIR AND EXISTS "${OPENBABEL_PLUGINDIR}")
        install(DIRECTORY "${OPENBABEL_PLUGINDIR}/"
                DESTINATION "${ATOMFORGE_PORTABLE_BIN_DIR}/openbabel")
    else()
        message(WARNING "OPENBABEL_PLUGINDIR not found; portable package may miss Open Babel plugins")
    endif()

    set(ATOMFORGE_PORTABLE_EXE_NAME "AtomForge${CMAKE_EXECUTABLE_SUFFIX}")
    get_filename_component(ATOMFORGE_CXX_COMPILER_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)

    # Collect library directories from pkg-config so the dependency scanner
    # can locate DLLs for glew, glfw, openbabel, symspg, etc.
    set(_pkg_lib_dirs "")
    foreach(_pc_libdirs GLFW3_LIBRARY_DIRS GLEW_LIBRARY_DIRS OPENBABEL_LIBRARY_DIRS SPGLIB_LIBRARY_DIRS)
        if(${_pc_libdirs})
            foreach(_d IN LISTS ${_pc_libdirs})
                get_filename_component(_d_bin "${_d}/../bin" ABSOLUTE)
                list(APPEND _pkg_lib_dirs "${_d}" "${_d_bin}")
            endforeach()
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _pkg_lib_dirs)
    set(ATOMFORGE_EXTRA_DLL_DIRS "${_pkg_lib_dirs}")

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/InstallPortableDependencies.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/InstallPortableDependencies.cmake"
        @ONLY
    )
    install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/InstallPortableDependencies.cmake")
else()
    # Standard installation (system-wide)
    if(WIN32)
        install(TARGETS AtomForge DESTINATION bin)
    else()
        install(TARGETS AtomForge DESTINATION bin)
        install(FILES README.md DESTINATION share/doc/atomforge)
    endif()
endif()

# ============================================================================
# CPack Configuration
# ============================================================================
set(CPACK_PACKAGE_NAME "AtomForge")
set(CPACK_PACKAGE_VERSION "0.2.0")
set(CPACK_PACKAGE_DESCRIPTION "Atomic structure builder for metallurgical simulation and atomistic modeling")

# Expose the version string to the source code.
target_compile_definitions(AtomForge PRIVATE ATOMFORGE_VERSION="${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_VENDOR "AtomForge")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/albert-hzbn/AtomForge")

if(WIN32)
    # Windows: ZIP archive for portable distribution
    set(CPACK_GENERATOR "ZIP")
    set(CPACK_ARCHIVE_COMPONENT_ALL TRUE)
    set(CPACK_COMPONENTS_ALL Unspecified)
else()
    # Linux: TGZ tarball for portable distribution
    set(CPACK_GENERATOR "TGZ")
    set(CPACK_ARCHIVE_COMPONENT_ALL TRUE)
    set(CPACK_COMPONENTS_ALL Unspecified)
endif()

include(CPack)
