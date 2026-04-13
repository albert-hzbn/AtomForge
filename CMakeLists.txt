cmake_minimum_required(VERSION 3.16)

project(AtomForge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ATOMFORGE_ENABLE_SPGLIB "Enable spglib/symspg support when available" ON)

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
    src/ui/PolyCrystalBuilderDialog.cpp
    src/ui/InterfaceBuilderDialog.cpp
    src/ui/CommonNeighbourAnalysis.cpp
    src/ui/EditMenuDialogs.cpp
    src/ui/FileBrowser.cpp
    src/ui/ImGuiSetup.cpp
    src/ui/LatticePlaneOverlay.cpp
    src/ui/VoronoiOverlay.cpp
    src/ui/MeasurementOverlay.cpp
    src/ui/PeriodicTableDialog.cpp
    src/ui/RadialDistributionAnalysis.cpp
    src/ui/StructureInfoDialog.cpp
    src/ui/TransformAtomsDialog.cpp
    src/util/PathUtils.cpp
    src/algorithms/CSLComputation.cpp
    src/algorithms/MeshLoader.cpp
    src/algorithms/NanoCrystalBuilder.cpp
    src/algorithms/PolyCrystalBuilder.cpp
    src/algorithms/BulkCrystalBuilder.cpp
    src/algorithms/InterfaceBuilder.cpp
    src/algorithms/VoronoiComputation.cpp
    src/io/StructureLoader.cpp
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/imgui_impl_glfw.cpp
    imgui/imgui_impl_opengl3.cpp
)

add_executable(AtomForge ${ATOMFORGE_SOURCES})

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
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui
)

target_compile_options(AtomForge PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wformat>
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

# ============================================================================
# Portable Build Option
# ============================================================================
option(BUILD_PORTABLE "Build portable package with bundled dependencies" OFF)

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
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_DESCRIPTION "Atomic structure builder for metallurgical simulation and atomistic modeling")
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
