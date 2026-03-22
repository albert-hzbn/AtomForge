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
    src/ElementData.cpp
    src/camera/Camera.cpp
    src/graphics/CylinderMesh.cpp
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
    src/ui/CommonNeighbourAnalysis.cpp
    src/ui/EditMenuDialogs.cpp
    src/ui/FileBrowser.cpp
    src/ui/ImGuiSetup.cpp
    src/ui/MeasurementOverlay.cpp
    src/ui/PeriodicTableDialog.cpp
    src/ui/RadialDistributionAnalysis.cpp
    src/ui/StructureInfoDialog.cpp
    src/ui/TransformAtomsDialog.cpp
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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/camera
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io
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
