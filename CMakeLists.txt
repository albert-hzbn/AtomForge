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
    src/ui/LatticePlaneOverlay.cpp
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

# ============================================================================
# Portable Build Option
# ============================================================================
option(BUILD_PORTABLE "Build portable package with bundled dependencies" OFF)

if(BUILD_PORTABLE)
    message(STATUS "Building portable package: dependencies will be bundled")
    
    if(WIN32)
        # Windows portable: copy required DLLs to build directory
        message(STATUS "Portable Windows build: DLLs will be copied to output directory")
        
        # Get the directory where DLLs are located
        get_target_property(GLFW3_LIB_DIRS PkgConfig::GLFW3 INTERFACE_LINK_DIRECTORIES)
        get_target_property(GLEW_LIB_DIRS PkgConfig::GLEW INTERFACE_LINK_DIRECTORIES)
        get_target_property(OPENBABEL_LIB_DIRS PkgConfig::OPENBABEL INTERFACE_LINK_DIRECTORIES)
        
        # Find common DLL paths (MSYS2 UCRT64)
        find_file(GLFW3_DLL glfw3.dll PATHS
            ${GLFW3_LIB_DIRS}
            ${CMAKE_PREFIX_PATH}/bin
            C:/msys2/ucrt64/bin
        )
        find_file(GLEW_DLL glew32.dll PATHS
            ${GLEW_LIB_DIRS}
            ${CMAKE_PREFIX_PATH}/bin
            C:/msys2/ucrt64/bin
        )
        find_file(OPENBABEL_DLL openbabel3.dll PATHS
            ${OPENBABEL_LIB_DIRS}
            ${CMAKE_PREFIX_PATH}/bin
            C:/msys2/ucrt64/bin
        )
        
        # Copy DLLs to build folder
        if(GLFW3_DLL)
            add_custom_command(TARGET AtomForge POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy "${GLFW3_DLL}" "$<TARGET_FILE_DIR:AtomForge>/"
                COMMENT "Copying GLFW3 DLL to output directory"
            )
        endif()
        if(GLEW_DLL)
            add_custom_command(TARGET AtomForge POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy "${GLEW_DLL}" "$<TARGET_FILE_DIR:AtomForge>/"
                COMMENT "Copying GLEW DLL to output directory"
            )
        endif()
        if(OPENBABEL_DLL)
            add_custom_command(TARGET AtomForge POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy "${OPENBABEL_DLL}" "$<TARGET_FILE_DIR:AtomForge>/"
                COMMENT "Copying OpenBabel DLL to output directory"
            )
        endif()
        
    elseif(UNIX AND NOT APPLE)
        # Linux portable: prepare for AppImage or tarball with bundled libs
        message(STATUS "Portable Linux build: prepare for AppImage or bundled tarball")
        
        # Set RPATH to load libraries from the same directory as executable
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
    if(WIN32)
        # Windows portable: output to a versioned folder
        install(TARGETS AtomForge DESTINATION "AtomForge")
        install(FILES README.md DESTINATION "AtomForge")
        # DLLs already copied by POST_BUILD commands above
        
    else()
        # Linux portable: create AppDir-like structure for easy redistribution
        install(TARGETS AtomForge DESTINATION "AtomForge/bin")
        install(FILES README.md DESTINATION "AtomForge")
    endif()
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
