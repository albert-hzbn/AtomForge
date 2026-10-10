# Native scientific analysis and simulation tools shared by the desktop and CLI.
# They need no Python runtime; spglib (optional) enables symmetry k-paths.
add_library(atomforge_science STATIC
    ${PROJECT_SOURCE_DIR}/src/science/Json.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ScienceCore.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ScienceData.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Analysis.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Potentials.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ManyBodyPotentials.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Simulation.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ReciprocalPath.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ScienceTools.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ResultPlots.cpp
    ${PROJECT_SOURCE_DIR}/src/science/AtomProperties.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Phonons.cpp
    ${PROJECT_SOURCE_DIR}/src/science/DislocationLines.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Clusters.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Diffraction.cpp
    ${PROJECT_SOURCE_DIR}/src/science/VaspElectronic.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Batch.cpp
    ${PROJECT_SOURCE_DIR}/src/science/LammpsExport.cpp
    ${PROJECT_SOURCE_DIR}/src/science/GifWriter.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ScienceCatalog.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Symmetry.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Workflows.cpp
    ${PROJECT_SOURCE_DIR}/src/science/TrajectoryStructure.cpp
    # Non-destructive structure-editing pipeline (modifiers, expressions, text syntax).
    ${PROJECT_SOURCE_DIR}/src/pipeline/Pipeline.cpp
    ${PROJECT_SOURCE_DIR}/src/pipeline/Modifiers.cpp
    ${PROJECT_SOURCE_DIR}/src/pipeline/Expression.cpp
    ${PROJECT_SOURCE_DIR}/src/pipeline/PipelineEditor.cpp
    # Nye-tensor algorithms shared with the desktop dislocation tools.
    ${PROJECT_SOURCE_DIR}/src/algorithms/NyeTensor.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/DislocationFit.cpp)
target_link_libraries(atomforge_science PUBLIC AtomForge::Core)
target_compile_options(atomforge_science PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wshadow>)

find_package(PkgConfig QUIET)

# Optional codecs for compressed VASP EIGENVAL inputs (.gz, .bz2, .xz):
# CMake package first, then pkg-config (MSYS2 and most Linux distributions).
function(atomforge_science_codec definition package target module)
    find_package(${package} QUIET)
    if(TARGET ${target})
        target_link_libraries(atomforge_science PRIVATE ${target})
    elseif(PkgConfig_FOUND)
        pkg_check_modules(ATOMFORGE_CODEC_${definition} QUIET IMPORTED_TARGET ${module})
        if(NOT ATOMFORGE_CODEC_${definition}_FOUND)
            return()
        endif()
        target_link_libraries(atomforge_science PRIVATE PkgConfig::ATOMFORGE_CODEC_${definition})
    else()
        return()
    endif()
    target_compile_definitions(atomforge_science PRIVATE ATOMFORGE_SCIENCE_${definition})
endfunction()
atomforge_science_codec(ZLIB ZLIB ZLIB::ZLIB zlib)
atomforge_science_codec(BZIP2 BZip2 BZip2::BZip2 bzip2)
atomforge_science_codec(LZMA LibLZMA LibLZMA::LibLZMA liblzma)

if(ATOMFORGE_ENABLE_SPGLIB)
    if(PkgConfig_FOUND)
        pkg_check_modules(ATOMFORGE_SCIENCE_SPGLIB QUIET IMPORTED_TARGET spglib)
        if(NOT ATOMFORGE_SCIENCE_SPGLIB_FOUND)
            pkg_check_modules(ATOMFORGE_SCIENCE_SPGLIB QUIET IMPORTED_TARGET symspg)
        endif()
    endif()
    if(ATOMFORGE_SCIENCE_SPGLIB_FOUND)
        target_compile_definitions(atomforge_science PRIVATE ATOMS_ENABLE_SPGLIB)
        target_link_libraries(atomforge_science PUBLIC PkgConfig::ATOMFORGE_SCIENCE_SPGLIB)
    else()
        message(STATUS "spglib not found: symmetry reciprocal-space paths are unavailable")
    endif()
endif()

# Shared library behind the Python package's atomforge.native module.
set_target_properties(atomforge_science PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(atomforge_science_native SHARED ${PROJECT_SOURCE_DIR}/src/science/NativeAPI.cpp)
target_link_libraries(atomforge_science_native PRIVATE atomforge_science)
set_target_properties(atomforge_science_native PROPERTIES PREFIX "" CXX_VISIBILITY_PRESET hidden)
if(NOT ATOMFORGE_BUILD_APP)
    install(TARGETS atomforge_science_native RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
endif()
