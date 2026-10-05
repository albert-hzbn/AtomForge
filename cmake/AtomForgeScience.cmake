# Native scientific analysis and simulation tools shared by the desktop and CLI.
# They need no Python runtime; spglib (optional) enables symmetry k-paths.
add_library(atomforge_science STATIC
    ${PROJECT_SOURCE_DIR}/src/science/Json.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ScienceCore.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ScienceData.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Analysis.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Potentials.cpp
    ${PROJECT_SOURCE_DIR}/src/science/Simulation.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ReciprocalPath.cpp
    ${PROJECT_SOURCE_DIR}/src/science/ScienceTools.cpp)
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
