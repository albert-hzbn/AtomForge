# The reusable core has no OpenGL, ImGui, or Open Babel dependencies.
find_package(Threads REQUIRED)
find_package(glm CONFIG QUIET)
if(NOT TARGET glm::glm)
    find_path(ATOMFORGE_GLM_INCLUDE_DIR glm/glm.hpp
        HINTS ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
    if(NOT ATOMFORGE_GLM_INCLUDE_DIR)
        message(FATAL_ERROR "GLM headers are required to build AtomForge core")
    endif()
endif()

add_library(atomforge_core STATIC
    ${PROJECT_SOURCE_DIR}/src/graphics/Picking.cpp
    ${PROJECT_SOURCE_DIR}/src/util/ElementData.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/AmorphousBuilder.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/SubstitutionalSolidSolutionBuilder.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/CommonNeighbourAnalysis.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/RadialDistributionAnalysis.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/AngularDistributionAnalysis.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/ShortRangeOrderAnalysis.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/InterstitialVoidAnalysis.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/VoronoiComputation.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/MeshLoader.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms/CellSculptorAlgo.cpp
)
add_library(AtomForge::Core ALIAS atomforge_core)
target_compile_features(atomforge_core PUBLIC cxx_std_17)
target_include_directories(atomforge_core PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(atomforge_core PUBLIC Threads::Threads)
if(TARGET glm::glm)
    target_link_libraries(atomforge_core PUBLIC glm::glm)
else()
    target_include_directories(atomforge_core SYSTEM PUBLIC ${ATOMFORGE_GLM_INCLUDE_DIR})
endif()
target_compile_options(atomforge_core PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wshadow>)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(atomforge_core PUBLIC OpenMP::OpenMP_CXX)
endif()
