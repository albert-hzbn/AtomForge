if(NOT BUILD_TESTING)
    return()
endif()

add_executable(atomforge_core_tests ${PROJECT_SOURCE_DIR}/tests/core_regressions.cpp)
target_link_libraries(atomforge_core_tests PRIVATE AtomForge::Core)
add_test(NAME core_regressions COMMAND atomforge_core_tests)

add_executable(atomforge_analysis_tests ${PROJECT_SOURCE_DIR}/tests/analysis_regressions.cpp)
target_link_libraries(atomforge_analysis_tests PRIVATE AtomForge::Core)
add_test(NAME analysis_regressions COMMAND atomforge_analysis_tests)

add_executable(atomforge_electronic_tests ${PROJECT_SOURCE_DIR}/tests/electronic_regressions.cpp)
target_link_libraries(atomforge_electronic_tests PRIVATE AtomForge::Core)
add_test(NAME electronic_regressions COMMAND atomforge_electronic_tests)

if(TARGET AtomForge)
    add_executable(atomforge_nanocrystal_metadata_tests
        ${PROJECT_SOURCE_DIR}/tests/nanocrystal_metadata.cpp
        ${PROJECT_SOURCE_DIR}/src/algorithms/NanoCrystalBuilder.cpp)
    target_include_directories(atomforge_nanocrystal_metadata_tests PRIVATE
        ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/util)
    target_link_libraries(atomforge_nanocrystal_metadata_tests PRIVATE AtomForge::Core)
    if(TARGET PkgConfig::SPGLIB)
        target_compile_definitions(atomforge_nanocrystal_metadata_tests PRIVATE ATOMS_ENABLE_SPGLIB)
        target_link_libraries(atomforge_nanocrystal_metadata_tests PRIVATE PkgConfig::SPGLIB)
    endif()
    add_test(NAME nanocrystal_metadata COMMAND atomforge_nanocrystal_metadata_tests)

    add_executable(atomforge_application_paths_tests
        ${PROJECT_SOURCE_DIR}/tests/application_paths.cpp
        ${PROJECT_SOURCE_DIR}/src/util/ApplicationPaths.cpp)
    target_include_directories(atomforge_application_paths_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_compile_definitions(atomforge_application_paths_tests PRIVATE
        ATOMFORGE_SOURCE_MANUAL="${ATOMFORGE_MANUAL_PDF}")
    add_test(NAME application_paths COMMAND atomforge_application_paths_tests)

    add_executable(atomforge_ui_layout_tests
        ${PROJECT_SOURCE_DIR}/tests/ui_layout.cpp
        ${PROJECT_SOURCE_DIR}/src/ui/ResponsiveLayout.cpp
        ${PROJECT_SOURCE_DIR}/src/ui/ImGuiSetup.cpp
        ${PROJECT_SOURCE_DIR}/imgui/imgui_impl_glfw.cpp
        ${PROJECT_SOURCE_DIR}/imgui/imgui_impl_opengl3.cpp
        ${PROJECT_SOURCE_DIR}/imgui/imgui.cpp
        ${PROJECT_SOURCE_DIR}/imgui/imgui_draw.cpp
        ${PROJECT_SOURCE_DIR}/imgui/imgui_tables.cpp
        ${PROJECT_SOURCE_DIR}/imgui/imgui_widgets.cpp)
    target_include_directories(atomforge_ui_layout_tests PRIVATE
        ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/imgui)
    target_link_libraries(atomforge_ui_layout_tests PRIVATE PkgConfig::GLFW3 OpenGL::GL ${CMAKE_DL_LIBS})
    add_test(NAME ui_layout COMMAND atomforge_ui_layout_tests)

    find_package(OpenGL REQUIRED)
    add_executable(atomforge_viewport_tests
        ${PROJECT_SOURCE_DIR}/tests/electronic_viewport.cpp
        ${PROJECT_SOURCE_DIR}/src/graphics/ElectronicViewport.cpp
        ${PROJECT_SOURCE_DIR}/src/graphics/Shader.cpp)
    target_link_libraries(atomforge_viewport_tests PRIVATE AtomForge::Core PkgConfig::GLFW3 PkgConfig::GLEW OpenGL::GL)
    add_test(NAME electronic_viewport COMMAND atomforge_viewport_tests)
    set_tests_properties(electronic_viewport PROPERTIES SKIP_RETURN_CODE 77)
endif()

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME python_electronic
        COMMAND ${CMAKE_COMMAND} -E env
            "ATOMFORGE_ELECTRONIC_LIBRARY=$<TARGET_FILE:atomforge_electronic>"
            ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/python/tests/test_electronic.py)
    foreach(suite test_atomforge test_notebook test_regressions)
        add_test(NAME python_${suite}
            COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/python/tests/${suite}.py)
    endforeach()
    if(TARGET AtomForge)
        add_test(NAME python_builders
            COMMAND ${CMAKE_COMMAND} -E env "ATOMFORGE_PATH=$<TARGET_FILE:AtomForge>"
                ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/python/tests/test_builders.py)
        set_tests_properties(python_builders PROPERTIES TIMEOUT 180)
        add_test(NAME cli_smoke
            COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/cli_smoke.py $<TARGET_FILE:AtomForge>)
        set_tests_properties(cli_smoke PROPERTIES TIMEOUT 180)
    endif()
endif()
