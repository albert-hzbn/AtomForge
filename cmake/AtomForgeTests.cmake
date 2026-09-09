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
        add_test(NAME cli_smoke
            COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/cli_smoke.py $<TARGET_FILE:AtomForge>)
        set_tests_properties(cli_smoke PROPERTIES TIMEOUT 180)
    endif()
endif()
