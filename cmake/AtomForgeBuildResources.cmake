# Serialize heavyweight compile/link work without reducing Release optimization.
# Ninja pools avoid launching waiting processes; the lock also protects Makefiles
# and concurrent build invocations using the same build directory.
set(_atomforge_build_lock "${CMAKE_BINARY_DIR}/atomforge-compiler.lock")
set(_atomforge_build_launcher
    "${CMAKE_COMMAND}"
    "-DATOMFORGE_BUILD_LOCK=${_atomforge_build_lock}"
    -P "${CMAKE_CURRENT_LIST_DIR}/RunBuildCommand.cmake" --)

if(CMAKE_GENERATOR MATCHES "Ninja|Makefiles")
    set(CMAKE_CXX_COMPILER_LAUNCHER
        ${_atomforge_build_launcher} ${CMAKE_CXX_COMPILER_LAUNCHER})
    # RULE_LAUNCH_LINK also covers static archivers and works with CMake 3.16.
    set_property(DIRECTORY PROPERTY RULE_LAUNCH_LINK
        "\"${CMAKE_COMMAND}\" \"-DATOMFORGE_BUILD_LOCK=${_atomforge_build_lock}\" -P \"${CMAKE_CURRENT_LIST_DIR}/RunBuildCommand.cmake\" --")
    message(STATUS "AtomForge: compile and link commands are serialized to bound build memory")
endif()

if(CMAKE_GENERATOR MATCHES "Ninja")
    set_property(GLOBAL APPEND PROPERTY JOB_POOLS atomforge_build=1)
    set(CMAKE_JOB_POOL_COMPILE atomforge_build)
    set(CMAKE_JOB_POOL_LINK atomforge_build)
endif()
