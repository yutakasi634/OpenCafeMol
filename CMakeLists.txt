cmake_minimum_required(VERSION 3.10)
enable_testing()
project(open_cafemol)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          17)

# for headder include with the path from project root
include_directories(${PROJECT_SOURCE_DIR})


# =============================================================================
# include CMake modules
include(CheckCXXCompilerFlag)

# =============================================================================
# Setting for OpenMM library

# search installed OpenMM lib, include and plugin path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
find_package(OpenMM REQUIRED)
message(STATUS "OpenMM include path : ${OPENMM_INCLUDE_DIR}")
message(STATUS "OpenMM library path : ${OPENMM_LIBRARY}")

# set option for openmm plugin directory
get_target_property(EXPECTED_OPENMM_PLUGIN_DIR OpenMM::OpenMM PLUGIN_DIRECTORY)
set(OPENMM_PLUGIN_DIR ${EXPECTED_OPENMM_PLUGIN_DIR} CACHE PATH "The plugin path for OpenMM")
message(STATUS "OpenMM plugin  path : ${OPENMM_PLUGIN_DIR}")

# =============================================================================
# Setting for toml11
set(TOML11_PRECOMPILE ON)

# =============================================================================
# warning flags

set(OPEN_AICG2PLUS_WARNING_FLAGS "")
check_cxx_compiler_flag("-Wall" OPEN_AICG2PLUS_COMPILER_HAS_WALL)
if(OPEN_AICG2PLUS_COMPILER_HAS_WALL)
    set(OPEN_AICG2PLUS_WARNING_FLAGS "${OPEN_AICG2PLUS_WARNING_FLAGS} -Wall")
endif()

message(STATUS "warning flags are ${OPEN_AICG2PLUS_WARNING_FLAGS}")

# =============================================================================
# optimization flags

set(OPEN_AICG2PLUS_OPTIMIZATION_FLAGS "")
option(DEBUG "compile in debug mode" OFF)
if(DEBUG)
    check_cxx_compiler_flag("-O0" OPEN_AICG2PLUS_COMPILER_HAS_O0)
    if(OPEN_AICG2PLUS_COMPILER_HAS_O0)
        set(OPEN_AICG2PLUS_OPTIMIZATION_FLAGS "${OPEN_AICG2PLUS_OPTIMIZATION_FLAGS} -O0")
    endif()
    check_cxx_compiler_flag("-g" OPEN_AICG2PLUS_COMPILER_HAS_G)
    if(OPEN_AICG2PLUS_COMPILER_HAS_O0)
        set(OPEN_AICG2PLUS_OPTIMIZATION_FLAGS "${OPEN_AICG2PLUS_OPTIMIZATION_FLAGS} -g")
    endif()
else()
    check_cxx_compiler_flag("-Ofast" OPEN_AICG2PLUS_COMPILER_HAS_OFAST)
    if(OPEN_AICG2PLUS_COMPILER_HAS_OFAST)
        set(OPEN_AICG2PLUS_OPTIMIZATION_FLAGS "${OPEN_AICG2PLUS_OPTIMIZATION_FLAGS} -Ofast")
    endif()
endif()

message(STATUS "optimization flags are ${OPEN_AICG2PLUS_OPTIMIZATION_FLAGS}")

# =============================================================================
# debugger flags

set(OPEN_AICG2PLUS_DEBUGGER_FLAGS "")
option(GDB "add flag for gdb" OFF)
if(GDB)
    check_cxx_compiler_flag("-ggdb" OPEN_AICG2PLUS_COMPILER_HAS_GGDB)
    if(OPEN_AICG2PLUS_COMPILER_HAS_GGDB)
        set(OPEN_AICG2PLUS_DEBUGGER_FLAGS "${OPEN_AICG2PLUS_DEBUGGER_FLAGS} -ggdb")
    endif()
endif()

message(STATUS "debugger flags are ${OPEN_AICG2PLUS_DEBUGGER_FLAGS}")

# =============================================================================
# compiler-dependent library configuration
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 8.0 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        link_libraries(stdc++fs)
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
    link_libraries(stdc++fs)
endif()

# =============================================================================
# Build the code

add_subdirectory(toml11)
add_subdirectory(fmt)
add_subdirectory(googletest)
add_subdirectory(src)
add_subdirectory(test)
