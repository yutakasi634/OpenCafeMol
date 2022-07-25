cmake_minimum_required(VERSION 3.10)
project(open_aicg2plus)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          17)

# for headder include with the path from project root
include_directories(${PROJECT_SOURCE_DIR})

# =============================================================================
# optimization flags

include(CheckCXXCompilerFlag)

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

# =============================================================================
# debugger flags

set(OPEN_AICG2PLUS_COMPILER_FLAGS "")
option(GDB "add flag for gdb")
if(GDB)
    check_cxx_compiler_flag("-ggdb" OPEN_AICG2PLUS_COMPILER_HAS_GGDB)
    if(OPEN_AICG2PLUS_COMPILER_HAS_GGDB)
        set(OPEN_AICG2PLUS_DEBUGGER_FLAGS "${OPEN_AICG2PLUS_COMPILER_FLAGS} -ggdb")
    endif()
endif()

# =============================================================================
# Build the code

add_subdirectory(src)
