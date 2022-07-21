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
option(DEBUG "disable optimization flags" OFF)
if(DEBUG)
    check_cxx_compiler_flag("-O0" OPEN_AICG2PLUS_COMPILER_HAS_O0)
    if(OPEN_AICG2PLUS_COMPILER_HAS_O0)
        set(OPEN_AICG2PLUS_OPTIMIZATION_FLAGS "-O0")
    endif()
else()
    check_cxx_compiler_flag("-Ofast" OPEN_AICG2PLUS_COMPILER_HAS_OFAST)
    if(OPEN_AICG2PLUS_COMPILER_HAS_OFAST)
        set(OPEN_AICG2PLUS_OPTIMIZATION_FLAGS "-Ofast")
    endif()
endif()

# =============================================================================
# Build the code

add_subdirectory(src)
