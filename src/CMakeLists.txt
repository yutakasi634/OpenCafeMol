# =============================================================================
# Setting for use OpenMM library

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
set(CMAKE_CXX_FLAGS "${OPEN_AICG2PLUS_WARNING_FLAGS} ${OPEN_AICG2PLUS_OPTIMIZATION_FLAGS} ${OPEN_AICG2PLUS_DEBUGGER_FLAGS} -DOPENMM_PLUGIN_DIR=${OPENMM_PLUGIN_DIR}")
message(STATUS "CXX compiler flags are ${CMAKE_CXX_FLAGS}")

add_executable(open_aicg2plus OpenAICG2Plus.cpp
    IntegratorGenerator.cpp
    Observer.cpp
    SystemGenerator.cpp
    Topology.cpp
    )
set_target_properties(open_aicg2plus PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

target_include_directories(open_aicg2plus PRIVATE OpenMM::OpenMM)
target_link_libraries     (open_aicg2plus PRIVATE OpenMM::OpenMM fmt)
