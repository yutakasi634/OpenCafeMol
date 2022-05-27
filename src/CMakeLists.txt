# search installed OpenMM lib, include and plugin path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
find_package(OpenMM REQUIRED)
# set environment variable for openmm plugin directory
get_target_property(OPENMM_PLUGIN_DIR OpenMM::OpenMM PLUGIN_DIRECTORY)

set(CMAKE_CXX_FLAGS "-DOPENMM_PLUGIN_DIR=${OPENMM_PLUGIN_DIR}")
add_executable(open_aicg2plus OpenAICG2Plus.cpp)
set_target_properties(open_aicg2plus PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

target_include_directories(open_aicg2plus PRIVATE OpenMM::OpenMM)
target_link_libraries     (open_aicg2plus PRIVATE OpenMM::OpenMM)
