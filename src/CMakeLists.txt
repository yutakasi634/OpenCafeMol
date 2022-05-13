# search installed OpenMM lib and include path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
find_package(OpenMM REQUIRED)

add_executable(open_aicg2plus OpenAICG2Plus.cpp)
set_target_properties(open_aicg2plus PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

target_include_directories(open_aicg2plus PRIVATE OpenMM::OpenMM)
target_link_libraries     (open_aicg2plus PRIVATE OpenMM::OpenMM)
