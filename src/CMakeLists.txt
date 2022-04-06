# add exexutable

add_executable(open_aicg2plus open_aicg2plus.cpp)
set_target_properties(open_aicg2plus PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    COMPILE_FLAGS "-lOpenMM")

target_include_directories(open_aicg2plus PUBLIC /usr/local/openmm/include)
target_link_directories(open_aicg2plus PUBLIC /usr/local/openmm/lib)
target_link_libraries(open_aicg2plus OpenMM)
