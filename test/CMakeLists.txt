add_executable(test_read_toml_ff_generator test_read_toml_ff_generator.cpp)

target_link_libraries(test_read_toml_ff_generator PRIVATE OpenMM::OpenMM fmt GTest::gtest_main)
set_target_properties(test_read_toml_ff_generator PROPERTIES
    COMPILE_FLAGS "${OPEN_AICG2PLUS_WARNING_FLAGS} ${OPEN_AICG2PLUS_OPTIMIZATION_FLAGS}")

include(GoogleTest)
gtest_discover_tests(test_read_toml_ff_generator)
