set(TEST_NAMES
    test_read_toml_ff_generator
    test_read_genesis_input
    )

include(GoogleTest)

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)

    target_link_libraries(${TEST_NAME} PRIVATE
        OpenMM::OpenMM
        toml11::toml11
        fmt
        GTest::gtest_main
        opencm
        )
    set_target_properties(${TEST_NAME} PROPERTIES
        COMPILE_FLAGS "${OPEN_AICG2PLUS_WARNING_FLAGS} ${OPEN_AICG2PLUS_OPTIMIZATION_FLAGS}")

    gtest_discover_tests(${TEST_NAME}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
endforeach(TEST_NAME)
