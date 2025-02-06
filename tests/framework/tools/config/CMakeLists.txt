# Add config test executable
add_executable(config_test
    ConfigTest.cc)

# Link test dependencies for config_test
target_link_libraries(config_test
    PRIVATE
        metada::framework::tools::config
        GTest::GTest
        GTest::Main
        gmock)

# Register test
gtest_discover_tests(config_test)

# Add format target for test files
AddFormatTarget(config_test ${CMAKE_CURRENT_SOURCE_DIR}) 