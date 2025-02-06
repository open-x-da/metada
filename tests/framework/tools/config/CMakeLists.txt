# Add config test executable
add_executable(config_test
    ConfigTest.cc)

# Link config library and setup testing
target_link_libraries(config_test
    PRIVATE
        metada::framework::tools::config)

# Use the AddTests macro to handle GTest/GMock setup and coverage
metada_add_test_with_coverage(config_test)


# Add format target for test files
AddFormatTarget(config_test ${CMAKE_CURRENT_SOURCE_DIR}) 