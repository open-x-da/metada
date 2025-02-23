# Create logging test executable
add_executable(logger_test
    LoggerTest.cpp)

# Link logger library and setup testing
target_link_libraries(logger_test
    PRIVATE
        metada::framework::common)

# Use the AddTests macro to handle GTest/GMock setup
metada_add_test_with_coverage(logger_test)

# Add format target for test files
AddFormatTarget(logger_test ${CMAKE_CURRENT_SOURCE_DIR}) 