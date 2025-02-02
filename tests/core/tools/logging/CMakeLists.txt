# Create logging test executable
add_executable(LoggerTest LoggerTest.cpp)

target_link_libraries(LoggerTest
    PRIVATE
        metada::core::tools::logging
        metada::backends::tools::logging
        GTest::gtest
        GTest::gtest_main
)

# Register test
gtest_discover_tests(LoggerTest)

# Add format target
AddFormatTarget(LoggerTest ${CMAKE_CURRENT_SOURCE_DIR})
