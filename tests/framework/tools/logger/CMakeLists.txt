# Create logging test executable
add_executable(LoggerTest LoggerTest.cpp)

target_link_libraries(LoggerTest
    PRIVATE
        metada::framework::tools::logger
        GTest::gtest_main
        GTest::gmock
)

# Register test
gtest_discover_tests(LoggerTest)

# Add format target
AddFormatTarget(LoggerTest ${CMAKE_CURRENT_SOURCE_DIR}) 