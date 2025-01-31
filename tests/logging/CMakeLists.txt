# Google Logging Test
if(Glog_FOUND AND USE_GLOG)
    add_executable(glog_test glog_test.cpp)
    target_link_libraries(glog_test
        PRIVATE
            glog::glog
    )

    # Add Google Test dependencies and configure test discovery
    AddTests(glog_test)

    # Add format target
    AddFormatTarget(glog_test ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# Create logging test executable
add_executable(logging_test logging_test.cpp)

target_link_libraries(logging_test
    PRIVATE
        metada::core
        GTest::gtest
        GTest::gtest_main
)

# Register test
gtest_discover_tests(logging_test)

# Add format target
AddFormatTarget(logging_test ${CMAKE_CURRENT_SOURCE_DIR})
