# Google Logging Test
add_executable(glog_test glog_test.cpp)
target_link_libraries(glog_test
    PRIVATE
        glog::glog
        GTest::gtest
        GTest::gtest_main
)

# Set C++17 standard for filesystem support
set_target_properties(glog_test PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# Only link filesystem library on older Apple systems
if(APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(glog_test PRIVATE c++fs)
endif()

gtest_discover_tests(glog_test)

# Add format target
add_format_target(glog_test ${CMAKE_CURRENT_SOURCE_DIR})
