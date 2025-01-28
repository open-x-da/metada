# Google Logging Test
add_executable(glog_test glog_test.cpp)
target_link_libraries(glog_test
    PRIVATE
        glog::glog
        GTest::gtest_main
)

# Set C++17 standard for filesystem support
# This is required because:
# 1. std::filesystem is used for managing file and directory operations
# 2. std::filesystem became standardized in C++17
# 3. Earlier versions only had experimental::filesystem which required separate linking
# Setting CXX_STANDARD_REQUIRED ensures the build fails if C++17 is not available
# rather than silently falling back to an older standard
set_target_properties(glog_test PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# On older Apple systems (Clang < 9.0), the filesystem library needs to be explicitly linked
# This is because std::filesystem was experimental before C++17 and required separate linking
# On newer systems, it's part of the standard library and doesn't need explicit linking
# Reference: https://en.cppreference.com/w/cpp/filesystem
if(APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(glog_test PRIVATE c++fs)
endif()

gtest_discover_tests(glog_test)

# Add format target
add_format_target(glog_test ${CMAKE_CURRENT_SOURCE_DIR})
