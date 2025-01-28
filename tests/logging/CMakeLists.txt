# Google Logging Test
add_executable(glog_test glog_test.cpp)
target_link_libraries(glog_test
    PRIVATE
        glog::glog
)

# Add Google Test dependencies and configure test discovery
AddTests(glog_test)

# Add format target
add_format_target(glog_test ${CMAKE_CURRENT_SOURCE_DIR})
