if(ENABLE_CUDA_ACCELERATION)
    # Device Query Test
    add_executable(device_query_test device_query_test.cpp)
    target_compile_features(device_query_test PRIVATE cuda_std_11)
    set_target_properties(device_query_test PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        LANGUAGE CUDA
    )
    target_link_libraries(device_query_test 
        PRIVATE 
            cuda_settings
            GTest::gtest_main
    )
    target_include_directories(device_query_test
        PRIVATE
            ${CMAKE_SOURCE_DIR}/src
    )
    gtest_discover_tests(device_query_test)

    # Add format target
    add_format_target(device_query_test ${CMAKE_CURRENT_SOURCE_DIR})

    # CUDA Test Suite
    add_executable(cuda_test cuda_test.cpp)
    target_link_libraries(cuda_test
        PRIVATE
            cuda_settings
            glog::glog
            GTest::gtest_main
    )
    gtest_discover_tests(cuda_test)

    # Add format target
    add_format_target(cuda_test ${CMAKE_CURRENT_SOURCE_DIR})
endif() 