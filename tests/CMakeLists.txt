# Only add CUDA tests if CUDA is available
if(USE_CUDA)
    # Device Query Test (CUDA)
    add_executable(device_query_test device_query_test.cpp)
    target_compile_features(device_query_test PRIVATE cuda_std_11)
    set_target_properties(device_query_test PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        LANGUAGE CUDA
    )
    target_link_libraries(device_query_test 
        PRIVATE 
            cuda_settings
            GTest::gtest 
            GTest::gtest_main
    )
    gtest_discover_tests(device_query_test)

    # CUDA Test Suite
    add_executable(cuda_test cuda_test.cpp)
    target_link_libraries(cuda_test
        PRIVATE
            cuda_settings
            glog::glog
            GTest::gtest
            GTest::gtest_main
    )
    gtest_discover_tests(cuda_test)
endif()

# Google Logging Test
add_executable(glog_test glog_test.cpp)
target_link_libraries(glog_test
    PRIVATE
        glog::glog
        GTest::gtest
        GTest::gtest_main
)
gtest_discover_tests(glog_test)

# Python Lorenz95 Test
add_executable(lorenz95_python_test lorenz95_python_test.cpp)
target_compile_definitions(lorenz95_python_test 
    PRIVATE 
        PY_SSIZE_T_CLEAN
)
target_link_libraries(lorenz95_python_test
    PRIVATE
        Python3::Python
        Python3::NumPy
        GTest::gtest
        GTest::gtest_main
)
gtest_discover_tests(lorenz95_python_test)

# Copy Python module to build directory
configure_file(lorenz95.py ${CMAKE_CURRENT_BINARY_DIR}/lorenz95.py COPYONLY) 