add_executable(device_query_test device_query_test.cpp)
target_compile_features(device_query_test PRIVATE cuda_std_11)
set_target_properties(device_query_test PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    LANGUAGE CUDA
)

target_link_libraries(device_query_test 
    PRIVATE 
    CUDA::cudart
    GTest::gtest 
    GTest::gtest_main
)

# Add Python Lorenz95 test
add_executable(lorenz95_python_test lorenz95_python_test.cpp)

# Set Python-specific compile definitions
target_compile_definitions(lorenz95_python_test 
    PRIVATE 
    PY_SSIZE_T_CLEAN
)

# Link libraries
target_link_libraries(lorenz95_python_test
    PRIVATE
    Python3::Python
    Python3::NumPy
    GTest::gtest
    GTest::gtest_main
)

# Register tests
gtest_discover_tests(device_query_test)
gtest_discover_tests(lorenz95_python_test)

# Copy Python module to build directory
configure_file(lorenz95.py ${CMAKE_CURRENT_BINARY_DIR}/lorenz95.py COPYONLY) 