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

# Add format target
add_format_target(lorenz95_python_test ${CMAKE_CURRENT_SOURCE_DIR}) 

# Copy Python module to build directory
configure_file(lorenz95.py ${CMAKE_CURRENT_BINARY_DIR}/lorenz95.py COPYONLY) 