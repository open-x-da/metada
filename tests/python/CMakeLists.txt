# Skip Python tests on Windows with MSYS2
# MSYS2 on Windows has compatibility issues with Python tests due to:
# 1. Path handling differences between Windows and MSYS2 environments
# 2. Potential conflicts between MSYS2's Python and system Python
# 3. Known issues with Python extension modules built under MSYS2
if(WIN32 AND MINGW)
    message(STATUS "Skipping Python tests on Windows with MSYS2")
    return()
endif()

# Python Lorenz95 Test
add_executable(lorenz95_python_test lorenz95_python_test.cpp)

# Define PY_SSIZE_T_CLEAN before including Python.h
# This ensures consistent behavior with Py_ssize_t for string lengths and array indices
# See: https://docs.python.org/3/c-api/intro.html#include-files
target_compile_definitions(lorenz95_python_test 
    PRIVATE 
        PY_SSIZE_T_CLEAN
)

# Link Python and NumPy libraries
target_link_libraries(lorenz95_python_test
    PRIVATE
        Python3::Python
        Python3::NumPy
)

# Add test setup and coverage
metada_add_test_with_coverage(lorenz95_python_test)

# Add format target
AddFormatTarget(lorenz95_python_test ${CMAKE_CURRENT_SOURCE_DIR}) 

# Copy Python module to build directory
configure_file(lorenz95.py ${CMAKE_CURRENT_BINARY_DIR}/lorenz95.py COPYONLY) 