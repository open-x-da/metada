cmake_minimum_required(VERSION 3.10)

# Forbid in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Please create a 'build' directory and run CMake from there.")
endif()

project(metada VERSION 0.1.0 LANGUAGES C CXX Fortran)

# Make CUDA optional
option(USE_CUDA "Enable CUDA support" ON)

# Find required packages
if(USE_CUDA)
    find_package(CUDAToolkit QUIET)
    if(NOT CUDAToolkit_FOUND)
        message(STATUS "CUDA not found - disabling CUDA support")
        set(USE_CUDA OFF)
    endif()
endif()

find_package(GTest REQUIRED)
find_package(Python3 COMPONENTS Development NumPy REQUIRED)
find_package(Doxygen OPTIONAL_COMPONENTS dot)
find_package(glog REQUIRED)

# Add custom CMake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

# Add source directory
add_subdirectory(src)

# Enable testing before adding subdirectory
enable_testing()
include(CTest)
add_subdirectory(tests)

# Configure documentation
include(Documentation)
configure_documentation()
