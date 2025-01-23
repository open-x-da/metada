cmake_minimum_required(VERSION 3.10)

# Forbid in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Please create a 'build' directory and run CMake from there.")
endif()

project(metada VERSION 0.1.0 LANGUAGES C CXX Fortran)

# Build options
option(USE_CUDA "Enable CUDA support" ON)
option(USE_DOXYGEN "Enable Doxygen documentation" ON)
option(USE_SPHINX "Enable Sphinx documentation" ON)

# Add custom CMake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(GTest REQUIRED)
find_package(Python3 COMPONENTS Development NumPy REQUIRED)
find_package(glog REQUIRED)

if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

# Add source directory
add_subdirectory(src)

# Configure CUDA
include(CUDA)
configure_cuda()

# Enable testing before adding subdirectory
enable_testing()
include(CTest)
add_subdirectory(tests)

# Configure documentation
include(Documentation)
configure_documentation()

# Include package configuration module
include(PackageConfig)

# Register required packages
register_required_package(GTest)
register_required_package(Python3)
register_required_package(glog)

# Register optional packages based on build options
if(USE_CUDA)
    register_optional_package(CUDAToolkit)
endif()

if(USE_DOXYGEN)
    register_optional_package(Doxygen)
endif()

if(USE_SPHINX)
    register_optional_package(Sphinx)
endif()

# Print configuration summary
include(PrintConfig)
print_configuration_summary()
