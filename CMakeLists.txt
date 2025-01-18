cmake_minimum_required(VERSION 3.10)

# Forbid in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Please create a 'build' directory and run CMake from there.")
endif()

project(cuda_learning VERSION 0.1.0 LANGUAGES C CXX Fortran)

# Find packages
find_package(CUDAToolkit REQUIRED)
find_package(GTest REQUIRED)
find_package(Python3 COMPONENTS Development NumPy REQUIRED)
find_package(Doxygen OPTIONAL_COMPONENTS dot)

# Add custom CMake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

# Main executable
add_executable(cuda_learning src/main.cpp)
target_link_libraries(cuda_learning PRIVATE CUDA::cudart)

# Enable testing before adding subdirectory
enable_testing()
include(CTest)
add_subdirectory(tests)

# Configure documentation
include(Documentation)
configure_documentation()
