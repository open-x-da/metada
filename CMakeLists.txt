# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Initial project declaration (required by CMake)
project(metada)

# Add custom CMake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Include project configuration module
include(ProjectSetup)

# Configure project with additional settings
metada_project_settings(metada
    VERSION 0.1.0
    LANGUAGES C CXX Fortran
)

# Find required packages
metada_find_package(GTest)
metada_find_package(Python3 COMPONENTS Development NumPy)
metada_find_package(glog)

# Configure CUDA support
include(CUDA)

# Add source directory
add_subdirectory(src)

# Add format target
add_format_target()

# Enable testing
include(CTest)
add_subdirectory(tests)

# Add documentation directory
add_subdirectory(docs)

# Print configuration summary
metada_project_summary()