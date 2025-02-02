# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Add custom CMake modules paths
list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts"
)

# Initial project declaration with version and languages
project(metada 
    VERSION 0.1.0
    LANGUAGES C CXX Fortran
)

# Include project configuration module
include(ProjectSetup)

# Configure build environment
metada_project_initialize()

# Find required packages
metada_find_package(Python3 COMPONENTS Development NumPy)
metada_find_package(GTest)
metada_find_package(Git)
metada_find_package(ClangFormat)
metada_find_package(Lcov)

metada_find_package(Glog OPTIONAL)

# Configure CUDA support
include(CUDA)

# Add source directory
add_subdirectory(src)

# Add test suite directory
include(CTest)
add_subdirectory(tests)

# Add applications directory
add_subdirectory(applications)

# Add documentation directory
add_subdirectory(docs)

# Print configuration summary
metada_project_summary()