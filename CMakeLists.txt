# Require CMake 3.30 or higher
cmake_minimum_required(VERSION 3.30)

# Add project-specific CMake modules
list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts"
)

# Define project metadata
project(metada
    VERSION 0.1.0
    LANGUAGES C CXX Fortran
)

# Load project configuration
include(ProjectSetup)
metada_project_initialize()

# Required dependencies
metada_find_package(Python3 COMPONENTS Development NumPy)
metada_find_package(GTest)
metada_find_package(Git)
metada_find_package(ClangFormat)
metada_find_package(Lcov)

# Optional dependencies
metada_find_package(yaml-cpp OPTIONAL QUIET)
metada_find_package(nlohmann_json OPTIONAL QUIET)
metada_find_package(Glog OPTIONAL QUIET)

# CUDA configuration
include(CUDA)

# Project structure
add_subdirectory(src)
include(CTest)
add_subdirectory(tests)
add_subdirectory(applications)
add_subdirectory(docs)

# Display build configuration
metada_project_summary()