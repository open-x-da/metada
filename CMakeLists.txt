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

# Find optional packages
metada_find_package(yaml-cpp OPTIONAL QUIET)
metada_find_package(nlohmann_json OPTIONAL QUIET)
metada_find_package(Glog OPTIONAL QUIET) 

# Ensure at least one config backend is available
if(NOT (yaml-cpp_FOUND OR nlohmann_json_FOUND))
    message(FATAL_ERROR "At least one configuration backend (yaml-cpp or nlohmann_json) must be available")
endif()

# Configure CUDA support
include(CUDA)

# Add subdirectories
add_subdirectory(src)
include(CTest)
add_subdirectory(tests)
add_subdirectory(applications)
add_subdirectory(docs)

# Print configuration summary
metada_project_summary()