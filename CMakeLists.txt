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

# If both are found, prioritize yaml-cpp
if(yaml-cpp_FOUND AND nlohmann_json_FOUND)
    set(CONFIG_BACKEND "YAML")
    message(STATUS "Both yaml-cpp and nlohmann_json found. Defaulting to yaml-cpp.")
elseif(yaml-cpp_FOUND)
    set(CONFIG_BACKEND "YAML")
    message(STATUS "Using yaml-cpp for configuration.")
elseif(nlohmann_json_FOUND)
    set(CONFIG_BACKEND "JSON")
    message(STATUS "Using nlohmann_json for configuration.")
else()
    message(FATAL_ERROR "Neither yaml-cpp nor nlohmann_json was found. At least one is required.")
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