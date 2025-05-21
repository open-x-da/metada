# Root CMakeLists.txt for the metada project
#
# This is the main CMake configuration file that:
# - Sets up project metadata and CMake requirements
# - Configures build options and dependencies
# - Defines the project structure
# - Enables testing and documentation

cmake_minimum_required(VERSION 3.30)

list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts"
)

project(metada
    VERSION 0.1.0
    LANGUAGES C CXX Fortran
)

include(ProjectSetup)
metada_project_initialize()

# Find required external dependencies
#metada_find_package(Python3 COMPONENTS Development NumPy)
metada_find_package(GTest CONFIG)
metada_find_package(Git)
metada_find_package(ClangFormat)


# Find optional external dependencies
metada_find_package(Lcov OPTIONAL QUIET)
metada_find_package(yaml-cpp OPTIONAL CONFIG QUIET)
metada_find_package(nlohmann_json OPTIONAL CONFIG QUIET)
metada_find_package(Glog OPTIONAL CONFIG QUIET)

# Configure CUDA support
include(CUDA)

# Add core components
add_subdirectory(src)
add_subdirectory(applications)

# Enable testing
include(CTest)
add_subdirectory(tests)

# Add documentation
add_subdirectory(docs)

# Enable precompiled headers for all targets
#include(EnablePrecompiledHeaders)

# Enable unity builds for all targets
#include(EnableUnityBuild)

# Show final configuration
metada_project_summary()
