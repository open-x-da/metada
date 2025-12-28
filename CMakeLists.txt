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
metada_find_package(Eigen3 CONFIG)

# Find optional external dependencies
metada_find_package(Lcov OPTIONAL QUIET)
metada_find_package(yaml-cpp OPTIONAL CONFIG QUIET)
metada_find_package(nlohmann_json OPTIONAL CONFIG QUIET)

metada_find_package(ng-log OPTIONAL CONFIG QUIET)

# Find xtl/xtensor libraries for array operations (optional)
option(WITH_XTL "Enable xtl/xtensor support for WRF backend array operations" OFF)
if(WITH_XTL)
    add_compile_definitions(METADA_WITH_XTL)
endif()

# Find Arrow/Parquet libraries for observation data I/O (optional)
option(WITH_ARROW "Enable Arrow/Parquet support for observation data I/O" OFF)
if(WITH_ARROW)
    metada_find_package(Arrow REQUIRED CONFIG)
    add_compile_definitions(METADA_WITH_ARROW)
endif()

# Find BUFR library for observation data I/O (optional)
option(WITH_BUFR "Enable BUFR support for observation data I/O" OFF)
if(WITH_BUFR)
    add_compile_definitions(METADA_WITH_BUFR)
endif()

# Configure WRFDA support (optional)
option(WITH_WRFDA "Enable WRFDA support for WRF backend" OFF)
if(WITH_WRFDA)
    add_compile_definitions(METADA_WITH_WRFDA)
    # Auto-detect MPI requirement from WRFDA configuration
    include(WRFDAConfig)
    detect_wrfda_mpi_config()
else()
    # Without WRFDA, default to serial execution
    set(METADA_USE_MPI OFF CACHE BOOL "MPI support (default: OFF)" FORCE)
    message(STATUS "METADA: WRFDA disabled, using serial execution")
endif()

# Configure MPI support (automatically determined from WRFDA if WITH_WRFDA is ON)
if(METADA_USE_MPI)
    # Find both CXX and Fortran MPI components (Fortran needed for WRFDA)
    find_package(MPI COMPONENTS CXX Fortran REQUIRED)
    add_compile_definitions(METADA_USE_MPI)
    message(STATUS "MPI support enabled (auto-detected from WRFDA)")
    message(STATUS "  MPI CXX compiler: ${MPI_CXX_COMPILER}")
    message(STATUS "  MPI CXX include dirs: ${MPI_CXX_INCLUDE_DIRS}")
    message(STATUS "  MPI CXX libraries: ${MPI_CXX_LIBRARIES}")
    message(STATUS "  MPI Fortran compiler: ${MPI_Fortran_COMPILER}")
    message(STATUS "  MPI Fortran include dirs: ${MPI_Fortran_INCLUDE_DIRS}")
    message(STATUS "  MPI Fortran libraries: ${MPI_Fortran_LIBRARIES}")
else()
    message(STATUS "MPI support disabled (serial execution only)")
endif()

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


# Show final configuration
metada_project_summary()
