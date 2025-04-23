# Findnetcdf-cxx4.cmake
#
# Finds the netCDF-C++ 4 library
#
# This will define the following variables:
#   netcdf-cxx4_FOUND        - True if the system has the netCDF-C++ 4 library
#   NETCDF_CXX4_INCLUDE_DIRS - The netCDF-C++ 4 include directory
#   NETCDF_CXX4_LIBRARIES    - The netCDF-C++ 4 libraries
#   NETCDF_CXX4_VERSION      - The netCDF-C++ 4 version (if available)
#
# and the following imported targets:
#   NetCDF::CXX4             - The netCDF-C++ 4 library target

# Find include directory
find_path(NETCDF_CXX4_INCLUDE_DIR
  NAMES ncFile.h
  PATH_SUFFIXES netcdf netcdf-cxx4 netcdf-cxx
  DOC "netCDF-C++ 4 include directory")
mark_as_advanced(NETCDF_CXX4_INCLUDE_DIR)

# Find library
find_library(NETCDF_CXX4_LIBRARY
  NAMES netcdf-cxx4 netcdf_c++4
  DOC "netCDF-C++ 4 library")
mark_as_advanced(NETCDF_CXX4_LIBRARY)

# Extract version information from the header file if available
if(NETCDF_CXX4_INCLUDE_DIR)
  # First try to extract version from ncException.h
  if(EXISTS "${NETCDF_CXX4_INCLUDE_DIR}/ncException.h")
    file(STRINGS "${NETCDF_CXX4_INCLUDE_DIR}/ncException.h" _netcdf_version_lines
      REGEX "#define[ \t]+NCXX4_VERSION_(MAJOR|MINOR|PATCH)")
    
    if(_netcdf_version_lines)
      string(REGEX REPLACE ".*NCXX4_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" _NETCDF_CXX4_VERSION_MAJOR "${_netcdf_version_lines}")
      string(REGEX REPLACE ".*NCXX4_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" _NETCDF_CXX4_VERSION_MINOR "${_netcdf_version_lines}")
      string(REGEX REPLACE ".*NCXX4_VERSION_PATCH[ \t]+([0-9]+).*" "\\1" _NETCDF_CXX4_VERSION_PATCH "${_netcdf_version_lines}")
      set(NETCDF_CXX4_VERSION "${_NETCDF_CXX4_VERSION_MAJOR}.${_NETCDF_CXX4_VERSION_MINOR}.${_NETCDF_CXX4_VERSION_PATCH}")
    endif()
  endif()
  
  # Try to get version using a different approach if not found yet
  if(NOT NETCDF_CXX4_VERSION)
    # Try to run netcdf-cxx4-config if available
    find_program(NETCDF_CXX4_CONFIG ncxx4-config HINTS "${NETCDF_CXX4_INCLUDE_DIR}/../bin")
    mark_as_advanced(NETCDF_CXX4_CONFIG)
    if(NETCDF_CXX4_CONFIG)
      execute_process(
        COMMAND ${NETCDF_CXX4_CONFIG} --version
        OUTPUT_VARIABLE NETCDF_CXX4_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    endif()
  endif()
  
  # Set a default version if still not found
  if(NOT NETCDF_CXX4_VERSION)
    set(NETCDF_CXX4_VERSION "unknown")
  endif()
endif()

# Handle standard FIND_PACKAGE arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(netcdf-cxx4
  REQUIRED_VARS
    NETCDF_CXX4_LIBRARY
    NETCDF_CXX4_INCLUDE_DIR
  VERSION_VAR NETCDF_CXX4_VERSION
)

# Set output variables
if(netcdf-cxx4_FOUND)
  set(NETCDF_CXX4_INCLUDE_DIRS ${NETCDF_CXX4_INCLUDE_DIR})
  set(NETCDF_CXX4_LIBRARIES ${NETCDF_CXX4_LIBRARY})
  
  # Set variables for package info system
  set(netcdf-cxx4_VERSION ${NETCDF_CXX4_VERSION})
  
  # Get the installation directory
  get_filename_component(_netcdf_cxx4_library_dir "${NETCDF_CXX4_LIBRARY}" DIRECTORY)
  get_filename_component(_netcdf_cxx4_install_dir "${_netcdf_cxx4_library_dir}" DIRECTORY)
  set(_netcdf_cxx4_library_dir "${_netcdf_cxx4_library_dir}" CACHE INTERNAL "Internal variable for netcdf-cxx4 library directory")
  set(_netcdf_cxx4_install_dir "${_netcdf_cxx4_install_dir}" CACHE INTERNAL "Internal variable for netcdf-cxx4 install directory")
  
  # Set the primary path variable for the package info system
  set(netcdf-cxx4_DIR "${_netcdf_cxx4_install_dir}" CACHE PATH "netcdf-cxx4 installation directory")
  
  # Create imported target
  if(NOT TARGET NetCDF::CXX4)
    add_library(NetCDF::CXX4 UNKNOWN IMPORTED)
    set_target_properties(NetCDF::CXX4 PROPERTIES
      IMPORTED_LOCATION "${NETCDF_CXX4_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_CXX4_INCLUDE_DIR}"
    )
  endif()
endif() 