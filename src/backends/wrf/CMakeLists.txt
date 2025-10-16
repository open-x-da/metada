# Find required dependencies
metada_find_package(xtl CONFIG)
metada_find_package(xtensor CONFIG)
metada_find_package(netcdf-cxx4)

# WRF Configuration Bridge Library
# This library provides C++ interface to WRF's Fortran configuration system
# Split into two libraries to ensure proper build order: Fortran first, then C++

# Step 1: Build Fortran module first
add_library(wrf_config_bridge_fortran STATIC
    wrf_config_bridge.F90
)

# Step 2: Build C++ wrapper that depends on Fortran module
add_library(wrf_config_bridge STATIC
    WRFConfigBridge.cpp
)

# Make C++ compilation depend on Fortran module being built first
add_dependencies(wrf_config_bridge wrf_config_bridge_fortran)

# Set include directories for the config bridge
target_include_directories(wrf_config_bridge
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Also set include directories for Fortran library
target_include_directories(wrf_config_bridge_fortran
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

# WRF linking configuration (reuse WRFDA settings)
option(WITH_WRFDA "Enable linking real WRFDA obs operator routines" ON)
set(WRFDA_MOD_DIR "" CACHE PATH "Path to WRFDA build directory with .mod files")

if(WITH_WRFDA AND WRFDA_MOD_DIR AND EXISTS ${WRFDA_MOD_DIR})
    message(STATUS "WRF config bridge: using WRFDA module dir ${WRFDA_MOD_DIR}")
    
    # Add WRFDA module directory for Fortran .mod files
    target_include_directories(wrf_config_bridge_fortran PUBLIC ${WRFDA_MOD_DIR})
    target_include_directories(wrf_config_bridge PUBLIC ${WRFDA_MOD_DIR})
    
    # Get WRF root directory
    get_filename_component(WRF_ROOT_DIR "${WRFDA_MOD_DIR}/../.." ABSOLUTE)
    message(STATUS "WRF config bridge: WRF root directory: ${WRF_ROOT_DIR}")
    
    # Define WRF library paths needed for module_configure
    # Note: Order matters for static linking - dependencies must come after dependents
    set(WRF_CONFIG_LIBRARIES
        "${WRFDA_MOD_DIR}/libwrfvar.a"
        "${WRF_ROOT_DIR}/external/io_netcdf/libwrfio_nf.a"
        "${WRF_ROOT_DIR}/external/io_grib1/libio_grib1.a"
        "${WRF_ROOT_DIR}/external/io_grib_share/libio_grib_share.a"
        "${WRF_ROOT_DIR}/external/io_int/libwrfio_int.a"
        "${WRF_ROOT_DIR}/frame/pack_utils.o"
        "${WRF_ROOT_DIR}/external/esmf_time_f90/libesmf_time.a"
        "${WRF_ROOT_DIR}/var/build/libesmf_time.a"
    )
    
    # Link WRF libraries to Fortran library
    foreach(LIB_PATH ${WRF_CONFIG_LIBRARIES})
        if(EXISTS ${LIB_PATH})
            message(STATUS "WRF config bridge: linking library: ${LIB_PATH}")
            target_link_libraries(wrf_config_bridge_fortran PUBLIC ${LIB_PATH})
        else()
            message(WARNING "WRF config bridge: library not found: ${LIB_PATH}")
        endif()
    endforeach()
    
    # Link NetCDF Fortran library (required by module_configure)
    find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff)
    if(NETCDF_FORTRAN_LIBRARY)
        message(STATUS "WRF config bridge: linking NetCDF Fortran: ${NETCDF_FORTRAN_LIBRARY}")
        target_link_libraries(wrf_config_bridge_fortran PUBLIC ${NETCDF_FORTRAN_LIBRARY})
    endif()
    
    if(WIN32)
        target_link_libraries(wrf_config_bridge_fortran PUBLIC ws2_32)
    endif()
    
    # Link C++ library to Fortran library
    target_link_libraries(wrf_config_bridge PUBLIC wrf_config_bridge_fortran)
else()
    message(WARNING "WRF config bridge: WRFDA not configured, WRF configuration initialization will not be available")
endif()

# Create the WRF backend library (header-only)
add_library(wrf_backends INTERFACE)
add_library(metada::backends::wrf ALIAS wrf_backends)

# Set include directories
target_include_directories(wrf_backends 
    INTERFACE 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link dependencies
target_link_libraries(wrf_backends 
    INTERFACE
        NetCDF::CXX4
        xtensor
        wrf_config_bridge
)

# Add format target
AddFormatTarget(wrf_backends ${CMAKE_CURRENT_SOURCE_DIR}) 
