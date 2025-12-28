# Find required dependencies
metada_find_package(netcdf-cxx4)

# Find optional xtl/xtensor dependencies (controlled by WITH_XTL option)
if(WITH_XTL)
    metada_find_package(xtl CONFIG)
    metada_find_package(xtensor CONFIG)
    message(STATUS "WRF backend: xtl/xtensor support enabled")
else()
    message(STATUS "WRF backend: xtl/xtensor support disabled")
endif()

# WRF Configuration Bridge Library
# This library provides C++ interface to WRF's Fortran configuration system
# Split into two libraries to ensure proper build order: Fortran first, then C++

# Step 1: Build Fortran module first
add_library(wrf_config_bridge_fortran STATIC
    wrf_config_bridge.F90
)

# Enable Fortran preprocessing (required for #ifndef directives)
# .F90 extension should enable this automatically, but we make it explicit
set_source_files_properties(wrf_config_bridge.F90 PROPERTIES
    Fortran_PREPROCESS ON
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
# WITH_WRFDA is now defined globally in root CMakeLists.txt
# Derive WRFDA_MOD_DIR from WRFDA_ROOT
if(WITH_WRFDA)
    if(WRFDA_ROOT)
        set(WRFDA_MOD_DIR "${WRFDA_ROOT}/var/build")
    else()
        message(FATAL_ERROR "WRF backend: WRFDA_ROOT must be set when WITH_WRFDA=ON")
    endif()
endif()

if(WITH_WRFDA AND EXISTS ${WRFDA_MOD_DIR})
    message(STATUS "WRF config bridge: using WRFDA module dir ${WRFDA_MOD_DIR}")
    
    # Add WRFDA module directory for Fortran .mod files
    target_include_directories(wrf_config_bridge_fortran PUBLIC ${WRFDA_MOD_DIR})
    target_include_directories(wrf_config_bridge PUBLIC ${WRFDA_MOD_DIR})
    
    # Get WRF root directory
    get_filename_component(WRF_ROOT_DIR "${WRFDA_MOD_DIR}/../.." ABSOLUTE)
    message(STATUS "WRF config bridge: WRF root directory: ${WRF_ROOT_DIR}")
    
    #============================================================================
    # Detect WRFDA's MPI configuration (STUBMPI)
    #============================================================================
    # STUBMPI is set during WRFDA compilation. We need to detect it to properly
    # compile METADA's Fortran code. We try multiple detection methods:
    # 1. Check user-provided option WRFDA_HAS_STUBMPI (manual override)
    # 2. Check WRFDA's configure.wrf file for STUBMPI flag (most reliable)
    # 3. Check if module_dm.mod exists (if not, likely STUBMPI)
    # 4. Use METADA_USE_MPI as a hint (but warn if mismatch)
    
    # Allow manual override via CMake option
    # Check if WRFDA_HAS_STUBMPI was set via command line
    if(DEFINED WRFDA_HAS_STUBMPI)
        if(WRFDA_HAS_STUBMPI)
            message(STATUS "WRF config bridge: WRFDA_HAS_STUBMPI manually set to TRUE")
            set(WRFDA_HAS_STUBMPI "TRUE")
        else()
            message(STATUS "WRF config bridge: WRFDA_HAS_STUBMPI manually set to FALSE")
            set(WRFDA_HAS_STUBMPI "FALSE")
        endif()
    else()
        set(WRFDA_HAS_STUBMPI "UNKNOWN")
        
        # Method 1: Check configure.wrf for STUBMPI flag (most reliable)
        set(WRF_CONFIGURE_FILE "${WRF_ROOT_DIR}/configure.wrf")
        if(EXISTS ${WRF_CONFIGURE_FILE})
            file(READ ${WRF_CONFIGURE_FILE} WRF_CONFIGURE_CONTENT)
            if(WRF_CONFIGURE_CONTENT MATCHES "-DSTUBMPI")
                message(STATUS "WRF config bridge: Found -DSTUBMPI in configure.wrf")
                set(WRFDA_HAS_STUBMPI "TRUE")
            else()
                if(WRF_CONFIGURE_CONTENT MATCHES "STUBMPI" AND NOT WRF_CONFIGURE_CONTENT MATCHES "DM_PARALLEL")
                    message(STATUS "WRF config bridge: Found STUBMPI (without DM_PARALLEL) in configure.wrf")
                    set(WRFDA_HAS_STUBMPI "TRUE")
                elseif(WRF_CONFIGURE_CONTENT MATCHES "DM_PARALLEL" AND NOT WRF_CONFIGURE_CONTENT MATCHES "STUBMPI")
                    message(STATUS "WRF config bridge: Found DM_PARALLEL without STUBMPI in configure.wrf")
                    set(WRFDA_HAS_STUBMPI "FALSE")
                endif()
            endif()
        endif()
        
        # Method 2: If still unknown, use METADA_USE_MPI as the source of truth
        # (METADA_USE_MPI is auto-detected from WRFDA configuration in root CMakeLists.txt)
        if(WRFDA_HAS_STUBMPI STREQUAL "UNKNOWN")
            if(METADA_USE_MPI)
                # METADA_USE_MPI is ON, so WRFDA should have MPI (no STUBMPI)
                message(STATUS "WRF config bridge: METADA_USE_MPI=ON, assuming WRFDA has MPI (no STUBMPI)")
                set(WRFDA_HAS_STUBMPI "FALSE")
            else()
                # METADA_USE_MPI is OFF, so WRFDA should have STUBMPI
                message(STATUS "WRF config bridge: METADA_USE_MPI=OFF, assuming WRFDA has STUBMPI")
                set(WRFDA_HAS_STUBMPI "TRUE")
            endif()
        endif()
    endif()
    
    # Set STUBMPI preprocessor definition for Fortran compilation
    # Note: We need to set this BEFORE the source file is processed
    # The error suggests WRFDA was built with STUBMPI (module_dm doesn't have real functions)
    # So we need to set STUBMPI to exclude the use statement
    # 
    # However, if METADA_USE_MPI is ON (auto-detected), WRFDA should have MPI.
    # The error might indicate a mismatch - let's be conservative and check the actual library
    if(WRFDA_HAS_STUBMPI STREQUAL "TRUE")
        # Set preprocessor definition for Fortran
        # Use both methods to ensure preprocessor sees it
        target_compile_definitions(wrf_config_bridge_fortran PRIVATE STUBMPI)
        # Also add directly to compiler flags for Fortran preprocessor
        target_compile_options(wrf_config_bridge_fortran PRIVATE "-DSTUBMPI")
        message(STATUS "WRF config bridge: Setting STUBMPI for Fortran compilation (WRFDA built without MPI)")
        message(STATUS "WRF config bridge: STUBMPI will be defined for Fortran preprocessor")
    else()
        # WRFDA should have MPI, but if we get a linking error, it might actually have STUBMPI
        # For now, don't set STUBMPI, but add a note
        message(STATUS "WRF config bridge: NOT setting STUBMPI (WRFDA built with MPI support)")
        message(STATUS "WRF config bridge: If you get 'wrf_set_dm_communicator not found' error,")
        message(STATUS "WRF config bridge: WRFDA might actually have STUBMPI. Set -DWRFDA_HAS_STUBMPI=TRUE to override.")
        # Verify consistency with METADA_USE_MPI (should match since it's auto-detected)
        if(NOT METADA_USE_MPI)
            message(WARNING "WRF config bridge: Inconsistency detected - WRFDA has MPI but METADA_USE_MPI is OFF")
            message(WARNING "WRF config bridge: This should not happen if auto-detection is working correctly")
        endif()
    endif()
    
    # Verify STUBMPI detection result
    message(STATUS "WRF config bridge: WRFDA_HAS_STUBMPI detection result: ${WRFDA_HAS_STUBMPI}")
    
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
    
    # Add RSL library when MPI is enabled (required for parallel communication)
    if(METADA_USE_MPI)
        list(APPEND WRF_CONFIG_LIBRARIES "${WRF_ROOT_DIR}/external/RSL_LITE/librsl_lite.a")
    endif()
    
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
    # Use NETCDF_ROOT for all NetCDF libraries (C, CXX, Fortran)
    if(NOT NETCDF_ROOT)
        set(NETCDF_ROOT "$ENV{NETCDF_ROOT}")
    endif()
    find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff HINTS ${NETCDF_ROOT} PATH_SUFFIXES lib lib64)
    if(NETCDF_FORTRAN_LIBRARY)
        message(STATUS "WRF config bridge: linking NetCDF Fortran: ${NETCDF_FORTRAN_LIBRARY}")
        target_link_libraries(wrf_config_bridge_fortran PUBLIC ${NETCDF_FORTRAN_LIBRARY})
    endif()
    
    if(WIN32)
        target_link_libraries(wrf_config_bridge_fortran PUBLIC ws2_32)
    endif()
    
    # Link C++ library to Fortran library
    target_link_libraries(wrf_config_bridge PUBLIC wrf_config_bridge_fortran)
    
    # Link MPI if enabled (auto-detected from WRFDA configuration)
    if(METADA_USE_MPI)
        target_link_libraries(wrf_config_bridge PUBLIC MPI::MPI_CXX)
        target_link_libraries(wrf_config_bridge_fortran PUBLIC MPI::MPI_Fortran)
        message(STATUS "WRF config bridge: Linked to MPI libraries (auto-detected from WRFDA)")
    endif()
else()
    message(WARNING "WRF config bridge: WRFDA not configured, WRF configuration initialization will not be available")
endif()

# WRFDA State Transforms Bridge Library
# This library provides C++ interface to WRFDA's state transformation routines
# (da_transform_xtoxa, da_transfer_xatowrf) for proven increment-to-state conversion

# Step 1: Build Fortran bridge module first
add_library(wrfda_state_transforms_bridge_fortran STATIC
    bridges/wrfda_state_transforms.F90
    bridges/wrfda_control_backend.F90
)

# Link WRFDA if available
if(WITH_WRFDA AND EXISTS ${WRFDA_MOD_DIR})
    # Convert to absolute path if relative
    get_filename_component(WRFDA_MOD_DIR_ABS "${WRFDA_MOD_DIR}" ABSOLUTE)
    target_include_directories(wrfda_state_transforms_bridge_fortran PUBLIC ${WRFDA_MOD_DIR_ABS})
    
    # Get WRF root directory (reuse from config bridge if available, otherwise calculate)
    if(NOT DEFINED WRF_ROOT_DIR)
        get_filename_component(WRF_ROOT_DIR "${WRFDA_MOD_DIR_ABS}/../.." ABSOLUTE)
        message(STATUS "WRFDA state transforms: calculated WRF root directory: ${WRF_ROOT_DIR}")
    endif()
    
    # Add WRF source directories for .inc files (needed by wrfda_control_backend.F90)
    # This allows includes like #include "da_dot.inc" instead of hardcoded paths
    set(WRF_DA_MINIMISATION_DIR "${WRF_ROOT_DIR}/var/da/da_minimisation")
    if(EXISTS ${WRF_DA_MINIMISATION_DIR})
        message(STATUS "WRFDA state transforms: adding WRF da_minimisation include dir: ${WRF_DA_MINIMISATION_DIR}")
        target_include_directories(wrfda_state_transforms_bridge_fortran PUBLIC ${WRF_DA_MINIMISATION_DIR})
    else()
        message(WARNING "WRFDA state transforms: WRF da_minimisation directory not found at ${WRF_DA_MINIMISATION_DIR}")
        message(WARNING "WRFDA state transforms: WRF_ROOT_DIR calculated as: ${WRF_ROOT_DIR}")
    endif()
    
    # Link same WRF libraries as config bridge
    foreach(LIB_PATH ${WRF_CONFIG_LIBRARIES})
        if(EXISTS ${LIB_PATH})
            target_link_libraries(wrfda_state_transforms_bridge_fortran PUBLIC ${LIB_PATH})
        endif()
    endforeach()
    
    # Link NetCDF Fortran
    if(NETCDF_FORTRAN_LIBRARY)
        target_link_libraries(wrfda_state_transforms_bridge_fortran PUBLIC ${NETCDF_FORTRAN_LIBRARY})
    endif()
    
    if(WIN32)
        target_link_libraries(wrfda_state_transforms_bridge_fortran PUBLIC ws2_32)
    endif()
endif()

# Step 2: Build C++ wrapper library
add_library(wrfda_state_transforms STATIC
    wrfda/WRFDAStateTransforms.cpp
)

# Make C++ compilation depend on Fortran module
add_dependencies(wrfda_state_transforms wrfda_state_transforms_bridge_fortran)

# Set include directories
target_include_directories(wrfda_state_transforms
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Link Fortran bridge
target_link_libraries(wrfda_state_transforms
    PUBLIC
        wrfda_state_transforms_bridge_fortran
)

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
        metada::base
        NetCDF::CXX4
        wrf_config_bridge
        wrfda_state_transforms
)

# Conditionally link xtensor if enabled
if(WITH_XTL)
    target_link_libraries(wrf_backends INTERFACE xtensor)
endif()

# Add format target
AddFormatTarget(wrf_backends ${CMAKE_CURRENT_SOURCE_DIR}) 
