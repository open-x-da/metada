# WRFDAConfig.cmake
# 
# Helper module for detecting WRFDA configuration, particularly MPI settings.
# This module provides functions to auto-detect WRFDA's build configuration
# and set METADA's MPI usage accordingly.

#============================================================================
# Function: detect_wrfda_mpi_config
#
# Automatically detects WRFDA's MPI configuration from configure.wrf and sets
# METADA_USE_MPI accordingly. This ensures METADA's MPI usage matches WRFDA's
# build configuration for compatibility.
#
# Parameters:
#   WRFDA_ROOT - Path to WRF root directory (required)
#
# Side effects:
#   Sets METADA_USE_MPI CACHE variable to ON or OFF based on WRFDA configuration
#
# Usage:
#   detect_wrfda_mpi_config()
#============================================================================
function(detect_wrfda_mpi_config)
    if(NOT WRFDA_ROOT)
        message(WARNING "WRFDA_ROOT not set. Cannot auto-detect MPI configuration.")
        message(WARNING "Defaulting to serial execution (METADA_USE_MPI=OFF)")
        set(METADA_USE_MPI OFF CACHE BOOL "MPI support (default: OFF)" FORCE)
        return()
    endif()
    
    set(WRFDA_MOD_DIR "${WRFDA_ROOT}/var/build")
    get_filename_component(WRF_ROOT_DIR "${WRFDA_MOD_DIR}/../.." ABSOLUTE)
    
    # Detect WRFDA's MPI configuration from configure.wrf
    set(WRF_CONFIGURE_FILE "${WRF_ROOT_DIR}/configure.wrf")
    set(WRFDA_HAS_MPI FALSE)
    
    if(EXISTS ${WRF_CONFIGURE_FILE})
        file(READ ${WRF_CONFIGURE_FILE} WRF_CONFIGURE_CONTENT)
        # Check if WRFDA was built with MPI (DM_PARALLEL without STUBMPI)
        if(WRF_CONFIGURE_CONTENT MATCHES "DM_PARALLEL" AND NOT WRF_CONFIGURE_CONTENT MATCHES "STUBMPI")
            set(WRFDA_HAS_MPI TRUE)
            message(STATUS "WRFDA: Detected MPI build (DM_PARALLEL without STUBMPI)")
        elseif(WRF_CONFIGURE_CONTENT MATCHES "-DSTUBMPI")
            set(WRFDA_HAS_MPI FALSE)
            message(STATUS "WRFDA: Detected serial build (STUBMPI)")
        else()
            # Default: assume MPI if DM_PARALLEL is mentioned
            if(WRF_CONFIGURE_CONTENT MATCHES "DM_PARALLEL")
                set(WRFDA_HAS_MPI TRUE)
                message(STATUS "WRFDA: Assuming MPI build (DM_PARALLEL found)")
            else()
                set(WRFDA_HAS_MPI FALSE)
                message(STATUS "WRFDA: Assuming serial build (no DM_PARALLEL)")
            endif()
        endif()
    else()
        message(WARNING "WRFDA configure.wrf not found at ${WRF_CONFIGURE_FILE}")
        message(WARNING "Cannot auto-detect MPI configuration. Defaulting to serial.")
        set(WRFDA_HAS_MPI FALSE)
    endif()
    
    # Set METADA_USE_MPI based on WRFDA's configuration
    if(WRFDA_HAS_MPI)
        set(METADA_USE_MPI ON CACHE BOOL "MPI support (auto-detected from WRFDA)" FORCE)
        message(STATUS "METADA: Auto-enabling MPI to match WRFDA configuration")
    else()
        set(METADA_USE_MPI OFF CACHE BOOL "MPI support (auto-detected from WRFDA)" FORCE)
        message(STATUS "METADA: Auto-disabling MPI to match WRFDA serial configuration")
    endif()
endfunction()

