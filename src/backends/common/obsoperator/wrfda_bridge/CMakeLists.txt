cmake_minimum_required(VERSION 3.20)

project(wrfda_bridge C Fortran)

add_library(wrfda_bridge STATIC)

set_target_properties(wrfda_bridge PROPERTIES
    OUTPUT_NAME wrfda_bridge
)

target_include_directories(wrfda_bridge
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
)

# WRFDA linking configuration
option(WITH_WRFDA "Enable linking real WRFDA obs operator routines" ON)
set(WRFDA_ROOT "" CACHE PATH "Path to WRF/var/da root directory (e.g., D:/linux/WRF/var/da)")

# Allow specifying the prebuilt WRFDA static library directly
set(WRFDA_VAR_LIB "" CACHE FILEPATH "Path to prebuilt libwrfvar.a (e.g., D:/linux/WRF/var/da/build/libwrfvar.a)")

if(WITH_WRFDA)
    # Only link WRFDA libraries if we have a complete build environment
    # This will be determined later in the conditional build logic

    # Real wrapper sources require WRFDA modules (.mod) to be found
    set(WRFDA_MOD_DIR "" CACHE PATH "Path to WRFDA build directory with .mod files (e.g., D:/linux/WRF/var/da/build)")
    if(WRFDA_MOD_DIR AND EXISTS ${WRFDA_MOD_DIR})
        message(STATUS "WRFDA bridge: using module dir ${WRFDA_MOD_DIR}")
        target_include_directories(wrfda_bridge PUBLIC ${WRFDA_MOD_DIR})
        
        # Add WRF inc directory for config_assigns.inc and other registry-generated includes
        get_filename_component(WRF_INC_DIR "${WRFDA_MOD_DIR}/../../inc" ABSOLUTE)
        if(EXISTS ${WRF_INC_DIR})
            message(STATUS "WRFDA bridge: adding WRF inc directory: ${WRF_INC_DIR}")
            target_include_directories(wrfda_bridge PUBLIC ${WRF_INC_DIR})
        else()
            message(WARNING "WRFDA bridge: WRF inc directory not found at ${WRF_INC_DIR}")
        endif()
        
        # Check if we have a complete WRFDA build environment
        set(WRFDA_COMPLETE_BUILD FALSE CACHE BOOL "Whether WRFDA has complete build environment")
        if(EXISTS "${WRFDA_MOD_DIR}/libwrfvar.a" AND EXISTS "${WRFDA_MOD_DIR}/da_control.o")
            # Try to find required dependencies
            find_package(netcdf-cxx4 QUIET)
            find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff)
            set(CRTM_LIBRARY_PATH "${WRFDA_MOD_DIR}/../external/crtm_2.3.0/libsrc/libCRTM.a")
            
            # Check for WRFDA libraries directly in WRF directory structure
            get_filename_component(WRF_ROOT_DIR "${WRFDA_MOD_DIR}/../.." ABSOLUTE)
            set(WRF_ESMF_TIME_LIB "${WRF_ROOT_DIR}/external/esmf_time_f90/libesmf_time.a")
            set(WRF_FFTPACK_LIB "${WRF_ROOT_DIR}/external/fftpack/fftpack5/libfftpack.a")
            set(WRF_IO_NETCDF_LIB "${WRF_ROOT_DIR}/external/io_netcdf/libwrfio_nf.a")
            set(WRF_CRTM_LIB "${WRF_ROOT_DIR}/var/external/crtm_2.3.0/libsrc/libCRTM.a")
            
            # Check if we have the essential WRFDA libraries
            if(netcdf-cxx4_FOUND AND NETCDF_FORTRAN_LIBRARY AND EXISTS ${CRTM_LIBRARY_PATH} AND EXISTS ${WRF_ESMF_TIME_LIB} AND EXISTS ${WRF_FFTPACK_LIB} AND EXISTS ${WRF_IO_NETCDF_LIB})
                set(WRFDA_COMPLETE_BUILD TRUE CACHE BOOL "Whether WRFDA has complete build environment" FORCE)
                message(STATUS "WRFDA bridge: Complete WRFDA build environment detected")
                message(STATUS "WRFDA bridge: Found all required WRFDA libraries in WRF directory structure")
            else()
                set(WRFDA_COMPLETE_BUILD FALSE CACHE BOOL "Whether WRFDA has complete build environment" FORCE)
                message(STATUS "WRFDA bridge: Incomplete WRFDA build environment - missing dependencies:")
                if(NOT netcdf-cxx4_FOUND)
                    message(STATUS "  - NetCDF-C++4 not found")
                endif()
                if(NOT NETCDF_FORTRAN_LIBRARY)
                    message(STATUS "  - NetCDF Fortran library not found")
                endif()
                if(NOT EXISTS ${CRTM_LIBRARY_PATH})
                    message(STATUS "  - CRTM library not found at ${CRTM_LIBRARY_PATH}")
                endif()
                if(NOT EXISTS ${WRF_ESMF_TIME_LIB})
                    message(STATUS "  - ESMF time library not found at ${WRF_ESMF_TIME_LIB}")
                endif()
                if(NOT EXISTS ${WRF_FFTPACK_LIB})
                    message(STATUS "  - FFT library not found at ${WRF_FFTPACK_LIB}")
                endif()
                if(NOT EXISTS ${WRF_IO_NETCDF_LIB})
                    message(STATUS "  - WRF I/O NetCDF library not found at ${WRF_IO_NETCDF_LIB}")
                endif()
                message(STATUS "WRFDA bridge: Using stub implementation")
            endif()
        endif()
        
        if(WRFDA_COMPLETE_BUILD)
            # Full WRFDA integration with all dependencies
            target_sources(wrfda_bridge PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/wrfda_dispatch.F90
            )
            
            # Set WRF root directory based on WRFDA_MOD_DIR
            get_filename_component(WRF_ROOT_DIR "${WRFDA_MOD_DIR}/../.." ABSOLUTE)
            message(STATUS "WRFDA bridge: WRF root directory: ${WRF_ROOT_DIR}")
            
            # Link WRFDA library and dependencies
            set(WRFDA_DA_CONTROL_OBJ "${WRFDA_MOD_DIR}/da_control.o")
            if (EXISTS ${WRFDA_DA_CONTROL_OBJ})
                message(STATUS "WRFDA bridge: linking extra object ${WRFDA_DA_CONTROL_OBJ}")
                target_link_libraries(wrfda_bridge PUBLIC ${WRFDA_DA_CONTROL_OBJ})
            endif()
            
            # Link the main WRFDA library
            set(WRFDA_LIBRARY_PATH "${WRFDA_MOD_DIR}/libwrfvar.a")
            if(EXISTS ${WRFDA_LIBRARY_PATH})
                message(STATUS "WRFDA bridge: linking WRFDA library: ${WRFDA_LIBRARY_PATH}")
                target_link_libraries(wrfda_bridge PUBLIC ${WRFDA_LIBRARY_PATH})
            endif()
            
            # Link all WRFDA-related libraries directly from WRF directory structure
            set(WRF_LIBRARIES
                "${WRF_ROOT_DIR}/external/esmf_time_f90/libesmf_time.a"
                "${WRF_ROOT_DIR}/external/fftpack/fftpack5/libfftpack.a"
                "${WRF_ROOT_DIR}/external/io_grib1/libio_grib1.a"
                "${WRF_ROOT_DIR}/external/io_grib_share/libio_grib_share.a"
                "${WRF_ROOT_DIR}/external/io_int/libwrfio_int.a"
                "${WRF_ROOT_DIR}/external/io_netcdf/libwrfio_nf.a"
                "${WRF_ROOT_DIR}/var/build/libesmf_time.a"
                "${WRF_ROOT_DIR}/var/external/bufr/libbufr.a"
                "${WRF_ROOT_DIR}/var/external/crtm_2.3.0/libsrc/libCRTM.a"
            )
            
            # Link additional WRF object files needed for linking
            # Note: module_domain.o, module_io.o, and module_machine.o are already in libwrfvar.a
            # Only include object files that are NOT in the main WRFDA library
            set(WRF_OBJECT_FILES
                "${WRF_ROOT_DIR}/frame/pack_utils.o"
                "${WRF_ROOT_DIR}/var/build/module_internal_header_util.o"
            )
            
            # Link each library if it exists
            foreach(LIB_PATH ${WRF_LIBRARIES})
                if(EXISTS ${LIB_PATH})
                    message(STATUS "WRFDA bridge: linking library: ${LIB_PATH}")
                    target_link_libraries(wrfda_bridge PUBLIC ${LIB_PATH})
                else()
                    message(WARNING "WRFDA bridge: library not found: ${LIB_PATH}")
                endif()
            endforeach()
            
            # Link additional WRF object files if they exist
            foreach(OBJ_PATH ${WRF_OBJECT_FILES})
                if(EXISTS ${OBJ_PATH})
                    message(STATUS "WRFDA bridge: linking object file: ${OBJ_PATH}")
                    target_link_libraries(wrfda_bridge PUBLIC ${OBJ_PATH})
                else()
                    message(WARNING "WRFDA bridge: object file not found: ${OBJ_PATH}")
                endif()
            endforeach()
            
            # Link NetCDF libraries (still use find_package for system libraries)
            target_link_libraries(wrfda_bridge PUBLIC NetCDF::CXX4)
            
            find_library(NETCDF_C_LIBRARY NAMES netcdf)
            if(NETCDF_C_LIBRARY)
                target_link_libraries(wrfda_bridge PUBLIC ${NETCDF_C_LIBRARY})
            endif()
            
            # WRFDA requires NetCDF Fortran library for Fortran I/O operations
            find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff)
            if(NETCDF_FORTRAN_LIBRARY)
                message(STATUS "WRFDA bridge: linking NetCDF Fortran library: ${NETCDF_FORTRAN_LIBRARY}")
                target_link_libraries(wrfda_bridge PUBLIC ${NETCDF_FORTRAN_LIBRARY})
            else()
                message(WARNING "WRFDA bridge: NetCDF Fortran library not found, WRFDA Fortran I/O may fail")
            endif()
            
            if(WIN32)
                target_link_libraries(wrfda_bridge PUBLIC ws2_32)
            endif()
            
        else()
            # Fast-fail: WRFDA dependencies are incomplete
            message(FATAL_ERROR "WRFDA bridge: Incomplete WRFDA build environment detected")
            message(FATAL_ERROR "WRFDA bridge: Missing required dependencies:")
            if(NOT netcdf-cxx4_FOUND)
                message(FATAL_ERROR "  - NetCDF-C++4 not found")
            endif()
            if(NOT NETCDF_FORTRAN_LIBRARY)
                message(FATAL_ERROR "  - NetCDF Fortran library not found")
            endif()
            if(NOT EXISTS ${CRTM_LIBRARY_PATH})
                message(FATAL_ERROR "  - CRTM library not found at ${CRTM_LIBRARY_PATH}")
            endif()
            if(NOT EXISTS ${WRF_ESMF_TIME_LIB})
                message(FATAL_ERROR "  - ESMF time library not found at ${WRF_ESMF_TIME_LIB}")
            endif()
            if(NOT EXISTS ${WRF_FFTPACK_LIB})
                message(FATAL_ERROR "  - FFT library not found at ${WRF_FFTPACK_LIB}")
            endif()
            if(NOT EXISTS ${WRF_IO_NETCDF_LIB})
                message(FATAL_ERROR "  - WRF I/O NetCDF library not found at ${WRF_IO_NETCDF_LIB}")
            endif()
            message(FATAL_ERROR "WRFDA bridge: Please ensure all WRFDA dependencies are properly installed and configured")
        endif()
        
    else()
        # WRFDA module directory not provided - continue without WRFDA integration
        message(STATUS "WRFDA bridge: WRFDA_MOD_DIR not set - skipping WRFDA integration")
        message(STATUS "WRFDA bridge: To enable WRFDA integration, set WRFDA_MOD_DIR to the WRFDA build directory")
        message(STATUS "WRFDA bridge: Example: -DWRFDA_MOD_DIR=/path/to/WRF/var/da/build")
        message(STATUS "WRFDA bridge: Building without WRFDA-specific features")
        # Add stub source file to ensure library can be created
        target_sources(wrfda_bridge PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/wrfda_bridge_stub.F90
        )
    endif()
else()
    # WRFDA is disabled - continue without WRFDA integration
    message(STATUS "WRFDA bridge: WITH_WRFDA is disabled - skipping WRFDA integration")
    message(STATUS "WRFDA bridge: To enable WRFDA integration, set -DWITH_WRFDA=ON")
    message(STATUS "WRFDA bridge: Building without WRFDA-specific features")
    # Add stub source file to ensure library can be created
    target_sources(wrfda_bridge PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/wrfda_bridge_stub.F90
    )
endif()

add_library(metada::bridges::wrfda ALIAS wrfda_bridge)


