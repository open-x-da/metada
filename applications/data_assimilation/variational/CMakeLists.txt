# Variational applications require WRFDA support
if(WITH_WRFDA)
    # TL/AD Checks Application
    add_executable(tl_ad_checks tl_ad_checks.cpp)

    # Incremental Variational Application
    add_executable(incremental_variational incremental_variational.cpp)

    target_link_libraries(tl_ad_checks
      PRIVATE
          metada::traits
          metada::framework::adapters
          metada::base
          metada::framework::runs
          metada::backends::common::utils
          metada::backends::common::io
          metada::backends::common::observation
          metada::backends::common::obsoperator
          metada::bridges::wrfda
          metada::framework::algorithms
          metada::backends::wrf
    )

    target_link_libraries(incremental_variational
      PRIVATE
          metada::traits
          metada::framework::adapters
          metada::base
          metada::framework::runs
          metada::backends::common::utils
          metada::backends::common::observation
          metada::backends::common::io
          metada::backends::common::obsoperator
          metada::bridges::wrfda
          metada::framework::algorithms
          metada::backends::wrf
    )
    
    # Link MPI libraries when MPI is enabled (required for WRFDA parallel operations)
    # WRFDA requires both CXX and Fortran MPI bindings (Fortran for libwrfvar.a and librsl_lite.a)
    if(METADA_USE_MPI)
        target_link_libraries(tl_ad_checks PRIVATE MPI::MPI_CXX MPI::MPI_Fortran)
        target_link_libraries(incremental_variational PRIVATE MPI::MPI_CXX MPI::MPI_Fortran)
    endif()

    # If WRFDA is enabled and we have a complete build environment, use linker group to resolve circular deps
    # Only link WRFDA libraries if we're not using the stub implementation
    if(WRFDA_VAR_LIB AND WRFDA_COMPLETE_BUILD)
        if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            # Explicitly add da_control.o if present and not archived in libwrfvar
            set(WRFDA_DA_CONTROL_OBJ "${WRFDA_MOD_DIR}/da_control.o")
            if (EXISTS ${WRFDA_DA_CONTROL_OBJ})
                target_link_libraries(tl_ad_checks PRIVATE ${WRFDA_DA_CONTROL_OBJ})
                target_link_libraries(incremental_variational PRIVATE ${WRFDA_DA_CONTROL_OBJ})
            endif()
            
            # Add module_quilt_outbuf_ops.o (required for quilt I/O operations)
            # This provides retrieve_pieces_of_field, store_piece_of_field, add_to_bufsize_for_field
            set(WRFDA_QUILT_OPS_OBJ "${WRFDA_MOD_DIR}/module_quilt_outbuf_ops.o")
            if (EXISTS ${WRFDA_QUILT_OPS_OBJ})
                target_link_libraries(tl_ad_checks PRIVATE ${WRFDA_QUILT_OPS_OBJ})
                target_link_libraries(incremental_variational PRIVATE ${WRFDA_QUILT_OPS_OBJ})
            endif()
            
            # Link NetCDF Fortran library (required by WRFDA)
            # Use NETCDF_ROOT for all NetCDF libraries (C, CXX, Fortran)
            if(NOT NETCDF_ROOT)
                set(NETCDF_ROOT "$ENV{NETCDF_ROOT}")
            endif()
            find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff HINTS ${NETCDF_ROOT} PATH_SUFFIXES lib lib64)
            if(NETCDF_FORTRAN_LIBRARY)
                target_link_libraries(tl_ad_checks PRIVATE ${NETCDF_FORTRAN_LIBRARY})
                target_link_libraries(incremental_variational PRIVATE ${NETCDF_FORTRAN_LIBRARY})
            endif()
            
            # Set CRTM library path (used by both executables)
            set(CRTM_LIBRARY_PATH "${WRFDA_MOD_DIR}/../external/crtm_2.3.0/libsrc/libCRTM.a")
            
            # Set RSL library path (required for MPI parallel communication in WRFDA)
            if(METADA_USE_MPI)
                get_filename_component(WRF_ROOT_DIR "${WRFDA_MOD_DIR}/../.." ABSOLUTE)
                set(RSL_LIBRARY_PATH "${WRF_ROOT_DIR}/external/RSL_LITE/librsl_lite.a")
            endif()
            
            # Configure tl_ad_checks linking
            # Use linker group to handle circular dependencies between WRFDA libraries
            # Note: --start-group/--end-group are GNU ld-specific (Linux only); macOS ld doesn't support these flags
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT APPLE)
                target_link_options(tl_ad_checks PRIVATE -Wl,--start-group)
            endif()
            target_link_libraries(tl_ad_checks PRIVATE ${WRFDA_VAR_LIB})
            # Link CRTM library AFTER libwrfvar.a (since libwrfvar depends on CRTM)
            if(EXISTS ${CRTM_LIBRARY_PATH})
                target_link_libraries(tl_ad_checks PRIVATE ${CRTM_LIBRARY_PATH})
            endif()
            # Link RSL library when MPI is enabled (required for parallel communication)
            if(METADA_USE_MPI AND EXISTS ${RSL_LIBRARY_PATH})
                target_link_libraries(tl_ad_checks PRIVATE ${RSL_LIBRARY_PATH})
            endif()
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT APPLE)
                target_link_options(tl_ad_checks PRIVATE -Wl,--end-group)
            endif()
            
            # Configure incremental_variational linking
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT APPLE)
                target_link_options(incremental_variational PRIVATE -Wl,--start-group)
            endif()
            target_link_libraries(incremental_variational PRIVATE ${WRFDA_VAR_LIB})
            # Link CRTM library AFTER libwrfvar.a (since libwrfvar depends on CRTM)
            if(EXISTS ${CRTM_LIBRARY_PATH})
                target_link_libraries(incremental_variational PRIVATE ${CRTM_LIBRARY_PATH})
            endif()
            # Link RSL library when MPI is enabled (required for parallel communication)
            if(METADA_USE_MPI AND EXISTS ${RSL_LIBRARY_PATH})
                target_link_libraries(incremental_variational PRIVATE ${RSL_LIBRARY_PATH})
            endif()
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT APPLE)
                target_link_options(incremental_variational PRIVATE -Wl,--end-group)
            endif()
        endif()
    elseif(WRFDA_VAR_LIB)
        message(STATUS "Variational application: Using WRFDA stub implementation (complete WRFDA environment not available)")
    endif()

    # Set output directory
    set_target_properties(tl_ad_checks
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    set_target_properties(incremental_variational
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    # Add format target to ensure consistent code style across the project
    AddFormatTarget(tl_ad_checks ${CMAKE_CURRENT_SOURCE_DIR}) 
    AddFormatTarget(incremental_variational ${CMAKE_CURRENT_SOURCE_DIR})
else()
    message(STATUS "Variational applications skipped: WRFDA support not available (WITH_WRFDA=OFF)")
endif()