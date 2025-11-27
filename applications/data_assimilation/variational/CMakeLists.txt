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
            
            # Link NetCDF Fortran library (required by WRFDA)
            find_library(NETCDF_FORTRAN_LIBRARY NAMES netcdff)
            if(NETCDF_FORTRAN_LIBRARY)
                target_link_libraries(tl_ad_checks PRIVATE ${NETCDF_FORTRAN_LIBRARY})
                target_link_libraries(incremental_variational PRIVATE ${NETCDF_FORTRAN_LIBRARY})
            endif()
            
            # Set CRTM library path (used by both executables)
            set(CRTM_LIBRARY_PATH "${WRFDA_MOD_DIR}/../external/crtm_2.3.0/libsrc/libCRTM.a")
            
            # Configure tl_ad_checks linking
            # Use linker group to handle circular dependencies between WRFDA libraries
            # Note: Library order matters - dependencies must come after dependents
            target_link_options(tl_ad_checks PRIVATE -Wl,--start-group)
            target_link_libraries(tl_ad_checks PRIVATE ${WRFDA_VAR_LIB})
            # Link CRTM library AFTER libwrfvar.a (since libwrfvar depends on CRTM)
            if(EXISTS ${CRTM_LIBRARY_PATH})
                target_link_libraries(tl_ad_checks PRIVATE ${CRTM_LIBRARY_PATH})
            endif()
            target_link_options(tl_ad_checks PRIVATE -Wl,--end-group)
            
            # Configure incremental_variational linking
            target_link_options(incremental_variational PRIVATE -Wl,--start-group)
            target_link_libraries(incremental_variational PRIVATE ${WRFDA_VAR_LIB})
            # Link CRTM library AFTER libwrfvar.a (since libwrfvar depends on CRTM)
            if(EXISTS ${CRTM_LIBRARY_PATH})
                target_link_libraries(incremental_variational PRIVATE ${CRTM_LIBRARY_PATH})
            endif()
            target_link_options(incremental_variational PRIVATE -Wl,--end-group)
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