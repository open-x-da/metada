# Enable Fortran support
enable_language(Fortran)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable with Fortran source
add_executable(obs_processing 
    obs_processing.cpp
)

# Link with required libraries
target_link_libraries(obs_processing
    PRIVATE
        metada::base
        metada::framework::adapters
        metada::framework::runs
        metada::traits
        metada::backends::common::utils
        metada::backends::common::io
        metada::backends::common::observation
)

# Link BUFR library if available (required for Fortran BUFR functions)
# Also ensure Fortran runtime is linked when using libraries with Fortran code
if(WITH_BUFR AND TARGET bufr::bufr_4)
    # Ensure Fortran language is enabled for this target
    enable_language(Fortran)
    
    target_link_libraries(obs_processing
        PRIVATE
            bufr::bufr_4
    )
    
    # CMake should automatically link Fortran runtime when Fortran is enabled,
    # but we ensure it's linked by enabling the language for this target
else()
    # Warn if BUFR is not available but the executable might need it
    message(WARNING "obs_processing may require BUFR support. Ensure WITH_BUFR is ON and bufr::bufr_4 target exists.")
endif()

# Set output directory
set_target_properties(obs_processing
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# readpb_config requires BUFR library - only build when available
if(WITH_BUFR AND TARGET bufr::bufr_4)
    add_executable(readpb_config readpb_config.f)

    target_link_libraries(readpb_config
        PRIVATE
            metada::backends::common::io
    )

    set_target_properties(readpb_config
        PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            Fortran_COMPILE_FLAGS "-ffixed-form"
    )
else()
    message(STATUS "readpb_config executable skipped: BUFR support not available")
endif()

# Add format target to ensure consistent code style
AddFormatTarget(obs_processing ${CMAKE_CURRENT_SOURCE_DIR})