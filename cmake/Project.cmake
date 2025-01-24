# Function to configure project settings
function(configure_project_settings name)
    cmake_parse_arguments(ARG
        ""
        "VERSION"
        "LANGUAGES"
        ${ARGN}
    )

    # Set default languages if not specified
    if(NOT ARG_LANGUAGES)
        set(ARG_LANGUAGES C CXX)
    endif()

    # Set default version if not specified
    if(NOT ARG_VERSION)
        set(ARG_VERSION 0.1.0)
    endif()

    # Update project settings
    project(${name}
        VERSION ${ARG_VERSION}
        LANGUAGES ${ARG_LANGUAGES}
    )

    # Include package configuration module
    include(PackageConfig)
endfunction()

# Macro to wrap project configuration (using macro to allow direct project() call)
macro(metada_project name)
    cmake_parse_arguments(ARG
        ""
        "VERSION"
        "LANGUAGES"
        ${ARGN}
    )

    # Set default languages if not specified
    if(NOT ARG_LANGUAGES)
        set(ARG_LANGUAGES C CXX)
    endif()

    # Set default version if not specified
    if(NOT ARG_VERSION)
        set(ARG_VERSION 0.1.0)
    endif()

    # Forbid in-source builds
    if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
        message(FATAL_ERROR "In-source builds are not allowed. Please create a 'build' directory and run CMake from there.")
    endif()

    # Call project command directly
    project(${name}
        VERSION ${ARG_VERSION}
        LANGUAGES ${ARG_LANGUAGES}
    )

    # Configure additional project settings
    configure_project_settings(${name})
endmacro() 