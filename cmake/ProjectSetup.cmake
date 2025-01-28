# Function to initialize project build environment
function(metada_project_initialize)
    # Forbid in-source builds
    if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
        message(FATAL_ERROR "In-source builds are not allowed. Please create a 'build' directory and run CMake from there.")
    endif()

    # Set build type if not specified
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    endif()

    # Include compiler flags configuration
    include(CompilerFlags)

    # Include code formatting configuration
    include(CodeFormat)

    # Include package configuration module
    include(package/Config)
    
    # Include printing utilities for configuration summary
    include(print/Config)
endfunction()