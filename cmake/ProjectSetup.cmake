# Function to initialize project build environment
function(metada_project_initialize)
    # Forbid in-source builds
    if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
        message(FATAL_ERROR "In-source builds are not allowed. Please create a 'build' directory and run CMake from there.")
    endif()

    # Set build type if not specified
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
    endif()

    # On older Apple systems (Clang < 9.0), the filesystem library needs to be explicitly linked
    # This is because std::filesystem was experimental before C++17 and required separate linking
    if(APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        link_libraries(c++fs)
    endif()

    # Add at the beginning of metada_project_initialize() function
    option(USE_GLOG "Enable Google logging backend" ON)
    
    # Option to control precompiled headers usage
    option(USE_PRECOMPILED_HEADERS "Enable precompiled headers for faster builds" ON)
    
    # Option to control unity builds
    option(USE_UNITY_BUILD "Enable unity builds for faster compilation" ON)
    option(UNITY_BUILD_BATCH_SIZE "Number of source files to batch in each unity source" 8)

    # Include compiler flags configuration
    include(CompilerFlags)

    # Include target property helpers
    include(TargetProperties)

    # Include code formatting configuration
    include(CodeFormat)

    # Include testing and coverage configuration
    include(TestingAndCoverage)
    
    # Include precompiled headers configuration
    include(PrecompiledHeaders)
    
    # Include unity build configuration
    include(UnityBuild)
    
    # Include package configuration module
    include(package/Config)
    
    # Include printing utilities for configuration summary
    include(print/Config)
endfunction()