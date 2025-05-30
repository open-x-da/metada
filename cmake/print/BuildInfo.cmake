# Print build configuration
function(print_build_info)
    message(STATUS "Build Information:")
    message(STATUS "  Build Configuration:")
    message(STATUS "    - Type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "    - Directory: ${CMAKE_BINARY_DIR}")
    message(STATUS "    - Install prefix: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "    - CMake Version: ${CMAKE_VERSION}")
    message(STATUS "    - CMake Path: ${CMAKE_COMMAND}")
    
    message(STATUS "  Build Acceleration:")
    if(USE_PRECOMPILED_HEADERS)
        message(STATUS "    - Precompiled headers: Enabled")
    else()
        message(STATUS "    - Precompiled headers: Disabled")
    endif()
    
    if(USE_UNITY_BUILD)
        message(STATUS "    - Unity builds: Enabled (batch size: ${UNITY_BUILD_BATCH_SIZE})")
    else()
        message(STATUS "    - Unity builds: Disabled")
    endif()
    
    message(STATUS "  Compiler Options:")
    print_compiler_flags()
endfunction()

function(print_compiler_flags)
    if(CMAKE_C_COMPILER)
        message(STATUS "    C Flags:")
        message(STATUS "      - Base flags: ${CMAKE_C_FLAGS}")
        message(STATUS "      - ${CMAKE_BUILD_TYPE} flags: ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
    endif()
    
    if(CMAKE_CXX_COMPILER)
        message(STATUS "    C++ Flags:")
        message(STATUS "      - Base flags: ${CMAKE_CXX_FLAGS}")
        message(STATUS "      - ${CMAKE_BUILD_TYPE} flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
    endif()
    
    if(CMAKE_Fortran_COMPILER)
        message(STATUS "    Fortran Flags:")
        message(STATUS "      - Base flags: ${CMAKE_Fortran_FLAGS}")
        message(STATUS "      - ${CMAKE_BUILD_TYPE} flags: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
    endif()

    if(CMAKE_CUDA_FLAGS)
        message(STATUS "    CUDA Flags:")
        message(STATUS "      - Base flags: ${CMAKE_CUDA_FLAGS}")
        message(STATUS "      - ${CMAKE_BUILD_TYPE} flags: ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE}}")
    endif()
endfunction() 