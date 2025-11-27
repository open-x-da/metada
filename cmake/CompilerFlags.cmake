# Configure default compiler flags based on compiler ID
function(configure_compiler_flags)
    # Configure C++20 as the project standard
    set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to use" FORCE)
    set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require C++ standard to be supported" FORCE)
    set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Disable compiler-specific extensions" FORCE)
       
    # Add module support for C++20
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fmodules-ts>")
    endif()
    
    # Log C++ standard being used
    message(STATUS "Using C++ standard: C++20")

    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra" 
            CACHE STRING "C compiler flags" FORCE)
        
        # Basic C++ compiler flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra" 
            CACHE STRING "C++ compiler flags" FORCE)
        
        # Add C++20 specific flags if needed
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            # Check for C++20 support
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
                message(WARNING "GCC ${CMAKE_CXX_COMPILER_VERSION} has limited support for C++20. GCC 10 or later required for full C++20 support.")
            endif()
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            # Check for C++20 support
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
                message(WARNING "Clang ${CMAKE_CXX_COMPILER_VERSION} has limited support for C++20. Clang 10 or later required for full C++20 support.")
            endif()
        endif()
        
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall" 
            CACHE STRING "Fortran compiler flags" FORCE)
    endif()
    
    # Configure debug-specific Fortran compiler flags
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
            set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ggdb -O0 -fdefault-real-8 -fconvert=big-endian -frecord-marker=4 -fbacktrace -fcheck=bounds,do,mem,pointer -ffpe-trap=invalid,zero,overflow"
                CACHE STRING "Fortran debug compiler flags" FORCE)
            message(STATUS "Added debug Fortran flags: -ggdb -O0 -fdefault-real-8 -fconvert=big-endian -frecord-marker=4 -fbacktrace -fcheck=bounds,do,mem,pointer -ffpe-trap=invalid,zero,overflow")
        endif()
    endif()
    
endfunction()

configure_compiler_flags()