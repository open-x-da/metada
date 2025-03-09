# Configure default compiler flags based on compiler ID
function(configure_compiler_flags)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra" 
            CACHE STRING "C compiler flags" FORCE)
        
        # Basic C++ compiler flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra" 
            CACHE STRING "C++ compiler flags" FORCE)
        
        # Add C++20 specific flags
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            # GCC needs specific flags for concepts before GCC 10
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts" 
                    CACHE STRING "C++ compiler flags" FORCE)
            endif()
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            # Clang may need specific flags for concepts
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts-ts" 
                    CACHE STRING "C++ compiler flags" FORCE)
            endif()
        endif()
        
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall" 
            CACHE STRING "Fortran compiler flags" FORCE)
    endif()
    
    # MSVC specific flags for C++20 concepts
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        # MSVC requires version 19.23 or later for C++20 concepts
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.23)
            message(WARNING "MSVC version ${CMAKE_CXX_COMPILER_VERSION} may not fully support C++20 concepts. Consider upgrading to MSVC 19.23 or later.")
        endif()
        
        # Add any MSVC-specific flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:__cplusplus /permissive-" 
            CACHE STRING "C++ compiler flags" FORCE)
    endif()
    
    # Print compiler information
    message(STATUS "C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
    message(STATUS "C++ Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS}")
    message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
endfunction()

configure_compiler_flags()