# Configure default compiler flags based on compiler ID
function(configure_compiler_flags)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra" 
            CACHE STRING "C compiler flags" FORCE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra" 
            CACHE STRING "C++ compiler flags" FORCE)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall" 
            CACHE STRING "Fortran compiler flags" FORCE)
    elseif(MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4" 
            CACHE STRING "C compiler flags" FORCE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4" 
            CACHE STRING "C++ compiler flags" FORCE)
    endif()
endfunction()

configure_compiler_flags()