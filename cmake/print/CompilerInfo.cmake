# Print compiler information
function(print_compiler_info)
    message(STATUS "Compiler Information:")
    
    if(CMAKE_C_COMPILER)
        message(STATUS "  - C compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
        message(STATUS "    Path: ${CMAKE_C_COMPILER}")
    endif()
    
    if(CMAKE_CXX_COMPILER)
        message(STATUS "  - C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
        message(STATUS "    Path: ${CMAKE_CXX_COMPILER}")
    endif()
    
    if(CMAKE_Fortran_COMPILER)
        # Get Fortran compiler version if not available through CMake variables
        if(NOT CMAKE_Fortran_COMPILER_VERSION)
            execute_process(
                COMMAND ${CMAKE_Fortran_COMPILER} --version
                OUTPUT_VARIABLE FORTRAN_VERSION_OUTPUT
                ERROR_VARIABLE FORTRAN_VERSION_ERROR
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            if(FORTRAN_VERSION_OUTPUT)
                string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" FORTRAN_VERSION "${FORTRAN_VERSION_OUTPUT}")
            endif()
        else()
            set(FORTRAN_VERSION ${CMAKE_Fortran_COMPILER_VERSION})
        endif()

        message(STATUS "  - Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${FORTRAN_VERSION}")
        message(STATUS "    Path: ${CMAKE_Fortran_COMPILER}")
    endif()

    if(CMAKE_CUDA_COMPILER)
        message(STATUS "  - CUDA compiler: ${CMAKE_CUDA_COMPILER_ID} ${CMAKE_CUDA_COMPILER_VERSION}")
        message(STATUS "    Path: ${CMAKE_CUDA_COMPILER}")
    endif()
endfunction() 