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
        message(STATUS "  - Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
        message(STATUS "    Path: ${CMAKE_Fortran_COMPILER}")
    endif()

    if(CMAKE_CUDA_COMPILER)
        message(STATUS "  - CUDA compiler: ${CMAKE_CUDA_COMPILER_ID} ${CMAKE_CUDA_COMPILER_VERSION}")
        message(STATUS "    Path: ${CMAKE_CUDA_COMPILER}")
    endif()
endfunction() 