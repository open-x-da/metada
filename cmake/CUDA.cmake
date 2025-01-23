# CUDA configuration function
function(configure_cuda)
    if(USE_CUDA)
        find_package(CUDAToolkit QUIET)
        if(NOT CUDAToolkit_FOUND)
            message(STATUS "CUDA not found - disabling CUDA support")
            set(USE_CUDA OFF PARENT_SCOPE)
            return()
        endif()

        # Create an interface library for CUDA settings
        add_library(cuda_settings INTERFACE)
        target_link_libraries(cuda_settings INTERFACE CUDA::cudart)
        target_compile_definitions(cuda_settings INTERFACE USE_CUDA)
        
        # Make it available to parent scope
        set(USE_CUDA ON PARENT_SCOPE)
    endif()
endfunction() 