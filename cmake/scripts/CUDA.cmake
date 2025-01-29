# Option to enable/disable CUDA support
option(ENABLE_CUDA_ACCELERATION "Enable GPU acceleration using CUDA (requires NVIDIA GPU and CUDA Toolkit)" ON)

# Function to configure CUDA support
function(configure_cuda_support)
    if(ENABLE_CUDA_ACCELERATION)
        metada_find_package(CUDAToolkit
            OPTIONAL
            QUIET
            CONDITION ENABLE_CUDA_ACCELERATION
        )

        if(NOT CUDAToolkit_FOUND)
            message(STATUS "CUDA acceleration disabled: CUDA Toolkit not found in PATH")
            set(ENABLE_CUDA_ACCELERATION OFF CACHE BOOL "Enable GPU acceleration using CUDA" FORCE)
            return()
        endif()

        # Create an interface library for CUDA settings
        add_library(cuda_settings INTERFACE)
        target_link_libraries(cuda_settings INTERFACE CUDA::cudart)
        target_compile_definitions(cuda_settings INTERFACE USE_CUDA)
        
        # Make it available to parent scope
        set(ENABLE_CUDA_ACCELERATION ON PARENT_SCOPE)
    endif()
endfunction()

# Configure CUDA support immediately when this module is included
configure_cuda_support() 