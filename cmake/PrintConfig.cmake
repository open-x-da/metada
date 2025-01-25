# Function to configure print support
function(configure_print_support)
    # Helper function to print a section header
    function(_print_section_header text)
        message(STATUS "")
        message(STATUS "${text}")
        string(LENGTH "${text}" text_length)
        string(REPEAT "-" ${text_length} separator)
        message(STATUS "${separator}")
        message(STATUS "")
    endfunction()

    # Helper function to print package info
    function(_print_package_info package_name)
        if(${package_name}_FOUND)
            set(version "${${package_name}_VERSION}")
            set(path "")
            
            # Convert package name to uppercase for compatibility
            string(TOUPPER ${package_name} PACKAGE_UPPER)
            
            # Try to get package path from different common variables
            if(DEFINED ${package_name}_DIR)
                set(path "${${package_name}_DIR}")
            elseif(DEFINED ${package_name}_ROOT)
                set(path "${${package_name}_ROOT}")
            # Special case for CUDA Toolkit path
            elseif(DEFINED ${package_name}_BIN_DIR)
                set(path "${${package_name}_BIN_DIR}")
            elseif(DEFINED ${PACKAGE_UPPER}_EXECUTABLE)
                set(path "${${PACKAGE_UPPER}_EXECUTABLE}")
            # Special case for Python3 executable path
            elseif(DEFINED _${package_name}_EXECUTABLE)
                set(path "${_${package_name}_EXECUTABLE}")
            elseif(DEFINED ${package_name}_PATH)
                set(path "${${package_name}_PATH}")
            endif()
            
            if(version AND path)
                message(STATUS "  - ${package_name}: Found (${version})")
                message(STATUS "    Path: ${path}")
            elseif(version)
                message(STATUS "  - ${package_name}: Found (${version})")
            elseif(path)
                message(STATUS "  - ${package_name}: Found")
                message(STATUS "    Path: ${path}")
            else()
                message(STATUS "  - ${package_name}: Found")
            endif()

            # Print component status if any components were requested
            get_property(components GLOBAL PROPERTY ${package_name}_COMPONENTS)
            if(components)
                message(STATUS "    Components:")
                foreach(component ${components})
                    if(${package_name}_${component}_FOUND)
                        message(STATUS "      - ${component}: Found")
                    else()
                        message(STATUS "      - ${component}: Not found")
                    endif()
                endforeach()
            endif()
        else()
            message(STATUS "  - ${package_name}: Not found")
        endif()
    endfunction()

    # Print configuration summary
    _print_section_header("Build Configuration")
    
    # Print compiler info with full paths
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
    
    # Print build information
    message(STATUS "Build Information:")
    message(STATUS "  Build Configuration:")
    message(STATUS "    - Type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "    - Directory: ${CMAKE_BINARY_DIR}")
    message(STATUS "    - Install prefix: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "  Compiler Options:")
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

    message(STATUS "  System Information:")
    message(STATUS "    - OS: ${CMAKE_HOST_SYSTEM_NAME} ${CMAKE_HOST_SYSTEM_VERSION}")
    message(STATUS "    - Architecture: ${CMAKE_HOST_SYSTEM_PROCESSOR}")

    message(STATUS "  CMake Information:")
    message(STATUS "    - Version: ${CMAKE_VERSION}")
    message(STATUS "    - Path: ${CMAKE_COMMAND}")
    
    # Print required packages
    _print_section_header("Required Packages")
    foreach(pkg IN LISTS METADA_REQUIRED_PACKAGES)
        _print_package_info(${pkg})
    endforeach()
    
    # Print optional packages
    _print_section_header("Optional Packages")
    foreach(pkg IN LISTS METADA_OPTIONAL_PACKAGES)
        _print_package_info(${pkg})
    endforeach()
endfunction()

# Configure print support immediately when this module is included
configure_print_support() 