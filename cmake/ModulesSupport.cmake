# Configure C++20 modules support
#
# This module provides functions to work with C++20 modules in the project.
# It sets up the necessary compiler flags and build rules for modules.

# Check for modules support in the compiler
function(check_modules_support)
    # Check compiler version for modules support
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.0)
            message(WARNING "GCC ${CMAKE_CXX_COMPILER_VERSION} has limited support for C++20 modules. GCC 11 or later recommended.")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0)
            message(WARNING "Clang ${CMAKE_CXX_COMPILER_VERSION} has limited support for C++20 modules. Clang 14 or later recommended.")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.29)
            message(WARNING "MSVC ${CMAKE_CXX_COMPILER_VERSION} has limited support for C++20 modules. MSVC 19.29 or later recommended.")
        endif()
    endif()

    # Log the status of modules support
    message(STATUS "Configuring C++20 modules support for ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
    
    # Setup GCC module compilation settings
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        # Add just -fmodules-ts flag globally (other flags may not be supported)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-fmodules-ts>)
    endif()
endfunction()

# Function to add module compiler flags
# This function detects if a target is an INTERFACE library and adds flags appropriately
function(add_module_compiler_flags target)
    # Check if the target is an INTERFACE library
    get_target_property(target_type ${target} TYPE)
    if(target_type STREQUAL "INTERFACE_LIBRARY")
        set(scope_keyword "INTERFACE")
    else()
        set(scope_keyword "PRIVATE")
    endif()
    
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        # We've already added the -fmodules-ts flag globally in check_modules_support
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${target} ${scope_keyword} 
            -fmodules 
            -std=c++20 
            -fbuiltin-module-map)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        target_compile_options(${target} ${scope_keyword} /experimental:module /std:c++20)
    endif()
endfunction()

# Add module extension properties based on compiler
function(set_module_file_properties files)
    foreach(file ${files})
        if(file MATCHES ".*\.cppm$|.*\.ixx$")
            set_source_files_properties(${file} PROPERTIES
                LANGUAGE CXX
            )
            
            if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
                set_source_files_properties(${file} PROPERTIES
                    COMPILE_OPTIONS "/interface"
                )
            elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                set_source_files_properties(${file} PROPERTIES
                    COMPILE_OPTIONS "-xc++-module"
                )
            endif()
        endif()
    endforeach()
endfunction()

# Configure module mapping
function(configure_module_mapping)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # For clang, create module.modulemap in the build directory
        file(WRITE ${CMAKE_BINARY_DIR}/module.modulemap "")
    endif()
endfunction()

# Run checks when the module is included
check_modules_support()
configure_module_mapping() 