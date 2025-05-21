# Configure precompiled headers functionality
#
# This module provides functions to easily add precompiled headers to targets.
# Requires CMake 3.16 or newer.

# Make sure we have the minimum required version for precompiled headers (3.16+)
if(CMAKE_VERSION VERSION_LESS 3.16)
    message(WARNING "CMake version ${CMAKE_VERSION} is too old for precompiled headers. Minimum required is 3.16.")
    function(metada_add_precompiled_headers)
        # No-op for older CMake versions
    endfunction()
    return()
endif()

# Function to add common precompiled headers to a target
# Usage: metada_add_precompiled_headers(target [PUBLIC|PRIVATE|INTERFACE] header1 [header2...])
function(metada_add_precompiled_headers target)
    # Parse arguments
    cmake_parse_arguments(ARGS "" "" "PUBLIC;PRIVATE;INTERFACE" ${ARGN})
    
    # Apply precompiled headers based on scope
    if(ARGS_PUBLIC)
        target_precompile_headers(${target} PUBLIC ${ARGS_PUBLIC})
    endif()
    
    if(ARGS_PRIVATE)
        target_precompile_headers(${target} PRIVATE ${ARGS_PRIVATE})
    endif()
    
    if(ARGS_INTERFACE)
        target_precompile_headers(${target} INTERFACE ${ARGS_INTERFACE})
    endif()
endfunction()

# Function to add standard C++ headers as precompiled headers to a target
function(metada_add_std_precompiled_headers target scope)
    target_precompile_headers(${target} ${scope}
        <vector>
        <string>
        <map>
        <unordered_map>
        <memory>
        <utility>
        <algorithm>
        <functional>
        <iostream>
    )
endfunction()

# Function to apply precompiled headers to all targets in the current directory
# Usage: metada_apply_precompiled_headers_to_all_targets([TARGETS target1 target2...] [PROJECT_HEADERS header1 header2...])
function(metada_apply_precompiled_headers_to_all_targets)
    if(NOT USE_PRECOMPILED_HEADERS)
        return()
    endif()
    
    # Parse arguments
    cmake_parse_arguments(ARGS "" "" "TARGETS;PROJECT_HEADERS" ${ARGN})
    
    # If no targets specified, get all targets in the current directory
    if(NOT ARGS_TARGETS)
        get_property(ARGS_TARGETS DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
    endif()
    
    foreach(target ${ARGS_TARGETS})
        # Check if target is a compilable type
        get_target_property(target_type ${target} TYPE)
        if(NOT (target_type STREQUAL "STATIC_LIBRARY" OR
                target_type STREQUAL "SHARED_LIBRARY" OR
                target_type STREQUAL "MODULE_LIBRARY" OR
                target_type STREQUAL "OBJECT_LIBRARY" OR
                target_type STREQUAL "EXECUTABLE"))
            continue()
        endif()
        
        # Add standard C++ headers
        metada_add_std_precompiled_headers(${target} PRIVATE)
        
        # Add project-specific headers if provided
        if(ARGS_PROJECT_HEADERS)
            target_precompile_headers(${target} PRIVATE ${ARGS_PROJECT_HEADERS})
        endif()
    endforeach()
endfunction() 