# Configure Unity Build functionality
#
# This module provides functions to easily enable Unity Builds for targets
# to reduce compilation time by combining source files.

# Function to enable unity build for a target
# Usage: metada_enable_unity_build(target [BATCH_SIZE num])
function(metada_enable_unity_build target)
    # Parse arguments
    cmake_parse_arguments(ARGS "" "BATCH_SIZE" "" ${ARGN})
    
    # Default batch size if not specified
    if(NOT ARGS_BATCH_SIZE)
        set(ARGS_BATCH_SIZE 8)
    endif()
    
    # Enable unity build
    set_target_properties(${target} PROPERTIES
        UNITY_BUILD ON
        UNITY_BUILD_BATCH_SIZE ${ARGS_BATCH_SIZE}
    )
    
    message(STATUS "Unity build enabled for target: ${target} (batch size: ${ARGS_BATCH_SIZE})")
endfunction()

# Function to enable unity build for all targets in a directory
# Usage: metada_enable_unity_build_for_directory([DIRECTORY dir] [BATCH_SIZE num] [RECURSIVE])
function(metada_enable_unity_build_for_directory)
    # Parse arguments
    cmake_parse_arguments(ARGS "RECURSIVE" "BATCH_SIZE;DIRECTORY" "" ${ARGN})
    
    # Default directory if not specified
    if(NOT ARGS_DIRECTORY)
        set(ARGS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    
    # Default batch size if not specified
    if(NOT ARGS_BATCH_SIZE)
        set(ARGS_BATCH_SIZE 8)
    endif()
    
    # Get all targets in the directory
    get_property(dir_targets DIRECTORY ${ARGS_DIRECTORY} PROPERTY BUILDSYSTEM_TARGETS)
    
    # Enable unity build for each compilable target
    foreach(target ${dir_targets})
        # Check if target is compilable
        get_target_property(target_type ${target} TYPE)
        if(target_type STREQUAL "STATIC_LIBRARY" OR
           target_type STREQUAL "SHARED_LIBRARY" OR
           target_type STREQUAL "MODULE_LIBRARY" OR
           target_type STREQUAL "OBJECT_LIBRARY")
            metada_enable_unity_build(${target} BATCH_SIZE ${ARGS_BATCH_SIZE})
        endif()
    endforeach()
    
    # Process subdirectories if RECURSIVE is specified
    if(ARGS_RECURSIVE)
        get_property(subdirectories DIRECTORY ${ARGS_DIRECTORY} PROPERTY SUBDIRECTORIES)
        foreach(subdir ${subdirectories})
            metada_enable_unity_build_for_directory(
                DIRECTORY ${subdir}
                BATCH_SIZE ${ARGS_BATCH_SIZE}
                RECURSIVE
            )
        endforeach()
    endif()
endfunction()

# Function to selectively disable unity build for specific source files
# Usage: metada_exclude_from_unity_build(target source1 [source2...])
function(metada_exclude_from_unity_build target)
    # Get the sources to exclude
    set(sources ${ARGN})
    
    # Skip if no sources provided
    if(NOT sources)
        return()
    endif()
    
    # Set properties for each source file
    foreach(source ${sources})
        set_source_files_properties(${source} PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
    endforeach()
    
    message(STATUS "Excluded ${CMAKE_CURRENT_LIST_LINE} source files from unity build for target: ${target}")
endfunction() 