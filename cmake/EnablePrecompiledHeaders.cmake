# Module to globally enable precompiled headers for the project
#
# Include this module in the main CMakeLists.txt after all targets are defined
# to automatically apply precompiled headers to appropriate targets.

include(PrecompiledHeaders)

# Skip if precompiled headers are disabled
if(NOT USE_PRECOMPILED_HEADERS)
    return()
endif()

message(STATUS "Enabling precompiled headers for faster compilation")

# Define common headers used across the project
set(METADA_COMMON_HEADERS
    <vector>
    <string>
    <memory>
    <utility>
    <algorithm>
    <functional>
    <iostream>
    <map>
    <unordered_map>
)

# Function to recursively find all non-interface library targets
function(find_all_compilable_targets dir result)
    set(targets)
    
    # Get all targets in this directory
    get_property(dir_targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
    foreach(target ${dir_targets})
        # Check if target is compilable
        get_target_property(target_type ${target} TYPE)
        if(target_type STREQUAL "STATIC_LIBRARY" OR
           target_type STREQUAL "SHARED_LIBRARY" OR
           target_type STREQUAL "MODULE_LIBRARY" OR
           target_type STREQUAL "OBJECT_LIBRARY" OR
           target_type STREQUAL "EXECUTABLE")
            list(APPEND targets ${target})
        endif()
    endforeach()
    
    # Get all subdirectories
    get_property(subdirs DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
    foreach(subdir ${subdirs})
        # Recursively get targets from subdirectories
        find_all_compilable_targets(${subdir} subdir_targets)
        list(APPEND targets ${subdir_targets})
    endforeach()
    
    # Return results
    set(${result} ${targets} PARENT_SCOPE)
endfunction()

# Find all compilable targets in the project
find_all_compilable_targets(${CMAKE_SOURCE_DIR} ALL_COMPILABLE_TARGETS)

# Apply precompiled headers to all compilable targets
foreach(target ${ALL_COMPILABLE_TARGETS})
    message(STATUS "Enabling precompiled headers for target: ${target}")
    target_precompile_headers(${target} PRIVATE ${METADA_COMMON_HEADERS})
endforeach() 