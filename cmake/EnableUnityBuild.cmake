# Module to globally enable unity builds for the project
#
# Include this module in the main CMakeLists.txt after all targets are defined
# to automatically apply unity builds to appropriate targets.

include(UnityBuild)

# Option to control unity build usage
option(USE_UNITY_BUILD "Enable unity builds for faster compilation" ON)
option(UNITY_BUILD_BATCH_SIZE "Number of source files to batch in each unity source" 8)

# Skip if unity builds are disabled
if(NOT USE_UNITY_BUILD)
    return()
endif()

message(STATUS "Enabling unity builds for faster compilation (batch size: ${UNITY_BUILD_BATCH_SIZE})")

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
           target_type STREQUAL "OBJECT_LIBRARY")
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

# Apply unity builds to all compilable targets
foreach(target ${ALL_COMPILABLE_TARGETS})
    message(STATUS "Enabling unity build for target: ${target}")
    metada_enable_unity_build(${target} BATCH_SIZE ${UNITY_BUILD_BATCH_SIZE})
endforeach() 