# Find all .gcda files in the build directory
file(GLOB_RECURSE GCDA_FILES "${CMAKE_BINARY_DIR}/*.gcda")

# Check if we're in a debug build configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(IS_DEBUG_BUILD TRUE)
else()
    set(IS_DEBUG_BUILD FALSE)
endif()

foreach(GCDA_FILE ${GCDA_FILES})
    if(IS_DEBUG_BUILD)
        message(STATUS "Found .gcda file: ${GCDA_FILE}")
    endif()
endforeach()

# Delete all .gcda files
foreach(GCDA_FILE ${GCDA_FILES})
    if(EXISTS "${GCDA_FILE}")
        file(REMOVE "${GCDA_FILE}")
        if(IS_DEBUG_BUILD)
            message(STATUS "Deleted: ${GCDA_FILE}")
        endif()
    endif()
endforeach()

if(IS_DEBUG_BUILD)
    message(STATUS "Cleaned ${CMAKE_BINARY_DIR} of .gcda files")
endif()