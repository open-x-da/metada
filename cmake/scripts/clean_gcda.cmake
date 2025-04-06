# Find all .gcda files in the build directory
file(GLOB_RECURSE GCDA_FILES "${CMAKE_BINARY_DIR}/*.gcda")

# Print found files for debugging
foreach(GCDA_FILE ${GCDA_FILES})
    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        message(STATUS "Found .gcda file: ${GCDA_FILE}")
    endif()
endforeach()

# Delete all .gcda files
foreach(GCDA_FILE ${GCDA_FILES})
    if(EXISTS "${GCDA_FILE}")
        file(REMOVE "${GCDA_FILE}")
        if(CMAKE_BUILD_TYPE MATCHES "Debug")
            message(STATUS "Deleted: ${GCDA_FILE}")
        endif()
    endif()
endforeach()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Cleaned ${CMAKE_BINARY_DIR} of .gcda files")
endif()