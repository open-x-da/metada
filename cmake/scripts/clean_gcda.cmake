# Find all .gcda files in the build directory
file(GLOB_RECURSE GCDA_FILES "${CMAKE_BINARY_DIR}/*.gcda")

# Print found files for debugging
foreach(GCDA_FILE ${GCDA_FILES})
    message(STATUS "Found .gcda file: ${GCDA_FILE}")
endforeach()

# Delete all .gcda files
foreach(GCDA_FILE ${GCDA_FILES})
    if(EXISTS "${GCDA_FILE}")
        file(REMOVE "${GCDA_FILE}")
        message(STATUS "Deleted: ${GCDA_FILE}")
    endif()
endforeach()

message(STATUS "Cleaned ${CMAKE_BINARY_DIR} of .gcda files")