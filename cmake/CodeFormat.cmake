# Find clang-format executable
find_program(CLANG_FORMAT_EXECUTABLE
    NAMES clang-format clang-format-16 clang-format-15 clang-format-14
    DOC "Path to clang-format executable"
)

# Function to add formatting target
function(add_format_target)
    if(NOT CLANG_FORMAT_EXECUTABLE)
        message(STATUS "clang-format not found, code formatting target will not be available")
        return()
    endif()

    # Get all source files
    file(GLOB_RECURSE ALL_SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/src/*.cpp
        ${CMAKE_SOURCE_DIR}/src/*.hpp
        ${CMAKE_SOURCE_DIR}/src/*.h
        ${CMAKE_SOURCE_DIR}/tests/*.cpp
        ${CMAKE_SOURCE_DIR}/tests/*.hpp
        ${CMAKE_SOURCE_DIR}/tests/*.h
    )

    # Create format target
    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -i --style=file ${ALL_SOURCE_FILES}
        COMMENT "Formatting source code with clang-format"
        VERBATIM
    )

    # Make the main target depend on format
    add_dependencies(metada format)
endfunction() 