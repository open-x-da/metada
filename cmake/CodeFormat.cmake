# Function to add formatting target
function(add_format_target target_name src_dir)
    if(NOT ClangFormat_EXECUTABLE)
        message(WARNING "clang-format not found, code formatting target will not be available")
        return()
    endif()

    # Get all source files
    file(GLOB_RECURSE ALL_SOURCE_FILES
        ${src_dir}/*.cpp
        ${src_dir}/*.cc
        ${src_dir}/*.hpp
        ${src_dir}/*.h
    )

    # Skip if no files to format
    list(LENGTH ALL_SOURCE_FILES SOURCE_COUNT)
    if(SOURCE_COUNT EQUAL 0)
        message(WARNING "No source files found in ${src_dir} for formatting")
        return()
    endif()

    # Create a timestamp file to track last format
    set(TIMESTAMP_FILE "${CMAKE_BINARY_DIR}/${target_name}_format_timestamp")

    # Create format target with timestamp check
    set(format_target_name "format_${target_name}")
    add_custom_command(
        OUTPUT ${TIMESTAMP_FILE}
        COMMAND ${ClangFormat_EXECUTABLE} -style=file -i ${ALL_SOURCE_FILES}
        COMMAND ${CMAKE_COMMAND} -E touch ${TIMESTAMP_FILE}
        DEPENDS ${ALL_SOURCE_FILES}
        COMMENT "Formatting source code with clang-format in ${src_dir}"
        VERBATIM
    )

    # Create a custom target that depends on the timestamp file.
    # This target will only trigger the formatting command when source files are modified
    # because the timestamp file's dependencies are set to the source files.
    add_custom_target(${format_target_name}
        DEPENDS ${TIMESTAMP_FILE}
    )

    # Make the target depend on format
    add_dependencies(${target_name} ${format_target_name})
endfunction() 