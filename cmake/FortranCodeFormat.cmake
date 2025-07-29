# FortranCodeFormat.cmake
#
# This module provides automatic Fortran code formatting using fprettify.
#
# Usage:
#   AddFortranFormatTarget(target_name src_dir [EXCLUDE_DIRS dir1;dir2;...])
#
# Parameters:
#   target_name: The CMake target name that will depend on formatting
#   src_dir: The source directory to search for Fortran files
#   EXCLUDE_DIRS: Optional list of directories to exclude from formatting
#
# Examples:
#   # Format all Fortran files in src_dir
#   AddFortranFormatTarget(my_target ${CMAKE_CURRENT_SOURCE_DIR})
#
#   # Exclude specific directory from formatting
#   AddFortranFormatTarget(my_target ${CMAKE_CURRENT_SOURCE_DIR} 
#       EXCLUDE_DIRS "fortran/src-region")
#
#   # Exclude multiple directories from formatting
#   AddFortranFormatTarget(my_target ${CMAKE_CURRENT_SOURCE_DIR} 
#       EXCLUDE_DIRS "fortran/src-region;fortran/src-globe;legacy")
#
# Supported Fortran file extensions:
#   .f90, .f95, .f03, .f08, .f, .for, .ftn, .fpp
#   .F90, .F95, .F03, .F08, .F, .FOR, .FTN, .FPP
#
# The function will:
#   1. Recursively search for all Fortran files in src_dir
#   2. Filter out files in excluded directories (if specified)
#   3. Create a custom target that formats files when they are modified
#   4. Make the target depend on the formatting target

# Function to add Fortran formatting target
function(AddFortranFormatTarget target_name src_dir)
    # Parse optional arguments for exclude directories
    cmake_parse_arguments(ARG
        ""
        ""
        "EXCLUDE_DIRS"
        ${ARGN}
    )
    if(NOT Fprettify_EXECUTABLE)
        message(WARNING "fprettify not found, Fortran code formatting target will not be available")
        return()
    endif()

    # Get all Fortran source files
    file(GLOB_RECURSE ALL_FORTRAN_FILES
        ${src_dir}/*.f90
        ${src_dir}/*.f95
        ${src_dir}/*.f03
        ${src_dir}/*.f08
        ${src_dir}/*.f
        ${src_dir}/*.for
        ${src_dir}/*.ftn
        ${src_dir}/*.fpp
        ${src_dir}/*.F90
        ${src_dir}/*.F95
        ${src_dir}/*.F03
        ${src_dir}/*.F08
        ${src_dir}/*.F
        ${src_dir}/*.FOR
        ${src_dir}/*.FTN
        ${src_dir}/*.FPP
    )

    # Filter out excluded directories
    if(ARG_EXCLUDE_DIRS)
        set(FILTERED_FORTRAN_FILES)
        set(EXCLUDED_FILES)
        foreach(file ${ALL_FORTRAN_FILES})
            set(should_exclude FALSE)
            foreach(exclude_dir ${ARG_EXCLUDE_DIRS})
                string(FIND "${file}" "${exclude_dir}" pos)
                if(NOT pos EQUAL -1)
                    set(should_exclude TRUE)
                    list(APPEND EXCLUDED_FILES ${file})
                    break()
                endif()
            endforeach()
            if(NOT should_exclude)
                list(APPEND FILTERED_FORTRAN_FILES ${file})
            endif()
        endforeach()
        set(ALL_FORTRAN_FILES ${FILTERED_FORTRAN_FILES})
        
        # Print excluded files for debugging
        if(EXCLUDED_FILES)
            message(STATUS "Fortran formatting: Excluded ${EXCLUDED_FILES} files from formatting")
        endif()
        
        # Print summary of files to be formatted
        list(LENGTH ALL_FORTRAN_FILES FORMAT_COUNT)
        list(LENGTH EXCLUDED_FILES EXCLUDED_COUNT)
        message(STATUS "Fortran formatting: Will format ${FORMAT_COUNT} files, excluded ${EXCLUDED_COUNT} files")
    endif()

    # Skip if no files to format
    list(LENGTH ALL_FORTRAN_FILES FORTRAN_COUNT)
    if(FORTRAN_COUNT EQUAL 0)
        message(WARNING "No Fortran source files found in ${src_dir} for formatting")
        return()
    endif()

    # Create a timestamp file to track last format
    set(FORTRAN_TIMESTAMP_FILE "${CMAKE_BINARY_DIR}/${target_name}_fortran_format_timestamp")

    # Create format target with timestamp check
    set(fortran_format_target_name "format_fortran_${target_name}")
    
    # Handle both executable and Python module forms of fprettify
    if(Fprettify_EXECUTABLE MATCHES ".*python.*-m.*fprettify")
        # Python module form
        add_custom_command(
            OUTPUT ${FORTRAN_TIMESTAMP_FILE}
            COMMAND ${Fprettify_EXECUTABLE} -r ${ALL_FORTRAN_FILES}
            COMMAND ${CMAKE_COMMAND} -E touch ${FORTRAN_TIMESTAMP_FILE}
            DEPENDS ${ALL_FORTRAN_FILES}
            COMMENT "Formatting Fortran source code with fprettify (Python module) in ${src_dir}"
            VERBATIM
        )
    else()
        # Executable form
        add_custom_command(
            OUTPUT ${FORTRAN_TIMESTAMP_FILE}
            COMMAND ${Fprettify_EXECUTABLE} -r ${ALL_FORTRAN_FILES}
            COMMAND ${CMAKE_COMMAND} -E touch ${FORTRAN_TIMESTAMP_FILE}
            DEPENDS ${ALL_FORTRAN_FILES}
            COMMENT "Formatting Fortran source code with fprettify in ${src_dir}"
            VERBATIM
        )
    endif()

    # Create a custom target that depends on the timestamp file.
    # This target will only trigger the formatting command when source files are modified
    # because the timestamp file's dependencies are set to the source files.
    add_custom_target(${fortran_format_target_name}
        DEPENDS ${FORTRAN_TIMESTAMP_FILE}
    )

    # Make the target depend on format
    add_dependencies(${target_name} ${fortran_format_target_name})
endfunction() 