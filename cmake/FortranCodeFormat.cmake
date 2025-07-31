# FortranCodeFormat.cmake
#
# This module provides automatic Fortran code formatting using fprettify.
#
# Usage:
#   AddFortranFormatTarget(target_name src_dir [EXCLUDE_DIRS dir1;dir2;...])
#   AddFortranFormatTargetSelective(target_name [FILES file1;file2;...] [DIRS dir1;dir2;...])
#
# Parameters:
#   target_name: The CMake target name that will depend on formatting
#   src_dir: The source directory to search for Fortran files
#   EXCLUDE_DIRS: Optional list of directories to exclude from formatting
#   FILES: Specific files to format (for AddFortranFormatTargetSelective)
#   DIRS: Specific directories to format (for AddFortranFormatTargetSelective)
#
# Examples:
#   # Format all Fortran files in src_dir
#   AddFortranFormatTarget(my_target ${CMAKE_CURRENT_SOURCE_DIR})
#
#   # Exclude specific directory from formatting
#   AddFortranFormatTarget(my_target ${CMAKE_CURRENT_SOURCE_DIR} 
#       EXCLUDE_DIRS "fortran/src-region")
#
#   # Format only specific files (FASTEST)
#   AddFortranFormatTargetSelective(my_target FILES "fortran/macom_logger.f90;fortran/macom_fortran_wrapper.f90")
#
#   # Format only specific directories
#   AddFortranFormatTargetSelective(my_target DIRS "fortran/wrapper;fortran/utils")
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
            COMMAND ${Fprettify_EXECUTABLE} -c ${CMAKE_SOURCE_DIR}/.fprettify.rc ${ALL_FORTRAN_FILES}
            COMMAND ${CMAKE_COMMAND} -E touch ${FORTRAN_TIMESTAMP_FILE}
            DEPENDS ${ALL_FORTRAN_FILES}
            COMMENT "Formatting Fortran source code with fprettify (Python module) in ${src_dir}"
            VERBATIM
        )
    else()
        # Executable form
        add_custom_command(
            OUTPUT ${FORTRAN_TIMESTAMP_FILE}
            COMMAND ${Fprettify_EXECUTABLE} -c ${CMAKE_SOURCE_DIR}/.fprettify.rc ${ALL_FORTRAN_FILES}
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

# Function to add selective Fortran formatting target (FASTER)
function(AddFortranFormatTargetSelective target_name)
    # Parse arguments for specific files and directories
    cmake_parse_arguments(ARG
        ""
        ""
        "FILES;DIRS"
        ${ARGN}
    )
    
    if(NOT Fprettify_EXECUTABLE)
        message(WARNING "fprettify not found, Fortran code formatting target will not be available")
        return()
    endif()

    set(FORTRAN_FILES_TO_FORMAT)

    # Add specific files if provided
    if(ARG_FILES)
        foreach(file ${ARG_FILES})
            # Handle both relative and absolute paths
            if(IS_ABSOLUTE ${file})
                # Absolute path provided
                if(EXISTS "${file}")
                    list(APPEND FORTRAN_FILES_TO_FORMAT "${file}")
                else()
                    message(WARNING "Fortran file not found: ${file}")
                endif()
            else()
                # Relative path - try both CMAKE_CURRENT_SOURCE_DIR and CMAKE_SOURCE_DIR
                if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
                    list(APPEND FORTRAN_FILES_TO_FORMAT "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
                elseif(EXISTS "${CMAKE_SOURCE_DIR}/${file}")
                    list(APPEND FORTRAN_FILES_TO_FORMAT "${CMAKE_SOURCE_DIR}/${file}")
                else()
                    message(WARNING "Fortran file not found: ${file} (tried ${CMAKE_CURRENT_SOURCE_DIR}/${file} and ${CMAKE_SOURCE_DIR}/${file})")
                endif()
            endif()
        endforeach()
    endif()

    # Add files from specific directories if provided
    if(ARG_DIRS)
        foreach(dir ${ARG_DIRS})
            file(GLOB_RECURSE DIR_FORTRAN_FILES
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.f90
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.f95
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.f03
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.f08
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.f
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.for
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.ftn
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.fpp
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.F95
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.F03
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.F08
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.F
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.FOR
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.FTN
                ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.FPP
            )
            list(APPEND FORTRAN_FILES_TO_FORMAT ${DIR_FORTRAN_FILES})
        endforeach()
    endif()

    # Skip if no files to format
    list(LENGTH FORTRAN_FILES_TO_FORMAT FORTRAN_COUNT)
    if(FORTRAN_COUNT EQUAL 0)
        message(WARNING "No Fortran source files specified for formatting")
        return()
    endif()

    # Remove duplicates
    list(REMOVE_DUPLICATES FORTRAN_FILES_TO_FORMAT)

    # Print summary
    list(LENGTH FORTRAN_FILES_TO_FORMAT FORMAT_COUNT)
    message(STATUS "Fortran formatting: Will format ${FORMAT_COUNT} specific files")

    # Create a timestamp file to track last format
    set(FORTRAN_TIMESTAMP_FILE "${CMAKE_BINARY_DIR}/${target_name}_fortran_format_selective_timestamp")

    # Create format target with timestamp check
    set(fortran_format_target_name "format_fortran_selective_${target_name}")
    
    # Handle both executable and Python module forms of fprettify
    if(Fprettify_EXECUTABLE MATCHES ".*python.*-m.*fprettify")
        # Python module form
        add_custom_command(
            OUTPUT ${FORTRAN_TIMESTAMP_FILE}
            COMMAND ${Fprettify_EXECUTABLE} -c ${CMAKE_SOURCE_DIR}/.fprettify.rc ${FORTRAN_FILES_TO_FORMAT}
            COMMAND ${CMAKE_COMMAND} -E touch ${FORTRAN_TIMESTAMP_FILE}
            DEPENDS ${FORTRAN_FILES_TO_FORMAT}
            COMMENT "Formatting specific Fortran files with fprettify (Python module)"
            VERBATIM
        )
    else()
        # Executable form
        add_custom_command(
            OUTPUT ${FORTRAN_TIMESTAMP_FILE}
            COMMAND ${Fprettify_EXECUTABLE} -c ${CMAKE_SOURCE_DIR}/.fprettify.rc ${FORTRAN_FILES_TO_FORMAT}
            COMMAND ${CMAKE_COMMAND} -E touch ${FORTRAN_TIMESTAMP_FILE}
            DEPENDS ${FORTRAN_FILES_TO_FORMAT}
            COMMENT "Formatting specific Fortran files with fprettify"
            VERBATIM
        )
    endif()

    # Create a custom target that depends on the timestamp file
    add_custom_target(${fortran_format_target_name}
        DEPENDS ${FORTRAN_TIMESTAMP_FILE}
    )

    # Make the target depend on format
    add_dependencies(${target_name} ${fortran_format_target_name})
endfunction() 