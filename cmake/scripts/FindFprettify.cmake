# FindFprettify.cmake
# Locates the fprettify Python package executable
#
# This module defines the following variables:
#   Fprettify_FOUND - True if fprettify executable is found
#   Fprettify_EXECUTABLE - Path to the fprettify executable
#   Fprettify_VERSION - Version of fprettify if available

# Try to find fprettify executable
find_program(Fprettify_EXECUTABLE
    NAMES fprettify
    PATHS
        ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin
        ${CMAKE_CURRENT_SOURCE_DIR}/.venv/bin
        $ENV{VIRTUAL_ENV}/bin
        $ENV{CONDA_PREFIX}/bin
        $ENV{CONDA_PREFIX}/Scripts
        $ENV{PYTHONPATH}
        $ENV{PATH}
    DOC "Path to fprettify executable"
)

# If not found in PATH, try to find it via Python
if(NOT Fprettify_EXECUTABLE)
    find_package(Python3 QUIET)
    if(Python3_FOUND)
        # Try to find fprettify via Python
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -m fprettify --version
            OUTPUT_VARIABLE FPRETTIFY_VERSION_OUTPUT
            ERROR_VARIABLE FPRETTIFY_VERSION_ERROR
            RESULT_VARIABLE FPRETTIFY_VERSION_RESULT
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        
        if(FPRETTIFY_VERSION_RESULT EQUAL 0)
            # fprettify is available via Python module
            set(Fprettify_EXECUTABLE "${Python3_EXECUTABLE} -m fprettify")
            set(Fprettify_VERSION ${FPRETTIFY_VERSION_OUTPUT})
        endif()
    endif()
endif()

# Handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Fprettify
    REQUIRED_VARS Fprettify_EXECUTABLE
    VERSION_VAR Fprettify_VERSION
)

# Mark as advanced
mark_as_advanced(Fprettify_EXECUTABLE Fprettify_VERSION) 