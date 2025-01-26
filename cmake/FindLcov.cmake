# Find lcov and genhtml executables
#
# This module defines the following variables:
#  Lcov_EXECUTABLE     - Path to lcov executable
#  Lcov_GENHTML       - Path to genhtml executable
#  Lcov_VERSION       - Version of lcov found
#  Lcov_FOUND         - True if lcov was found

# Look for lcov executable
find_program(Lcov_EXECUTABLE
    NAMES 
        lcov
    PATHS
        # Unix-like systems
        /usr/bin
        /usr/local/bin
        /opt/local/bin
        # Windows - MSYS2/Cygwin paths
        "C:/msys64/mingw64/bin"
        "C:/cygwin64/bin"
        # macOS - Homebrew
        /usr/local/opt/lcov/bin
        /opt/homebrew/bin
    DOC "Path to lcov executable"
)

# Look for genhtml executable
find_program(Lcov_GENHTML
    NAMES 
        genhtml
    PATHS
        # Unix-like systems
        /usr/bin
        /usr/local/bin
        /opt/local/bin
        # Windows - MSYS2/Cygwin paths
        "C:/msys64/mingw64/bin"
        "C:/cygwin64/bin"
        # macOS - Homebrew
        /usr/local/opt/lcov/bin
        /opt/homebrew/bin
    DOC "Path to genhtml executable"
)

if(Lcov_EXECUTABLE)
    # Get version information with error logging
    if(WIN32)
        execute_process(
            COMMAND perl "${Lcov_EXECUTABLE}" --version
            OUTPUT_VARIABLE Lcov_VERSION_OUTPUT
            ERROR_VARIABLE Lcov_VERSION_ERROR
            RESULT_VARIABLE Lcov_VERSION_RESULT
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    else()
        execute_process(
            COMMAND "${Lcov_EXECUTABLE}" --version
            OUTPUT_VARIABLE Lcov_VERSION_OUTPUT
            ERROR_VARIABLE Lcov_VERSION_ERROR
            RESULT_VARIABLE Lcov_VERSION_RESULT
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()

    if(NOT Lcov_VERSION_RESULT EQUAL 0)
        message(WARNING "Failed to get lcov version: ${Lcov_VERSION_ERROR}")
    else()
        if(Lcov_VERSION_OUTPUT MATCHES "LCOV version ([0-9]+\\.[0-9]+)")
            set(Lcov_VERSION "${CMAKE_MATCH_1}")
        elseif(Lcov_VERSION_OUTPUT MATCHES "LCOV version ([^\\s]+)")
            # Handle non-standard version formats like "rcinfo-cache-5767-g08b7c8fcc"
            set(Lcov_VERSION "${CMAKE_MATCH_1}")
        else()
            message(WARNING "Failed to parse lcov version from: ${Lcov_VERSION_OUTPUT}")
        endif()
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lcov
    REQUIRED_VARS 
        Lcov_EXECUTABLE
        Lcov_GENHTML
    VERSION_VAR Lcov_VERSION
)

# Mark advanced variables
mark_as_advanced(
    Lcov_EXECUTABLE
    Lcov_GENHTML
    Lcov_VERSION
) 