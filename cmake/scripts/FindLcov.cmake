# Find lcov and genhtml executables
#
# This module defines the following variables:
#  Lcov_EXECUTABLE     - Path to lcov executable
#  Lcov_GENHTML       - Path to genhtml executable
#  Lcov_VERSION       - Version of lcov found
#  Lcov_FOUND         - True if lcov was found
#
# The module respects the following variables if set:
#  LCOV_ROOT_DIR      - Custom path to look for lcov

# Allow users to specify custom search paths
if(DEFINED ENV{LCOV_ROOT_DIR})
    set(LCOV_ROOT_DIR $ENV{LCOV_ROOT_DIR})
endif()

# Set up executable names
set(LCOV_NAMES lcov)
set(GENHTML_NAMES genhtml)

# Platform-specific search paths
if(WIN32)
    # Windows-specific paths
    set(PLATFORM_PATHS
        # Windows - MSYS2/Cygwin paths
        "C:/msys64/mingw64/bin"
        "C:/msys64/usr/bin"
        "C:/cygwin64/bin"
        "C:/cygwin/bin"
    )
elseif(APPLE)
    # macOS-specific paths
    set(PLATFORM_PATHS
        # macOS - Homebrew
        /usr/local/opt/lcov/bin
        /opt/homebrew/bin
        /opt/homebrew/opt/lcov/bin
        # macOS standard paths
        /usr/local/bin
        /usr/bin
        /opt/local/bin
    )
else()
    # Linux/Unix specific paths
    set(PLATFORM_PATHS
        # Linux standard system directories
        /usr/bin
        /usr/local/bin
        /bin
        /opt/bin
        /opt/local/bin
        # Package manager specific directories
        /snap/bin
        # Spack installations (generic patterns)
        $ENV{SPACK_ROOT}/var/spack/environments/*/view/bin
        $ENV{HOME}/spack/var/spack/environments/*/view/bin
        /*/spack/var/spack/environments/*/view/bin
        # Common custom Spack installation locations
        /opt/spack/var/spack/environments/*/view/bin
        /mnt/*/spack/var/spack/environments/*/view/bin
    )
    
    # Special handling for WSL - explicitly avoid Windows paths
    if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux" AND EXISTS "/proc/version")
        file(READ "/proc/version" PROC_VERSION)
        if(PROC_VERSION MATCHES "Microsoft" OR PROC_VERSION MATCHES "WSL")
            message(STATUS "WSL detected - prioritizing Linux paths for lcov")
        endif()
    endif()
endif()

# Look for lcov executable
find_program(Lcov_EXECUTABLE
    NAMES ${LCOV_NAMES}
    PATHS
        # User-specified directories first
        ${LCOV_ROOT_DIR}
        # Platform-specific paths
        ${PLATFORM_PATHS}
    DOC "Path to lcov executable"
)

# Look for genhtml executable
find_program(Lcov_GENHTML
    NAMES ${GENHTML_NAMES}
    PATHS
        # User-specified directories first
        ${LCOV_ROOT_DIR}
        # Platform-specific paths
        ${PLATFORM_PATHS}
    DOC "Path to genhtml executable"
)

# If not found in specific locations, check using the PATH environment variable
if(NOT Lcov_EXECUTABLE)
    find_program(Lcov_EXECUTABLE NAMES ${LCOV_NAMES})
endif()

if(NOT Lcov_GENHTML)
    find_program(Lcov_GENHTML NAMES ${GENHTML_NAMES})
endif()

# Special handling for Spack environments
if((NOT Lcov_EXECUTABLE OR NOT Lcov_GENHTML) AND DEFINED ENV{SPACK_ENV})
    if(NOT Lcov_EXECUTABLE)
        find_program(Lcov_EXECUTABLE
            NAMES ${LCOV_NAMES}
            PATHS $ENV{SPACK_ENV}/.spack-env/view/bin
            NO_DEFAULT_PATH
        )
    endif()
    
    if(NOT Lcov_GENHTML)
        find_program(Lcov_GENHTML
            NAMES ${GENHTML_NAMES}
            PATHS $ENV{SPACK_ENV}/.spack-env/view/bin
            NO_DEFAULT_PATH
        )
    endif()
endif()

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

# Provide a clear status message about coverage availability
if(Lcov_FOUND)
    message(STATUS "Coverage testing: ENABLED (lcov ${Lcov_VERSION} found)")
else()
    message(STATUS "Coverage testing: DISABLED (lcov/genhtml not found)")
endif()

# Mark advanced variables
mark_as_advanced(
    Lcov_EXECUTABLE
    Lcov_GENHTML
    Lcov_VERSION
) 