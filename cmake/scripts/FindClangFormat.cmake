# Find clang-format executable
#
# This module defines the following variables:
#  ClangFormat_EXECUTABLE - Path to clang-format executable
#  ClangFormat_VERSION   - Version of clang-format found
#  ClangFormat_FOUND     - True if clang-format was found
#
# The module respects the following variables if set:
#  CLANG_FORMAT_ROOT_DIR - Custom path to look for clang-format

# Allow users to specify custom search paths
if(DEFINED ENV{CLANG_FORMAT_ROOT_DIR})
    set(CLANG_FORMAT_ROOT_DIR $ENV{CLANG_FORMAT_ROOT_DIR})
endif()

# Set up platform-specific executable names and default locations
set(CLANG_FORMAT_NAMES clang-format clang-format-17 clang-format-16 clang-format-15 clang-format-14 clang-format-13 clang-format-12 clang-format-11 clang-format-10)

# Platform-specific search paths
if(WIN32)
    # Windows-specific executable names and paths
    list(APPEND CLANG_FORMAT_NAMES clang-format.exe)
    set(PLATFORM_PATHS
        # Windows - LLVM standard install paths
        "C:/Program Files/LLVM/bin"
        "C:/Program Files (x86)/LLVM/bin"
        # Windows - MSYS2 paths
        "C:/msys64/mingw64/bin"
        "C:/msys64/usr/bin"
        # Windows - Visual Studio bundled clang
        "$ENV{ProgramFiles\(x86\)}/Microsoft Visual Studio/*/Professional/VC/Tools/Llvm/bin"
        "$ENV{ProgramFiles\(x86\)}/Microsoft Visual Studio/*/Enterprise/VC/Tools/Llvm/bin"
        "$ENV{ProgramFiles\(x86\)}/Microsoft Visual Studio/*/Community/VC/Tools/Llvm/bin"
    )
elseif(APPLE)
    # macOS-specific paths
    set(PLATFORM_PATHS
        # macOS - Homebrew
        /usr/local/opt/llvm/bin
        # macOS - Command Line Tools
        /Library/Developer/CommandLineTools/usr/bin
        # macOS standard paths
        /usr/local/bin
        /usr/bin
        /opt/homebrew/bin
        /opt/homebrew/opt/llvm/bin
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
        /usr/lib/llvm/*/bin
        /usr/lib64/llvm/bin
        # Spack installations (generic patterns)
        $ENV{SPACK_ROOT}/var/spack/environments/*/view/bin
        $ENV{HOME}/spack/var/spack/environments/*/view/bin
        /*/spack/var/spack/environments/*/view/bin
    )
    
    # Special handling for WSL - explicitly avoid Windows paths
    if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux" AND EXISTS "/proc/version")
        file(READ "/proc/version" PROC_VERSION)
        if(PROC_VERSION MATCHES "Microsoft" OR PROC_VERSION MATCHES "WSL")
            message(STATUS "WSL detected - prioritizing Linux paths for clang-format")
        endif()
    endif()
endif()

# Look for clang-format in platform-specific paths
find_program(ClangFormat_EXECUTABLE
    NAMES ${CLANG_FORMAT_NAMES}
    PATHS
        # User-specified directories first (platform independent)
        ${CLANG_FORMAT_ROOT_DIR}
        # Platform-specific paths
        ${PLATFORM_PATHS}
    DOC "Path to clang-format executable"
)

# If not found in platform-specific locations, check using the PATH environment variable
if(NOT ClangFormat_EXECUTABLE)
    find_program(ClangFormat_EXECUTABLE
        NAMES ${CLANG_FORMAT_NAMES}
    )
endif()

# Special handling for Spack environments (platform independent)
if(NOT ClangFormat_EXECUTABLE AND DEFINED ENV{SPACK_ENV})
    find_program(ClangFormat_EXECUTABLE
        NAMES ${CLANG_FORMAT_NAMES}
        PATHS $ENV{SPACK_ENV}/.spack-env/view/bin
        NO_DEFAULT_PATH
    )
endif()

if(ClangFormat_EXECUTABLE)
    # Get version information
    execute_process(
        COMMAND ${ClangFormat_EXECUTABLE} --version
        OUTPUT_VARIABLE ClangFormat_VERSION_OUTPUT
        ERROR_VARIABLE ClangFormat_VERSION_ERROR
        RESULT_VARIABLE ClangFormat_VERSION_RESULT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(NOT ClangFormat_VERSION_RESULT EQUAL 0)
        message(WARNING "Failed to get clang-format version: ${ClangFormat_VERSION_ERROR}")
    else()
        # Extract version from output (format: "clang-format version X.Y.Z")
        if(ClangFormat_VERSION_OUTPUT MATCHES "version ([0-9]+\\.[0-9]+\\.[0-9]+)")
            set(ClangFormat_VERSION "${CMAKE_MATCH_1}")
        elseif(ClangFormat_VERSION_OUTPUT MATCHES "version ([0-9]+)")
            set(ClangFormat_VERSION "${CMAKE_MATCH_1}.0.0")
        endif()
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ClangFormat
    REQUIRED_VARS ClangFormat_EXECUTABLE
    VERSION_VAR ClangFormat_VERSION
)

# Mark advanced variables
mark_as_advanced(
    ClangFormat_EXECUTABLE
    ClangFormat_VERSION
) 