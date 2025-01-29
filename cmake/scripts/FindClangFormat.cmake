# Find clang-format executable
#
# This module defines the following variables:
#  ClangFormat_EXECUTABLE - Path to clang-format executable
#  ClangFormat_VERSION   - Version of clang-format found
#  ClangFormat_FOUND     - True if clang-format was found

# Look for clang-format in standard paths and known install locations
find_program(ClangFormat_EXECUTABLE
    NAMES 
        # Generic names
        clang-format
        # Versioned names (newer to older)
        clang-format-17
        clang-format-16 
        clang-format-15
        clang-format-14
        # Windows specific names
        clang-format.exe
    PATHS
        # Unix-like systems
        /usr/bin
        /usr/local/bin
        /opt/local/bin
        /opt/bin
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
        # macOS - Homebrew
        /usr/local/opt/llvm/bin
        # macOS - Command Line Tools
        /Library/Developer/CommandLineTools/usr/bin
    DOC "Path to clang-format executable"
)

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