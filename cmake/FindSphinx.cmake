find_program(SPHINX_EXECUTABLE
    NAMES sphinx-build sphinx-build.exe
    DOC "Path to sphinx-build executable"
)

if(SPHINX_EXECUTABLE)
    # Get version information
    execute_process(
        COMMAND ${SPHINX_EXECUTABLE} --version
        OUTPUT_VARIABLE SPHINX_VERSION_OUTPUT
        ERROR_VARIABLE SPHINX_VERSION_ERROR
        RESULT_VARIABLE SPHINX_VERSION_RESULT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(SPHINX_VERSION_RESULT EQUAL 0)
        # Extract version from output (format: "sphinx-build X.Y.Z")
        string(REGEX REPLACE "^sphinx-build ([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1"
            SPHINX_VERSION "${SPHINX_VERSION_OUTPUT}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sphinx
    REQUIRED_VARS SPHINX_EXECUTABLE
    VERSION_VAR SPHINX_VERSION
)

mark_as_advanced(SPHINX_EXECUTABLE SPHINX_VERSION) 