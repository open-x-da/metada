# Options for documentation generation
option(BUILD_DOXYGEN "Build Doxygen API documentation. Requires Doxygen to be installed." ON)
option(BUILD_SPHINX "Build Sphinx documentation. Requires Sphinx to be installed." ON)

# Function to configure documentation support
function(configure_documentation_support)
    # Create main docs target first
    add_custom_target(docs)

    # Configure Doxygen
    if(BUILD_DOXYGEN)
        metada_find_package(Doxygen 
            OPTIONAL 
            CONDITION BUILD_DOXYGEN 
            QUIET 
            COMPONENTS dot
        )

        if(NOT Doxygen_FOUND)
            message(STATUS "Doxygen documentation generation disabled: Doxygen executable not found in PATH")
            set(BUILD_DOXYGEN OFF CACHE BOOL "Build Doxygen API documentation" FORCE)
        elseif(BUILD_DOXYGEN)
            # Configure Doxygen
            set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/docs/doxygen)
            set(DOXYGEN_GENERATE_HTML YES)
            set(DOXYGEN_GENERATE_LATEX NO)
            
            # Optional dot tool configuration
            if(Doxygen_DOT_FOUND)
                set(DOXYGEN_HAVE_DOT YES)
                set(DOXYGEN_UML_LOOK YES)
                set(DOXYGEN_CALL_GRAPH YES)
                set(DOXYGEN_CALLER_GRAPH YES)
            endif()

            # Add Doxygen target
            doxygen_add_docs(doxygen
                ${PROJECT_SOURCE_DIR}/src
                ${PROJECT_SOURCE_DIR}/tests
                COMMENT "Generating documentation with Doxygen"
            )
            add_dependencies(docs doxygen)
        endif()
    endif()

    # Configure Sphinx
    if(BUILD_SPHINX)
        metada_find_package(Sphinx
            OPTIONAL
            CONDITION BUILD_SPHINX
            QUIET
            COMPONENTS build
        )

        if(NOT Sphinx_FOUND)
            message(STATUS "Sphinx documentation generation disabled: Sphinx executable not found in PATH")
            set(BUILD_SPHINX OFF CACHE BOOL "Build Sphinx documentation" FORCE)
        elseif(BUILD_SPHINX)
            # Add Sphinx target
            add_custom_target(sphinx
                ${SPHINX_EXECUTABLE} -b html
                ${CMAKE_CURRENT_SOURCE_DIR}/docs/sphinx
                ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx
                COMMENT "Generating documentation with Sphinx"
                VERBATIM
            )
            # We need to ensure Python tests are run before building Sphinx documentation
            # because the tests generate Python API documentation that Sphinx needs
            if(TARGET lorenz95_python_test)
                add_dependencies(sphinx lorenz95_python_test)
            endif()
            add_dependencies(docs sphinx)
        endif()
    endif()
endfunction()

# Configure documentation support immediately when this module is included
configure_documentation_support() 