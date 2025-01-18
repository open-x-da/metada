# Find documentation tools
find_package(Doxygen OPTIONAL_COMPONENTS dot)
find_package(Sphinx OPTIONAL_COMPONENTS build)

function(configure_documentation)
    # Create main docs target first
    add_custom_target(docs)

    if(DOXYGEN_FOUND)
        # Configure Doxygen
        set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/docs/doxygen)
        set(DOXYGEN_GENERATE_HTML YES)
        set(DOXYGEN_GENERATE_LATEX NO)
        
        # Optional dot tool configuration
        if(DOXYGEN_DOT_FOUND)
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

    if(SPHINX_FOUND)
        # Add Sphinx target
        add_custom_target(sphinx
            ${SPHINX_EXECUTABLE} -b html
            ${CMAKE_CURRENT_SOURCE_DIR}/docs/sphinx
            ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx
            COMMENT "Generating documentation with Sphinx"
            VERBATIM
        )
        add_dependencies(sphinx lorenz95_python_test)
        add_dependencies(docs sphinx)
    endif()
endfunction() 