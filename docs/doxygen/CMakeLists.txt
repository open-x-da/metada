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
        # Download and configure doxygen-awesome-css theme
        # This provides a modern, responsive theme for the Doxygen documentation
        include(FetchContent)
        FetchContent_Declare(doxygen-awesome-css
            GIT_REPOSITORY
                https://github.com/jothepro/doxygen-awesome-css.git
            GIT_TAG
                v2.3.4
        )
        FetchContent_MakeAvailable(doxygen-awesome-css)
        # Configure Doxygen
        set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/docs/doxygen)
        set(DOXYGEN_GENERATE_HTML         YES)
        set(DOXYGEN_GENERATE_TREEVIEW     YES)
        set(DOXYGEN_GENERATE_LATEX         NO)
        set(DOXYGEN_HTML_EXTRA_STYLESHEET
            ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome.css)
        
        if(Doxygen_dot_FOUND)
            set(DOXYGEN_HAVE_DOT YES)
            set(DOXYGEN_UML_LOOK YES)
            set(DOXYGEN_CALL_GRAPH YES)
            set(DOXYGEN_CALLER_GRAPH YES)
            set(DOXYGEN_DOT_IMAGE_FORMAT svg)
            set(DOXYGEN_DOT_TRANSPARENT YES)
        endif()

        doxygen_add_docs(doxygen
            ${PROJECT_SOURCE_DIR}/src
            ${PROJECT_SOURCE_DIR}/tests
            COMMENT "Generating documentation with Doxygen"
        )
        add_dependencies(docs doxygen)
    endif()
endif() 