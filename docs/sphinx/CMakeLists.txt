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
        add_custom_target(sphinx
            ${SPHINX_EXECUTABLE} -b html
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating documentation with Sphinx"
            VERBATIM
        )
        if(TARGET lorenz95_python_test)
            add_dependencies(sphinx lorenz95_python_test)
        endif()
        add_dependencies(docs sphinx)
    endif()
endif() 