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
        # Add dependencies on Python test targets that generate API documentation
        # This ensures documentation is generated before Sphinx build
        foreach(test_target ${PYTHON_TEST_TARGETS})
            if(TARGET ${test_target})
                add_dependencies(sphinx ${test_target})
            endif()
        endforeach()
    endif()
endif() 