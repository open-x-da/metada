# Documentation options
option(BUILD_DOXYGEN "Build Doxygen API documentation" ON)
option(BUILD_SPHINX "Build Sphinx documentation" ON)

# Add subdirectories for each documentation type to check availability
add_subdirectory(doxygen)
add_subdirectory(sphinx)

# Create main docs target only if at least one documentation system is available
if((BUILD_DOXYGEN AND Doxygen_FOUND) OR (BUILD_SPHINX AND Sphinx_FOUND))
    add_custom_target(docs)
    
    # Add dependencies if documentation systems are available
    if(BUILD_DOXYGEN AND Doxygen_FOUND)
        add_dependencies(docs doxygen)
    endif()
    
    if(BUILD_SPHINX AND Sphinx_FOUND)
        add_dependencies(docs sphinx)
    endif()
endif() 