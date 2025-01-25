# Documentation options
option(BUILD_DOXYGEN "Build Doxygen API documentation" ON)
option(BUILD_SPHINX "Build Sphinx documentation" ON)

# Create main docs target
add_custom_target(docs)

# Add subdirectories for each documentation type
add_subdirectory(doxygen)
add_subdirectory(sphinx) 