# Add subdirectories
add_subdirectory(core)
add_subdirectory(backends)
add_subdirectory(applications)

# Create the main library
add_library(metada_core INTERFACE)
target_link_libraries(metada_core
    INTERFACE
        metada_core_logging
        metada_backends_logging
)

# Set include directories
target_include_directories(metada_core
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Install rules for libraries only
install(TARGETS metada_core
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
) 