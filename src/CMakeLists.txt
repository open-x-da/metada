# Add subdirectories
add_subdirectory(core)
add_subdirectory(backends)
add_subdirectory(applications)

# Create the main library
add_library(metada_core INTERFACE)
add_library(metada::core ALIAS metada_core)

target_link_libraries(metada_core
    INTERFACE
        metada::core::tools::logging
        metada::backends_logging
)

# Set include directories
target_include_directories(metada_core
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)