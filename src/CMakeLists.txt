# Add subdirectories
add_subdirectory(framework)
add_subdirectory(backends)

# Create the main library
add_library(metada_core INTERFACE)
add_library(metada::core ALIAS metada_core)

target_link_libraries(metada_core
    INTERFACE
        metada::framework::tools::logger
        metada::backends::tools::logger
)

# Set include directories
target_include_directories(metada_core
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)