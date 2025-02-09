# Create core interface library
add_library(core_interface INTERFACE)
add_library(metada::framework::core ALIAS core_interface)

# Set include directories
target_include_directories(core_interface
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link required dependencies
target_link_libraries(core_interface
    INTERFACE
        metada::framework::tools::config
        metada::framework::tools::logger
)

# Add format target
AddFormatTarget(core_interface ${CMAKE_CURRENT_SOURCE_DIR})