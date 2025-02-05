# Create config interface library
add_library(config_interface INTERFACE)
add_library(metada::framework::tools::config ALIAS config_interface)

# Set include directories
target_include_directories(config_interface
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add format target
AddFormatTarget(config_interface ${CMAKE_CURRENT_SOURCE_DIR})