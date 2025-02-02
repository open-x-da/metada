# Create core logging interface library and set up includes
add_library(metada_core_tools_logger INTERFACE)
add_library(metada::core::tools::logger ALIAS metada_core_tools_logger)

# Set include directories
target_include_directories(metada_core_tools_logger
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add format target
AddFormatTarget(metada_core_tools_logger ${CMAKE_CURRENT_SOURCE_DIR})