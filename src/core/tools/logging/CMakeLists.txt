# Create core logging interface library and set up includes
add_library(metada_core_logging INTERFACE)
add_library(metada::core::tools::logging ALIAS metada_core_logging)

# Set include directories
target_include_directories(metada_core_logging
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add format target
AddFormatTarget(metada_core_logging ${CMAKE_CURRENT_SOURCE_DIR})