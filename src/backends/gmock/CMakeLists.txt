# Create JSON config backend library
add_library(gmock_config INTERFACE)
add_library(metada::backends::gmock ALIAS gmock_config)

# Set include directories
target_include_directories(gmock_config
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Add format target
AddFormatTarget(gmock_config
    ${CMAKE_CURRENT_SOURCE_DIR})
