# Create Console logger backend library
add_library(logger_console INTERFACE)

# Set include directories
target_include_directories(logger_console
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Add format target
AddFormatTarget(logger_console
    ${CMAKE_CURRENT_SOURCE_DIR})