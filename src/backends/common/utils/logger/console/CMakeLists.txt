# Create Console logger backend library
add_library(logger_console INTERFACE)

add_library(metada::backends::common::utils::logger::console
    ALIAS logger_console)

# Set include directories
target_include_directories(logger_console
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link dependencies
target_link_libraries(logger_console
    INTERFACE
        metada::framework::common)

# Add format target
AddFormatTarget(logger_console
    ${CMAKE_CURRENT_SOURCE_DIR})