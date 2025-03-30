# Create JSON config backend library
add_library(gmock_utils INTERFACE)
add_library(metada::backends::gmock ALIAS gmock_utils)

# Set include directories
target_include_directories(gmock_utils
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Add format target
AddFormatTarget(gmock_utils
    ${CMAKE_CURRENT_SOURCE_DIR})
