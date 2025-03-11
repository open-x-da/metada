# Create common interface library
add_library(common INTERFACE)
add_library(metada::framework::common ALIAS common)

# Set include directories
target_include_directories(common
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Add subdirectories for utilities
add_subdirectory(utils/logger)

# Add format target
AddFormatTarget(common ${CMAKE_CURRENT_SOURCE_DIR})
