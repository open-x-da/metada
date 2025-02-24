# Create core interface library
add_library(runs INTERFACE)
add_library(metada::framework::runs ALIAS runs)

# Set include directories
target_include_directories(runs
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link required dependencies
target_link_libraries(runs
    INTERFACE
        metada::framework::common
)

# Add format target
AddFormatTarget(runs ${CMAKE_CURRENT_SOURCE_DIR})