# Create common observation library
add_library(common_observation INTERFACE)
add_library(metada::backends::common::observation ALIAS common_observation)

# Set include directories
target_include_directories(common_observation
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link dependencies
target_link_libraries(common_observation
    INTERFACE
        metada::framework::adapters)

# Add format target
AddFormatTarget(common_observation ${CMAKE_CURRENT_SOURCE_DIR}) 