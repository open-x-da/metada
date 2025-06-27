# Create common observation library
add_library(observation INTERFACE)
add_library(metada::backends::common::observation ALIAS observation)

# Set include directories
target_include_directories(observation
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_sources(observation
    INTERFACE
        GridObservation.hpp)

# Link dependencies
target_link_libraries(observation
    INTERFACE
        metada::framework::adapters)

# Add format target
AddFormatTarget(observation ${CMAKE_CURRENT_SOURCE_DIR}) 