# Define the scientific representation library
add_library(adapters INTERFACE)

# Create alias target
add_library(metada::framework::adapters ALIAS adapters)

# Set up include directories
target_include_directories(adapters
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

AddFormatTarget(adapters ${CMAKE_CURRENT_SOURCE_DIR})