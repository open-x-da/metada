# Define the scientific representation library
add_library(interfaces INTERFACE)

# Create alias target
add_library(metada::framework::interfaces ALIAS interfaces)

# Set up include directories
target_include_directories(interfaces
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(interfaces
    INTERFACE
        metada::framework::common
)

AddFormatTarget(interfaces ${CMAKE_CURRENT_SOURCE_DIR})