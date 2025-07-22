# Create the Simple backend library (header-only)
add_library(simple_backends INTERFACE)
add_library(metada::backends::simple ALIAS simple_backends)

# Set include directories
target_include_directories(simple_backends 
    INTERFACE 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link dependencies
target_link_libraries(simple_backends
    INTERFACE
        metada::base
)

# Add format target for the header files
AddFormatTarget(simple_backends ${CMAKE_CURRENT_SOURCE_DIR}) 