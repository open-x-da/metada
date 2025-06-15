# Create the Simple backend library
add_library(simple_backends INTERFACE)
add_library(metada::backends::simple ALIAS simple_backends)

# Set include directories
target_include_directories(simple_backends 
    INTERFACE 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Set headers as interface sources to show them in IDE
target_sources(simple_backends
    INTERFACE
        SimpleGeometry.hpp
        SimpleGeometryIterator.hpp
)

# Add format target
AddFormatTarget(simple_backends ${CMAKE_CURRENT_SOURCE_DIR}) 