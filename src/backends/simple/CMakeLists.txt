# Create the Simple backend library
add_library(simple_backends)
add_library(metada::backends::simple ALIAS simple_backends)

# Set include directories
target_include_directories(simple_backends 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Set sources
target_sources(simple_backends
    PUBLIC
        SimpleGeometry.hpp
        SimpleGeometryIterator.hpp
        SimpleState.hpp
        SimpleModel.hpp
)

target_link_libraries(simple_backends
    PRIVATE
        metada::base
)

# Add format target
AddFormatTarget(simple_backends ${CMAKE_CURRENT_SOURCE_DIR}) 