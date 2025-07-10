# Create the Simple backend library (header-only)
add_library(simple_backends INTERFACE)
add_library(metada::backends::simple ALIAS simple_backends)

# Set include directories
target_include_directories(simple_backends 
    INTERFACE 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Set sources (for IDEs)
target_sources(simple_backends
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SimpleGeometry.hpp>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SimpleGeometryIterator.hpp>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SimpleState.hpp>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SimpleModel.hpp>
)

target_link_libraries(simple_backends
    INTERFACE
        metada::base
)

# Add format target
AddFormatTarget(simple_backends ${CMAKE_CURRENT_SOURCE_DIR}) 