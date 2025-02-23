# Create common interface library
add_library(common INTERFACE)
add_library(metada::framework::common ALIAS common)

# Set include directories
target_include_directories(common
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Add all headers recursively
file(GLOB_RECURSE PUBLIC_HEADERS 
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)

# Set public headers property
set_target_properties(common PROPERTIES
    PUBLIC_HEADER "${PUBLIC_HEADERS}"
)

# Add format target
AddFormatTarget(common ${CMAKE_CURRENT_SOURCE_DIR})
