# Lite backend implementation
# Simple concrete backends for testing and development

# Create the Lite backend library (header-only)
add_library(lite_backends INTERFACE)
add_library(metada::backends::lite ALIAS lite_backends)

# Set include directories
target_include_directories(lite_backends 
    INTERFACE 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link dependencies
target_link_libraries(lite_backends
    INTERFACE
        metada::base
)

# Add format target for the header files
AddFormatTarget(lite_backends ${CMAKE_CURRENT_SOURCE_DIR})

set(LITE_BACKEND_HEADERS
    LiteState.hpp
    LiteGeometry.hpp
    LiteModel.hpp
    LiteObs.hpp
    LiteObsOperator.hpp
    LiteBEC.hpp
)

# Install headers
install(FILES ${LITE_BACKEND_HEADERS}
        DESTINATION include/metada/backends/lite) 