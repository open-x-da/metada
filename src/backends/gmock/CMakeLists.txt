# Create gmock backend library
add_library(gmock_backends INTERFACE)
add_library(metada::backends::gmock ALIAS gmock_backends)

# Set include directories
target_include_directories(gmock_backends
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Set headers as interface sources to show them in IDE
target_sources(gmock_backends
    INTERFACE
        MockConfig.hpp
        MockGeometry.hpp
        MockGridPoint.hpp
        MockGeometryIterator.hpp
        MockLogger.hpp
        MockModel.hpp
        MockState.hpp
        MockObservation.hpp)

# Add format target
AddFormatTarget(gmock_backends ${CMAKE_CURRENT_SOURCE_DIR})
