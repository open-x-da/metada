# Create gmock backend library
add_library(gmock_backends INTERFACE)
add_library(metada::backends::gmock ALIAS gmock_backends)

# Set include directories
target_include_directories(gmock_backends
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(gmock_backends
    INTERFACE
        GTest::gtest_main
        GTest::gmock)

# Add format target
AddFormatTarget(gmock_backends ${CMAKE_CURRENT_SOURCE_DIR})
