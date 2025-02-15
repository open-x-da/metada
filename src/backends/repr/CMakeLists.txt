# Create Mock State backend library
add_library(repr_mock INTERFACE)
add_library(metada::backends::repr::mock
    ALIAS repr_mock)

# Set include directories
target_include_directories(repr_mock
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

# Link dependencies
target_link_libraries(repr_mock
    INTERFACE
        metada::framework::repr
        metada::framework::tools::config)

# Add format target
AddFormatTarget(repr_mock
    ${CMAKE_CURRENT_SOURCE_DIR})