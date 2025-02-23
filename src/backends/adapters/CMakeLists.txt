# Create Mock State backend library
add_library(adapters_backend INTERFACE)
add_library(metada::backends::adapters_backend
    ALIAS adapters_backend)

# Set include directories
target_include_directories(adapters_backend
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

# Link dependencies
target_link_libraries(adapters_backend
    INTERFACE
        metada::framework::adapters
        metada::framework::tools::config)

# Add format target
AddFormatTarget(adapters_backend
    ${CMAKE_CURRENT_SOURCE_DIR})