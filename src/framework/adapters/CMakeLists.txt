# Define the scientific representation library
add_library(adapters INTERFACE)

# Create alias target
add_library(metada::framework::adapters ALIAS adapters)

# Set up include directories
target_include_directories(adapters
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(adapters
    INTERFACE
        metada::framework::tools::config
        metada::framework::tools::logger
)

AddFormatTarget(adapters ${CMAKE_CURRENT_SOURCE_DIR})