# Define the scientific representation library as interface (header-only)
add_library(adapters INTERFACE)

# Create alias target
add_library(metada::framework::adapters ALIAS adapters)

# Set up include directories for interface library
target_include_directories(adapters
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/utils/config>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/utils/logger>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/io>
        )

AddFormatTarget(adapters ${CMAKE_CURRENT_SOURCE_DIR})