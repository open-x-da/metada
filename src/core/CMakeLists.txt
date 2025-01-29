# Create core logging interface library
add_library(metada_core_logging INTERFACE)

target_include_directories(metada_core_logging
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
) 