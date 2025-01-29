# Create backend logging library
add_library(metada_backends_logging INTERFACE)

target_link_libraries(metada_backends_logging
    INTERFACE
        glog::glog
        metada_core_logging
)

target_include_directories(metada_backends_logging
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
) 