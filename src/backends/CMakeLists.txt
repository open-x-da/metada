# Create backend logging library
add_library(metada_backends_logging INTERFACE)

if(Glog_FOUND AND USE_GLOG)
    target_compile_definitions(metada_backends_logging
        INTERFACE
            USE_GLOG
    )
    
    target_link_libraries(metada_backends_logging
        INTERFACE
            glog::glog
            metada_core_logging
    )
else()
    target_link_libraries(metada_backends_logging
        INTERFACE
            metada_core_logging
    )
endif()

target_include_directories(metada_backends_logging
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
) 

# Add format target to ensure consistent code style across the project
add_format_target(metada_backends_logging ${CMAKE_CURRENT_SOURCE_DIR})