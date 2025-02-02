# Create backend logger library
add_library(metada_backends_tools_logger INTERFACE)
add_library(metada::backends::tools::logger ALIAS metada_backends_tools_logger)

if(Glog_FOUND AND USE_GLOG)
    target_compile_definitions(metada_backends_tools_logger
        INTERFACE
            USE_GLOG
    )
    
    target_link_libraries(metada_backends_tools_logger
        INTERFACE
            glog::glog
            metada::core::tools::logger
    )
else()
    target_link_libraries(metada_backends_tools_logger
        INTERFACE
            metada::core::tools::logger
    )
endif()

target_include_directories(metada_backends_tools_logger
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
) 

# Add format target to ensure consistent code style across the project
AddFormatTarget(metada_backends_tools_logger ${CMAKE_CURRENT_SOURCE_DIR})