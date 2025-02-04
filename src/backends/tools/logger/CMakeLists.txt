# Create backend logger library
add_library(metada_backends_tools_logger INTERFACE)
add_library(metada::backends::tools::logger ALIAS metada_backends_tools_logger)

# Always link with the framework logger
target_link_libraries(metada_backends_tools_logger
    INTERFACE
        metada::framework::tools::logger
)

target_include_directories(metada_backends_tools_logger
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add glog if available and enabled
if(Glog_FOUND AND USE_GLOG)   
    target_link_libraries(metada_backends_tools_logger
        INTERFACE
            glog::glog
    )
endif()

# Add format target to ensure consistent code style across the project
AddFormatTarget(metada_backends_tools_logger ${CMAKE_CURRENT_SOURCE_DIR})