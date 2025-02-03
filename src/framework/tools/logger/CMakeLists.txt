# Create framework logging interface library and set up includes
add_library(metada_framework_tools_logger INTERFACE)
add_library(metada::framework::tools::logger ALIAS metada_framework_tools_logger)

target_link_libraries(metada_framework_tools_logger
    INTERFACE
        metada::backends::tools::logger
)

if(USE_GLOG)
    target_compile_definitions(metada_framework_tools_logger
        INTERFACE
            USE_GLOG
    )
endif()

# Set include directories
target_include_directories(metada_framework_tools_logger
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add format target
AddFormatTarget(metada_framework_tools_logger ${CMAKE_CURRENT_SOURCE_DIR})