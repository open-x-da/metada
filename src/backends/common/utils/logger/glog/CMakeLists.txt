# Create Glog logger backend library
add_library(logger_glog INTERFACE)

# Set include directories
target_include_directories(logger_glog
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link dependencies
target_link_libraries(logger_glog
    INTERFACE
        metada::framework::adapters
        glog::glog)

# Add format target
AddFormatTarget(logger_glog
    ${CMAKE_CURRENT_SOURCE_DIR})