# Create Glog logger backend library
add_library(logger_glog
    GlogLogger.cpp)
    
add_library(metada::backends::tools::logger::glog
    ALIAS logger_glog)

# Set include directories
target_include_directories(logger_glog
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

# Link dependencies
target_link_libraries(logger_glog
    PUBLIC
        metada::framework::tools::logger
        glog::glog)

# Add format target
AddFormatTarget(logger_glog
    ${CMAKE_CURRENT_SOURCE_DIR})