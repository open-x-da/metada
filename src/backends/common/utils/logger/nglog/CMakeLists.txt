# Create nglog logger backend library
add_library(logger_nglog INTERFACE)

# Set include directories
target_include_directories(logger_nglog
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link dependencies
target_link_libraries(logger_nglog
    INTERFACE
        metada::framework::adapters
        ng-log::ng-log)

# Add format target
AddFormatTarget(logger_nglog ${CMAKE_CURRENT_SOURCE_DIR}) 