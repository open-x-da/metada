# Create logger library with source files
add_library(logger
    LogStream.cpp
    ILogger.cpp
)

# Create alias target for easy referencing
add_library(metada::framework::common::utils::logger ALIAS logger)

# Set include directories
target_include_directories(logger
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
)

# Link against any dependencies
target_link_libraries(logger
    PUBLIC
        metada::framework::common
)

# Set C++20 properties for the target
target_compile_features(logger PUBLIC cxx_std_20)

# Add format target
AddFormatTarget(logger ${CMAKE_CURRENT_SOURCE_DIR}) 