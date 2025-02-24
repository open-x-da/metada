# Create JSON config backend library
add_library(config_json
    JsonConfig.cpp)
add_library(metada::backends::common::utils::config::json
    ALIAS config_json)

# Set include directories
target_include_directories(config_json
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link dependencies
target_link_libraries(config_json
    PUBLIC
        metada::framework::common
        nlohmann_json::nlohmann_json)

# Add format target
AddFormatTarget(config_json
    ${CMAKE_CURRENT_SOURCE_DIR})
