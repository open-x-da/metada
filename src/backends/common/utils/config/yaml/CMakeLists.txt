# Create YAML config backend library
add_library(config_yaml
    YamlConfig.cpp)

# Set include directories
target_include_directories(config_yaml
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link dependencies
target_link_libraries(config_yaml
    PUBLIC
        metada::framework::adapters
        yaml-cpp::yaml-cpp)

# Add format target
AddFormatTarget(config_yaml
    ${CMAKE_CURRENT_SOURCE_DIR})
