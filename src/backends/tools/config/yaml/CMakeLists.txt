# Create YAML config backend library
add_library(config_yaml
    YamlConfig.cpp)
add_library(metada::backends::tools::config::yaml
    ALIAS config_yaml)

# Set include directories
target_include_directories(config_yaml
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

# Link dependencies
target_link_libraries(config_yaml
    PUBLIC
        metada::framework::tools::config
        yaml-cpp::yaml-cpp)

# Add format target
AddFormatTarget(config_yaml
    ${CMAKE_CURRENT_SOURCE_DIR})
