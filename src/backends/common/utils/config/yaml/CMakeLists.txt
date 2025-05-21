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

# Add precompiled headers if enabled
if(USE_PRECOMPILED_HEADERS)
    # Add standard C++ headers 
    metada_add_std_precompiled_headers(config_yaml PRIVATE)
    
    # Add project-specific headers
    target_precompile_headers(config_yaml
        PRIVATE 
        <yaml-cpp/yaml.h>
    )
endif()

# Add format target
AddFormatTarget(config_yaml
    ${CMAKE_CURRENT_SOURCE_DIR})
