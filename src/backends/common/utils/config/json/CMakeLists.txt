# Create JSON config backend library
add_library(config_json
    JsonConfig.cpp)
add_library(metada::backends::common::utils::config
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

# Add precompiled headers if enabled
if(USE_PRECOMPILED_HEADERS)
    # Add standard C++ headers 
    metada_add_std_precompiled_headers(config_json PRIVATE)
    
    # Add project-specific headers
    target_precompile_headers(config_json
        PRIVATE 
        <nlohmann/json.hpp>
    )
endif()

# Add format target
AddFormatTarget(config_json
    ${CMAKE_CURRENT_SOURCE_DIR})
