if(CONFIG_BACKEND STREQUAL "JSON")
    # Create JSON config backend library
    add_library(config_json
        JsonConfig.cc)
    add_library(metada::backends::tools::config::json 
        ALIAS config_json)

    # Set include directories
    target_include_directories(config_json
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR})

    # Link dependencies
    target_link_libraries(config_json
        PUBLIC
            metada::framework::tools::config
            nlohmann_json::nlohmann_json)

    # Add format target
    AddFormatTarget(config_json 
        ${CMAKE_CURRENT_SOURCE_DIR})
endif()