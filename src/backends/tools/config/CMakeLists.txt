# Add subdirectories for each backend
add_subdirectory(yaml)
add_subdirectory(json)

# Create the main config backends library
add_library(config_backends INTERFACE)
add_library(metada::backends::tools::config ALIAS config_backends)

# Set include directories for the interface library
target_include_directories(config_backends
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Configure compile definitions based on selected backend
target_compile_definitions(config_backends
    INTERFACE
        $<$<STREQUAL:${CONFIG_BACKEND},YAML>:USE_YAML_CONFIG>
        $<$<STREQUAL:${CONFIG_BACKEND},JSON>:USE_JSON_CONFIG>
)

# Link with available backends using nested IF generator expressions
target_link_libraries(config_backends
    INTERFACE
        $<IF:$<STREQUAL:${CONFIG_BACKEND},YAML>,
            metada::backends::tools::config::yaml,
            $<IF:$<STREQUAL:${CONFIG_BACKEND},JSON>,
                metada::backends::tools::config::json,
                $<ERROR:Invalid CONFIG_BACKEND value: ${CONFIG_BACKEND}. Must be either YAML or JSON>>>
) 