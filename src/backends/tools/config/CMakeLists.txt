# Add subdirectories for each backend
add_subdirectory(yaml)
add_subdirectory(json)

# Create the main config backends library
add_library(config_backends INTERFACE)
add_library(metada::backends::tools::config ALIAS config_backends)

# Link with available backends
target_link_libraries(config_backends
    INTERFACE
        $<$<BOOL:${yaml-cpp_FOUND}>:metada::backends::tools::config::yaml>
        $<$<BOOL:${nlohmann_json_FOUND}>:metada::backends::tools::config::json>
) 