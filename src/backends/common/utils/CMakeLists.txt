# Common utilities library
#
# Provides configuration and logging utilities with selectable backends

# Create the common utils library (only public interface)
add_library(common_utils INTERFACE)
add_library(metada::backends::common::utils ALIAS common_utils)

# Add include directories for utility headers
target_include_directories(common_utils
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
)

# First build configuration utilities (required by logger)
add_subdirectory(config)

# Then build logger utilities that depend on config
add_subdirectory(logger)

# Link both internal components to the public interface
target_link_libraries(common_utils
    INTERFACE
        config_backend
        logger_backend
)

# Add a message showing what backends are used with configuration details
message(STATUS "Selected backend configurations:")
message(STATUS "  - Config backend: ${CONFIG_BACKEND} (Configurable via -DCONFIG_BACKEND=<YAML|JSON>)")
message(STATUS "  - Logger backend: ${LOGGER_BACKEND} (Configurable via -DLOGGER_BACKEND=<NGLOG|CONSOLE>)")