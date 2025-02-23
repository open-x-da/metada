# Define available config backends and set default
if(nlohmann_json_FOUND)
    set(CONFIG_BACKEND "JSON" CACHE STRING "Select configuration backend")
elseif(yaml-cpp_FOUND)
    set(CONFIG_BACKEND "YAML" CACHE STRING "Select configuration backend") 
else()
    message(FATAL_ERROR "No configuration backend available. Install yaml-cpp or nlohmann_json")
endif()

# Set available options based on found packages
set(CONFIG_BACKENDS "")
if(yaml-cpp_FOUND)
    list(APPEND CONFIG_BACKENDS "YAML")
endif()
if(nlohmann_json_FOUND)
    list(APPEND CONFIG_BACKENDS "JSON")
endif()

# Set the allowed values for CONFIG_BACKEND
set_property(CACHE CONFIG_BACKEND PROPERTY STRINGS ${CONFIG_BACKENDS})

# Add selected backend subdirectory
if(CONFIG_BACKEND IN_LIST CONFIG_BACKENDS)
    string(TOLOWER ${CONFIG_BACKEND} CONFIG_BACKEND_LOWER)
    add_subdirectory(${CONFIG_BACKEND_LOWER})
else()
    message(FATAL_ERROR "Invalid CONFIG_BACKEND value: ${CONFIG_BACKEND}. Must be one of: ${CONFIG_BACKENDS}")
endif()

# Create the main config backends library
add_library(config_backends INTERFACE)
add_library(metada::backends::common::utils::config ALIAS config_backends)

# Set include directories for the interface library
target_include_directories(config_backends
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Configure compile definitions and link libraries based on selected backend
target_compile_definitions(config_backends
    INTERFACE
        USE_${CONFIG_BACKEND}_CONFIG
)

target_link_libraries(config_backends
    INTERFACE
        metada::backends::common::utils::config::${CONFIG_BACKEND_LOWER}
)
