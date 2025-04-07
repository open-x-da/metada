# Configuration backend selection
#
# Supports JSON (nlohmann_json) or YAML (yaml-cpp) configuration formats.
# At least one backend must be available.

# Set default backend based on availability
if(NOT DEFINED CONFIG_BACKEND)
    if(nlohmann_json_FOUND AND yaml-cpp_FOUND)
        set(CONFIG_BACKEND "YAML" CACHE STRING "Select configuration backend (YAML or JSON)")
        set_property(CACHE CONFIG_BACKEND PROPERTY STRINGS "YAML" "JSON")
    elseif(nlohmann_json_FOUND)
        set(CONFIG_BACKEND "JSON" CACHE STRING "Select configuration backend")
    elseif(yaml-cpp_FOUND)
        set(CONFIG_BACKEND "YAML" CACHE STRING "Select configuration backend")
    else()
        message(FATAL_ERROR "No configuration backend available. Install yaml-cpp or nlohmann_json")
    endif()
else()
    # Validate user-selected backend
    string(TOUPPER ${CONFIG_BACKEND} CONFIG_BACKEND_UPPER)
    if(CONFIG_BACKEND_UPPER STREQUAL "JSON" AND NOT nlohmann_json_FOUND)
        message(FATAL_ERROR "JSON configuration backend requested but nlohmann_json is not available")
    elseif(CONFIG_BACKEND_UPPER STREQUAL "YAML" AND NOT yaml-cpp_FOUND)
        message(FATAL_ERROR "YAML configuration backend requested but yaml-cpp is not available")
    elseif(NOT CONFIG_BACKEND_UPPER STREQUAL "JSON" AND NOT CONFIG_BACKEND_UPPER STREQUAL "YAML")
        message(FATAL_ERROR "Invalid configuration backend: ${CONFIG_BACKEND}. Valid options are: YAML, JSON")
    endif()
    
    # Normalize backend name
    set(CONFIG_BACKEND ${CONFIG_BACKEND_UPPER})
endif()

# Create the config library (internal implementation)
add_library(config_backend INTERFACE)

# Include utils.hpp and other common headers
target_include_directories(config_backend
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
)

# Add the selected backend as a dependency
string(TOLOWER ${CONFIG_BACKEND} CONFIG_BACKEND_LOWER)
add_subdirectory(${CONFIG_BACKEND_LOWER})

# Link the selected backend to the interface library
target_link_libraries(config_backend
    INTERFACE
        config_${CONFIG_BACKEND_LOWER}
)