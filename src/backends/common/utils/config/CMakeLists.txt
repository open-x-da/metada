# Configuration backend selection
#
# Supports JSON (nlohmann_json) or YAML (yaml-cpp) configuration formats.
# At least one backend must be available.

if(nlohmann_json_FOUND)
    set(CONFIG_BACKEND "JSON" CACHE STRING "Select configuration backend")
elseif(yaml-cpp_FOUND)
    set(CONFIG_BACKEND "YAML" CACHE STRING "Select configuration backend")
else()
    message(FATAL_ERROR "No configuration backend available. Install yaml-cpp or nlohmann_json")
endif()

string(TOLOWER ${CONFIG_BACKEND} CONFIG_BACKEND_LOWER)
add_subdirectory(${CONFIG_BACKEND_LOWER})