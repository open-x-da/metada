#ifndef METADA_BACKENDS_TOOLS_CONFIG_CONFIGBACKENDSELECTOR_HPP_
#define METADA_BACKENDS_TOOLS_CONFIG_CONFIGBACKENDSELECTOR_HPP_

#include "ConfigTraits.hpp"  // Generic ConfigTraits

#ifdef USE_YAML_CONFIG
#include "yaml/YamlConfig.hpp"
#elif defined(USE_JSON_CONFIG)
#include "json/JsonConfig.hpp"
#endif

namespace metada {
namespace framework {
namespace tools {
namespace config {

/**
 * @brief Configuration backend selector traits
 * 
 * This specialization of ConfigTraits for void provides compile-time selection
 * of the appropriate configuration backend based on build configuration.
 * 
 * The selection is controlled by the following preprocessor definitions:
 * - USE_YAML_CONFIG: Selects the YAML configuration backend which provides
 *   configuration storage and access using YAML format
 * - USE_JSON_CONFIG: Selects the JSON configuration backend which provides
 *   configuration storage and access using JSON format
 *
 * Only one backend can be selected at a time through CMake configuration.
 * The selected backend determines the format used for:
 * - Loading configuration from files/strings
 * - Accessing configuration values
 * - Saving configuration to files
 *
 * @see YamlConfig YAML configuration backend using yaml-cpp library
 * @see JsonConfig JSON configuration backend using nlohmann::json library
 * @see ConfigTraits Generic configuration traits template
 */
template <>
struct ConfigTraits<void> {
#ifdef USE_YAML_CONFIG
    /** @brief Selected configuration backend type (YAML implementation) */
    using ConfigBackend = backends::tools::config::yaml::YamlConfig;
#elif defined(USE_JSON_CONFIG)
    /** @brief Selected configuration backend type (JSON implementation) */
    using ConfigBackend = backends::tools::config::json::JsonConfig;
#endif
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_CONFIG_CONFIGBACKENDSELECTOR_HPP_
