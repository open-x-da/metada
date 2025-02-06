#ifndef METADA_BACKENDS_TOOLS_CONFIG_CONFIG_BACKEND_SELECTOR_H_
#define METADA_BACKENDS_TOOLS_CONFIG_CONFIG_BACKEND_SELECTOR_H_

#include "ConfigTraits.h"  // Generic ConfigTraits

#ifdef USE_YAML_CONFIG
#include "yaml/YamlConfig.h"
#elif defined(USE_JSON_CONFIG)
#include "json/JsonConfig.h"
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
 * - USE_YAML_CONFIG: Selects the YAML configuration backend
 * - USE_JSON_CONFIG: Selects the JSON configuration backend
 *
 * Only one backend can be selected at a time. The selection is determined at
 * build time through CMake configuration.
 *
 * @see YamlConfig YAML configuration backend implementation
 * @see JsonConfig JSON configuration backend implementation
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

#endif  // METADA_BACKENDS_TOOLS_CONFIG_CONFIG_BACKEND_SELECTOR_H_
