#ifndef METADA_BACKENDS_TOOLS_CONFIG_JSON_JSON_CONFIG_TRAITS_H_
#define METADA_BACKENDS_TOOLS_CONFIG_JSON_JSON_CONFIG_TRAITS_H_

#include "ConfigTraits.h"
#include "JsonConfig.h"

namespace metada {
namespace framework {
namespace tools {
namespace config {

/**
 * @brief Configuration traits specialization for JSON backend
 *
 * This template specialization defines the configuration backend type to be
 * used when working with JSON configuration files. It maps the void type
 * parameter to the JsonConfig implementation.
 *
 * The ConfigTraits template is used throughout the framework to provide a
 * consistent way to specify which configuration backend should be used. This
 * specialization ensures that when no specific backend is requested (void
 * parameter), the JSON backend will be used as the default.
 *
 * @see JsonConfig The JSON configuration backend implementation
 * @see ConfigTraits The base configuration traits template
 */
template <>
struct ConfigTraits<void> {
  /** @brief Type definition for the JSON configuration backend */
  using ConfigBackend = backends::tools::config::json::JsonConfig;
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_CONFIG_JSON_JSON_CONFIG_TRAITS_H_