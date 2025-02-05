#ifndef METADA_BACKENDS_TOOLS_CONFIG_YAML_YAML_CONFIG_TRAITS_H_
#define METADA_BACKENDS_TOOLS_CONFIG_YAML_YAML_CONFIG_TRAITS_H_

#include "ConfigTraits.h"
#include "YamlConfig.h"

namespace metada {
namespace framework {
namespace tools {
namespace config {

/**
 * @brief Configuration traits specialization for YAML backend
 *
 * This template specialization defines the configuration backend type to be
 * used when working with YAML configuration files. It maps the void type
 * parameter to the YamlConfig implementation.
 *
 * The ConfigTraits template is used throughout the framework to provide a
 * consistent way to specify which configuration backend should be used. This
 * specialization ensures that when no specific backend is requested (void
 * parameter), the YAML backend will be used as the default.
 *
 * @see YamlConfig The YAML configuration backend implementation
 * @see ConfigTraits The base configuration traits template
 */
template <>
struct ConfigTraits<void> {
  /** @brief Type definition for the YAML configuration backend */
  using ConfigBackend = backends::tools::config::yaml::YamlConfig;
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_CONFIG_YAML_YAML_CONFIG_TRAITS_H_