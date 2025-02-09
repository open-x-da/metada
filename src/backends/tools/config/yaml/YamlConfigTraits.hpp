#ifndef METADA_BACKENDS_TOOLS_CONFIG_YAML_YAMLCONFIGTRAITS_HPP_
#define METADA_BACKENDS_TOOLS_CONFIG_YAML_YAMLCONFIGTRAITS_HPP_

#include "ConfigTraits.hpp"  // Generic ConfigTraits from framework
#include "YamlConfig.hpp"

/**
 * @file YamlConfigTraits.hpp
 * @brief Traits specialization for YAML configuration backend
 *
 * This header provides a specialization of the ConfigTraits template for the
 * YAML configuration backend. It maps the generic configuration interface to
 * the YAML-specific implementation provided by YamlConfig.
 */
template <>
struct metada::framework::tools::config::ConfigTraits<void> {
  using ConfigBackend = metada::backends::tools::config::yaml::YamlConfig;
};

#endif  // METADA_BACKENDS_TOOLS_CONFIG_YAML_YAMLCONFIGTRAITS_HPP_
