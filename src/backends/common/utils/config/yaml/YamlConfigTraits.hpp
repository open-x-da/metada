#pragma once

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
struct metada::framework::common::utils::config::ConfigTraits<void> {
  using ConfigBackend = metada::backends::common::utils::config::yaml::YamlConfig;
};
