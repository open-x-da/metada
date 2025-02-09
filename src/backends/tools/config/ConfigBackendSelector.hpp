#ifndef METADA_BACKENDS_TOOLS_CONFIG_CONFIGBACKENDSELECTOR_HPP_
#define METADA_BACKENDS_TOOLS_CONFIG_CONFIGBACKENDSELECTOR_HPP_

/**
 * @file ConfigBackendSelector.hpp
 * @brief Header file for selecting the configuration backend based on compile flags
 *
 * This header provides conditional inclusion of either YAML or JSON configuration
 * backend traits based on compile-time flags (USE_YAML_CONFIG or USE_JSON_CONFIG).
 * It allows the configuration system to be built with different backends without
 * code changes.
 */

#ifdef USE_YAML_CONFIG
#include "yaml/YamlConfigTraits.hpp"
#elif defined(USE_JSON_CONFIG)
#include "json/JsonConfigTraits.hpp"
#endif

#endif  // METADA_BACKENDS_TOOLS_CONFIG_CONFIGBACKENDSELECTOR_HPP_
