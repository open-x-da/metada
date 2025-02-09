#ifndef METADA_BACKENDS_TOOLS_CONFIG_JSON_JSONCONFIGTRAITS_HPP_
#define METADA_BACKENDS_TOOLS_CONFIG_JSON_JSONCONFIGTRAITS_HPP_

#include "ConfigTraits.hpp"  // Generic ConfigTraits from framework
#include "JsonConfig.hpp"

/**
 * @file JsonConfigTraits.hpp
 * @brief Traits specialization for JSON configuration backend
 *
 * This header provides a specialization of the ConfigTraits template for the
 * JSON configuration backend. It maps the generic configuration interface to
 * the JSON-specific implementation provided by JsonConfig.
 */
template <>
struct metada::framework::tools::config::ConfigTraits<void> {
  using ConfigBackend = metada::backends::tools::config::json::JsonConfig;
};

#endif  // METADA_BACKENDS_TOOLS_CONFIG_JSON_JSONCONFIGTRAITS_HPP_