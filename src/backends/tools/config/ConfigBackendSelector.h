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

template <>
struct ConfigTraits<void> {
#ifdef USE_YAML_CONFIG
    using ConfigBackend = backends::tools::config::yaml::YamlConfig;
#elif defined(USE_JSON_CONFIG)
    using ConfigBackend = backends::tools::config::json::JsonConfig;
#endif
};

}  // namespace config
}  // namespace tools
}  // namespace framework
}  // namespace metada

#endif  // METADA_BACKENDS_TOOLS_CONFIG_CONFIG_BACKEND_SELECTOR_H_
