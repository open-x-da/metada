#pragma once
#include "BackendTraits.hpp"

// Include backend headers based on CMake configuration
#ifdef CONFIG_BACKEND_JSON
#include "../backends/common/utils/config/json/JsonConfig.hpp"
#endif

#ifdef CONFIG_BACKEND_YAML
#include "../backends/common/utils/config/yaml/YamlConfig.hpp"
#endif

#ifdef LOGGER_BACKEND_CONSOLE
#include "../backends/common/utils/logger/console/ConsoleLogger.hpp"
#endif

#ifdef LOGGER_BACKEND_NGLOG
#include "../backends/common/utils/logger/nglog/NgLogger.hpp"
#endif

#include "../backends/lite/LiteBEC.hpp"
#include "../backends/lite/LiteGeometry.hpp"
#include "../backends/lite/LiteModel.hpp"
#include "../backends/lite/LiteObs.hpp"
#include "../backends/lite/LiteObsOperator.hpp"
#include "../backends/lite/LiteState.hpp"

namespace metada::traits {

/**
 * @brief Lite backend tag for concrete testing
 */
struct LiteBackendTag {};

/**
 * @brief Backend traits specialization for lite backends
 */
template <>
struct BackendTraits<LiteBackendTag> {
  // Config backend selection
#ifdef CONFIG_BACKEND_JSON
  using ConfigBackend = backends::config::JsonConfig;
#elif defined(CONFIG_BACKEND_YAML)
  using ConfigBackend = backends::config::YamlConfig;
#else
  using ConfigBackend = backends::config::YamlConfig; // Default
#endif

  // Logger backend selection
#ifdef LOGGER_BACKEND_NGLOG
  using LoggerBackend = backends::logger::NgLogger<ConfigBackend>;
#elif defined(LOGGER_BACKEND_CONSOLE)
  using LoggerBackend = backends::logger::ConsoleLogger<ConfigBackend>;
#else
  using LoggerBackend = backends::logger::NgLogger<ConfigBackend>; // Default
#endif

  using StateBackend = metada::backends::lite::LiteState;
  using GeometryBackend = metada::backends::lite::LiteGeometry;
  using ObservationBackend = metada::backends::lite::LiteObs;
  using ObsOperatorBackend = metada::backends::lite::LiteObsOperator;
  using ModelBackend = metada::backends::lite::LiteModel;
  using BackgroundErrorCovarianceBackend =
      metada::backends::lite::LiteBEC;
};

}  // namespace metada::traits 