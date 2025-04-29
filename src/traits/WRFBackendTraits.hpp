#pragma once
#include <BackendTraits.hpp>

// Include backend headers based on CMake configuration
#ifdef CONFIG_BACKEND_JSON
#include "../backends/common/utils/config/json/JsonConfig.hpp"
#endif

#ifdef CONFIG_BACKEND_YAML
#include "../backends/common/utils/config/yaml/YamlConfig.hpp"
#endif

#ifdef LOGGER_BACKEND_GLOG
#include "../backends/common/utils/logger/glog/GoogleLogger.hpp"
#endif

#ifdef LOGGER_BACKEND_CONSOLE
#include "../backends/common/utils/logger/console/ConsoleLogger.hpp"
#endif

#include "../backends/wrf/WRFGeometry.hpp"
#include "../backends/wrf/WRFGeometryIterator.hpp"
#include "../backends/wrf/WRFState.hpp"
#include "../backends/wrf/WRFModel.hpp"

namespace metada::traits {

struct WRFBackendTag {};

template<>
struct BackendTraits<WRFBackendTag> {
      // Select ConfigBackend based on CMake configuration
#ifdef CONFIG_BACKEND_JSON
  using ConfigBackend = backends::config::JsonConfig;
#elif defined(CONFIG_BACKEND_YAML)
  using ConfigBackend = backends::config::YamlConfig;
#else
  using ConfigBackend = backends::config::YamlConfig; // Default
#endif

  // Select LoggerBackend based on CMake configuration
#ifdef LOGGER_BACKEND_GLOG
  using LoggerBackend = backends::logger::GoogleLogger<ConfigBackend>;
#elif defined(LOGGER_BACKEND_CONSOLE)
  using LoggerBackend = backends::logger::ConsoleLogger<ConfigBackend>;
#else
  using LoggerBackend = backends::logger::GoogleLogger<ConfigBackend>; // Default
#endif

  using GeometryBackend = backends::wrf::WRFGeometry<ConfigBackend>;
  using GeometryIteratorBackend =
      backends::wrf::WRFGeometryIterator<ConfigBackend>;
  using StateBackend = backends::wrf::WRFState<ConfigBackend>;
  using ModelBackend = backends::wrf::WRFModel<ConfigBackend>;
};

}  // namespace metada::traits