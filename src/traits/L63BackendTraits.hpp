#pragma once

#include "BackendTraits.hpp"

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

#ifdef LOGGER_BACKEND_NGLOG
#include "../backends/common/utils/logger/nglog/NgLogger.hpp"
#endif

#include "../backends/common/observation/GridObservation.hpp"
#include "../backends/common/obsoperator/IdentityObsOperator.hpp"

namespace metada::traits {

struct L63BackendTag {};

template<>
struct BackendTraits<L63BackendTag> {
  // Select ConfigBackend based on CMake configuration
#ifdef CONFIG_BACKEND_JSON
  using ConfigBackend = backends::config::JsonConfig;
#elif defined(CONFIG_BACKEND_YAML)
  using ConfigBackend = backends::config::YamlConfig;
#else
  using ConfigBackend = backends::config::YamlConfig; // Default
#endif

  // Select LoggerBackend based on CMake configuration
#ifdef LOGGER_BACKEND_NGLOG
  using LoggerBackend = backends::logger::NgLogger<ConfigBackend>;
#elif defined(LOGGER_BACKEND_GLOG)
  using LoggerBackend = backends::logger::GoogleLogger<ConfigBackend>;
#elif defined(LOGGER_BACKEND_CONSOLE)
  using LoggerBackend = backends::logger::ConsoleLogger<ConfigBackend>;
#else
  using LoggerBackend = backends::logger::NgLogger<ConfigBackend>; // Default
#endif

  // TODO: Add L63-specific backend implementations when available
  // For now, using placeholder types that need to be implemented
  
  /** @brief L63 geometry backend implementation (placeholder) */
  using GeometryBackend = void; // TODO: Implement L63Geometry
  
  /** @brief L63 geometry iterator backend implementation (placeholder) */
  using GeometryIteratorBackend = void; // TODO: Implement L63GeometryIterator
  
  /** @brief L63 state vector backend implementation (placeholder) */
  using StateBackend = void; // TODO: Implement L63State
  
  /** @brief L63 ensemble backend implementation (placeholder) */
  using EnsembleBackend = void; // TODO: Implement L63Ensemble
  
  /** @brief L63 observation backend implementation */
  using ObservationBackend = backends::common::observation::GridObservation;
  
  /** @brief L63 observation iterator backend implementation */
  using ObservationIteratorBackend = metada::backends::common::observation::GridObservationIterator;
  
  /** @brief Identity observation operator backend implementation */
  using ObsOperatorBackend = backends::common::obsoperator::IdentityObsOperator<StateBackend, ObservationBackend>;
  
  /** @brief L63 model backend implementation (placeholder) */
  using ModelBackend = void; // TODO: Implement L63Model
};

} // namespace metada::traits
