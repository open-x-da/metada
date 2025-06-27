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

#include "../backends/simple/SimpleGeometry.hpp"
#include "../backends/simple/SimpleGeometryIterator.hpp"
#include "../backends/simple/SimpleState.hpp"
#include "../backends/simple/SimpleEnsemble.hpp"
#include "../backends/common/observation/GridObservation.hpp"
#include "../backends/common/obsoperator/IdentityObsOperator.hpp"
#include "../backends/simple/SimpleModel.hpp"

namespace metada::traits {

/**
 * @brief Tag type for simple backend implementations
 * 
 * @details This empty struct serves as a tag type that can be used with
 * framework adapters to select simple implementations for testing and development
 * purposes. It provides a lightweight alternative to more complex backend
 * implementations.
 */
struct SimpleBackendTag {};

/**
 * @brief BackendTraits specialization for simple implementations
 * 
 * @details This specialization maps the SimpleBackendTag to concrete simple
 * implementation types for each backend component. These implementations provide
 * basic functionality suitable for testing and development scenarios where
 * full-featured backends are not required.
 */
template<>
struct BackendTraits<SimpleBackendTag> {
  /** @brief Configuration backend implementation */
  #ifdef CONFIG_BACKEND_JSON
  using ConfigBackend = backends::config::JsonConfig;
  #elif defined(CONFIG_BACKEND_YAML)
  using ConfigBackend = backends::config::YamlConfig;
  #else
  using ConfigBackend = backends::config::YamlConfig; // Default
  #endif

  /** @brief Logger backend implementation */
  #ifdef LOGGER_BACKEND_NGLOG
  using LoggerBackend = backends::logger::NgLogger<ConfigBackend>;
  #elif defined(LOGGER_BACKEND_GLOG)
  using LoggerBackend = backends::logger::GoogleLogger<ConfigBackend>;
  #elif defined(LOGGER_BACKEND_CONSOLE)
  using LoggerBackend = backends::logger::ConsoleLogger<ConfigBackend>;
  #else
  using LoggerBackend = backends::logger::NgLogger<ConfigBackend>; // Default
  #endif

  /** @brief Simple geometry backend implementation */
  using GeometryBackend = backends::simple::SimpleGeometry;
  
  /** @brief Simple geometry iterator backend implementation */
  using GeometryIteratorBackend = backends::simple::SimpleGeometryIterator<typename GeometryBackend::const_iterator>;
  
  /** @brief Simple state vector backend implementation */
  using StateBackend = backends::simple::SimpleState;
  
  /** @brief Simple ensemble backend implementation */
  using EnsembleBackend = backends::simple::SimpleEnsemble;
  
  /** @brief Simple observation backend implementation */
  using ObservationBackend = backends::common::observation::GridObservation;
  
  /** @brief Identity observation operator backend implementation */
  using ObsOperatorBackend = backends::common::obsoperator::IdentityObsOperator<StateBackend>;
  
  /** @brief Simple model backend implementation */
  using ModelBackend = backends::simple::SimpleModel;
};

} // namespace metada::traits