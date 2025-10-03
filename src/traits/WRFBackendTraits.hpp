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

#include "../backends/common/io/BufrObsIO.hpp"

#include "../backends/wrf/WRFGeometry.hpp"
#include "../backends/wrf/WRFGeometryIterator.hpp"
#include "../backends/wrf/WRFState.hpp"
#include "../backends/wrf/WRFModel.hpp"
#include "../backends/wrf/WRFBackgroundErrorCovariance.hpp"
#include "../backends/wrf/WRFObservation.hpp"
#include "../backends/wrf/WRFObsOperator.hpp"
#include "../backends/common/observation/GridObservation.hpp"
#include "../backends/common/observation/PrepBUFRObservation.hpp"
#include "../backends/common/observation/PrepBUFRObservationIterator.hpp"

namespace metada::traits {

/**
 * @brief Tag type for WRF backend implementations
 * 
 * @details This empty struct serves as a tag type that can be used with
 * framework adapters to select WRF-specific implementations. It enables the use
 * of WRF model components and data structures within the framework.
 */
struct WRFBackendTag {};

/**
 * @brief BackendTraits specialization for WRF implementations
 * 
 * @details This specialization maps the WRFBackendTag to concrete WRF-specific
 * implementation types for each backend component. These implementations provide
 * integration with the Weather Research and Forecasting (WRF) model system,
 * including its geometry, state vector, and model components.
 */
template<>
struct BackendTraits<WRFBackendTag> {
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
  #elif defined(LOGGER_BACKEND_CONSOLE)
  using LoggerBackend = backends::logger::ConsoleLogger<ConfigBackend>;
  #else
  using LoggerBackend = backends::logger::NgLogger<ConfigBackend>; // Default
  #endif

  using ObsIOBackend = backends::io::BufrObsIO<ConfigBackend>;

  /** @brief WRF-specific geometry backend implementation */
  using GeometryBackend = backends::wrf::WRFGeometry<ConfigBackend>;
  
  /** @brief WRF-specific geometry iterator backend implementation */
  using GeometryIteratorBackend = backends::wrf::WRFGeometryIterator<ConfigBackend>;
  
  /** @brief WRF-specific state vector backend implementation */
  using StateBackend = backends::wrf::WRFState<ConfigBackend, GeometryBackend>;
  
  /** @brief WRF-specific ensemble backend implementation (placeholder) */
  using EnsembleBackend = void; // TODO: Implement WRFEnsemble
  
  /** @brief WRF observation backend implementation */
  using ObservationBackend = backends::wrf::WRFObservation<GeometryBackend>;
  
  /** @brief WRF observation iterator backend implementation */
  using ObservationIteratorBackend = metada::backends::common::observation::PrepBUFRObservationIterator;
  
  /** @brief WRF observation operator backend implementation */
  using ObsOperatorBackend = backends::wrf::WRFObsOperator<StateBackend, ObservationBackend>;
  
  /** @brief WRF-specific model backend implementation */
  using ModelBackend = backends::wrf::WRFModel<ConfigBackend, StateBackend>;
  
  /** @brief WRF-specific background error covariance backend implementation */
  using BackgroundErrorCovarianceBackend = backends::wrf::WRFBackgroundErrorCovariance<ConfigBackend>;
};

}  // namespace metada::traits