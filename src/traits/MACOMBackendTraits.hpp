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

// Include MACOM backend headers first to provide context
#include "../backends/macom/MACOMGeometry.hpp"
#include "../backends/macom/MACOMGeometryIterator.hpp"
#include "../backends/macom/MACOMState.hpp"
#include "../backends/macom/MACOMModel.hpp"

// Now include the MACOM-specific observation-related headers
#include "../backends/macom/MACOMObservation.hpp"
#include "../backends/macom/MACOMObsOperator.hpp"

// Keep common headers for fallback compatibility
#include "../backends/common/observation/GridObservation.hpp"
#include "../backends/common/obsoperator/IdentityObsOperator.hpp"

#ifdef USE_MPI
#include "../backends/macom/MACOMParallel.hpp"
#endif

namespace metada::traits {

/**
 * @brief Tag type for MACOM backend implementations.
 * @details This empty struct serves as a tag for selecting the MACOM backend.
 */
struct MACOMBackendTag {};

/**
 * @brief BackendTraits specialization for MACOM implementations.
 * @details This specialization maps the MACOMBackendTag to concrete MACOM
 * implementation types for each backend component.
 */
template<>
struct BackendTraits<MACOMBackendTag> {
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
    using LoggerBackend = backends::logger::GoogleLogger<ConfigBackend>; // Default
#endif

    /** @brief MACOM geometry backend implementation */
    using GeometryBackend = backends::macom::MACOMGeometry<ConfigBackend>;

    /** @brief MACOM geometry iterator backend implementation */
    using GeometryIteratorBackend = backends::macom::MACOMGeometryIterator<ConfigBackend>;

    /** @brief MACOM state vector backend implementation */
    using StateBackend = backends::macom::MACOMState<ConfigBackend, GeometryBackend>;

    /** @brief Observation backend implementation for MACOM */
    using ObservationBackend = backends::macom::MACOMObservation;
    
    /** @brief Observation operator backend implementation for MACOM */
    using ObsOperatorBackend = backends::macom::MACOMObsOperator<StateBackend, ObservationBackend>;

    /** @brief MACOM model backend implementation */
    using ModelBackend = backends::macom::MACOMModel<ConfigBackend, StateBackend>;

#ifdef USE_MPI
    /** @brief MACOM parallel backend implementation */
    using ParallelBackend = backends::macom::MACOMParallel;
#endif
};

}     // namespace metada::traits