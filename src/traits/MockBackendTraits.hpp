/**
 * @file MockBackendTraits.hpp
 * @brief Specialization of BackendTraits for mock implementations
 * @ingroup traits
 * @author Metada Framework Team
 *
 * @details
 * This file provides a specialization of the BackendTraits template for mock
 * implementations used in testing. The MockBackendTag serves as a tag type
 * that can be used with framework adapters to inject mock implementations
 * for unit and integration testing.
 */

#pragma once

#include "BackendTraits.hpp"
#include "../backends/gmock/MockConfig.hpp"
#include "../backends/gmock/MockLogger.hpp"
#include "../backends/gmock/MockGeometry.hpp"
#include "../backends/gmock/MockGeometryIterator.hpp"
#include "../backends/gmock/MockState.hpp"
#include "../backends/gmock/MockModel.hpp"
#include "../backends/gmock/MockObservation.hpp"
#include "../backends/gmock/MockObsOperator.hpp"
#include "../backends/gmock/MockObsIO.hpp"

namespace metada::traits {

/**
 * @brief Tag type for mock backend implementations
 * 
 * @details This empty struct serves as a tag type that can be used with
 * framework adapters to select mock implementations for testing purposes.
 * It allows test code to use the same adapter interfaces with mock backends
 * that can be controlled using Google Mock expectations.
 */
struct MockBackendTag {};

/**
 * @brief BackendTraits specialization for mock implementations
 * 
 * @details This specialization maps the MockBackendTag to concrete mock
 * implementation types for each backend component. These mock implementations
 * are based on Google Mock and provide test-friendly interfaces for setting
 * expectations and verifying interactions in unit and integration tests.
 */
template<>
struct BackendTraits<MockBackendTag> {
  /** @brief Mock implementation of configuration backend */
  using ConfigBackend = backends::gmock::MockConfig;
  
  /** @brief Mock implementation of logging backend */
  using LoggerBackend = backends::gmock::MockLogger<ConfigBackend>;
  
  /** @brief Mock implementation of geometry backend */
  using GeometryBackend = backends::gmock::MockGeometry<ConfigBackend>;
  
  /** @brief Mock implementation of geometry iterator backend */
  using GeometryIteratorBackend = backends::gmock::MockGeometryIterator;
  
  /** @brief Mock implementation of state vector backend */
  using StateBackend = backends::gmock::MockState<ConfigBackend>;
  
  /** @brief Mock implementation of model backend */
  using ModelBackend = backends::gmock::MockModel<ConfigBackend, StateBackend>;
  
  /** @brief Mock implementation of observation backend */
  using ObservationBackend = backends::gmock::MockObservation<ConfigBackend>;

  /** @brief Mock implementation of observation operator backend */
  using ObsOperatorBackend = backends::gmock::MockObsOperator<ConfigBackend, StateBackend, ObservationBackend>;
  
  /** @brief Mock implementation of observation I/O backend */
  using ObsIOBackend = backends::gmock::MockObsIO<ConfigBackend>;
};

} // namespace metada::traits
