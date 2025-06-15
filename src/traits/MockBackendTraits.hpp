/**
 * @file MockBackendTraits.hpp
 * @brief Specialization of BackendTraits for mock implementations used in testing
 * @ingroup traits
 * @author Metada Framework Team
 *
 * @details
 * This file provides a specialization of the BackendTraits template for mock
 * implementations used in testing. The MockBackendTag serves as a tag type
 * that can be used with framework adapters to inject mock implementations
 * for unit and integration testing. These mock implementations are based on
 * Google Mock and provide test-friendly interfaces for setting expectations
 * and verifying interactions.
 */

#pragma once

#include "BackendTraits.hpp"
#include "../backends/gmock/MockConfig.hpp"
#include "../backends/gmock/MockEnsemble.hpp"
#include "../backends/gmock/MockLogger.hpp"
#include "../backends/gmock/MockGeometry.hpp"
#include "../backends/gmock/MockGeometryIterator.hpp"
#include "../backends/gmock/MockState.hpp"
#include "../backends/gmock/MockModel.hpp"
#include "../backends/gmock/MockObservation.hpp"
#include "../backends/gmock/MockObsOperator.hpp"

namespace metada::traits {

/**
 * @brief Tag type for mock backend implementations
 * 
 * @details This empty struct serves as a tag type that can be used with
 * framework adapters to select mock implementations for testing purposes.
 * It enables test code to use the same adapter interfaces with mock backends
 * that can be controlled using Google Mock expectations, facilitating unit
 * and integration testing of framework components.
 */
struct MockBackendTag {};

/**
 * @brief BackendTraits specialization for mock implementations
 * 
 * @details This specialization maps the MockBackendTag to concrete mock
 * implementation types for each backend component. These mock implementations
 * are based on Google Mock and provide test-friendly interfaces for setting
 * expectations and verifying interactions in unit and integration tests.
 * Each mock component can be configured to return predefined values or
 * simulate specific behaviors for testing scenarios.
 */
template<>
struct BackendTraits<MockBackendTag> {
  /** @brief Mock configuration backend implementation for testing */
  using ConfigBackend = backends::gmock::MockConfig;
  
  /** @brief Mock logging backend implementation for testing */
  using LoggerBackend = backends::gmock::MockLogger<ConfigBackend>;
  
  /** @brief Mock geometry backend implementation for testing */
  using GeometryBackend = backends::gmock::MockGeometry;
  
  /** @brief Mock geometry iterator backend implementation for testing */
  using GeometryIteratorBackend = backends::gmock::MockGeometryIterator;
  
  /** @brief Mock state vector backend implementation for testing */
  using StateBackend = backends::gmock::MockState<ConfigBackend, GeometryBackend>;
  
  /** @brief Mock model backend implementation for testing */
  using ModelBackend = backends::gmock::MockModel<ConfigBackend, StateBackend>;
  
  /** @brief Mock observation backend implementation for testing */
  using ObservationBackend = backends::gmock::MockObservation<ConfigBackend>;

  /** @brief Mock observation operator backend implementation for testing */
  using ObsOperatorBackend = backends::gmock::MockObsOperator<ConfigBackend, StateBackend, ObservationBackend>;

  /** @brief Mock ensemble backend implementation for testing */
  using EnsembleBackend = backends::gmock::MockEnsemble<ConfigBackend, GeometryBackend>;
};

} // namespace metada::traits
