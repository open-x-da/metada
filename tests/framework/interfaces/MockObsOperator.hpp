#pragma once

#include <gmock/gmock.h>

#include <string>
#include <vector>

#include "IObsOperator.hpp"
#include "utils/config/IConfig.hpp"

namespace metada::tests {

using framework::IIncrement;
using framework::IObservation;
using framework::IObsOperator;
using framework::IState;

/**
 * @brief Mock implementation of IObsOperator for testing
 *
 * This mock class provides test doubles for all methods in the IObsOperator
 * interface, allowing tests to verify interactions with observation operators
 * without requiring real implementations.
 */
class MockObsOperator : public IObsOperator {
 public:
  // Core operations
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));

  // Forward operator
  MOCK_METHOD(void, apply, (const IState& state, IObservation& obs),
              (const, override));

  // Tangent linear and adjoint
  MOCK_METHOD(void, applyTangentLinear,
              (const IIncrement& dx, IObservation& dy), (const, override));
  MOCK_METHOD(void, applyAdjoint, (const IObservation& dy, IIncrement& dx),
              (const, override));

  // Metadata
  MOCK_METHOD(const std::vector<std::string>&, getRequiredStateVars, (),
              (const, override));
  MOCK_METHOD(const std::vector<std::string>&, getRequiredObsVars, (),
              (const, override));
};

}  // namespace metada::tests