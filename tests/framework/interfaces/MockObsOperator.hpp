#pragma once

#include <gmock/gmock.h>

#include <string>
#include <vector>

#include "IObsOperator.hpp"
#include "utils/config/IConfig.hpp"

namespace metada::tests {

using framework::IObservation;
using framework::IObsOperator;
using framework::IState;

/**
 * @brief Mock implementation of IObsOperator for testing
 *
 * Provides mock methods for all observation operator operations:
 * - Forward operator (H): state -> observation
 * - Tangent linear (H'): increment -> observation
 * - Adjoint (H'*): observation -> increment
 */
class MockObsOperator : public IObsOperator {
 public:
  // Core operations
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));

  // Forward operator: H(x)
  MOCK_METHOD(void, apply, (const IState& state, IObservation& obs),
              (const, override));

  // Tangent linear: H'(x)δx
  MOCK_METHOD(void, applyTangentLinear, (const void* dx, IObservation& dy),
              (const, override));

  // Adjoint: H'*(x)δy
  MOCK_METHOD(void, applyAdjoint, (const IObservation& dy, void* dx),
              (const, override));

  // Metadata
  MOCK_METHOD(const std::vector<std::string>&, getRequiredStateVars, (),
              (const, override));
  MOCK_METHOD(const std::vector<std::string>&, getRequiredObsVars, (),
              (const, override));
};

}  // namespace metada::tests