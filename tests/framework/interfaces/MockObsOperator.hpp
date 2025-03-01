#pragma once

#include <gmock/gmock.h>

#include <string>
#include <vector>

#include "IObsOperator.hpp"

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
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, finalize, (), (override));

  // Forward operator
  MOCK_METHOD(void, apply, (const IState&, IObservation&), (const, override));

  // Tangent linear operator
  MOCK_METHOD(void, applyTangentLinear, (const IIncrement&, IObservation&),
              (const, override));

  // Adjoint operator
  MOCK_METHOD(void, applyAdjoint, (const IObservation&, IIncrement&),
              (const, override));

  // Error handling
  MOCK_METHOD(void, setObservationError, (const IObservation&), (override));
  MOCK_METHOD(double, getObservationError, (const IObservation&),
              (const, override));

  // Configuration
  MOCK_METHOD(void, setParameter, (const std::string&, double), (override));
  MOCK_METHOD(double, getParameter, (const std::string&), (const, override));

  // Required variables
  MOCK_METHOD(const std::vector<std::string>&, getRequiredStateVariables, (),
              (const, override));
  MOCK_METHOD(const std::vector<std::string>&, getRequiredObsVariables, (),
              (const, override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));
};

}  // namespace metada::tests