/**
 * @file MockTimeStepper.hpp
 * @brief Mock implementation of ITimeStepper interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the ITimeStepper interface using
 * Google Mock. It allows testing code that depends on ITimeStepper by providing
 * mock implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * @see ITimeStepper
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>

#include "ITimeStepper.hpp"

namespace metada::tests {

using framework::IState;
using framework::ITimeStepper;

/**
 * @brief Mock implementation of ITimeStepper for testing
 *
 * Provides mock methods for all time stepper operations including:
 * - Basic stepping operations
 * - Time step configuration
 * - Adaptive time stepping
 * - Time management
 */
class MockTimeStepper : public ITimeStepper {
 public:
  // Core stepping operations
  MOCK_METHOD(void, step, (const IState& currentState, IState& nextState),
              (override));

  // Time step configuration
  MOCK_METHOD(double, getTimeStep, (), (const, override));
  MOCK_METHOD(void, setTimeStep, (double dt), (override));
  MOCK_METHOD(double, getMaxStableTimeStep, (const IState& state), (override));

  // Adaptive time stepping
  MOCK_METHOD(bool, usesAdaptiveTimeStep, (), (const, override));
  MOCK_METHOD(void, setAdaptiveTimeStep, (bool enable), (override));

  // Time management
  MOCK_METHOD(double, getCurrentTime, (), (const, override));
  MOCK_METHOD(void, setCurrentTime, (double time), (override));
  MOCK_METHOD(std::string, getTimeSteppingScheme, (), (const, override));

  // These methods don't exist in the parent interface, so we shouldn't mark
  // them override
  MOCK_METHOD(std::vector<std::string>, getSupportedSchemes, (), (const));
  MOCK_METHOD(void, setTimeSteppingScheme, (const std::string& scheme));
};

}  // namespace metada::tests