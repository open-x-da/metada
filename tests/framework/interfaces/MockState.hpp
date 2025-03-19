/**
 * @file MockState.hpp
 * @brief Mock implementation of IState interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IState interface using Google
 * Mock. It allows testing code that depends on IState by providing mock
 * implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 *
 * @see IState
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include "IState.hpp"

namespace metada::tests {

using framework::Config;
using framework::IConfig;
using framework::IState;

/**
 * @brief Mock implementation of IState for testing
 *
 * @details
 * Provides mock methods for all IState interface operations, organized into
 * the following categories:
 *
 * @par Core Operations
 * - initialize() - Initialize state from configuration
 * - reset() - Reset state to initial values
 * - validate() - Validate state consistency
 * - isInitialized() - Check initialization status
 *
 * @par Copy/Move Operations
 * - copyFrom() - Copy state from another instance
 * - moveFrom() - Move state from another instance
 * - equals() - Compare equality with another state
 *
 * @par Data Access
 * - getData() - Get raw pointer to data
 * - getData() const - Get const raw pointer to data
 *
 * @par Metadata Operations
 * - setMetadata() - Set metadata key-value pair
 * - getMetadata() - Get metadata value by key
 *
 * @par State Information
 * - getVariableNames() - Get names of state variables
 * - getDimensions() - Get dimensions of state space
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
class MockState : public IState {
 public:
  // Disable default constructor
  MockState() = delete;

  // Destructor
  ~MockState() override = default;

  // Copy constructor
  MockState(const MockState& other) = delete;

  // Copy assignment operator
  MockState& operator=(const MockState& other) = delete;

  // Move constructor
  MockState(MockState&& other) noexcept : config_(other.config_) {
    // Explicit Move Constructor (Even Without Data Members)
  }

  // Move assignment operator
  MockState& operator=(MockState&& other) {
    // Explicit Move Assignment (Even Without Data Members)
    return *this;
  }

  // Constructor that initializes state from config
  explicit MockState(const IConfig& config) : config_(config) { initialize(); }

  // Clone operation
  MOCK_METHOD(std::unique_ptr<IState>, clone, (), (const, override));

  // Core state operations
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, zero, (), (override));

  // Compare operations
  MOCK_METHOD(bool, equals, (const IState& other), (const, override));

  // Data access
  MOCK_METHOD(void*, getData, (), (override));
  MOCK_METHOD(const void*, getData, (), (const, override));

  // State information queries
  MOCK_METHOD(const std::vector<std::string>&, getVariableNames, (),
              (const, override));
  MOCK_METHOD(bool, hasVariable, (const std::string& name), (const, override));
  MOCK_METHOD(const std::vector<size_t>&, getDimensions,
              (const std::string& name), (const, override));

  // Arithmetic operations
  MOCK_METHOD(void, add, (const IState& other), (override));
  MOCK_METHOD(void, subtract, (const IState& other), (override));
  MOCK_METHOD(void, multiply, (double scalar), (override));
  MOCK_METHOD(double, dot, (const IState& other), (const, override));
  MOCK_METHOD(double, norm, (), (const, override));

  // Get the config
  const IConfig& config() const { return config_; }

 private:
  const IConfig& config_;
};

}  // namespace metada::tests