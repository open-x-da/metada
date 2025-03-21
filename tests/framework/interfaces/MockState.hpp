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

#include <unordered_map>

#include "IState.hpp"

namespace metada::tests {

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
 * - zero() - Set all values to zero
 *
 * @par Comparison Operations
 * - equals() - Compare equality with another state
 *
 * @par Arithmetic Operations
 * - add() - Add another state to this one
 * - subtract() - Subtract another state from this one
 * - multiply() - Multiply this state by a scalar
 * - dot() - Calculate dot product with another state
 * - norm() - Calculate the norm of this state
 *
 * @par Data Access
 * - getData() - Get raw pointer to data
 * - getData() const - Get const raw pointer to data
 *
 * @par State Information
 * - getVariableNames() - Get names of state variables
 * - hasVariable() - Check if state contains a specific variable
 * - getDimensions() - Get dimensions of state space for a variable
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
  MockState& operator=([[maybe_unused]] MockState&& other) {
    // Explicit Move Assignment (Even Without Data Members)
    return *this;
  }

  // Constructor that initializes state from config
  explicit MockState(const IConfig& config) : config_(config) { initialize(); }

  // Clone operation
  std::unique_ptr<MockState> clone() const {
    auto cloned = std::make_unique<MockState>(config_);
    return cloned;
  }

  // Core state operations
  MOCK_METHOD(void, initialize, (), (override));
  MOCK_METHOD(void, zero, (), (override));

  // Compare operations
  MOCK_METHOD(bool, equals, (const IState& other), (const, override));

  // Arithmetic operations
  MOCK_METHOD(void, add, (const IState& other), (override));
  MOCK_METHOD(void, subtract, (const IState& other), (override));
  MOCK_METHOD(void, multiply, (double scalar), (override));
  MOCK_METHOD(double, dot, (const IState& other), (const, override));
  MOCK_METHOD(double, norm, (), (const, override));

  // Data access

  void* getData() { return data_.empty() ? nullptr : data_.data(); }

  const void* getData() const { return data_.empty() ? nullptr : data_.data(); }

  // State information queries
  const std::vector<std::string>& getVariableNames() const {
    return variableNames_;
  }

  const std::vector<size_t>& getDimensions(const std::string& name) const {
    auto it = dimensions_.find(name);
    if (it == dimensions_.end()) {
      throw std::runtime_error("Variable " + name +
                               " not found in state dimensions");
    }
    return it->second;
  }

  // Helper method to set dimensions for a variable
  void setDimensions(const std::string& name, const std::vector<size_t>& dims) {
    if (std::find(variableNames_.begin(), variableNames_.end(), name) ==
        variableNames_.end()) {
      throw std::runtime_error(
          "Cannot set dimensions for non-existent variable: " + name);
    }
    dimensions_[name] = dims;
  }

  // Test helper methods
  void setVariables(const std::vector<std::string>& variables) {
    variableNames_ = variables;
  }

  void setData(const std::vector<double>& data) { data_ = data; }

  // Get the config
  const IConfig& config() const { return config_; }

 private:
  const IConfig& config_;
  std::vector<std::string> variableNames_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
  std::vector<double> data_;
};

}  // namespace metada::tests