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
 * - Clone operations
 * - Move semantics
 * - Data access and manipulation
 *
 * @see IState
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <cmath>
#include <iostream>
#include <unordered_map>

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of IState for testing
 *
 * @details
 * Provides mock methods for all IState interface operations, organized into
 * the following categories:
 *
 * @par Core Operations
 * - initialize() - Initialize state from configuration and geometry
 * - zero() - Set all values to zero
 * - clone() - Create a deep copy of the state
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
 * - setData() - Set the internal data vector (helper method for testing)
 *
 * @par State Information
 * - getVariableNames() - Get names of state variables
 * - hasVariable() - Check if state contains a specific variable
 * - size() - Get total size of the state vector
 * - setVariables() - Set variable names (helper method for testing)
 *
 * @par Memory Management
 * - Move constructor and assignment operator
 * - Deleted copy constructor and assignment operator (non-copyable)
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
template <typename ConfigBackend, typename GeometryBackend>
class MockState {
 public:
  // Disable default constructor
  MockState() = delete;

  // Destructor
  ~MockState() = default;

  // Copy constructor
  MockState(const MockState& other) = delete;

  // Copy assignment operator
  MockState& operator=(const MockState& other) = delete;

  // Move constructor
  MockState(MockState&& other) noexcept
      : config_(other.config_),
        geometry_(other.geometry_),
        variableNames_(std::move(other.variableNames_)),
        data_(std::move(other.data_)) {
    other.variableNames_.clear();
    other.data_.clear();
  }

  // Move assignment operator
  MockState& operator=(MockState&& other) noexcept {
    if (this != &other) {
      variableNames_ = std::move(other.variableNames_);
      data_ = std::move(other.data_);
      other.variableNames_.clear();
      other.data_.clear();
    }
    return *this;
  }

  // Constructor that initializes state from config and geometry
  MockState(const ConfigBackend& config, const GeometryBackend& geometry)
      : config_(config), geometry_(geometry) {
    initialize();
  }

  // Clone operation
  std::unique_ptr<MockState> clone() const {
    auto cloned = std::make_unique<MockState>(config_, geometry_);
    cloned->setData(data_);                // Copy the data
    cloned->setVariables(variableNames_);  // Copy the variable names
    return cloned;
  }

  // Core state operations
  MOCK_METHOD(void, initialize, ());
  MOCK_METHOD(void, zero, ());

  // Compare operations
  MOCK_METHOD(bool, equals, (const MockState& other), (const));

  // Arithmetic operations
  MOCK_METHOD(void, add, (const MockState& other));
  MOCK_METHOD(void, subtract, (const MockState& other));
  MOCK_METHOD(void, multiply, (double scalar));

  // Implement dot product calculation
  double dot(const MockState& other) const {
    if (data_.size() != other.data_.size()) return 0.0;
    double result = 0.0;
    for (size_t i = 0; i < data_.size(); ++i) {
      result += data_[i] * other.data_[i];
    }
    return result;
  }

  // Implement norm calculation
  double norm() const {
    double result = 0.0;
    for (size_t i = 0; i < data_.size(); ++i) {
      result += data_[i] * data_[i];
    }
    return std::sqrt(result);
  }

  // File I/O operations
  MOCK_METHOD(void, saveToFile, (const std::string& filename), (const));

  // Data access
  void* getData() { return data_.empty() ? nullptr : data_.data(); }

  const void* getData() const { return data_.empty() ? nullptr : data_.data(); }

  // State information queries
  const std::vector<std::string>& getVariableNames() const {
    return variableNames_;
  }

  /**
   * @brief Get the total size of the state vector
   * @return Total number of elements in the state vector
   */
  size_t size() const { return data_.size(); }

  // Test helper methods
  void setVariables(const std::vector<std::string>& variables) {
    variableNames_ = variables;
  }

  void setData(const std::vector<double>& data) { data_ = data; }

  // Get the config
  const ConfigBackend& config() const { return config_; }
  // Get the geometry
  const GeometryBackend& geometry() const { return geometry_; }

 private:
  const ConfigBackend& config_;
  const GeometryBackend& geometry_;
  std::vector<std::string> variableNames_;
  std::vector<double> data_;
};

// Output operator for MockState
template <typename ConfigBackend, typename GeometryBackend>
inline std::ostream& operator<<(
    std::ostream& os, const MockState<ConfigBackend, GeometryBackend>& state) {
  os << "MockState(size=" << state.size() << ")";
  return os;
}

}  // namespace metada::backends::gmock