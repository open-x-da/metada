/**
 * @file MockObservation.hpp
 * @brief Mock implementation of IObservation interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IObservation interface using
 * Google Mock. It allows testing code that depends on IObservation by providing
 * mock implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 *
 * @see IObservation
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>
#include <vector>

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of IObservation for testing
 *
 * @details
 * Provides mock methods for all IObservation interface operations, organized
 * into the following categories:
 *
 * @par Lifecycle Management
 * - initialize() - Initialize observation from configuration
 * - reset() - Reset observation to initial values
 * - validate() - Validate observation consistency
 * - isValid() - Check if observation is valid
 * - isInitialized() - Check initialization status
 *
 * @par Copy/Move Operations
 * - copyFrom() - Copy observation from another instance
 * - moveFrom() - Move observation from another instance
 * - equals() - Compare equality with another observation
 *
 * @par Data Access
 * - getData() - Get raw pointer to data
 * - getData() const - Get const raw pointer to data
 * - getUncertainty() - Get raw pointer to uncertainty data
 * - getUncertainty() const - Get const raw pointer to uncertainty data
 * - getSize() - Get size of observation data
 *
 * @par Metadata Operations
 * - setMetadata() - Set metadata key-value pair
 * - getMetadata() - Get metadata value by key
 * - hasMetadata() - Check if metadata key exists
 *
 * @par Observation Information
 * - getVariableNames() - Get names of observation variables
 * - hasVariable() - Check if variable exists
 * - getDimensions() - Get dimensions of observation space
 *
 * @par Spatiotemporal Metadata
 * - setLocations() - Set spatial locations for observations
 * - setTimes() - Set timestamps for observations
 * - getLocations() - Get spatial locations
 * - getTimes() - Get timestamps
 *
 * @par Arithmetic Operations
 * - add() - Add another observation
 * - subtract() - Subtract another observation
 * - multiply() - Multiply by scalar
 *
 * @par Quality Control
 * - setQualityFlags() - Set quality control flags
 * - getQualityFlags() - Get quality control flags
 * - setConfidenceValues() - Set confidence values
 * - getConfidenceValues() - Get confidence values
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
template <typename ConfigBackend>
class MockObservation {
 public:
  // Disable default constructor
  MockObservation() = delete;

  // Destructor
  ~MockObservation() = default;

  // Copy constructor
  MockObservation(const MockObservation& other) = delete;

  /**
   * @brief Copy assignment operator
   */
  MockObservation& operator=(const MockObservation& other) = delete;

  /**
   * @brief Move constructor
   */
  MockObservation(MockObservation&& other) noexcept : config_(other.config_) {
    // Explicit Move Constructor (Even Without Data Members)
  }

  /**
   * @brief Move assignment operator
   */
  MockObservation& operator=(
      [[maybe_unused]] MockObservation&& other) noexcept {
    // Explicit Move Assignment (Even Without Data Members)
    return *this;
  }

  /**
   * @brief Constructor that initializes observation from config
   */
  explicit MockObservation(const ConfigBackend& config) : config_(config) {
    initialize();
  }

  // Clone operation
  std::unique_ptr<MockObservation> clone() const {
    auto cloned = std::make_unique<MockObservation>(config_);
    return cloned;
  }

  // Lifecycle management
  MOCK_METHOD(void, initialize, ());

  MOCK_METHOD(void, applyQC, ());
  MOCK_METHOD(bool, equals, (const MockObservation& other), (const));

  // Arithmetic operations
  MOCK_METHOD(void, add, (const MockObservation& other));
  MOCK_METHOD(void, subtract, (const MockObservation& other));
  MOCK_METHOD(void, multiply, (double scalar));

  void loadFromFile([[maybe_unused]] const std::string& filename) {
    // TODO: Implement
  }

  void saveToFile([[maybe_unused]] const std::string& filename) const {
    // TODO: Implement
  }

  // Get the data
  void* getData() { return &data_; }

  const void* getData() const { return &data_; }

  // Get the variable names
  const std::vector<std::string>& getVariableNames() const {
    return variableNames_;
  }

  // Get the dimensions
  const std::vector<size_t>& getDimensions(const std::string& name) const {
    return dimensions_.at(name);
  }

  // Get the covariance matrix
  const std::vector<double>& getCovariance() const { return covariance_; }

  // Test helper methods
  void setVariables(const std::vector<std::string>& variables) {
    variableNames_ = variables;
  }

  void setDimensions(const std::string& name, const std::vector<size_t>& dims) {
    dimensions_[name] = dims;
  }

  void setData(const std::vector<double>& data) { data_ = data; }

  void setCovariance(const std::vector<double>& cov) { covariance_ = cov; }

 private:
  const ConfigBackend& config_;
  std::vector<std::string> variableNames_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
  std::vector<double> data_;
  std::vector<double> covariance_;
};

}  // namespace metada::backends::gmock