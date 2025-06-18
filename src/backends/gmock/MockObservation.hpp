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
 * - Hierarchical data organization: type → variable → data
 *
 * @see IObservation
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>
#include <unordered_map>
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
 * - getTypeNames() - Get names of observation types
 * - getVariableNames(typeName) - Get variables for a specific type
 * - hasVariable() - Check if variable exists
 * - getSize(typeName, varName) - Get size for specific type/variable
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

  // File I/O operations - these match the concept requirements
  void loadFromFile([[maybe_unused]] const std::string& filename) {
    // TODO: Implement
  }

  void saveToFile([[maybe_unused]] const std::string& filename) const {
    // TODO: Implement
  }

  // Get the data
  void* getData() {
    if (data_.empty()) {
      return nullptr;
    }
    // Return pointer to the first type's first variable's data
    return data_.begin()->second.begin()->second.data();
  }

  const void* getData() const {
    if (data_.empty()) {
      return nullptr;
    }
    // Return pointer to the first type's first variable's data
    return data_.begin()->second.begin()->second.data();
  }

  // Get the observation type names
  std::vector<std::string> getTypeNames() const {
    std::vector<std::string> types;
    for (const auto& [type_name, _] : data_) {
      types.push_back(type_name);
    }
    return types;
  }

  // Get the variable names for a specific type
  std::vector<std::string> getVariableNames(const std::string& typeName) const {
    std::vector<std::string> variables;
    auto it = data_.find(typeName);
    if (it != data_.end()) {
      for (const auto& [var_name, _] : it->second) {
        variables.push_back(var_name);
      }
    }
    return variables;
  }

  // Get the size for a specific type/variable
  size_t getSize(const std::string& typeName,
                 const std::string& varName) const {
    std::string key = typeName + "." + varName;
    auto it = sizes_.find(key);
    if (it != sizes_.end()) {
      return it->second;
    }
    return 0;
  }

  // Get the covariance matrix
  const std::vector<double>& getCovariance() const { return covariance_; }

  // Get error value for a specific type/variable
  float getError(const std::string& typeName,
                 const std::string& varName) const {
    std::string key = typeName + "." + varName;
    auto it = errors_.find(key);
    if (it != errors_.end()) {
      return it->second;
    }
    return 0.0f;
  }

  // Get missing value for a specific type/variable
  float getMissingValue(const std::string& typeName,
                        const std::string& varName) const {
    std::string key = typeName + "." + varName;
    auto it = missing_values_.find(key);
    if (it != missing_values_.end()) {
      return it->second;
    }
    return -999.0f;
  }

  // Test helper methods
  void setTypeNames(const std::vector<std::string>& types) {
    typeNames_ = types;
  }

  void setVariableNames(const std::string& typeName,
                        const std::vector<std::string>& variables) {
    variableNames_[typeName] = variables;
  }

  void setSize(const std::string& typeName, const std::string& varName,
               size_t size) {
    std::string key = typeName + "." + varName;
    sizes_[key] = size;
  }

  void setData(const std::vector<double>& data) {
    // Store data in the first available type/variable or create default
    if (data_.empty()) {
      data_["default"]["default"] = data;
      setSize("default", "default", data.size());
    } else {
      // Use the first type and variable
      auto& first_type = data_.begin()->second;
      auto& first_var = first_type.begin()->second;
      first_var = data;
      setSize(data_.begin()->first, first_type.begin()->first, data.size());
    }
  }

  void setCovariance(const std::vector<double>& cov) { covariance_ = cov; }

  void setError(const std::string& typeName, const std::string& varName,
                float error) {
    std::string key = typeName + "." + varName;
    errors_[key] = error;
  }

  void setMissingValue(const std::string& typeName, const std::string& varName,
                       float missingValue) {
    std::string key = typeName + "." + varName;
    missing_values_[key] = missingValue;
  }

  // Helper method to set up hierarchical data structure
  void setupHierarchicalData(const std::string& typeName,
                             const std::string& varName,
                             const std::vector<double>& data) {
    data_[typeName][varName] = data;
    setSize(typeName, varName, data.size());
  }

 private:
  const ConfigBackend& config_;
  std::vector<std::string> typeNames_;
  std::unordered_map<std::string, std::vector<std::string>> variableNames_;
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::vector<double>>>
      data_;
  std::unordered_map<std::string, size_t> sizes_;
  std::unordered_map<std::string, float> errors_;
  std::unordered_map<std::string, float> missing_values_;
  std::vector<double> covariance_;
};

}  // namespace metada::backends::gmock