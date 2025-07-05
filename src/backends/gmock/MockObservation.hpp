/**
 * @file MockObservation.hpp
 * @brief Mock implementation of observation backend for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the observation backend interface
 * using Google Mock. It allows testing code that depends on observation
 * backends by providing mock implementations of all interface methods that can
 * be configured with expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 * - Point-based observations with location information
 * - Iteration capabilities
 *
 * @see Observation
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>
#include <unordered_map>
#include <vector>

#include "Location.hpp"
#include "MockObservationIterator.hpp"

namespace metada::backends::gmock {

using metada::framework::CoordinateSystem;
using metada::framework::Location;

// Use the MockObservationPoint from MockObservationIterator.hpp
using ObservationPoint = MockObservationPoint;

/**
 * @brief Mock implementation of observation backend for testing
 *
 * @details
 * Provides mock methods for all observation backend interface operations,
 * organized into the following categories:
 *
 * @par Lifecycle Management
 * - initialize() - Initialize observation from configuration
 * - applyQC() - Apply quality control
 * - equals() - Compare equality with another observation
 *
 * @par Copy/Move Operations
 * - clone() - Clone observation
 *
 * @par Data Access and Iteration
 * - begin()/end() - Iteration support
 * - size() - Get number of observations
 * - operator[] - Direct indexing
 * - getData() - Get raw pointer to data
 * - getData<T>() - Template data access
 *
 * @par Variable Information
 * - getTypeNames() - Get names of observation types
 * - getVariableNames(typeName) - Get variables for a specific type
 *
 * @par Arithmetic Operations
 * - add() - Add another observation
 * - subtract() - Subtract another observation
 * - multiply() - Multiply by scalar
 *
 * @par File I/O Operations
 * - loadFromFile() - Load observations from file
 * - saveToFile() - Save observations to file
 *
 * @par Geographic Filtering
 * - getObservationsInBox() - Get observations in geographic bounding box
 * - getObservationsInVerticalRange() - Get observations in vertical range
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
template <typename ConfigBackend>
class MockObservation {
 public:
  // Type aliases for concept compliance
  using iterator_type = MockObservationIterator;
  using value_type = ObservationPoint;

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
  MockObservation(MockObservation&& other) noexcept
      : config_(other.config_),
        observations_(std::move(other.observations_)),
        type_variable_map_(std::move(other.type_variable_map_)),
        covariance_(std::move(other.covariance_)) {}

  /**
   * @brief Move assignment operator
   */
  MockObservation& operator=(MockObservation&& other) noexcept {
    if (this != &other) {
      observations_ = std::move(other.observations_);
      type_variable_map_ = std::move(other.type_variable_map_);
      covariance_ = std::move(other.covariance_);
    }
    return *this;
  }

  /**
   * @brief Constructor that initializes observation from config
   */
  explicit MockObservation(const ConfigBackend& config) : config_(config) {
    // Set up default type/variable mapping for testing
    type_variable_map_["obs_A"]["temperature"] = {0, 1, 2};
    type_variable_map_["obs_A"]["pressure"] = {0, 1, 2};

    // Add some default observations for testing
    addObservation(Location(45.0, -120.0, 1000.0), 25.5, 0.5);
    addObservation(Location(46.0, -121.0, 850.0), 15.2, 0.3);
    addObservation(Location(47.0, -122.0, 500.0), -5.8, 0.7);

    // Set default covariance
    covariance_ = {1.0, 0.0, 0.0, 1.0};
  }

  // Iteration capabilities
  MockObservationIterator begin() const {
    return MockObservationIterator(&observations_, 0);
  }

  MockObservationIterator end() const {
    return MockObservationIterator(&observations_, observations_.size());
  }

  size_t size() const { return observations_.size(); }

  const ObservationPoint& operator[](size_t index) const {
    return observations_[index];
  }

  // Clone operation
  std::unique_ptr<MockObservation> clone() const {
    auto cloned = std::make_unique<MockObservation>(config_);
    cloned->observations_ = observations_;
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

  // File I/O operations
  MOCK_METHOD(void, loadFromFile,
              (const std::string& filename, double error,
               double missing_value));
  MOCK_METHOD(void, saveToFile, (const std::string& filename), (const));

  // Geographic filtering
  MOCK_METHOD(std::vector<ObservationPoint>, getObservationsInBox,
              (double min_lat, double max_lat, double min_lon, double max_lon),
              (const));
  MOCK_METHOD(std::vector<ObservationPoint>, getObservationsInVerticalRange,
              (double min_level, double max_level), (const));

  // Get the data
  void* getData() {
    if (observations_.empty()) {
      return nullptr;
    }
    return const_cast<void*>(static_cast<const void*>(observations_.data()));
  }

  const void* getData() const {
    if (observations_.empty()) {
      return nullptr;
    }
    return static_cast<const void*>(observations_.data());
  }

  // Template data access for backward compatibility
  template <typename T>
  T getData() const {
    if constexpr (std::is_same_v<T, std::vector<double>>) {
      std::vector<double> values;
      values.reserve(observations_.size());
      for (const auto& obs : observations_) {
        values.push_back(obs.value);
      }
      return values;
    }
    return T{};
  }

  // Get the observation type names
  std::vector<std::string> getTypeNames() const {
    std::vector<std::string> types;
    for (const auto& [type_name, _] : type_variable_map_) {
      types.push_back(type_name);
    }
    return types;
  }

  // Get the variable names for a specific type
  std::vector<std::string> getVariableNames(const std::string& typeName) const {
    std::vector<std::string> variables;
    auto it = type_variable_map_.find(typeName);
    if (it != type_variable_map_.end()) {
      for (const auto& [var_name, _] : it->second) {
        variables.push_back(var_name);
      }
    }
    return variables;
  }

  // Get the observation error variances
  std::vector<double> getCovariance() const {
    std::vector<double> errors;
    errors.reserve(observations_.size());
    for (const auto& obs : observations_) {
      errors.push_back(obs.error * obs.error);  // Convert error to variance
    }
    return errors;
  }

  // Test helper methods
  void setObservations(const std::vector<ObservationPoint>& obs) {
    observations_ = obs;
  }

  void addObservation(const Location& location, double value, double error) {
    observations_.emplace_back(location, value, error);
  }

  void setCovariance(const std::vector<double>& cov) { covariance_ = cov; }

  void setTypeVariableMap(
      const std::unordered_map<
          std::string, std::unordered_map<std::string, std::vector<size_t>>>&
          map) {
    type_variable_map_ = map;
  }

 private:
  const ConfigBackend& config_;
  std::vector<ObservationPoint> observations_;

  // Mapping from type/variable to observation indices (for backward
  // compatibility)
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::vector<size_t>>>
      type_variable_map_;

  std::vector<double> covariance_;
};

}  // namespace metada::backends::gmock