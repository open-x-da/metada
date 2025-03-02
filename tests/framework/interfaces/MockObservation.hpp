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

#include "IObservation.hpp"
#include "utils/config/Config.hpp"

namespace metada::tests {

using metada::framework::Config;
using metada::framework::IConfig;
using metada::framework::IObservation;

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
class MockObservation : public IObservation {
 private:
  const Config<MockConfig>& config_;

 public:
  // Disable default constructor
  MockObservation() = delete;

  /**
   * @brief Constructor that initializes observation from config
   */
  explicit MockObservation(const Config<MockConfig>& config) : config_(config) {
    initialize(config.backend());
  }

  /**
   * @brief Copy constructor
   */
  MockObservation(const MockObservation& other) : config_(other.config_) {}

  /**
   * @brief Move constructor
   */
  MockObservation(MockObservation&& other) noexcept : config_(other.config_) {}

  /**
   * @brief Copy assignment operator
   */
  MockObservation& operator=([[maybe_unused]] const MockObservation& other) {
    return *this;
  }

  /**
   * @brief Move assignment operator
   */
  MockObservation& operator=(
      [[maybe_unused]] MockObservation&& other) noexcept {
    return *this;
  }

  // Lifecycle management
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(void, reset, (), (override));
  MOCK_METHOD(void, validate, (), (const, override));
  MOCK_METHOD(bool, isValid, (), (const, override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));

  // Copy/Move operations
  MOCK_METHOD(void, copyFrom, (const IObservation& other), (override));
  MOCK_METHOD(void, moveFrom, (IObservation && other), (override));
  MOCK_METHOD(bool, equals, (const IObservation& other), (const, override));

  // Data access
  MOCK_METHOD(void*, getData, (), (override));
  MOCK_METHOD(const void*, getData, (), (const, override));
  MOCK_METHOD(void*, getUncertainty, (), (override));
  MOCK_METHOD(const void*, getUncertainty, (), (const, override));
  MOCK_METHOD(size_t, getSize, (), (const, override));

  // Metadata operations
  MOCK_METHOD(void, setMetadata,
              (const std::string& key, const std::string& value), (override));
  MOCK_METHOD(std::string, getMetadata, (const std::string& key),
              (const, override));
  MOCK_METHOD(bool, hasMetadata, (const std::string& key), (const, override));

  // Observation information
  MOCK_METHOD(const std::vector<std::string>&, getVariableNames, (),
              (const, override));
  MOCK_METHOD(bool, hasVariable, (const std::string& name), (const, override));
  MOCK_METHOD(const std::vector<size_t>&, getDimensions, (), (const, override));

  // Spatiotemporal metadata
  MOCK_METHOD(void, setLocations,
              (const std::vector<std::vector<double>>& locations), (override));
  MOCK_METHOD(void, setTimes, (const std::vector<double>& timestamps),
              (override));
  MOCK_METHOD(const std::vector<std::vector<double>>&, getLocations, (),
              (const, override));
  MOCK_METHOD(const std::vector<double>&, getTimes, (), (const, override));

  // Arithmetic operations
  MOCK_METHOD(void, add, (const IObservation& other), (override));
  MOCK_METHOD(void, subtract, (const IObservation& other), (override));
  MOCK_METHOD(void, multiply, (double scalar), (override));

  // Quality control
  MOCK_METHOD(void, setQualityFlags, (const std::vector<int>& flags),
              (override));
  MOCK_METHOD(const std::vector<int>&, getQualityFlags, (), (const, override));
  MOCK_METHOD(void, setConfidenceValues, (const std::vector<double>& values),
              (override));
  MOCK_METHOD(const std::vector<double>&, getConfidenceValues, (),
              (const, override));

  // Get the config
  const Config<MockConfig>& config() const { return config_; }
};

}  // namespace metada::tests