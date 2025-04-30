/**
 * @file MockModel.hpp
 * @brief Mock implementation of IModel interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IModel interface using Google
 * Mock. It allows testing code that depends on IModel by providing mock
 * implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 *
 * @see IModel
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>

#include "DateTime.hpp"

namespace metada::backends::gmock {

// Using declaration for DateTime
using core::DateTime;

/**
 * @brief Mock implementation of IModel for testing
 *
 * @details
 * Provides mock methods for all IModel interface operations, organized
 * into the following categories:
 *
 * @par Lifecycle Management
 * - initialize() - Initialize model from configuration
 * - reset() - Reset model to initial state
 * - finalize() - Finalize model
 * - isInitialized() - Check if model is initialized
 *
 * @par Parameter Management
 * - getParameter() - Get a parameter value by name
 * - setParameter() - Set a parameter value
 *
 * @par Model Execution
 * - run() - Run the model from initial to final state
 *
 * @tparam ConfigBackend The backend configuration type
 * @tparam StateBackend The backend state type
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
template <typename ConfigBackend, typename StateBackend>
class MockModel {
 public:
  /**
   * @brief Default constructor is deleted
   *
   * MockModel objects must be constructed with a Config object.
   */
  MockModel() = delete;

  /**
   * @brief Destructor
   */
  ~MockModel() = default;

  /**
   * @brief Copy constructor is deleted
   */
  MockModel(const MockModel&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  MockModel& operator=(const MockModel&) = delete;

  /**
   * @brief Move constructor
   *
   * @param other The other MockModel object to move from
   */
  MockModel(MockModel&& other) noexcept : config_(std::move(other.config_)) {}

  /**
   * @brief Move assignment operator
   *
   * @param other The other MockModel object to move from
   * @return Reference to this MockModel object
   */
  MockModel& operator=(MockModel&& other) noexcept {
    config_ = std::move(other.config_);
    return *this;
  }

  /**
   * @brief Constructor taking a configuration object
   *
   * @param config Configuration object with model parameters
   */
  explicit MockModel(const ConfigBackend& config) : config_(config) {}

  // Lifecycle management
  MOCK_METHOD(void, initialize, (const ConfigBackend& config));
  MOCK_METHOD(void, reset, ());
  MOCK_METHOD(void, finalize, ());
  MOCK_METHOD(bool, isInitialized, (), (const));

  // Parameter management
  MOCK_METHOD(std::string, getParameter, (const std::string& name), (const));
  MOCK_METHOD(void, setParameter,
              (const std::string& name, const std::string& value));

  // Model execution
  MOCK_METHOD(void, run,
              (const StateBackend& initialState, StateBackend& finalState,
               const DateTime& startTime, const DateTime& endTime));

 private:
  const ConfigBackend& config_; /**< Reference to the configuration object */
};

}  // namespace metada::backends::gmock