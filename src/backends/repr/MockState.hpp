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

#ifndef METADA_FRAMEWORK_REPR_TESTS_MOCK_STATE_HPP_
#define METADA_FRAMEWORK_REPR_TESTS_MOCK_STATE_HPP_

#include <gmock/gmock.h>

#include "IState.hpp"

// Forward declarations with namespace forwarding
namespace metada::framework::tools::config {
class IConfig;
template <typename T>
class Config;
}  // namespace metada::framework::tools::config

namespace metada::framework::tools::config::tests {
class MockConfig;
}  // namespace metada::framework::tools::config::tests

namespace metada::backends::repr {

// Use namespace aliases for shorter references
namespace config = metada::framework::tools::config;

using ::testing::_;
using ::testing::Return;

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
class MockState : public framework::repr::IState {
 private:
  const config::Config<config::tests::MockConfig>& config_;

 public:
  // Disable default constructor
  MockState() = delete;

  // Constructor that initializes state from config
  MockState(const config::Config<config::tests::MockConfig>& config)
      : config_(config) {
    // Set default behavior for initialize
    ON_CALL(*this, initialize(_)).WillByDefault(Return());
    initialize(config.backend());
  }

  MockState(const MockState& other) = delete;
  MockState& operator=(const MockState& other) = delete;

  MockState(MockState&& other) noexcept = delete;
  MockState& operator=(MockState&& other) noexcept = delete;

  // Core state operations
  MOCK_METHOD(void, initialize, (const config::IConfig& config), (override));
  MOCK_METHOD(void, reset, (), (override));
  MOCK_METHOD(void, validate, (), (const, override));

  // Copy operations
  MOCK_METHOD(void, copyFrom, (const framework::repr::IState& other),
              (override));
  MOCK_METHOD(void, moveFrom, (framework::repr::IState && other), (override));
  MOCK_METHOD(bool, equals, (const framework::repr::IState& other),
              (const, override));

  // Data access
  MOCK_METHOD(void*, getData, (), (override));
  MOCK_METHOD(const void*, getData, (), (const, override));

  // Metadata operations
  MOCK_METHOD(void, setMetadata,
              (const std::string& key, const std::string& value), (override));
  MOCK_METHOD(std::string, getMetadata, (const std::string& key),
              (const, override));

  // State information queries
  MOCK_METHOD(const std::vector<std::string>&, getVariableNames, (),
              (const, override));
  MOCK_METHOD(const std::vector<size_t>&, getDimensions, (), (const, override));

  // Get the config
  const config::Config<config::tests::MockConfig>& config() const {
    return config_;
  }
};

}  // namespace metada::backends::repr

#endif  // METADA_FRAMEWORK_REPR_TESTS_MOCK_STATE_HPP_