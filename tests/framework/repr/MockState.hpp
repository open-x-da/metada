/**
 * @file MockState.hpp
 * @brief Mock implementation of IState interface for testing
 * @ingroup tests
 *
 * @details
 * This mock class provides a test double for the IState interface using Google
 * Mock. It allows testing code that depends on IState by providing mock
 * implementations of all interface methods that can be configured with
 * expectations and behaviors.
 */

#ifndef METADA_TESTS_FRAMEWORK_REPR_MOCK_STATE_HPP_
#define METADA_TESTS_FRAMEWORK_REPR_MOCK_STATE_HPP_

#include <gmock/gmock.h>

#include "Config.hpp"
#include "IConfig.hpp"
#include "IState.hpp"
#include "MockConfig.hpp"

using metada::framework::tools::config::Config;
using metada::framework::tools::config::tests::MockConfig;

namespace metada {
namespace framework {
namespace repr {
namespace tests {

using ::testing::_;
using ::testing::Return;
/**
 * @brief Mock implementation of IState for testing
 *
 * @details
 * Provides mock methods for all IState interface operations:
 *
 * Core operations:
 * - initialize() - Initialize state from configuration
 * - reset() - Reset state to initial values
 * - validate() - Validate state consistency
 * - isInitialized() - Check initialization status
 *
 * Copy/Move operations:
 * - copyFrom() - Copy state from another instance
 * - moveFrom() - Move state from another instance
 * - equals() - Compare equality with another state
 *
 * Data access:
 * - getData() - Get raw pointer to data
 * - getData() const - Get const raw pointer to data
 *
 * Metadata operations:
 * - setMetadata() - Set metadata key-value pair
 * - getMetadata() - Get metadata value by key
 *
 * State information:
 * - getVariableNames() - Get names of state variables
 * - getDimensions() - Get dimensions of state space
 */
class MockState : public IState {
 private:
  const Config<MockConfig>& config_;

 public:
  // Disable default constructor
  MockState() = delete;

  // Constructor that initializes state from config
  MockState(const Config<MockConfig>& config) : config_(config) {
    // Set default behavior for initialize
    ON_CALL(*this, initialize(_)).WillByDefault(Return());
    initialize(config.backend());
  }

  MockState(const MockState& other) = delete;
  MockState& operator=(const MockState& other) = delete;

  MockState(MockState&& other) noexcept = delete;
  MockState& operator=(MockState&& other) noexcept = delete;

  // Core state operations
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(void, reset, (), (override));
  MOCK_METHOD(void, validate, (), (const, override));

  // Copy operations
  MOCK_METHOD(void, copyFrom, (const IState& other), (override));
  MOCK_METHOD(void, moveFrom, (IState && other), (override));
  MOCK_METHOD(bool, equals, (const IState& other), (const, override));

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
  const Config<MockConfig>& config() const { return config_; }
};

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_REPR_MOCK_STATE_HPP_