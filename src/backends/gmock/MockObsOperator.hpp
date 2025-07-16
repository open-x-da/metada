#pragma once

#include <gmock/gmock.h>

#include <string>
#include <vector>

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of ObsOperatorBackend for testing
 *
 * @tparam ConfigBackend Type of the configuration backend
 * @tparam StateBackend Type of the state backend
 * @tparam ObsBackend Type of the observation backend
 *
 * This class provides mock methods for all observation operator operations:
 * - Initialization and configuration
 * - Forward operator (H): mapping from model state to observation space
 * - Access to required state and observation variables
 *
 * It is primarily intended for unit testing components that depend on
 * observation operators.
 */
template <typename ConfigBackend, typename StateBackend, typename ObsBackend>
class MockObsOperator {
 public:
  /**
   * @brief Default constructor is deleted
   */
  MockObsOperator() = delete;

  /**
   * @brief Destructor
   */
  ~MockObsOperator() = default;

  /**
   * @brief Copy constructor is deleted
   * @param other Source object
   */
  MockObsOperator(const MockObsOperator& other) = delete;

  /**
   * @brief Copy assignment operator is deleted
   * @param other Source object
   * @return Reference to this object
   */
  MockObsOperator& operator=(const MockObsOperator& other) = delete;

  /**
   * @brief Move constructor
   * @param other Source object
   */
  MockObsOperator(MockObsOperator&& other) noexcept : config_(other.config_) {}

  /**
   * @brief Move assignment operator
   * @param other Source object
   * @return Reference to this object
   */
  MockObsOperator& operator=(
      [[maybe_unused]] MockObsOperator&& other) noexcept {
    return *this;
  }

  /**
   * @brief Constructor with configuration
   * @param config Configuration backend reference
   */
  MockObsOperator(const ConfigBackend& config) : config_(config) {
    initialize(config);
  }

  /**
   * @brief Initialize the observation operator
   * @param config Configuration backend reference
   */
  MOCK_METHOD(void, initialize, (const ConfigBackend& config));

  /**
   * @brief Check if the observation operator is initialized
   * @return True if initialized, false otherwise
   */
  MOCK_METHOD(bool, isInitialized, (), (const));

  /**
   * @brief Apply the forward observation operator H(x)
   * @param state Input state
   * @param obs Output observation
   */
  MOCK_METHOD(std::vector<double>, apply,
              (const StateBackend& state, const ObsBackend& obs), (const));

  /**
   * @brief Get the list of required state variables
   * @return Vector of state variable names
   */
  MOCK_METHOD(const std::vector<std::string>&, getRequiredStateVars, (),
              (const));

  /**
   * @brief Get the list of required observation variables
   * @return Vector of observation variable names
   */
  MOCK_METHOD(const std::vector<std::string>&, getRequiredObsVars, (), (const));

  /**
   * @brief Apply tangent linear observation operator H(dx)
   * @param state_increment Input state increment
   * @param reference_state Reference state for linearization
   * @param obs Reference observation
   * @return Vector of observation space increments
   */
  MOCK_METHOD(std::vector<double>, applyTangentLinear,
              (const StateBackend& state_increment,
               const StateBackend& reference_state, const ObsBackend& obs),
              (const));

  /**
   * @brief Apply adjoint observation operator H^T(dy)
   * @param obs_increment Input observation increment
   * @param reference_state Reference state for adjoint
   * @param state_increment Output state increment
   * @param obs Reference observation
   */
  MOCK_METHOD(void, applyAdjoint,
              (const std::vector<double>& obs_increment,
               const StateBackend& reference_state,
               StateBackend& state_increment, const ObsBackend& obs),
              (const));

  /**
   * @brief Check if tangent linear and adjoint operators are available
   * @return True if linearization is supported
   */
  MOCK_METHOD(bool, supportsLinearization, (), (const));

  /**
   * @brief Check if the observation operator is linear
   * @return True if the operator is linear
   */
  MOCK_METHOD(bool, isLinear, (), (const));

 private:
  /** @brief Reference to the configuration backend */
  const ConfigBackend& config_;
};

}  // namespace metada::backends::gmock