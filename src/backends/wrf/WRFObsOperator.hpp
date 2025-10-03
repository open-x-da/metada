#pragma once

#include <memory>
#include <string>
#include <vector>

#include "WRFDAObsOperator.hpp"

namespace metada::backends::wrf {

/**
 * @brief WRF observation operator backend that bridges framework adapters to WRFDA implementation
 *
 * This class serves as the WRF observation operator backend, providing a clean interface
 * between the framework adapter and the WRFDA-specific implementation. It follows the
 * same pattern as other WRF backend classes (WRFObservation, WRFState, WRFGeometry).
 *
 * Key responsibilities:
 * 1. Implement framework ObsOperatorBackendImpl concept requirements
 * 2. Delegate operations to WRFDAObsOperator implementation
 * 3. Provide WRF-specific configuration and initialization
 * 4. Handle linearization and adjoint operations
 * 5. Bridge between framework types and WRFDA types
 *
 * @tparam StateBackend The WRF state backend type
 * @tparam ObsBackend The WRF observation backend type
 */
template <typename StateBackend, typename ObsBackend>
class WRFObsOperator {
 public:
  // Delete default constructor and copying (required by framework concept)
  WRFObsOperator() = delete;
  WRFObsOperator(const WRFObsOperator&) = delete;
  WRFObsOperator& operator=(const WRFObsOperator&) = delete;

  /**
   * @brief Constructor with configuration
   * @param config Configuration object containing operator settings
   *
   * @details Initializes the WRF observation operator with the provided configuration.
   * The configuration is passed through to the underlying WRFDA implementation.
   */
  template <typename ConfigBackend>
  explicit WRFObsOperator(const ConfigBackend& config)
      : wrfda_impl_(config), initialized_(false) {
    // WRFDAObsOperator handles its own initialization in constructor
    initialized_ = wrfda_impl_.isInitialized();
  }

  /**
   * @brief Move constructor
   * @param other WRFObsOperator to move from
   */
  WRFObsOperator(WRFObsOperator&& other) noexcept
      : wrfda_impl_(std::move(other.wrfda_impl_)),
        initialized_(other.initialized_) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   * @param other WRFObsOperator to move from
   * @return Reference to this object
   */
  WRFObsOperator& operator=(WRFObsOperator&& other) noexcept {
    if (this != &other) {
      wrfda_impl_ = std::move(other.wrfda_impl_);
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Initialize the observation operator
   * @param config Configuration object
   *
   * @details Initializes the underlying WRFDA implementation. This method
   * is required by the framework concept but may be redundant since
   * WRFDAObsOperator initializes itself in its constructor.
   */
  template <typename ConfigBackend>
  void initialize(const ConfigBackend& config) {
    wrfda_impl_.initialize(config);
    initialized_ = wrfda_impl_.isInitialized();
  }

  /**
   * @brief Check if the operator is initialized
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const {
    return initialized_ && wrfda_impl_.isInitialized();
  }

  /**
   * @brief Apply the forward observation operator: H(x)
   * @param state Model state to transform
   * @param obs Observations defining observation locations
   * @return Vector of simulated observations
   *
   * @details Maps model state to observation space by applying the WRFDA
   * observation operator. This computes H(x) where H is the observation
   * operator and x is the model state.
   */
  std::vector<double> apply(const StateBackend& state,
                            const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }
    return wrfda_impl_.apply(state, obs);
  }

  /**
   * @brief Apply the tangent linear observation operator: H'(δx)
   * @param state_increment State increment to transform
   * @param reference_state Reference state around which to linearize
   * @param obs Observations defining observation locations
   * @return Vector of observation increments
   *
   * @details Applies the tangent linear of the observation operator to a
   * state increment. This is used in variational data assimilation for
   * computing the linearized observation operator.
   */
  std::vector<double> applyTangentLinear(
      const StateBackend& state_increment,
      const StateBackend& reference_state,
      const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }
    return wrfda_impl_.applyTangentLinear(state_increment, reference_state, obs);
  }

  /**
   * @brief Apply the adjoint observation operator: H^T(δy)
   * @param obs_increment Observation space increment
   * @param reference_state Reference state around which adjoint is computed
   * @param result_state State to accumulate adjoint result in
   * @param obs Observations defining observation locations
   *
   * @details Applies the adjoint of the observation operator to map from
   * observation space back to state space. This is essential for computing
   * gradients in variational data assimilation.
   */
  void applyAdjoint(const std::vector<double>& obs_increment,
                    const StateBackend& reference_state,
                    StateBackend& result_state,
                    const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFObsOperator not initialized");
    }
    wrfda_impl_.applyAdjoint(obs_increment, reference_state, result_state, obs);
  }

  /**
   * @brief Get required state variables
   * @return Reference to vector of required state variable names
   */
  const std::vector<std::string>& getRequiredStateVars() const {
    return wrfda_impl_.getRequiredStateVars();
  }

  /**
   * @brief Get required observation variables
   * @return Reference to vector of required observation variable names
   */
  const std::vector<std::string>& getRequiredObsVars() const {
    return wrfda_impl_.getRequiredObsVars();
  }

  /**
   * @brief Check if linearization is supported
   * @return True if tangent linear and adjoint operators are available
   */
  bool supportsLinearization() const {
    return wrfda_impl_.supportsLinearization();
  }

  /**
   * @brief Check if the operator is linear
   * @return True if the observation operator is linear
   *
   * @details Linear operators have the property that H(x+dx) = H(x) + H(dx).
   * For linear operators, the tangent linear and forward operators are identical.
   */
  bool isLinear() const {
    return wrfda_impl_.isLinear();
  }

  /**
   * @brief Get access to underlying WRFDA implementation
   * @return Reference to WRFDAObsOperator implementation
   *
   * @details Provides access to the underlying WRFDA implementation for
   * advanced use cases. This should be used sparingly and only when
   * WRF-specific functionality is needed.
   */
  const common::obsoperator::WRFDAObsOperator<StateBackend, ObsBackend>&
  wrfdaImpl() const {
    return wrfda_impl_;
  }

 private:
  /// Underlying WRFDA observation operator implementation
  common::obsoperator::WRFDAObsOperator<StateBackend, ObsBackend> wrfda_impl_;
  
  /// Initialization status flag
  bool initialized_;
};

}  // namespace metada::backends::wrf
