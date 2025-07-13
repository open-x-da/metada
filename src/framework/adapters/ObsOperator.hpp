#pragma once

#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "ObsOperatorConcepts.hpp"
#include "ObservationConcepts.hpp"
#include "StateConcepts.hpp"

namespace metada::framework {

// Forward declaration
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

template <typename BackendTag>
  requires ObservationBackendType<BackendTag>
class Observation;

template <typename BackendTag>
  requires StateBackendType<BackendTag>
class State;

/**
 * @brief Adapter class for observation operators in data assimilation systems
 *
 * @details This class provides a type-safe interface for observation operators
 * that map between model state space and observation space. It wraps a backend
 * implementation that satisfies the ObsOperatorBackendType concept, providing
 * consistent access patterns across different backend implementations.
 *
 * The ObsOperator class inherits from NonCopyable to prevent unintended
 * copying, but supports move semantics for efficient resource management.
 *
 * Features:
 * - Forward observation operator (H): Maps state to observation space
 * - Type-safe data access through templates
 * - Delegation to backend implementation
 * - Comprehensive error handling
 * - Move semantics for efficient resource management
 *
 * @tparam BackendTag The backend tag type that must satisfy
 * ObsOperatorBackendType concept
 *
 * @see ObsOperatorBackendType
 * @see NonCopyable
 */
template <typename BackendTag>
  requires ObsOperatorBackendType<BackendTag>
class ObsOperator : public NonCopyable {
 public:
  /** @brief Type alias for the backend implementation type */
  using ObsOperatorBackend =
      typename traits::BackendTraits<BackendTag>::ObsOperatorBackend;

  /** @brief Default constructor is deleted since we need a backend */
  ObsOperator() = delete;

  /**
   * @brief Constructor that initializes observation operator with configuration
   *
   * @details Creates and initializes the observation operator backend using the
   * provided configuration object. The backend handles its own initialization
   * in its constructor.
   *
   * @tparam ConfigBackend The configuration backend type
   * @param config Configuration object containing initialization parameters
   */
  template <typename ConfigBackend>
  explicit ObsOperator(const Config<ConfigBackend>& config)
      : backend_(config.backend()) {
    logger_.Info() << "ObsOperator constructed";
  }

  /**
   * @brief Move constructor
   *
   * @details Transfers ownership of the backend from another observation
   * operator.
   *
   * @param other The observation operator to move from
   */
  ObsOperator(ObsOperator&& other) noexcept
      : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment operator
   *
   * @details Transfers ownership of the backend from another observation
   * operator.
   *
   * @param other The observation operator to move from
   * @return Reference to this observation operator after assignment
   */
  ObsOperator& operator=(ObsOperator&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Get const reference to the backend implementation
   *
   * @return Const reference to the backend implementation
   */
  const ObsOperatorBackend& backend() const { return backend_; }

  /**
   * @brief Check if the observation operator is initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return backend_.isInitialized(); }

  /**
   * @brief Apply forward observation operator: H(x)
   *
   * @details Maps state to observation space (y = H(x)) by transforming model
   * state variables into simulated observations that can be compared with
   * actual observations.
   *
   * @tparam StateBackend The state backend type
   * @tparam ObsBackend The observation backend type
   * @param state Model state to transform
   * @param obs Output observation to store the result
   * @throws std::runtime_error If the observation operator is not initialized
   */
  template <typename StateBackend, typename ObsBackend>
  std::vector<double> apply(const State<StateBackend>& state,
                            const Observation<ObsBackend>& obs) const {
    logger_.Debug() << "Applying observation operator";

    return backend_.apply(state.backend(), obs.backend());
  }

  /**
   * @brief Get the state variables required by this observation operator
   *
   * @return Const reference to vector of required state variable names
   */
  const std::vector<std::string>& getRequiredStateVars() const {
    return backend_.getRequiredStateVars();
  }

  /**
   * @brief Get the observation variables required by this observation operator
   *
   * @return Const reference to vector of required observation variable names
   */
  const std::vector<std::string>& getRequiredObsVars() const {
    return backend_.getRequiredObsVars();
  }

  /**
   * @brief Apply tangent linear observation operator: H dx
   *
   * @details Applies the tangent linear of the observation operator to a state
   * increment. This is used in variational data assimilation for computing the
   * linearized observation operator around a reference trajectory.
   *
   * @tparam StateBackend The state backend type
   * @tparam ObsBackend The observation backend type
   * @param state_increment State increment to transform
   * @param reference_state Reference state around which to linearize
   * @param obs Reference observations for context
   * @return Vector containing the transformed increment in observation space
   * @throws std::runtime_error If the observation operator is not initialized
   */
  template <typename StateBackend, typename ObsBackend>
  std::vector<double> applyTangentLinear(
      const Increment<State<StateBackend>>& state_increment,
      const State<StateBackend>& reference_state,
      const Observation<ObsBackend>& obs) const {
    logger_.Debug() << "Applying tangent linear observation operator";

    return backend_.applyTangentLinear(state_increment.entity().backend(),
                                       reference_state.backend(),
                                       obs.backend());
  }

  /**
   * @brief Apply adjoint observation operator: H^T delta_y
   *
   * @details Applies the adjoint of the observation operator to an observation
   * increment. This maps from observation space back to state space and is
   * essential for computing gradients in variational data assimilation.
   *
   * @param obs_increment Observation space increment
   * @param reference_state Reference state around which the adjoint is computed
   * @return State increment containing the adjoint transformation result
   * @throws std::runtime_error If the observation operator is not initialized
   */
  Increment<BackendTag> applyAdjoint(
      const std::vector<double>& obs_increment,
      const State<BackendTag>& reference_state) const;

  /**
   * @brief Check if tangent linear and adjoint operators are available
   *
   * @return True if tangent linear/adjoint operators are supported
   */
  bool supportsLinearization() const {
    return backend_.supportsLinearization();
  }

  /**
   * @brief Check if the observation operator is linear
   *
   * @details Linear observation operators have the property that H(x+dx) = H(x)
   * + H(dx). For linear operators, the tangent linear and forward operators are
   * identical.
   *
   * @return True if the observation operator is linear
   */
  bool isLinear() const { return backend_.isLinear(); }

 private:
  ObsOperatorBackend backend_; /**< Backend implementation */
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

// Implementation of applyAdjoint for ObsOperator
template <typename BackendTag>
  requires ObsOperatorBackendType<BackendTag>
Increment<BackendTag> ObsOperator<BackendTag>::applyAdjoint(
    const std::vector<double>& obs_increment,
    const State<BackendTag>& reference_state) const {
  // Create an increment from the reference state
  auto increment = Increment<BackendTag>::createFromEntity(reference_state);
  increment.zero();
  backend_.applyAdjoint(obs_increment, reference_state.backend(),
                        increment.state().backend());
  return increment;
}

}  // namespace metada::framework