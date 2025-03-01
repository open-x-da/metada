#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "IObsOperator.hpp"
#include "Increment.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Main observation operator class template providing a generic interface
 *
 * This class template provides a static interface for observation operators
 * using a backend specified by the ObsOperatorBackend template parameter. The
 * backend must implement the IObsOperator interface.
 *
 * The observation operator is responsible for:
 * - Mapping model state variables to observation space (H operator)
 * - Computing tangent linear approximations (H' operator)
 * - Computing adjoint operations (H'* operator)
 * - Managing observation errors and uncertainties
 *
 * @tparam ObsOperatorBackend The observation operator backend type
 * @tparam StateType The state type used by this operator
 * @tparam ObsType The observation type used by this operator
 */
template <typename ObsOperatorBackend, typename StateType, typename ObsType>
class ObsOperator {
 private:
  ObsOperatorBackend& backend_;  ///< Reference to backend implementation

 public:
  /** @brief Default constructor - deleted to prevent usage without backend */
  ObsOperator() = delete;

  /** @brief Constructor with backend reference */
  explicit ObsOperator(ObsOperatorBackend& backend) : backend_(backend) {}

  /** @brief Get direct access to the backend instance */
  ObsOperatorBackend& backend() { return backend_; }

  /** @brief Get const access to the backend instance */
  const ObsOperatorBackend& backend() const { return backend_; }

  // Core operations with fluent interface
  ObsOperator& initialize() {
    backend_.initialize();
    return *this;
  }

  ObsOperator& finalize() {
    backend_.finalize();
    return *this;
  }

  // Forward operator: model state -> observation space
  void apply(const State<StateType>& state,
             Observation<ObsType>& observation) const {
    if (!isInitialized()) {
      throw std::runtime_error("ObsOperator not initialized");
    }

    if (!state.isInitialized()) {
      throw std::runtime_error("State not initialized");
    }

    backend_.apply(state.backend(), observation.backend());
  }

  // Tangent linear operator: increment -> observation space
  void applyTangentLinear(const Increment<StateType>& increment,
                          Observation<ObsType>& observation) const {
    if (!isInitialized()) {
      throw std::runtime_error("ObsOperator not initialized");
    }

    if (!increment.isInitialized()) {
      throw std::runtime_error("Increment not initialized");
    }

    backend_.applyTangentLinear(increment.backend(), observation.backend());
  }

  // Adjoint operator: observation -> increment space
  void applyAdjoint(const Observation<ObsType>& observation,
                    Increment<StateType>& increment) const {
    if (!isInitialized()) {
      throw std::runtime_error("ObsOperator not initialized");
    }

    if (!observation.isValid()) {
      throw std::runtime_error("Observation not valid");
    }

    backend_.applyAdjoint(observation.backend(), increment.backend());
  }

  // Error handling with fluent interface
  ObsOperator& setObservationError(const Observation<ObsType>& obs) {
    if (!obs.isValid()) {
      throw std::runtime_error("Cannot set error for invalid observation");
    }

    backend_.setObservationError(obs.backend());
    return *this;
  }

  double getObservationError(const Observation<ObsType>& obs) const {
    if (!obs.isValid()) {
      throw std::runtime_error("Cannot get error for invalid observation");
    }

    return backend_.getObservationError(obs.backend());
  }

  // Configuration with fluent interface
  ObsOperator& setParameter(const std::string& name, double value) {
    backend_.setParameter(name, value);
    return *this;
  }

  double getParameter(const std::string& name) const {
    return backend_.getParameter(name);
  }

  // Required variables
  const std::vector<std::string>& getRequiredStateVariables() const {
    return backend_.getRequiredStateVariables();
  }

  const std::vector<std::string>& getRequiredObsVariables() const {
    return backend_.getRequiredObsVariables();
  }

  bool isInitialized() const { return backend_.isInitialized(); }

  // Validate compatibility between state and observation
  bool validateCompatibility(const State<StateType>& state,
                             const Observation<ObsType>& obs) const {
    // Check if state contains all required variables
    for (const auto& var : getRequiredStateVariables()) {
      if (!state.hasVariable(var)) {
        return false;
      }
    }

    // Check if observation contains all required variables
    for (const auto& var : getRequiredObsVariables()) {
      if (!obs.hasAttribute("variable") ||
          obs.getAttribute("variable") != var) {
        return false;
      }
    }

    return true;
  }
};

}  // namespace metada::framework