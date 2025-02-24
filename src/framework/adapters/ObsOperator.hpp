#pragma once

#include <memory>
#include <vector>

#include "IObsOperator.hpp"
#include "Increment.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace metada::framework::repr {

/**
 * @brief Main observation operator class template providing a generic interface
 *
 * This class template provides a static interface for observation operators
 * using a backend specified by the ObsOperatorBackend template parameter. The
 * backend must implement the IObsOperator interface.
 *
 * @tparam ObsOperatorBackend The observation operator backend type
 * @tparam StateType The state type used by this operator
 * @tparam ObsType The observation type used by this operator
 */
template <typename ObsOperatorBackend, typename StateType, typename ObsType>
class ObsOperator {
 private:
  ObsOperatorBackend
      backend_;  ///< Instance of the observation operator backend

 public:
  /** @brief Default constructor */
  ObsOperator() = default;

  /** @brief Get direct access to the backend instance */
  ObsOperatorBackend& backend() { return backend_; }

  /** @brief Get const access to the backend instance */
  const ObsOperatorBackend& backend() const { return backend_; }

  // Core operations
  void initialize() { backend_.initialize(); }
  void finalize() { backend_.finalize(); }

  // Forward operator: model state -> observation space
  void apply(const State<StateType>& state,
             Observation<ObsType>& observation) const {
    backend_.apply(state.backend(), observation.backend());
  }

  // Tangent linear operator: increment -> observation space
  void applyTangentLinear(const Increment<StateType>& increment,
                          Observation<ObsType>& observation) const {
    backend_.applyTangentLinear(increment.backend(), observation.backend());
  }

  // Adjoint operator: observation -> increment space
  void applyAdjoint(const Observation<ObsType>& observation,
                    Increment<StateType>& increment) const {
    backend_.applyAdjoint(observation.backend(), increment.backend());
  }

  // Error handling
  void setObservationError(const Observation<ObsType>& obs) {
    backend_.setObservationError(obs.backend());
  }

  double getObservationError(const Observation<ObsType>& obs) const {
    return backend_.getObservationError(obs.backend());
  }

  // Configuration
  void setParameter(const std::string& name, double value) {
    backend_.setParameter(name, value);
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
};

}  // namespace metada::framework::repr