#pragma once

#include <memory>
#include <string>

#include "IModel.hpp"
#include "State.hpp"

namespace metada::framework::repr {

/**
 * @brief Main model class template providing a generic interface to model
 * implementations
 *
 * This class template provides a static interface for model operations using a
 * backend specified by the ModelBackend template parameter. The backend must
 * implement the IModel interface.
 *
 * Key features:
 * - Model lifecycle management
 * - State handling with type safety
 * - Parameter configuration
 * - Model metadata access
 *
 * @tparam ModelBackend The model backend type that implements IModel
 * @tparam StateType The state type used by this model
 */
template <typename ModelBackend, typename StateType>
class Model {
 private:
  ModelBackend backend_;  ///< Instance of the model backend

 public:
  /** @brief Default constructor */
  Model() = default;

  /** @brief Get direct access to the backend instance */
  ModelBackend& backend() { return backend_; }

  /** @brief Get const access to the backend instance */
  const ModelBackend& backend() const { return backend_; }

  // Core model operations
  void initialize() { backend_.initialize(); }
  void step(double dt) { backend_.step(dt); }
  void finalize() { backend_.finalize(); }

  // State management
  void setState(const State<StateType>& state) {
    const StateType& typedState = state.backend();
    backend_.setState(typedState);
  }

  State<StateType> getState() const {
    const IState& backendState = backend_.getState();
    const StateType& typedState = dynamic_cast<const StateType&>(backendState);
    return State<StateType>(typedState);
  }

  // Model configuration
  void setParameter(const std::string& name, double value) {
    backend_.setParameter(name, value);
  }

  double getParameter(const std::string& name) const {
    return backend_.getParameter(name);
  }

  // Model metadata
  std::string getName() const { return backend_.getName(); }
  std::string getVersion() const { return backend_.getVersion(); }

  // Model state
  double getCurrentTime() const { return backend_.getCurrentTime(); }
  bool isInitialized() const { return backend_.isInitialized(); }
};

}  // namespace metada::framework::repr