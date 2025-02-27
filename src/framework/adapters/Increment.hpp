#pragma once

#include <memory>
#include <string>
#include <vector>

#include "../interfaces/IIncrement.hpp"

namespace metada::framework {

// Forward declaration
template <typename T>
class State;

/**
 * @brief Main increment class template providing a generic interface to
 * increment implementations
 *
 * This class template provides a static interface for handling state
 * increments/perturbations using a backend specified by the IncrementBackend
 * template parameter. The backend must implement the IIncrement interface.
 *
 * @tparam IncrementBackend The increment backend type that implements
 * IIncrement
 */
template <typename IncrementBackend>
class Increment {
 private:
  IncrementBackend backend_;  ///< Instance of the increment backend

 public:
  /** @brief Default constructor */
  Increment() = default;

  /** @brief Get direct access to the backend instance */
  IncrementBackend& backend() { return backend_; }

  /** @brief Get const access to the backend instance */
  const IncrementBackend& backend() const { return backend_; }

  // Core increment operations
  void initialize() { backend_.initialize(); }
  void zero() { backend_.zero(); }
  void scale(double alpha) { backend_.scale(alpha); }

  // Linear algebra operations
  void axpy(double alpha, const Increment& other) {
    backend_.axpy(alpha, other.backend());
  }

  double dot(const Increment& other) const {
    return backend_.dot(other.backend());
  }

  double norm() const { return backend_.norm(); }

  // State operations - forward declare these methods
  template <typename StateType>
  void addToState(State<StateType>& state) const;

  template <typename StateType>
  void differenceFromStates(const State<StateType>& state1,
                            const State<StateType>& state2);

  /**
   * @brief Create an increment from the difference of two states
   *
   * Factory method that creates a new increment representing state1 - state2.
   *
   * @param state1 First state
   * @param state2 Second state
   * @return A new increment containing the difference
   */
  template <typename StateType>
  static Increment createFromDifference(const State<StateType>& state1,
                                        const State<StateType>& state2);

  // Data access
  template <typename T>
  T& getData() {
    return *static_cast<T*>(backend_.getData());
  }

  template <typename T>
  const T& getData() const {
    return *static_cast<const T*>(backend_.getData());
  }

  // Metadata operations
  void setMetadata(const std::string& key, const std::string& value) {
    backend_.setMetadata(key, value);
  }

  std::string getMetadata(const std::string& key) const {
    return backend_.getMetadata(key);
  }

  // Increment information
  const std::vector<size_t>& getDimensions() const {
    return backend_.getDimensions();
  }

  bool isInitialized() const { return backend_.isInitialized(); }
};

}  // namespace metada::framework

// Include State.hpp after Increment class definition
#include "State.hpp"

namespace metada::framework {

// Implementation of methods that depend on State
template <typename IncrementBackend>
template <typename StateType>
void Increment<IncrementBackend>::addToState(State<StateType>& state) const {
  backend_.addToState(state.backend());
}

template <typename IncrementBackend>
template <typename StateType>
void Increment<IncrementBackend>::differenceFromStates(
    const State<StateType>& state1, const State<StateType>& state2) {
  backend_.differenceFromStates(state1.backend(), state2.backend());
}

template <typename IncrementBackend>
template <typename StateType>
Increment<IncrementBackend> Increment<IncrementBackend>::createFromDifference(
    const State<StateType>& state1, const State<StateType>& state2) {
  Increment result;
  result.initialize();
  result.differenceFromStates(state1, state2);
  return result;
}

}  // namespace metada::framework