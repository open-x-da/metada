/**
 * @file State.hpp
 * @brief Template class providing a generic interface to state implementations
 * @ingroup repr
 * @author Metada Framework Team
 *
 * @details
 * This header provides a template class that wraps state backend
 * implementations and provides a unified interface for state operations. The
 * State class delegates operations to the backend while providing type
 * safety and a consistent API.
 *
 * The State class template is designed to:
 * - Provide a generic interface to different state backend implementations
 * - Support initialization from configuration objects
 * - Enable core state operations (reset, validate, etc)
 * - Implement proper copy/move semantics
 * - Ensure type-safe data access
 * - Support metadata management
 * - Allow state information queries
 *
 * @see IState
 * @see Config
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "utils/NonCopyable.hpp"

namespace metada::framework {

// Forward declaration
template <typename T>
class Config;

/**
 * @brief Main state class template providing a generic interface to state
 * implementations
 *
 * @details
 * This class template wraps a state backend implementation and provides a
 * type-safe interface for all state operations. It delegates operations to
 * the backend while adding:
 * - Type safety through templates
 * - Proper copy/move semantics
 * - Consistent interface across different backends
 *
 * The backend must implement the IState interface to provide the core
 * functionality, while this wrapper adds type safety and convenience.
 *
 * Example usage:
 * @code
 * Config<T> config;
 * State<Backend> state(config);
 * state.reset();
 * auto& data = state.getData<double>();
 * @endcode
 *
 * @tparam Backend The state backend type that implements IState interface
 *
 * @see IState
 */
template <typename Backend>
class State : private NonCopyable {
 public:
  // Constructors
  State() = delete;  // Disable default constructor since we need a backend

  // Destructor
  ~State() = default;

  /**
   * @brief Constructor that initializes state with configuration
   *
   * @tparam T The configuration backend type
   * @param[in] config Configuration object containing initialization parameters
   * @throws std::runtime_error If backend initialization fails
   */
  template <typename T>
  explicit State(const Config<T>& config)
      : backend_(config.backend()), initialized_(true) {}

  /**
   * @brief Move constructor
   * @param other State instance to move from
   */
  State(State&& other) noexcept
      : backend_(std::move(other.backend_)), initialized_(other.initialized_) {
    // Reset the moved-from object's state
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   * @param other State instance to move from
   * @return Reference to this instance
   */
  State& operator=(State&& other) noexcept {
    if (this != &other) {
      // Move the backend state
      backend_ = std::move(other.backend_);

      // Transfer initialization state
      initialized_ = other.initialized_;

      // Reset the moved-from object's state
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Clone the state
   * @return A new state with the same configuration
   */
  State clone() const { return State(std::move(*backend_.clone())); }

  // Core state operations
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Set all values to zero
   * @return Reference to this state
   */
  State& zero() {
    backend_.zero();
    return *this;
  }

  /**
   * @brief Calculate the dot product with another state
   * @param other State to calculate dot product with
   * @return Resulting dot product value
   */
  double dot(const State& other) const { return backend_.dot(other.backend_); }

  /**
   * @brief Calculate the L2 norm of this state
   * @return L2 norm value
   */
  double norm() const { return backend_.norm(); }

  /**
   * @brief Equality operator
   * @param other State instance to compare with
   * @return true if states are equal, false otherwise
   */
  bool operator==(const State& other) const {
    return backend_.equals(other.backend_) &&
           initialized_ == other.initialized_;
  }

  /**
   * @brief Inequality operator
   * @param other State instance to compare with
   * @return true if states are not equal, false otherwise
   */
  bool operator!=(const State& other) const { return !(*this == other); }

  // Data access
  /**
   * @brief Get typed access to underlying data
   * @tparam T Type of data to access
   * @return Reference to data of type T
   */
  template <typename T>
  T& getData() {
    return *static_cast<T*>(backend_.getData());
  }

  /**
   * @brief Get const typed access to underlying data
   * @tparam T Type of data to access
   * @return Const reference to data of type T
   */
  template <typename T>
  const T& getData() const {
    return *static_cast<const T*>(backend_.getData());
  }

  // State information
  /**
   * @brief Get names of state variables
   * @return Const reference to vector of variable names
   */
  const std::vector<std::string>& getVariableNames() const {
    return backend_.getVariableNames();
  }

  /**
   * @brief Check if the state contains a specific variable
   * @param name Name of the variable to check
   * @return true if the variable exists in the state, false otherwise
   */
  bool hasVariable(const std::string& name) const {
    const auto& variables = getVariableNames();
    return std::find(variables.begin(), variables.end(), name) !=
           variables.end();
  }

  /**
   * @brief Get dimensions of state space
   * @return Const reference to vector of dimension sizes
   */
  const std::vector<size_t>& getDimensions(const std::string& name) const {
    return backend_.getDimensions(name);
  }

  /**
   * @brief Get direct access to the backend instance
   * @return Reference to backend implementation
   */
  Backend& backend() { return backend_; }

  /**
   * @brief Get const access to the backend instance
   * @return Const reference to backend implementation
   */
  const Backend& backend() const { return backend_; }

  /**
   * @brief Addition operator
   * @param other State to add
   * @return New state containing the sum
   */
  State operator+(const State& other) const {
    State result(std::move(*backend_.clone()));  // Create copy of this state
    result.backend_.add(other.backend_);
    return std::move(result);  // Explicitly mark for move
  }

  /**
   * @brief Subtraction operator
   * @param other State to subtract
   * @return New state containing the difference
   */
  State operator-(const State& other) const {
    State result(std::move(*backend_.clone()));  // Create copy of this state
    result.backend_.subtract(other.backend_);
    return std::move(result);  // Explicitly mark for move
  }

  /**
   * @brief Addition assignment operator
   * @param other State to add
   * @return Reference to this state
   */
  State& operator+=(const State& other) {
    backend_.add(other.backend_);
    return *this;
  }

  /**
   * @brief Subtraction assignment operator
   * @param other State to subtract
   * @return Reference to this state
   */
  State& operator-=(const State& other) {
    backend_.subtract(other.backend_);
    return *this;
  }

  /**
   * @brief Multiplication operator
   * @param scalar Value to multiply by
   * @return New state containing the product
   */
  State operator*(double scalar) const {
    State result(std::move(*backend_.clone()));  // Create copy of this state
    result.backend_.multiply(scalar);
    return std::move(result);  // Explicitly mark for move
  }

  /**
   * @brief Non-member multiplication operator (scalar * State) (scalar on
   * left)
   * @param scalar Value to multiply by
   * @param state State to multiply
   * @return New state containing the product
   */
  friend State operator*(double scalar, const State& state) {
    return state * scalar;  // Reuse the other operator
  }

  /**
   * @brief Multiplication assignment operator
   * @param scalar Value to multiply by
   * @return Reference to this state
   */
  State& operator*=(double scalar) {
    backend_.multiply(scalar);
    return *this;
  }

  /**
   * @brief Create an increment representing the difference between this state
   * and another
   *
   * @tparam IncrementType The increment type to create
   * @param other The state to subtract from this one
   * @return An increment representing (this - other)
   */
  template <typename IncrementType>
  IncrementType createIncrementTo(const State& other) const;

  /**
   * @brief Apply an increment to this state
   *
   * @tparam IncrementType The increment type to apply
   * @param increment The increment to apply
   * @return Reference to this state after applying the increment
   */
  template <typename IncrementType>
  State& applyIncrement(const IncrementType& increment);

 private:
  // Private constructor to be used internally by clone
  explicit State(Backend&& backend)
      : backend_(std::move(backend)), initialized_(true) {}

  Backend backend_;          ///< Instance of the state backend
  bool initialized_{false};  ///< Initialization flag
};

}  // namespace metada::framework

// Include Increment.hpp after State class definition to resolve circular
// dependency
#include "Increment.hpp"

namespace metada::framework {

// Implementation of methods that depend on Increment
template <typename Backend>
template <typename IncrementType>
IncrementType State<Backend>::createIncrementTo(const State& other) const {
  // Use the factory method in Increment
  return IncrementType::createFromDifference(*this, other);
}

template <typename Backend>
template <typename IncrementType>
State<Backend>& State<Backend>::applyIncrement(const IncrementType& increment) {
  // Use the applyTo method in Increment
  increment.applyTo(*this);
  return *this;
}

}  // namespace metada::framework