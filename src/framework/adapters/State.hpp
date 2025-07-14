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
 * - Enable core state operations (zero, dot, norm, etc.)
 * - Implement proper move semantics (non-copyable)
 * - Ensure type-safe data access
 * - Support arithmetic operations (+, -, *, etc.)
 * - Allow state information queries
 * - Support increment creation and application
 *
 * @see Config
 * @see Increment
 * @see StateBackendType
 */

#pragma once

#include <concepts>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "Geometry.hpp"
#include "Location.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "StateConcepts.hpp"

namespace metada::framework {

/**
 * @brief Forward declaration of Config class
 */
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Forward declaration of Increment class
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class Increment;

/**
 * @brief Main state class template providing a generic interface to state
 * implementations
 *
 * @details
 * This class template wraps a state backend implementation and provides a
 * type-safe interface for all state operations. It delegates operations to
 * the backend while adding:
 * - Type safety through templates
 * - Proper move semantics (non-copyable)
 * - Consistent interface across different backends
 * - Support for arithmetic operations
 * - Integration with increment operations
 *
 * The backend tag must satisfy the StateBackendType concept, which ensures
 * it provides valid backend implementation types through BackendTraits.
 *
 * Example usage:
 * @code
 * Config<T> config;
 * Geometry<T> geometry(config);
 * State<Backend> state(config, geometry);
 * state.zero();
 * auto& data = state.getData<double>();
 * @endcode
 *
 * @tparam BackendTag The tag type that defines the state backend through
 * BackendTraits
 *
 * @see StateBackendType
 * @see Config
 * @see Increment
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class State : private NonCopyable {
 public:
  using StateBackend = typename traits::BackendTraits<BackendTag>::StateBackend;
  using Geometry = Geometry<BackendTag>;

  // Constructors
  State() = delete;  // Disable default constructor since we need a backend

  // Destructor
  ~State() = default;

  /**
   * @brief Constructor that initializes state with configuration and geometry
   *
   * @param[in] config Configuration object containing initialization parameters
   * @param[in] geometry Geometry object providing grid information
   * @throws std::runtime_error If backend initialization fails
   */
  State(const Config<BackendTag>& config, const Geometry& geometry)
      : backend_(config.backend(), geometry.backend()),
        geometry_(&geometry),
        initialized_(true) {
    logger_.Info() << "State constructed";
  }

  /**
   * @brief Move constructor
   * @param other State instance to move from
   */
  State(State&& other) noexcept
      : backend_(std::move(other.backend_)),
        geometry_(other.geometry_),
        initialized_(other.initialized_) {
    // Reset the moved-from object's state
    other.geometry_ = nullptr;
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
      // Transfer geometry reference
      geometry_ = other.geometry_;
      // Transfer initialization state
      initialized_ = other.initialized_;
      // Reset the moved-from object's state
      other.geometry_ = nullptr;
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
  /**
   * @brief Check if the state is initialized
   * @return true if initialized, false otherwise
   */
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

  /**
   * @brief Get pointer to underlying data array
   * @tparam T Type of data to access
   * @return Pointer to data array of type T
   */
  template <typename T>
  T* getDataPtr() {
    return static_cast<T*>(backend_.getData());
  }

  /**
   * @brief Get const pointer to underlying data array
   * @tparam T Type of data to access
   * @return Const pointer to data array of type T
   */
  template <typename T>
  const T* getDataPtr() const {
    return static_cast<const T*>(backend_.getData());
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
   * @brief Get the total size of the state vector
   * @return Total number of elements in the state vector
   */
  size_t size() const { return backend_.size(); }

  /**
   * @brief Save state data to file
   * @param filename Path to save state file
   */
  void saveToFile(const std::string& filename) const {
    backend_.saveToFile(filename);
  }

  // Arithmetic operations
  /**
   * @brief Addition operator
   * @param other State to add
   * @return New state containing the sum
   */
  State operator+(const State& other) const {
    State result(std::move(*backend_.clone()));  // Create copy of this state
    result.backend_.add(other.backend_);
    return result;  // NRVO will handle this
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
   * @brief Addition assignment operator for increment
   * @param increment Increment to add
   * @return Reference to this state
   */
  State& operator+=(const Increment<BackendTag>& increment) {
    backend_.add(increment.state().backend());
    return *this;
  }

  /**
   * @brief Subtraction operator
   * @param other State to subtract
   * @return New state containing the difference
   */
  State operator-(const State& other) const {
    State result(std::move(*backend_.clone()));  // Create copy of this state
    result.backend_.subtract(other.backend_);
    return result;  // NRVO will handle this
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
    return result;  // NRVO will handle this
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
   * @brief Get direct access to the backend instance
   * @return Reference to backend implementation
   */
  StateBackend& backend() { return backend_; }

  /**
   * @brief Get const access to the backend instance
   * @return Const reference to backend implementation
   */
  const StateBackend& backend() const { return backend_; }

  /**
   * @brief Get the associated geometry (if any)
   * @return Pointer to the geometry instance, or nullptr if not set
   */
  const Geometry* geometry() const { return geometry_; }

  double& at(const Location& loc) { return backend_.at(loc); }
  const double& at(const Location& loc) const { return backend_.at(loc); }

 private:
  /**
   * @brief Private constructor to be used internally by clone
   * @param backend Backend instance to move from
   */
  explicit State(StateBackend&& backend)
      : backend_(std::move(backend)), geometry_(nullptr), initialized_(true) {}

  StateBackend backend_;
  const Geometry* geometry_ = nullptr;
  bool initialized_{false};
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

// Output operator
template <typename BackendTag>
  requires StateBackendType<BackendTag>
inline std::ostream& operator<<(std::ostream& os,
                                const State<BackendTag>& state) {
  return os << state.backend();
}

}  // namespace metada::framework