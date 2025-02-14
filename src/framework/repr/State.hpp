/**
 * @file State.hpp
 * @brief Template class providing a generic interface to state implementations
 * @ingroup repr
 *
 * @details
 * This header provides a template class that wraps state backend
 * implementations and provides a unified interface for state operations. The
 * State class delegates operations to the backend while providing type
 * safety and a consistent API.
 *
 * Key features:
 * - Generic interface to different state backend implementations
 * - Support for initialization from configuration
 * - Core state operations (reset, validate, etc)
 * - Copy/move semantics
 * - Data access with type safety
 * - Metadata management
 * - State information queries
 */

#ifndef METADA_SRC_FRAMEWORK_REPR_STATE_HPP_
#define METADA_SRC_FRAMEWORK_REPR_STATE_HPP_

#include <memory>
#include <string>
#include <vector>

#include "Config.hpp"
#include "IState.hpp"

using metada::framework::tools::config::Config;
namespace metada {
namespace framework {
namespace repr {

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
 * @tparam StateBackend The state backend type that implements IState interface
 */
template <typename StateBackend>
class State {
 private:
  StateBackend& backend_;    ///< Instance of the state backend
  bool initialized_{false};  ///< Initialization flag

 public:
  // Constructors
  State() = delete;  // Disable default constructor since we need a backend

  /**
   * @brief Basic constructor with backend only
   * @param backend Reference to backend implementation
   */
  explicit State(StateBackend& backend)
      : backend_(backend), initialized_(true) {}

  /**
   * @brief Copy constructor
   * @param other State instance to copy from
   */
  explicit State(const State& other)
      : backend_(other.backend_), initialized_(other.initialized_) {}

  /**
   * @brief Constructor that initializes state with backend and configuration
   *
   * @details
   * This constructor takes both a backend implementation and a configuration
   * object. It initializes the backend with the provided configuration and sets
   * the initialized flag to true upon successful initialization.
   *
   * @tparam T The configuration backend type
   * @param[in] backend Reference to the state backend implementation
   * @param[in] config Configuration object containing initialization parameters
   * @throws std::runtime_error If backend initialization fails
   */
  template <typename T>
  explicit State(StateBackend& backend, const Config<T>& config)
      : backend_(backend), initialized_(false) {
    backend_.initialize(config.backend());
    initialized_ = true;
  }

  // Core state operations
  void reset() { backend_.reset(); }
  void validate() const { backend_.validate(); }
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Move constructor
   * @param other State instance to move from
   */
  State(State&& other) noexcept
      : backend_(other.backend_), initialized_(other.initialized_) {
    // No need for moveFrom since we're just taking the reference
  }

  /**
   * @brief Copy assignment operator
   * @param other State instance to copy from
   * @return Reference to this instance
   */
  State& operator=(const State& other) {
    if (this != &other) {
      backend_.copyFrom(other.backend_);
      initialized_ = other.initialized_;
    }
    return *this;
  }

  /**
   * @brief Move assignment operator
   * @param other State instance to move from
   * @return Reference to this instance
   */
  State& operator=(State&& other) noexcept {
    if (this != &other) {
      backend_.moveFrom(std::move(other.backend_));
      initialized_ = other.initialized_;
    }
    return *this;
  }

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

  // Metadata operations
  /**
   * @brief Set metadata value for given key
   * @param key Metadata key
   * @param value Value to set
   */
  void setMetadata(const std::string& key, const std::string& value) {
    backend_.setMetadata(key, value);
  }

  /**
   * @brief Get metadata value for given key
   * @param key Metadata key
   * @return Value associated with key
   */
  std::string getMetadata(const std::string& key) const {
    return backend_.getMetadata(key);
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
   * @brief Get dimensions of state space
   * @return Const reference to vector of dimension sizes
   */
  const std::vector<size_t>& getDimensions() const {
    return backend_.getDimensions();
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
};

}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_SRC_FRAMEWORK_REPR_STATE_HPP_