#pragma once

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "NonCopyable.hpp"
#include "ObservationConcepts.hpp"

namespace metada::framework {

/**
 * @brief Forward declaration of Config class
 *
 * @tparam BackendTag The backend tag type that must satisfy ConfigBackendType
 * concept
 */
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Adapter class for observation implementations in data assimilation
 * systems
 *
 * @details This class provides a type-safe interface for handling observational
 * data in data assimilation workflows. It wraps a backend implementation that
 * satisfies the ObservationBackendType concept, providing consistent access
 * patterns and operations across different backend implementations.
 *
 * The Observation class inherits from NonCopyable to prevent unintended
 * copying, but provides explicit clone functionality for when copies are needed
 * in operations like ensemble generation, observation perturbation, and
 * observation-space calculations in data assimilation algorithms.
 *
 * Features:
 * - Type-safe data access through templates
 * - Delegation to backend implementation
 * - Comprehensive error handling
 * - Arithmetic operations for observation-space calculations
 * - Explicit cloning mechanism for controlled copying
 * - Move semantics for efficient resource management
 *
 * @tparam BackendTag The backend tag type that must satisfy
 * ObservationBackendType concept
 *
 * @see ObservationBackendType
 * @see NonCopyable
 */
template <typename BackendTag>
  requires ObservationBackendType<BackendTag>
class Observation : private NonCopyable {
 public:
  /** @brief Type alias for the backend implementation type */
  using ObservationBackend =
      typename traits::BackendTraits<BackendTag>::ObservationBackend;

  /** @brief Default constructor is deleted since we need a backend */
  Observation() = delete;

  /** @brief Destructor */
  ~Observation() = default;

  /**
   * @brief Constructor that initializes observation with configuration
   *
   * @param config Configuration object containing initialization parameters
   */
  explicit Observation(const Config<BackendTag>& config)
      : backend_(config.backend()), initialized_(true) {}

  /**
   * @brief Move constructor
   *
   * @details Transfers ownership of the backend from another observation
   * and updates initialization state.
   *
   * @param other The observation to move from
   */
  Observation(Observation&& other) noexcept
      : backend_(std::move(other.backend_)), initialized_(other.initialized_) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   *
   * @details Transfers ownership of the backend from another observation
   * and updates initialization state.
   *
   * @param other The observation to move from
   * @return Reference to this observation after assignment
   */
  Observation& operator=(Observation&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Clone the observation
   *
   * @details Creates a deep copy of this observation, including all data and
   * metadata. This is the preferred way to create copies of observations when
   * needed.
   *
   * @return A new observation with the same data and configuration
   */
  Observation clone() const {
    return Observation(std::move(*backend_.clone()));
  }

  /**
   * @brief Check if the observation is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply quality control procedures to the observation
   *
   * @details Delegates to the backend implementation to perform quality control
   * checks and filtering on the observation data.
   */
  void applyQC() { backend_.applyQC(); }

  /**
   * @brief Load observation data from a file
   *
   * @details Delegates to the backend implementation to read observation data
   * from the specified file and initialize this observation object.
   *
   * @param filename Path to the file containing observation data
   */
  void loadFromFile(const std::string& filename) {
    backend_.loadFromFile(filename);
    initialized_ = true;
  }

  /**
   * @brief Save observation data to a file
   *
   * @details Delegates to the backend implementation to write the current
   * observation data to the specified file.
   *
   * @param filename Path where observation data will be saved
   */
  void saveToFile(const std::string& filename) const {
    backend_.saveToFile(filename);
  }

  /**
   * @brief Equality operator
   *
   * @details Compares this observation with another for equality by delegating
   * to the backend implementation and checking initialization status.
   *
   * @param other The observation to compare with
   * @return True if the observations are equal, false otherwise
   */
  bool operator==(const Observation& other) const {
    return backend_.equals(other.backend_) &&
           initialized_ == other.initialized_;
  }

  /**
   * @brief Inequality operator
   *
   * @details Compares this observation with another for inequality by negating
   * the result of the equality operator.
   *
   * @param other The observation to compare with
   * @return True if the observations are not equal, false otherwise
   */
  bool operator!=(const Observation& other) const { return !(*this == other); }

  /**
   * @brief Get typed access to the underlying data
   *
   * @details Provides access to the observation data cast to the specified
   * type. The caller is responsible for ensuring the type is correct.
   *
   * @tparam T The expected data type
   * @return Reference to the data cast to the specified type
   */
  template <typename T>
  T& getData() {
    return *static_cast<T*>(backend_.getData());
  }

  /**
   * @brief Get const typed access to the underlying data
   *
   * @details Provides read-only access to the observation data cast to the
   * specified type. The caller is responsible for ensuring the type is correct.
   *
   * @tparam T The expected data type
   * @return Const reference to the data cast to the specified type
   */
  template <typename T>
  const T& getData() const {
    return *static_cast<const T*>(backend_.getData());
  }

  /**
   * @brief Get the names of variables in this observation
   *
   * @details Retrieves the list of variable names contained in this
   * observation.
   *
   * @return Const reference to vector of variable names
   */
  const std::vector<std::string>& getVariableNames() const {
    return backend_.getVariableNames();
  }

  /**
   * @brief Check if a specific variable exists in this observation
   *
   * @details Searches for the specified variable name in the list of variables
   * contained in this observation.
   *
   * @param name Name of the variable to check
   * @return True if the variable exists, false otherwise
   */
  bool hasVariable(const std::string& name) const {
    const auto& variables = getVariableNames();
    return std::find(variables.begin(), variables.end(), name) !=
           variables.end();
  }

  /**
   * @brief Get the dimensions of a specific variable in this observation
   *
   * @details Retrieves the dimension sizes for the specified variable.
   *
   * @param name Name of the variable to get dimensions for
   * @return Const reference to vector of dimension sizes
   */
  const std::vector<size_t>& getDimensions(const std::string& name) const {
    return backend_.getDimensions(name);
  }

  /**
   * @brief Addition operator
   *
   * @details Creates a new observation that is the sum of this observation and
   * another.
   *
   * @param other The observation to add to this one
   * @return A new observation containing the result of the addition
   */
  Observation operator+(const Observation& other) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.add(other.backend_);
    return result;
  }

  /**
   * @brief Compound addition operator
   *
   * @details Adds another observation to this one in-place.
   *
   * @param other The observation to add to this one
   * @return Reference to this observation after addition
   */
  Observation& operator+=(const Observation& other) {
    backend_.add(other.backend_);
    return *this;
  }

  /**
   * @brief Subtraction operator
   *
   * @details Creates a new observation that is the difference between this
   * observation and another.
   *
   * @param other The observation to subtract from this one
   * @return A new observation containing the result of the subtraction
   */
  Observation operator-(const Observation& other) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.subtract(other.backend_);
    return result;
  }

  /**
   * @brief Compound subtraction operator
   *
   * @details Subtracts another observation from this one in-place.
   *
   * @param other The observation to subtract from this one
   * @return Reference to this observation after subtraction
   */
  Observation& operator-=(const Observation& other) {
    backend_.subtract(other.backend_);
    return *this;
  }

  /**
   * @brief Multiplication operator
   *
   * @details Creates a new observation that is this observation multiplied by a
   * scalar.
   *
   * @param scalar The scalar value to multiply this observation by
   * @return A new observation containing the result of the multiplication
   */
  Observation operator*(double scalar) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.multiply(scalar);
    return result;
  }

  /**
   * @brief Compound multiplication operator
   *
   * @details Multiplies this observation by a scalar in-place.
   *
   * @param scalar The scalar value to multiply this observation by
   * @return Reference to this observation after multiplication
   */
  Observation& operator*=(double scalar) {
    backend_.multiply(scalar);
    return *this;
  }

  /**
   * @brief Multiplication operator with scalar
   *
   * @details Creates a new observation that is the result of multiplying
   * a scalar by an observation. This enables the commutative property
   * of scalar multiplication.
   *
   * @param scalar The scalar value to multiply the observation by
   * @param obs The observation to multiply
   * @return A new observation containing the result of the multiplication
   */
  friend Observation operator*(double scalar, const Observation& obs) {
    return obs * scalar;
  }

  /**
   * @brief Access to backend instance
   *
   * @details Provides direct access to the underlying backend implementation.
   * This should be used with caution as it bypasses the adapter's interface.
   *
   * @return Reference to the backend instance
   */
  ObservationBackend& backend() { return backend_; }

  /**
   * @brief Const access to backend instance
   *
   * @details Provides read-only access to the underlying backend
   * implementation. This should be used with caution as it bypasses the
   * adapter's interface.
   *
   * @return Const reference to the backend instance
   */
  const ObservationBackend& backend() const { return backend_; }

 private:
  /**
   * @brief Private constructor to be used internally by clone
   *
   * @details Creates an observation by taking ownership of an existing backend
   * instance. This is primarily used by the clone method.
   *
   * @param backend Backend instance to move from
   */
  explicit Observation(ObservationBackend&& backend)
      : backend_(std::move(backend)), initialized_(true) {}

  ObservationBackend backend_;  ///< Backend implementation instance
  bool initialized_{false};     ///< Initialization flag
};

}  // namespace metada::framework