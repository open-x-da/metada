#pragma once

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "IObservation.hpp"
#include "utils/NonCopyable.hpp"
#include "utils/config/Config.hpp"

namespace metada::framework {

// Concepts to check backend requirements
template <typename T>
concept HasClone = requires(const T& t) {
  { t.clone() } -> std::convertible_to<std::unique_ptr<T>>;
};

template <typename T>
concept HasGetData = requires(T& t, const T& ct) {
  { t.getData() } -> std::same_as<void*>;
  { ct.getData() } -> std::same_as<const void*>;
};

template <typename T>
concept HasGetVariableNames = requires(const T& t) {
  { t.getVariableNames() } -> std::same_as<const std::vector<std::string>&>;
};

template <typename T>
concept HasGetDimensions = requires(const T& t) {
  {
    t.getDimensions(std::string{})
  } -> std::same_as<const std::vector<size_t>&>;
};

// Combined concept for all backend requirements
template <typename T>
concept ObservationBackend = HasClone<T> && HasGetData<T> &&
                             HasGetVariableNames<T> && HasGetDimensions<T>;

/**
 * @brief Adapter class for observation implementations in data assimilation
 * systems
 *
 * This template class provides a type-safe interface for handling observational
 * data using a backend that implements the IObservation interface.
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
 * @tparam Backend The backend implementation type
 */
template <typename Backend>
  requires ObservationBackend<Backend>
class Observation : private NonCopyable {
 public:
  /** @brief Default constructor is deleted since we need a backend */
  Observation() = delete;

  /** @brief Destructor */
  ~Observation() = default;

  /**
   * @brief Constructor that initializes observation with configuration
   *
   * @tparam T The configuration backend type
   * @param config Configuration object containing initialization parameters
   */
  template <typename T>
  explicit Observation(const Config<T>& config)
      : backend_(config.backend()), initialized_(true) {}

  /**
   * @brief Move constructor
   *
   * Transfers ownership of the backend from another observation
   * and updates initialization state.
   */
  Observation(Observation&& other) noexcept
      : backend_(std::move(other.backend_)), initialized_(other.initialized_) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   *
   * Transfers ownership of the backend from another observation
   * and updates initialization state.
   *
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
   * @return A new observation with the same data and configuration
   */
  Observation clone() const {
    return Observation(std::move(*backend_.clone()));
  }

  /**
   * @brief Check if the observation is properly initialized
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply quality control procedures to the observation
   */
  void applyQC() { backend_.applyQC(); }

  /**
   * @brief Load observation data from a file
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
   * @param filename Path where observation data will be saved
   */
  void saveToFile(const std::string& filename) const {
    backend_.saveToFile(filename);
  }

  /**
   * @brief Equality operator
   *
   * @param other The observation to compare with
   * @return bool True if the observations are equal, false otherwise
   */
  bool operator==(const Observation& other) const {
    return backend_.equals(other.backend_) &&
           initialized_ == other.initialized_;
  }

  /**
   * @brief Inequality operator
   *
   * @param other The observation to compare with
   * @return bool True if the observations are not equal, false otherwise
   */
  bool operator!=(const Observation& other) const { return !(*this == other); }

  /**
   * @brief Get typed access to the underlying data
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
   * @return Const reference to vector of variable names
   */
  const std::vector<std::string>& getVariableNames() const {
    return backend_.getVariableNames();
  }

  /**
   * @brief Check if a specific variable exists in this observation
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
   * @brief Get the dimensions of this observation
   *
   * @return Const reference to vector of dimension sizes
   */
  const std::vector<size_t>& getDimensions(const std::string& name) const {
    return backend_.getDimensions(name);
  }

  /**
   * @brief Addition operator
   * @param other The observation to add to this one
   * @return Observation The result of the addition
   */
  Observation operator+(const Observation& other) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.add(other.backend_);
    return result;
  }

  /**
   * @brief Compound addition operator
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
   * @param other The observation to subtract from this one
   * @return Observation The result of the subtraction
   */
  Observation operator-(const Observation& other) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.subtract(other.backend_);
    return result;
  }

  /**
   * @brief Compound subtraction operator
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
   * @param scalar The scalar value to multiply this observation by
   * @return Observation The result of the multiplication
   */
  Observation operator*(double scalar) const {
    Observation result(std::move(*backend_.clone()));
    result.backend_.multiply(scalar);
    return result;
  }

  /**
   * @brief Compound multiplication operator
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
   * @param scalar The scalar value to multiply this observation by
   * @param obs The observation to multiply
   * @return Observation The result of the multiplication
   */
  friend Observation operator*(double scalar, const Observation& obs) {
    return obs * scalar;
  }

  /**
   * @brief Access to backend instance
   * @return Backend& Reference to the backend instance
   */
  Backend& backend() { return backend_; }

  /**
   * @brief Access to backend instance
   * @return const Backend& Reference to the backend instance
   */
  const Backend& backend() const { return backend_; }

 private:
  /**
   * @brief Private constructor to be used internally by clone
   * @param backend Backend instance to move from
   */
  explicit Observation(Backend&& backend)
      : backend_(std::move(backend)), initialized_(true) {}

  Backend backend_;          ///< Backend implementation instance
  bool initialized_{false};  ///< Initialization flag
};

}  // namespace metada::framework