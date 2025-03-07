#pragma once

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "IObservation.hpp"
#include "utils/config/Config.hpp"

namespace metada::framework {

/**
 * @brief Adapter class for observation implementations in data assimilation
 * systems
 *
 * This template class provides a type-safe interface for handling observational
 * data using a backend that implements the IObservation interface.
 *
 * Unlike the previous design, observations are now copyable to allow for
 * operations like ensemble generation, observation perturbation, and
 * observation-space calculations in data assimilation algorithms.
 *
 * Features:
 * - Type-safe data access through templates
 * - Delegation to backend implementation
 * - Comprehensive error handling
 * - Arithmetic operations for observation-space calculations
 *
 * @tparam Backend The backend implementation type
 */
template <typename Backend>
class Observation {
 private:
  Backend backend_;          ///< Backend implementation instance
  bool initialized_{false};  ///< Initialization flag

 public:
  /** @brief Default constructor is deleted since we need a backend */
  Observation() = delete;

  /**
   * @brief Constructor that initializes observation with configuration
   *
   * @tparam T The configuration backend type
   * @param config Configuration object containing initialization parameters
   */
  template <typename T>
  explicit Observation(const Config<T>& config)
      : backend_(config), initialized_(true) {}

  /**
   * @brief Constructor with backend rvalue reference
   *
   * This constructor takes ownership of the backend instance.
   * Note: The backend must be already initialized with a config.
   */
  explicit Observation(Backend&& backend)
      : backend_(std::move(backend)), initialized_(true) {}

  /**
   * @brief Constructor with backend lvalue reference
   *
   * This constructor is particularly useful for shared_ptr backends or
   * when working with mock objects that can't be moved.
   * Note: The backend must be already initialized with a config.
   */
  explicit Observation(Backend& backend)
      : backend_(backend), initialized_(true) {}

  /**
   * @brief Copy constructor
   */
  Observation(const Observation& other) : backend_(other.backend_.config()) {
    backend_.copyFrom(other.backend_);
    initialized_ = other.initialized_;
  }

  /**
   * @brief Move constructor
   */
  Observation(Observation&& other) noexcept
      : backend_(other.backend_.config()) {
    backend_.moveFrom(std::move(other.backend_));
    initialized_ = other.initialized_;
    other.initialized_ = false;
  }

  /**
   * @brief Copy assignment operator
   */
  Observation& operator=(const Observation& other) {
    if (this != &other) {
      backend_.copyFrom(other.backend_);
      initialized_ = other.initialized_;
    }
    return *this;
  }

  /**
   * @brief Move assignment operator
   */
  Observation& operator=(Observation&& other) noexcept {
    if (this != &other) {
      backend_.moveFrom(std::move(other.backend_));
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /** @brief Access to backend instance */
  Backend& backend() { return backend_; }
  const Backend& backend() const { return backend_; }

  // Lifecycle management
  Observation& reset() {
    backend_.reset();
    return *this;
  }

  void validate() const { backend_.validate(); }

  bool isInitialized() const { return initialized_; }

  // Copy/Move operations
  Observation& copyFrom(const Observation& other) {
    backend_.copyFrom(other.backend_);
    initialized_ = other.initialized_;
    return *this;
  }

  Observation& moveFrom(Observation&& other) {
    backend_.moveFrom(std::move(other.backend_));
    initialized_ = other.initialized_;
    other.initialized_ = false;
    return *this;
  }

  bool equals(const Observation& other) const {
    return backend_.equals(other.backend_);
  }

  /**
   * @brief Addition assignment operator
   *
   * @param other The observation to add to this one
   * @return Observation& Reference to this observation after addition
   * @throws std::runtime_error If either observation is invalid
   */
  Observation& operator+=(const Observation& other) {
    backend_.add(other.backend_);
    return *this;
  }

  /**
   * @brief Subtraction assignment operator
   *
   * @param other The observation to subtract from this one
   * @return Observation& Reference to this observation after subtraction
   * @throws std::runtime_error If either observation is invalid
   */
  Observation& operator-=(const Observation& other) {
    backend_.subtract(other.backend_);
    return *this;
  }

  /**
   * @brief Multiplication assignment operator
   *
   * @param scalar The scalar value to multiply this observation by
   * @return Observation& Reference to this observation after multiplication
   * @throws std::runtime_error If the observation is invalid
   */
  Observation& operator*=(double scalar) {
    backend_.multiply(scalar);
    return *this;
  }

  /**
   * @brief Equality operator
   *
   * @param other The observation to compare with
   * @return bool True if the observations are equal, false otherwise
   */
  bool operator==(const Observation& other) const { return equals(other); }

  /**
   * @brief Inequality operator
   *
   * @param other The observation to compare with
   * @return bool True if the observations are not equal, false otherwise
   */
  bool operator!=(const Observation& other) const { return !equals(other); }

  // Type-safe data access
  template <typename T>
  T& getData() {
    return *static_cast<T*>(backend_.getData());
  }

  template <typename T>
  const T& getData() const {
    return *static_cast<const T*>(backend_.getData());
  }

  template <typename T>
  T& getUncertainty() {
    return *static_cast<T*>(backend_.getUncertainty());
  }

  template <typename T>
  const T& getUncertainty() const {
    return *static_cast<const T*>(backend_.getUncertainty());
  }

  size_t getSize() const { return backend_.getSize(); }

  // Metadata operations
  Observation& setMetadata(const std::string& key, const std::string& value) {
    backend_.setMetadata(key, value);
    return *this;
  }

  std::string getMetadata(const std::string& key) const {
    return backend_.getMetadata(key);
  }

  bool hasMetadata(const std::string& key) const {
    return backend_.hasMetadata(key);
  }

  // Observation information
  const std::vector<std::string>& getVariableNames() const {
    return backend_.getVariableNames();
  }

  bool hasVariable(const std::string& name) const {
    return backend_.hasVariable(name);
  }

  const std::vector<size_t>& getDimensions() const {
    return backend_.getDimensions();
  }

  // Spatiotemporal metadata
  Observation& setLocations(const std::vector<std::vector<double>>& locations) {
    backend_.setLocations(locations);
    return *this;
  }

  Observation& setTimes(const std::vector<double>& timestamps) {
    backend_.setTimes(timestamps);
    return *this;
  }

  const std::vector<std::vector<double>>& getLocations() const {
    return backend_.getLocations();
  }

  const std::vector<double>& getTimes() const { return backend_.getTimes(); }

  // Arithmetic operations
  Observation& add(const Observation& other) {
    backend_.add(other.backend_);
    return *this;
  }

  Observation& subtract(const Observation& other) {
    backend_.subtract(other.backend_);
    return *this;
  }

  Observation& multiply(double scalar) {
    backend_.multiply(scalar);
    return *this;
  }

  Observation operator+(const Observation& other) const {
    Observation result(*this);
    result.add(other);
    return result;
  }

  Observation operator-(const Observation& other) const {
    Observation result(*this);
    result.subtract(other);
    return result;
  }

  Observation operator*(double scalar) const {
    Observation result(*this);
    result.multiply(scalar);
    return result;
  }

  // Quality control
  Observation& setQualityFlags(const std::vector<int>& flags) {
    backend_.setQualityFlags(flags);
    return *this;
  }

  const std::vector<int>& getQualityFlags() const {
    return backend_.getQualityFlags();
  }

  Observation& setConfidenceValues(const std::vector<double>& values) {
    backend_.setConfidenceValues(values);
    return *this;
  }

  const std::vector<double>& getConfidenceValues() const {
    return backend_.getConfidenceValues();
  }
};

// Free function for scalar multiplication
template <typename Backend>
Observation<Backend> operator*(double scalar, const Observation<Backend>& obs) {
  return obs * scalar;
}

}  // namespace metada::framework