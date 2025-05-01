/**
 * @file Model.hpp
 * @brief Adapter implementing the IModel interface for various backends
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains a template implementation of the IModel interface
 * that works with various backend model implementations. It serves as an
 * adapter between the core interfaces and concrete model implementations.
 */

#pragma once

#include <stdexcept>
#include <string>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "ModelConcepts.hpp"
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
 * @brief Forward declaration of State class
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class State;

/**
 * @brief Main model class template providing a generic interface to model
 * implementations
 *
 * @details
 * This class adapts the backend model implementation to the IModel interface.
 * It provides capability-based query methods for optional capabilities using
 * a component-based design that leverages C++20 features.
 *
 * The backend tag must satisfy the ModelBackendType concept, which ensures
 * it provides valid backend implementation types through BackendTraits.
 *
 * @tparam BackendTag The tag type that defines the backend through
 * BackendTraits
 */
template <typename BackendTag>
  requires ModelBackendType<BackendTag>
class Model : private NonCopyable {
 public:
  using ModelBackend = typename traits::BackendTraits<BackendTag>::ModelBackend;

  /**
   * @brief Default constructor is deleted
   *
   * Model objects must be constructed with a Config object.
   */
  Model() = delete;

  /**
   * @brief Destructor
   */
  ~Model() = default;  // Default destructor

  /**
   * @brief Construct a new Model with the given configuration
   *
   * @param config Configuration object with model parameters
   *
   * @note The model must be initialized with initialize() before use, and when
   * no longer needed, resources should be released with finalize().
   */
  explicit Model(const Config<BackendTag>& config)
      : backend_(config.backend()) {}

  /**
   * @brief Move constructor
   *
   * @param other The other Model object to move from
   */
  Model(Model&& other) noexcept : backend_(std::move(other.backend_)) {}

  /**
   * @brief Move assignment operator
   *
   * @param other The other Model object to move from
   * @return Reference to this Model object
   */
  Model& operator=(Model&& other) noexcept {
    if (this != &other) {
      backend_ = std::move(other.backend_);
    }
    return *this;
  }

  /**
   * @brief Initialize the model with the given configuration
   *
   * This method must be called before using the model.
   *
   * @param config Configuration object with model parameters
   * @throws std::runtime_error if initialization fails
   */
  void initialize(const Config<BackendTag>& config) {
    if (backend_.isInitialized()) {
      return;  // Already initialized
    }

    try {
      backend_.initialize(config.backend());
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model initialization failed: ") +
                               e.what());
    }
  }

  /**
   * @brief Reset the model to its initial state
   *
   * This method resets the model state without releasing resources.
   *
   * @throws std::runtime_error if reset fails
   */
  void reset() {
    if (!backend_.isInitialized()) {
      throw std::runtime_error("Model not initialized");
    }

    try {
      backend_.reset();
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model reset failed: ") + e.what());
    }
  }

  /**
   * @brief Finalize the model, releasing all resources
   *
   * @throws std::runtime_error if finalization fails
   */
  void finalize() {
    if (!backend_.isInitialized()) {
      return;  // No need to finalize if not initialized
    }

    try {
      backend_.finalize();
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model finalization failed: ") +
                               e.what());
    }
  }

  /**
   * @brief Check if the model is initialized
   *
   * @return true if the model is initialized, false otherwise
   */
  bool isInitialized() const { return backend_.isInitialized(); }

  /**
   * @brief Type-safe run method that works directly with State templates
   *
   * @tparam StateType The type of state used by the model
   * @param initialState The initial state of the model
   * @param finalState The final state after model run (output parameter)
   * @throws std::runtime_error if the model run fails
   */
  void run(const State<BackendTag>& initialState,
           State<BackendTag>& finalState) {
    if (!backend_.isInitialized()) {
      throw std::runtime_error("Model not initialized");
    }

    try {
      backend_.run(initialState.backend(), finalState.backend());
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model run failed: ") + e.what());
    }
  }

  /**
   * @brief Get the backend instance
   *
   * This method provides direct access to the backend model implementation.
   * Use with caution as it bypasses the adapter's type safety and error
   * handling.
   *
   * @return Reference to the backend instance
   */
  ModelBackend& backend() { return backend_; }

  /**
   * @brief Get the backend instance (const version)
   *
   * @return Const reference to the backend instance
   */
  const ModelBackend& backend() const { return backend_; }

 private:
  ModelBackend backend_;
};

}  // namespace metada::framework