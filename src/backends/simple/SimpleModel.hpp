/**
 * @file SimpleModel.hpp
 * @brief Simple identity model implementation
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides a simple identity model implementation for testing and
 * development purposes. The model simply copies the input state to the output
 * state without any transformation, making it useful for testing data
 * assimilation algorithms.
 */

#pragma once

#include <stdexcept>
#include <string>

#include "SimpleState.hpp"

namespace metada::backends::simple {

/**
 * @brief Simple identity model backend implementation
 *
 * @details
 * This class implements a basic identity model that:
 * - Copies input state to output state without modification
 * - Supports initialization and configuration
 * - Provides proper resource management
 * - Satisfies the ModelBackendImpl concept requirements
 *
 * The identity model is useful for:
 * - Testing data assimilation algorithms
 * - Development and debugging
 * - Cases where no model evolution is needed
 */
class SimpleModel {
 public:
  // Delete default constructor
  SimpleModel() = delete;

  // Delete copy constructor and assignment
  SimpleModel(const SimpleModel&) = delete;
  SimpleModel& operator=(const SimpleModel&) = delete;

  /**
   * @brief Constructor that initializes from configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  explicit SimpleModel(const ConfigBackend& config) : initialized_(false) {
    initialize(config);
  }

  /**
   * @brief Move constructor
   * @param other Model to move from
   */
  SimpleModel(SimpleModel&& other) noexcept : initialized_(other.initialized_) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   * @param other Model to move from
   * @return Reference to this model
   */
  SimpleModel& operator=(SimpleModel&& other) noexcept {
    if (this != &other) {
      initialized_ = other.initialized_;
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Destructor
   */
  ~SimpleModel() = default;

  /**
   * @brief Initialize the model with configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  void initialize([[maybe_unused]] const ConfigBackend& config) {
    if (initialized_) {
      return;  // Already initialized
    }

    try {
      // Read model parameters from config if needed
      // For identity model, we don't need any special parameters
      initialized_ = true;
    } catch (const std::exception& e) {
      throw std::runtime_error("SimpleModel initialization failed: " +
                               std::string(e.what()));
    }
  }

  /**
   * @brief Check if the model is initialized
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Reset the model to initial state
   */
  void reset() {
    if (!initialized_) {
      throw std::runtime_error("SimpleModel not initialized");
    }
    // For identity model, reset is a no-op
  }

  /**
   * @brief Finalize the model, releasing resources
   */
  void finalize() {
    if (!initialized_) {
      return;  // No need to finalize if not initialized
    }
    initialized_ = false;
  }

  /**
   * @brief Run the model (identity operation)
   *
   * @details
   * The identity model simply copies the input state to the output state.
   * This is useful for testing data assimilation algorithms where no
   * model evolution is needed.
   *
   * @param initialState Input state
   * @param finalState Output state (will be identical to input)
   */
  void run(const SimpleState& initialState, SimpleState& finalState) {
    if (!initialized_) {
      throw std::runtime_error("SimpleModel not initialized");
    }

    try {
      // Identity operation: copy input to output
      // Get data from input state
      const auto input_data =
          static_cast<const double*>(initialState.getData());
      size_t state_size = initialState.size();

      // Copy data to output state
      auto output_data = static_cast<double*>(finalState.getData());
      for (size_t i = 0; i < state_size; ++i) {
        output_data[i] = input_data[i];
      }
    } catch (const std::exception& e) {
      throw std::runtime_error("SimpleModel run failed: " +
                               std::string(e.what()));
    }
  }

 private:
  bool initialized_;  ///< Initialization status
};

}  // namespace metada::backends::simple