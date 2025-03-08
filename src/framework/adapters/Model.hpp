/**
 * @file Model.hpp
 * @brief Template class providing a generic interface to model implementations
 * @ingroup repr
 * @author Metada Framework Team
 *
 * @details
 * This header provides a template class that wraps model backend
 * implementations and provides a unified interface for model operations. The
 * Model class delegates operations to the backend while providing type
 * safety and a consistent API.
 *
 * The Model class template is designed to:
 * - Provide a generic interface to different model backend implementations
 * - Support initialization from configuration objects
 * - Enable core model operations (run, etc)
 * - Support capability interfaces for specialized functionality (time stepping,
 * AI prediction, etc)
 * - Implement proper exception handling and error reporting
 *
 * @see IModel
 * @see ITimeStepper
 * @see IAIPredictor
 * @see IBatchProcessor
 * @see Config
 */

#pragma once

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "../common/utils/NonCopyable.hpp"
#include "../interfaces/IAIPredictor.hpp"
#include "../interfaces/IBatchProcessor.hpp"
#include "../interfaces/IHardwareAccelerator.hpp"
#include "../interfaces/IModel.hpp"
#include "../interfaces/ITimeStepper.hpp"
#include "model/AIPredictorImpl.hpp"
#include "model/BatchProcessorImpl.hpp"
#include "model/HardwareAcceleratorImpl.hpp"
#include "model/TimeStepperImpl.hpp"

namespace metada::framework {

// Forward declarations
template <typename T>
class Config;

template <typename T>
class State;

/**
 * @brief Main model class template providing a generic interface to model
 * implementations
 *
 * @details
 * This class template wraps a model backend implementation and provides a
 * type-safe interface for all model operations. It delegates operations to
 * the backend through composition, following the adapter pattern.
 *
 * Key features:
 * - Type safety through templates
 * - Composition-based design (not inheritance)
 * - Capability interface implementations as needed
 * - Proper error handling and diagnostics
 * - Consistent interface across different backends
 *
 * The backend must implement the necessary operations to provide the core
 * functionality required by the model.
 *
 * @tparam Backend The backend type implementing the model
 * @tparam ConfigType The type of configuration object used
 */
template <typename Backend,
          typename ConfigType = Config<typename Backend::config_type>>
class Model : private NonCopyable {
 private:
  Backend backend_;          ///< Instance of the model backend
  bool initialized_{false};  ///< Initialization state
  ConfigType configCopy_;    ///< Copy of the configuration

  // Capability interface implementations
  std::optional<TimeStepperImpl<Backend>> timeStepper_;
  std::optional<AIPredictorImpl<Backend>> aiPredictor_;
  std::optional<BatchProcessorImpl<Backend>> batchProcessor_;
  std::optional<HardwareAcceleratorImpl<Backend>> hardwareAccelerator_;

 public:
  /**
   * @brief Default constructor is deleted
   *
   * Model objects must be constructed with a Config object.
   */
  Model() = delete;

  /**
   * @brief Construct a new Model with the given configuration
   *
   * @param config Configuration object with model parameters
   */
  explicit Model(const ConfigType& config)
      : backend_(config), configCopy_(config) {
    initializeCapabilities();
  }

  /**
   * @brief Initialize the model with the given configuration
   *
   * @param config Configuration object with model parameters
   * @throws std::runtime_error if initialization fails
   */
  void initialize(const IConfig& config) {
    try {
      backend_.initialize(config);
      initialized_ = true;
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model initialization failed: ") +
                               e.what());
    }
  }

  /**
   * @brief Reset the model to its initial state
   */
  void reset() {
    try {
      backend_.reset();
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model reset failed: ") + e.what());
    }
  }

  /**
   * @brief Check if the model is initialized
   *
   * @return true if the model is initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Get the configuration object used to create this model
   *
   * @return Const reference to the configuration object
   */
  const ConfigType& getConfig() const { return configCopy_; }

  /**
   * @brief Get a parameter from the model
   *
   * @param name The name of the parameter
   * @return The value of the parameter as a string
   * @throws std::out_of_range if the parameter doesn't exist
   */
  std::string getParameter(const std::string& name) const {
    try {
      return backend_.getParameter(name);
    } catch (const std::exception&) {
      throw std::out_of_range("Parameter not found: " + name);
    }
  }

  /**
   * @brief Set a parameter in the model
   *
   * @param name The name of the parameter
   * @param value The value of the parameter
   * @throws std::invalid_argument if the parameter is invalid
   */
  void setParameter(const std::string& name, const std::string& value) {
    try {
      backend_.setParameter(name, value);
      // We don't need to update a local parameters map anymore since we have
      // the full config
    } catch (const std::exception& e) {
      throw std::invalid_argument(std::string("Invalid parameter: ") +
                                  e.what());
    }
  }

  /**
   * @brief Run the model from initial state to final state
   *
   * This is the primary execution method.
   *
   * @param initialState The initial state of the model
   * @param finalState The final state after model execution (output parameter)
   * @param startTime The start time for the model run
   * @param endTime The end time for the model run
   * @throws std::runtime_error if the model run fails
   */
  void run(const IState& initialState, IState& finalState, double startTime,
           double endTime) {
    if (!initialized_) {
      throw std::runtime_error("Model not initialized");
    }

    try {
      backend_.run(initialState, finalState, startTime, endTime);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model run failed: ") + e.what());
    }
  }

  /**
   * @brief Get the model's time stepping capability
   *
   * @return Pointer to the model's time stepping capability, or nullptr if not
   * supported
   */
  ITimeStepper* getTimeStepper() {
    return timeStepper_ ? &*timeStepper_ : nullptr;
  }

  /**
   * @brief Get the model's AI prediction capability
   *
   * @return Pointer to the model's AI prediction capability, or nullptr if not
   * supported
   */
  IAIPredictor* getAIPredictor() {
    return aiPredictor_ ? &*aiPredictor_ : nullptr;
  }

  /**
   * @brief Get the model's batch processing capability
   *
   * @return Pointer to the model's batch processing capability, or nullptr if
   * not supported
   */
  IBatchProcessor* getBatchProcessor() {
    return batchProcessor_ ? &*batchProcessor_ : nullptr;
  }

  /**
   * @brief Get the model's hardware acceleration capability
   *
   * @return Pointer to the model's hardware acceleration capability, or nullptr
   * if not supported
   */
  IHardwareAccelerator* getHardwareAccelerator() {
    return hardwareAccelerator_ ? &*hardwareAccelerator_ : nullptr;
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
  Backend& getBackend() { return backend_; }

  /**
   * @brief Get the backend instance (const version)
   *
   * @return Const reference to the backend instance
   */
  const Backend& getBackend() const { return backend_; }

  /**
   * @brief Type-safe run method that works directly with State templates
   *
   * @tparam StateType The type of state used by the model
   * @param initialState The initial state of the model
   * @param finalState The final state after model run (output parameter)
   * @param startTime The start time for the model run
   * @param endTime The end time for the model run
   * @throws std::runtime_error if the model run fails
   */
  template <typename StateType>
  void run(const State<StateType>& initialState, State<StateType>& finalState,
           double startTime, double endTime) {
    if (!initialized_) {
      throw std::runtime_error("Model not initialized");
    }

    try {
      backend_.run(initialState.getBackend(), finalState.getBackend(),
                   startTime, endTime);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model run failed: ") + e.what());
    }
  }

 private:
  /**
   * @brief Initialize the capabilities based on what the backend supports
   */
  void initializeCapabilities() {
    // Check if backend supports time stepping
    if constexpr (requires(Backend& b) {
                    b.step(std::declval<IState&>(), std::declval<IState&>());
                  }) {
      timeStepper_.emplace(backend_);
    }

    // Check if backend supports AI prediction
    if constexpr (requires(Backend& b) {
                    b.predict(std::declval<IState&>(), std::declval<IState&>(),
                              double());
                  }) {
      aiPredictor_.emplace(backend_);
    }

    // Check if backend supports batch processing
    if constexpr (requires(Backend& b) {
                    b.processBatch(
                        std::declval<std::vector<std::shared_ptr<IState>>&>(),
                        std::declval<std::vector<std::shared_ptr<IState>>&>(),
                        double(), double());
                  }) {
      batchProcessor_.emplace(backend_);
    }

    // Check if backend supports hardware acceleration
    if constexpr (requires(Backend& b) {
                    b.setDevice(std::declval<std::string>());
                  }) {
      hardwareAccelerator_.emplace(backend_);
    }
  }
};

}  // namespace metada::framework