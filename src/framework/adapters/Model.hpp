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

#include <concepts>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "../common/utils/NonCopyable.hpp"
#include "../interfaces/IAIPredictor.hpp"
#include "../interfaces/IBatchProcessor.hpp"
#include "../interfaces/IHardwareAccelerator.hpp"
#include "../interfaces/IModel.hpp"
#include "../interfaces/ITimeStepper.hpp"
// #include "model/AIPredictorImpl.hpp"
// #include "model/BatchProcessorImpl.hpp"
//  #include "model/HardwareAcceleratorImpl.hpp"
// #include "model/TimeStepperImpl.hpp"

namespace metada::framework {

// Forward declarations
template <typename T>
class Config;

template <typename T>
class State;

// C++20 concepts for capability detection
// Time stepping capability - check for all required methods
template <typename T>
concept HasTimeSteppingCapability =
    requires(T& backend, IState& state1, IState& state2, double dt, bool enable,
             std::string scheme) {
      backend.step(state1, state2);
      { backend.getTimeStep() } -> std::convertible_to<double>;
      backend.setTimeStep(dt);
      { backend.getMaxStableTimeStep(state1) } -> std::convertible_to<double>;
      { backend.usesAdaptiveTimeStep() } -> std::convertible_to<bool>;
      backend.setAdaptiveTimeStep(enable);
      { backend.getCurrentTime() } -> std::convertible_to<double>;
      backend.setCurrentTime(dt);
      { backend.getTimeSteppingScheme() } -> std::convertible_to<std::string>;
    };

// AI prediction capability - check for all required methods
template <typename T>
concept HasAIPredictionCapability =
    requires(T& backend, IState& state1, IState& state2, double leadTime,
             std::vector<std::unique_ptr<IState>>& predictions,
             std::vector<double>& leadTimes, std::string weightsPath,
             std::string device) {
      backend.predict(state1, state2, leadTime);
      backend.predictMultiple(state1, predictions, leadTimes);
      { backend.supportsLeadTime(leadTime) } -> std::convertible_to<bool>;
      {
        backend.getSupportedLeadTimes()
      } -> std::convertible_to<std::vector<double>>;
      { backend.getArchitecture() } -> std::convertible_to<std::string>;
      { backend.getWeightsVersion() } -> std::convertible_to<std::string>;
      { backend.isDeterministic() } -> std::convertible_to<bool>;
      { backend.isReady() } -> std::convertible_to<bool>;
      backend.loadWeights(weightsPath);
      backend.setDevice(device);
      { backend.getDevice() } -> std::convertible_to<std::string>;
    };

// Batch processing capability - check for all required methods
template <typename T>
concept HasBatchProcessingCapability =
    requires(T& backend, std::vector<std::shared_ptr<IState>>& input,
             std::vector<std::shared_ptr<IState>>& output, double startTime,
             double endTime, size_t batchSize, bool enable) {
      backend.processBatch(input, output, startTime, endTime);
      { backend.getMaxBatchSize() } -> std::convertible_to<size_t>;
      backend.setBatchSize(batchSize);
      { backend.getBatchSize() } -> std::convertible_to<size_t>;
      { backend.supportsVariableBatchSize() } -> std::convertible_to<bool>;
      { backend.isBatchProcessingSupported() } -> std::convertible_to<bool>;
      { backend.getOptimalBatchSize() } -> std::convertible_to<size_t>;
      { backend.supportsAutomaticBatching() } -> std::convertible_to<bool>;
      backend.setAutomaticBatching(enable);
    };

// Hardware acceleration capability - check for all required methods
template <typename T>
concept HasHardwareAccelerationCapability =
    requires(T& backend, std::string device, size_t memSize, bool enable) {
      backend.setDevice(device);
      { backend.getDevice() } -> std::convertible_to<std::string>;
      {
        backend.getAvailableDevices()
      } -> std::convertible_to<std::vector<std::string>>;
      { backend.supportsDeviceType(device) } -> std::convertible_to<bool>;
      { backend.getAvailableMemory() } -> std::convertible_to<size_t>;
      backend.setMemoryLimit(memSize);
      { backend.getMemoryLimit() } -> std::convertible_to<size_t>;
      backend.setAsynchronousExecution(enable);
      { backend.isAsynchronousExecutionEnabled() } -> std::convertible_to<bool>;
      backend.synchronize();
    };

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
template <typename Backend>
class Model : private NonCopyable {
 private:
  Backend backend_;          ///< Instance of the model backend
  bool initialized_{false};  ///< Initialization state

  // Capability interface implementations
  // std::optional<TimeStepperImpl<Backend>> timeStepper_;
  // std::optional<AIPredictorImpl<Backend>> aiPredictor_;
  // std::optional<BatchProcessorImpl<Backend>> batchProcessor_;
  // std::optional<HardwareAcceleratorImpl<Backend>> hardwareAccelerator_;

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
  template <typename T>
  explicit Model(const Config<T>& config) : backend_(config) {
    initializeCapabilities();
  }

  /**
   * @brief Initialize the model with the given configuration
   *
   * @param config Configuration object with model parameters
   * @throws std::runtime_error if initialization fails
   */
  template <typename T>
  void initialize(const Config<T>& config) {
    try {
      backend_.initialize(config.backend());
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
  // ITimeStepper* getTimeStepper() {
  //   return timeStepper_ ? &*timeStepper_ : nullptr;
  // }

  /**
   * @brief Get the model's AI prediction capability
   *
   * @return Pointer to the model's AI prediction capability, or nullptr if not
   * supported
   */
  // IAIPredictor* getAIPredictor() {
  //   return aiPredictor_ ? &*aiPredictor_ : nullptr;
  // }

  /**
   * @brief Get the model's batch processing capability
   *
   * @return Pointer to the model's batch processing capability, or nullptr if
   * not supported
   */
  // IBatchProcessor* getBatchProcessor() {
  //   return batchProcessor_ ? &*batchProcessor_ : nullptr;
  // }

  /**
   * @brief Get the model's hardware acceleration capability
   *
   * @return Pointer to the model's hardware acceleration capability, or nullptr
   * if not supported
   */
  // IHardwareAccelerator* getHardwareAccelerator() {
  //   return hardwareAccelerator_ ? &*hardwareAccelerator_ : nullptr;
  // }

  /**
   * @brief Get the backend instance
   *
   * This method provides direct access to the backend model implementation.
   * Use with caution as it bypasses the adapter's type safety and error
   * handling.
   *
   * @return Reference to the backend instance
   */
  Backend& backend() { return backend_; }

  /**
   * @brief Get the backend instance (const version)
   *
   * @return Const reference to the backend instance
   */
  const Backend& backend() const { return backend_; }

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
      backend_.run(initialState.backend(), finalState.backend(), startTime,
                   endTime);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Model run failed: ") + e.what());
    }
  }

 private:
  /**
   * @brief Initialize the capabilities based on what the backend supports
   *
   * Uses C++20 concepts to check if the backend supports various capabilities.
   */
  void initializeCapabilities() {
    // Check if backend supports time stepping
    // if constexpr (HasTimeSteppingCapability<Backend>) {
    //   timeStepper_.emplace(backend_);
    // }

    // Check if backend supports AI prediction
    // if constexpr (HasAIPredictionCapability<Backend>) {
    //   aiPredictor_.emplace(backend_);
    // }

    // Check if backend supports batch processing
    // if constexpr (HasBatchProcessingCapability<Backend>) {
    //   batchProcessor_.emplace(backend_);
    // }

    // Check if backend supports hardware acceleration
    // The backend must provide ALL of the required methods for hardware
    // acceleration
    // if constexpr (HasHardwareAccelerationCapability<Backend>) {
    //   // Make sure the backend has all the methods required by
    //   // HardwareAcceleratorImpl
    //   if constexpr (requires(Backend& b, std::string device, size_t memSize,
    //                             bool enable) {
    //                       b.setDevice(device);
    //                       { b.getDevice() } ->
    //                       std::convertible_to<std::string>;
    //                       {
    //                         b.supportsDeviceType(device)
    //                       } -> std::convertible_to<bool>;
    //                       {
    //                         b.getAvailableDevices()
    //                       } -> std::convertible_to<std::vector<std::string>>;
    //                       { b.getAvailableMemory() } ->
    //                       std::convertible_to<size_t>;
    //                       b.setMemoryLimit(memSize);
    //                       { b.getMemoryLimit() } ->
    //                       std::convertible_to<size_t>;
    //                       b.setAsynchronousExecution(enable);
    //                       {
    //                         b.isAsynchronousExecutionEnabled()
    //                       } -> std::convertible_to<bool>;
    //                       b.synchronize();
    //                     }) {
    //     hardwareAccelerator_.emplace(backend_);
    //   }
    // }
  }
};

}  // namespace metada::framework