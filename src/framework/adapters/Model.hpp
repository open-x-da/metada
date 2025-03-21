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

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "../common/utils/NonCopyable.hpp"
#include "../interfaces/IAIPredictor.hpp"
#include "../interfaces/IBatchProcessor.hpp"
#include "../interfaces/IHardwareAccelerator.hpp"
#include "../interfaces/ITimeStepper.hpp"
#include "State.hpp"
#include "model/AIPredictorImpl.hpp"
#include "model/BatchProcessorImpl.hpp"
#include "model/HardwareAcceleratorImpl.hpp"
#include "model/TimeStepperImpl.hpp"

namespace metada::framework {

// Forward declarations
template <typename T>
class Config;

// C++20 concepts for capability detection
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

template <typename T>
concept HasHardwareAccelerationCapability =
    requires(T& backend, std::string device, size_t memSize, bool enable) {
      backend.setDevice(device);
      { backend.getDevice() } -> std::convertible_to<std::string>;
      { backend.supportsDeviceType(device) } -> std::convertible_to<bool>;
      {
        backend.getAvailableDevices()
      } -> std::convertible_to<std::vector<std::string>>;
      { backend.getAvailableMemory() } -> std::convertible_to<size_t>;
      backend.setMemoryLimit(memSize);
      { backend.getMemoryLimit() } -> std::convertible_to<size_t>;
      backend.setAsynchronousExecution(enable);
      { backend.isAsynchronousExecutionEnabled() } -> std::convertible_to<bool>;
      backend.synchronize();
    };

// Component implementation for TimeSteppingCapability
template <typename Backend, bool HasCapability = false>
struct TimeStepperComponent {
  void initialize(Backend&) {}
  void reset() {}
  ITimeStepper* get() const { return nullptr; }
};

template <typename Backend>
struct TimeStepperComponent<Backend, true> {
  std::optional<TimeStepperImpl<Backend>> impl_;

  void initialize(Backend& backend) { impl_.emplace(backend); }
  void reset() { impl_.reset(); }
  ITimeStepper* get() const { return impl_ ? &*impl_ : nullptr; }
};

// Component implementation for AIPredictionCapability
template <typename Backend, bool HasCapability = false>
struct AIPredictorComponent {
  void initialize(Backend&) {}
  void reset() {}
  IAIPredictor* get() const { return nullptr; }
};

template <typename Backend>
struct AIPredictorComponent<Backend, true> {
  std::optional<AIPredictorImpl<Backend>> impl_;

  void initialize(Backend& backend) { impl_.emplace(backend); }
  void reset() { impl_.reset(); }
  IAIPredictor* get() const { return impl_ ? &*impl_ : nullptr; }
};

// Component implementation for BatchProcessingCapability
template <typename Backend, bool HasCapability = false>
struct BatchProcessorComponent {
  void initialize(Backend&) {}
  void reset() {}
  IBatchProcessor* get() const { return nullptr; }
};

template <typename Backend>
struct BatchProcessorComponent<Backend, true> {
  std::optional<BatchProcessorImpl<Backend>> impl_;

  void initialize(Backend& backend) { impl_.emplace(backend); }
  void reset() { impl_.reset(); }
  IBatchProcessor* get() const { return impl_ ? &*impl_ : nullptr; }
};

// Component implementation for HardwareAccelerationCapability
template <typename Backend, bool HasCapability = false>
struct HardwareAcceleratorComponent {
  void initialize(Backend&) {}
  void reset() {}
  IHardwareAccelerator* get() const { return nullptr; }
};

template <typename Backend>
struct HardwareAcceleratorComponent<Backend, true> {
  std::optional<HardwareAcceleratorImpl<Backend>> impl_;

  void initialize(Backend& backend) { impl_.emplace(backend); }
  void reset() { impl_.reset(); }
  IHardwareAccelerator* get() const { return impl_ ? &*impl_ : nullptr; }
};

/**
 * @brief Main model class template providing a generic interface to model
 * implementations
 *
 * @details
 * This class adapts the backend model implementation to the IModel interface.
 * It provides capability-based query methods for optional capabilities using
 * a component-based design that leverages C++20 features.
 *
 * @tparam Backend The backend model type
 */
template <typename Backend>
class Model : private NonCopyable {
 private:
  Backend backend_;          ///< Instance of the model backend
  bool initialized_{false};  ///< Initialization state

  // Capabilities implemented as components with zero overhead for unsupported
  // capabilities
  [[no_unique_address]] TimeStepperComponent<
      Backend, HasTimeSteppingCapability<Backend>> timeStepper_;
  [[no_unique_address]] AIPredictorComponent<
      Backend, HasAIPredictionCapability<Backend>> aiPredictor_;
  [[no_unique_address]] BatchProcessorComponent<
      Backend, HasBatchProcessingCapability<Backend>> batchProcessor_;
  [[no_unique_address]] HardwareAcceleratorComponent<
      Backend, HasHardwareAccelerationCapability<Backend>> hardwareAccelerator_;

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
   *
   * @note The model must be initialized with initialize() before use, and when
   * no longer needed, resources should be released with finalize().
   */
  template <typename T>
  explicit Model(const Config<T>& config) : backend_(config) {
    initializeCapabilities();
  }

  /**
   * @brief Initialize the model with the given configuration
   *
   * This method must be called before using the model.
   *
   * @param config Configuration object with model parameters
   * @throws std::runtime_error if initialization fails
   */
  template <typename T>
  void initialize(const Config<T>& config) {
    if (initialized_) {
      return;  // Already initialized
    }

    try {
      backend_.initialize(config.backend());
      initializeCapabilities();
      initialized_ = true;
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
    if (!initialized_) {
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
    if (!initialized_) {
      return;  // No need to finalize if not initialized
    }

    try {
      backend_.finalize();
      // Reset all capabilities
      timeStepper_.reset();
      aiPredictor_.reset();
      batchProcessor_.reset();
      hardwareAccelerator_.reset();
      initialized_ = false;
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
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Get a model parameter value
   *
   * @param name The parameter name
   * @return The parameter value as a string
   * @throws std::runtime_error if the parameter does not exist or cannot be
   * retrieved
   */
  std::string getParameter(const std::string& name) const {
    try {
      return backend_.getParameter(name);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Failed to get parameter '") + name +
                               "': " + e.what());
    }
  }

  /**
   * @brief Set a model parameter value
   *
   * @param name The parameter name
   * @param value The parameter value as a string
   * @throws std::runtime_error if the parameter does not exist or cannot be set
   */
  void setParameter(const std::string& name, const std::string& value) {
    try {
      backend_.setParameter(name, value);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Failed to set parameter '") + name +
                               "' to '" + value + "': " + e.what());
    }
  }

  /**
   * @brief Get the model's time stepping capability
   *
   * @return Pointer to the model's time stepping capability, or nullptr
   * if not supported
   */
  ITimeStepper* getTimeStepper() { return timeStepper_.get(); }

  /**
   * @brief Get the model's AI prediction capability
   *
   * @return Pointer to the model's AI prediction capability, or nullptr
   * if not supported
   */
  IAIPredictor* getAIPredictor() { return aiPredictor_.get(); }

  /**
   * @brief Get the model's batch processing capability
   *
   * @return Pointer to the model's batch processing capability, or nullptr
   * if not supported
   */
  IBatchProcessor* getBatchProcessor() { return batchProcessor_.get(); }

  /**
   * @brief Get the model's hardware acceleration capability
   *
   * @return Pointer to the model's hardware acceleration capability, or nullptr
   * if not supported
   */
  IHardwareAccelerator* getHardwareAccelerator() {
    return hardwareAccelerator_.get();
  }

  /**
   * @brief Run the model from initialState to finalState
   *
   * @param initialState The initial state of the model
   * @param finalState The final state after model run (output parameter)
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

 private:
  /**
   * @brief Initialize the capabilities based on what the backend supports
   *
   * Uses C++20 concepts to check if the backend supports various capabilities.
   */
  void initializeCapabilities() {
    // Initialize all component capabilities
    // The component will handle whether the initialization actually happens
    timeStepper_.initialize(backend_);
    aiPredictor_.initialize(backend_);
    batchProcessor_.initialize(backend_);
    hardwareAccelerator_.initialize(backend_);
  }
};

}  // namespace metada::framework