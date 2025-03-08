/**
 * @file IAIPredictor.hpp
 * @brief Interface defining the capability for AI-based prediction models
 * @ingroup repr
 *
 * @details
 * This header provides the capability interface for AI-based prediction models,
 * such as GraphCast, Pangu-Weather, FourCastNet, etc. Models that implement
 * this interface can make predictions using pre-trained neural networks or
 * other machine learning techniques.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

namespace metada::framework {

// Forward declarations
class IState;

/**
 * @brief Capability interface for AI-based prediction models
 *
 * @details
 * This interface defines the capabilities required for AI-based prediction
 * models. Models that implement this capability can make predictions using
 * pre-trained neural networks or other machine learning techniques.
 *
 * This capability is typically used by:
 * - GraphCast, Pangu-Weather, FourCastNet, etc.
 * - Any model using machine learning for prediction
 * - Hybrid models that incorporate AI components
 *
 * Implementations should ensure:
 * - Proper loading and management of pre-trained weights
 * - Efficient inference on appropriate hardware (CPU, GPU, TPU)
 * - Correct handling of input/output formats
 */
class IAIPredictor {
 public:
  /**
   * @brief Virtual destructor
   */
  virtual ~IAIPredictor() = default;

  /**
   * @brief Load pre-trained weights from a file
   *
   * @param weightsPath Path to the pre-trained weights file
   * @throws std::runtime_error if loading fails
   */
  virtual void loadWeights(const std::string& weightsPath) = 0;

  /**
   * @brief Set the computational device for the model
   *
   * @param device The device to use (e.g., "cpu", "cuda:0", "tpu:0")
   * @throws std::runtime_error if the device is not available
   */
  virtual void setDevice(const std::string& device) = 0;

  /**
   * @brief Get the current computational device
   *
   * @return The current device (e.g., "cpu", "cuda:0", "tpu:0")
   */
  virtual std::string getDevice() const = 0;

  /**
   * @brief Make a single-leadtime prediction from an initial state
   *
   * @param initialState The initial state
   * @param prediction The predicted state (output parameter)
   * @param leadTime The lead time for the prediction in seconds
   * @throws std::runtime_error if prediction fails
   */
  virtual void predict(const IState& initialState, IState& prediction,
                       double leadTime) = 0;

  /**
   * @brief Make multi-leadtime predictions from an initial state
   *
   * @param initialState The initial state
   * @param predictions Vector of pointers to predicted states (output
   * parameter)
   * @param leadTimes Vector of lead times for the predictions in seconds
   * @throws std::runtime_error if prediction fails
   */
  virtual void predictMultiple(
      const IState& initialState,
      std::vector<std::unique_ptr<IState>>& predictions,
      const std::vector<double>& leadTimes) = 0;

  /**
   * @brief Check if the model supports a specific lead time
   *
   * @param leadTime The lead time to check in seconds
   * @return true if the model supports the lead time, false otherwise
   */
  virtual bool supportsLeadTime(double leadTime) const = 0;

  /**
   * @brief Get the supported lead times for the model
   *
   * @return Vector of supported lead times in seconds
   */
  virtual std::vector<double> getSupportedLeadTimes() const = 0;

  /**
   * @brief Get the model architecture name
   *
   * @return The name of the model architecture (e.g., "GraphCast",
   * "Pangu-Weather")
   */
  virtual std::string getArchitecture() const = 0;

  /**
   * @brief Get the version of the pre-trained weights
   *
   * @return The version of the pre-trained weights
   */
  virtual std::string getWeightsVersion() const = 0;

  /**
   * @brief Check if the model is deterministic
   *
   * @return true if the model is deterministic, false if it's stochastic
   */
  virtual bool isDeterministic() const = 0;

  /**
   * @brief Check if the model has been properly loaded and is ready for
   * predictions
   *
   * @return true if the model is ready, false otherwise
   */
  virtual bool isReady() const = 0;
};

}  // namespace metada::framework