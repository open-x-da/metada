/**
 * @file AIPredictorImpl.hpp
 * @brief Implementation of the AI prediction capability for models
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains the implementation of the IAIPredictor interface
 * that delegates to a backend model implementation. It provides AI prediction
 * functionality for machine learning models like GraphCast and Pangu-Weather.
 */

#pragma once

#include <memory>
#include <vector>

#include "../../interfaces/IAIPredictor.hpp"

namespace metada::framework {

/**
 * @brief Implementation of the AI prediction capability
 *
 * @details
 * This class implements the IAIPredictor interface by delegating
 * to the backend model implementation. It provides AI prediction
 * functionality for machine learning models like GraphCast and Pangu-Weather.
 *
 * @tparam Backend The backend type implementing the model
 */
template <typename Backend>
class AIPredictorImpl : public IAIPredictor {
 private:
  Backend& backend_;  ///< Reference to the backend model

 public:
  /**
   * @brief Construct a new AIPredictorImpl
   *
   * @param backend Reference to the backend model
   */
  explicit AIPredictorImpl(Backend& backend) : backend_(backend) {}

  /**
   * @brief Load pre-trained weights from a file
   *
   * @param weightsPath Path to the pre-trained weights file
   */
  void loadWeights(const std::string& weightsPath) override {
    backend_.loadWeights(weightsPath);
  }

  /**
   * @brief Set the computational device for the model
   *
   * @param device The device to use (e.g., "cpu", "cuda:0", "tpu:0")
   */
  void setDevice(const std::string& device) override {
    backend_.setDevice(device);
  }

  /**
   * @brief Get the current computational device
   *
   * @return The current device (e.g., "cpu", "cuda:0", "tpu:0")
   */
  std::string getDevice() const override { return backend_.getDevice(); }

  /**
   * @brief Make a single-leadtime prediction from an initial state
   *
   * @param initialState The initial state
   * @param prediction The predicted state (output parameter)
   * @param leadTime The lead time for the prediction in seconds
   */
  void predict(const IState& initialState, IState& prediction,
               double leadTime) override {
    backend_.predict(initialState, prediction, leadTime);
  }

  /**
   * @brief Make multi-leadtime predictions from an initial state
   *
   * @param initialState The initial state
   * @param predictions Vector of pointers to predicted states (output
   * parameter)
   * @param leadTimes Vector of lead times for the predictions in seconds
   */
  void predictMultiple(const IState& initialState,
                       std::vector<std::unique_ptr<IState>>& predictions,
                       const std::vector<double>& leadTimes) override {
    backend_.predictMultiple(initialState, predictions, leadTimes);
  }

  /**
   * @brief Check if the model supports a specific lead time
   *
   * @param leadTime The lead time to check in seconds
   * @return true if the model supports the lead time, false otherwise
   */
  bool supportsLeadTime(double leadTime) const override {
    return backend_.supportsLeadTime(leadTime);
  }

  /**
   * @brief Get the supported lead times for the model
   *
   * @return Vector of supported lead times in seconds
   */
  std::vector<double> getSupportedLeadTimes() const override {
    return backend_.getSupportedLeadTimes();
  }

  /**
   * @brief Get the model architecture name
   *
   * @return The name of the model architecture
   */
  std::string getArchitecture() const override {
    return backend_.getArchitecture();
  }

  /**
   * @brief Get the version of the pre-trained weights
   *
   * @return The version of the pre-trained weights
   */
  std::string getWeightsVersion() const override {
    return backend_.getWeightsVersion();
  }

  /**
   * @brief Check if the model is deterministic
   *
   * @return true if the model is deterministic, false if it's stochastic
   */
  bool isDeterministic() const override { return backend_.isDeterministic(); }

  /**
   * @brief Check if the model has been properly loaded and is ready for
   * predictions
   *
   * @return true if the model is ready, false otherwise
   */
  bool isReady() const override { return backend_.isReady(); }
};

}  // namespace metada::framework