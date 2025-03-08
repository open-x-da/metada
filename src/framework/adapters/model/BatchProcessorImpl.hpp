/**
 * @file BatchProcessorImpl.hpp
 * @brief Implementation of the batch processing capability for models
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains the implementation of the IBatchProcessor interface
 * that delegates to a backend model implementation. It provides batch
 * processing functionality for both AI models and physical models that support
 * it.
 */

#pragma once

#include <memory>
#include <vector>

#include "../../interfaces/IBatchProcessor.hpp"

namespace metada::framework {

/**
 * @brief Implementation of the batch processing capability
 *
 * @details
 * This class implements the IBatchProcessor interface by delegating
 * to the backend model implementation. It provides batch processing
 * functionality for both AI models and physical models that support it.
 *
 * @tparam Backend The backend type implementing the model
 */
template <typename Backend>
class BatchProcessorImpl : public IBatchProcessor {
 private:
  Backend& backend_;  ///< Reference to the backend model

 public:
  /**
   * @brief Construct a new BatchProcessorImpl
   *
   * @param backend Reference to the backend model
   */
  explicit BatchProcessorImpl(Backend& backend) : backend_(backend) {}

  /**
   * @brief Get the maximum batch size supported by the model
   *
   * @return The maximum batch size
   */
  size_t getMaxBatchSize() const override { return backend_.getMaxBatchSize(); }

  /**
   * @brief Set the batch size for the model
   *
   * @param batchSize The batch size to use
   */
  void setBatchSize(size_t batchSize) override {
    backend_.setBatchSize(batchSize);
  }

  /**
   * @brief Get the current batch size
   *
   * @return The current batch size
   */
  size_t getBatchSize() const override { return backend_.getBatchSize(); }

  /**
   * @brief Process a batch of states at once
   *
   * @param inputBatch Vector of input states
   * @param outputBatch Vector of output states (output parameter)
   * @param startTime The start time for the model run
   * @param endTime The end time for the model run
   */
  void processBatch(const std::vector<std::shared_ptr<IState>>& inputBatch,
                    std::vector<std::shared_ptr<IState>>& outputBatch,
                    double startTime, double endTime) override {
    backend_.processBatch(inputBatch, outputBatch, startTime, endTime);
  }

  /**
   * @brief Check if the model supports variable batch sizes
   *
   * @return true if the model supports variable batch sizes, false otherwise
   */
  bool supportsVariableBatchSize() const override {
    return backend_.supportsVariableBatchSize();
  }

  /**
   * @brief Check if the current hardware configuration supports batch
   * processing
   *
   * @return true if batch processing is supported, false otherwise
   */
  bool isBatchProcessingSupported() const override {
    return backend_.isBatchProcessingSupported();
  }

  /**
   * @brief Get the optimal batch size for the current hardware configuration
   *
   * @return The optimal batch size
   */
  size_t getOptimalBatchSize() const override {
    return backend_.getOptimalBatchSize();
  }

  /**
   * @brief Check if the model supports automatic batching
   *
   * @return true if automatic batching is supported, false otherwise
   */
  bool supportsAutomaticBatching() const override {
    return backend_.supportsAutomaticBatching();
  }

  /**
   * @brief Enable or disable automatic batching
   *
   * @param enable true to enable automatic batching, false to disable
   */
  void setAutomaticBatching(bool enable) override {
    backend_.setAutomaticBatching(enable);
  }
};

}  // namespace metada::framework