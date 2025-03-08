/**
 * @file IBatchProcessor.hpp
 * @brief Interface defining the capability for batch processing in models
 * @ingroup repr
 *
 * @details
 * This header provides the capability interface for models that support
 * batch processing, allowing them to process multiple states in a single
 * operation. This is particularly useful for AI models, but can also be
 * used by traditional models that support parallel execution.
 */

#pragma once

#include <memory>
#include <vector>

namespace metada::framework {

// Forward declarations
class IState;

/**
 * @brief Capability interface for batch processing
 *
 * @details
 * This interface defines the capabilities required for models that
 * support batch processing. Models that implement this capability
 * can process multiple states in a single operation, which can
 * significantly improve performance.
 *
 * This capability is typically used by:
 * - AI-based models that can process batches efficiently
 * - Traditional models that support parallel execution
 * - Any model that can benefit from vectorized operations
 *
 * Implementations should ensure:
 * - Efficient memory usage for batched operations
 * - Proper parallelization and hardware utilization
 * - Consistent results compared to single-state processing
 */
class IBatchProcessor {
 public:
  /**
   * @brief Virtual destructor
   */
  virtual ~IBatchProcessor() = default;

  /**
   * @brief Get the maximum batch size supported by the model
   *
   * @return The maximum batch size
   */
  virtual size_t getMaxBatchSize() const = 0;

  /**
   * @brief Set the batch size for the model
   *
   * @param batchSize The batch size to use
   * @throws std::invalid_argument if the batch size is invalid
   */
  virtual void setBatchSize(size_t batchSize) = 0;

  /**
   * @brief Get the current batch size
   *
   * @return The current batch size
   */
  virtual size_t getBatchSize() const = 0;

  /**
   * @brief Process a batch of states at once
   *
   * @param inputBatch Vector of input states
   * @param outputBatch Vector of output states (output parameter)
   * @param startTime The start time for the model run
   * @param endTime The end time for the model run
   * @throws std::runtime_error if batch processing fails
   */
  virtual void processBatch(
      const std::vector<std::shared_ptr<IState>>& inputBatch,
      std::vector<std::shared_ptr<IState>>& outputBatch, double startTime,
      double endTime) = 0;

  /**
   * @brief Check if the model supports variable batch sizes
   *
   * @return true if the model supports variable batch sizes, false otherwise
   */
  virtual bool supportsVariableBatchSize() const = 0;

  /**
   * @brief Check if the current hardware configuration supports batch
   * processing
   *
   * @return true if batch processing is supported, false otherwise
   */
  virtual bool isBatchProcessingSupported() const = 0;

  /**
   * @brief Get the optimal batch size for the current hardware configuration
   *
   * @return The optimal batch size
   */
  virtual size_t getOptimalBatchSize() const = 0;

  /**
   * @brief Check if the model supports automatic batching
   *
   * Automatic batching means the model can automatically group individual
   * requests into batches for improved performance.
   *
   * @return true if automatic batching is supported, false otherwise
   */
  virtual bool supportsAutomaticBatching() const = 0;

  /**
   * @brief Enable or disable automatic batching
   *
   * @param enable true to enable automatic batching, false to disable
   */
  virtual void setAutomaticBatching(bool enable) = 0;
};

}  // namespace metada::framework