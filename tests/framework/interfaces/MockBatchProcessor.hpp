/**
 * @file MockBatchProcessor.hpp
 * @brief Mock implementation of IBatchProcessor interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IBatchProcessor interface
 * using Google Mock. It allows testing code that depends on IBatchProcessor by
 * providing mock implementations of all interface methods that can be
 * configured with expectations and behaviors.
 *
 * @see IBatchProcessor
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <memory>
#include <vector>

#include "IBatchProcessor.hpp"

namespace metada::tests {

using framework::IBatchProcessor;
using framework::IState;

/**
 * @brief Mock implementation of IBatchProcessor for testing
 *
 * Provides mock methods for all batch processing operations including:
 * - Batch size configuration
 * - Batch processing execution
 * - Batch processing capabilities
 */
class MockBatchProcessor : public IBatchProcessor {
 public:
  // Batch size configuration
  MOCK_METHOD(size_t, getMaxBatchSize, (), (const, override));
  MOCK_METHOD(void, setBatchSize, (size_t batchSize), (override));
  MOCK_METHOD(size_t, getBatchSize, (), (const, override));

  // Batch processing execution
  MOCK_METHOD(void, processBatch,
              (const std::vector<std::shared_ptr<IState>>& inputBatch,
               std::vector<std::shared_ptr<IState>>& outputBatch,
               double startTime, double endTime),
              (override));

  // Batch processing capabilities
  MOCK_METHOD(bool, supportsVariableBatchSize, (), (const, override));
  MOCK_METHOD(bool, isBatchProcessingSupported, (), (const, override));
  MOCK_METHOD(size_t, getOptimalBatchSize, (), (const, override));
  MOCK_METHOD(bool, supportsAutomaticBatching, (), (const, override));
  MOCK_METHOD(void, setAutomaticBatching, (bool enable), (override));
};

}  // namespace metada::tests