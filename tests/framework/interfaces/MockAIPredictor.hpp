/**
 * @file MockAIPredictor.hpp
 * @brief Mock implementation of IAIPredictor interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IAIPredictor interface using
 * Google Mock. It allows testing code that depends on IAIPredictor by providing
 * mock implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * @see IAIPredictor
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <memory>
#include <string>
#include <vector>

#include "IAIPredictor.hpp"

namespace metada::tests {

using framework::IAIPredictor;
using framework::IState;

/**
 * @brief Mock implementation of IAIPredictor for testing
 *
 * Provides mock methods for all AI prediction operations including:
 * - Weight management
 * - Device configuration
 * - Prediction methods
 * - Model metadata
 */
class MockAIPredictor : public IAIPredictor {
 public:
  // Weight management
  MOCK_METHOD(void, loadWeights, (const std::string& weightsPath), (override));

  // Device configuration
  MOCK_METHOD(void, setDevice, (const std::string& device), (override));
  MOCK_METHOD(std::string, getDevice, (), (const, override));

  // Prediction methods
  MOCK_METHOD(void, predict,
              (const IState& initialState, IState& prediction, double leadTime),
              (override));

  MOCK_METHOD(void, predictMultiple,
              (const IState& initialState,
               std::vector<std::unique_ptr<IState>>& predictions,
               const std::vector<double>& leadTimes),
              (override));

  // Model capabilities
  MOCK_METHOD(bool, supportsLeadTime, (double leadTime), (const, override));
  MOCK_METHOD(std::vector<double>, getSupportedLeadTimes, (),
              (const, override));

  // Model metadata
  MOCK_METHOD(std::string, getArchitecture, (), (const, override));
  MOCK_METHOD(std::string, getWeightsVersion, (), (const, override));
  MOCK_METHOD(bool, isDeterministic, (), (const, override));
  MOCK_METHOD(bool, isReady, (), (const, override));
};

}  // namespace metada::tests