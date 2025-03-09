/**
 * @file MockModel.hpp
 * @brief Mock implementation of IModel interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IModel interface using Google
 * Mock. It allows testing code that depends on IModel by providing mock
 * implementations of all interface methods that can be configured with
 * expectations and behaviors.
 *
 * The mock implementation supports:
 * - Setting expectations on method calls
 * - Configuring return values and behaviors
 * - Verifying interaction patterns
 * - Testing error conditions
 *
 * @see IModel
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>

#include "IAIPredictor.hpp"
#include "IBatchProcessor.hpp"
#include "IHardwareAccelerator.hpp"
#include "IModel.hpp"
#include "utils/config/Config.hpp"

namespace metada::tests {

using framework::Config;
using framework::IAIPredictor;
using framework::IBatchProcessor;
using framework::IConfig;
using framework::IHardwareAccelerator;
using framework::IModel;
using framework::IState;
using framework::ITimeStepper;

/**
 * @brief Mock implementation of IModel for testing
 *
 * @details
 * Provides mock methods for all IModel interface operations, organized
 * into the following categories:
 *
 * @par Lifecycle Management
 * - initialize() - Initialize model from configuration
 * - reset() - Reset model to initial state
 * - isInitialized() - Check if model is initialized
 *
 * @par Parameter Management
 * - getParameter() - Get a parameter value by name
 * - setParameter() - Set a parameter value
 *
 * @par Model Execution
 * - run() - Run the model from initial to final state
 *
 * @par Capability Access
 * - getTimeStepper() - Get time stepping capability
 * - getAIPredictor() - Get AI prediction capability
 * - getBatchProcessor() - Get batch processing capability
 * - getHardwareAccelerator() - Get hardware acceleration capability
 *
 * @note All mock methods use Google Mock's MOCK_METHOD macro to enable
 * setting expectations and verifying calls.
 */
class MockModel : public IModel {
 public:
  /**
   * @brief Constructor taking a configuration object
   *
   * @param config Configuration object with model parameters
   */
  explicit MockModel(const Config<MockConfig>& config)
      : IModel(config.backend()) {}

  // Lifecycle management
  MOCK_METHOD(void, initialize, (const IConfig& config), (override));
  MOCK_METHOD(void, reset, (), (override));
  MOCK_METHOD(bool, isInitialized, (), (const, override));

  // Parameter management
  MOCK_METHOD(std::string, getParameter, (const std::string& name),
              (const, override));
  MOCK_METHOD(void, setParameter,
              (const std::string& name, const std::string& value), (override));

  // Model execution
  MOCK_METHOD(void, run,
              (const IState& initialState, IState& finalState, double startTime,
               double endTime),
              (override));

  // Capability access
  MOCK_METHOD(ITimeStepper*, getTimeStepper, (), (override));
  MOCK_METHOD(IAIPredictor*, getAIPredictor, (), (override));
  MOCK_METHOD(IBatchProcessor*, getBatchProcessor, (), (override));
  MOCK_METHOD(IHardwareAccelerator*, getHardwareAccelerator, (), (override));
};

}  // namespace metada::tests