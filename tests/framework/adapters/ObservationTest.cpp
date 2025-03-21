/**
 * @file ObservationTest.cpp
 * @brief Unit tests for the Observation class template
 * @ingroup tests
 * @author Metada Framework Team
 *
 * This test suite verifies the functionality of the Observation class template,
 * which provides a generic interface for observation implementations.
 *
 * The tests cover:
 * - Construction and initialization
 * - Core operations (reset, validate, data access)
 * - Copy/move semantics
 * - Arithmetic operations (+, -, *, +=, -=, *=)
 * - Comparison operators (==, !=)
 * - Metadata management
 * - Location and time data handling
 * - Quality control flags
 * - Error handling for invalid states
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 *
 * @see Observation
 * @see IObservation
 * @see MockObservation
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "MockConfig.hpp"
#include "MockLogger.hpp"
#include "MockObservation.hpp"
#include "MockState.hpp"
#include "Observation.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Const;
using ::testing::Return;
using ::testing::ReturnRef;

using metada::framework::Config;
using metada::framework::Logger;
using metada::framework::Observation;
using metada::framework::runs::ApplicationContext;

using Traits = AppTraits<MockLogger, MockConfig, MockState, MockObservation>;

/**
 * @brief Test fixture for Observation class tests
 *
 * Provides common setup and test data for Observation tests including:
 * - Mock objects and application context
 * - Sample observation data and metadata
 * - Helper methods for observation creation and access
 */
class ObservationTest : public ::testing::Test {
 protected:
  /** @brief Application context instance */
  std::unique_ptr<ApplicationContext<Traits>> context_;

  /** @brief Variable names for test data */
  std::vector<std::string> variableNames_;

  /** @brief Variable dimensions */
  std::vector<size_t> dimensions_;

  /** @brief Sample location coordinates */
  std::vector<std::vector<double>> locations_;

  /** @brief Sample timestamps */
  std::vector<double> times_;

  /** @brief Quality control flags */
  std::vector<int> qualityFlags_;

  /** @brief Confidence values */
  std::vector<double> confidenceValues_;

  /** @brief Sample temperature value */
  double temperature_;

  /** @brief Sample uncertainty value */
  double uncertainty_;

  /** @brief First test observation */
  std::unique_ptr<Observation<Traits::ObservationType>> obs1_;

  /** @brief Second test observation */
  std::unique_ptr<Observation<Traits::ObservationType>> obs2_;

  /**
   * @brief Set up test fixture
   *
   * Initializes test data, creates observation instances, and sets up mock
   * expectations
   */
  void SetUp() override {
    // Create a new application context for this test
    context_ = std::make_unique<ApplicationContext<Traits>>("ObservationTest");

    // Initialize test data
    variableNames_ = {"temperature", "pressure"};
    dimensions_ = {1, 1};
    locations_ = {{45.0, -120.0, 100.0}, {46.0, -121.0, 200.0}};
    times_ = {1609459200.0, 1609545600.0};  // Example timestamps (Unix time)
    qualityFlags_ = {0, 1};                 // 0 = good, 1 = suspect
    confidenceValues_ = {0.95, 0.85};
    temperature_ = 25.5;
    uncertainty_ = 0.5;

    // Create observations
    obs1_ = std::make_unique<Observation<Traits::ObservationType>>(getConfig());
    obs2_ = std::make_unique<Observation<Traits::ObservationType>>(getConfig());
  }

  /**
   * @brief Clean up test fixture
   *
   * Releases test resources and cleans up test data
   */
  void TearDown() override {
    // Clean up observations
    obs1_.reset();
    obs2_.reset();

    // Then clean up other resources
    variableNames_.clear();
    dimensions_.clear();
    locations_.clear();
    times_.clear();
    qualityFlags_.clear();
    confidenceValues_.clear();
    context_.reset();
  }

  /**
   * @brief Get logger instance from context
   * @return Reference to logger
   */
  Logger<Traits::LoggerType>& getLogger() { return context_->getLogger(); }

  /**
   * @brief Get configuration instance from context
   * @return Reference to configuration
   */
  Config<Traits::ConfigType>& getConfig() { return context_->getConfig(); }

  /**
   * @brief Create new observation instance
   * @return New observation initialized with test config
   */
  Observation<Traits::ObservationType> createObservation() {
    return Observation<Traits::ObservationType>(getConfig());
  }
};

}  // namespace metada::tests