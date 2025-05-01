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

using metada::framework::ApplicationContext;
using metada::framework::Config;
using metada::framework::Logger;
using metada::framework::Observation;

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

/**
 * @brief Test observation construction and movement
 *
 * Verifies:
 * - Configuration constructor initialization
 * - Move constructor behavior
 * - Move assignment behavior
 * - Initialization state preservation
 * - Clone functionality
 */
TEST_F(ObservationTest, ConstructionAndMovement) {
  // Test configuration constructor
  auto observation = createObservation();
  EXPECT_TRUE(observation.isInitialized());

  // Test move constructor
  auto observation_moved(std::move(observation));
  EXPECT_FALSE(observation.isInitialized());  // Source should be invalidated
  EXPECT_TRUE(observation_moved.isInitialized());

  // Test move assignment
  observation = std::move(observation_moved);
  EXPECT_FALSE(
      observation_moved.isInitialized());  // Source should be invalidated
  EXPECT_TRUE(observation.isInitialized());

  // Test clone functionality - use ON_CALL since we can't directly mock clone
  auto cloned_observation = observation.clone();
  EXPECT_TRUE(cloned_observation.isInitialized());
  EXPECT_NE(&observation.backend(),
            &cloned_observation.backend());  // Should be different instances
}

/**
 * @brief Test input/output operations
 *
 * Verifies:
 * - loadFromFile functionality
 * - saveToFile functionality
 * - Initialization state after loading
 */
TEST_F(ObservationTest, InputOutputOperations) {
  const std::string test_filename = "test_observation.dat";

  // Test loadFromFile - can't mock these directly, test the behavior
  obs1_->loadFromFile(test_filename);
  EXPECT_TRUE(obs1_->isInitialized());

  // Test saveToFile
  obs1_->saveToFile(test_filename);
}

/**
 * @brief Test quality control operations
 *
 * Verifies:
 * - applyQC functionality
 */
TEST_F(ObservationTest, QualityControl) {
  EXPECT_CALL(obs1_->backend(), applyQC()).Times(1);
  obs1_->applyQC();
}

/**
 * @brief Test comparison operations
 *
 * Verifies:
 * - equals method
 * - Equality operator (==)
 * - Inequality operator (!=)
 * - Comparing observations with different initialization states
 */
TEST_F(ObservationTest, ComparisonOperations) {
  // Test equality and inequality operators
  EXPECT_CALL(obs1_->backend(), equals(testing::Ref(obs2_->backend())))
      .WillOnce(Return(true))
      .WillOnce(Return(true))
      .WillOnce(Return(false))
      .WillOnce(Return(false));

  EXPECT_TRUE(*obs1_ == *obs2_);
  EXPECT_FALSE(*obs1_ != *obs2_);

  EXPECT_FALSE(*obs1_ == *obs2_);
  EXPECT_TRUE(*obs1_ != *obs2_);

  // Test comparison with different initialization states
  auto uninit_obs = createObservation();
  uninit_obs = Observation<Traits::ObservationType>(
      std::move(*obs1_));  // Move to invalidate obs1_

  // Reset obs1_ for further tests
  obs1_.reset(new Observation<Traits::ObservationType>(getConfig()));

  EXPECT_FALSE(*obs1_ == uninit_obs);
  EXPECT_TRUE(*obs1_ != uninit_obs);
}

/**
 * @brief Test data access and variable information
 *
 * Verifies:
 * - getData (const and non-const)
 * - getVariableNames
 * - hasVariable (for existing and non-existing variables)
 * - getDimensions
 */
TEST_F(ObservationTest, DataAccessAndInformation) {
  obs1_->backend().setVariables(variableNames_);
  obs1_->backend().setDimensions("temperature", dimensions_);
  obs1_->backend().setData(confidenceValues_);

  // Test variable names access
  const auto& vars = obs1_->getVariableNames();
  EXPECT_EQ(vars, variableNames_);

  // Test dimensions access
  const auto& dims = obs1_->getDimensions("temperature");
  EXPECT_EQ(dims, dimensions_);

  // Test data access
  const double* data = &obs1_->getData<double>();
  EXPECT_NE(data, nullptr);  // Verify we got a valid pointer
  // Verify actual data values
  for (size_t i = 0; i < confidenceValues_.size(); ++i) {
    EXPECT_DOUBLE_EQ(data[i], confidenceValues_[i]);
  }

  // Test hasVariable
  EXPECT_TRUE(obs1_->hasVariable("temperature"));
  EXPECT_FALSE(obs1_->hasVariable("nonexistent_variable"));
}

/**
 * @brief Test arithmetic operations
 *
 * Verifies:
 * - Binary operators (+, -, *)
 * - Assignment operators (+=, -=, *=)
 * - Scalar multiplication (both left and right)
 */
TEST_F(ObservationTest, ArithmeticOperations) {
  // Test addition operators
  EXPECT_CALL(obs1_->backend(), add(testing::Ref(obs2_->backend()))).Times(1);

  auto result = *obs1_ + *obs2_;  // Binary operator
  *obs1_ += *obs2_;               // Assignment operator

  // Test subtraction operators
  EXPECT_CALL(obs1_->backend(), subtract(testing::Ref(obs2_->backend())))
      .Times(1);

  result = *obs1_ - *obs2_;  // Binary operator
  *obs1_ -= *obs2_;          // Assignment operator

  // Test multiplication operators
  EXPECT_CALL(obs1_->backend(), multiply(2.0)).Times(1);

  result = *obs1_ * 2.0;  // Right scalar multiplication
  result = 2.0 * *obs1_;  // Left scalar multiplication (friend operator)
  *obs1_ *= 2.0;          // Assignment operator
}

/**
 * @brief Test backend access
 *
 * Verifies:
 * - Non-const backend access
 * - Const backend access
 */
TEST_F(ObservationTest, BackendAccess) {
  // Test non-const backend access
  auto& backend = obs1_->backend();
  EXPECT_EQ(&backend, &(obs1_->backend()));

  // Test const backend access
  const auto& const_obs = *obs1_;
  const auto& const_backend = const_obs.backend();
  EXPECT_EQ(&const_backend, &(const_obs.backend()));
}

}  // namespace metada::tests