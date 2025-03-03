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

  /** @brief Mock configuration */
  Config<MockConfig> mockConfig_;

  /** @brief Mock observation backend */
  MockObservation mockBackend_{mockConfig_};

  /** @brief Test observation instance */
  Observation<MockObservation> observation_{mockBackend_};

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

    // Setup common expectations for test data access
    ON_CALL(obs1_->backend(), getData()).WillByDefault(Return(&temperature_));
    ON_CALL(obs1_->backend(), getUncertainty())
        .WillByDefault(Return(&uncertainty_));
    ON_CALL(obs1_->backend(), getSize()).WillByDefault(Return(2));
    ON_CALL(obs1_->backend(), isValid()).WillByDefault(Return(true));
    ON_CALL(obs1_->backend(), isInitialized()).WillByDefault(Return(true));
    ON_CALL(obs1_->backend(), getLocations())
        .WillByDefault(ReturnRef(locations_));
    ON_CALL(obs1_->backend(), getTimes()).WillByDefault(ReturnRef(times_));
    ON_CALL(obs1_->backend(), getQualityFlags())
        .WillByDefault(ReturnRef(qualityFlags_));

    // Default metadata behavior - can be overridden in specific tests
    ON_CALL(obs1_->backend(), hasMetadata("key1")).WillByDefault(Return(true));
    ON_CALL(obs1_->backend(), getMetadata("key1"))
        .WillByDefault(Return("value1"));
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
 * @brief Test observation construction
 *
 * Verifies:
 * - Configuration constructor initialization
 * - Copy constructor behavior
 * - Move constructor behavior
 * - Initialization state preservation
 */
TEST_F(ObservationTest, ConstructorTests) {
  // Test configuration constructor - use the helper method
  auto observation = createObservation();
  EXPECT_TRUE(observation.isInitialized());

  // For the copy constructor test, we need to use ON_CALL instead of
  // EXPECT_CALL because we can't set expectations on state_copy.backend()
  // before it exists
  ON_CALL(observation.backend(), copyFrom(_)).WillByDefault(Return());

  // Test copy constructor
  Observation<Traits::ObservationType> observation_copy(observation);
  EXPECT_CALL(observation_copy.backend(), equals(_)).WillOnce(Return(true));
  EXPECT_TRUE(observation_copy == observation);
  EXPECT_TRUE(
      observation_copy
          .isInitialized());  // Copy should preserve initialization state

  // For the move constructor test, we need to use ON_CALL
  ON_CALL(observation_copy.backend(), moveFrom(_)).WillByDefault(Return());

  // Test move constructor
  Observation<Traits::ObservationType> observation_moved(
      std::move(observation_copy));
  // Verify the moved-from state is no longer initialized
  EXPECT_FALSE(observation_copy.isInitialized());
  EXPECT_TRUE(observation_moved.isInitialized());
}

/**
 * @brief Test core observation operations
 *
 * Verifies:
 * - reset() functionality
 * - validate() functionality
 */
TEST_F(ObservationTest, CoreObservationOperations) {
  // Use the pre-created observation object
  EXPECT_CALL(obs1_->backend(), reset()).Times(1);
  obs1_->reset();

  EXPECT_CALL(obs1_->backend(), validate()).Times(1);
  obs1_->validate();
}

/**
 * @brief Test data access
 *
 * Verifies getData() returns correct typed data
 */
TEST_F(ObservationTest, GetDataReturnsCorrectValue) {
  EXPECT_EQ(obs1_->getData<double>(), temperature_);
}

/**
 * @brief Test uncertainty access
 *
 * Verifies getUncertainty() returns correct uncertainty value
 */
TEST_F(ObservationTest, GetUncertaintyReturnsCorrectValue) {
  EXPECT_EQ(obs1_->getUncertainty<double>(), uncertainty_);
}

/**
 * @brief Test location data access
 *
 * Verifies getLocations() returns correct location coordinates
 */
TEST_F(ObservationTest, LocationsAreCorrectlyReturned) {
  // Expect the getLocations method to be called
  EXPECT_CALL(obs1_->backend(), getLocations()).Times(1);

  // Get the locations
  const auto& locs = obs1_->getLocations();

  // Verify the location data
  EXPECT_EQ(locs.size(), 2);
  EXPECT_EQ(locs[0][0], 45.0);
  EXPECT_EQ(locs[0][1], -120.0);
  EXPECT_EQ(locs[0][2], 100.0);
  EXPECT_EQ(locs[1][0], 46.0);
  EXPECT_EQ(locs[1][1], -121.0);
  EXPECT_EQ(locs[1][2], 200.0);
}

/**
 * @brief Test location setting
 *
 * Verifies setLocations() properly updates backend
 */
TEST_F(ObservationTest, SetLocationsCallsBackend) {
  // Use the pre-created observation with mockBackend_
  std::vector<std::vector<double>> newLocations{{50.0, -110.0, 200.0},
                                                {51.0, -111.0, 300.0}};

  // Set expectations
  EXPECT_CALL(obs1_->backend(),
              setLocations(testing::ElementsAreArray(newLocations)))
      .Times(1);

  // Test the method
  obs1_->setLocations(newLocations);
}

/**
 * @brief Test method chaining
 *
 * Verifies fluent interface works for multiple operations
 */
TEST_F(ObservationTest, FluentInterfaceWorks) {
  // Use the pre-created observation with mockBackend_
  std::vector<std::vector<double>> newLocations{{50.0, -110.0, 200.0},
                                                {51.0, -111.0, 300.0}};
  std::vector<double> newTimes{1234567892.0, 1234567893.0};
  std::vector<int> newFlags{3, 4};

  // Set expectations for each method call
  EXPECT_CALL(obs1_->backend(),
              setLocations(testing::ElementsAreArray(newLocations)))
      .Times(1);
  EXPECT_CALL(obs1_->backend(), setTimes(testing::ElementsAreArray(newTimes)))
      .Times(1);
  EXPECT_CALL(obs1_->backend(),
              setQualityFlags(testing::ElementsAreArray(newFlags)))
      .Times(1);

  // Test chaining of methods
  obs1_->setLocations(newLocations)
      .setTimes(newTimes)
      .setQualityFlags(newFlags);
}

/**
 * @brief Test metadata operations
 *
 * Verifies:
 * - Setting metadata
 * - Getting metadata
 * - Checking metadata existence
 */
TEST_F(ObservationTest, MetadataOperationsWork) {
  // Use the pre-created observation with mockBackend_
  const std::string key = "key1";
  const std::string value = "value1";

  // Set expectations for metadata operations
  EXPECT_CALL(obs1_->backend(), setMetadata(key, value)).Times(1);
  EXPECT_CALL(obs1_->backend(), getMetadata(key)).WillOnce(Return(value));
  EXPECT_CALL(obs1_->backend(), hasMetadata(key)).WillOnce(Return(true));

  // Test the methods
  obs1_->setMetadata(key, value);
  EXPECT_EQ(obs1_->getMetadata(key), value);
  EXPECT_TRUE(obs1_->hasMetadata(key));
}

/**
 * @brief Test invalid data access
 *
 * Verifies exception throwing for invalid observation access
 */
TEST_F(ObservationTest, ThrowsOnInvalidAccess) {
  // Create a new observation with invalid state for this specific test
  auto invalidObs = createObservation();

  // Override the default valid state for this test
  EXPECT_CALL(invalidObs.backend(), isValid()).WillRepeatedly(Return(false));

  // Test that accessing data from an invalid observation throws an exception
  EXPECT_THROW(invalidObs.getData<double>(), std::runtime_error);
}

/**
 * @brief Test copy assignment
 *
 * Verifies:
 * - Backend copyFrom() called correctly
 * - Initialization state preserved
 */
TEST_F(ObservationTest, CopyAssignment) {
  // Use the pre-created observation objects
  EXPECT_CALL(obs2_->backend(), copyFrom(testing::Ref(obs1_->backend())))
      .Times(1);
  *obs2_ = *obs1_;

  // Verify initialization of observation is preserved
  EXPECT_TRUE(obs2_->isInitialized());
}

/**
 * @brief Test move assignment
 *
 * Verifies:
 * - Proper observation transfer
 * - Source invalidation
 * - Destination validation
 */
TEST_F(ObservationTest, MoveAssignment) {
  // Create a new observation for move testing
  auto temp_observation = createObservation();
  EXPECT_TRUE(temp_observation.isInitialized());

  // Set up expectations for move assignment
  EXPECT_CALL(obs1_->backend(), moveFrom(_)).Times(1);

  // Test move assignment
  *obs1_ = std::move(temp_observation);

  // Verify observations after move
  EXPECT_FALSE(temp_observation.isInitialized());
  EXPECT_TRUE(obs1_->isInitialized());
}

/**
 * @brief Test arithmetic operations
 *
 * Verifies add(), subtract(), multiply() functionality
 */
TEST_F(ObservationTest, ArithmeticOperationsWork) {
  // Create observations for this test
  auto result = createObservation();

  // Setup expectations
  EXPECT_CALL(result.backend(), add(_)).Times(1);
  EXPECT_CALL(result.backend(), subtract(_)).Times(1);
  EXPECT_CALL(result.backend(), multiply(2.0)).Times(1);

  // Test arithmetic operations
  result.add(*obs1_);
  result.subtract(*obs1_);
  result.multiply(2.0);
}

/**
 * @brief Test operator overloads
 *
 * Verifies +, -, * operator functionality
 */
TEST_F(ObservationTest, OperatorOverloadsWork) {
  // Create a result observation for testing
  auto result = createObservation();

  // Test operator overloads
  EXPECT_NO_THROW(result = *obs1_ + *obs2_);
  EXPECT_NO_THROW(auto result1 = *obs1_ + *obs2_);
  EXPECT_NO_THROW(result = *obs1_ - *obs2_);
  EXPECT_NO_THROW(auto result2 = *obs1_ - *obs2_);
  EXPECT_NO_THROW(result = *obs1_ * 2.0);
  EXPECT_NO_THROW(auto result3 = *obs1_ * 2.0);
  EXPECT_NO_THROW(result = 2.0 * *obs1_);
  EXPECT_NO_THROW(auto result4 = 2.0 * *obs1_);
}

/**
 * @brief Test arithmetic assignment operators
 *
 * Verifies +=, -=, *= operator functionality
 */
TEST_F(ObservationTest, ArithmeticAssignmentOperatorsWork) {
  // Set expectations for the backend operations
  EXPECT_CALL(obs1_->backend(), isValid())
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(obs2_->backend(), isValid())
      .Times(2)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(obs1_->backend(), add(_)).Times(1);
  EXPECT_CALL(obs1_->backend(), subtract(_)).Times(1);
  EXPECT_CALL(obs1_->backend(), multiply(2.0)).Times(1);

  // Test the operators
  *obs1_ += *obs2_;
  *obs1_ -= *obs2_;
  *obs1_ *= 2.0;

  // No need for assertions as we're verifying the calls to the backend
}

/**
 * @brief Test comparison operators
 *
 * Verifies == and != operator functionality
 */
TEST_F(ObservationTest, ComparisonOperatorsWork) {
  // Create observations for comparison
  auto obsA = createObservation();
  auto obsB = createObservation();

  // Set up expectations for the equals method
  EXPECT_CALL(obsA.backend(), equals(_))
      .WillOnce(Return(true))
      .WillOnce(Return(true))
      .WillOnce(Return(false))
      .WillOnce(Return(false));

  // Test the operators
  EXPECT_TRUE(obsA == obsA);   // Should be equal to itself
  EXPECT_FALSE(obsA != obsA);  // Should not be unequal to itself

  EXPECT_FALSE(obsA == obsB);  // Should not be equal to a different observation
  EXPECT_TRUE(obsA != obsB);   // Should be unequal to a different observation
}

}  // namespace metada::tests