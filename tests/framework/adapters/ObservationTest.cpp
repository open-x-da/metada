/**
 * @file ObservationTest.cpp
 * @brief Unit tests for the Observation class template
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the Observation class template,
 * which provides a generic interface for observation implementations. The tests
 * cover:
 *
 * Core functionality:
 * - Initialization and construction
 * - Observation operations (reset, validate)
 * - Data access and type safety
 *
 * Advanced features:
 * - Metadata management
 * - Spatiotemporal metadata
 * - Quality control
 * - Copy/move semantics
 * - Arithmetic operations
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
 * @details
 * Provides common test data and setup/teardown for Observation tests:
 * - Variable names and dimensions
 * - Locations and timestamps
 * - Quality flags and confidence values
 * - Test data for observations and uncertainties
 */
class ObservationTest : public ::testing::Test {
 protected:
  /** @brief Application context for testing */
  std::unique_ptr<ApplicationContext<Traits>> context_;

  // Test data
  /** @brief Names of observation variables */
  std::vector<std::string> variableNames_;

  /** @brief Dimensions of observation variables */
  std::vector<size_t> dimensions_;

  /** @brief Spatial locations of observations */
  std::vector<std::vector<double>> locations_;

  /** @brief Timestamps of observations */
  std::vector<double> times_;

  /** @brief Quality flags for observations */
  std::vector<int> qualityFlags_;

  /** @brief Confidence values for observations */
  std::vector<double> confidenceValues_;

  /** @brief Sample temperature observation value */
  double temperature_;

  /** @brief Sample uncertainty value */
  double uncertainty_;

  // Mock config and observation
  /** @brief Mock configuration object */
  Config<MockConfig> mockConfig_;

  /** @brief Mock observation backend */
  MockObservation mockBackend_{mockConfig_};

  /** @brief Observation adapter using mock backend */
  Observation<MockObservation> observation_{mockBackend_};

  // Observations for tests
  /** @brief First test observation instance */
  std::unique_ptr<Observation<Traits::ObservationType>> obs1_;

  /** @brief Second test observation instance */
  std::unique_ptr<Observation<Traits::ObservationType>> obs2_;

  /**
   * @brief Set up test data before each test
   *
   * Initializes:
   * - Mock expectations for common operations
   * - Test data access
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
   * @brief Clean up test data after each test
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
   * @brief Get reference to the logger from context
   * @return Reference to the logger instance
   */
  Logger<Traits::LoggerType>& getLogger() { return context_->getLogger(); }

  /**
   * @brief Get reference to the config from context
   * @return Reference to the configuration instance
   */
  Config<Traits::ConfigType>& getConfig() { return context_->getConfig(); }

  /**
   * @brief Create a fresh Observation object for tests that need a new instance
   * @return A new Observation object initialized with the test config
   */
  Observation<Traits::ObservationType> createObservation() {
    return Observation<Traits::ObservationType>(getConfig());
  }
};

/**
 * @brief Test initialization with config
 *
 * @details
 * Verifies that the constructor properly initializes the observation
 * by calling initialize on the backend with the provided config.
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
 * @details
 * Verifies basic observation manipulation:
 * - reset() properly resets observation
 * - validate() checks observation consistency
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
 * @details
 * Verifies that getData returns the correct typed data from the backend.
 */
TEST_F(ObservationTest, GetDataReturnsCorrectValue) {
  EXPECT_EQ(obs1_->getData<double>(), temperature_);
}

/**
 * @brief Test uncertainty access
 *
 * @details
 * Verifies that getUncertainty returns the correct typed uncertainty data.
 */
TEST_F(ObservationTest, GetUncertaintyReturnsCorrectValue) {
  EXPECT_EQ(obs1_->getUncertainty<double>(), uncertainty_);
}

/**
 * @brief Test location data access
 *
 * @details
 * Verifies that getLocations returns the correct location data.
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
 * @brief Test setting locations
 *
 * @details
 * Verifies that setLocations properly calls the backend.
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
 * @brief Test fluent interface
 *
 * @details
 * Verifies that method chaining works correctly.
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
 * @details
 * Verifies that metadata can be set and retrieved correctly.
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
 * @brief Test exception on invalid access
 *
 * @details
 * Verifies that accessing data from an invalid observation throws an exception.
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
 * @details
 * Verifies proper observation copying between instances:
 * - Correct delegation to backend copyFrom()
 * - Proper observation transfer
 * - Preservation of initialization state
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
 * @details
 * Verifies proper move semantics:
 * - Correct observation transfer
 * - Source observation invalidation
 * - Destination observation validation
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
 * @details
 * Verifies that arithmetic operations work correctly.
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
 * @details
 * Verifies that operator overloads work correctly.
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
 * @details
 * Verifies that the arithmetic assignment operators (+=, -=, *=) work
 * correctly.
 */
TEST_F(ObservationTest, ArithmeticAssignmentOperatorsWork) {
  // Set expectations for the backend operations
  EXPECT_CALL(obs1_->backend(), add(testing::Ref(obs2_->backend()))).Times(1);
  EXPECT_CALL(obs1_->backend(), subtract(testing::Ref(obs2_->backend())))
      .Times(1);
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
 * @details
 * Verifies that the equality (==) and inequality (!=) operators work correctly.
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