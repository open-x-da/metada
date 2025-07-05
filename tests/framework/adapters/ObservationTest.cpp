/**
 * @file ObservationTest.cpp
 * @brief Unit tests for the Observation adapter class template
 * @ingroup tests
 * @author Metada Framework Team
 *
 * This test suite verifies the functionality of the Observation adapter class
 * template, which provides a generic interface for observation implementations
 * in data assimilation systems.
 *
 * The tests cover:
 * - Construction and initialization
 * - Move semantics and cloning
 * - Data access and iteration capabilities
 * - Arithmetic operations (+, -, *, +=, -=, *=)
 * - Comparison operators (==, !=)
 * - Input/output operations (file loading and saving)
 * - Quality control application
 * - Backend access (const and non-const)
 * - Geographic filtering operations
 * - Point-based observations with location information
 *
 * The test suite uses Google Test/Mock framework for mocking backend
 * implementations and verifying adapter behavior through assertions.
 *
 * @see framework::Observation
 * @see traits::BackendTraits
 * @see traits::MockBackendTag
 * @see backends::gmock::MockObservation
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <vector>

#include "Config.hpp"
#include "Location.hpp"
#include "MockBackendTraits.hpp"
#include "MockObservation.hpp"
#include "MockObservationIterator.hpp"
#include "Observation.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Const;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::CoordinateSystem;
using framework::Location;
using framework::Observation;

/**
 * @brief Test fixture for Observation class tests
 *
 * Provides common setup and test data for Observation tests including:
 * - Mock objects and application context
 * - Sample observation data with location information
 * - Helper methods for observation creation and access
 */
class ObservationTest : public ::testing::Test {
 protected:
  /** @brief Sample observation locations */
  std::vector<std::pair<double, double>> locations_;  // lat, lon pairs

  /** @brief Sample vertical levels */
  std::vector<double> levels_;

  /** @brief Sample observation values */
  std::vector<double> values_;

  /** @brief Sample observation errors */
  std::vector<double> errors_;

  /** @brief First test observation */
  std::unique_ptr<Observation<traits::MockBackendTag>> obs1_;

  /** @brief Second test observation */
  std::unique_ptr<Observation<traits::MockBackendTag>> obs2_;

  /** @brief Configuration instance */
  std::unique_ptr<Config<traits::MockBackendTag>> config_;

  /**
   * @brief Set up test fixture
   *
   * Initializes test data, creates observation instances, and sets up mock
   * expectations
   */
  void SetUp() override {
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    auto config_file = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file);
    obs1_ = std::make_unique<Observation<traits::MockBackendTag>>(*config_);
    obs2_ = std::make_unique<Observation<traits::MockBackendTag>>(*config_);

    // Initialize test data
    initializeTestData();
  }

  /**
   * @brief Clean up test fixture
   *
   * Releases test resources and cleans up test data
   */
  void TearDown() override {
    // Clean up observations
    resetObservations();

    // Then clean up other resources
    locations_.clear();
    levels_.clear();
    values_.clear();
    errors_.clear();
    config_.reset();
  }

  /**
   * @brief Create new observation instance
   * @return New observation initialized with test config
   */
  Observation<traits::MockBackendTag> createObservation() {
    return Observation<traits::MockBackendTag>(*config_);
  }

  // Helper function to initialize test data
  void initializeTestData() {
    // Sample locations (lat, lon)
    locations_ = {{45.0, -120.0}, {46.0, -121.0}, {47.0, -122.0}};

    // Sample vertical levels (pressure in hPa)
    levels_ = {1000.0, 850.0, 500.0};

    // Sample observation values
    values_ = {25.5, 15.2, -5.8};

    // Sample observation errors
    errors_ = {0.5, 0.3, 0.7};
  }

  // Helper function to create observations
  void createObservations() {
    obs1_ = std::make_unique<Observation<traits::MockBackendTag>>(*config_);
    obs2_ = std::make_unique<Observation<traits::MockBackendTag>>(*config_);
  }

  // Helper function to reset observations
  void resetObservations() {
    obs1_.reset();
    obs2_.reset();
  }

  // Helper function to verify data access
  void verifyDataAccess() {
    const auto data = obs1_->getData<std::vector<double>>();
    EXPECT_FALSE(data.empty());  // Verify we got non-empty data
    // Verify actual data values
    for (size_t i = 0; i < data.size() && i < values_.size(); ++i) {
      EXPECT_DOUBLE_EQ(data[i], values_[i]);
    }
  }

  // Helper function to verify iteration
  void verifyIteration() {
    size_t count = 0;
    for (const auto& obs : *obs1_) {
      EXPECT_TRUE(obs.is_valid);
      if (obs.location.getCoordinateSystem() == CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs.location.getGeographicCoords();
        EXPECT_GE(lat, -90.0);
        EXPECT_LE(lat, 90.0);
        EXPECT_GE(lon, -180.0);
        EXPECT_LE(lon, 180.0);
      }
      count++;
    }
    EXPECT_EQ(count, obs1_->size());
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
  const double error = 0.1;
  const double missing_value = -999.0;

  // Test loadFromFile
  EXPECT_CALL(obs1_->backend(),
              loadFromFile(test_filename, error, missing_value))
      .Times(1);
  obs1_->loadFromFile(test_filename, error, missing_value);
  EXPECT_TRUE(obs1_->isInitialized());

  // Test saveToFile
  EXPECT_CALL(obs1_->backend(), saveToFile(test_filename)).Times(1);
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
  uninit_obs = Observation<traits::MockBackendTag>(
      std::move(*obs1_));  // Move to invalidate obs1_

  // Reset obs1_ for further tests
  obs1_.reset(new Observation<traits::MockBackendTag>(*config_));

  EXPECT_FALSE(*obs1_ == uninit_obs);
  EXPECT_TRUE(*obs1_ != uninit_obs);
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