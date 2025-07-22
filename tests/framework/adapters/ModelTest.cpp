/**
 * @file ModelTest.cpp
 * @brief Unit tests for the Model class template
 * @ingroup tests
 * @author Metada Framework Team
 *
 * This test suite verifies the functionality of the Model class template,
 * which provides a generic interface for model implementations.
 *
 * The tests cover:
 * - Construction and initialization
 * - Parameter management
 * - Model execution
 * - Error handling
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>

#include "Config.hpp"
#include "Logger.hpp"
#include "MockBackendTraits.hpp"
#include "Model.hpp"
#include "State.hpp"

namespace metada::tests {

using namespace metada::framework;

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Eq;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::Throw;

using framework::Config;
using framework::Geometry;
using framework::Model;
using framework::State;

/**
 * @brief Test fixture for Model class template
 */
class ModelTest : public ::testing::Test {
 protected:
  std::string config_file_;
  std::unique_ptr<Config<traits::MockBackendTag>> config_;
  std::unique_ptr<Geometry<traits::MockBackendTag>> geometry_;

  // Test objects
  std::unique_ptr<Model<traits::MockBackendTag>> mockModel_;

  void SetUp() override {
    // Create and configure the mock model
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file_);

    // Initialize the Logger singleton before creating any objects that use it
    Logger<traits::MockBackendTag>::Init(*config_);

    geometry_ = std::make_unique<Geometry<traits::MockBackendTag>>(*config_);
    mockModel_ = std::make_unique<Model<traits::MockBackendTag>>(*config_);
  }

  void TearDown() override {
    // Clean up
    config_.reset();
    geometry_.reset();
    mockModel_.reset();
    // Reset the Logger singleton after tests
    Logger<traits::MockBackendTag>::Reset();
  }
};

/**
 * @brief Test construction and initialization
 */
TEST_F(ModelTest, ConstructAndInitialize) {
  // Create model
  Model<traits::MockBackendTag> model(*config_);

  // Verify the model is not initialized yet
  EXPECT_FALSE(model.isInitialized());

  // Setup expectations
  EXPECT_CALL(model.backend(), initialize(_));

  // Initialize the model
  model.initialize(*config_);

  // Setup expectations
  EXPECT_CALL(model.backend(), isInitialized()).WillOnce(Return(true));

  // Verify the model is now initialized
  EXPECT_TRUE(model.isInitialized());
}

/**
 * @brief Test model execution
 */
TEST_F(ModelTest, ModelExecution) {
  // Create model and mock states
  Model<traits::MockBackendTag> model(*config_);
  State<traits::MockBackendTag> initialState(*config_, *geometry_);
  State<traits::MockBackendTag> finalState(*config_, *geometry_);

  // Setup expectations for initialization
  EXPECT_CALL(model.backend(), initialize(_));

  // Initialize the model
  model.initialize(*config_);

  // Setup expectations for isInitialized
  EXPECT_CALL(model.backend(), isInitialized()).WillOnce(Return(true));

  // Setup expectations for run
  EXPECT_CALL(model.backend(), run(_, _));

  // Run the model
  model.run(initialState, finalState);
}

/**
 * @brief Test error handling on initialization
 */
TEST_F(ModelTest, InitializationError) {
  // Create model
  Model<traits::MockBackendTag> model(*config_);

  // Setup expectations to throw on initialize
  EXPECT_CALL(model.backend(), initialize(_))
      .WillOnce(Throw(std::runtime_error("Initialization error")));

  // Expect exception on initialize
  EXPECT_THROW(model.initialize(*config_), std::runtime_error);
  EXPECT_FALSE(model.isInitialized());
}

/**
 * @brief Test error handling on run
 */
TEST_F(ModelTest, RunError) {
  // Create model and mock states
  Model<traits::MockBackendTag> model(*config_);
  State<traits::MockBackendTag> initialState(*config_, *geometry_);
  State<traits::MockBackendTag> finalState(*config_, *geometry_);

  // Initialize the model
  EXPECT_CALL(model.backend(), initialize(_));
  model.initialize(*config_);

  // Setup expectations for isInitialized
  EXPECT_CALL(model.backend(), isInitialized()).WillOnce(Return(true));

  // Setup expectations to throw on run
  EXPECT_CALL(model.backend(), run(_, _))
      .WillOnce(Throw(std::runtime_error("Run error")));

  // Expect exception on run
  EXPECT_THROW(model.run(initialState, finalState), std::runtime_error);
}

/**
 * @brief Test model finalization
 */
TEST_F(ModelTest, Finalization) {
  // Create model
  Model<traits::MockBackendTag> model(*config_);

  // First initialize the model
  EXPECT_CALL(model.backend(), initialize(_));
  model.initialize(*config_);

  // Setup expectations for isInitialized
  EXPECT_CALL(model.backend(), isInitialized())
      .WillOnce(Return(true))
      .WillOnce(Return(true));

  // Verify the model is now initialized
  EXPECT_TRUE(model.isInitialized());

  // Setup expectations for finalize
  EXPECT_CALL(model.backend(), finalize());

  // Finalize the model
  model.finalize();

  // Setup expectations for isInitialized
  EXPECT_CALL(model.backend(), isInitialized()).WillOnce(Return(false));

  // Verify the model is no longer initialized
  EXPECT_FALSE(model.isInitialized());
}

/**
 * @brief Test error handling during finalization
 */
TEST_F(ModelTest, FinalizationError) {
  // Create model
  Model<traits::MockBackendTag> model(*config_);

  // First initialize the model
  EXPECT_CALL(model.backend(), initialize(_));
  model.initialize(*config_);

  // Setup expectations for isInitialized
  EXPECT_CALL(model.backend(), isInitialized()).WillOnce(Return(true));

  // Setup expectations for finalize to throw
  EXPECT_CALL(model.backend(), finalize())
      .WillOnce(Throw(std::runtime_error("Finalization error")));

  // Expect exception on finalize
  EXPECT_THROW(model.finalize(), std::runtime_error);
}

}  // namespace metada::tests