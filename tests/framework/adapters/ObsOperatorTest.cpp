#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "Config.hpp"
#include "MockBackendTraits.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::Observation;
using framework::ObsOperator;
using framework::State;

/**
 * @brief Test fixture for ObsOperator class
 *
 * Tests the observation operator class that provides mapping
 * between model state space and observation space in data assimilation.
 * The ObsOperator transforms model state variables into simulated observations
 * that can be compared with actual observations.
 */
class ObsOperatorTest : public ::testing::Test {
 protected:
  // Test data
  std::vector<std::string> state_vars_{"temperature", "pressure"};
  std::vector<std::string> obs_vars_{"radiance"};

  // Test vectors for expected return values
  std::vector<std::string> state_vars_return_;
  std::vector<std::string> obs_vars_return_;

  std::unique_ptr<Config<traits::MockBackendTag>> config_;

  void SetUp() override {
    // Setup return values for metadata methods
    state_vars_return_ = state_vars_;
    obs_vars_return_ = obs_vars_;

    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    auto config_file = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file);
  }

  void TearDown() override { config_.reset(); }
};

/**
 * @brief Test that ObsOperator correctly follows NonCopyable pattern
 *
 * Verifies that ObsOperator adheres to the NonCopyable design pattern
 * by ensuring it cannot be copied but can be moved. This is important
 * for resource management and ownership semantics.
 */
TEST_F(ObsOperatorTest, NonCopyableButMovable) {
  using ObsOperatorType = ObsOperator<traits::MockBackendTag>;

  // Static assertions to verify type traits
  static_assert(!std::is_default_constructible_v<ObsOperatorType>,
                "ObsOperator should not be default constructible");
  static_assert(!std::is_copy_constructible_v<ObsOperatorType>,
                "ObsOperator should not be copy constructible");
  static_assert(!std::is_copy_assignable_v<ObsOperatorType>,
                "ObsOperator should not be copy assignable");
  static_assert(std::is_move_constructible_v<ObsOperatorType>,
                "ObsOperator should be move constructible");
  static_assert(std::is_move_assignable_v<ObsOperatorType>,
                "ObsOperator should be move assignable");
}

/**
 * @brief Test constructor with Config
 *
 * Verifies that ObsOperator can be properly constructed with a Config object
 * and that the backend is correctly initialized during construction.
 */
TEST_F(ObsOperatorTest, ConstructorWithConfig) {
  // Create ObsOperator using config
  ObsOperator<traits::MockBackendTag> obsOp(*config_);

  // Verify backend was initialized
  auto& backend = obsOp.backend();
  EXPECT_CALL(backend, isInitialized()).WillOnce(Return(true));
  EXPECT_TRUE(obsOp.isInitialized());
}

/**
 * @brief Test move operations
 *
 * Verifies that ObsOperator correctly implements move semantics,
 * ensuring that resources are properly transferred during move operations
 * and that the moved-from object is left in a valid but uninitialized state.
 */
TEST_F(ObsOperatorTest, MoveOperations) {
  // Create source operator
  ObsOperator<traits::MockBackendTag> obsOp1(*config_);
  auto& backend1 = obsOp1.backend();

  // Test move constructor
  ObsOperator<traits::MockBackendTag> obsOp2(std::move(obsOp1));
  auto& backend2 = obsOp2.backend();

  // Verify backend is initialized
  EXPECT_CALL(backend2, isInitialized()).WillOnce(Return(true));
  EXPECT_TRUE(obsOp2.isInitialized());
  EXPECT_CALL(backend1, isInitialized()).WillOnce(Return(false));
  EXPECT_FALSE(obsOp1.isInitialized());

  // Test move assignment
  ObsOperator<traits::MockBackendTag> obsOp3(*config_);
  obsOp3 = std::move(obsOp2);

  auto& backend3 = obsOp3.backend();
  EXPECT_CALL(backend3, isInitialized()).WillOnce(Return(true));
  EXPECT_TRUE(obsOp3.isInitialized());
  EXPECT_CALL(backend2, isInitialized()).WillOnce(Return(false));
  EXPECT_FALSE(obsOp2.isInitialized());
}

/**
 * @brief Test required variables access
 *
 * Verifies that ObsOperator correctly exposes the required state and
 * observation variables by delegating to the backend implementation.
 * These methods are essential for determining compatibility between
 * State and Observation objects with this operator.
 */
TEST_F(ObsOperatorTest, RequiredVariables) {
  // Create ObsOperator using config
  ObsOperator<traits::MockBackendTag> obsOp(*config_);
  auto& backend = obsOp.backend();

  // Verify required variables access
  EXPECT_CALL(backend, getRequiredStateVars())
      .WillOnce(ReturnRef(state_vars_return_));
  EXPECT_CALL(backend, getRequiredObsVars())
      .WillOnce(ReturnRef(obs_vars_return_));

  // Test required variables accessors
  const auto& state_vars = obsOp.getRequiredStateVars();
  const auto& obs_vars = obsOp.getRequiredObsVars();

  EXPECT_EQ(state_vars, state_vars_);
  EXPECT_EQ(obs_vars, obs_vars_);
}

/**
 * @brief Test apply method
 *
 * Verifies that the apply method correctly delegates to the backend
 * implementation. This method transforms a State object into simulated
 * observations that are stored in the provided Observation object.
 */
TEST_F(ObsOperatorTest, Apply) {
  // Create ObsOperator using config
  ObsOperator<traits::MockBackendTag> obsOp(*config_);
  auto& backend = obsOp.backend();

  State<traits::MockBackendTag> state(*config_);
  Observation<traits::MockBackendTag> obs(*config_);

  EXPECT_CALL(backend, apply(_, _)).Times(1);
  obsOp.apply(state, obs);
}

}  // namespace metada::tests
