#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "Increment.hpp"
#include "MockConfig.hpp"
#include "MockLogger.hpp"
#include "MockObsOperator.hpp"
#include "MockObservation.hpp"
#include "MockState.hpp"
#include "ObsOperator.hpp"
#include "State.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::Increment;
using framework::Logger;
using framework::Observation;
using framework::ObsOperator;
using framework::State;
using framework::runs::ApplicationContext;

// Define test traits without MockIncrement
using Traits = AppTraits<MockLogger, MockConfig, MockState, MockObservation,
                         MockObsOperator>;

/**
 * @brief Test fixture for ObsOperator class
 *
 * Tests the observation operator class that provides mapping
 * between model space and observation space in data assimilation.
 */
class ObsOperatorTest : public ::testing::Test {
 protected:
  /** @brief Application context instance */
  std::unique_ptr<ApplicationContext<Traits>> context_;

  // Test data
  std::vector<std::string> state_vars_{"temperature", "pressure"};
  std::vector<std::string> obs_vars_{"radiance"};

  // Mock configuration
  Config<MockConfig> mock_config_;

  // Test vectors for expected return values
  std::vector<std::string> state_vars_return_;
  std::vector<std::string> obs_vars_return_;

  void SetUp() override {
    context_ = std::make_unique<ApplicationContext<Traits>>("ObsOperatorTest");

    // Setup return values for metadata methods
    state_vars_return_ = state_vars_;
    obs_vars_return_ = obs_vars_;
  }

  void TearDown() override { context_.reset(); }

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
};

/**
 * @brief Test that ObsOperator correctly follows NonCopyable pattern
 */
TEST_F(ObsOperatorTest, NonCopyableButMovable) {
  using ObsOperatorType = ObsOperator<Traits::ObsOperatorType>;

  // Verify type traits
  EXPECT_FALSE(std::is_copy_constructible_v<ObsOperatorType>);
  EXPECT_FALSE(std::is_copy_assignable_v<ObsOperatorType>);
  EXPECT_TRUE(std::is_move_constructible_v<ObsOperatorType>);
  EXPECT_TRUE(std::is_move_assignable_v<ObsOperatorType>);
}

/**
 * @brief Test ObsOperator initialization constraints
 */
TEST_F(ObsOperatorTest, CannotInitializeWithoutConfig) {
  using ObsOperatorType = ObsOperator<Traits::ObsOperatorType>;

  // Verify that default construction is not possible
  EXPECT_FALSE(std::is_default_constructible_v<ObsOperatorType>);
}

/**
 * @brief Test initialization from configuration
 */
TEST_F(ObsOperatorTest, CanInitializeFromConfig) {
  // Setup expectations for initialization
  EXPECT_CALL(getConfig().backend(), Get("obs_operator"))
      .WillOnce(Return(framework::ConfigValue(true)));

  // Create ObsOperator using config
  ObsOperator<Traits::ObsOperatorType> obsOp(getConfig());

  // Verify initialization succeeded
  EXPECT_TRUE(obsOp.isInitialized());

  // Verify we can't initialize again
  EXPECT_THROW(obsOp.initialize(getConfig()), std::runtime_error);
}

/**
 * @brief Test initialization failure handling
 */
TEST_F(ObsOperatorTest, InitializationFailsWithInvalidConfig) {
  // Setup expectations for failed initialization
  EXPECT_CALL(getConfig().backend(), Get("obs_operator"))
      .WillOnce(Return(framework::ConfigValue(false)));

  // Expect exception when constructing with invalid config
  EXPECT_THROW(
      { ObsOperator<Traits::ObsOperatorType> obsOp(getConfig()); },
      std::runtime_error);
}

/**
 * @brief Test forward operator
 */
TEST_F(ObsOperatorTest, ApplyForwardOperator) {
  // Create test instances
  ObsOperator<Traits::ObsOperatorType> obsOp(getConfig());
  State<Traits::StateType> state(getConfig());
  Observation<Traits::ObservationType> obs(getConfig());

  // Test forward operator (H)
  EXPECT_CALL(obsOp.backend(), apply(_, _)).Times(1);
  obsOp.apply(state, obs);
}

/**
 * @brief Test tangent linear operator - using State differences
 */
TEST_F(ObsOperatorTest, ApplyTangentLinearOperator) {
  // Create test instances
  ObsOperator<Traits::ObsOperatorType> obsOp(getConfig());
  State<Traits::StateType> state1(getConfig());
  State<Traits::StateType> state2(getConfig());
  Observation<Traits::ObservationType> dy(getConfig());

  // Create increments from states
  using IncrementType = Increment<State<Traits::StateType>>;
  IncrementType dx = IncrementType::createFromDifference(state1, state2);
  // Test tangent linear operator (H')
  EXPECT_CALL(obsOp.backend(), applyTangentLinear(_, _)).Times(1);

  obsOp.applyTangentLinear(dx, dy);
}

/**
 * @brief Test adjoint operator - using State differences
 */
TEST_F(ObsOperatorTest, ApplyAdjointOperator) {
  // Create test instances
  ObsOperator<Traits::ObsOperatorType> obsOp(getConfig());
  State<Traits::StateType> state1(getConfig());
  State<Traits::StateType> state2(getConfig());
  Observation<Traits::ObservationType> dy(getConfig());

  // Create increments from states
  using IncrementType = Increment<State<Traits::StateType>>;
  IncrementType dx = IncrementType::createFromDifference(state1, state2);

  // Test adjoint operator (H'*)
  EXPECT_CALL(obsOp.backend(), applyAdjoint(_, _)).Times(1);

  obsOp.applyAdjoint(dy, dx);
}

/**
 * @brief Test accessor methods for required variables
 */
TEST_F(ObsOperatorTest, RequiredVariables) {
  ObsOperator<Traits::ObsOperatorType> obsOp(getConfig());

  // Setup expectations for required variables
  EXPECT_CALL(obsOp.backend(), getRequiredStateVars())
      .WillOnce(ReturnRef(state_vars_return_));
  EXPECT_CALL(obsOp.backend(), getRequiredObsVars())
      .WillOnce(ReturnRef(obs_vars_return_));

  // Test getters
  EXPECT_EQ(obsOp.getRequiredStateVars(), state_vars_);
  EXPECT_EQ(obsOp.getRequiredObsVars(), obs_vars_);
}

/**
 * @brief Test exception handling for uninitialized state
 */
TEST_F(ObsOperatorTest, ThrowsWhenNotInitialized) {
  // Create test instances
  ObsOperator<Traits::ObsOperatorType> obsOp1(getConfig());
  State<Traits::StateType> state(getConfig());
  Observation<Traits::ObservationType> obs(getConfig());
  State<Traits::StateType> state1(getConfig());
  State<Traits::StateType> state2(getConfig());

  // Create increment from states
  using IncrementType = Increment<State<Traits::StateType>>;
  IncrementType dx = IncrementType::createFromDifference(state1, state2);

  // Create uninitialized operator through move
  ObsOperator<Traits::ObsOperatorType> uninitializedOp(std::move(obsOp1));

  // Test exception throwing
  EXPECT_THROW(uninitializedOp.apply(state, obs), std::runtime_error);
  EXPECT_THROW(uninitializedOp.applyTangentLinear(dx, obs), std::runtime_error);
  EXPECT_THROW(uninitializedOp.applyAdjoint(obs, dx), std::runtime_error);
}

/**
 * @brief Test move operations
 */
TEST_F(ObsOperatorTest, MoveOperations) {
  // Create source operator
  ObsOperator<Traits::ObsOperatorType> obsOp1(getConfig());
  EXPECT_TRUE(obsOp1.isInitialized());

  // Test move constructor
  ObsOperator<Traits::ObsOperatorType> obsOp2(std::move(obsOp1));
  EXPECT_TRUE(obsOp2.isInitialized());
  EXPECT_FALSE(obsOp1.isInitialized());  // NOLINT: accessing moved-from object

  // Test move assignment
  ObsOperator<Traits::ObsOperatorType> obsOp3(getConfig());
  obsOp3 = std::move(obsOp2);
  EXPECT_TRUE(obsOp3.isInitialized());
  EXPECT_FALSE(obsOp2.isInitialized());  // NOLINT: accessing moved-from object
}

}  // namespace metada::tests
