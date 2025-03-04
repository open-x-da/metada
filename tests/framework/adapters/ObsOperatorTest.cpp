#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "MockConfig.hpp"
#include "MockIncrement.hpp"
#include "MockLogger.hpp"
#include "MockObsOperator.hpp"
#include "MockObservation.hpp"
#include "MockState.hpp"
#include "ObsOperator.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Increment;
using framework::Observation;
using framework::ObsOperator;
using framework::State;
using framework::runs::ApplicationContext;

using Traits = AppTraits<MockLogger, MockConfig, MockState, MockIncrement,
                         MockObservation, MockObsOperator>;

class ObsOperatorTest : public ::testing::Test {
 protected:
  /** @brief Application context instance */
  std::unique_ptr<ApplicationContext<Traits>> context_;

  // Test data
  std::vector<std::string> state_vars_;
  std::vector<std::string> obs_vars_;

  // Adapters
  std::unique_ptr<State<Traits::StateType>> stateAdapter_;
  std::unique_ptr<Observation<Traits::ObservationType>> obsAdapter_;
  std::unique_ptr<Increment<Traits::IncrementType>> incrementAdapter_;

  // System under test
  std::unique_ptr<ObsOperator<Traits::ObsOperatorType>> obsOperator_;

  void SetUp() override {
    context_ = std::make_unique<ApplicationContext<Traits>>("ObsOperatorTest");

    // Setup test data
    state_vars_ = {"temperature", "pressure"};
    obs_vars_ = {"radiance"};
  }

  void TearDown() override {
    context_.reset();
    state_vars_.clear();
    obs_vars_.clear();
  }
};

TEST_F(ObsOperatorTest, NonCopyableButMovable) {
  using ObsOperatorType = ObsOperator<Traits::ObsOperatorType>;

  // Verify type traits
  EXPECT_FALSE(std::is_copy_constructible_v<ObsOperatorType>);
  EXPECT_FALSE(std::is_copy_assignable_v<ObsOperatorType>);
  EXPECT_TRUE(std::is_move_constructible_v<ObsOperatorType>);
  EXPECT_TRUE(std::is_move_assignable_v<ObsOperatorType>);
}

TEST_F(ObsOperatorTest, CannotInitializeWithoutConfig) {
  using ObsOperatorType = ObsOperator<Traits::ObsOperatorType>;

  // Verify that default construction is not possible
  EXPECT_FALSE(std::is_default_constructible_v<ObsOperatorType>);

  // Verify that construction without config is not possible
  EXPECT_FALSE((std::is_constructible_v<ObsOperatorType>));
}

TEST_F(ObsOperatorTest, CanInitializeFromConfig) {
  // Create mock config
  Config<MockConfig> config;

  ObsOperator<Traits::ObsOperatorType> obsOp(config);

  // Verify initialization succeeded
  EXPECT_TRUE(obsOp.isInitialized());

  // Verify we can't initialize again
  EXPECT_THROW(obsOp.initialize(config), std::runtime_error);
}

TEST_F(ObsOperatorTest, ApplyOperations) {
  Config<MockConfig> config;
  ObsOperator<Traits::ObsOperatorType> obsOp(config);

  // Create mock state and observation
  State<Traits::StateType> state(config);
  Observation<Traits::ObservationType> obs(config);

  // Test apply operation
  EXPECT_CALL(obsOp.backend(),
              apply(testing::Ref(state.backend()), testing::Ref(obs.backend())))
      .Times(1);
  obsOp.apply(state, obs);

  // Test tangent linear operation
  Increment<Traits::IncrementType> dx(config);
  EXPECT_CALL(obsOp.backend(), applyTangentLinear(testing::Ref(dx.backend()),
                                                  testing::Ref(obs.backend())))
      .Times(1);
  obsOp.applyTangentLinear(dx, obs);

  // Test adjoint operation
  EXPECT_CALL(obsOp.backend(), applyAdjoint(testing::Ref(obs.backend()),
                                            testing::Ref(dx.backend())))
      .Times(1);
  obsOp.applyAdjoint(obs, dx);
}

TEST_F(ObsOperatorTest, RequiredVariables) {
  Config<MockConfig> config;
  ObsOperator<Traits::ObsOperatorType> obsOp(config);

  // Setup test data
  std::vector<std::string> stateVars = {"temperature", "pressure"};
  std::vector<std::string> obsVars = {"radiance"};

  // Test getRequiredStateVars
  EXPECT_CALL(obsOp.backend(), getRequiredStateVars())
      .WillOnce(ReturnRef(stateVars));
  EXPECT_EQ(obsOp.getRequiredStateVars(), stateVars);

  // Test getRequiredObsVars
  EXPECT_CALL(obsOp.backend(), getRequiredObsVars())
      .WillOnce(ReturnRef(obsVars));
  EXPECT_EQ(obsOp.getRequiredObsVars(), obsVars);
}

TEST_F(ObsOperatorTest, ThrowsWhenNotInitialized) {
  Config<MockConfig> config;
  ObsOperator<Traits::ObsOperatorType> obsOp(config);

  // Force uninitialized state
  obsOp = ObsOperator<Traits::ObsOperatorType>(std::move(obsOp));

  State<Traits::StateType> state(config);
  Observation<Traits::ObservationType> obs(config);
  Increment<Traits::IncrementType> dx(config);

  EXPECT_THROW(obsOp.apply(state, obs), std::runtime_error);
  EXPECT_THROW(obsOp.applyTangentLinear(dx, obs), std::runtime_error);
  EXPECT_THROW(obsOp.applyAdjoint(obs, dx), std::runtime_error);
}

TEST_F(ObsOperatorTest, MoveOperations) {
  Config<MockConfig> config;
  ObsOperator<Traits::ObsOperatorType> obsOp1(config);

  // Test move constructor
  ObsOperator<Traits::ObsOperatorType> obsOp2(std::move(obsOp1));
  EXPECT_TRUE(obsOp2.isInitialized());
  EXPECT_FALSE(obsOp1.isInitialized());  // NOLINT: accessing moved-from object

  // Test move assignment
  ObsOperator<Traits::ObsOperatorType> obsOp3(config);
  obsOp3 = std::move(obsOp2);
  EXPECT_TRUE(obsOp3.isInitialized());
  EXPECT_FALSE(obsOp2.isInitialized());  // NOLINT: accessing moved-from object
}

}  // namespace metada::tests
