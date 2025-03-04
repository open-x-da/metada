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
  MockConfig config;

  ObsOperator<Traits::ObsOperatorType> obsOp(config);

  // Verify initialization succeeded
  EXPECT_TRUE(obsOp.isInitialized());

  // Verify we can't initialize again
  EXPECT_THROW(obsOp.initialize(config), std::runtime_error);
}

}  // namespace metada::tests
