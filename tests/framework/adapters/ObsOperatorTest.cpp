#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "MockIncrement.hpp"
#include "MockObsOperator.hpp"
#include "MockObservation.hpp"
#include "MockState.hpp"
#include "ObsOperator.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Increment;
using framework::Observation;
using framework::ObsOperator;
using framework::State;

class ObsOperatorTest : public ::testing::Test {
 protected:
  // Test data
  std::vector<std::string> stateVars_{"temperature", "pressure", "humidity"};
  std::vector<std::string> obsVars_{"satellite_radiance", "brightness_temp"};

  void SetUp() override {
    // Setup mock state
    EXPECT_CALL(mockState_, isInitialized()).WillRepeatedly(Return(true));
    EXPECT_CALL(mockState_, hasVariable(_)).WillRepeatedly(Return(true));

    // Setup mock observation
    EXPECT_CALL(mockObs_, isValid()).WillRepeatedly(Return(true));
    EXPECT_CALL(mockObs_, hasAttribute(_)).WillRepeatedly(Return(true));
    EXPECT_CALL(mockObs_, getAttribute(_))
        .WillRepeatedly(Return("satellite_radiance"));

    // Setup mock increment
    EXPECT_CALL(mockIncrement_, isInitialized()).WillRepeatedly(Return(true));

    // Setup mock operator
    EXPECT_CALL(mockOperator_, getRequiredStateVariables())
        .WillRepeatedly(ReturnRef(stateVars_));
    EXPECT_CALL(mockOperator_, getRequiredObsVariables())
        .WillRepeatedly(ReturnRef(obsVars_));
    EXPECT_CALL(mockOperator_, isInitialized()).WillRepeatedly(Return(true));

    // Create adapters with mocks
    stateAdapter_ = std::make_unique<State<MockState>>(mockState_);
    obsAdapter_ = std::make_unique<Observation<MockObservation>>(mockObs_);
    incrementAdapter_ =
        std::make_unique<Increment<MockIncrement>>(mockIncrement_);

    // Create operator with mock backend
    obsOperator_ = std::make_unique<
        ObsOperator<MockObsOperator, MockState, MockObservation>>(
        mockOperator_);
  }

  // Mocks
  MockState mockState_;
  MockObservation mockObs_;
  MockIncrement mockIncrement_;
  MockObsOperator mockOperator_;

  // Adapters
  std::unique_ptr<State<MockState>> stateAdapter_;
  std::unique_ptr<Observation<MockObservation>> obsAdapter_;
  std::unique_ptr<Increment<MockIncrement>> incrementAdapter_;

  // System under test
  std::unique_ptr<ObsOperator<MockObsOperator, MockState, MockObservation>>
      obsOperator_;
};

TEST_F(ObsOperatorTest, InitializeCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, initialize()).Times(1);

  // Call method under test
  obsOperator_->initialize();
}

TEST_F(ObsOperatorTest, FinalizeCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, finalize()).Times(1);

  // Call method under test
  obsOperator_->finalize();
}

TEST_F(ObsOperatorTest, ApplyCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, apply(_, _)).Times(1);

  // Call method under test
  obsOperator_->apply(*stateAdapter_, *obsAdapter_);
}

TEST_F(ObsOperatorTest, ApplyTangentLinearCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, applyTangentLinear(_, _)).Times(1);

  // Call method under test
  obsOperator_->applyTangentLinear(*incrementAdapter_, *obsAdapter_);
}

TEST_F(ObsOperatorTest, ApplyAdjointCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, applyAdjoint(_, _)).Times(1);

  // Call method under test
  obsOperator_->applyAdjoint(*obsAdapter_, *incrementAdapter_);
}

TEST_F(ObsOperatorTest, SetObservationErrorCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, setObservationError(_)).Times(1);

  // Call method under test
  obsOperator_->setObservationError(*obsAdapter_);
}

TEST_F(ObsOperatorTest, GetObservationErrorCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, getObservationError(_)).WillOnce(Return(0.5));

  // Call method under test
  double error = obsOperator_->getObservationError(*obsAdapter_);

  // Verify result
  ASSERT_DOUBLE_EQ(error, 0.5);
}

TEST_F(ObsOperatorTest, SetParameterCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, setParameter("bias_correction", 0.1)).Times(1);

  // Call method under test
  obsOperator_->setParameter("bias_correction", 0.1);
}

TEST_F(ObsOperatorTest, GetParameterCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, getParameter("bias_correction"))
      .WillOnce(Return(0.1));

  // Call method under test
  double value = obsOperator_->getParameter("bias_correction");

  // Verify result
  ASSERT_DOUBLE_EQ(value, 0.1);
}

TEST_F(ObsOperatorTest, GetRequiredStateVariablesCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, getRequiredStateVariables())
      .WillOnce(ReturnRef(stateVars_));

  // Call method under test
  const auto& vars = obsOperator_->getRequiredStateVariables();

  // Verify result
  ASSERT_EQ(vars.size(), 3);
  ASSERT_EQ(vars[0], "temperature");
  ASSERT_EQ(vars[1], "pressure");
  ASSERT_EQ(vars[2], "humidity");
}

TEST_F(ObsOperatorTest, GetRequiredObsVariablesCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, getRequiredObsVariables())
      .WillOnce(ReturnRef(obsVars_));

  // Call method under test
  const auto& vars = obsOperator_->getRequiredObsVariables();

  // Verify result
  ASSERT_EQ(vars.size(), 2);
  ASSERT_EQ(vars[0], "satellite_radiance");
  ASSERT_EQ(vars[1], "brightness_temp");
}

TEST_F(ObsOperatorTest, IsInitializedCallsBackend) {
  // Setup expectation
  EXPECT_CALL(mockOperator_, isInitialized()).WillOnce(Return(true));

  // Call method under test
  bool initialized = obsOperator_->isInitialized();

  // Verify result
  ASSERT_TRUE(initialized);
}

TEST_F(ObsOperatorTest, ValidateCompatibilityChecksStateAndObs) {
  // Test successful validation
  EXPECT_TRUE(
      obsOperator_->validateCompatibility(*stateAdapter_, *obsAdapter_));

  // Test failed validation with missing state variable
  EXPECT_CALL(mockState_, hasVariable("temperature")).WillOnce(Return(false));
  EXPECT_FALSE(
      obsOperator_->validateCompatibility(*stateAdapter_, *obsAdapter_));

  // Reset expectations
  testing::Mock::VerifyAndClearExpectations(&mockState_);
  EXPECT_CALL(mockState_, hasVariable(_)).WillRepeatedly(Return(true));

  // Test failed validation with missing observation attribute
  EXPECT_CALL(mockObs_, hasAttribute("variable")).WillOnce(Return(false));
  EXPECT_FALSE(
      obsOperator_->validateCompatibility(*stateAdapter_, *obsAdapter_));
}

TEST_F(ObsOperatorTest, FluentInterfaceWorks) {
  // Setup expectations
  EXPECT_CALL(mockOperator_, initialize()).Times(1);
  EXPECT_CALL(mockOperator_, setParameter("param1", 1.0)).Times(1);
  EXPECT_CALL(mockOperator_, setParameter("param2", 2.0)).Times(1);

  // Test method chaining
  obsOperator_->initialize()
      .setParameter("param1", 1.0)
      .setParameter("param2", 2.0);
}

TEST_F(ObsOperatorTest, ThrowsWhenApplyingWithUninitializedOperator) {
  // Setup expectations
  EXPECT_CALL(mockOperator_, isInitialized()).WillOnce(Return(false));

  // Test exception
  EXPECT_THROW(obsOperator_->apply(*stateAdapter_, *obsAdapter_),
               std::runtime_error);
}

TEST_F(ObsOperatorTest, ThrowsWhenApplyingWithUninitializedState) {
  // Setup expectations
  EXPECT_CALL(mockState_, isInitialized()).WillOnce(Return(false));

  // Test exception
  EXPECT_THROW(obsOperator_->apply(*stateAdapter_, *obsAdapter_),
               std::runtime_error);
}

TEST_F(ObsOperatorTest, ThrowsWhenApplyingAdjointWithInvalidObs) {
  // Setup expectations
  EXPECT_CALL(mockObs_, isValid()).WillOnce(Return(false));

  // Test exception
  EXPECT_THROW(obsOperator_->applyAdjoint(*obsAdapter_, *incrementAdapter_),
               std::runtime_error);
}

}  // namespace metada::tests
