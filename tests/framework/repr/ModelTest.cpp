#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MockIncrement.hpp"
#include "MockModel.hpp"
#include "MockState.hpp"
#include "Model.hpp"

namespace metada {
namespace framework {
namespace repr {
namespace tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

class ModelTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup common test data
    required_vars_ = {"temperature", "pressure"};
  }

  std::vector<std::string> required_vars_;
};

TEST_F(ModelTest, InitializeDelegatestoBackend) {
  MockModel mock_backend;
  Model<MockModel> model;

  EXPECT_CALL(model.backend(), initialize()).Times(1);

  model.initialize();
}

TEST_F(ModelTest, FinalizeDelegatestoBackend) {
  MockModel mock_backend;
  Model<MockModel> model;

  EXPECT_CALL(model.backend(), finalize()).Times(1);

  model.finalize();
}

TEST_F(ModelTest, StepOperations) {
  MockModel mock_backend;
  Model<MockModel> model;
  State<MockState> state1, state2;

  EXPECT_CALL(model.backend(), step(_, _)).Times(1);

  model.step(state1, state2);
}

TEST_F(ModelTest, TangentLinearOperations) {
  MockModel mock_backend;
  Model<MockModel> model;
  Increment<MockIncrement> inc1, inc2;

  EXPECT_CALL(model.backend(), stepTL(_, _)).Times(1);

  model.stepTL(inc1, inc2);
}

TEST_F(ModelTest, AdjointOperations) {
  MockModel mock_backend;
  Model<MockModel> model;
  Increment<MockIncrement> inc1, inc2;

  EXPECT_CALL(model.backend(), stepAD(_, _)).Times(1);

  model.stepAD(inc2, inc1);
}

TEST_F(ModelTest, ParameterOperations) {
  MockModel mock_backend;
  Model<MockModel> model;

  EXPECT_CALL(model.backend(), setParameter("dt", 0.1)).Times(1);
  EXPECT_CALL(model.backend(), getParameter("dt")).WillOnce(Return(0.1));

  model.setParameter("dt", 0.1);
  EXPECT_DOUBLE_EQ(model.getParameter("dt"), 0.1);
}

TEST_F(ModelTest, ModelInformation) {
  MockModel mock_backend;
  Model<MockModel> model;

  EXPECT_CALL(model.backend(), getRequiredStateVariables())
      .WillOnce(ReturnRef(required_vars_));
  EXPECT_CALL(model.backend(), isInitialized()).WillOnce(Return(true));

  EXPECT_EQ(model.getRequiredStateVariables(), required_vars_);
  EXPECT_TRUE(model.isInitialized());
}

TEST_F(ModelTest, BackendAccessors) {
  Model<MockModel> model;

  // Test that we can get non-const access to backend
  MockModel& backend = model.backend();
  EXPECT_CALL(backend, initialize()).Times(1);
  backend.initialize();

  // Test that we can get const access to backend
  const Model<MockModel>& const_model = model;
  const MockModel& const_backend = const_model.backend();
  EXPECT_CALL(const_backend, isInitialized()).WillOnce(Return(true));
  const_backend.isInitialized();
}

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada