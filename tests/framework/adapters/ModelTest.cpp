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
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;

class ModelTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(ModelTest, InitializeDelegatestoBackend) {
  Model<MockModel, MockState> model;
  EXPECT_CALL(model.backend(), initialize()).Times(1);
  model.initialize();
}

TEST_F(ModelTest, FinalizeDelegatestoBackend) {
  Model<MockModel, MockState> model;
  EXPECT_CALL(model.backend(), finalize()).Times(1);
  model.finalize();
}

TEST_F(ModelTest, StepDelegatestoBackend) {
  Model<MockModel, MockState> model;
  EXPECT_CALL(model.backend(), step(0.1)).Times(1);
  model.step(0.1);
}

TEST_F(ModelTest, StateManagement) {
  Model<MockModel, MockState> model;

  // Create a State wrapper around MockState
  MockState mockStateBackend;
  State<MockState> mockState(mockStateBackend);

  // Set expectations for setState
  EXPECT_CALL(model.backend(), setState(testing::Ref(mockState.backend())))
      .Times(1);
  model.setState(mockState);

  // Set expectations for getState
  MockState backendState;
  EXPECT_CALL(model.backend(), getState()).WillOnce(ReturnRef(backendState));

  State<MockState> retrievedState = model.getState();
  EXPECT_EQ(&retrievedState.backend(), &backendState);
}

TEST_F(ModelTest, ParameterOperations) {
  Model<MockModel, MockState> model;

  EXPECT_CALL(model.backend(), setParameter("dt", 0.1)).Times(1);
  EXPECT_CALL(model.backend(), getParameter("dt")).WillOnce(Return(0.1));

  model.setParameter("dt", 0.1);
  EXPECT_DOUBLE_EQ(model.getParameter("dt"), 0.1);
}

TEST_F(ModelTest, ModelMetadata) {
  Model<MockModel, MockState> model;

  EXPECT_CALL(model.backend(), getName()).WillOnce(Return("TestModel"));
  EXPECT_CALL(model.backend(), getVersion()).WillOnce(Return("1.0.0"));
  EXPECT_CALL(model.backend(), getCurrentTime()).WillOnce(Return(10.5));
  EXPECT_CALL(model.backend(), isInitialized()).WillOnce(Return(true));

  EXPECT_EQ(model.getName(), "TestModel");
  EXPECT_EQ(model.getVersion(), "1.0.0");
  EXPECT_DOUBLE_EQ(model.getCurrentTime(), 10.5);
  EXPECT_TRUE(model.isInitialized());
}

TEST_F(ModelTest, BackendAccessors) {
  Model<MockModel, MockState> model;

  // Test that we can get non-const access to backend
  MockModel& backend = model.backend();
  EXPECT_CALL(backend, initialize()).Times(1);
  backend.initialize();

  // Test that we can get const access to backend
  const Model<MockModel, MockState>& const_model = model;
  const MockModel& const_backend = const_model.backend();
  EXPECT_CALL(const_backend, isInitialized()).WillOnce(Return(true));
  const_backend.isInitialized();
}

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada