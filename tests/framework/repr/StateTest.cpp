#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MockState.hpp"
#include "State.hpp"

namespace metada {
namespace framework {
namespace repr {
namespace tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

class StateTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Setup common test data
    variable_names_ = {"temperature", "pressure"};
    dimensions_ = {10, 20};
    test_data_ = new double[10];  // Example data
    for (int i = 0; i < 10; i++) test_data_[i] = i;
  }

  void TearDown() override { delete[] test_data_; }

  std::vector<std::string> variable_names_;
  std::vector<size_t> dimensions_;
  double* test_data_;
};

TEST_F(StateTest, InitializeDelegatestoBackend) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), initialize()).Times(1);

  state.initialize();
}

TEST_F(StateTest, ResetDelegatestoBackend) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), reset()).Times(1);

  state.reset();
}

TEST_F(StateTest, ValidateDelegatestoBackend) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), validate()).Times(1);

  state.validate();
}

TEST_F(StateTest, GetDataReturnsTypedPointer) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), getData()).WillOnce(Return(test_data_));

  double* data = &state.getData<double>();
  EXPECT_EQ(data, test_data_);
}

TEST_F(StateTest, MetadataOperations) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), setMetadata("key1", "value1")).Times(1);
  EXPECT_CALL(state.backend(), getMetadata("key1")).WillOnce(Return("value1"));

  state.setMetadata("key1", "value1");
  EXPECT_EQ(state.getMetadata("key1"), "value1");
}

TEST_F(StateTest, StateInformation) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), getVariableNames())
      .WillOnce(ReturnRef(variable_names_));
  EXPECT_CALL(state.backend(), getDimensions())
      .WillOnce(ReturnRef(dimensions_));

  EXPECT_EQ(state.getVariableNames(), variable_names_);
  EXPECT_EQ(state.getDimensions(), dimensions_);
}

TEST_F(StateTest, BackendAccessors) {
  State<MockState> state;

  // Test that we can get non-const access to backend
  MockState& backend = state.backend();
  EXPECT_CALL(backend, initialize()).Times(1);
  backend.initialize();

  // Test that we can get const access to backend
  const State<MockState>& const_state = state;
  const MockState& const_backend = const_state.backend();
  EXPECT_CALL(const_backend, validate()).Times(1);
  const_backend.validate();
}

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada