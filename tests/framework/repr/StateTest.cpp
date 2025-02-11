/**
 * @file StateTest.cpp
 * @brief Unit tests for the State class template
 */

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

/**
 * @brief Test fixture for State class tests
 */
class StateTest : public ::testing::Test {
 protected:
  /**
   * @brief Set up test data before each test
   */
  void SetUp() override {
    // Setup common test data
    variable_names_ = {"temperature", "pressure"};
    dimensions_ = {10, 20};
    test_data_ = new double[10];  // Example data
    for (int i = 0; i < 10; i++) test_data_[i] = i;
  }

  /**
   * @brief Clean up test data after each test
   */
  void TearDown() override { delete[] test_data_; }

  std::vector<std::string> variable_names_;  ///< Test variable names
  std::vector<size_t> dimensions_;           ///< Test dimensions
  double* test_data_;                        ///< Test data array
};

/**
 * @brief Test that initialize() delegates to backend
 */
TEST_F(StateTest, InitializeDelegatestoBackend) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), initialize()).Times(1);

  state.initialize();
}

/**
 * @brief Test that reset() delegates to backend
 */
TEST_F(StateTest, ResetDelegatestoBackend) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), reset()).Times(1);

  state.reset();
}

/**
 * @brief Test that validate() delegates to backend
 */
TEST_F(StateTest, ValidateDelegatestoBackend) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), validate()).Times(1);

  state.validate();
}

/**
 * @brief Test that getData() returns correctly typed pointer
 */
TEST_F(StateTest, GetDataReturnsTypedPointer) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), getData()).WillOnce(Return(test_data_));

  double* data = &state.getData<double>();
  EXPECT_EQ(data, test_data_);
}

/**
 * @brief Test metadata get/set operations
 */
TEST_F(StateTest, MetadataOperations) {
  MockState mock_backend;
  State<MockState> state;

  EXPECT_CALL(state.backend(), setMetadata("key1", "value1")).Times(1);
  EXPECT_CALL(state.backend(), getMetadata("key1")).WillOnce(Return("value1"));

  state.setMetadata("key1", "value1");
  EXPECT_EQ(state.getMetadata("key1"), "value1");
}

/**
 * @brief Test state information accessors
 */
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

/**
 * @brief Test backend accessor methods
 */
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