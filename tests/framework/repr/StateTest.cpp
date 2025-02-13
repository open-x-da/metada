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
 * @brief Test that getData() returns correctly typed pointer
 */
TEST_F(StateTest, GetDataReturnsTypedPointer) {
  MockState mock_backend;
  State<MockState> state(mock_backend);

  EXPECT_CALL(mock_backend, getData()).WillOnce(Return(test_data_));

  double* data = &state.getData<double>();
  EXPECT_EQ(data, test_data_);
}

/**
 * @brief Test metadata get/set operations
 */
TEST_F(StateTest, MetadataOperations) {
  MockState mock_backend;
  State<MockState> state(mock_backend);

  EXPECT_CALL(mock_backend, setMetadata("key1", "value1")).Times(1);
  EXPECT_CALL(mock_backend, getMetadata("key1")).WillOnce(Return("value1"));

  state.setMetadata("key1", "value1");
  EXPECT_EQ(state.getMetadata("key1"), "value1");
}

/**
 * @brief Test state information accessors
 */
TEST_F(StateTest, StateInformation) {
  MockState mock_backend;
  State<MockState> state(mock_backend);

  EXPECT_CALL(mock_backend, getVariableNames())
      .WillOnce(ReturnRef(variable_names_));
  EXPECT_CALL(mock_backend, getDimensions()).WillOnce(ReturnRef(dimensions_));

  EXPECT_EQ(state.getVariableNames(), variable_names_);
  EXPECT_EQ(state.getDimensions(), dimensions_);
}

/**
 * @brief Test copy operations
 */
TEST_F(StateTest, CopyOperations) {
  MockState mock_backend1, mock_backend2;
  State<MockState> state1(mock_backend1);
  State<MockState> state2(mock_backend2);

  // Test copy assignment
  EXPECT_CALL(mock_backend2, copyFrom(testing::Ref(mock_backend1))).Times(1);
  state2 = state1;

  // Test equality comparison
  EXPECT_CALL(mock_backend1, equals(testing::Ref(mock_backend2)))
      .WillOnce(Return(true))
      .WillOnce(Return(false));

  EXPECT_TRUE(state1 == state2);
  EXPECT_FALSE(state1 == state2);
}

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada