/**
 * @file StateTest.cpp
 * @brief Unit tests for the State class template
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the State class template,
 * which provides a generic interface for state implementations. The tests
 * cover:
 *
 * Core functionality:
 * - Initialization and construction
 * - State operations (reset, validate)
 * - Data access and type safety
 *
 * Advanced features:
 * - Metadata management
 * - State information queries
 * - Copy/move semantics
 * - Equality comparison
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 *
 * @see State
 * @see IState
 * @see MockState
 */

#ifndef METADA_TESTS_FRAMEWORK_REPR_STATE_TEST_HPP_
#define METADA_TESTS_FRAMEWORK_REPR_STATE_TEST_HPP_

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <optional>

#include "Config.hpp"
#include "MockConfig.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada {
namespace framework {
namespace repr {
namespace tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using metada::framework::tools::config::Config;
using metada::framework::tools::config::tests::MockConfig;

/**
 * @brief Test fixture for State class tests
 *
 * @details
 * Provides common test data and setup/teardown for State tests:
 * - Variable names and dimensions
 * - Test data array
 * - Proper cleanup of allocated resources
 */
class StateTest : public ::testing::Test {
 protected:
  Config<MockConfig> config_;
  /**
   * @brief Set up test data before each test
   *
   * Initializes:
   * - Variable names array
   * - Dimensions array
   * - Test data array with sequential values
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
   *
   * Frees allocated test data array
   */
  void TearDown() override { delete[] test_data_; }

  std::vector<std::string> variable_names_;  ///< Test variable names
  std::vector<size_t> dimensions_;           ///< Test dimensions
  double* test_data_;                        ///< Test data array
};

/**
 * @brief Test all available constructors
 *
 * @details
 * Verifies proper initialization and state transfer for:
 * - Configuration constructor
 * - Copy constructor with state validation
 * - Move constructor with proper state transfer
 */
TEST_F(StateTest, ConstructorTests) {
  // Test configuration constructor
  State<MockState> state(config_);
  EXPECT_TRUE(state.isInitialized());

  // Test copy constructor
  State<MockState> state_copy(state);
  EXPECT_CALL(state_copy.backend(), equals(_)).WillOnce(Return(true));
  EXPECT_TRUE(state_copy == state);

  // Test move constructor
  State<MockState> state_moved(std::move(state_copy));
  // Verify the moved-from state is no longer initialized
  EXPECT_FALSE(state_copy.isInitialized());
  EXPECT_TRUE(state_moved.isInitialized());
}

/**
 * @brief Test core state operations
 *
 * @details
 * Verifies basic state manipulation:
 * - reset() properly resets state
 * - validate() checks state consistency
 */
TEST_F(StateTest, CoreStateOperations) {
  State<MockState> state(config_);

  // Test reset
  EXPECT_CALL(state.backend(), reset()).Times(1);
  state.reset();

  // Test validate
  EXPECT_CALL(state.backend(), validate()).Times(1);
  state.validate();
}

/**
 * @brief Test copy assignment
 *
 * @details
 * Verifies proper state copying between instances:
 * - Correct delegation to backend copyFrom()
 * - Proper state transfer
 */
TEST_F(StateTest, CopyAssignment) {
  State<MockState> state1(config_);
  State<MockState> state2(config_);

  // Test copy assignment
  EXPECT_CALL(state2.backend(), copyFrom(testing::Ref(state1.backend())))
      .Times(1);
  state2 = state1;
}

/**
 * @brief Test move assignment
 *
 * @details
 * Verifies proper move semantics:
 * - Correct state transfer
 * - Source state invalidation
 * - Destination state validation
 */
TEST_F(StateTest, MoveAssignment) {
  State<MockState> state1(config_);
  State<MockState> state2(config_);

  // Set up expectations for move assignment
  EXPECT_CALL(state2.backend(), moveFrom(_)).Times(1);

  // Test move assignment
  state2 = std::move(state1);

  // Verify states after move
  EXPECT_FALSE(state1.isInitialized());
  EXPECT_TRUE(state2.isInitialized());
}

/**
 * @brief Test equality and inequality comparison
 *
 * @details
 * Verifies comparison operators:
 * - Equality (==) with both true/false cases
 * - Inequality (!=) with both true/false cases
 * - Proper delegation to backend equals()
 */
TEST_F(StateTest, EqualityComparison) {
  State<MockState> state1(config_);
  State<MockState> state2(config_);

  EXPECT_CALL(state1.backend(), equals(testing::Ref(state2.backend())))
      .WillOnce(Return(true))
      .WillOnce(Return(false));
  EXPECT_TRUE(state1 == state2);
  EXPECT_FALSE(state1 == state2);

  EXPECT_CALL(state1.backend(), equals(testing::Ref(state2.backend())))
      .WillOnce(Return(true))
      .WillOnce(Return(false));
  EXPECT_FALSE(state1 != state2);
  EXPECT_TRUE(state1 != state2);
}

/**
 * @brief Test that getData() returns correctly typed pointer
 *
 * @details
 * Verifies type-safe data access:
 * - Proper type casting
 * - Correct pointer return
 * - Backend delegation
 */
TEST_F(StateTest, GetDataReturnsTypedPointer) {
  State<MockState> state(config_);

  EXPECT_CALL(state.backend(), getData()).WillOnce(Return(test_data_));

  double* data = &state.getData<double>();
  EXPECT_EQ(data, test_data_);
}

/**
 * @brief Test metadata get/set operations
 *
 * @details
 * Verifies metadata management:
 * - Setting key-value pairs
 * - Retrieving values by key
 * - Backend delegation for storage
 */
TEST_F(StateTest, MetadataOperations) {
  State<MockState> state(config_);

  EXPECT_CALL(state.backend(), setMetadata("key1", "value1")).Times(1);
  EXPECT_CALL(state.backend(), getMetadata("key1")).WillOnce(Return("value1"));

  state.setMetadata("key1", "value1");
  EXPECT_EQ(state.getMetadata("key1"), "value1");
}

/**
 * @brief Test state information accessors
 *
 * @details
 * Verifies state query operations:
 * - Variable name retrieval
 * - Dimension information access
 * - Backend delegation
 */
TEST_F(StateTest, StateInformation) {
  State<MockState> state(config_);

  EXPECT_CALL(state.backend(), getVariableNames())
      .WillOnce(ReturnRef(variable_names_));
  EXPECT_CALL(state.backend(), getDimensions())
      .WillOnce(ReturnRef(dimensions_));

  EXPECT_EQ(state.getVariableNames(), variable_names_);
  EXPECT_EQ(state.getDimensions(), dimensions_);
}

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_REPR_STATE_TEST_HPP_