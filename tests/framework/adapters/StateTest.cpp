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
 * - State operations (zero, dot, norm)
 * - Data access and type safety
 *
 * Advanced features:
 * - State information queries
 * - Move semantics
 * - Equality comparison
 * - Arithmetic operations
 * - Increment operations
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 *
 * @see State
 * @see StateBackendType
 * @see MockState
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>

#include "Config.hpp"
#include "MockBackendTraits.hpp"
#include "State.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::ByMove;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::Geometry;
using framework::State;

/**
 * @brief Test fixture for State class tests
 *
 * @details
 * Provides common test data and setup/teardown for State tests:
 * - Configuration file path
 * - Mock configuration object
 * - Test state objects
 * - Proper cleanup of allocated resources
 */
class StateTest : public ::testing::Test {
 protected:
  std::string config_file_;
  std::unique_ptr<Config<traits::MockBackendTag>> config_;
  std::unique_ptr<Geometry<traits::MockBackendTag>> geometry_;

  // Common State objects for tests
  std::unique_ptr<State<traits::MockBackendTag>> state1_;  ///< First test state
  std::unique_ptr<State<traits::MockBackendTag>>
      state2_;  ///< Second test state

  /**
   * @brief Set up test data before each test
   *
   * Initializes:
   * - Configuration file path
   * - Mock configuration object
   * - Geometry object
   * - Common State objects for testing
   */
  void SetUp() override {
    // Get the directory where the test file is located
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file_);
    geometry_ = std::make_unique<Geometry<traits::MockBackendTag>>(*config_);

    // Create common State objects
    state1_ =
        std::make_unique<State<traits::MockBackendTag>>(*config_, *geometry_);
    state2_ =
        std::make_unique<State<traits::MockBackendTag>>(*config_, *geometry_);
  }

  /**
   * @brief Clean up test data after each test
   */
  void TearDown() override {
    config_.reset();
    geometry_.reset();
    state1_.reset();
    state2_.reset();
  }

  /**
   * @brief Create a fresh State object for tests that need a new instance
   * @return A new State object initialized with the test config
   */
  State<traits::MockBackendTag> createState() {
    // Since Config is non-copyable, we pass a reference to the config
    return State<traits::MockBackendTag>(*config_, *geometry_);
  }
};

/**
 * @brief Test construction and initialization
 */
TEST_F(StateTest, Construction) {
  // Test configuration constructor - use the helper method
  auto state = createState();
  EXPECT_TRUE(state.isInitialized());
}

/**
 * @brief Test clone operation
 */
TEST_F(StateTest, Clone) {
  // Call the clone method
  auto state_clone = state1_->clone();

  // Verify the clone was successful
  EXPECT_TRUE(state_clone.isInitialized());
}

/**
 * @brief Test move operations
 */
TEST_F(StateTest, MoveOperations) {
  // First create two separate states to test with
  auto state_source = createState();
  EXPECT_TRUE(state_source.isInitialized());

  // Test move construction
  State<traits::MockBackendTag> state_moved(std::move(state_source));
  EXPECT_TRUE(state_moved.isInitialized());
  EXPECT_FALSE(
      state_source.isInitialized());  // Source should be in moved-from state

  // Create another state to test move assignment
  auto state_source2 = createState();
  EXPECT_TRUE(state_source2.isInitialized());

  // Create a target for move assignment
  auto state_target = createState();

  // Test move assignment
  state_target = std::move(state_source2);
  EXPECT_TRUE(state_target.isInitialized());
  EXPECT_FALSE(
      state_source2.isInitialized());  // Source should be in moved-from state
}

/**
 * @brief Test equality comparison
 */
TEST_F(StateTest, EqualityComparison) {
  EXPECT_CALL(state1_->backend(), equals(testing::Ref(state2_->backend())))
      .WillOnce(Return(true))
      .WillOnce(Return(false));
  EXPECT_TRUE(*state1_ == *state2_);
  EXPECT_FALSE(*state1_ == *state2_);

  EXPECT_CALL(state1_->backend(), equals(testing::Ref(state2_->backend())))
      .WillOnce(Return(true))
      .WillOnce(Return(false));
  EXPECT_FALSE(*state1_ != *state2_);
  EXPECT_TRUE(*state1_ != *state2_);
}

/**
 * @brief Test data access
 */
TEST_F(StateTest, DataAccess) {
  // Initialize test data in the state
  std::vector<double> expected_data = {1.0, 2.0, 3.0, 4.0, 5.0};
  state1_->backend().setData(expected_data);

  // Test non-const data access
  auto data_ptr = state1_->template getDataPtr<double>();
  EXPECT_NE(data_ptr, nullptr);  // Verify we got a valid pointer
  // Verify actual data values
  for (size_t i = 0; i < expected_data.size(); ++i) {
    EXPECT_DOUBLE_EQ(data_ptr[i], expected_data[i]);
  }

  // Test const data access
  const auto& const_state = *state1_;
  const auto const_data_ptr = const_state.template getDataPtr<double>();
  EXPECT_NE(const_data_ptr, nullptr);  // Verify we got a valid pointer
  // Verify actual data values
  for (size_t i = 0; i < expected_data.size(); ++i) {
    EXPECT_DOUBLE_EQ(const_data_ptr[i], expected_data[i]);
  }

  // Test data modification through non-const access
  data_ptr[0] = 10.0;
  EXPECT_DOUBLE_EQ(state1_->template getDataPtr<double>()[0], 10.0);
  EXPECT_DOUBLE_EQ(data_ptr[0],
                   10.0);  // Verify both access paths show the change

  // Test empty data case
  state1_->backend().setData({});
  EXPECT_EQ(state1_->size(), 0);
}

/**
 * @brief Test state information methods
 */
TEST_F(StateTest, StateInformation) {
  // Setup test variables and dimensions
  std::vector<std::string> expected_vars = {"temperature", "pressure",
                                            "humidity"};
  state1_->backend().setVariables(expected_vars);

  // Test getVariableNames
  const auto& vars = state1_->getVariableNames();
  EXPECT_EQ(vars, expected_vars);
  EXPECT_EQ(vars.size(), 3);
  EXPECT_EQ(vars[0], "temperature");
  EXPECT_EQ(vars[1], "pressure");
  EXPECT_EQ(vars[2], "humidity");

  // Test hasVariable
  EXPECT_TRUE(state1_->hasVariable("temperature"));
  EXPECT_TRUE(state1_->hasVariable("pressure"));
  EXPECT_TRUE(state1_->hasVariable("humidity"));
  EXPECT_FALSE(state1_->hasVariable("nonexistent_variable"));

  // Test size method
  EXPECT_CALL(state1_->backend(), size()).WillOnce(Return(6000));  // 10*20*30
  size_t state_size = state1_->size();
  EXPECT_EQ(state_size, 6000);
}

/**
 * @brief Test state operations (zero, dot, norm)
 */
TEST_F(StateTest, StateOperations) {
  // Test zero operation
  EXPECT_CALL(state1_->backend(), zero()).Times(1);
  State<traits::MockBackendTag>& result = state1_->zero();
  EXPECT_EQ(&result, state1_.get());

  // Test dot product
  EXPECT_CALL(state1_->backend(), dot(testing::Ref(state2_->backend())))
      .WillOnce(Return(42.0));
  double dot_product = state1_->dot(*state2_);
  EXPECT_DOUBLE_EQ(dot_product, 42.0);

  // Test norm calculation
  EXPECT_CALL(state1_->backend(), norm()).WillOnce(Return(10.0));
  double norm = state1_->norm();
  EXPECT_DOUBLE_EQ(norm, 10.0);
}

/**
 * @brief Test arithmetic operations
 */
TEST_F(StateTest, ArithmeticOperations) {
  // Create a result state for testing
  auto result = createState();

  // Test addition operator
  {
    EXPECT_NO_THROW(result = *state1_ + *state2_);
  }

  // Test subtraction operator
  {
    EXPECT_NO_THROW(result = *state1_ - *state2_);
  }

  // Test addition assignment (doesn't require cloning)
  EXPECT_CALL(state1_->backend(), add(testing::Ref(state2_->backend())))
      .Times(1);
  *state1_ += *state2_;

  // Test subtraction assignment (doesn't require cloning)
  EXPECT_CALL(state1_->backend(), subtract(testing::Ref(state2_->backend())))
      .Times(1);
  *state1_ -= *state2_;

  // Test multiplication by scalar
  {
    EXPECT_NO_THROW(result = *state1_ * 2.0);
  }

  // Test scalar * state multiplication (uses same codepath)
  {
    EXPECT_NO_THROW(result = 3.0 * *state1_);
  }

  // Test multiplication assignment (doesn't require cloning)
  EXPECT_CALL(state1_->backend(), multiply(2.0)).Times(1);
  *state1_ *= 2.0;
}

/**
 * @brief Test arithmetic error handling
 */
TEST_F(StateTest, ArithmeticErrors) {
  // Test addition with incompatible states
  EXPECT_CALL(state1_->backend(), add(testing::Ref(state2_->backend())))
      .WillOnce(testing::Throw(std::runtime_error("Incompatible states")));
  EXPECT_THROW(*state1_ += *state2_, std::runtime_error);

  // Test subtraction with incompatible states
  EXPECT_CALL(state1_->backend(), subtract(testing::Ref(state2_->backend())))
      .WillOnce(testing::Throw(std::runtime_error("Incompatible states")));
  EXPECT_THROW(*state1_ -= *state2_, std::runtime_error);
}

/**
 * @brief Test increment creation and application
 */
TEST_F(StateTest, IncrementOperations) {
  // We need a mock increment type for testing
  struct MockIncrement {
    // Keep track of which states were used to create the increment
    const State<traits::MockBackendTag>* state_a;
    const State<traits::MockBackendTag>* state_b;

    static MockIncrement createFromDifference(
        const State<traits::MockBackendTag>& a,
        const State<traits::MockBackendTag>& b) {
      return MockIncrement{&a, &b};
    }

    void applyTo(State<traits::MockBackendTag>& state) const {
      // Verify we got the expected states in createFromDifference
      EXPECT_NE(state_a, nullptr);
      EXPECT_NE(state_b, nullptr);

      // Call zero method on the state to verify the call
      state.zero();
    }
  };

  // Test creating an increment
  auto increment = state1_->createIncrementTo<MockIncrement>(*state2_);

  // Verify the increment was created with the correct states
  EXPECT_EQ(increment.state_a, state1_.get());
  EXPECT_EQ(increment.state_b, state2_.get());

  // Test applying the increment - should call zero() on state1_
  EXPECT_CALL(state1_->backend(), zero()).Times(1);
  state1_->applyIncrement(increment);
}

}  // namespace metada::tests