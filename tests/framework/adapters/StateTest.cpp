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
 * - Copy/move semantics
 * - Equality comparison
 * - Arithmetic operations
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 *
 * @see State
 * @see IState
 * @see MockState
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "MockConfig.hpp"
#include "MockLogger.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::ByMove;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::Logger;
using framework::State;
using framework::runs::ApplicationContext;

using Traits = AppTraits<MockLogger, MockConfig, MockState>;

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
  std::unique_ptr<ApplicationContext<Traits>>
      context_;  ///< Application context for testing
  std::vector<std::string> variable_names_;  ///< Test variable names
  std::vector<size_t> dimensions_;           ///< Test dimensions
  std::vector<double> test_data_;            ///< Test data array

  // Common State objects for tests
  std::unique_ptr<State<Traits::StateType>> state1_;  ///< First test state
  std::unique_ptr<State<Traits::StateType>> state2_;  ///< Second test state

  /**
   * @brief Set up test data before each test
   *
   * Initializes:
   * - Application context
   * - Variable names array
   * - Dimensions array
   * - Test data array with sequential values
   * - Common State objects
   */
  void SetUp() override {
    // Create a new application context for testing
    context_ = std::make_unique<ApplicationContext<Traits>>("StateTest");

    // Setup common test data
    variable_names_ = {"temperature", "pressure"};
    dimensions_ = {10, 20};
    test_data_.resize(10);  // Example data
    for (int i = 0; i < 10; i++) test_data_[i] = i;

    // Create common State objects
    state1_ = std::make_unique<State<Traits::StateType>>(getConfig());
    state2_ = std::make_unique<State<Traits::StateType>>(getConfig());

    // Setup common expectations for test data access
    ON_CALL(state1_->backend(), getData())
        .WillByDefault(Return(test_data_.data()));
    ON_CALL(state1_->backend(), getVariableNames())
        .WillByDefault(ReturnRef(variable_names_));
    ON_CALL(state1_->backend(), getDimensions(testing::Ref(variable_names_[0])))
        .WillByDefault(ReturnRef(dimensions_));
  }

  /**
   * @brief Clean up test data after each test
   */
  void TearDown() override {
    // Clean up State objects first
    state1_.reset();
    state2_.reset();

    // Then clean up other resources
    test_data_.clear();
    dimensions_.clear();
    variable_names_.clear();
    context_.reset();
  }

  /**
   * @brief Get reference to the logger from context
   */
  Logger<Traits::LoggerType>& getLogger() { return context_->getLogger(); }

  /**
   * @brief Get reference to the config from context
   */
  Config<Traits::ConfigType>& getConfig() { return context_->getConfig(); }

  /**
   * @brief Create a fresh State object for tests that need a new instance
   * @return A new State object initialized with the test config
   */
  State<Traits::StateType> createState() {
    // Since Config is non-copyable, we pass a reference to the config
    return State<Traits::StateType>(getConfig());
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
  // Set expectation on the backend
  EXPECT_CALL(state1_->backend(), clone()).Times(1);

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
  State<Traits::StateType> state_moved(std::move(state_source));
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
  // We've already set up the default behavior in SetUp()
  double* data = &state1_->getData<double>();
  EXPECT_EQ(data, test_data_.data());
}

/**
 * @brief Test state information methods
 */
TEST_F(StateTest, StateInformation) {
  // Setup expectations for getVariableNames
  EXPECT_CALL(state1_->backend(), getVariableNames())
      .WillRepeatedly(ReturnRef(variable_names_));

  // Test getVariableNames method
  const auto& vars = state1_->getVariableNames();
  EXPECT_EQ(vars, variable_names_);
  EXPECT_EQ(vars.size(), 2);
  EXPECT_EQ(vars[0], "temperature");
  EXPECT_EQ(vars[1], "pressure");

  // Setup expectations for getDimensions
  EXPECT_CALL(state1_->backend(), getDimensions("temperature"))
      .WillRepeatedly(ReturnRef(dimensions_));

  // Test getDimensions method
  const auto& dims = state1_->getDimensions("temperature");
  EXPECT_EQ(dims, dimensions_);
  EXPECT_EQ(dims.size(), 2);
  EXPECT_EQ(dims[0], 10);
  EXPECT_EQ(dims[1], 20);

  // Test hasVariable method
  EXPECT_TRUE(state1_->hasVariable("temperature"));
  EXPECT_TRUE(state1_->hasVariable("pressure"));
  EXPECT_FALSE(state1_->hasVariable("nonexistent_variable"));
}

/**
 * @brief Test state operations (zero, dot, norm)
 */
TEST_F(StateTest, StateOperations) {
  // Test zero operation
  EXPECT_CALL(state1_->backend(), zero()).Times(1);
  State<Traits::StateType>& result = state1_->zero();
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

  // Setup for binary arithmetic operations that require cloning
  {
    EXPECT_CALL(state1_->backend(), clone()).Times(1);

    // The add operation is called on the clone, not the original
    EXPECT_CALL(result.backend(), add(testing::Ref(state2_->backend())))
        .Times(1);

    result = *state1_ + *state2_;
  }

  // Test subtraction operator
  {
    EXPECT_CALL(state1_->backend(), clone()).Times(1);

    EXPECT_CALL(result.backend(), subtract(testing::Ref(state2_->backend())))
        .Times(1);

    result = *state1_ - *state2_;
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
    EXPECT_CALL(state1_->backend(), clone()).Times(1);
    EXPECT_CALL(result.backend(), multiply(2.0)).Times(1);

    result = *state1_ * 2.0;
  }

  // Test scalar * state multiplication (uses same codepath)
  {
    EXPECT_CALL(state1_->backend(), clone()).Times(1);
    EXPECT_CALL(result.backend(), multiply(3.0)).Times(1);

    result = 3.0 * *state1_;
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
    const State<Traits::StateType>* state_a;
    const State<Traits::StateType>* state_b;

    static MockIncrement createFromDifference(
        const State<Traits::StateType>& a, const State<Traits::StateType>& b) {
      return MockIncrement{&a, &b};
    }

    void applyTo(State<Traits::StateType>& state) const {
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