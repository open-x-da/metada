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
    ON_CALL(state1_->backend(), getDimensions())
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
 * @brief Test all available constructors
 *
 * @details
 * Verifies proper initialization and state transfer for:
 * - Configuration constructor
 * - Copy constructor with state validation
 * - Move constructor with proper state transfer
 */
TEST_F(StateTest, ConstructorTests) {
  // Test configuration constructor - use the helper method
  auto state = createState();
  EXPECT_TRUE(state.isInitialized());

  // For the copy constructor test, we need to use ON_CALL instead of
  // EXPECT_CALL because we can't set expectations on state_copy.backend()
  // before it exists
  ON_CALL(state.backend(), copyFrom(_)).WillByDefault(Return());

  // Test copy constructor
  State<Traits::StateType> state_copy(state);
  EXPECT_CALL(state_copy.backend(), equals(_)).WillOnce(Return(true));
  EXPECT_TRUE(state_copy == state);
  EXPECT_TRUE(
      state_copy.isInitialized());  // Copy should preserve initialization state

  // For the move constructor test, we need to use ON_CALL
  ON_CALL(state_copy.backend(), moveFrom(_)).WillByDefault(Return());

  // Test move constructor
  State<Traits::StateType> state_moved(std::move(state_copy));
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
  // Use the pre-created state object
  EXPECT_CALL(state1_->backend(), reset()).Times(1);
  state1_->reset();

  EXPECT_CALL(state1_->backend(), validate()).Times(1);
  state1_->validate();
}

/**
 * @brief Test copy assignment
 *
 * @details
 * Verifies proper state copying between instances:
 * - Correct delegation to backend copyFrom()
 * - Proper state transfer
 * - Preservation of initialization state
 */
TEST_F(StateTest, CopyAssignment) {
  // Use the pre-created state objects
  EXPECT_CALL(state2_->backend(), copyFrom(testing::Ref(state1_->backend())))
      .Times(1);
  *state2_ = *state1_;

  // Verify initialization state is preserved
  EXPECT_TRUE(state2_->isInitialized());
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
  // Create a new state for move testing
  auto temp_state = createState();
  EXPECT_TRUE(temp_state.isInitialized());

  // Set up expectations for move assignment
  EXPECT_CALL(state1_->backend(), moveFrom(_)).Times(1);

  // Test move assignment
  *state1_ = std::move(temp_state);

  // Verify states after move
  EXPECT_FALSE(temp_state.isInitialized());
  EXPECT_TRUE(state1_->isInitialized());
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
 * @brief Test that getData() returns correctly typed pointer
 *
 * @details
 * Verifies type-safe data access:
 * - Proper type casting
 * - Correct pointer return
 * - Backend delegation
 */
TEST_F(StateTest, GetDataReturnsTypedPointer) {
  // We've already set up the default behavior in SetUp()
  double* data = &state1_->getData<double>();
  EXPECT_EQ(data, test_data_.data());
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
  EXPECT_CALL(state1_->backend(), setMetadata("key1", "value1")).Times(1);
  EXPECT_CALL(state1_->backend(), getMetadata("key1"))
      .WillOnce(Return("value1"));

  state1_->setMetadata("key1", "value1");
  EXPECT_EQ(state1_->getMetadata("key1"), "value1");
}

/**
 * @brief Test state information accessors
 *
 * @details
 * Verifies state query operations:
 * - Variable name retrieval
 * - Variable existence check
 * - Dimension information access
 * - Backend delegation
 */
TEST_F(StateTest, StateInformation) {
  // We've already set up the default behavior in SetUp()
  EXPECT_EQ(state1_->getVariableNames(), variable_names_);
  EXPECT_EQ(state1_->getDimensions(), dimensions_);

  // Test hasVariable method
  EXPECT_TRUE(state1_->hasVariable("temperature"));
  EXPECT_TRUE(state1_->hasVariable("pressure"));
  EXPECT_FALSE(state1_->hasVariable("nonexistent_variable"));
}

/**
 * @brief Test arithmetic operators
 *
 * @details
 * Verifies state arithmetic operations:
 * - Addition operator (+)
 * - Subtraction operator (-)
 * - Addition assignment (+=)
 * - Subtraction assignment (-=)
 */
TEST_F(StateTest, ArithmeticOperations) {
  // Create a result state for testing
  auto result = createState();

  EXPECT_NO_THROW(result = *state1_ + *state2_);
  EXPECT_NO_THROW(auto result1 = *state1_ + *state2_);
  EXPECT_NO_THROW(result = *state1_ - *state2_);
  EXPECT_NO_THROW(auto result2 = *state1_ - *state2_);

  // Test addition assignment
  EXPECT_CALL(state1_->backend(), add(testing::Ref(state2_->backend())))
      .Times(1);
  *state1_ += *state2_;

  // Test subtraction assignment
  EXPECT_CALL(state1_->backend(), subtract(testing::Ref(state2_->backend())))
      .Times(1);
  *state1_ -= *state2_;

  // Test multiplication
  EXPECT_NO_THROW(result = *state1_ * 2.0);
  EXPECT_NO_THROW(result = 2.0 * *state1_);

  // Test multiplication assignment
  EXPECT_CALL(state1_->backend(), multiply(2.0)).Times(1);
  *state1_ *= 2.0;
}

/**
 * @brief Test arithmetic error conditions
 *
 * @details
 * Verifies proper error handling:
 * - Incompatible state dimensions
 * - Invalid state data
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
 * @brief Test MockState copy and move operations
 *
 * @details
 * Verifies that MockState properly supports:
 * - Copy construction
 * - Copy assignment
 * - Move construction
 * - Move assignment
 *
 * This test ensures that our mock implementation correctly supports
 * the copy/move semantics we've implemented in the State classes.
 */
TEST_F(StateTest, MockStateCopyAndMove) {
  // Get a reference to the backend of state1_
  auto& mockState1 = state1_->backend();

  // Test copy construction of MockState
  MockState mockState2(mockState1);

  // Test copy assignment of MockState
  MockState mockState3(getConfig());
  mockState3 = mockState1;

  // Test move construction of MockState
  MockState mockState4(std::move(mockState3));

  // Test move assignment of MockState
  MockState mockState5(getConfig());
  mockState5 = std::move(mockState4);

  // Verify that all operations completed without errors
  // Note: We don't need to verify specific behavior since the mock
  // implementation doesn't actually copy/move data, just maintains
  // the reference to the config
  SUCCEED() << "MockState copy and move operations completed successfully";
}

/**
 * @brief Test the zero() method that sets all state values to zero
 *
 * Verifies:
 * - The backend zero() method is called correctly
 * - The method returns a reference to itself for method chaining
 */
TEST_F(StateTest, ZeroMethod) {
  // Create a test state
  auto state = createState();

  // Expect the zero method to be called
  EXPECT_CALL(state.backend(), zero()).Times(1);

  // Call the zero method and verify it returns reference to self
  State<Traits::StateType>& result = state.zero();
  EXPECT_EQ(&result, &state);
}

}  // namespace metada::tests