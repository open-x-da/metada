#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "AppTraits.hpp"
#include "ApplicationContext.hpp"
#include "Increment.hpp"
#include "MockConfig.hpp"
#include "MockLogger.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using framework::Config;
using framework::Increment;
using framework::Logger;
using framework::State;
using framework::runs::ApplicationContext;

using Traits = AppTraits<MockLogger, MockConfig, MockState>;

class IncrementTest : public ::testing::Test {
 protected:
  std::unique_ptr<ApplicationContext<Traits>> context_;
  // Test data
  std::vector<size_t> dimensions_;

  // Mock entity for testing
  std::unique_ptr<State<Traits::StateType>> entity1_;
  std::unique_ptr<State<Traits::StateType>> entity2_;

  void SetUp() override {
    context_ = std::make_unique<ApplicationContext<Traits>>("IncrementTest");
    dimensions_ = {10, 20};

    entity1_ = std::make_unique<State<Traits::StateType>>(getConfig());
    entity2_ = std::make_unique<State<Traits::StateType>>(getConfig());
  }

  void TearDown() override {
    context_.reset();
    dimensions_.clear();
    entity1_.reset();
    entity2_.reset();
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
   * @brief Create a test increment from the difference between entity1_ and
   * entity2_
   */
  Increment<State<Traits::StateType>> createTestIncrement() {
    return Increment<State<Traits::StateType>>(*entity1_, *entity2_);
  }
};

/**
 * @brief Test creating an increment from two entities
 *
 * @details
 * Verifies that an increment can be properly created as the difference
 * between two state entities.
 */
TEST_F(IncrementTest, CreateFromTwoEntities) {
  // Create increment from the difference of two entities
  Increment<State<Traits::StateType>> increment(*entity2_, *entity1_);

  // Verify the increment was created successfully
  EXPECT_TRUE(increment.isInitialized());
}

/**
 * @brief Test creating an increment using the factory method
 */
// TEST_F(IncrementTest, CreateFromDifference) {
//  Create increment using the factory method
// auto increment = Increment<State<Traits::StateType>>::createFromDifference(
//     *entity2_, *entity1_);

// Verify the increment was created successfully
// EXPECT_TRUE(increment.isInitialized());
//}

/**
 * @brief Test zero operation
 */
// TEST_F(IncrementTest, ZeroOperation) {
//   auto increment = createTestIncrement();

// Expect the zero operation to be called on the backend
// EXPECT_CALL(entity1_->backend(), zero()).Times(1);

//   increment.zero();
// }

/**
 * @brief Test scale operation
 */
TEST_F(IncrementTest, ScaleOperation) {
  auto increment = createTestIncrement();

  // Expect the scale operation to be called on the backend
  // EXPECT_CALL(entity1_->backend(), operator*=(2.5)).Times(1);

  increment.scale(2.5);
}

/**
 * @brief Test axpy operation
 */
TEST_F(IncrementTest, AxpyOperation) {
  auto increment1 = createTestIncrement();
  auto increment2 = createTestIncrement();
  // Expect the axpy operation to be implemented correctly
  // EXPECT_CALL(entity2_->backend(), operator*=(3.0)).Times(1);
  // EXPECT_CALL(entity1_->backend(), operator+=(entity2_->backend())).Times(1);

  increment1.axpy(3.0, increment2);
}

/**
 * @brief Test dot product operation
 */
// TEST_F(IncrementTest, DotProductOperation) {
//   auto increment1 = createTestIncrement();
//   auto increment2 = createTestIncrement();

// Expect the dot product to be calculated on the backend
// EXPECT_CALL(entity1_->backend(), dot(entity2_->backend()))
//     .WillOnce(Return(42.0));

//   EXPECT_DOUBLE_EQ(increment1.dot(increment2), 42.0);
// }

/**
 * @brief Test norm operation
 */
// TEST_F(IncrementTest, NormOperation) {
//   auto increment = createTestIncrement();

// Expect the norm to be calculated on the backend
// EXPECT_CALL(entity1_->backend(), norm()).WillOnce(Return(5.0));

//   EXPECT_DOUBLE_EQ(increment.norm(), 5.0);
// }

/**
 * @brief Test applying an increment to an entity
 */
TEST_F(IncrementTest, ApplyToOperation) {
  auto increment = createTestIncrement();

  // Expect the increment to be applied to entity2
  // EXPECT_CALL(entity2_->backend(), operator+=(entity1_->backend())).Times(1);

  increment.applyTo(*entity2_);
}

/**
 * @brief Test operator overloads
 */
TEST_F(IncrementTest, OperatorOverloads) {
  auto increment1 = createTestIncrement();
  auto increment2 = createTestIncrement();

  // Test addition operator
  // EXPECT_CALL(entity1_->backend(), operator+=(entity2_->backend())).Times(1);
  auto result1 = increment1 + increment2;

  // Test multiplication operator
  // EXPECT_CALL(entity1_->backend(), operator*=(3.5)).Times(1);
  auto result2 = increment1 * 3.5;

  // Test non-member multiplication operator
  // EXPECT_CALL(entity1_->backend(), operator*=(4.0)).Times(1);
  auto result3 = 4.0 * increment1;
}

}  // namespace metada::tests