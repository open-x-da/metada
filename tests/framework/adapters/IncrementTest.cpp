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
  };

  /**
   * @brief Get reference to the logger from context
   */
  Logger<Traits::LoggerType>& getLogger() { return context_->getLogger(); }

  /**
   * @brief Get reference to the config from context
   */
  Config<Traits::ConfigType>& getConfig() { return context_->getConfig(); }
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

}  // namespace metada::tests