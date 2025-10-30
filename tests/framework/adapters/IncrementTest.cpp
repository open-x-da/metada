#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "Geometry.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "MockBackendTraits.hpp"
#include "MockConfig.hpp"
#include "MockLogger.hpp"
#include "MockState.hpp"
#include "State.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using Config = framework::Config<traits::MockBackendTag>;
using Logger = framework::Logger<traits::MockBackendTag>;
using State = framework::State<traits::MockBackendTag>;
using Increment = metada::framework::Increment<traits::MockBackendTag>;

class IncrementTest : public ::testing::Test {
 protected:
  std::vector<size_t> dimensions_;
  std::unique_ptr<metada::framework::Config<traits::MockBackendTag>> config_;
  std::unique_ptr<metada::framework::Geometry<traits::MockBackendTag>>
      geometry_;
  std::unique_ptr<metada::framework::State<traits::MockBackendTag>> entity1_;
  std::unique_ptr<metada::framework::State<traits::MockBackendTag>> entity2_;

  void SetUp() override {
    dimensions_ = {10, 20};
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    auto config_file = (test_dir / "test_config.yaml").string();
    config_ =
        std::make_unique<metada::framework::Config<traits::MockBackendTag>>(
            config_file);
    metada::framework::Logger<traits::MockBackendTag>::Init(*config_);
    geometry_ =
        std::make_unique<metada::framework::Geometry<traits::MockBackendTag>>(
            *config_);
    entity1_ =
        std::make_unique<metada::framework::State<traits::MockBackendTag>>(
            *config_, *geometry_);
    entity2_ =
        std::make_unique<metada::framework::State<traits::MockBackendTag>>(
            *config_, *geometry_);
  }

  void TearDown() override {
    dimensions_.clear();
    entity1_.reset();
    entity2_.reset();
    geometry_.reset();
    config_.reset();
    metada::framework::Logger<traits::MockBackendTag>::Reset();
  }
};

TEST_F(IncrementTest, RandomizeOperation) {
  // Create increment from geometry
  Increment increment =
      Increment::createFromGeometry(entity1_->geometry()->backend());
  increment.randomize();

  // Get data as vector
  auto data = increment.getData<std::vector<double>>();
  for (size_t i = 0; i < data.size(); ++i) {
    EXPECT_GE(data[i], -0.5);
    EXPECT_LE(data[i], 0.5);
  }
}

TEST_F(IncrementTest, DotProductImplementation) {
  auto& state1 = *entity1_;
  auto& state2 = *entity2_;
  state1.backend().setData({1.0, 2.0, 3.0, 4.0, 5.0});
  state2.backend().setData({5.0, 4.0, 3.0, 2.0, 1.0});

  // Create increments from geometry and transfer state data
  Increment inc1 = Increment::createFromGeometry(state1.geometry()->backend());
  Increment inc2 = Increment::createFromGeometry(state2.geometry()->backend());
  inc1.incrementBackend().transferFromState(state1.backend());
  inc2.incrementBackend().transferFromState(state2.backend());

  double expected = 1.0 * 5.0 + 2.0 * 4.0 + 3.0 * 3.0 + 4.0 * 2.0 + 5.0 * 1.0;
  EXPECT_DOUBLE_EQ(inc1.dot(inc2), expected);
}

}  // namespace metada::tests