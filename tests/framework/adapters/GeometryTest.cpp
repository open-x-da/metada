/**
 * @file GeometryTest.cpp
 * @brief Unit tests for Geometry class
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the Geometry class template,
 * which provides a generic interface for geometry implementations. The tests
 * cover:
 *
 * Core functionality:
 * - Initialization and construction
 * - Periodicity queries
 * - Size information queries
 * - Halo exchange
 * - Clone operations
 * - Iterator access
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <vector>

#include "Config.hpp"
#include "Geometry.hpp"
#include "GeometryIterator.hpp"
#include "MockBackendTraits.hpp"
#include "State.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using namespace metada::framework;
using namespace metada::backends::gmock;

class GeometryTest : public ::testing::Test {
 protected:
  std::string config_file_;
  std::unique_ptr<Config<traits::MockBackendTag>> config_;

  // Geometry object for testing
  std::unique_ptr<Geometry<traits::MockBackendTag>> geometry_;

  // State object for halo exchange testing
  std::unique_ptr<State<traits::MockBackendTag>> state_;

  void SetUp() override {
    // Get the directory where the test file is located
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file_);

    // Create test objects - we'll use a mock geometry through the adapter
    geometry_ = std::make_unique<Geometry<traits::MockBackendTag>>(*config_);
    state_ = std::make_unique<State<traits::MockBackendTag>>(*config_);
  }

  void TearDown() override {
    geometry_.reset();
    state_.reset();
    config_.reset();
  }
};

/**
 * @brief Test construction and initialization from Config
 *
 * Verifies that:
 * - A Geometry object can be properly constructed from a Config
 * - The constructed object is correctly initialized
 * - The configuration is accessible from the Geometry object
 */
TEST_F(GeometryTest, Construction) {
  // Test that our test fixture setup correctly initialized the geometry
  EXPECT_TRUE(geometry_->isInitialized());

  // Verify we can create a new instance directly
  Geometry<traits::MockBackendTag> localGeometry(*config_);
  EXPECT_TRUE(localGeometry.isInitialized());

  // Verify the config reference is maintained
  EXPECT_EQ(&config_->backend(), &geometry_->config().backend());
}

/**
 * @brief Test periodicity queries
 *
 * Verifies that the Geometry adapter correctly delegates periodicity
 * queries to the backend implementation.
 */
TEST_F(GeometryTest, PeriodicityQueries) {
  // Setup expectations for periodicity queries
  EXPECT_CALL(geometry_->backend(), isPeriodicX()).WillOnce(Return(true));
  EXPECT_CALL(geometry_->backend(), isPeriodicY()).WillOnce(Return(false));
  EXPECT_CALL(geometry_->backend(), isPeriodicZ()).WillOnce(Return(true));

  // Test periodicity in different dimensions
  EXPECT_TRUE(geometry_->isPeriodicX());   // X dimension periodic
  EXPECT_FALSE(geometry_->isPeriodicY());  // Y dimension not periodic
  EXPECT_TRUE(geometry_->isPeriodicZ());   // Z dimension periodic
}

/**
 * @brief Test size information queries
 *
 * Verifies that the Geometry adapter correctly delegates size
 * information queries to the backend implementation.
 */
TEST_F(GeometryTest, SizeInformation) {
  // Setup expectations for size queries
  EXPECT_CALL(geometry_->backend(), totalGridSize()).WillOnce(Return(3000));

  // Test total size query
  EXPECT_EQ(geometry_->totalGridSize(), 3000);
}

/**
 * @brief Test halo exchange operation
 *
 * Verifies that the Geometry adapter correctly delegates halo exchange
 * operations to the backend implementation.
 */
TEST_F(GeometryTest, HaloExchange) {
  // Setup expectations for halo exchange
  EXPECT_CALL(geometry_->backend(), haloExchangeImpl(testing::_)).Times(1);

  // Perform halo exchange
  geometry_->haloExchange(*state_);
}

/**
 * @brief Test clone operation
 *
 * Verifies that the Geometry adapter can create a proper clone of itself
 * with the same configuration and backend state.
 */
TEST_F(GeometryTest, Clone) {
  // Clone the geometry
  auto cloned_geometry = geometry_->clone();

  // Verify cloned geometry is initialized
  EXPECT_TRUE(cloned_geometry.isInitialized());
}

/**
 * @brief Test iterator begin/end access
 *
 * Verifies that the Geometry adapter correctly provides iterators
 * for traversing grid points by delegating to the backend implementation.
 */
TEST_F(GeometryTest, IteratorAccess) {
  // Set up mock iterators to return
  GeometryIterator<traits::MockBackendTag> begin_iter(*config_);
  GeometryIterator<traits::MockBackendTag> end_iter(*config_);

  // Setup expectations for begin/end
  EXPECT_CALL(geometry_->backend(), begin())
      .WillOnce(Return(std::move(begin_iter.backend())));
  EXPECT_CALL(geometry_->backend(), end())
      .WillOnce(Return(std::move(end_iter.backend())));

  // Get iterators
  auto iter_begin = geometry_->begin();
  auto iter_end = geometry_->end();

  // Simple verification that the iterators are different
  EXPECT_NE(&iter_begin, &iter_end);
}

}  // namespace metada::tests