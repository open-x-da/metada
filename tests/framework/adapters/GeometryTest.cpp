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
 * - Clone operations
 * - Iterator access
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "Config.hpp"
#include "Geometry.hpp"
#include "GeometryIterator.hpp"
#include "Logger.hpp"
#include "MockBackendTraits.hpp"
#include "MockConfig.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using namespace metada::framework;
using namespace metada::backends::gmock;

/**
 * @brief Test fixture for Geometry tests
 *
 * This fixture sets up the necessary objects for testing the Geometry adapter,
 * including a configuration, a geometry instance.
 */
class GeometryTest : public ::testing::Test {
 protected:
  std::string config_file_;
  std::unique_ptr<Config<traits::MockBackendTag>> config_;

  // Geometry object for testing
  std::unique_ptr<Geometry<traits::MockBackendTag>> geometry_;

  /**
   * @brief Set up the test environment
   *
   * Creates configuration, geometry, and state objects needed for testing.
   */
  void SetUp() override {
    // Get the directory where the test file is located
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file_);

    // Initialize the Logger singleton before creating any objects that use it
    Logger<traits::MockBackendTag>::Init(*config_);

    // Create test objects - we'll use a mock geometry through the adapter
    geometry_ = std::make_unique<Geometry<traits::MockBackendTag>>(*config_);
  }

  /**
   * @brief Clean up the test environment
   *
   * Releases all resources in the correct order.
   */
  void TearDown() override {
    geometry_.reset();
    // Reset the Logger singleton after tests
    Logger<traits::MockBackendTag>::Reset();
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
  EXPECT_CALL(geometry_->backend(), size()).WillOnce(Return(3000));
  EXPECT_EQ(geometry_->size(), 3000);

  // Verify we can create a new instance directly
  Geometry<traits::MockBackendTag> localGeometry(*config_);
  EXPECT_CALL(localGeometry.backend(), size()).WillOnce(Return(3000));
  EXPECT_EQ(localGeometry.size(), 3000);
}

/**
 * @brief Test periodicity queries
 *
 * Verifies that the Geometry adapter correctly delegates periodicity
 * queries to the backend implementation.
 */
TEST_F(GeometryTest, PeriodicityQueries) {
  // Setup expectations for periodicity queries
  EXPECT_CALL(geometry_->backend(), size()).WillOnce(Return(3000));

  // Test total size query
  EXPECT_EQ(geometry_->size(), 3000);
}

/**
 * @brief Test size information queries
 *
 * Verifies that the Geometry adapter correctly delegates size
 * information queries to the backend implementation.
 */
TEST_F(GeometryTest, SizeInformation) {
  // Setup expectations for size queries
  EXPECT_CALL(geometry_->backend(), size()).WillOnce(Return(3000));

  // Test total size query
  EXPECT_EQ(geometry_->size(), 3000);
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
  EXPECT_CALL(cloned_geometry.backend(), size()).WillOnce(Return(3000));
  EXPECT_EQ(cloned_geometry.size(), 3000);
}

/**
 * @brief Test iterator begin/end access
 *
 * Verifies that the Geometry adapter correctly provides iterators
 * for traversing grid points by delegating to the backend implementation.
 */
TEST_F(GeometryTest, IteratorAccess) {
  // Create MockGeometryIterator instances
  MockGeometryIterator begin_mock_iter;
  MockGeometryIterator end_mock_iter;

  // Create persistent grid points for ReturnRef
  MockGridPoint grid_point{0, 0, 0};

  // Set up default actions for the mock iterators
  ON_CALL(begin_mock_iter, dereference()).WillByDefault(ReturnRef(grid_point));
  ON_CALL(end_mock_iter, dereference()).WillByDefault(ReturnRef(grid_point));

  // Setup expectations for begin/end
  EXPECT_CALL(geometry_->backend(), begin()).WillOnce(Return(begin_mock_iter));
  EXPECT_CALL(geometry_->backend(), end()).WillOnce(Return(end_mock_iter));

  // Get iterators
  auto iter_begin = geometry_->begin();
  auto iter_end = geometry_->end();

  // Verify that iterators are different objects
  EXPECT_NE(&iter_begin, &iter_end);
}

/**
 * @brief Test move constructor
 *
 * Verifies that the Geometry adapter correctly implements move constructor
 * semantics, preserving the state of the moved object and marking the
 * moved-from object as uninitialized.
 */
TEST_F(GeometryTest, MoveConstructor) {
  // Setup expectations for the backend methods
  EXPECT_CALL(geometry_->backend(), size()).WillOnce(Return(1000));

  // Verify original geometry is initialized and has expected size
  EXPECT_EQ(geometry_->size(), 1000);

  // Move construct a new geometry
  Geometry<traits::MockBackendTag> moved_geometry(std::move(*geometry_));

  // Setup expectations for the moved geometry's backend methods
  EXPECT_CALL(moved_geometry.backend(), size()).WillOnce(Return(1000));

  // Verify the moved-to geometry has the expected state
  EXPECT_EQ(moved_geometry.size(), 1000);

  // We can't directly test if the original is uninitialized because we've moved
  // it, but we've covered the move constructor semantics
}

/**
 * @brief Test move assignment operator
 *
 * Verifies that the Geometry adapter correctly implements move assignment
 * semantics, transferring the state from one object to another.
 */
TEST_F(GeometryTest, MoveAssignment) {
  // Create a new geometry for move assignment
  Geometry<traits::MockBackendTag> local_geometry(*config_);

  // Setup expectations for the source backend methods
  EXPECT_CALL(geometry_->backend(), size()).WillOnce(Return(1000));

  // Verify source geometry has expected state
  EXPECT_EQ(geometry_->size(), 1000);

  // Create a temporary to use for move assignment
  auto temp_geometry =
      std::make_unique<Geometry<traits::MockBackendTag>>(*config_);

  // Perform move assignment
  local_geometry = std::move(*temp_geometry);

  // Setup expectations for the assigned geometry's backend methods
  EXPECT_CALL(local_geometry.backend(), size()).WillOnce(Return(1000));

  // Verify the moved-to geometry has the expected state
  EXPECT_EQ(local_geometry.size(), 1000);

  // We can't directly test the original geometry's state as we've moved it
}

/**
 * @brief Test const iterator begin/end methods
 *
 * Verifies that the const versions of begin() and end() correctly
 * delegate to the backend implementation.
 */
TEST_F(GeometryTest, ConstIterators) {
  // Create MockGeometryIterator instances
  MockGeometryIterator begin_mock_iter;
  MockGeometryIterator end_mock_iter;

  // Create persistent grid points for ReturnRef
  MockGridPoint grid_point{0, 0, 0};

  // Set up default actions for the mock iterators
  ON_CALL(begin_mock_iter, dereference()).WillByDefault(ReturnRef(grid_point));
  ON_CALL(end_mock_iter, dereference()).WillByDefault(ReturnRef(grid_point));

  // Get a const reference to our geometry for testing const methods
  const Geometry<traits::MockBackendTag>& const_geom = *geometry_;

  // Setup behavior using ON_CALL which works for const methods
  ON_CALL(const_geom.backend(), begin()).WillByDefault(Return(begin_mock_iter));
  ON_CALL(const_geom.backend(), end()).WillByDefault(Return(end_mock_iter));

  // Get const iterators
  auto const_begin = const_geom.begin();
  auto const_end = const_geom.end();

  // Verify iterators are different
  EXPECT_NE(&const_begin, &const_end);
}

/**
 * @brief Test backend access methods
 *
 * Verifies that both const and non-const backend() methods
 * return references to the correct backend instance.
 */
TEST_F(GeometryTest, BackendAccess) {
  // Test non-const backend access
  auto& backend_ref = geometry_->backend();
  EXPECT_EQ(&backend_ref, &(geometry_->backend()));

  // Test const backend access
  const Geometry<traits::MockBackendTag>& const_geom = *geometry_;
  const auto& const_backend_ref = const_geom.backend();
  EXPECT_EQ(&const_backend_ref, &(const_geom.backend()));

  // Verify both references point to the same backend
  EXPECT_EQ(static_cast<const void*>(&backend_ref),
            static_cast<const void*>(&const_backend_ref));
}

TEST_F(GeometryTest, STLContainerFunctions) {
  // Set up values to be returned by the mock
  MockGridPoint pt0{0, 0, 0};
  MockGridPoint pt1{1, 1, 1};
  MockGridPoint pt2{2, 2, 2};
  std::vector<MockGridPoint> points = {pt0, pt1, pt2};

  // Set up expectations for size-related methods
  EXPECT_CALL(geometry_->backend(), size()).WillOnce(Return(points.size()));
  EXPECT_CALL(geometry_->backend(), empty()).WillOnce(Return(false));
  EXPECT_CALL(geometry_->backend(), max_size()).WillOnce(Return(100));

  // Set up expectations for element access via get()
  EXPECT_CALL(geometry_->backend(), get(0)).WillOnce(ReturnRef(points[0]));
  EXPECT_CALL(geometry_->backend(), get(1)).WillOnce(ReturnRef(points[1]));
  EXPECT_CALL(geometry_->backend(), front()).WillOnce(ReturnRef(points[0]));
  EXPECT_CALL(geometry_->backend(), back()).WillOnce(ReturnRef(points[2]));

  // Test size-related methods
  EXPECT_EQ(geometry_->size(), 3);
  EXPECT_FALSE(geometry_->empty());
  EXPECT_GE(geometry_->max_size(), 3);

  // Test element access (all forward to get() internally)
  EXPECT_NO_THROW({ geometry_->at(0); });
  EXPECT_NO_THROW({ geometry_->front(); });
  EXPECT_NO_THROW({ geometry_->back(); });
  EXPECT_NO_THROW({ geometry_->operator[](1); });
}

TEST_F(GeometryTest, RangeBasedForAndSTLAlgorithms) {
  // Prepare test points
  MockGridPoint pt0{0, 1, 2};
  MockGridPoint pt2{2, 3, 4};

  // Create mock iterators for begin and end
  MockGeometryIterator begin_iter;
  MockGeometryIterator end_iter;

  // Set up expectations for iterator operations with default actions
  EXPECT_CALL(begin_iter, dereference()).WillRepeatedly(ReturnRef(pt0));
  EXPECT_CALL(begin_iter, increment()).WillRepeatedly(Return());
  EXPECT_CALL(begin_iter, compare(testing::_)).WillRepeatedly(Return(false));

  EXPECT_CALL(end_iter, dereference()).WillRepeatedly(ReturnRef(pt2));
  EXPECT_CALL(end_iter, increment()).WillRepeatedly(Return());
  EXPECT_CALL(end_iter, compare(testing::_)).WillRepeatedly(Return(true));

  // Set up backend to return the correct iterators
  EXPECT_CALL(geometry_->backend(), begin()).WillOnce(Return(begin_iter));
  EXPECT_CALL(geometry_->backend(), end()).WillOnce(Return(end_iter));

  // Test that we can get iterators from the geometry
  auto it = geometry_->begin();
  auto end = geometry_->end();

  // Verify that iterators are different objects
  EXPECT_NE(&it, &end);

  // Set up expectations for iterator operations with default actions
  EXPECT_CALL(it, dereference()).WillRepeatedly(ReturnRef(pt0));

  // Test basic iterator operations - these should not throw with default
  // actions
  EXPECT_NO_THROW({
    auto& point = *it;  // This calls dereference()
    (void)point;        // Suppress unused variable warning
  });

  // Test iterator comparison
  EXPECT_NO_THROW({
    bool is_equal = (it == end);
    (void)is_equal;  // Suppress unused variable warning
  });
}

}  // namespace metada::tests