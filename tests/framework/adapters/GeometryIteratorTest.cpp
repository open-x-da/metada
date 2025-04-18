/**
 * @file GeometryIteratorTest.cpp
 * @brief Unit tests for GeometryIterator class
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the GeometryIterator class
 * template, which provides a generic interface for geometry iterator
 * implementations. The tests cover:
 *
 * - Basic iterator operations (dereference, increment, comparison)
 * - Iteration through sequences of grid points
 * - Iterator semantics and behavior
 *
 * The GeometryIterator is designed as a forward iterator that wraps
 * backend-specific iterator implementations while providing a consistent
 * interface across different backends. This test suite ensures that the adapter
 * correctly delegates operations to the backend implementation.
 *
 * The test suite uses Google Test/Mock framework for mocking and assertions.
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <vector>

#include "Config.hpp"
#include "GeometryIterator.hpp"
#include "MockBackendTraits.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using namespace metada::framework;
using namespace metada::backends::gmock;

/**
 * @brief Test fixture for GeometryIterator class tests
 *
 * Tests basic functionality of the iterator without requiring a complete
 * Geometry object. The fixture provides:
 * - A mock configuration
 * - A mock iterator instance
 * - Test grid points for iteration sequences
 */
class GeometryIteratorTest : public ::testing::Test {
 protected:
  std::string config_file_;
  std::unique_ptr<Config<traits::MockBackendTag>> config_;
  // Mock iterator for testing
  std::unique_ptr<GeometryIterator<traits::MockBackendTag>> iter_;
  // Test grid points for iteration sequence
  std::vector<MockGridPoint> points_;

  void SetUp() override {
    // Create and configure the mock model
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file_);

    // Create mock iterator and test points
    iter_ =
        std::make_unique<GeometryIterator<traits::MockBackendTag>>(*config_);
    points_ = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
  }

  void TearDown() override {
    iter_.reset();
    config_.reset();
  }
};

/**
 * @brief Test basic iterator operations (dereference, increment, comparison)
 *
 * Verifies that the GeometryIterator correctly delegates basic operations
 * to the backend implementation:
 * - Dereference operator (*) returns the current grid point
 * - Equality/inequality operators compare iterator positions
 * - Pre-increment and post-increment operators advance the iterator
 */
TEST_F(GeometryIteratorTest, BasicOperations) {
  // Set up test point
  MockGridPoint point{1, 2, 3};

  // Set up expectations on the backend iterator
  EXPECT_CALL(iter_->backend(), dereference()).WillOnce(Return(point));
  EXPECT_CALL(iter_->backend(), compare(testing::_))
      .WillOnce(Return(true))
      .WillOnce(Return(true));
  EXPECT_CALL(iter_->backend(), increment()).Times(2);

  // Test dereference
  MockGridPoint result = **iter_;
  EXPECT_EQ(result.x, point.x);
  EXPECT_EQ(result.y, point.y);
  EXPECT_EQ(result.z, point.z);

  // Test comparison
  GeometryIterator<traits::MockBackendTag> otherIter(*config_);
  EXPECT_TRUE(*iter_ == otherIter);
  EXPECT_FALSE(*iter_ != otherIter);

  // Test increment operations
  ++(*iter_);  // Pre-increment
  (*iter_)++;  // Post-increment
}

/**
 * @brief Test iteration through a sequence of points
 *
 * Verifies that the GeometryIterator can be used to traverse a sequence
 * of grid points, simulating how it would be used in a range-based for loop
 * or standard algorithm. This test ensures:
 * - Proper iteration semantics (begin to end)
 * - Correct dereferencing of points during iteration
 * - Proper termination when reaching the end
 */
TEST_F(GeometryIteratorTest, IterationSequence) {
  // Create a second iterator to act as the end sentinel
  auto endIter =
      std::make_unique<GeometryIterator<traits::MockBackendTag>>(*config_);

  // Setup expectations - mock the sequence of points
  // First expect dereference calls to return our test points in sequence
  {
    ::testing::InSequence seq;
    for (const auto& point : points_) {
      EXPECT_CALL(iter_->backend(), dereference()).WillOnce(Return(point));
    }
  }

  // Configure compare behavior - return false until we reach the end
  size_t compareCallCount = 0;
  EXPECT_CALL(iter_->backend(), compare(testing::_))
      .WillRepeatedly([this, &compareCallCount](const auto&) {
        // Return false for first points_.size() calls, then true
        return (compareCallCount++ >= points_.size());
      });

  // Configure increment to be called for each point
  EXPECT_CALL(iter_->backend(), increment()).Times(points_.size());

  // Simulate iterating through sequence
  std::size_t pointIndex = 0;
  while (*iter_ != *endIter && pointIndex < points_.size()) {
    MockGridPoint point = **iter_;
    EXPECT_EQ(point.x, points_[pointIndex].x);
    EXPECT_EQ(point.y, points_[pointIndex].y);
    EXPECT_EQ(point.z, points_[pointIndex].z);
    ++(*iter_);
    pointIndex++;
  }

  EXPECT_EQ(pointIndex, points_.size());
}

}  // namespace metada::tests