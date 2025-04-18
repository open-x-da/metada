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

#include "GeometryIterator.hpp"
#include "MockBackendTraits.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using namespace metada::framework;
using namespace metada::backends::gmock;

// Test fixture helper to create mock iterators - since MockGeometryIterator
// has a deleted default constructor and can't be copied, use shared_ptr
class MockGeometryIteratorFactory {
 public:
  // Factory method to create a mock iterator
  static std::shared_ptr<MockGeometryIterator> createMockIterator() {
    // Create a mock iterator with a null context
    return std::make_shared<MockGeometryIterator>(nullptr);
  }
};

/**
 * @brief Test fixture for GeometryIterator class tests
 *
 * Tests basic functionality of the iterator without requiring a complete
 * Geometry object. The fixture provides:
 * - Mock iterator instances
 * - Test grid points for iteration sequences
 */
class GeometryIteratorTest : public ::testing::Test {
 protected:
  // Mock iterator for testing - using a factory approach since default ctor is
  // deleted
  std::shared_ptr<MockGeometryIterator> mock_iter_;
  std::unique_ptr<GeometryIterator<traits::MockBackendTag>> iter_;
  // Test grid points for iteration sequence
  std::vector<MockGridPoint> points_;

  void SetUp() override {
    // Create mock iterator and test points using a factory approach
    mock_iter_ = MockGeometryIteratorFactory::createMockIterator();
    iter_ =
        std::make_unique<GeometryIterator<traits::MockBackendTag>>(*mock_iter_);
    points_ = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
  }

  void TearDown() override {
    iter_.reset();
    // No explicit cleanup needed for shared_ptr
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
  EXPECT_CALL(*mock_iter_, dereference()).WillOnce(Return(point));
  EXPECT_CALL(*mock_iter_, compare(testing::_))
      .WillOnce(Return(true))
      .WillOnce(Return(false));
  EXPECT_CALL(*mock_iter_, increment()).Times(2);

  // Test dereference
  MockGridPoint result = **iter_;
  EXPECT_EQ(result.x, point.x);
  EXPECT_EQ(result.y, point.y);
  EXPECT_EQ(result.z, point.z);

  // Test comparison
  std::shared_ptr<MockGeometryIterator> otherMockIter =
      MockGeometryIteratorFactory::createMockIterator();
  GeometryIterator<traits::MockBackendTag> otherIter(*otherMockIter);
  EXPECT_TRUE(*iter_ == otherIter);
  EXPECT_TRUE(*iter_ != otherIter);

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
  std::shared_ptr<MockGeometryIterator> endMockIter =
      MockGeometryIteratorFactory::createMockIterator();
  GeometryIterator<traits::MockBackendTag> endIter(*endMockIter);

  // Setup expectations - mock the sequence of points
  // First expect dereference calls to return our test points in sequence
  {
    ::testing::InSequence seq;
    for (const auto& point : points_) {
      EXPECT_CALL(*mock_iter_, dereference()).WillOnce(Return(point));
    }
  }

  // Configure compare behavior - return false until we reach the end
  size_t compareCallCount = 0;
  EXPECT_CALL(*mock_iter_, compare(testing::_))
      .WillRepeatedly([this, &compareCallCount](const auto&) {
        // Return false for first points_.size() calls, then true
        return (compareCallCount++ >= points_.size());
      });

  // Configure increment to be called for each point
  EXPECT_CALL(*mock_iter_, increment()).Times(points_.size());

  // Simulate iterating through sequence
  std::size_t pointIndex = 0;
  while (*iter_ != endIter && pointIndex < points_.size()) {
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