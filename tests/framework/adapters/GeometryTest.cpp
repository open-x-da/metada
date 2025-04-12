/**
 * @file GeometryTest.cpp
 * @brief Unit tests for Geometry and GeometryIterator classes
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the Geometry class template and
 * GeometryIterator, which provide generic interfaces for geometry
 * implementations. The tests cover:
 *
 * Core functionality:
 * - Initialization and construction
 * - Grid iteration
 * - Periodicity queries
 * - Halo exchange
 * - GeometryIterator operations
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

// Define the backend tag for testing
using TestBackendTag = metada::traits::MockBackendTag;

/**
 * @brief Test fixture for GeometryIterator class tests
 *
 * Tests basic functionality of the iterator without requiring a complete
 * Geometry object
 */
class GeometryIteratorTest : public ::testing::Test {
 protected:
  // Mock iterator for testing
  std::unique_ptr<MockGeometryIterator<MockConfig>> iter_;
  // Test grid points for iteration sequence
  std::vector<MockGridPoint> points_;

  void SetUp() override {
    // Create mock iterator and test points
    iter_ = std::make_unique<MockGeometryIterator<MockConfig>>();
    points_ = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
  }

  void TearDown() override { iter_.reset(); }
};

/**
 * @brief Test basic iterator operations (dereference, increment, comparison)
 */
TEST_F(GeometryIteratorTest, BasicOperations) {
  // Set up test point
  MockGridPoint point{1, 2, 3};

  // Set up expectations
  EXPECT_CALL(*iter_, dereference()).WillOnce(Return(point));
  EXPECT_CALL(*iter_, compare(testing::_)).WillOnce(Return(true));
  EXPECT_CALL(*iter_, increment()).Times(2);

  // Test dereference
  MockGridPoint result = **iter_;
  EXPECT_EQ(result.x, point.x);
  EXPECT_EQ(result.y, point.y);
  EXPECT_EQ(result.z, point.z);

  // Test comparison
  MockGeometryIterator<MockConfig> otherIter;
  EXPECT_TRUE(*iter_ == otherIter);
  EXPECT_FALSE(*iter_ != otherIter);

  // Test increment operations
  ++(*iter_);  // Pre-increment
  (*iter_)++;  // Post-increment
}

/**
 * @brief Test iteration through a sequence of points
 */
TEST_F(GeometryIteratorTest, IterationSequence) {
  // Create a second iterator to act as the end sentinel
  auto endIter = std::make_unique<MockGeometryIterator<MockConfig>>();

  // Configure the iterator sequence
  iter_->mockIteratorSequence(points_);

  // Setup end comparison
  EXPECT_CALL(*iter_, compare(testing::_)).WillRepeatedly(Return(false));
  EXPECT_CALL(*endIter, compare(testing::_)).WillRepeatedly(Return(true));

  // Simulate iterating through sequence
  std::size_t pointIndex = 0;
  while (!(*iter_ == *endIter) && pointIndex < points_.size()) {
    MockGridPoint point = **iter_;
    EXPECT_EQ(point.x, points_[pointIndex].x);
    EXPECT_EQ(point.y, points_[pointIndex].y);
    EXPECT_EQ(point.z, points_[pointIndex].z);
    ++(*iter_);
    pointIndex++;
  }

  EXPECT_EQ(pointIndex, points_.size());
}

// Note: In a complete implementation, we would also have a GeometryTest fixture
// similar to StateTest, with tests for the full Geometry class functionalities:
// - Construction from Config
// - Periodicity queries
// - Size information
// - Halo exchange with State
// - Clone operation
// - Iterator begin()/end() access
//
// These tests would require the MockBackendTag to fully satisfy the
// GeometryBackendType concept including StateBackend support for haloExchange.

}  // namespace metada::tests

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}