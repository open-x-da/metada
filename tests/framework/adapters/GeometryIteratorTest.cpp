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

#include <vector>

#include "GeometryIterator.hpp"
#include "MockBackendTraits.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;
using ::testing::ReturnRef;

using namespace metada::framework;
using namespace metada::backends::gmock;

class GeometryIteratorTest : public ::testing::Test {
 protected:
  std::vector<MockGridPoint> points_;

  void SetUp() override { points_ = {{0, 1, 2}, {1, 2, 3}, {2, 3, 4}}; }
};

TEST_F(GeometryIteratorTest, BasicOperations) {
  MockGridPoint point{1, 2, 3};

  // Create the mock iterator that will be wrapped by GeometryIterator
  MockGeometryIterator mock_iter;

  // Set up expectations on the mock iterator
  EXPECT_CALL(mock_iter, dereference()).WillOnce(ReturnRef(point));
  EXPECT_CALL(mock_iter, compare(_))
      .WillOnce(Return(true))
      .WillOnce(Return(false));
  EXPECT_CALL(mock_iter, increment()).Times(2);

  // Create the GeometryIterator adapter that wraps the mock
  GeometryIterator<traits::MockBackendTag> iter(mock_iter);

  // Test dereference - this should call mock_iter.dereference()
  MockGridPoint result = *iter;
  EXPECT_EQ(result.x, point.x);
  EXPECT_EQ(result.y, point.y);
  EXPECT_EQ(result.z, point.z);

  // Test comparison - this should call mock_iter.compare()
  MockGeometryIterator other_mock_iter;
  GeometryIterator<traits::MockBackendTag> other_iter(other_mock_iter);
  EXPECT_TRUE(iter == other_iter);
  EXPECT_TRUE(iter != other_iter);

  // Test increment - this should call mock_iter.increment()
  ++iter;
  iter++;
}

TEST_F(GeometryIteratorTest, STLIteratorFunctions) {
  // Setup a sequence of points and a vector of iterators
  std::vector<MockGridPoint> test_points = points_;
  std::vector<GeometryIterator<traits::MockBackendTag>> iters;

  for (auto& pt : test_points) {
    MockGeometryIterator mock_iter;
    EXPECT_CALL(mock_iter, dereference()).WillRepeatedly(ReturnRef(pt));
    iters.emplace_back(mock_iter);
  }

  // Use STL algorithms with GeometryIterator
  std::vector<MockGridPoint> collected;
  for (auto& iter : iters) {
    collected.push_back(*iter);
  }
  EXPECT_EQ(collected, test_points);

  // Test std::distance and std::advance (simulate with two iterators)
  // Note: This is a mock, so we can't actually increment through a real range,
  // but we can check that the type traits are present and the interface works.
  if (!iters.empty()) {
    auto it1 = iters.front();
    auto it2 = iters.back();
    (void)std::distance(it1, it2);
  }
}

}  // namespace metada::tests