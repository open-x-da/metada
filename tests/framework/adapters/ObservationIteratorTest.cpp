/**
 * @file ObservationIteratorTest.cpp
 * @brief Unit tests for ObservationIterator class
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This test suite verifies the functionality of the ObservationIterator class
 * template, which provides a generic interface for observation iterator
 * implementations. The tests cover:
 *
 * - Basic iterator operations (dereference, increment, comparison)
 * - Iteration through sequences of observation points
 * - Iterator semantics and behavior
 * - Copy and move operations
 * - Backend access methods
 * - STL compatibility and iterator traits
 *
 * The ObservationIterator is designed as a forward iterator that wraps
 * backend-specific iterator implementations while providing a consistent
 * interface across different backends. This test suite ensures that the adapter
 * correctly delegates operations to the backend implementation.
 *
 * The test suite uses a concrete MockObservationIterator for testing,
 * allowing verification of actual iteration behavior.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "MockBackendTraits.hpp"
#include "MockObservationIterator.hpp"
#include "ObservationIterator.hpp"

namespace metada::tests {

using namespace metada::framework;
using namespace metada::backends::gmock;

class ObservationIteratorTest : public ::testing::Test {
 protected:
  std::vector<MockObservationPoint> points_;

  void SetUp() override {
    // Create test observation points
    points_ = {MockObservationPoint(Location(45.0, -120.0, 1000.0), 25.5, 0.5),
               MockObservationPoint(Location(46.0, -121.0, 850.0), 15.2, 0.3),
               MockObservationPoint(Location(47.0, -122.0, 500.0), -5.8, 0.7)};
  }
};

/**
 * @brief Test basic iterator operations
 *
 * Verifies that the ObservationIterator correctly delegates basic operations
 * to the underlying backend iterator implementation.
 */
TEST_F(ObservationIteratorTest, BasicOperations) {
  // Create backend iterator pointing to first element
  MockObservationIterator backend_iter(&points_, 0);

  // Create the ObservationIterator adapter that wraps the backend
  ObservationIterator<traits::MockBackendTag> iter(backend_iter);

  // Test dereference - should return first point
  MockObservationPoint result = *iter;
  EXPECT_EQ(result.location, points_[0].location);
  EXPECT_EQ(result.value, points_[0].value);
  EXPECT_EQ(result.error, points_[0].error);

  // Test arrow operator
  auto* ptr = iter.operator->();
  EXPECT_EQ(ptr->value, points_[0].value);

  // Test pre-increment
  ++iter;
  MockObservationPoint second_result = *iter;
  EXPECT_EQ(second_result.value, points_[1].value);

  // Test post-increment
  auto iter_copy = iter++;
  MockObservationPoint third_result = *iter;
  EXPECT_EQ(third_result.value, points_[2].value);

  // The copy should still point to second element
  MockObservationPoint copy_result = *iter_copy;
  EXPECT_EQ(copy_result.value, points_[1].value);
}

/**
 * @brief Test default constructor
 *
 * Verifies that the default constructor creates a valid iterator.
 */
TEST_F(ObservationIteratorTest, DefaultConstructor) {
  ObservationIterator<traits::MockBackendTag> iter;

  // Default constructed iterator should be valid
  // We can't test much without actual data, but we verify it compiles
  EXPECT_TRUE(true);  // Basic existence test
}

/**
 * @brief Test copy constructor and assignment
 *
 * Verifies that the ObservationIterator can be copied correctly
 * and that copies behave independently.
 */
TEST_F(ObservationIteratorTest, CopyOperations) {
  MockObservationIterator backend_iter(&points_, 0);

  // Create original iterator
  ObservationIterator<traits::MockBackendTag> original(backend_iter);

  // Test copy constructor
  ObservationIterator<traits::MockBackendTag> copy_constructed(original);

  // Test copy assignment
  ObservationIterator<traits::MockBackendTag> copy_assigned;
  copy_assigned = original;

  // Verify all iterators point to same data initially
  EXPECT_EQ((*original).value, points_[0].value);
  EXPECT_EQ((*copy_constructed).value, points_[0].value);
  EXPECT_EQ((*copy_assigned).value, points_[0].value);

  // Advance original and verify copies are independent
  ++original;
  EXPECT_EQ((*original).value, points_[1].value);
  EXPECT_EQ((*copy_constructed).value,
            points_[0].value);                          // Should still be first
  EXPECT_EQ((*copy_assigned).value, points_[0].value);  // Should still be first
}

/**
 * @brief Test move constructor and assignment
 *
 * Verifies that the ObservationIterator can be moved correctly.
 */
TEST_F(ObservationIteratorTest, MoveOperations) {
  MockObservationIterator backend_iter(&points_, 0);

  // Create original iterator
  ObservationIterator<traits::MockBackendTag> original(backend_iter);

  // Test move constructor
  ObservationIterator<traits::MockBackendTag> move_constructed(
      std::move(original));

  // Test move assignment
  MockObservationIterator backend_iter2(&points_, 1);
  ObservationIterator<traits::MockBackendTag> temp(backend_iter2);
  ObservationIterator<traits::MockBackendTag> move_assigned;
  move_assigned = std::move(temp);

  // Verify moved-to objects work
  EXPECT_EQ((*move_constructed).value, points_[0].value);
  EXPECT_EQ((*move_assigned).value, points_[1].value);
}

/**
 * @brief Test backend access methods
 *
 * Verifies that the backend() methods provide correct access
 * to the underlying backend iterator.
 */
TEST_F(ObservationIteratorTest, BackendAccess) {
  MockObservationIterator backend_iter(&points_, 0);
  ObservationIterator<traits::MockBackendTag> iter(backend_iter);

  // Test non-const backend access
  auto& backend_ref = iter.backend();
  EXPECT_EQ((*backend_ref).value, points_[0].value);

  // Test const backend access
  const auto& const_iter = iter;
  const auto& const_backend_ref = const_iter.backend();
  EXPECT_EQ((*const_backend_ref).value, points_[0].value);
}

/**
 * @brief Test iterator type traits
 *
 * Verifies that the ObservationIterator has the correct iterator traits
 * for use with STL algorithms.
 */
TEST_F(ObservationIteratorTest, IteratorTraits) {
  using IteratorType = ObservationIterator<traits::MockBackendTag>;

  // Test iterator category
  static_assert(std::is_same_v<
                typename std::iterator_traits<IteratorType>::iterator_category,
                std::forward_iterator_tag>);

  // Test value type
  static_assert(
      std::is_same_v<typename std::iterator_traits<IteratorType>::value_type,
                     MockObservationPoint>);

  // Test that pointer and reference types are defined
  using pointer = typename std::iterator_traits<IteratorType>::pointer;
  using reference = typename std::iterator_traits<IteratorType>::reference;

  EXPECT_TRUE(true);  // Compilation test
}

/**
 * @brief Test iteration through observation sequence
 *
 * Verifies that the ObservationIterator can be used to iterate
 * through a sequence of observation points.
 */
TEST_F(ObservationIteratorTest, IterationSequence) {
  // Create iterators for begin and end
  MockObservationIterator begin_backend(&points_, 0);
  MockObservationIterator end_backend(&points_, points_.size());

  ObservationIterator<traits::MockBackendTag> begin_iter(begin_backend);
  ObservationIterator<traits::MockBackendTag> end_iter(end_backend);

  // Collect values through iteration
  std::vector<double> collected_values;
  auto current = begin_iter;

  for (size_t i = 0; i < points_.size(); ++i) {
    collected_values.push_back((*current).value);
    ++current;
  }

  // Verify all values were collected correctly
  EXPECT_EQ(collected_values.size(), points_.size());
  for (size_t i = 0; i < points_.size(); ++i) {
    EXPECT_EQ(collected_values[i], points_[i].value);
  }
}

/**
 * @brief Test comparison operators
 *
 * Verifies that equality and inequality operators work correctly.
 */
TEST_F(ObservationIteratorTest, ComparisonOperators) {
  MockObservationIterator backend_iter1(&points_, 0);
  MockObservationIterator backend_iter2(&points_, 0);
  MockObservationIterator backend_iter3(&points_, 1);

  ObservationIterator<traits::MockBackendTag> iter1(backend_iter1);
  ObservationIterator<traits::MockBackendTag> iter2(backend_iter2);
  ObservationIterator<traits::MockBackendTag> iter3(backend_iter3);

  // Test equality - same position
  EXPECT_TRUE(iter1 == iter2);
  EXPECT_FALSE(iter1 != iter2);

  // Test inequality - different positions
  EXPECT_FALSE(iter1 == iter3);
  EXPECT_TRUE(iter1 != iter3);
}

/**
 * @brief Test STL algorithm compatibility
 *
 * Verifies that ObservationIterator works with STL algorithms.
 */
TEST_F(ObservationIteratorTest, STLCompatibility) {
  // Create a vector of iterators pointing to different elements
  std::vector<ObservationIterator<traits::MockBackendTag>> iterators;

  for (size_t i = 0; i < points_.size(); ++i) {
    MockObservationIterator backend_iter(&points_, i);
    iterators.emplace_back(backend_iter);
  }

  // Use STL algorithm to find iterator pointing to specific value
  auto found =
      std::find_if(iterators.begin(), iterators.end(),
                   [](const auto& iter) { return (*iter).value == 15.2; });

  EXPECT_NE(found, iterators.end());
  EXPECT_EQ((**found).value, 15.2);
}

}  // namespace metada::tests