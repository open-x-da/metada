/**
 * @file MockObservationIterator.hpp
 * @brief Simple concrete observation iterator for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This class provides a simple concrete implementation of an observation
 * iterator for testing purposes. It implements the standard forward iterator
 * interface and can iterate over a vector of MockObservationPoint objects.
 *
 * The implementation supports:
 * - Standard forward iterator operations
 * - Iteration over observation point collections
 * - Copy and move semantics
 * - Comparison operations
 *
 * @see ObservationIterator
 * @see MockObservationPoint
 */

#pragma once

#include <iterator>
#include <vector>

#include "Location.hpp"

namespace metada::backends::gmock {

using metada::framework::Location;

/**
 * @brief Mock observation point for iterator testing
 */
struct MockObservationPoint {
  Location location;
  double value{0.0};
  double error{1.0};
  bool is_valid{true};

  MockObservationPoint() : location(0, 0, 0) {}
  MockObservationPoint(const Location& loc, double val, double err)
      : location(loc), value(val), error(err), is_valid(true) {}

  bool operator==(const MockObservationPoint& other) const {
    return location == other.location && value == other.value &&
           error == other.error && is_valid == other.is_valid;
  }
};

/**
 * @brief Simple concrete observation iterator for testing
 *
 * @details
 * Provides a concrete implementation of forward iterator operations for
 * iterating over MockObservationPoint collections:
 *
 * @par Iterator Operations
 * - operator*() - Dereference iterator
 * - operator->() - Arrow operator
 * - operator++() - Pre-increment
 * - operator++(int) - Post-increment
 * - operator==() - Equality comparison
 * - operator!=() - Inequality comparison
 *
 * @par Standard Iterator Traits
 * - iterator_category - Forward iterator tag
 * - value_type - MockObservationPoint
 * - difference_type - ptrdiff_t
 * - pointer - MockObservationPoint*
 * - reference - MockObservationPoint&
 *
 * @note This is a concrete iterator implementation that can be copied,
 * moved, and used with STL algorithms.
 */
class MockObservationIterator {
 public:
  // Standard iterator type traits
  using iterator_category = std::forward_iterator_tag;
  using value_type = MockObservationPoint;
  using difference_type = std::ptrdiff_t;
  using pointer = MockObservationPoint*;
  using reference = MockObservationPoint&;

  // Default constructor
  MockObservationIterator() = default;

  // Copy constructor
  MockObservationIterator(const MockObservationIterator&) = default;

  // Copy assignment
  MockObservationIterator& operator=(const MockObservationIterator&) = default;

  // Move constructor
  MockObservationIterator(MockObservationIterator&&) noexcept = default;

  // Move assignment
  MockObservationIterator& operator=(MockObservationIterator&&) noexcept =
      default;

  // Destructor
  ~MockObservationIterator() = default;

  // Simple concrete iterator implementation for testing
  MockObservationIterator(const std::vector<MockObservationPoint>* data,
                          size_t index)
      : data_(data), index_(index) {}

  // Iterator operations
  reference operator*() const {
    return const_cast<reference>((*data_)[index_]);
  }

  pointer operator->() const {
    return &const_cast<reference>((*data_)[index_]);
  }

  MockObservationIterator& operator++() {
    ++index_;
    return *this;
  }

  MockObservationIterator operator++(int) {
    MockObservationIterator tmp = *this;
    ++index_;
    return tmp;
  }

  bool operator==(const MockObservationIterator& other) const {
    return data_ == other.data_ && index_ == other.index_;
  }

  bool operator!=(const MockObservationIterator& other) const {
    return !(*this == other);
  }

 private:
  const std::vector<MockObservationPoint>* data_{nullptr};
  size_t index_{0};
};

}  // namespace metada::backends::gmock