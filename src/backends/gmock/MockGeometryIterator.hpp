/**
 * @file MockGeometryIterator.hpp
 * @brief Mock implementation of geometry iterator backend for testing
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This file provides a Google Mock implementation of a geometry iterator
 * backend for use in unit tests. It implements all interfaces required by the
 * GeometryIteratorBackendType concept.
 *
 * The mock implementation allows tests to:
 * - Set expectations on iterator operations
 * - Verify interactions with the iterator
 * - Return controlled test values
 * - Simulate iteration behavior without requiring a real implementation
 *
 * @see GeometryIterator
 * @see GeometryIteratorBackendType
 * @see MockGeometry
 */

#pragma once

#include <gmock/gmock.h>

#include <iterator>
#include <memory>
#include <vector>

#include "MockGridPoint.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of a geometry iterator backend for testing
 *
 * @details
 * This class provides a Google Mock implementation of a geometry iterator
 * backend for use in unit tests. It implements all methods required by the
 * GeometryIteratorBackendType concept and allows tests to set expectations on
 * iterator behavior.
 *
 * The mock iterator supports:
 * - Dereference operations to access grid points
 * - Increment operations to advance the iterator
 * - Comparison operations to check iterator positions
 * - Proper copy and move semantics
 *
 * This mock is used in tests to verify that the GeometryIterator adapter
 * correctly delegates operations to the backend implementation.
 *
 * @see MockGeometry
 * @see GeometryIterator
 * @see GeometryIteratorBackendType
 */
class MockGeometryIterator {
 public:
  // STL-style type aliases
  using iterator_category = std::forward_iterator_tag;
  using value_type = MockGridPoint;
  using difference_type = std::ptrdiff_t;
  using pointer = MockGridPoint*;
  using reference = MockGridPoint&;

  /** @brief Default constructor (for end/sentinels) */
  MockGeometryIterator() = default;

  /**
   * @brief Constructor that takes a mock geometry or other context
   * @param context Pointer to context object (can be null for testing)
   */
  explicit MockGeometryIterator([[maybe_unused]] void* context) {}

  /**
   * @brief Copy constructor
   *
   * @details
   * Google Mock's MOCK_METHOD members are not copyable by default, which would
   * normally make this class non-copyable. However, STL iterators and the
   * framework concepts require that iterators be copyable. To satisfy this,
   * we provide a copy constructor with an empty body. This allows the iterator
   * to be copied as required by STL algorithms and the framework, but does NOT
   * copy the internal Google Mock state (which is fine for test iterators).
   *
   * This is a common workaround for using Google Mock with STL-style iterators.
   */
  MockGeometryIterator(const MockGeometryIterator&) {}

  /**
   * @brief Copy assignment operator
   *
   * @details
   * Same rationale as the copy constructor: this is required for STL and
   * framework compatibility, but does not copy any Google Mock state.
   */
  MockGeometryIterator& operator=(const MockGeometryIterator&) { return *this; }

  /**
   * @brief Move constructor
   * @param other Iterator to move from
   */
  MockGeometryIterator(MockGeometryIterator&& other) noexcept
      : MockGeometryIterator(other) {}  // Delegate to copy constructor

  /**
   * @brief Move assignment
   * @param other Iterator to move from
   * @return Reference to this iterator after assignment
   */
  MockGeometryIterator& operator=(MockGeometryIterator&& other) noexcept {
    if (this != &other) {
      *this = other;  // Delegate to copy assignment
    }
    return *this;
  }

  /**
   * @brief Destructor
   */
  ~MockGeometryIterator() = default;

  /** @brief Mock method for dereferencing the iterator */
  MOCK_METHOD(reference, dereference, (), (const));

  /** @brief Mock method for incrementing the iterator */
  MOCK_METHOD(void, increment, ());

  /** @brief Mock method for comparing iterators */
  MOCK_METHOD(bool, compare, (const MockGeometryIterator&), (const));

  /**
   * @brief Dereference operator to access the current grid point
   * @return The current grid point
   */
  reference operator*() const { return dereference(); }

  /**
   * @brief Arrow operator (returns pointer to a static value for testability)
   * @return Pointer to the current grid point
   */
  pointer operator->() const { return &dereference(); }

  /**
   * @brief Pre-increment operator to advance to the next grid point
   * @return Reference to this iterator after advancement
   */
  MockGeometryIterator& operator++() {
    increment();
    return *this;
  }

  /**
   * @brief Post-increment operator to advance to the next grid point
   * @return Copy of this iterator after advancement
   */
  MockGeometryIterator operator++(int) {
    MockGeometryIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  /**
   * @brief Equality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators point to the same position
   */
  bool operator==(const MockGeometryIterator& other) const {
    return compare(other);
  }

  /**
   * @brief Inequality comparison operator
   * @param other Iterator to compare with
   * @return True if iterators point to different positions
   */
  bool operator!=(const MockGeometryIterator& other) const {
    return !compare(other);
  }
};

}  // namespace metada::backends::gmock