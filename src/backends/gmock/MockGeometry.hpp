/**
 * @file MockGeometry.hpp
 * @brief Mock implementation of geometry backend for testing
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This file provides a Google Mock implementation of a geometry backend
 * for use in unit tests. It implements all interfaces required by the
 * GeometryBackendType concept.
 *
 * The mock implementation allows tests to:
 * - Set expectations on method calls
 * - Verify interactions with the geometry backend
 * - Return controlled test values
 * - Simulate backend behavior without requiring a real implementation
 *
 * @see GeometryBackendType
 * @see MockGeometryIterator
 */

#pragma once

#include <gmock/gmock.h>

#include <cstddef>
#include <vector>

#include "MockGeometryIterator.hpp"
#include "MockGridPoint.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of a geometry backend for testing
 *
 * @details
 * This class provides a Google Mock implementation of a geometry backend
 * for use in unit tests. It implements all methods required by the
 * GeometryBackendType concept and allows tests to set expectations on
 * method calls.
 *
 * The mock geometry supports:
 * - Iteration through grid points via begin/end methods
 * - Size information queries
 * - Halo exchange operations
 * - Proper move semantics (non-copyable)
 * - Cloning capability
 *
 * @see MockGeometryIterator
 * @see GeometryBackendType
 */
class MockGeometry {
 public:
  /** @brief Type of elements stored in the geometry */
  using value_type = MockGridPoint;
  /** @brief Reference to an element */
  using reference = value_type&;
  /** @brief Const reference to an element */
  using const_reference = const value_type&;
  /** @brief Pointer to an element */
  using pointer = value_type*;
  /** @brief Const pointer to an element */
  using const_pointer = const value_type*;
  /** @brief Type used for size and indexing */
  using size_type = std::size_t;
  /** @brief Type used for iterator differences */
  using difference_type = std::ptrdiff_t;
  /** @brief Iterator type for traversing grid points */
  using iterator = MockGeometryIterator;
  /** @brief Const iterator type for traversing grid points */
  using const_iterator = MockGeometryIterator;

  /**
   * @brief Constructor that takes a mock config
   * @tparam ConfigBackend The mock configuration backend type
   * @param config Mock configuration object
   * @details Initializes the geometry with a set of sample grid points
   */
  template <typename ConfigBackend>
  explicit MockGeometry(const ConfigBackend& /*config*/)
      : gridPoints_({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}}) {}

  /** @brief Default constructor is deleted to ensure proper initialization */
  MockGeometry() = delete;

  /** @brief Copy constructor is deleted (non-copyable) */
  MockGeometry(const MockGeometry&) = delete;

  /** @brief Copy assignment operator is deleted (non-copyable) */
  MockGeometry& operator=(const MockGeometry&) = delete;

  /**
   * @brief Move constructor
   * @param other Geometry to move from
   * @details Moves the grid points from the source geometry
   */
  MockGeometry(MockGeometry&& other) noexcept
      : gridPoints_(std::move(other.gridPoints_)) {}

  /**
   * @brief Move assignment operator
   * @param other Geometry to move from
   * @return Reference to this geometry after assignment
   * @details Moves the grid points from the source geometry if not
   * self-assignment
   */
  MockGeometry& operator=(MockGeometry&& other) noexcept {
    if (this != &other) {
      gridPoints_ = std::move(other.gridPoints_);
    }
    return *this;
  }

  /**
   * @brief Create a clone of this geometry
   * @return A new MockGeometry instance with the same grid points
   * @details Creates a deep copy of the geometry with identical grid points
   */
  MockGeometry clone() const { return MockGeometry(this->gridPoints_); }

  /** @brief Mock method for getting iterator to first grid point */
  MOCK_METHOD(iterator, begin, ());

  /** @brief Mock method for getting iterator past the last grid point */
  MOCK_METHOD(iterator, end, ());

  /** @brief Mock method for getting const iterator to first grid point */
  MOCK_METHOD(const_iterator, begin, (), (const));

  /** @brief Mock method for getting const iterator past the last grid point */
  MOCK_METHOD(const_iterator, end, (), (const));

  /** @brief Mock method for getting const iterator to first grid point */
  MOCK_METHOD(const_iterator, cbegin, (), (const));
  /** @brief Mock method for getting const iterator past the last grid point */
  MOCK_METHOD(const_iterator, cend, (), (const));

  /** @brief Mock method for getting the number of grid points */
  MOCK_METHOD(size_type, size, (), (const));
  /** @brief Mock method for checking if geometry is empty */
  MOCK_METHOD(bool, empty, (), (const));
  /** @brief Mock method for getting maximum possible size */
  MOCK_METHOD(size_type, max_size, (), (const));

  /** @brief Mock method for accessing grid point by index */
  MOCK_METHOD(reference, get, (size_type idx), ());
  /** @brief Mock method for accessing grid point by index (const version) */
  MOCK_METHOD(const_reference, get, (size_type idx), (const));

  /** @brief Access grid point by index */
  reference operator[](size_type idx) { return get(idx); }
  /** @brief Access grid point by index (const version) */
  const_reference operator[](size_type idx) const { return get(idx); }
  /** @brief Access grid point by index with bounds checking */
  reference at(size_type idx) { return get(idx); }
  /** @brief Access grid point by index with bounds checking (const version) */
  const_reference at(size_type idx) const { return get(idx); }
  /** @brief Mock method for accessing first grid point */
  MOCK_METHOD(reference, front, (), ());
  /** @brief Mock method for accessing first grid point (const version) */
  MOCK_METHOD(const_reference, front, (), (const));
  /** @brief Mock method for accessing last grid point */
  MOCK_METHOD(reference, back, (), ());
  /** @brief Mock method for accessing last grid point (const version) */
  MOCK_METHOD(const_reference, back, (), (const));

 private:
  /** @brief Sample grid points for testing */
  std::vector<MockGridPoint> gridPoints_;

  /**
   * @brief Private constructor for cloning
   * @param points Vector of grid points to initialize with
   * @details Used internally by clone() to create a copy with the same points
   */
  MockGeometry(const std::vector<MockGridPoint>& points)
      : gridPoints_(points) {}
};

}  // namespace metada::backends::gmock