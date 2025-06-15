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
#include <memory>
#include <stdexcept>
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
 * - Periodicity queries in X, Y, and Z dimensions
 * - Size information queries
 * - Halo exchange operations
 * - Initialization status checks
 * - Proper move semantics (non-copyable)
 * - Cloning capability
 *
 * @tparam ConfigBackend The mock configuration backend type
 *
 * @see MockGeometryIterator
 * @see GeometryBackendType
 */
template <typename ConfigBackend>
class MockGeometry {
 public:
  // STL-style type aliases
  using value_type = MockGridPoint;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using iterator = MockGeometryIterator;
  using const_iterator = MockGeometryIterator;

  /**
   * @brief Constructor that takes a mock config
   * @param config Mock configuration object
   */
  explicit MockGeometry(const ConfigBackend& config)
      : config_(config), gridPoints_({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}}) {}

  /** @brief Default constructor is deleted to ensure proper initialization */
  MockGeometry() = delete;

  /** @brief Copy constructor is deleted (non-copyable) */
  MockGeometry(const MockGeometry&) = delete;

  /** @brief Copy assignment operator is deleted (non-copyable) */
  MockGeometry& operator=(const MockGeometry&) = delete;

  /**
   * @brief Move constructor
   * @param other Geometry to move from
   */
  MockGeometry(MockGeometry&& other) noexcept
      : config_(other.config_), gridPoints_(std::move(other.gridPoints_)) {}

  /**
   * @brief Move assignment operator
   * @param other Geometry to move from
   * @return Reference to this geometry after assignment
   */
  MockGeometry& operator=(MockGeometry&& other) noexcept {
    if (this != &other) {
      gridPoints_ = std::move(other.gridPoints_);
    }
    return *this;
  }

  /**
   * @brief Create a clone of this geometry
   * @return A new MockGeometry instance with the same configuration
   */
  MockGeometry clone() const { return MockGeometry(config_); }

  /** @brief Mock method for getting iterator to first grid point */
  MOCK_METHOD(iterator, begin, ());

  /** @brief Mock method for getting iterator past the last grid point */
  MOCK_METHOD(iterator, end, ());

  /** @brief Mock method for getting const iterator to first grid point */
  MOCK_METHOD(const_iterator, begin, (), (const));

  /** @brief Mock method for getting const iterator past the last grid point */
  MOCK_METHOD(const_iterator, end, (), (const));

  // STL-style cbegin/cend (mocked)
  MOCK_METHOD(const_iterator, cbegin, (), (const));
  MOCK_METHOD(const_iterator, cend, (), (const));

  // Size information
  MOCK_METHOD(size_type, size, (), (const));
  MOCK_METHOD(bool, empty, (), (const));
  MOCK_METHOD(size_type, max_size, (), (const));

  // Element access (mocked)
  MOCK_METHOD(reference, get, (size_type idx), ());
  MOCK_METHOD(const_reference, get, (size_type idx), (const));

  reference operator[](size_type idx) { return get(idx); }
  const_reference operator[](size_type idx) const { return get(idx); }
  reference at(size_type idx) { return get(idx); }
  const_reference at(size_type idx) const { return get(idx); }
  MOCK_METHOD(reference, front, (), ());
  MOCK_METHOD(const_reference, front, (), (const));
  MOCK_METHOD(reference, back, (), ());
  MOCK_METHOD(const_reference, back, (), (const));

 private:
  /** @brief Reference to the configuration backend */
  const ConfigBackend& config_;

  /** @brief Sample grid points for testing */
  std::vector<MockGridPoint> gridPoints_;
};

}  // namespace metada::backends::gmock