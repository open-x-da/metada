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

#include <memory>
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
  /** @brief Iterator type for traversing grid points */
  using iterator = MockGeometryIterator;

  /** @brief Const iterator type for traversing grid points */
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
      : config_(std::move(other.config_)),
        gridPoints_(std::move(other.gridPoints_)) {}

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

  /**
   * @brief Get the total number of grid points
   * @return Total number of grid points in the geometry
   */
  MOCK_METHOD(std::size_t, totalGridSize, (), (const));

  /**
   * @brief Check if the geometry is periodic in X dimension
   * @return True if periodic in X, false otherwise
   */
  MOCK_METHOD(bool, isPeriodicX, (), (const));

  /**
   * @brief Check if the geometry is periodic in Y dimension
   * @return True if periodic in Y, false otherwise
   */
  MOCK_METHOD(bool, isPeriodicY, (), (const));

  /**
   * @brief Check if the geometry is periodic in Z dimension
   * @return True if periodic in Z, false otherwise
   */
  MOCK_METHOD(bool, isPeriodicZ, (), (const));

  /**
   * @brief Check if the geometry is properly initialized
   * @return True if initialized, false otherwise
   */
  MOCK_METHOD(bool, isInitialized, (), (const));

  /**
   * @brief Perform halo exchange operation on a state
   * @param state The state object to perform halo exchange on
   * @tparam StateBackend The backend type of the state
   */
  template <typename StateBackend>
  void haloExchange(StateBackend& state) const {
    // Call the mock method with void* to avoid template issues in mocking
    haloExchangeImpl(static_cast<void*>(&state));
  }

  /**
   * @brief Implementation method for halo exchange (for mocking)
   * @param state Pointer to the state object (cast to void*)
   */
  MOCK_METHOD(void, haloExchangeImpl, (void* state), (const));

 private:
  /** @brief Reference to the configuration backend */
  const ConfigBackend& config_;

  /** @brief Sample grid points for testing */
  std::vector<MockGridPoint> gridPoints_;
};

}  // namespace metada::backends::gmock