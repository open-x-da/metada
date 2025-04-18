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

#include <memory>
#include <vector>
#include <gmock/gmock.h>

#include "MockGridPoint.hpp"
#include "MockGeometryIterator.hpp"

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
 * - Proper move semantics (non-copyable)
 *
 * @tparam ConfigBackend The mock configuration backend type
 *
 * @see MockGeometryIterator
 * @see GeometryBackendType
 */
template <typename ConfigBackend>
class MockGeometry {
public:
    using iterator = MockGeometryIterator;
    using const_iterator = MockGeometryIterator;

    /**
     * @brief Constructor that takes a mock config
     * @param config Mock configuration object
     */
    explicit MockGeometry(const ConfigBackend& config) : config_(config), gridPoints_({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}}) {}
    
    // Deleted default constructor
    MockGeometry() = delete;
    
    // Delete copy constructor and assignment operator
    MockGeometry(const MockGeometry&) = delete;
    MockGeometry& operator=(const MockGeometry&) = delete;

    // Move constructor and assignment operator
    MockGeometry(MockGeometry&& other) noexcept : config_(std::move(other.config_)), gridPoints_(std::move(other.gridPoints_)) {}
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
    MockGeometry clone() const {
        return MockGeometry(config_);
    }

    // Mock methods required by GeometryBackendType concept
    MOCK_METHOD(iterator, begin, ());
    MOCK_METHOD(iterator, end, ());
    MOCK_METHOD(const_iterator, begin, (), (const));
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
    const ConfigBackend& config_;
    std::vector<MockGridPoint> gridPoints_;
};

}    // namespace metada::backends::gmock 