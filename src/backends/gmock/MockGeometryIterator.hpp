/**
 * @file MockGeometryIterator.hpp
 * @brief Mock implementation of geometry iterator backend for testing
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This file provides a Google Mock implementation of a geometry iterator backend
 * for use in unit tests. It implements all interfaces required by the
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
#include <memory>
#include <vector>

#include "MockGridPoint.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of a geometry iterator backend for testing
 *
 * @details
 * This class provides a Google Mock implementation of a geometry iterator backend
 * for use in unit tests. It implements all methods required by the
 * GeometryIteratorBackendType concept and allows tests to set expectations on
 * iterator behavior.
 *
 * The mock iterator supports:
 * - Dereference operations to access grid points
 * - Increment operations to advance the iterator
 * - Comparison operations to check iterator positions
 * - Proper copy and move semantics
 *
 * @see MockGeometry
 * @see GeometryIterator
 * @see GeometryIteratorBackendType
 */
class MockGeometryIterator {
public:
    using GridPoint = MockGridPoint;
    MockGeometryIterator() = delete;
    
    // Constructor that takes a mock geometry or other context
    explicit MockGeometryIterator([[maybe_unused]] void* context) {}
    
    // No copy or move operations - Google Mock objects can't be copied or moved
    MockGeometryIterator(const MockGeometryIterator&) {};
    MockGeometryIterator& operator=(const MockGeometryIterator&) { return *this; }
    MockGeometryIterator(MockGeometryIterator&&) {};
    MockGeometryIterator& operator=(MockGeometryIterator&&) { return *this; }
    
    MOCK_METHOD(GridPoint, dereference, (), (const));
    MOCK_METHOD(void, increment, ());
    MOCK_METHOD(bool, compare, (const MockGeometryIterator&), (const));
    
    /**
     * @brief Dereference operator to access the current grid point
     * @return The current grid point
     */
    GridPoint operator*() const { 
        return dereference(); 
    }
    
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
     * @return Pointer to this iterator after advancement
     * @note This is a simplified implementation for testing purposes
     */
    MockGeometryIterator* operator++(int) {
        increment();
        return this;
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

} // namespace metada::backends::gmock 