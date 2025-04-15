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

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of a grid point for testing
 * 
 * @details
 * This simple struct represents a grid point in the mock geometry.
 * It provides equality comparison operators to facilitate testing
 * of iterator behavior and grid point access.
 */
struct MockGridPoint {
    int x, y, z;
    
    bool operator==(const MockGridPoint& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
    
    bool operator!=(const MockGridPoint& other) const {
        return !(*this == other);
    }
};

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
 * - Proper move semantics (non-copyable)
 *
 * @tparam ConfigBackend The mock configuration backend type
 *
 * @see MockGeometry
 * @see GeometryIterator
 * @see GeometryIteratorBackendType
 */
template <typename ConfigBackend>
class MockGeometryIterator {
public:
    MockGeometryIterator() = delete;

    // Delete copy and move constructors and assignment operators
    MockGeometryIterator(const MockGeometryIterator&) = delete;
    MockGeometryIterator& operator=(const MockGeometryIterator&) = delete;

    MockGeometryIterator(MockGeometryIterator&& other) noexcept : config_(std::move(other.config_)) {}
    MockGeometryIterator& operator=(MockGeometryIterator&& other) noexcept {
        if (this != &other) {
            config_ = std::move(other.config_);
        }
        return *this;
    }
    
    /** 
     * @brief Constructor from a config 
     * @param config Reference to the configuration backend
     */
    MockGeometryIterator(const ConfigBackend& config) : config_(config) {}

    // GMock objects are not copyable or movable, so define our own comparison behavior
    // without trying to mock operator== directly
    MOCK_METHOD(MockGridPoint, dereference, (), (const));
    MOCK_METHOD(void, increment, ());
    MOCK_METHOD(bool, compare, (const MockGeometryIterator&), (const));
    
    /**
     * @brief Dereference operator to access the current grid point
     * @return The current grid point
     */
    MockGridPoint operator*() const { 
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

private:
    /** @brief Reference to the configuration backend */
    const ConfigBackend& config_;
};

} // namespace metada::backends::gmock 