/**
 * @file MockGeometryIterator.hpp
 * @brief Mock implementation of geometry iterator backend for testing
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This file provides a Google Mock implementation of a geometry iterator backend
 * for use in unit tests. It implements all interfaces required by the
 * GeometryIteratorBackendImpl concept.
 */

#pragma once

#include <gmock/gmock.h>
#include <memory>
#include <vector>

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of a grid point for testing
 * 
 * This simple struct represents a grid point in the mock geometry
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
 * GeometryIteratorBackendImpl concept and allows tests to set expectations on
 * iterator behavior.
 *
 * @tparam ConfigBackend The mock configuration backend type
 */
template <typename ConfigBackend>
class MockGeometryIterator {
public:
    MockGeometryIterator() = default;
    
    // GMock objects are not copyable or movable, so define our own comparison behavior
    // without trying to mock operator== directly
    MOCK_METHOD(MockGridPoint, dereference, (), (const));
    MOCK_METHOD(void, increment, ());
    MOCK_METHOD(bool, compare, (const MockGeometryIterator*), (const));
    
    // Operator implementations that delegate to the mockable methods
    MockGridPoint operator*() const { 
        return dereference(); 
    }
    
    MockGeometryIterator& operator++() {
        increment();
        return *this;
    }
    
    // Post-increment - return pointer to this object (safe approximation)
    MockGeometryIterator* operator++(int) {
        increment();
        return this;
    }
    
    bool operator==(const MockGeometryIterator& other) const {
        return compare(&other);
    }
    
    bool operator!=(const MockGeometryIterator& other) const {
        return !compare(&other);
    }
    
    // Helper method to mock a sequence of points for iteration
    void mockIteratorSequence(const std::vector<MockGridPoint>& points) {
        for (size_t i = 0; i < points.size(); ++i) {
            EXPECT_CALL(*this, dereference())
                .WillRepeatedly(testing::Return(points[i]));
            
            if (i < points.size() - 1) {
                EXPECT_CALL(*this, increment())
                    .WillOnce(testing::DoDefault());
            }
        }
    }
};

} // namespace metada::backends::gmock 