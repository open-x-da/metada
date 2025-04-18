#pragma once

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

}  // namespace metada::backends::gmock 