/**
 * @file MockGeometry.hpp
 * @brief Mock implementation of geometry backend for testing
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This file provides a Google Mock implementation of a geometry backend
 * for use in unit tests. It implements all interfaces required by the
 * GeometryBackendImpl concept.
 */

#pragma once

#include <memory>
#include <vector>
#include <gmock/gmock.h>

#include "MockGeometryIterator.hpp"

namespace metada::backends::gmock {

/**
 * @brief Mock implementation of a geometry backend for testing
 *
 * @details
 * This class provides a Google Mock implementation of a geometry backend
 * for use in unit tests. It implements all methods required by the
 * GeometryBackendImpl concept and allows tests to set expectations on
 * method calls.
 *
 * @tparam ConfigBackend The mock configuration backend type
 */
template <typename ConfigBackend>
class MockGeometry {
public:
    using iterator = MockGeometryIterator<ConfigBackend>;
    using const_iterator = MockGeometryIterator<ConfigBackend>;

    /**
     * @brief Constructor that takes a mock config
     * @param config Mock configuration object
     */
    explicit MockGeometry([[maybe_unused]] const ConfigBackend& config) {}
    
    // Deleted default constructor
    MockGeometry() = delete;
    
    // Delete copy constructor and assignment operator
    MockGeometry(const MockGeometry&) = delete;
    MockGeometry& operator=(const MockGeometry&) = delete;

    // Move constructor and assignment
    MockGeometry(MockGeometry&&) noexcept = default;
    MockGeometry& operator=(MockGeometry&&) noexcept = default;

    // Mock methods required by GeometryBackendImpl concept
    MOCK_METHOD(iterator, begin, ());
    MOCK_METHOD(iterator, end, ());
    MOCK_METHOD(const_iterator, begin, (), (const));
    MOCK_METHOD(const_iterator, end, (), (const));
    
    MOCK_METHOD(std::size_t, totalGridSize, (), (const));
    
    MOCK_METHOD(bool, isPeriodicX, (), (const));
    MOCK_METHOD(bool, isPeriodicY, (), (const));
    MOCK_METHOD(bool, isPeriodicZ, (), (const));
    
    MOCK_METHOD(bool, isInitialized, (), (const));
     
    MOCK_METHOD(MockGeometry, clone, (), (const));
    
    template <typename StateBackend>
    void haloExchange(StateBackend& state) {
        // Call the mock method with void* to avoid template issues in mocking
        haloExchangeImpl(static_cast<void*>(&state));
    }
    
    MOCK_METHOD(void, haloExchangeImpl, (void* state));
};

} // namespace metada::backends::gmock 