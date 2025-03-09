/**
 * @file MockHardwareAccelerator.hpp
 * @brief Mock implementation of IHardwareAccelerator interface for testing
 * @ingroup tests
 * @author Metada Framework Team
 *
 * @details
 * This mock class provides a test double for the IHardwareAccelerator interface
 * using Google Mock. It allows testing code that depends on
 * IHardwareAccelerator by providing mock implementations of all interface
 * methods.
 *
 * @see IHardwareAccelerator
 * @see testing::Mock
 */

#pragma once

#include <gmock/gmock.h>

#include <string>
#include <vector>

#include "IHardwareAccelerator.hpp"

namespace metada::tests {

using framework::IHardwareAccelerator;

/**
 * @brief Mock implementation of IHardwareAccelerator for testing
 *
 * Provides mock methods for all hardware acceleration operations including:
 * - Device selection and management
 * - Memory management
 * - Execution control
 */
class MockHardwareAccelerator : public IHardwareAccelerator {
 public:
  // Device operations
  MOCK_METHOD(void, setDevice, (const std::string& device), (override));
  MOCK_METHOD(std::string, getDevice, (), (const, override));
  MOCK_METHOD(bool, supportsDeviceType, (const std::string& deviceType),
              (const, override));
  MOCK_METHOD(std::vector<std::string>, getAvailableDevices, (),
              (const, override));

  // Memory management
  MOCK_METHOD(size_t, getAvailableMemory, (), (const, override));
  MOCK_METHOD(void, setMemoryLimit, (size_t memoryLimitBytes), (override));
  MOCK_METHOD(size_t, getMemoryLimit, (), (const, override));

  // Execution control
  MOCK_METHOD(void, setAsynchronousExecution, (bool enable), (override));
  MOCK_METHOD(bool, isAsynchronousExecutionEnabled, (), (const, override));
  MOCK_METHOD(void, synchronize, (), (override));
};

}  // namespace metada::tests