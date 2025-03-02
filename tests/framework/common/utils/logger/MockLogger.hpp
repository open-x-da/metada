#pragma once

#include <gmock/gmock.h>

#include "utils/logger/ILogger.hpp"

namespace metada::tests {

using framework::ILogger;

/**
 * @brief Mock logger backend for testing
 *
 * Implements ILogger interface using Google Mock to provide mock logging
 * methods. Used to verify Logger's interaction with its backend implementation.
 *
 * This class inherits non-copyable behavior from ILogger but supports move
 * semantics to allow transferring ownership when needed.
 */
class MockLogger : public ILogger {
 public:
  /**
   * @brief Default constructor
   */
  MockLogger() = default;

  /**
   * @brief Move constructor - explicitly defined for Google Mock compatibility
   */
  MockLogger(MockLogger&&) noexcept {}

  /**
   * @brief Move assignment - explicitly defined for Google Mock compatibility
   */
  MockLogger& operator=(MockLogger&&) noexcept { return *this; }

  // Mock methods
  MOCK_METHOD(void, Info, (const std::string& message), (override));
  MOCK_METHOD(void, Warning, (const std::string& message), (override));
  MOCK_METHOD(void, Error, (const std::string& message), (override));
  MOCK_METHOD(void, Debug, (const std::string& message), (override));

  /**
   * @brief Mock initialization method
   * @param app_name Application name (unused in mock)
   */
  static void Init([[maybe_unused]] const std::string& app_name) {}

  /**
   * @brief Mock shutdown method
   */
  static void Shutdown() {}
};

}  // namespace metada::tests
