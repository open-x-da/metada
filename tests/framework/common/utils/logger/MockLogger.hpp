#pragma once

#include <gmock/gmock.h>

#include <string>

#include "utils/logger/ILogger.hpp"
#include "utils/logger/LogStream.hpp"

namespace metada::tests {

using framework::ILogger;
using framework::LogStream;

/**
 * @brief Mock implementation of ILogger for testing purposes
 *
 * This class provides a mock implementation of the ILogger interface using
 * Google Mock. It can be used to verify that logging calls are made with the
 * expected parameters and at the expected severity levels during tests.
 *
 * The mock provides methods to set expectations on logging calls and verify
 * that they were satisfied. It supports both the traditional string-based API
 * and the newer stream-based API.
 *
 * Example usage:
 * @code
 * MockLogger mock_logger;
 *
 * // Test traditional API
 * EXPECT_CALL(mock_logger, Info("Starting test"));
 * EXPECT_CALL(mock_logger, Error("Test failed"));
 *
 * // Test stream API (needs special handling)
 * // For this, you typically verify Info/Error etc. gets called with the
 * // final formatted string
 * EXPECT_CALL(mock_logger, Info("User 123 logged in from 192.168.1.1"));
 *
 * // Then pass the mock to code under test
 * TestFunction(mock_logger);
 * @endcode
 */
class MockLogger : public ILogger {
 public:
  /**
   * @brief Default constructor
   */
  MockLogger() = default;

  /**
   * @brief Move constructor (no-op for mock)
   */
  MockLogger(MockLogger&&) noexcept {}

  /**
   * @brief Move assignment (no-op for mock)
   */
  MockLogger& operator=(MockLogger&&) noexcept { return *this; }

  // Mock methods for string-based API
  MOCK_METHOD(void, Info, (const std::string& message), (override));
  MOCK_METHOD(void, Warning, (const std::string& message), (override));
  MOCK_METHOD(void, Error, (const std::string& message), (override));
  MOCK_METHOD(void, Debug, (const std::string& message), (override));

  // Methods for stream-based API
  // Note: These use the base class implementations which ultimately call the
  // mocked methods
  LogStream InfoStream() override { return framework::ILogger::InfoStream(); }
  LogStream WarningStream() override {
    return framework::ILogger::WarningStream();
  }
  LogStream ErrorStream() override { return framework::ILogger::ErrorStream(); }
  LogStream DebugStream() override { return framework::ILogger::DebugStream(); }

  /**
   * @brief No-op initialization for testing
   */
  static void Init([[maybe_unused]] const std::string& app_name) {}

  /**
   * @brief No-op shutdown for testing
   */
  static void Shutdown() {}
};

}  // namespace metada::tests
