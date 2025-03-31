#pragma once

#include <gmock/gmock.h>

#include <string>

#include "common/utils/logger/LogStream.hpp"

namespace metada::backends::gmock {

using framework::LogLevel;
using framework::LogStream;

/**
 * @brief Mock implementation of logger for testing purposes
 *
 * This class provides a mock implementation of the logger interface using
 * Google Mock. It can be used to verify that logging calls are made with the
 * expected parameters and at the expected severity levels during tests.
 *
 * The mock provides methods to set expectations on logging calls and verify
 * that they were satisfied.
 *
 * Example usage:
 * @code
 * MockLogger mock_logger;
 *
 * // Test message logging at different severity levels
 * EXPECT_CALL(mock_logger, LogMessage(LogLevel::Info, "Starting test"));
 * EXPECT_CALL(mock_logger, LogMessage(LogLevel::Error, "Test failed"));
 *
 * // When using stream API, the eventual message string is what gets logged
 * // So you set expectations on the LogMessage method with the final string
 * EXPECT_CALL(mock_logger, LogMessage(LogLevel::Info, "User 123 logged in from
 * 192.168.1.1"));
 *
 * // Then pass the mock to code under test
 * TestFunction(mock_logger);
 * @endcode
 */
class MockLogger {
 public:
  /**
   * @brief Disabled default constructor
   */
  MockLogger() = default;

  /**
   * @brief Default destructor
   */
  ~MockLogger() = default;
  /**
   * @brief Copy constructor (no-op for mock)
   */
  MockLogger(const MockLogger&) = delete;

  /**
   * @brief Copy assignment (no-op for mock)
   */
  MockLogger& operator=(const MockLogger&) = delete;

  /**
   * @brief Move constructor (no-op for mock)
   */
  MockLogger(MockLogger&&) noexcept {}

  /**
   * @brief Move assignment (no-op for mock)
   */
  MockLogger& operator=(MockLogger&&) noexcept { return *this; }

    /**
   * @brief Mock method for logging messages
   * @param[in] level Log level
   * @param[in] message Log message
   */
  MOCK_METHOD(void, LogMessage, (LogLevel level, const std::string& message));

  /**
   * @brief No-op initialization for testing
   */
  static void Init([[maybe_unused]] const std::string& app_name) {}

  /**
   * @brief No-op shutdown for testing
   */
  static void Shutdown() {}
};

}  // namespace metada::backends::gmock
