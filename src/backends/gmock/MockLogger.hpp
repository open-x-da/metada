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
 * @details This class provides a mock implementation of the logger interface using
 * Google Mock. It can be used to verify that logging calls are made with the
 * expected parameters and at the expected severity levels during tests.
 *
 * The mock implements the LoggerBackend concept and provides methods to set 
 * expectations on logging calls and verify that they were satisfied. It's designed
 * to work with the LogStream class for stream-based logging in tests.
 *
 * Example usage:
 * @code
 * MockLogger<MockConfig> mock_logger(config);
 *
 * // Test message logging at different severity levels
 * EXPECT_CALL(mock_logger, LogMessage(LogLevel::Info, "Starting test"));
 * EXPECT_CALL(mock_logger, LogMessage(LogLevel::Error, "Test failed"));
 *
 * // When using stream API, the eventual message string is what gets logged
 * // So you set expectations on the LogMessage method with the final string
 * EXPECT_CALL(mock_logger, LogMessage(LogLevel::Info, "User 123 logged in from 192.168.1.1"));
 *
 * // Then pass the mock to code under test
 * TestFunction(mock_logger);
 * @endcode
 * 
 * @tparam ConfigBackend The configuration backend type required by the logger
 */
template <typename ConfigBackend>
class MockLogger {
 public:
  /**
   * @brief Default constructor is disabled
   * 
   * @details Logger backends should always be initialized with a configuration.
   */
  MockLogger() = delete;

  /**
   * @brief Destructor that calls Shutdown
   * 
   * @details Ensures proper cleanup of any logger resources when the mock is destroyed.
   */
  ~MockLogger() { Shutdown(); }

  /**
   * @brief Copy constructor is disabled
   * 
   * @details Logger backends are not intended to be copied.
   */
  MockLogger(const MockLogger&) = delete;

  /**
   * @brief Copy assignment operator is disabled
   * 
   * @details Logger backends are not intended to be copied.
   */
  MockLogger& operator=(const MockLogger&) = delete;

  /**
   * @brief Move constructor - explicitly defined for Google Mock compatibility
   * 
   * @details Google Mock requires move operations to be defined for proper test setup.
   */
  MockLogger(MockLogger&&) noexcept {}

  /**
   * @brief Move assignment - explicitly defined for Google Mock compatibility
   * 
   * @details Google Mock requires move operations to be defined for proper test setup.
   * 
   * @return Reference to this MockLogger instance
   */
  MockLogger& operator=(MockLogger&&) noexcept { return *this; }

  /**
   * @brief Constructor that takes a configuration backend
   * 
   * @details Initializes the mock logger with the provided configuration and
   * calls the static Init method with a default application name.
   *
   * @param config The configuration backend instance
   */
  explicit MockLogger([[maybe_unused]] const ConfigBackend& config) {
    Init("MockLogger");
  }

  /**
   * @brief Mock method for logging messages
   * 
   * @details This is the primary method that will be mocked in tests. It receives
   * log messages with their associated severity level. When using the stream API,
   * this method receives the final formatted message string.
   *
   * @param level The severity level of the log message
   * @param message The formatted log message string
   */
  MOCK_METHOD(void, LogMessage, (LogLevel level, const std::string& message));

  /**
   * @brief No-op initialization for testing
   * 
   * @details In a real logger, this would initialize the logging system.
   * For the mock, it's a no-op that can be used to verify initialization calls.
   *
   * @param app_name The name of the application for logging identification
   */
  static void Init([[maybe_unused]] const std::string& app_name) {}

  /**
   * @brief No-op shutdown for testing
   * 
   * @details In a real logger, this would clean up logging resources.
   * For the mock, it's a no-op that can be used to verify shutdown calls.
   */
  static void Shutdown() {}
};

}  // namespace metada::backends::gmock
