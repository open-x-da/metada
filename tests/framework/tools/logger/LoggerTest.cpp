/**
 * @file LoggerTest.cpp
 * @brief Unit tests for the Logger class template implementation
 *
 * This test suite verifies that the Logger class template correctly forwards
 * logging calls to the underlying logger backend implementation. It uses a mock
 * logger backend to verify the interaction between the Logger class and its
 * backend.
 *
 * The tests ensure that:
 * - Each logging method (Info, Warning, Error, Debug) properly delegates to the backend
 * - Messages are passed through unmodified
 * - Methods are called exactly once per invocation
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ILogger.h"
#include "Logger.h"

/**
 * @brief Mock logger class for testing Logger implementation
 *
 * This mock class implements the ILogger interface using Google Mock to provide
 * mock implementations of the logging methods for testing purposes.
 */
class MockLogger : public metada::framework::tools::logger::ILogger {
 public:
  /** @copydoc ILogger::Info */
  MOCK_METHOD(void, Info, (const std::string& message), (override));

  /** @copydoc ILogger::Warning */
  MOCK_METHOD(void, Warning, (const std::string& message), (override));

  /** @copydoc ILogger::Error */
  MOCK_METHOD(void, Error, (const std::string& message), (override));

  /** @copydoc ILogger::Debug */
  MOCK_METHOD(void, Debug, (const std::string& message), (override));

  /**
   * @brief Mock initialization method
   * @param app_name Unused application name parameter
   */
  static void Init([[maybe_unused]] const std::string& app_name) {}

  /**
   * @brief Mock shutdown method
   */
  static void Shutdown() {}
};

/**
 * @brief Test fixture for Logger tests
 *
 * This fixture class provides a common setup for testing the Logger template
 * with a mock logging backend.
 */
class LoggerTest : public ::testing::Test {
 protected:
  /** @brief Logger instance using MockLogger as backend for testing */
  metada::framework::tools::logger::Logger<MockLogger> logger;
};

/**
 * @brief Test that Info() calls underlying implementation
 *
 * Verifies that calling Info() on the Logger forwards the message
 * to the backend's Info() method exactly once with the correct message.
 *
 * Test Steps:
 * 1. Set expectation for Info() call with specific message
 * 2. Call logger.Info() with test message
 * 3. Verify expectation is met on test completion
 */
TEST_F(LoggerTest, InfoCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Info("test message")).Times(1);
  logger.Info("test message");
}

/**
 * @brief Test that Warning() calls underlying implementation
 *
 * Verifies that calling Warning() on the Logger forwards the message
 * to the backend's Warning() method exactly once with the correct message.
 *
 * Test Steps:
 * 1. Set expectation for Warning() call with specific message
 * 2. Call logger.Warning() with test message
 * 3. Verify expectation is met on test completion
 */
TEST_F(LoggerTest, WarningCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Warning("test message")).Times(1);
  logger.Warning("test message");
}

/**
 * @brief Test that Error() calls underlying implementation
 *
 * Verifies that calling Error() on the Logger forwards the message
 * to the backend's Error() method exactly once with the correct message.
 *
 * Test Steps:
 * 1. Set expectation for Error() call with specific message
 * 2. Call logger.Error() with test message
 * 3. Verify expectation is met on test completion
 */
TEST_F(LoggerTest, ErrorCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Error("test message")).Times(1);
  logger.Error("test message");
}

/**
 * @brief Test that Debug() calls underlying implementation
 *
 * Verifies that calling Debug() on the Logger forwards the message
 * to the backend's Debug() method exactly once with the correct message.
 *
 * Test Steps:
 * 1. Set expectation for Debug() call with specific message
 * 2. Call logger.Debug() with test message
 * 3. Verify expectation is met on test completion
 */
TEST_F(LoggerTest, DebugCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Debug("test message")).Times(1);
  logger.Debug("test message");
}
