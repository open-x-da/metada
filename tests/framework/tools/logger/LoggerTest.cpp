/**
 * @file LoggerTest.cpp
 * @brief Unit tests for the Logger class template
 *
 * This test suite verifies the Logger class template correctly delegates
 * logging operations to its backend implementation. It uses Google Mock to
 * create a mock logger backend and validate the interactions between Logger and
 * backend.
 *
 * Key test areas:
 * - Proper delegation of logging methods (Info, Warning, Error, Debug)
 * - Message passing integrity
 * - Call frequency verification
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ILogger.hpp"
#include "Logger.hpp"

/**
 * @brief Mock logger backend for testing
 *
 * Implements ILogger interface using Google Mock to provide mock logging
 * methods. Used to verify Logger's interaction with its backend implementation.
 */
class MockLogger : public metada::framework::tools::logger::ILogger {
 public:
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

/**
 * @brief Test fixture for Logger class tests
 *
 * Provides a test environment with a Logger instance using MockLogger backend.
 * Enables consistent testing of Logger's delegation behavior.
 */
class LoggerTest : public ::testing::Test {
 protected:
  /** @brief Logger instance with mock backend */
  metada::framework::tools::logger::Logger<MockLogger> logger;
};

/**
 * @brief Verify Info() method delegation
 *
 * Tests that Logger::Info() properly delegates to backend's Info() method:
 * - Correct message forwarding
 * - Single invocation per call
 */
TEST_F(LoggerTest, InfoCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Info("test message")).Times(1);
  logger.Info("test message");
}

/**
 * @brief Verify Warning() method delegation
 *
 * Tests that Logger::Warning() properly delegates to backend's Warning()
 * method:
 * - Correct message forwarding
 * - Single invocation per call
 */
TEST_F(LoggerTest, WarningCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Warning("test message")).Times(1);
  logger.Warning("test message");
}

/**
 * @brief Verify Error() method delegation
 *
 * Tests that Logger::Error() properly delegates to backend's Error() method:
 * - Correct message forwarding
 * - Single invocation per call
 */
TEST_F(LoggerTest, ErrorCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Error("test message")).Times(1);
  logger.Error("test message");
}

/**
 * @brief Verify Debug() method delegation
 *
 * Tests that Logger::Debug() properly delegates to backend's Debug() method:
 * - Correct message forwarding
 * - Single invocation per call
 */
TEST_F(LoggerTest, DebugCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Debug("test message")).Times(1);
  logger.Debug("test message");
}
