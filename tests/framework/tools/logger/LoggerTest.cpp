/**
 * @file LoggerTest.cpp
 * @brief Unit tests for Logger class
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ILogger.h"
#include "Logger.h"

/**
 * @brief Mock logger class for testing Logger implementation
 */
class MockLogger : public metada::ILogger {
 public:
  MOCK_METHOD(void, Info, (const std::string& message), (override));
  MOCK_METHOD(void, Warning, (const std::string& message), (override));
  MOCK_METHOD(void, Error, (const std::string& message), (override));
  MOCK_METHOD(void, Debug, (const std::string& message), (override));

  static void Init(const std::string&) {}
  static void Shutdown() {}
};

/**
 * @brief Test fixture for Logger tests
 */
class LoggerTest : public ::testing::Test {
 protected:
  metada::Logger<MockLogger> logger;
};

/**
 * @brief Test that Info() calls underlying implementation
 */
TEST_F(LoggerTest, InfoCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Info("test message")).Times(1);
  logger.Info("test message");
}

/**
 * @brief Test that Warning() calls underlying implementation
 */
TEST_F(LoggerTest, WarningCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Warning("test message")).Times(1);
  logger.Warning("test message");
}

/**
 * @brief Test that Error() calls underlying implementation
 */
TEST_F(LoggerTest, ErrorCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Error("test message")).Times(1);
  logger.Error("test message");
}

/**
 * @brief Test that Debug() calls underlying implementation
 */
TEST_F(LoggerTest, DebugCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Debug("test message")).Times(1);
  logger.Debug("test message");
}
