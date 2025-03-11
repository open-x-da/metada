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
 * - Move semantics (non-copyable behavior)
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "MockLogger.hpp"
#include "utils/logger/LogStream.hpp"
#include "utils/logger/Logger.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;

using framework::Logger;

/**
 * @brief Test fixture for Logger class tests
 *
 * Provides a test environment with a Logger instance using MockLogger backend.
 * Enables consistent testing of Logger's delegation behavior.
 */
class LoggerTest : public ::testing::Test {
 protected:
  /** @brief Logger instance with mock backend */
  Logger<MockLogger> logger;

  void SetUp() override { MockLogger::Init("test"); }

  void TearDown() override { MockLogger::Shutdown(); }
};

/**
 * @brief Verify Info method delegation
 *
 * Tests that Logger::Info() stream properly delegates to backend's LogMessage
 * method:
 * - Correct message forwarding
 * - Correct log level
 * - Single invocation per call
 */
TEST_F(LoggerTest, InfoStreamCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), LogMessage(LogLevel::Info, "test message"))
      .Times(1);
  logger.Info() << "test message";
}

/**
 * @brief Verify Warning method delegation
 *
 * Tests that Logger::Warning() stream properly delegates to backend's
 * LogMessage method:
 * - Correct message forwarding
 * - Correct log level
 * - Single invocation per call
 */
TEST_F(LoggerTest, WarningStreamCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), LogMessage(LogLevel::Warning, "test message"))
      .Times(1);
  logger.Warning() << "test message";
}

/**
 * @brief Verify Error method delegation
 *
 * Tests that Logger::Error() stream properly delegates to backend's LogMessage
 * method:
 * - Correct message forwarding
 * - Correct log level
 * - Single invocation per call
 */
TEST_F(LoggerTest, ErrorStreamCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), LogMessage(LogLevel::Error, "test message"))
      .Times(1);
  logger.Error() << "test message";
}

/**
 * @brief Verify Debug method delegation
 *
 * Tests that Logger::Debug() stream properly delegates to backend's LogMessage
 * method:
 * - Correct message forwarding
 * - Correct log level
 * - Single invocation per call
 */
TEST_F(LoggerTest, DebugStreamCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), LogMessage(LogLevel::Debug, "test message"))
      .Times(1);
  logger.Debug() << "test message";
}

/**
 * @brief Verify move semantics and non-copyable behavior
 *
 * Tests that Logger properly supports:
 * - Move construction
 * - Move assignment
 * - Non-copyable behavior (verified at compile time)
 */
TEST_F(LoggerTest, MoveSemantics) {
  // Setup expectations for the original logger
  EXPECT_CALL(logger.backend(), LogMessage(LogLevel::Info, "original logger"))
      .Times(1);
  logger.Info() << "original logger";

  // Test move construction
  Logger<MockLogger> moved_logger(std::move(logger));

  // Setup expectations for the moved logger
  EXPECT_CALL(moved_logger.backend(),
              LogMessage(LogLevel::Info, "moved logger"))
      .Times(1);
  moved_logger.Info() << "moved logger";

  // Create a new logger for move assignment test
  Logger<MockLogger> another_logger;
  EXPECT_CALL(another_logger.backend(),
              LogMessage(LogLevel::Info, "another logger"))
      .Times(1);
  another_logger.Info() << "another logger";

  // Test move assignment
  Logger<MockLogger> assigned_logger;
  assigned_logger = std::move(another_logger);

  // Setup expectations for the assigned logger
  EXPECT_CALL(assigned_logger.backend(),
              LogMessage(LogLevel::Info, "assigned logger"))
      .Times(1);
  assigned_logger.Info() << "assigned logger";

  // Note: The following would not compile due to deleted copy operations:
  // Logger<MockLogger> copy_logger(logger); // Copy construction - deleted
  // Logger<MockLogger> assign_logger;
  // assign_logger = logger; // Copy assignment - deleted
}

/**
 * @brief Test complex stream-based logging
 *
 * Verifies that the stream-based API correctly composes messages with
 * multiple data types and forwards them properly.
 */
TEST_F(LoggerTest, StreamBasedLoggingWithMultipleTypes) {
  // The backend's LogMessage method should be called with the composed message
  EXPECT_CALL(logger.backend(), LogMessage(LogLevel::Info, "test 42 message"))
      .Times(1);

  // Use the stream-based API to compose a message
  logger.Info() << "test " << 42 << " message";
}

/**
 * @brief Test logging with variable types
 *
 * Verifies that the stream-based API can handle different variable types
 * including numeric types, boolean values, and strings.
 */
TEST_F(LoggerTest, StreamBasedLoggingWithVariableTypes) {
  // Test with different types
  EXPECT_CALL(logger.backend(),
              LogMessage(LogLevel::Warning, "warning 3.14 value"))
      .Times(1);
  logger.Warning() << "warning " << 3.14 << " value";

  EXPECT_CALL(logger.backend(),
              LogMessage(LogLevel::Error, "error: value=true"))
      .Times(1);
  logger.Error() << "error: value=" << true;

  EXPECT_CALL(logger.backend(),
              LogMessage(LogLevel::Debug, "debug complex = (1,2)"))
      .Times(1);
  std::stringstream ss;
  ss << "(" << 1 << "," << 2 << ")";
  logger.Debug() << "debug complex = " << ss.str();
}

}  // namespace metada::tests
