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

#include "MockLogger.hpp"
#include "utils/logger/Logger.hpp"

namespace metada {
namespace framework {
namespace common {
namespace utils {
namespace logger {
namespace tests {

using ::testing::Return;

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

}  // namespace tests
}  // namespace logger
}  // namespace utils
}  // namespace common
}  // namespace framework
}  // namespace metada
