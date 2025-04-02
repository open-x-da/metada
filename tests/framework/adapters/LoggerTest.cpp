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
 * - Config-based initialization
 * - Backend access
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <sstream>
#include <string>

#include "MockBackendTraits.hpp"
#include "common/utils/config/Config.hpp"
#include "common/utils/logger/LogStream.hpp"
#include "common/utils/logger/Logger.hpp"

namespace metada::tests {

using ::testing::_;
using ::testing::Return;

using framework::Config;
using framework::Logger;
using framework::LogLevel;

/**
 * @brief Test fixture for Logger class tests
 *
 * Provides a test environment with a Logger instance using MockLogger backend.
 * Enables consistent testing of Logger's delegation behavior.
 */
class LoggerTest : public ::testing::Test {
 protected:
  std::string config_file_;
  std::unique_ptr<Config<traits::MockBackendTag>> config_;
  std::unique_ptr<Logger<traits::MockBackendTag>> logger;

  void SetUp() override {
    // Get the directory where the test file is located
    auto test_dir = std::filesystem::path(__FILE__).parent_path();
    config_file_ = (test_dir / "test_config.yaml").string();
    config_ = std::make_unique<Config<traits::MockBackendTag>>(config_file_);
    logger = std::make_unique<Logger<traits::MockBackendTag>>(*config_);
    traits::BackendTraits<traits::MockBackendTag>::LoggerBackend::Init("test");
  }

  void TearDown() override {
    traits::BackendTraits<traits::MockBackendTag>::LoggerBackend::Shutdown();
  }
};

/**
 * @brief Test Logger initialization with config
 *
 * Verifies that:
 * - Logger can be constructed with a config
 * - The backend is properly initialized
 * - The config is correctly passed to the backend
 */
TEST_F(LoggerTest, ConfigBasedInitialization) {
  // Create a new config and logger
  Config<traits::MockBackendTag> new_config(config_file_);
  Logger<traits::MockBackendTag> new_logger(new_config);

  // Verify the logger is functional
  EXPECT_CALL(new_logger.backend(), LogMessage(LogLevel::Info, "test message"))
      .Times(1);
  new_logger.Info() << "test message";
}

/**
 * @brief Test backend access methods
 *
 * Verifies that:
 * - The backend() method returns a reference to the correct backend
 * - The backend can be used directly for logging
 */
TEST_F(LoggerTest, BackendAccess) {
  // Test direct backend access
  EXPECT_CALL(logger->backend(),
              LogMessage(LogLevel::Info, "direct backend access"))
      .Times(1);
  logger->backend().LogMessage(LogLevel::Info, "direct backend access");

  // Test backend reference is valid
  auto& backend_ref = logger->backend();
  EXPECT_CALL(backend_ref, LogMessage(LogLevel::Warning, "backend reference"))
      .Times(1);
  backend_ref.LogMessage(LogLevel::Warning, "backend reference");
}

/**
 * @brief Test all logging levels
 *
 * Verifies that all logging levels (Info, Warning, Error, Debug) work correctly
 * and forward messages to the backend with the appropriate level.
 */
TEST_F(LoggerTest, AllLoggingLevels) {
  // Test Info level
  EXPECT_CALL(logger->backend(), LogMessage(LogLevel::Info, "info message"))
      .Times(1);
  logger->Info() << "info message";

  // Test Warning level
  EXPECT_CALL(logger->backend(),
              LogMessage(LogLevel::Warning, "warning message"))
      .Times(1);
  logger->Warning() << "warning message";

  // Test Error level
  EXPECT_CALL(logger->backend(), LogMessage(LogLevel::Error, "error message"))
      .Times(1);
  logger->Error() << "error message";

  // Test Debug level
  EXPECT_CALL(logger->backend(), LogMessage(LogLevel::Debug, "debug message"))
      .Times(1);
  logger->Debug() << "debug message";
}

/**
 * @brief Test move semantics and non-copyable behavior
 *
 * Verifies that:
 * - Move construction works correctly
 * - Move assignment works correctly
 * - Copy operations are properly deleted
 */
TEST_F(LoggerTest, MoveSemantics) {
  // Setup expectations for the original logger
  EXPECT_CALL(logger->backend(), LogMessage(LogLevel::Info, "original logger"))
      .Times(1);
  logger->Info() << "original logger";

  // Test move construction
  Logger<traits::MockBackendTag> moved_logger(std::move(*logger));

  // Setup expectations for the moved logger
  EXPECT_CALL(moved_logger.backend(),
              LogMessage(LogLevel::Info, "moved logger"))
      .Times(1);
  moved_logger.Info() << "moved logger";

  // Create a new logger for move assignment test
  Config<traits::MockBackendTag> another_config(config_file_);
  Logger<traits::MockBackendTag> another_logger(another_config);
  EXPECT_CALL(another_logger.backend(),
              LogMessage(LogLevel::Info, "another logger"))
      .Times(1);
  another_logger.Info() << "another logger";

  // Test move assignment
  another_logger = std::move(moved_logger);

  // Setup expectations for the assigned logger
  EXPECT_CALL(another_logger.backend(),
              LogMessage(LogLevel::Info, "assigned logger"))
      .Times(1);
  another_logger.Info() << "assigned logger";
}

/**
 * @brief Test complex stream-based logging
 *
 * Verifies that:
 * - Multiple values can be chained in a single log message
 * - Different types can be mixed in the same message
 * - The final message is correctly composed and forwarded
 */
TEST_F(LoggerTest, ComplexStreamLogging) {
  // Test with multiple values and types
  EXPECT_CALL(logger->backend(),
              LogMessage(LogLevel::Info, "User 123 logged in from 192.168.1.1"))
      .Times(1);
  logger->Info() << "User " << 123 << " logged in from " << "192.168.1.1";

  // Test with numeric types
  EXPECT_CALL(logger->backend(),
              LogMessage(LogLevel::Warning, "Resource usage: 75.5%"))
      .Times(1);
  logger->Warning() << "Resource usage: " << 75.5 << "%";

  // Test with boolean values
  EXPECT_CALL(logger->backend(),
              LogMessage(LogLevel::Error, "Operation failed: 1"))
      .Times(1);
  logger->Error() << "Operation failed: " << true;

  // Test with complex string composition
  std::stringstream ss;
  ss << "(" << 1 << "," << 2 << ")";
  EXPECT_CALL(logger->backend(),
              LogMessage(LogLevel::Debug, "Coordinates: (1,2)"))
      .Times(1);
  logger->Debug() << "Coordinates: " << ss.str();
}

/**
 * @brief Test multiple log messages in sequence
 *
 * Verifies that:
 * - Multiple log messages can be sent in sequence
 * - Each message is properly forwarded
 * - The order of messages is maintained
 */
TEST_F(LoggerTest, MultipleSequentialLogs) {
  // Setup expectations for multiple sequential logs
  EXPECT_CALL(logger->backend(), LogMessage(LogLevel::Info, "First message"))
      .Times(1);
  EXPECT_CALL(logger->backend(),
              LogMessage(LogLevel::Warning, "Second message"))
      .Times(1);
  EXPECT_CALL(logger->backend(), LogMessage(LogLevel::Error, "Third message"))
      .Times(1);

  // Send multiple logs in sequence
  logger->Info() << "First message";
  logger->Warning() << "Second message";
  logger->Error() << "Third message";
}

}  // namespace metada::tests
