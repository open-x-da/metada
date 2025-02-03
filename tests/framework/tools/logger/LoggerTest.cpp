#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ILogger.h"
#include "LoggerTestInterface.h"

namespace mlogger = metada::framework::tools::logger;

namespace mtest {

// Mock logger implementing ILogger interface
class MockLogger : public mlogger::ILogger {
 public:
  MOCK_METHOD(void, Info, (const std::string& message), (override));
  MOCK_METHOD(void, Warning, (const std::string& message), (override));
  MOCK_METHOD(void, Error, (const std::string& message), (override));
  MOCK_METHOD(void, Debug, (const std::string& message), (override));

  static void Init(const std::string&) {}
  static void Shutdown() {}
};

class LoggerTest : public ::testing::Test {
 protected:
  mlogger::LoggerTestInterface<MockLogger> logger;
};

TEST_F(LoggerTest, InfoCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Info("test message")).Times(1);
  logger.Info("test message");
}

TEST_F(LoggerTest, WarningCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Warning("test message")).Times(1);
  logger.Warning("test message");
}

TEST_F(LoggerTest, ErrorCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Error("test message")).Times(1);
  logger.Error("test message");
}

TEST_F(LoggerTest, DebugCallsUnderlyingImplementation) {
  EXPECT_CALL(logger.backend(), Debug("test message")).Times(1);
  logger.Debug("test message");
}

}  // namespace mtest
