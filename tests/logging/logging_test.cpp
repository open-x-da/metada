#include <gtest/gtest.h>

#include "backends/logging/GoogleLogger.h"
#include "core/logging/Logger.h"

using Logger = metada::logging::Logger<metada::logging::GoogleLogger>;

class LoggingTest : public ::testing::Test {
 protected:
  void SetUp() override { metada::logging::GoogleLogger::Init("logging_test"); }

  void TearDown() override { metada::logging::GoogleLogger::Shutdown(); }

  Logger logger;
};

TEST_F(LoggingTest, BasicLogging) {
  EXPECT_NO_THROW(logger.Info("Info message"));
  EXPECT_NO_THROW(logger.Warning("Warning message"));
  EXPECT_NO_THROW(logger.Error("Error message"));
  EXPECT_NO_THROW(logger.Debug("Debug message"));
}