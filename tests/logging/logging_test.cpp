#include <gtest/gtest.h>

#include "core/tools/logging/Logger.h"

#ifdef USE_GLOG
#include "backends/tools/logging/google/GoogleLogger.h"
using LoggerBackend = metada::logging::GoogleLogger;
#else
#include "backends/tools/logging/default/DefaultLogger.h"
using LoggerBackend = metada::logging::DefaultLogger;
#endif

using Logger = metada::logging::Logger<LoggerBackend>;

class LoggingTest : public ::testing::Test {
 protected:
  void SetUp() override { LoggerBackend::Init("logging_test"); }

  void TearDown() override { LoggerBackend::Shutdown(); }

  Logger logger;
};

TEST_F(LoggingTest, BasicLogging) {
  EXPECT_NO_THROW(logger.Info("Info message"));
  EXPECT_NO_THROW(logger.Warning("Warning message"));
  EXPECT_NO_THROW(logger.Error("Error message"));
  EXPECT_NO_THROW(logger.Debug("Debug message"));
}