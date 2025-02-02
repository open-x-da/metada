#include <gtest/gtest.h>

#include "Logger.h"

#ifdef USE_GLOG
#include "google/GoogleLogger.h"
using LoggerBackend = metada::logging::GoogleLogger;
#else
#include "default/DefaultLogger.h"
using LoggerBackend = metada::logging::DefaultLogger;
#endif

namespace metada {
namespace logging {
namespace test {

using Logger = metada::logging::Logger<LoggerBackend>;

class LoggerTest : public ::testing::Test {
 protected:
  void SetUp() override { LoggerBackend::Init("LoggerTest"); }

  void TearDown() override { LoggerBackend::Shutdown(); }

  Logger logger;
};

TEST_F(LoggerTest, BasicLogging) {
  EXPECT_NO_THROW(logger.Info("Info message"));
  EXPECT_NO_THROW(logger.Warning("Warning message"));
  EXPECT_NO_THROW(logger.Error("Error message"));
  EXPECT_NO_THROW(logger.Debug("Debug message"));
}

}  // namespace test
}  // namespace logging
}  // namespace metada
