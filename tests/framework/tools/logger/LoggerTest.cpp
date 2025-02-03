#include <gtest/gtest.h>

#include "Logger.h"

#ifdef USE_GLOG
#include "google/GoogleLogger.h"
namespace mlogger = metada::backends::tools::logger::google_logger;
#else
#include "default/DefaultLogger.h"
namespace mlogger = metada::backends::tools::logger::default_logger;
#endif

namespace mtest {
using namespace metada::framework::tools::logger;
using Backend = mlogger::GoogleLogger;
using Logger = Logger<Backend>;

class LoggerTest : public ::testing::Test {
 protected:
  void SetUp() override { Backend::Init("LoggerTest"); }

  void TearDown() override { Backend::Shutdown(); }

  Logger logger;
};

TEST_F(LoggerTest, BasicLogging) {
  EXPECT_NO_THROW(logger.Info("Info message"));
  EXPECT_NO_THROW(logger.Warning("Warning message"));
  EXPECT_NO_THROW(logger.Error("Error message"));
  EXPECT_NO_THROW(logger.Debug("Debug message"));
}

}  // namespace mtest
