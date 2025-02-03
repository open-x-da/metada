#include "backends/tools/logger/google/GoogleLogger.h"
#include "framework/tools/logger/Logger.h"

using Logger = metada::framework::tools::logger::Logger<
    metada::backends::tools::logger::google_logger::GoogleLogger>;

int main() {
  Logger logger;
  metada::backends::tools::logger::google_logger::GoogleLogger::Init("my_app");

  logger.Info("Application started");
  // ... your code ...
  logger.Info("Application finished");

  metada::backends::tools::logger::google_logger::GoogleLogger::Shutdown();
  return 0;
}
