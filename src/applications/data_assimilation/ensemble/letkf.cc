#include "backends/tools/logger/google/GoogleLogger.h"
#include "core/tools/logger/Logger.h"

using Logger = metada::logger::Logger<metada::logger::GoogleLogger>;

int main() {
  Logger logger;
  metada::logger::GoogleLogger::Init("my_app");

  logger.Info("Application started");
  // ... your code ...
  logger.Info("Application finished");

  metada::logger::GoogleLogger::Shutdown();
  return 0;
}
