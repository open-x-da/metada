#include "backends/logging/google/GoogleLogger.h"
#include "core/logging/Logger.h"

using Logger = metada::logging::Logger<metada::logging::GoogleLogger>;

int main() {
  Logger logger;
  metada::logging::GoogleLogger::Init("my_app");

  logger.Info("Application started");
  // ... your code ...
  logger.Info("Application finished");

  metada::logging::GoogleLogger::Shutdown();
  return 0;
}
