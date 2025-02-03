#include "Logger.h"
#include "default/DefaultLogger.h"

namespace {
namespace mlogger = metada::framework::tools::logger;
using LoggerBackend =
    metada::backends::tools::logger::default_logger::DefaultLogger;
using Logger = mlogger::Logger<LoggerBackend>;

class LoggerInitializer {
 public:
  explicit LoggerInitializer(const std::string& app_name) {
    LoggerBackend::Init(app_name);
  }
  ~LoggerInitializer() { LoggerBackend::Shutdown(); }
};
}  // namespace

int main() {
  // Initialize logger with RAII
  LoggerInitializer log_init("letkf_app");
  Logger logger;

  logger.Info("LETKF application starting...");

  try {
    // Add your LETKF implementation here
    logger.Debug("Initializing LETKF parameters");

    // Example logging points you might add:
    logger.Info("Loading ensemble members");
    logger.Info("Processing observations");
    logger.Info("Computing analysis");

    logger.Info("LETKF application completed successfully");
    return 0;
  } catch (const std::exception& e) {
    logger.Error("LETKF application failed: " + std::string(e.what()));
    return 1;
  }
}
