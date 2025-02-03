#include "Logger.h"
#include "LoggerConfig.h"

namespace {
namespace mlogger = metada::framework::tools::logger;
using Logger = mlogger::Logger<metada::logger::LoggerBackend>;

class LoggerInitializer {
 public:
  explicit LoggerInitializer(const std::string& app_name) {
    metada::logger::LoggerBackend::Init(app_name);
  }
  ~LoggerInitializer() { metada::logger::LoggerBackend::Shutdown(); }
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
