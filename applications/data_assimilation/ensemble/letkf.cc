#include "Logger.h"
#include "google/GoogleLogger.h"

namespace {
// Logger type alias for cleaner code
namespace glog = metada::backends::tools::logger::google_logger;
namespace mlog = metada::framework::tools::logger;

using Logger = mlog::Logger<glog::GoogleLogger>;

// RAII wrapper for logger initialization
class LoggerInitializer {
 public:
  explicit LoggerInitializer(const std::string& app_name) {
    glog::GoogleLogger::Init(app_name);
  }
  ~LoggerInitializer() { glog::GoogleLogger::Shutdown(); }
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
