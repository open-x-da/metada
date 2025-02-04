/**
 * @file letkf.cc
 * @brief Local Ensemble Transform Kalman Filter (LETKF) application
 */

#include "Logger.h"
#include "LoggerTraits.h"
#include "console/ConsoleLogger.h"
// #include "google/GoogleLogger.h"

/**
 * @brief RAII wrapper for logger initialization and shutdown
 */
class LoggerInitializer {
 public:
  /**
   * @brief Initialize logger with application name
   * @param app_name Name of the application
   */
  explicit LoggerInitializer(const std::string& app_name) {
    metada::framework::tools::logger::LoggerTraits<void>::LoggerBackend::Init(
        app_name);
  }
  /**
   * @brief Shutdown logger on destruction
   */
  ~LoggerInitializer() {
    metada::framework::tools::logger::LoggerTraits<
        void>::LoggerBackend::Shutdown();
  }
};

/**
 * @brief Main entry point for LETKF application
 * @return 0 on success, 1 on failure
 */
int main() {
  // Initialize logger with RAII
  LoggerInitializer log_init("letkf_app");

  metada::framework::tools::logger::Logger<
      metada::framework::tools::logger::LoggerTraits<void>::LoggerBackend>
      logger;

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
