/**
 * @file letkf.cc
 * @brief Local Ensemble Transform Kalman Filter (LETKF) application
 */

#include "Config.h"
#include "JsonConfigTraits.h"
#include "Logger.h"
#include "console/ConsoleLoggerTraits.h"

using namespace metada::framework::tools::config;

// Use Config with appropriate backend traits
using JsonConfigType = Config<ConfigTraits<void>::ConfigBackend>;

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
 * @brief Load LETKF configuration from file
 * @param logger Logger instance for output
 * @param config_file Path to configuration file
 * @return true if configuration loaded successfully
 */
template <typename Logger>
bool LoadConfiguration(Logger& logger, const std::string& config_file) {
  try {
    JsonConfigType config;
    if (!config.LoadFromFile(config_file)) {
      logger.Error("Failed to load configuration from: " + config_file);
      return false;
    }
    logger.Info("Loaded configuration from: " + config_file);
    return true;
  } catch (const std::exception& e) {
    logger.Error("Error loading configuration: " + std::string(e.what()));
    return false;
  }
}

bool isYamlFile(const std::string& filename) {
  size_t pos = filename.find_last_of('.');
  if (pos == std::string::npos) return false;
  std::string ext = filename.substr(pos);
  return ext == ".yaml" || ext == ".yml";
}

/**
 * @brief Main entry point for LETKF application
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return 0 on success, 1 on failure
 */
int main(int argc, char* argv[]) {
  // Initialize logger with RAII
  LoggerInitializer log_init("letkf_app");

  metada::framework::tools::logger::Logger<
      metada::framework::tools::logger::LoggerTraits<void>::LoggerBackend>
      logger;

  logger.Info("LETKF application starting...");

  try {
    // Check command line arguments
    if (argc != 2) {
      logger.Error("Usage: letkf <config_file>");
      return 1;
    }

    // Load configuration
    std::string config_file = argv[1];
    if (!LoadConfiguration(logger, config_file)) {
      return 1;
    }

    // Example: Read LETKF parameters from configuration
    int ensemble_size;
    double inflation_factor;
    std::vector<std::string> state_variables;

    JsonConfigType config;
    if (!config.LoadFromFile(config_file)) {
      logger.Error("Failed to load configuration file: " + config_file);
      return 1;
    }

    try {
      ensemble_size = std::get<int>(config.Get("letkf.ensemble_size", 32));
      inflation_factor =
          std::get<double>(config.Get("letkf.inflation_factor", 1.1));
      state_variables = std::get<std::vector<std::string>>(config.Get(
          "letkf.state_variables",
          std::vector<std::string>{"temperature", "pressure", "humidity"}));
    } catch (const std::exception& e) {
      logger.Error("Error reading configuration values: " +
                   std::string(e.what()));
      return 1;
    }

    // Log configuration
    logger.Info("LETKF Configuration:");
    logger.Info("  - Ensemble Size: " + std::to_string(ensemble_size));
    logger.Info("  - Inflation Factor: " + std::to_string(inflation_factor));
    logger.Info("  - State Variables:");
    for (const auto& var : state_variables) {
      logger.Info("    * " + var);
    }

    // Add your LETKF implementation here
    logger.Debug("Initializing LETKF parameters");
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
