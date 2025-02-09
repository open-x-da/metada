/**
 * @file letkf.cpp
 * @brief Local Ensemble Transform Kalman Filter (LETKF) application
 *
 * This application implements the Local Ensemble Transform Kalman Filter
 * algorithm for data assimilation. It reads configuration from a YAML/JSON file
 * and performs ensemble-based state estimation.
 */

#include "Config.hpp"
#include "ConfigBackendSelector.hpp"
#include "Logger.hpp"
#include "console/ConsoleLoggerTraits.hpp"

using namespace metada::framework::tools::config;

// Use Config with appropriate backend traits
using ConfigType = Config<ConfigTraits<void>::ConfigBackend>;

/**
 * @brief RAII wrapper for logger initialization and shutdown
 *
 * Provides automatic initialization and cleanup of the logger system
 * using RAII principles.
 */
class LoggerInitializer {
 public:
  /**
   * @brief Initialize logger with application name
   * @param app_name Name of the application for logging identification
   */
  explicit LoggerInitializer(const std::string& app_name) {
    metada::framework::tools::logger::LoggerTraits<void>::LoggerBackend::Init(
        app_name);
  }
  /**
   * @brief Shutdown logger on destruction
   *
   * Ensures proper cleanup of logger resources when object goes out of scope
   */
  ~LoggerInitializer() {
    metada::framework::tools::logger::LoggerTraits<
        void>::LoggerBackend::Shutdown();
  }
};

/**
 * @brief Load LETKF configuration from file
 *
 * Attempts to load and parse configuration settings from the specified file.
 * Handles both YAML and JSON formats.
 *
 * @param logger Logger instance for status and error reporting
 * @param config_file Path to configuration file
 * @return true if configuration loaded successfully, false otherwise
 */
template <typename Logger>
bool LoadConfiguration(Logger& logger, const std::string& config_file) {
  try {
    ConfigType config;
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

/**
 * @brief Check if a filename has a YAML extension
 *
 * @param filename Name of file to check
 * @return true if file has .yaml or .yml extension, false otherwise
 */
bool isYamlFile(const std::string& filename) {
  size_t pos = filename.find_last_of('.');
  if (pos == std::string::npos) return false;
  std::string ext = filename.substr(pos);
  return ext == ".yaml" || ext == ".yml";
}

/**
 * @brief Main entry point for LETKF application
 *
 * Initializes the application, loads configuration, sets up LETKF parameters
 * and executes the data assimilation algorithm.
 *
 * Expected configuration format:
 * ```yaml
 * letkf:
 *   ensemble_size: 32
 *   inflation_factor: 1.1
 *   state_variables:
 *     - temperature
 *     - pressure
 *     - humidity
 * ```
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments, expects config file path as first
 * argument
 * @return 0 on successful execution, 1 on error
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

    ConfigType config;
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
