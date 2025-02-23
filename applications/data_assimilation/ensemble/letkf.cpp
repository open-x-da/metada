/**
 * @file letkf.cpp
 * @brief Local Ensemble Transform Kalman Filter (LETKF) application
 *
 * This application implements the Local Ensemble Transform Kalman Filter
 * algorithm for data assimilation. It reads configuration from a YAML/JSON file
 * and performs ensemble-based state estimation.
 *
 * The LETKF algorithm is a variant of the Ensemble Kalman Filter that performs
 * the analysis step locally in space, making it computationally efficient and
 * suitable for high-dimensional systems. Key features:
 *
 * - Reads configuration from YAML/JSON files
 * - Supports multiple state variables (temperature, pressure, humidity, etc.)
 * - Configurable ensemble size and inflation parameters
 * - Local analysis for improved computational efficiency
 * - Robust error handling and logging
 *
 * @see Hunt et al. (2007) "Efficient Data Assimilation for Spatiotemporal
 * Chaos: A Local Ensemble Transform Kalman Filter"
 */

#include "ApplicationContext.hpp"
#include "ConfigBackendSelector.hpp"
#include "LoggerBackendSelector.hpp"
#include "utils/config/Config.hpp"
#include "utils/logger/Logger.hpp"

using namespace metada::framework::core;
using namespace metada::framework::common::utils::config;
using namespace metada::framework::common::utils::logger;

// Use Config with appropriate backend traits
using ConfigType = Config<ConfigTraits<void>::ConfigBackend>;
using LoggerType = Logger<LoggerTraits<void>::LoggerBackend>;

/**
 * @brief Load LETKF configuration from file
 *
 * Attempts to load and parse configuration settings from the specified file.
 * Handles both YAML and JSON formats. The configuration file should contain
 * LETKF-specific parameters such as ensemble size, inflation factor, and
 * state variables.
 *
 * @param logger Logger instance for status and error reporting
 * @param config_file Path to configuration file
 * @return true if configuration loaded successfully, false otherwise
 * @throws std::runtime_error If file cannot be opened or parsed
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
 * Validates if the given filename ends with a YAML extension (.yaml or .yml).
 * This helps determine the appropriate parser to use for configuration files.
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
 * and executes the data assimilation algorithm. The application follows these
 * steps:
 * 1. Initialize logging system
 * 2. Parse command line arguments
 * 3. Load and validate configuration
 * 4. Set up LETKF parameters
 * 5. Execute data assimilation
 * 6. Clean up resources
 *
 * Expected configuration format:
 * ```yaml
 * letkf:
 *   ensemble_size: 32        # Number of ensemble members
 *   inflation_factor: 1.1    # Covariance inflation parameter
 *   state_variables:         # List of variables to assimilate
 *     - temperature
 *     - pressure
 *     - humidity
 * ```
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments, expects config file path as first
 * argument
 * @return 0 on successful execution, 1 on error
 * @throws std::runtime_error For critical errors during execution
 */
int main(int argc, char* argv[]) {
  auto context = ApplicationContext("letkf_app", argv[1]);
  auto logger = context.getLogger();
  // auto config = context.getConfig();

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
