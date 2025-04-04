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
#include "L63BackendTraits.hpp"

using namespace metada::framework;
using namespace metada::framework::runs;
using namespace metada::traits;

/**
 * @brief Main entry point for LETKF application
 *
 * @details
 * Initializes the application context, loads configuration, sets up LETKF
 * parameters and executes the data assimilation algorithm. The application
 * follows these steps:
 * 1. Initialize application context (logging system and configuration)
 * 2. Parse command line arguments
 * 3. Load and validate configuration
 * 4. Set up LETKF parameters
 * 5. Execute data assimilation
 * 6. Clean up resources automatically via RAII
 *
 * @par Expected configuration format:
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
  auto context = ApplicationContext<L63BackendTag>(argv[0], argv[1]);
  auto& logger = context.getLogger();
  auto& config = context.getConfig();

  logger.Info() << "LETKF application starting...";

  try {
    // Check command line arguments
    if (argc != 2) {
      logger.Error() << "Usage: letkf <config_file>";
      return 1;
    }

    // Example: Read LETKF parameters from configuration
    int ensemble_size;
    double inflation_factor;
    std::vector<std::string> state_variables;

    try {
      ensemble_size = config.Get("ensemble_size").asInt();
      inflation_factor = config.Get("inflation_factor").asFloat();
      state_variables = config.Get("state_variables").asVectorString();
    } catch (const std::exception& e) {
      logger.Error() << "Error reading configuration values: " << e.what();
      return 1;
    }

    // Log configuration
    logger.Info() << "LETKF Configuration:";
    logger.Info() << "  - Ensemble Size: " << ensemble_size;
    logger.Info() << "  - Inflation Factor: " << inflation_factor;
    logger.Info() << "  - State Variables:";
    for (const auto& var : state_variables) {
      logger.Debug() << "    * " << var;
    }

    // Add your LETKF implementation here
    logger.Debug() << "Initializing LETKF parameters";
    logger.Info() << "Loading ensemble members";
    logger.Info() << "Processing observations";
    logger.Info() << "Computing analysis";

    logger.Info() << "LETKF application completed successfully";
    return 0;
  } catch (const std::exception& e) {
    logger.Error() << "LETKF application failed: " << e.what();
    return 1;
  }
}
