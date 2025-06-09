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
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "LETKF.hpp"
#include "MockBackendTraits.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::MockBackendTag;

int main(int argc, char** argv) {
  auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
  auto& logger = context.getLogger();
  auto& config = context.getConfig();

  logger.Info() << "LETKF application starting...";

  try {
    // Check command line arguments
    if (argc != 2) {
      logger.Error() << "Usage: letkf <config_file>";
      return 1;
    }

    // Read LETKF parameters from configuration
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
    logger.Info() << "  - State Variables: " << state_variables.size();
    for (const auto& var : state_variables) {
      logger.Debug() << "    * " << var;
    }

    // Initialize geometry from config
    logger.Info() << "Initializing geometry";
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Create Ensemble, Observation, ObsOperator using adapters
    fwk::Ensemble<BackendTag> ensemble(config, geometry);
    fwk::Observation<BackendTag> observation(config);
    fwk::ObsOperator<BackendTag> obs_operator(config);

    // Create LETKF and run analysis
    fwk::LETKF<BackendTag> letkf(ensemble, observation, obs_operator,
                                 inflation_factor);
    letkf.analyse();

    logger.Info() << "LETKF application completed successfully";
    return 0;
  } catch (const std::exception& e) {
    logger.Error() << "LETKF application failed: " << e.what();
    return 1;
  }
}
