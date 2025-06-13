/**
 * @file main.cpp
 * @brief Driver for the Local Ensemble Transform Kalman Filter (LETKF)
 * application
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
  // Initialize application context
  auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
  auto& logger = context.getLogger();
  auto& config = context.getConfig();

  logger.Info() << "LETKF application starting...";

  try {
    // Validate command line arguments
    if (argc != 2) {
      logger.Error() << "Usage: letkf <config_file>";
      return 1;
    }

    // Read configuration
    int ensemble_size = config.Get("ensemble_size").asInt();
    double inflation_factor = config.Get("inflation_factor").asFloat();
    auto state_variables = config.Get("state_variables").asVectorString();

    // Log configuration
    logger.Info() << "LETKF Configuration:";
    logger.Info() << "  - Ensemble Size: " << ensemble_size;
    logger.Info() << "  - Inflation Factor: " << inflation_factor;
    logger.Info() << "  - State Variables: " << state_variables.size();

    // Initialize components
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));
    fwk::Ensemble<BackendTag> ensemble(config, geometry);
    fwk::Observation<BackendTag> observation(config);
    fwk::ObsOperator<BackendTag> obs_operator(config);

    // Run LETKF analysis
    fwk::LETKF<BackendTag> letkf(ensemble, observation, obs_operator,
                                 inflation_factor);
    letkf.Analyse();

    logger.Info() << "LETKF application completed successfully";
    return 0;
  } catch (const std::exception& e) {
    logger.Error() << "LETKF application failed: " << e.what();
    return 1;
  }
}
