/**
 * @file main.cpp
 * @brief Driver for the Ensemble Transform Kalman Filter (ETKF)
 * application
 */

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "ETKF.hpp"
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "SimpleBackendTraits.hpp"
#include "State.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::SimpleBackendTag;

int main(int argc, char** argv) {
  // Initialize application context
  auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
  auto& logger = context.getLogger();
  auto& config = context.getConfig();

  logger.Info() << "ETKF application starting...";

  try {
    // Validate command line arguments
    if (argc != 2) {
      logger.Error() << "Usage: etkf <config_file>";
      return 1;
    }

    // Initialize components
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    fwk::Ensemble<BackendTag> ensemble(config.GetSubsection("ensemble"),
                                       geometry);
    ensemble.RecomputeMean();
    ensemble.RecomputePerturbations();

    // fwk::Observation<BackendTag> observation(config);
    // fwk::ObsOperator<BackendTag> obs_operator(config);

    // Run ETKF analysis
    // fwk::ETKF<BackendTag> etkf(ensemble, observation, obs_operator,
    //                           inflation_factor);
    // etkf.Analyse();

    logger.Info() << "ETKF application completed successfully";
    return 0;
  } catch (const std::exception& e) {
    logger.Error() << "ETKF application failed: " << e.what();
    return 1;
  }
}
