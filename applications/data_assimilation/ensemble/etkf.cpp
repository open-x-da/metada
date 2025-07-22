/**
 * @file main.cpp
 * @brief Driver program for the Ensemble Transform Kalman Filter (ETKF)
 * @details This application implements the ETKF data assimilation algorithm.
 *          It reads configuration from a file, initializes required components
 *          like ensemble members, observations and observation operators, and
 *          performs the analysis step.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: etkf <config_file>
 * @return 0 on success, 1 on failure
 */

#include "ETKF.hpp"

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "LiteBackendTraits.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::LiteBackendTag;

int main(int argc, char** argv) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: etkf <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "ETKF application starting...";

    // Initialize components
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    fwk::Ensemble<BackendTag> ensemble(config.GetSubsection("ensemble"),
                                       geometry);

    fwk::Observation<BackendTag> observations(
        config.GetSubsection("observations"));

    fwk::ObsOperator<BackendTag> obs_operator(
        config.GetSubsection("obs_operator"));

    // Run ETKF analysis
    fwk::ETKF<BackendTag> etkf(ensemble, observations, obs_operator,
                               config.GetSubsection("analysis"));
    etkf.Analyse();
    etkf.saveEnsemble();

    logger.Info() << "ETKF application completed successfully";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ETKF application failed: " << e.what() << std::endl;
    return 1;
  }
}
