/**
 * @file enkf.cpp
 * @brief Driver program for the Ensemble Kalman Filter (EnKF) data assimilation
 * algorithm
 * @details This application implements the ensemble Kalman filter data
 * assimilation algorithm. It reads configuration from a file, initializes
 * required components like ensemble members, observations and observation
 *          operators, and performs the analysis step.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: enkf <config_file>
 * @return 0 on success, 1 on failure
 */

#include "EnKF.hpp"

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "ControlVariableBackend.hpp"
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "IdentityControlVariableBackend.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "SimpleBackendTraits.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::SimpleBackendTag;

int main(int argc, char** argv) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: enkf <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info()
        << "Ensemble Kalman Filter data assimilation application starting...";

    // Get EnKF parameters from configuration
    auto analysis_config = config.GetSubsection("analysis");
    double inflation_factor = analysis_config.Get("inflation").asFloat();
    std::string inflation_method =
        analysis_config.Get("inflation_method").asString();

    logger.Info() << "EnKF parameters:";
    logger.Info() << "  Inflation factor: " << inflation_factor;
    logger.Info() << "  Inflation method: " << inflation_method;

    // Initialize geometry
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Initialize ensemble
    fwk::Ensemble<BackendTag> ensemble(config.GetSubsection("ensemble"),
                                       geometry);
    logger.Info() << "Ensemble initialized with " << ensemble.Size()
                  << " members";

    // Initialize observations
    fwk::Observation<BackendTag> obs(config.GetSubsection("observations"));
    logger.Info() << "Observations loaded: " << obs.size() << " observations";

    // Create identity control backend for ensemble methods
    auto control_backend =
        std::make_shared<fwk::IdentityControlVariableBackend<BackendTag>>();

    // Initialize observation operator
    fwk::ObsOperator<BackendTag> obs_op(config.GetSubsection("obs_operator"),
                                        *control_backend);
    logger.Info() << "Observation operator initialized";

    // Create EnKF algorithm instance
    fwk::EnKF<BackendTag> enkf(ensemble, obs, obs_op, config);

    logger.Info() << "Starting EnKF analysis...";

    // Perform EnKF analysis
    enkf.Analyse();

    // Get analysis results
    auto results = enkf.getAnalysisResults();

    // Log results
    logger.Info() << "=== EnKF Analysis Results ===";
    logger.Info() << "Innovation norm: " << results.innovation_norm;
    logger.Info() << "Analysis increment norm: "
                  << results.analysis_increment_norm;
    logger.Info() << "Background spread: " << results.background_spread;
    logger.Info() << "Analysis spread: " << results.analysis_spread;
    logger.Info() << "Max Kalman gain: " << results.max_kalman_gain;
    logger.Info() << "Min Kalman gain: " << results.min_kalman_gain;
    logger.Info() << "Condition number: " << results.condition_number;
    logger.Info() << "Ensemble size: " << results.ensemble_size;
    logger.Info() << "Observation count: " << results.observation_count;
    logger.Info() << "Inflation method: " << results.inflation_method;
    logger.Info() << "Inflation factor: " << results.inflation_factor;

    // Log spread statistics
    logger.Info() << "=== Spread Statistics ===";
    double spread_ratio = results.analysis_spread / results.background_spread;
    logger.Info() << "Background spread: " << results.background_spread;
    logger.Info() << "Analysis spread: " << results.analysis_spread;
    logger.Info() << "Spread ratio (analysis/background): " << spread_ratio;

    // Log Kalman gain statistics
    logger.Info() << "=== Kalman Gain Statistics ===";
    double gain_range = results.max_kalman_gain - results.min_kalman_gain;
    double gain_mean =
        (results.max_kalman_gain + results.min_kalman_gain) / 2.0;
    logger.Info() << "Max Kalman gain: " << results.max_kalman_gain;
    logger.Info() << "Min Kalman gain: " << results.min_kalman_gain;
    logger.Info() << "Gain range: " << gain_range;
    logger.Info() << "Gain mean: " << gain_mean;

    // Log numerical stability
    logger.Info() << "=== Numerical Stability ===";
    logger.Info() << "Condition number: " << results.condition_number;
    if (results.condition_number > 1e12) {
      logger.Warning()
          << "High condition number detected - potential numerical instability";
    } else if (results.condition_number > 1e8) {
      logger.Warning()
          << "Moderate condition number - monitor for numerical issues";
    } else {
      logger.Info() << "Good numerical stability";
    }

    // Save analysis results
    enkf.saveEnsemble();

    logger.Info()
        << "Ensemble Kalman Filter data assimilation application completed "
           "successfully";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "EnKF application failed: " << e.what() << std::endl;
    return 1;
  }
}