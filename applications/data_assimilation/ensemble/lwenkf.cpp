/**
 * @file lwenkf.cpp
 * @brief Driver program for the Local Weighted Ensemble Kalman Filter (LWEnKF)
 * data assimilation algorithm
 * @details This application implements the local weighted ensemble Kalman
 * filter data assimilation algorithm. It reads configuration from a file,
 * initializes required components like ensemble members, observations and
 * observation operators, and performs the analysis step with localization and
 * weighting features.
 *
 * The LWEnKF extends the traditional EnKF by:
 * - Using distance-based localization to reduce spurious correlations
 * - Applying adaptive weighting to ensemble members
 * - Better handling of observation error correlations
 * - Improved treatment of model error and observation bias
 *
 * @author Metada Development Team
 * @date 2024
 */

#include "LWEnKF.hpp"

#include <iostream>
#include <memory>

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

int main(int argc, char* argv[]) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: lwenkf <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Local Weighted Ensemble Kalman Filter data assimilation "
                     "application starting...";

    // Get analysis configuration
    auto analysis_config = config.GetSubsection("analysis");
    double inflation_factor = analysis_config.Get("inflation").asFloat();
    double localization_radius =
        analysis_config.Get("localization_radius").asFloat();
    std::string inflation_method =
        analysis_config.Get("inflation_method").asString();
    std::string localization_function =
        analysis_config.Get("localization_function").asString();
    std::string weighting_scheme =
        analysis_config.Get("weighting_scheme").asString();

    logger.Info() << "LWEnKF parameters:";
    logger.Info() << "  Inflation factor: " << inflation_factor;
    logger.Info() << "  Localization radius: " << localization_radius;
    logger.Info() << "  Inflation method: " << inflation_method;
    logger.Info() << "  Localization function: " << localization_function;
    logger.Info() << "  Weighting scheme: " << weighting_scheme;

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

    // Create LWEnKF object
    fwk::LWEnKF<BackendTag> lwenkf(ensemble, obs, obs_op, config);
    logger.Info() << "LWEnKF object created successfully";

    // Perform analysis
    logger.Info() << "Starting LWEnKF analysis...";
    lwenkf.Analyse();

    // Get analysis results
    auto results = lwenkf.getAnalysisResults();

    // Log analysis results
    logger.Info() << "=== LWEnKF Analysis Results ===";
    logger.Info() << "Innovation norm: " << results.innovation_norm;
    logger.Info() << "Analysis increment norm: "
                  << results.analysis_increment_norm;
    logger.Info() << "Background spread: " << results.background_spread;
    logger.Info() << "Analysis spread: " << results.analysis_spread;

    logger.Info() << "=== Localization and Weighting ===";
    logger.Info() << "Localization radius: " << results.localization_radius;
    logger.Info() << "Localization function: " << results.localization_function;
    logger.Info() << "Weighting scheme: " << results.weighting_scheme;
    logger.Info() << "Max weight: " << results.max_weight;
    logger.Info() << "Min weight: " << results.min_weight;
    logger.Info() << "Weight variance: " << results.weight_variance;

    logger.Info() << "=== Kalman Gain Statistics ===";
    double gain_range = results.max_kalman_gain - results.min_kalman_gain;
    double gain_mean =
        (results.max_kalman_gain + results.min_kalman_gain) / 2.0;
    logger.Info() << "Max Kalman gain: " << results.max_kalman_gain;
    logger.Info() << "Min Kalman gain: " << results.min_kalman_gain;
    logger.Info() << "Gain range: " << gain_range;
    logger.Info() << "Gain mean: " << gain_mean;

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

    logger.Info() << "=== Ensemble Statistics ===";
    logger.Info() << "Ensemble size: " << results.ensemble_size;
    logger.Info() << "Observation count: " << results.observation_count;
    logger.Info() << "Inflation method: " << results.inflation_method;
    logger.Info() << "Inflation factor: " << results.inflation_factor;

    // Save ensemble
    lwenkf.saveEnsemble();

    logger.Info() << "Local Weighted Ensemble Kalman Filter data assimilation "
                     "application completed successfully";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Error in LWEnKF application: " << e.what() << std::endl;
    return 1;
  }
}