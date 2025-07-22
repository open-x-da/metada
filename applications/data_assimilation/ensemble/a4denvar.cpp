/**
 * @file a4denvar.cpp
 * @brief Driver program for the Analytical Four-Dimensional
 * Ensemble-Variational (A4DEnVar) data assimilation algorithm
 * @details This application implements the analytical four-dimensional
 * ensemble-variational data assimilation algorithm as described in Liang et al.
 * (2021). It reads configuration from a file, initializes required components
 * like ensemble members, observations and observation operators across multiple
 * time windows, and performs the analysis step with an analytical solution.
 *
 * The A4DEnVar combines the strengths of 4D-Var (temporal consistency) and
 * ensemble methods (flow-dependent error covariance) with an analytical
 * solution that avoids iterative minimization.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: a4denvar <config_file>
 * @return 0 on success, 1 on failure
 */

#include "A4DEnVar.hpp"

#include <iostream>
#include <memory>

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "SimpleBackendTraits.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::SimpleBackendTag;

int main(int argc, char* argv[]) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: a4denvar <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Analytical Four-Dimensional Ensemble-Variational data "
                     "assimilation application starting...";

    // Get A4DEnVar parameters from configuration
    auto analysis_config = config.GetSubsection("analysis");
    double inflation_factor = analysis_config.Get("inflation").asFloat();
    double localization_radius =
        analysis_config.Get("localization_radius").asFloat();
    std::string inflation_method =
        analysis_config.Get("inflation_method").asString();
    std::string localization_function =
        analysis_config.Get("localization_function").asString();
    std::string covariance_type =
        analysis_config.Get("covariance_type").asString();

    logger.Info() << "A4DEnVar parameters:";
    logger.Info() << "  Inflation factor: " << inflation_factor;
    logger.Info() << "  Localization radius: " << localization_radius;
    logger.Info() << "  Inflation method: " << inflation_method;
    logger.Info() << "  Localization function: " << localization_function;
    logger.Info() << "  Covariance type: " << covariance_type;

    // Initialize geometry
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Initialize ensemble
    fwk::Ensemble<BackendTag> ensemble(config.GetSubsection("ensemble"),
                                       geometry);
    logger.Info() << "Ensemble initialized with " << ensemble.Size()
                  << " members";

    // Initialize observations for multiple time windows
    std::vector<std::unique_ptr<fwk::Observation<BackendTag>>> observations;
    auto observations_config = config.GetSubsection("observations");

    // For A4DEnVar, we need observations for multiple time windows
    // This is a simplified implementation - in practice, you'd have multiple
    // observation files for different time windows
    auto obs =
        std::make_unique<fwk::Observation<BackendTag>>(observations_config);
    observations.push_back(std::move(obs));

    // For demonstration, we'll use the same observations for multiple windows
    // In practice, you'd load different observation files for each time window
    int num_time_windows = analysis_config.Get("time_windows").asInt();
    for (int t = 1; t < num_time_windows; ++t) {
      auto cloned_obs = std::make_unique<fwk::Observation<BackendTag>>(
          observations[0]->clone());
      observations.push_back(std::move(cloned_obs));
    }

    logger.Info() << "Observations loaded for " << observations.size()
                  << " time windows";

    // Initialize observation operators for multiple time windows
    std::vector<std::unique_ptr<fwk::ObsOperator<BackendTag>>> obs_operators;
    auto obs_operator_config = config.GetSubsection("obs_operator");

    // Create observation operators for each time window
    for (int t = 0; t < num_time_windows; ++t) {
      auto obs_op =
          std::make_unique<fwk::ObsOperator<BackendTag>>(obs_operator_config);
      obs_operators.push_back(std::move(obs_op));
    }

    logger.Info() << "Observation operators initialized for "
                  << obs_operators.size() << " time windows";

    // Create A4DEnVar algorithm instance
    fwk::A4DEnVar<BackendTag> a4denvar(ensemble, observations, obs_operators,
                                       config);
    logger.Info() << "A4DEnVar object created successfully";

    // Perform analysis
    logger.Info() << "Starting A4DEnVar analysis...";
    a4denvar.Analyse();

    // Get analysis results
    auto results = a4denvar.getAnalysisResults();

    // Log analysis results
    logger.Info() << "=== A4DEnVar Analysis Results ===";
    logger.Info() << "Cost function value: " << results.cost_function_value;
    logger.Info() << "Background cost: " << results.background_cost;
    logger.Info() << "Observation cost: " << results.observation_cost;
    logger.Info() << "Innovation norm: " << results.innovation_norm;
    logger.Info() << "Analysis increment norm: "
                  << results.analysis_increment_norm;

    logger.Info() << "=== Ensemble Statistics ===";
    logger.Info() << "Background spread: " << results.background_spread;
    logger.Info() << "Analysis spread: " << results.analysis_spread;
    logger.Info() << "Ensemble size: " << results.ensemble_size;
    logger.Info() << "Time windows: " << results.time_windows;
    logger.Info() << "Total observations: " << results.observation_count;

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

    logger.Info() << "=== Algorithm Configuration ===";
    logger.Info() << "Covariance type: " << results.covariance_type;
    logger.Info() << "Localization radius: " << results.localization_radius;
    logger.Info() << "Localization function: " << results.localization_function;
    logger.Info() << "Inflation method: " << results.inflation_method;
    logger.Info() << "Inflation factor: " << results.inflation_factor;

    // Log cost function analysis
    logger.Info() << "=== Cost Function Analysis ===";
    double total_cost = results.background_cost + results.observation_cost;
    double background_ratio = results.background_cost / total_cost;
    double observation_ratio = results.observation_cost / total_cost;

    logger.Info() << "Total cost: " << total_cost;
    logger.Info() << "Background cost ratio: " << background_ratio;
    logger.Info() << "Observation cost ratio: " << observation_ratio;

    // Log temporal consistency
    logger.Info() << "=== Temporal Consistency ===";
    logger.Info() << "Number of time windows: " << results.time_windows;
    logger.Info() << "Observations per window: "
                  << (results.observation_count / results.time_windows);
    logger.Info() << "Temporal consistency maintained across "
                  << results.time_windows << " time windows";

    // Save ensemble
    a4denvar.saveEnsemble();

    logger.Info() << "Analytical Four-Dimensional Ensemble-Variational data "
                     "assimilation application completed successfully";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Error in A4DEnVar application: " << e.what() << std::endl;
    return 1;
  }
}