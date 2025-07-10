/**
 * @file variational.cpp
 * @brief Driver program for Variational Data Assimilation (4DVAR/3DVAR/FGAT)
 * @details This application implements variational data assimilation algorithms
 *          including 4DVAR, 3DVAR, and FGAT. It reads configuration from a
 * file, initializes required components like background state, observations,
 *          observation operators, model, and background error covariance, and
 *          performs the variational analysis.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: variational <config_file>
 * @return 0 on success, 1 on failure
 */

#include "Variational.hpp"

#include "ApplicationContext.hpp"
#include "BackgroundErrorCovariance.hpp"
#include "Config.hpp"
#include "Geometry.hpp"
#include "Model.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"
#include "WRFBackendTraits.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::WRFBackendTag;

int main(int argc, char** argv) {
  // Initialize application context
  auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
  auto& logger = context.getLogger();
  auto& config = context.getConfig();

  logger.Info() << "Variational data assimilation application starting...";

  try {
    // Validate command line arguments
    if (argc != 2) {
      logger.Error() << "Usage: variational <config_file>";
      return 1;
    }

    // Get variational type from configuration
    std::string var_type = config.Get("variational_type").asString();
    logger.Info() << "Variational method: " << var_type;

    // Initialize geometry
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Initialize background state
    fwk::State<BackendTag> background(config.GetSubsection("background"),
                                      geometry);

    // Initialize model
    fwk::Model<BackendTag> model(config.GetSubsection("model"));

    // Initialize observations
    std::vector<fwk::Observation<BackendTag>> observations;
    auto obs_configs = config.GetSubsectionsFromVector("observations");
    observations.reserve(obs_configs.size());

    for (const auto& obs_config : obs_configs) {
      observations.emplace_back(obs_config);
    }

    if (observations.empty()) {
      logger.Error() << "No observations specified in configuration";
      return 1;
    }

    logger.Info() << "Loaded " << observations.size() << " observation time(s)";

    // Initialize observation operators
    std::vector<fwk::ObsOperator<BackendTag>> obs_operators;
    auto obs_op_configs = config.GetSubsectionsFromVector("obs_operators");
    obs_operators.reserve(obs_op_configs.size());

    for (const auto& obs_op_config : obs_op_configs) {
      obs_operators.emplace_back(obs_op_config);
    }

    // Validate observation and operator counts match
    if (observations.size() != obs_operators.size()) {
      logger.Error() << "Number of observations (" << observations.size()
                     << ") must match number of observation operators ("
                     << obs_operators.size() << ")";
      return 1;
    }

    // Initialize background error covariance
    fwk::BackgroundErrorCovariance<BackendTag> bg_error_cov(
        config.GetSubsection("background_covariance"));

    // Create variational algorithm instance
    fwk::Variational<BackendTag> variational(
        config, background, observations, obs_operators, model, bg_error_cov);

    logger.Info() << "Starting " << var_type << " analysis...";

    // Perform variational analysis
    auto results = variational.analyze();

    // Log results
    logger.Info() << "=== Analysis Results ===";
    logger.Info() << "Variational type: " << results.variational_type;
    logger.Info() << "Final cost: " << results.final_cost;
    logger.Info() << "Cost reduction: " << results.cost_reduction;
    logger.Info() << "Background cost: " << results.background_cost;
    logger.Info() << "Observation cost: " << results.observation_cost;
    logger.Info() << "Iterations: " << results.iterations;
    logger.Info() << "Converged: " << (results.converged ? "Yes" : "No");
    logger.Info() << "Convergence reason: " << results.convergence_reason;

    // Save analysis results
    variational.saveAnalysis(results);

    // Compute and display innovation statistics
    auto innovation_stats =
        variational.computeInnovationStatistics(results.analysis_state);

    logger.Info() << "=== Innovation Statistics ===";
    for (size_t i = 0; i < innovation_stats.size(); ++i) {
      logger.Info() << "Time " << i << ": mean=" << innovation_stats[i].first
                    << ", rms=" << innovation_stats[i].second;
    }

    // Optional: Perform gradient test if requested
    if (config.Get("perform_gradient_test").asBool()) {
      logger.Info() << "Performing gradient test...";
      bool gradient_test_passed = variational.performGradientTest(background);
      logger.Info() << "Gradient test: "
                    << (gradient_test_passed ? "PASSED" : "FAILED");
    }

    logger.Info()
        << "Variational data assimilation application completed successfully";
    return 0;

  } catch (const std::exception& e) {
    logger.Error() << "Variational application failed: " << e.what();
    return 1;
  }
}