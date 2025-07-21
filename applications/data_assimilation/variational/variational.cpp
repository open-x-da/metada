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
#include "SimpleBackendTraits.hpp"
#include "State.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::SimpleBackendTag;

int main(int argc, char** argv) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: variational <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info() << "Variational data assimilation application starting...";

    // Get variational type from configuration
    std::string var_type = config.Get("variational_type").asString();
    logger.Info() << "Variational method: " << var_type;

    // Initialize geometry
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Initialize background state
    fwk::State<BackendTag> background(config.GetSubsection("background"),
                                      geometry);
    logger.Info() << "Background state: " << background;

    // Initialize model
    fwk::Model<BackendTag> model(config.GetSubsection("model"));

    // Initialize observations
    // For 3DVAR, we have a single observation time with multiple types
    std::vector<fwk::Observation<BackendTag>> observations;
    observations.emplace_back(config.GetSubsection("observations"));

    logger.Info() << "Loaded observations for " << var_type << " analysis";

    // Initialize observation operators
    // For 3DVAR, we need a single observation operator matching the observation
    // structure
    std::vector<fwk::ObsOperator<BackendTag>> obs_operators;
    obs_operators.emplace_back(config.GetSubsection("obs_operator"));

    logger.Info() << "Loaded observation operator for " << var_type
                  << " analysis";

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
      bool gradient_test_passed = variational.performGradientTest(background);
      logger.Info() << "Gradient test "
                    << (gradient_test_passed ? "PASSED" : "FAILED");
    }

    logger.Info()
        << "Variational data assimilation application completed successfully";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Variational application failed: " << e.what() << std::endl;
    return 1;
  }
}