/**
 * @file incremental_variational.cpp
 * @brief Driver program for Incremental Variational Data Assimilation
 * (4DVAR/3DVAR/FGAT)
 * @details This application implements incremental variational data
 * assimilation algorithms including 4DVAR, 3DVAR, and FGAT using the
 * incremental formulation where the cost function works with analysis
 * increments δx instead of total states x. This offers better numerical
 * conditioning and is the standard approach used in operational systems like
 * WRFDA.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: incremental_variational <config_file>
 * @return 0 on success, 1 on failure
 */

#include "ApplicationContext.hpp"
#include "BackgroundErrorCovariance.hpp"
#include "Config.hpp"
#include "Geometry.hpp"
#include "IncrementalVariational.hpp"
#include "Model.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"
#include "WRFBackendTraits.hpp"

namespace fwk = metada::framework;
using BackendTag = metada::traits::WRFBackendTag;

int main(int argc, char** argv) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: incremental_variational <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info()
        << "Incremental Variational data assimilation application starting...";

    // Get variational type from configuration
    std::string var_type = config.Get("variational_type").asString();
    logger.Info() << "Variational method: " << var_type
                  << " (incremental formulation)";

    // Initialize observations first (data-driven workflow)
    // For 3DVAR, we have a single observation time with multiple types
    std::vector<fwk::Observation<BackendTag>> observations;
    observations.emplace_back(config.GetSubsection("observation"));

    logger.Info() << "Loaded observations for " << var_type << " analysis";
    logger.Info() << "Observations: " << observations.front();

    // Initialize geometry
    fwk::Geometry<BackendTag> geometry(config.GetSubsection("geometry"));

    // Initialize background state
    fwk::State<BackendTag> background(config.GetSubsection("background"),
                                      geometry);
    logger.Info() << "Background state: " << background;

    // Initialize model
    fwk::Model<BackendTag> model(config.GetSubsection("model"));

    // Initialize observation operators
    // For 3DVAR, we need a single observation operator matching the observation
    // structure
    std::vector<fwk::ObsOperator<BackendTag>> obs_operators;
    obs_operators.emplace_back(config.GetSubsection("obs_operator"),
                               observations.front());

    logger.Info() << "Loaded observation operator for " << var_type
                  << " analysis";

    // Initialize background error covariance
    fwk::BackgroundErrorCovariance<BackendTag> bg_error_cov(
        config.GetSubsection("background_covariance"));

    // Create incremental variational algorithm instance
    fwk::IncrementalVariational<BackendTag> incremental_variational(
        config, background, observations, obs_operators, model, bg_error_cov);

    logger.Info() << "Starting " << var_type << " incremental analysis...";

    // Perform incremental variational analysis
    auto results = incremental_variational.analyze();

    // Log results
    logger.Info() << "=== Incremental Analysis Results ===";
    logger.Info() << "Variational type: " << results.variational_type;
    logger.Info() << "Final cost: " << results.final_cost;
    logger.Info() << "Cost reduction: " << results.cost_reduction;
    logger.Info() << "Background cost: " << results.background_cost;
    logger.Info() << "Observation cost: " << results.observation_cost;
    logger.Info() << "Iterations: " << results.iterations;
    logger.Info() << "Converged: " << (results.converged ? "Yes" : "No");
    logger.Info() << "Convergence reason: " << results.convergence_reason;

    // Save analysis results
    incremental_variational.saveAnalysis(results);

    // Compute and display innovation statistics
    auto innovation_stats = incremental_variational.computeInnovationStatistics(
        results.analysis_state);

    logger.Info() << "=== Innovation Statistics ===";
    for (size_t i = 0; i < innovation_stats.size(); ++i) {
      logger.Info() << "Time " << i << ": mean=" << innovation_stats[i].first
                    << ", rms=" << innovation_stats[i].second;
    }

    // Display pre-computed innovation vectors
    const auto& innovations = incremental_variational.getInnovations();
    logger.Info() << "=== Pre-computed Innovation Vectors ===";
    for (size_t i = 0; i < innovations.size(); ++i) {
      logger.Info() << "Time window " << i << ": " << innovations[i].size()
                    << " innovations";
      if (!innovations[i].empty()) {
        double mean_innovation = 0.0;
        for (double val : innovations[i]) {
          mean_innovation += val;
        }
        mean_innovation /= innovations[i].size();
        logger.Info() << "  Mean innovation: " << mean_innovation;
      }
    }

    // Optional: Perform gradient test if requested
    if (config.Get("perform_gradient_test").asBool()) {
      // Create a small test increment
      auto test_increment =
          fwk::Increment<BackendTag>::createFromEntity(background);
      test_increment.randomize();
      test_increment *= 0.01;  // Small perturbation

      bool gradient_test_passed =
          incremental_variational.performGradientTest(test_increment);
      logger.Info() << "Incremental gradient test "
                    << (gradient_test_passed ? "PASSED" : "FAILED");
    }

    logger.Info() << "Incremental Variational data assimilation application "
                     "completed successfully";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Incremental Variational application failed: " << e.what()
              << std::endl;
    return 1;
  }
}
