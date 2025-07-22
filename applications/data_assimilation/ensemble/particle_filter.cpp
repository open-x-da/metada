/**
 * @file particle_filter.cpp
 * @brief Driver program for the Particle Filter (PF) data assimilation
 * algorithm
 * @details This application implements the particle filter data assimilation
 *          algorithm. It reads configuration from a file, initializes required
 *          components like ensemble members, observations and observation
 *          operators, and performs the analysis step.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *            Expected format: particle_filter <config_file>
 * @return 0 on success, 1 on failure
 */

#include "ApplicationContext.hpp"
#include "Config.hpp"
#include "Ensemble.hpp"
#include "Geometry.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "ParticleFilter.hpp"
// #include "SimpleBackendTraits.hpp"
#include "MACOMBackendTraits.hpp"

namespace fwk = metada::framework;
// using BackendTag = metada::traits::SimpleBackendTag;
using BackendTag = metada::traits::MACOMBackendTag;

int main(int argc, char** argv) {
  try {
    // Validate command line arguments
    if (argc != 2) {
      std::cerr << "Usage: particle_filter <config_file>" << std::endl;
      return 1;
    }

    // Initialize application context
    auto context = fwk::ApplicationContext<BackendTag>(argc, argv);
    auto& logger = context.getLogger();
    auto& config = context.getConfig();

    logger.Info()
        << "Particle Filter data assimilation application starting...";

    // Get particle filter parameters from configuration
    auto analysis_config = config.GetSubsection("analysis");
    int resampling_threshold =
        analysis_config.Get("resampling_threshold").asInt();
    std::string resampling_method =
        analysis_config.Get("resampling_method").asString();
    bool jittering_enabled = analysis_config.Get("jittering_enabled").asBool();
    double jittering_std = analysis_config.Get("jittering_std").asFloat();

    logger.Info() << "Particle Filter parameters:";
    logger.Info() << "  Resampling threshold: " << resampling_threshold;
    logger.Info() << "  Resampling method: " << resampling_method;
    logger.Info() << "  Jittering enabled: "
                  << (jittering_enabled ? "Yes" : "No");
    if (jittering_enabled) {
      logger.Info() << "  Jittering std: " << jittering_std;
    }

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

    // Initialize observation operator
    fwk::ObsOperator<BackendTag> obs_op(config.GetSubsection("obs_operator"));
    logger.Info() << "Observation operator initialized";

    // Create particle filter algorithm instance
    fwk::ParticleFilter<BackendTag> particle_filter(ensemble, obs, obs_op,
                                                    config);

    logger.Info() << "Starting Particle Filter analysis...";

    // Perform particle filter analysis
    particle_filter.Analyse();

    // Get analysis results
    auto results = particle_filter.getAnalysisResults();

    // Log results
    logger.Info() << "=== Particle Filter Analysis Results ===";
    logger.Info() << "Effective sample size: " << results.effective_sample_size;
    logger.Info() << "Max weight: " << results.max_weight;
    logger.Info() << "Min weight: " << results.min_weight;
    logger.Info() << "Weight variance: " << results.weight_variance;
    logger.Info() << "Resampling performed: "
                  << (results.resampling_performed ? "Yes" : "No");
    logger.Info() << "Resampling method: " << results.resampling_method;
    logger.Info() << "Resampling threshold: " << results.resampling_threshold;

    // Log weight statistics
    logger.Info() << "=== Weight Statistics ===";
    const auto& weights = results.weights;
    double mean_weight = 1.0 / weights.size();
    double weight_std = 0.0;
    for (const auto& weight : weights) {
      double diff = weight - mean_weight;
      weight_std += diff * diff;
    }
    weight_std = std::sqrt(weight_std / weights.size());

    logger.Info() << "Mean weight: " << mean_weight;
    logger.Info() << "Weight std: " << weight_std;
    logger.Info() << "Coefficient of variation: " << (weight_std / mean_weight);

    // Log likelihood statistics
    logger.Info() << "=== Likelihood Statistics ===";
    const auto& likelihoods = results.likelihood_values;
    double max_likelihood =
        *std::max_element(likelihoods.begin(), likelihoods.end());
    double min_likelihood =
        *std::min_element(likelihoods.begin(), likelihoods.end());
    double mean_likelihood =
        std::accumulate(likelihoods.begin(), likelihoods.end(), 0.0) /
        likelihoods.size();

    logger.Info() << "Max likelihood: " << max_likelihood;
    logger.Info() << "Min likelihood: " << min_likelihood;
    logger.Info() << "Mean likelihood: " << mean_likelihood;

    // Save analysis results
    particle_filter.saveEnsemble();

    // Optional: Save detailed analysis results
    if (analysis_config.Get("save_detailed_results").asBool()) {
      logger.Info() << "Saving detailed analysis results...";
      // Implementation would save detailed results to file
    }

    logger.Info() << "Particle Filter data assimilation application completed "
                     "successfully";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Particle Filter application failed: " << e.what()
              << std::endl;
    return 1;
  }
}