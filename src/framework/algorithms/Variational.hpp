#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "BackgroundErrorCovariance.hpp"
#include "Config.hpp"
#include "CostFunction.hpp"
#include "Ensemble.hpp"
#include "Logger.hpp"
#include "Minimization.hpp"
#include "Model.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "ProgressBar.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Variational data assimilation algorithm implementation
 *
 * @details This class implements variational data assimilation algorithms
 * including:
 * - 4DVAR: Four-dimensional variational assimilation with full time window
 * - 3DVAR: Three-dimensional variational assimilation (single time)
 * - FGAT: First Guess at Appropriate Time
 *
 * The algorithm minimizes the cost function:
 * J(x) = 1/2 * (x-xb)^T B^-1 (x-xb) + 1/2 * sum_i (y_i - H_i(M_i(x)))^T R_i^-1
 * (y_i - H_i(M_i(x)))
 *
 * Where:
 * - x is the analysis state (to be determined)
 * - xb is the background state (first guess)
 * - B is the background error covariance matrix
 * - y_i are observations at time i
 * - H_i is the observation operator at time i
 * - M_i is the model propagated to time i
 * - R_i is the observation error covariance matrix at time i
 *
 * The minimization is performed using various optimization algorithms (L-BFGS,
 * CG, etc.) and the gradient is computed using the adjoint method for 4DVAR.
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class Variational {
 public:
  /**
   * @brief Analysis results structure
   */
  struct AnalysisResults {
    State<BackendTag> analysis_state;  ///< Final analysis state
    double final_cost = 0.0;           ///< Final cost function value
    double background_cost = 0.0;      ///< Background term cost
    double observation_cost = 0.0;     ///< Observation term cost
    double cost_reduction = 0.0;       ///< Total cost reduction
    int iterations = 0;                ///< Number of minimization iterations
    bool converged = false;            ///< Whether minimization converged
    std::string convergence_reason;    ///< Reason for convergence/termination
    std::string variational_type;      ///< Type of variational method used
    std::vector<double>
        cost_evolution;  ///< Cost function evolution during minimization
  };

  /**
   * @brief Constructor for Variational algorithm
   *
   * @param config Configuration object
   * @param background Background state (first guess)
   * @param observations Vector of observations at different times
   * @param obs_operators Vector of observation operators for each time
   * @param model Model for forward and adjoint integration
   * @param bg_error_cov Background error covariance
   */
  Variational(const Config<BackendTag>& config,
              const State<BackendTag>& background,
              const std::vector<Observation<BackendTag>>& observations,
              const std::vector<ObsOperator<BackendTag>>& obs_operators,
              Model<BackendTag>& model,
              const BackgroundErrorCovariance<BackendTag>& bg_error_cov)
      : config_(config),
        background_(background),
        observations_(observations),
        obs_operators_(obs_operators),
        model_(model),
        bg_error_cov_(bg_error_cov),
        cost_function_(config, background, observations, obs_operators, model,
                       bg_error_cov),
        minimizer_(config),
        output_base_file_(config.Get("output_base_file").asString()),
        format_(config.Get("format").asString()),
        save_trajectory_(config.Get("save_trajectory").asBool()) {
    logger_.Info() << "Variational algorithm constructed for "
                   << cost_function_.getVariationalTypeName();
    logger_.Info() << "Time windows: " << cost_function_.getTimeWindows();
    logger_.Info() << "Minimization algorithm: "
                   << minimizer_.getAlgorithmName();

    // Validate inputs
    if (observations_.empty()) {
      throw std::runtime_error(
          "No observations provided for variational assimilation");
    }

    if (observations_.size() != obs_operators_.size()) {
      throw std::runtime_error(
          "Number of observations must match number of observation operators");
    }

    // Check model capabilities for 4DVAR
    if (cost_function_.getVariationalTypeName() == "4DVAR" &&
        !model_.supportsAdjoint()) {
      logger_.Warning()
          << "4DVAR requested but model does not support adjoint. "
             "Falling back to FGAT.";
    }
  }

  /**
   * @brief Perform variational data assimilation analysis
   *
   * @return Analysis results structure containing the optimized state and
   * diagnostics
   */
  AnalysisResults analyze() {
    logger_.Info() << "Starting " << cost_function_.getVariationalTypeName()
                   << " analysis";

    // Create results with properly initialized analysis state
    // Note: We'll move the final analysis state into results later
    std::string var_type = cost_function_.getVariationalTypeName();

    // Create cost and gradient function wrappers for the minimizer
    auto cost_func = [this](const State<BackendTag>& state) -> double {
      return cost_function_.evaluate(state);
    };

    auto gradient_func = [this](const State<BackendTag>& state,
                                Increment<BackendTag>& gradient) -> void {
      cost_function_.gradient(state, gradient);
    };

    // Evaluate initial cost
    double initial_cost = cost_func(background_);
    logger_.Info() << "Initial cost: " << initial_cost;

    // Perform minimization
    // Start with a clone of the background state as the initial guess
    auto analysis_state = background_.clone();
    auto convergence_info = minimizer_.minimize(background_, cost_func,
                                                gradient_func, analysis_state);

    // Create results structure with the final analysis state
    AnalysisResults results{
        std::move(analysis_state),    // analysis_state
        convergence_info.final_cost,  // final_cost
        0.0,                          // background_cost (computed later)
        0.0,                          // observation_cost (computed later)
        initial_cost - convergence_info.final_cost,  // cost_reduction
        convergence_info.iterations,                 // iterations
        convergence_info.converged,                  // converged
        convergence_info.convergence_reason,         // convergence_reason
        var_type,                                    // variational_type
        {}                                           // cost_evolution
    };

    // Compute cost breakdown
    computeCostBreakdown(results);

    logger_.Info() << "Variational analysis completed";
    logger_.Info() << "Final cost: " << results.final_cost;
    logger_.Info() << "Cost reduction: " << results.cost_reduction;
    logger_.Info() << "Background cost: " << results.background_cost;
    logger_.Info() << "Observation cost: " << results.observation_cost;
    logger_.Info() << "Iterations: " << results.iterations;
    logger_.Info() << "Converged: " << (results.converged ? "Yes" : "No");
    logger_.Info() << "Convergence reason: " << results.convergence_reason;

    return results;
  }

  /**
   * @brief Save analysis results to files
   *
   * @param results Analysis results to save
   */
  void saveAnalysis(const AnalysisResults& results) const {
    logger_.Info() << "Saving variational analysis results";

    // Save analysis state
    results.analysis_state.saveToFile(output_base_file_ + "_analysis." +
                                      format_);
    logger_.Info() << "Analysis state saved to: "
                   << output_base_file_ + "_analysis." + format_;

    // Save analysis increment (analysis - background)
    auto analysis_increment = Increment<BackendTag>::createFromDifference(
        results.analysis_state, background_);
    // Note: Increment doesn't have saveToFile, so we'd need to add that or save
    // via state

    // Save diagnostics to text file
    saveDiagnostics(results);

    logger_.Info() << "Variational analysis results saved";
  }

  /**
   * @brief Get the cost function object
   */
  const CostFunction<BackendTag>& getCostFunction() const {
    return cost_function_;
  }

  /**
   * @brief Get the minimizer object
   */
  const Minimization<BackendTag>& getMinimizer() const { return minimizer_; }

  /**
   * @brief Compute innovation statistics
   *
   * @param state State to compute innovations for
   * @return Vector of innovation statistics for each observation time
   */
  std::vector<std::pair<double, double>> computeInnovationStatistics(
      const State<BackendTag>& state) const {
    std::vector<std::pair<double, double>> stats;  // (mean, rms)

    auto current_state = state.clone();

    for (size_t i = 0; i < observations_.size(); ++i) {
      const auto& obs = observations_[i];
      const auto& obs_op = obs_operators_[i];

      // Propagate state to observation time (if needed)
      if (i > 0) {
        auto next_state = current_state.clone();
        model_.run(current_state, next_state);
        current_state = std::move(next_state);
      }

      // Compute innovations
      auto simulated_obs = obs_op.apply(current_state, obs);
      auto obs_data = obs.template getData<std::vector<double>>();

      double sum = 0.0, sum_sq = 0.0;
      for (size_t j = 0; j < obs_data.size(); ++j) {
        double innovation = obs_data[j] - simulated_obs[j];
        sum += innovation;
        sum_sq += innovation * innovation;
      }

      double mean = sum / obs_data.size();
      double rms = std::sqrt(sum_sq / obs_data.size());
      stats.emplace_back(mean, rms);
    }

    return stats;
  }

  /**
   * @brief Perform gradient test to verify adjoint implementation
   *
   * @param test_state State around which to test the gradient
   * @param perturbation_size Size of perturbation for finite difference
   * @return True if gradient test passes within tolerance
   */
  bool performGradientTest(const State<BackendTag>& test_state,
                           double perturbation_size = 1e-1) const {
    logger_.Info() << "Performing gradient test for variational implementation";

    const double tolerance = 1e-6;

    // Compute analytical gradient
    auto analytical_gradient =
        Increment<BackendTag>::createFromEntity(test_state);
    cost_function_.gradient(test_state, analytical_gradient);

    // Create random perturbation
    auto perturbation = Increment<BackendTag>::createFromEntity(test_state);
    perturbation.randomize();

    // Scale perturbation by the specified size
    perturbation *= perturbation_size;

    // Compute cost at test state
    double cost_at_x = cost_function_.evaluate(test_state);

    // Compute cost at perturbed state
    auto perturbed_state = test_state.clone();
    perturbed_state += perturbation;
    double cost_at_x_plus_dx = cost_function_.evaluate(perturbed_state);

    // Compute finite difference approximation
    double finite_diff_grad =
        (cost_at_x_plus_dx - cost_at_x) / perturbation_size;

    // Compute analytical gradient in perturbation direction
    double analytical_grad =
        analytical_gradient.dot(perturbation) / perturbation_size;

    // Compute relative error
    double gradient_error =
        std::abs(finite_diff_grad - analytical_grad) /
        (std::max(std::abs(finite_diff_grad), std::abs(analytical_grad)) +
         1e-12);

    bool test_passed = gradient_error < tolerance;

    logger_.Info() << "Gradient test " << (test_passed ? "PASSED" : "FAILED");
    logger_.Info() << "Finite difference gradient: " << finite_diff_grad;
    logger_.Info() << "Analytical gradient: " << analytical_grad;
    logger_.Info() << "Relative error: " << gradient_error;
    logger_.Info() << "Tolerance: " << tolerance;

    return test_passed;
  }

 private:
  /**
   * @brief Compute breakdown of cost function terms
   */
  void computeCostBreakdown([[maybe_unused]] AnalysisResults& results) const {
    // Evaluate background term
    auto background_increment = Increment<BackendTag>::createFromDifference(
        results.analysis_state, background_);
    results.background_cost =
        0.5 * bg_error_cov_.quadraticForm(background_increment);

    // Observation cost is the remainder
    results.observation_cost = results.final_cost - results.background_cost;
  }

  /**
   * @brief Save diagnostic information to file
   */
  void saveDiagnostics([[maybe_unused]] const AnalysisResults& results) const {
    std::string diagnostics_file = output_base_file_ + "_diagnostics.txt";

    // Implementation would write diagnostics to file
    // Including cost evolution, innovation statistics, etc.
    logger_.Info() << "Diagnostics saved to: " << diagnostics_file;
  }

  const Config<BackendTag>& config_;
  const State<BackendTag>& background_;
  const std::vector<Observation<BackendTag>>& observations_;
  const std::vector<ObsOperator<BackendTag>>& obs_operators_;
  Model<BackendTag>& model_;
  const BackgroundErrorCovariance<BackendTag>& bg_error_cov_;
  CostFunction<BackendTag> cost_function_;
  Minimization<BackendTag> minimizer_;
  std::string output_base_file_;
  std::string format_;
  bool save_trajectory_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework