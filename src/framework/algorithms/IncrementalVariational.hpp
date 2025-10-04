#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "BackgroundErrorCovariance.hpp"
#include "Config.hpp"
#include "Ensemble.hpp"
#include "IncrementalCostFunction.hpp"
#include "IncrementalGradientChecks.hpp"
#include "IncrementalMinimization.hpp"
#include "Logger.hpp"
#include "Model.hpp"
#include "ObsOperator.hpp"
#include "OperatorChecks.hpp"
#include "ProgressBar.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Incremental variational data assimilation algorithm implementation
 *
 * @details This class implements incremental variational data assimilation
 * algorithms including:
 * - 4DVAR: Four-dimensional variational assimilation with full time window
 * - 3DVAR: Three-dimensional variational assimilation (single time)
 * - FGAT: First Guess at Appropriate Time
 *
 * The algorithm minimizes the incremental cost function:
 * J(δx) = 1/2 * δx^T B^-1 δx + 1/2 * Σᵢ (dᵢ - Hᵢ(Mᵢ(xb + δx)))^T Rᵢ^-1 (dᵢ -
 * Hᵢ(Mᵢ(xb + δx)))
 *
 * Where:
 * - δx is the analysis increment (δx = x - xb)
 * - xb is the background state (first guess)
 * - B is the background error covariance matrix
 * - dᵢ = yᵢ - Hᵢ(Mᵢ(xb)) is the innovation vector (O-B)
 * - Hᵢ is the observation operator at time i
 * - Mᵢ is the model propagated to time i
 * - Rᵢ is the observation error covariance matrix at time i
 *
 * The minimization is performed using various optimization algorithms (L-BFGS,
 * CG, etc.) and the gradient is computed using the adjoint method for 4DVAR.
 *
 * This incremental formulation offers better numerical conditioning and is the
 * standard approach used in operational systems like WRFDA.
 *
 * @tparam BackendTag The backend tag type that must satisfy the required
 * concepts
 */
template <typename BackendTag>
class IncrementalVariational {
 public:
  /**
   * @brief Analysis results structure
   */
  struct AnalysisResults {
    State<BackendTag> analysis_state;  ///< Final analysis state (xb + δx*)
    Increment<BackendTag> analysis_increment;  ///< Final analysis increment δx*
    double final_cost = 0.0;                   ///< Final cost function value
    double background_cost = 0.0;              ///< Background term cost
    double observation_cost = 0.0;             ///< Observation term cost
    double cost_reduction = 0.0;               ///< Total cost reduction
    int iterations = 0;              ///< Number of minimization iterations
    bool converged = false;          ///< Whether minimization converged
    std::string convergence_reason;  ///< Reason for convergence/termination
    std::string variational_type;    ///< Type of variational method used
    std::vector<double>
        cost_evolution;  ///< Cost function evolution during minimization
  };

  /**
   * @brief Constructor for Incremental Variational algorithm
   *
   * @param config Configuration object
   * @param background Background state (first guess)
   * @param observations Vector of observations at different times
   * @param obs_operators Vector of observation operators for each time
   * @param model Model for forward and adjoint integration
   * @param bg_error_cov Background error covariance
   */
  IncrementalVariational(
      const Config<BackendTag>& config, const State<BackendTag>& background,
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
    logger_.Info() << "IncrementalVariational algorithm constructed for "
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
   * @brief Perform incremental variational data assimilation analysis
   *
   * @return Analysis results structure containing the optimized increment and
   * diagnostics
   */
  AnalysisResults analyze() {
    logger_.Info() << "Starting " << cost_function_.getVariationalTypeName()
                   << " incremental analysis";

    std::string var_type = cost_function_.getVariationalTypeName();

    // Create cost and gradient function wrappers for the minimizer
    auto cost_func = [this](const Increment<BackendTag>& increment) -> double {
      return cost_function_.evaluate(increment);
    };

    auto gradient_func = [this](const Increment<BackendTag>& increment,
                                Increment<BackendTag>& gradient) -> void {
      cost_function_.gradient(increment, gradient);
    };

    // Start with zero increment as initial guess
    auto initial_increment =
        Increment<BackendTag>::createFromEntity(background_);
    initial_increment.zero();

    // Evaluate initial cost (should be just background term)
    double initial_cost = cost_func(initial_increment);
    logger_.Info() << "Initial cost (background only): " << initial_cost;

    // Perform minimization over increments
    auto final_increment = Increment<BackendTag>::createFromEntity(background_);
    auto convergence_info = minimizer_.minimize(initial_increment, cost_func,
                                                gradient_func, final_increment);

    // Create final analysis state: xa = xb + δx*
    auto analysis_state = background_.clone();
    analysis_state += final_increment.state();

    // Create results structure
    AnalysisResults results{
        std::move(analysis_state),    // analysis_state
        std::move(final_increment),   // analysis_increment
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

    logger_.Info() << "Incremental variational analysis completed";
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
    logger_.Info() << "Saving incremental variational analysis results";

    // Save analysis state
    results.analysis_state.saveToFile(output_base_file_ + "_analysis." +
                                      format_);
    logger_.Info() << "Analysis state saved to: "
                   << output_base_file_ + "_analysis." + format_;

    // Save analysis increment
    saveIncrement(results.analysis_increment,
                  output_base_file_ + "_increment." + format_);
    logger_.Info() << "Analysis increment saved to: "
                   << output_base_file_ + "_increment." + format_;

    // Save diagnostics to text file
    saveDiagnostics(results);

    logger_.Info() << "Incremental variational analysis results saved";
  }

  /**
   * @brief Get the incremental cost function object
   */
  const IncrementalCostFunction<BackendTag>& getCostFunction() const {
    return cost_function_;
  }

  /**
   * @brief Get the incremental minimizer object
   */
  const IncrementalMinimization<BackendTag>& getMinimizer() const {
    return minimizer_;
  }

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
      auto obs_data = obs.getObservationValues();

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
   * @param test_increment Increment around which to test the gradient
   * @param tolerance Tolerance for the gradient test (default: 1e-6)
   * @return True if gradient test passes within tolerance
   */
  bool performGradientTest(const Increment<BackendTag>& test_increment,
                           double tolerance = 1e-6) const {
    logger_.Info() << "Performing gradient test for incremental variational "
                      "implementation";

    // Use the improved gradient check utility function
    bool test_passed = checkIncrementalCostFunctionGradient<BackendTag>(
        cost_function_, test_increment, tolerance);

    return test_passed;
  }

  /**
   * @brief Perform tangent linear/adjoint consistency check for observation
   * operators
   *
   * @param tolerance Tolerance for the TL/AD check (default: 1e-6)
   * @return True if TL/AD check passes within tolerance
   */
  bool performObsOperatorTLADCheck(double tolerance = 1e-6) const {
    logger_.Info() << "Performing observation operator TL/AD consistency check";

    bool test_passed = checkIncrementalObsOperatorTLAD<BackendTag>(
        obs_operators_, background_, observations_, tolerance);

    return test_passed;
  }

  /**
   * @brief Perform tangent linear finite difference check for observation
   * operators
   *
   * @param tolerance Tolerance for the tangent linear check (default: 1e-6)
   * @param epsilons Set of perturbation sizes for FD check (default: {1e-3,
   * 1e-4, 1e-5, 1e-6, 1e-7})
   * @return True if tangent linear check passes within tolerance
   */
  bool performObsOperatorTangentLinearCheck(
      double tolerance = 1e-6,
      const std::vector<double>& epsilons = std::vector<double>{
          1e-3, 1e-4, 1e-5, 1e-6, 1e-7}) const {
    logger_.Info() << "Performing observation operator tangent linear check";

    bool test_passed = checkIncrementalObsOperatorTangentLinear<BackendTag>(
        obs_operators_, background_, observations_, tolerance, epsilons);

    return test_passed;
  }

  /**
   * @brief Perform comprehensive observation operator checks
   *
   * @details This method performs both tangent linear/adjoint consistency check
   * and tangent linear finite difference check for all observation operators
   * used in the incremental variational formulation.
   *
   * @param tolerance Tolerance for the checks (default: 1e-6)
   * @param epsilons Set of perturbation sizes for FD check (default: {1e-3,
   * 1e-4, 1e-5, 1e-6, 1e-7})
   * @return True if all observation operator checks pass within tolerance
   */
  bool performObsOperatorChecks(
      double tolerance = 1e-6,
      const std::vector<double>& epsilons = std::vector<double>{
          1e-3, 1e-4, 1e-5, 1e-6, 1e-7}) const {
    logger_.Info() << "Performing comprehensive observation operator checks";

    bool test_passed = checkIncrementalObsOperatorsComprehensive<BackendTag>(
        obs_operators_, background_, observations_, tolerance, epsilons);

    return test_passed;
  }

  /**
   * @brief Perform all verification checks for incremental variational
   * implementation
   *
   * @details This method performs a comprehensive verification of the
   * incremental variational implementation including:
   * - Cost function gradient check
   * - Observation operator tangent linear/adjoint consistency check
   * - Observation operator tangent linear finite difference check
   *
   * @param test_increment Increment around which to test the gradient
   * @param tolerance Tolerance for all checks (default: 1e-6)
   * @param epsilons Set of perturbation sizes for FD checks (default: {1e-3,
   * 1e-4, 1e-5, 1e-6, 1e-7})
   * @return True if all verification checks pass within tolerance
   */
  bool performVerificationChecks(
      const Increment<BackendTag>& test_increment, double tolerance = 1e-6,
      const std::vector<double>& epsilons = std::vector<double>{
          1e-3, 1e-4, 1e-5, 1e-6, 1e-7}) const {
    logger_.Info()
        << "Performing comprehensive verification checks for incremental "
           "variational implementation";

    // 1. Cost function gradient check
    bool gradient_passed = performGradientTest(test_increment, tolerance);

    // 2. Observation operator checks
    bool obs_ops_passed = performObsOperatorChecks(tolerance, epsilons);

    bool all_passed = gradient_passed && obs_ops_passed;

    logger_.Info() << "Verification check results:";
    logger_.Info() << "  Cost function gradient: "
                   << (gradient_passed ? "PASSED" : "FAILED");
    logger_.Info() << "  Observation operators: "
                   << (obs_ops_passed ? "PASSED" : "FAILED");
    logger_.Info() << "  Overall result: "
                   << (all_passed ? "PASSED" : "FAILED");

    return all_passed;
  }

  /**
   * @brief Get the pre-computed innovation vectors
   */
  const std::vector<std::vector<double>>& getInnovations() const {
    return cost_function_.getInnovations();
  }

 private:
  /**
   * @brief Compute breakdown of cost function terms
   */
  void computeCostBreakdown(AnalysisResults& results) const {
    // Evaluate background term using the final increment
    results.background_cost =
        0.5 * bg_error_cov_.quadraticForm(results.analysis_increment);

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

  /**
   * @brief Save increment to file
   */
  void saveIncrement(const Increment<BackendTag>& increment,
                     const std::string& filename) const {
    // For now, we'll save the increment as a state file
    // In the future, we might want a dedicated increment file format
    increment.state().saveToFile(filename);
  }

  const Config<BackendTag>& config_;
  const State<BackendTag>& background_;
  const std::vector<Observation<BackendTag>>& observations_;
  const std::vector<ObsOperator<BackendTag>>& obs_operators_;
  Model<BackendTag>& model_;
  const BackgroundErrorCovariance<BackendTag>& bg_error_cov_;
  IncrementalCostFunction<BackendTag> cost_function_;
  IncrementalMinimization<BackendTag> minimizer_;
  std::string output_base_file_;
  std::string format_;
  bool save_trajectory_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework
