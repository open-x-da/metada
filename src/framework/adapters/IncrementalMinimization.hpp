#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "BackendTraits.hpp"
#include "ConfigConcepts.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "NonCopyable.hpp"
#include "State.hpp"

// Framework-independent optimizers
#include "ConjugateGradientOptimizer.hpp"
#include "LBFGSOptimizer.hpp"
#include "OptimizerBase.hpp"
#include "SteepestDescentOptimizer.hpp"

namespace metada::framework {

// Forward declarations
template <typename BackendTag>
  requires ConfigBackendType<BackendTag>
class Config;

/**
 * @brief Framework adapter for incremental minimization algorithms
 *
 * @details This class bridges between the framework types
 * (Increment<BackendTag>) and the pure mathematical optimization algorithms in
 * base::optimization. It converts framework types to std::vector<double>, calls
 * the appropriate optimizer, and converts results back to framework types.
 *
 * The key difference from the regular Minimization class is that this works
 * with increments δx instead of total states x, implementing the incremental
 * variational formulation.
 *
 * @tparam BackendTag The backend tag type
 */
template <typename BackendTag>
  requires StateBackendType<BackendTag>
class IncrementalMinimization : public NonCopyable {
 public:
  /**
   * @brief Cost function type definition for incremental formulation
   *
   * Function that takes an increment and returns the cost function value
   */
  using CostFunction = std::function<double(const Increment<BackendTag>&)>;

  /**
   * @brief Gradient function type definition for incremental formulation
   *
   * Function that takes an increment and computes the gradient into an
   * increment
   */
  using GradientFunction =
      std::function<void(const Increment<BackendTag>&, Increment<BackendTag>&)>;

  /**
   * @brief Convergence information structure
   */
  struct ConvergenceInfo {
    bool converged = false;          ///< Whether convergence was achieved
    int iterations = 0;              ///< Number of iterations performed
    double initial_cost = 0.0;       ///< Initial cost function value
    double final_cost = 0.0;         ///< Final cost function value
    double gradient_norm = 0.0;      ///< Final gradient norm
    double cost_reduction = 0.0;     ///< Total cost reduction achieved
    std::string convergence_reason;  ///< Reason for convergence/termination

    // Performance metrics
    int cost_evaluations = 0;  ///< Total number of cost function evaluations
    int gradient_evaluations = 0;  ///< Total number of gradient evaluations
  };

  /** @brief Default constructor is deleted */
  IncrementalMinimization() = delete;

  /**
   * @brief Constructor with configuration
   *
   * @param config Configuration object containing minimization parameters
   */
  explicit IncrementalMinimization(const Config<BackendTag>& config)
      : optimizer_(createOptimizer(config)) {
    logger_.Info() << "IncrementalMinimization adapter constructed with "
                   << optimizer_->getName() << " algorithm";
    logger_.Info() << "Max iterations: " << optimizer_->getMaxIterations();
    logger_.Info() << "Tolerance: " << optimizer_->getTolerance();
    logger_.Info() << "Gradient tolerance: "
                   << optimizer_->getGradientTolerance();
  }

  /**
   * @brief Move constructor
   */
  IncrementalMinimization(IncrementalMinimization&& other) noexcept = default;

  /**
   * @brief Move assignment operator
   */
  IncrementalMinimization& operator=(IncrementalMinimization&& other) noexcept =
      default;

  /**
   * @brief Destructor
   */
  ~IncrementalMinimization() = default;

  /**
   * @brief Minimize the cost function over increments δx
   *
   * @param initial_increment Initial guess for the increment δx₀
   * @param cost_func Cost function J(δx)
   * @param gradient_func Gradient function ∇J(δx)
   * @param final_increment Output: final optimized increment δx*
   * @return Convergence information
   */
  ConvergenceInfo minimize(const Increment<BackendTag>& initial_increment,
                           const CostFunction& cost_func,
                           const GradientFunction& gradient_func,
                           Increment<BackendTag>& final_increment) const {
    logger_.Info() << "Starting incremental minimization with "
                   << optimizer_->getName();

    // Convert initial increment to std::vector<double>
    auto initial_data =
        initial_increment.template getData<std::vector<double>>();
    std::vector<double> final_data(initial_data.size());

    // Create cost and gradient function wrappers for the optimizer
    auto cost_wrapper = [&cost_func, &initial_increment](
                            const std::vector<double>& x) -> double {
      // Create a temporary increment from the vector data
      auto temp_increment = Increment<BackendTag>(initial_increment.state());
      auto temp_data = temp_increment.template getData<std::vector<double>>();
      std::copy(x.begin(), x.end(), temp_data.begin());
      return cost_func(temp_increment);
    };

    auto gradient_wrapper =
        [&gradient_func, &initial_increment](
            const std::vector<double>& x) -> std::vector<double> {
      // Create a temporary increment from the vector data
      auto temp_increment = Increment<BackendTag>(initial_increment.state());
      auto temp_data = temp_increment.template getData<std::vector<double>>();
      std::copy(x.begin(), x.end(), temp_data.begin());

      // Create a temporary gradient increment
      auto temp_gradient = Increment<BackendTag>(initial_increment.state());
      temp_gradient.zero();

      // Compute gradient
      gradient_func(temp_increment, temp_gradient);

      // Extract gradient data
      auto gradient_data =
          temp_gradient.template getData<std::vector<double>>();
      return gradient_data;
    };

    // Perform minimization
    auto convergence_info = optimizer_->minimize(initial_data, cost_wrapper,
                                                 gradient_wrapper, final_data);

    // Convert final result back to increment
    final_increment = Increment<BackendTag>(initial_increment.state());
    auto final_increment_data =
        final_increment.template getData<std::vector<double>>();
    std::copy(final_data.begin(), final_data.end(),
              final_increment_data.begin());

    // Convert convergence info to our format
    ConvergenceInfo result;
    result.converged = convergence_info.converged;
    result.iterations = convergence_info.iterations;
    result.initial_cost = convergence_info.initial_cost;
    result.final_cost = convergence_info.final_cost;
    result.gradient_norm = convergence_info.gradient_norm;
    result.cost_reduction =
        convergence_info.initial_cost - convergence_info.final_cost;
    result.convergence_reason = convergence_info.convergence_reason;
    result.cost_evaluations = convergence_info.cost_evaluations;
    result.gradient_evaluations = convergence_info.gradient_evaluations;

    logger_.Info() << "Incremental minimization completed";
    logger_.Info() << "Converged: " << (result.converged ? "Yes" : "No");
    logger_.Info() << "Iterations: " << result.iterations;
    logger_.Info() << "Cost reduction: " << result.cost_reduction;
    logger_.Info() << "Final gradient norm: " << result.gradient_norm;

    return result;
  }

  /**
   * @brief Get the optimizer name
   */
  std::string getAlgorithmName() const { return optimizer_->getName(); }

  /**
   * @brief Get the maximum number of iterations
   */
  int getMaxIterations() const { return optimizer_->getMaxIterations(); }

  /**
   * @brief Get the cost function tolerance
   */
  double getTolerance() const { return optimizer_->getTolerance(); }

  /**
   * @brief Get the gradient tolerance
   */
  double getGradientTolerance() const {
    return optimizer_->getGradientTolerance();
  }

 private:
  /**
   * @brief Create the appropriate optimizer based on configuration
   */
  std::unique_ptr<base::optimization::OptimizerBase> createOptimizer(
      const Config<BackendTag>& config) {
    std::string algorithm = config.Get("minimization_algorithm").asString();
    int max_iterations = config.Get("max_iterations").asInt();
    double tolerance = config.Get("cost_tolerance").asFloat();
    double gradient_tolerance = config.Get("gradient_tolerance").asFloat();

    if (algorithm == "L-BFGS") {
      return std::make_unique<base::optimization::LBFGSOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    } else if (algorithm == "CG") {
      return std::make_unique<base::optimization::ConjugateGradientOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    } else if (algorithm == "SteepestDescent") {
      return std::make_unique<base::optimization::SteepestDescentOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    } else {
      logger_.Warning() << "Unknown minimization algorithm: " << algorithm
                        << ". Defaulting to L-BFGS.";
      return std::make_unique<base::optimization::LBFGSOptimizer>(
          max_iterations, tolerance, gradient_tolerance);
    }
  }

  std::unique_ptr<base::optimization::OptimizerBase> optimizer_;
  Logger<BackendTag>& logger_ = Logger<BackendTag>::Instance();
};

}  // namespace metada::framework
