#pragma once

#include <cmath>
#include <functional>

#include "OptimizerBase.hpp"

namespace metada::base::optimization {

/**
 * @brief Framework-independent Steepest Descent optimizer
 *
 * @details Simple gradient descent optimization algorithm.
 * Uses the negative gradient as the search direction at each iteration.
 *
 * This implementation is completely independent of the metada framework
 * and can be used in any C++ project.
 */
class SteepestDescentOptimizer : public OptimizerBase {
 public:
  /**
   * @brief Constructor with parameters
   */
  explicit SteepestDescentOptimizer(int max_iterations = 100,
                                    double tolerance = 1e-8,
                                    double gradient_tolerance = 1e-6,
                                    bool line_search_enabled = true)
      : max_iterations_(max_iterations),
        tolerance_(tolerance),
        gradient_tolerance_(gradient_tolerance),
        line_search_enabled_(line_search_enabled) {}

  /**
   * @brief Minimize the cost function using Steepest Descent
   */
  Result minimize(const std::vector<double>& initial_guess,
                  const CostFunction& cost_func,
                  const GradientFunction& gradient_func,
                  std::vector<double>& solution) override {
    Result result;
    solution = initial_guess;

    // Initial evaluation
    double current_cost = cost_func(solution);
    result.initial_cost = current_cost;
    result.cost_evaluations = 1;
    result.gradient_evaluations = 0;

    for (int iter = 0; iter < max_iterations_; ++iter) {
      double previous_cost = current_cost;

      // Compute gradient
      auto gradient = gradient_func(solution);
      double gradient_norm = computeNorm(gradient);
      result.gradient_evaluations++;

      // Search direction is negative gradient
      auto search_direction = gradient;
      scaleVector(search_direction, -1.0);

      // Line search
      double step_size =
          lineSearch(solution, search_direction, cost_func, gradient_func,
                     current_cost, line_search_enabled_);

      if (step_size <= 0.0) {
        result.convergence_reason = "Line search failed";
        break;
      }

      // Update solution
      std::vector<double> new_solution = solution;
      addScaledVector(new_solution, search_direction, step_size);

      double new_cost = cost_func(new_solution);
      result.cost_evaluations++;

      double cost_change = std::abs(new_cost - current_cost);

      solution = std::move(new_solution);
      current_cost = new_cost;

      result.iterations++;

      // Log iteration information if logger is set
      if (iteration_logger_) {
        iteration_logger_(result.iterations, previous_cost, current_cost,
                          gradient_norm, step_size);
      }

      // Check convergence
      if (checkConvergence(cost_change, gradient_norm, result.iterations,
                           max_iterations_, tolerance_, gradient_tolerance_,
                           result.convergence_reason)) {
        result.converged = (result.iterations < max_iterations_);
        result.gradient_norm = gradient_norm;
        break;
      }
    }

    result.final_cost = current_cost;
    result.cost_reduction = result.initial_cost - result.final_cost;

    return result;
  }

  std::string getName() const override { return "Steepest Descent"; }
  int getMaxIterations() const override { return max_iterations_; }
  void setMaxIterations(int max_iterations) override {
    max_iterations_ = max_iterations;
  }
  double getTolerance() const override { return tolerance_; }
  void setTolerance(double tolerance) override { tolerance_ = tolerance; }
  double getGradientTolerance() const override { return gradient_tolerance_; }
  void setGradientTolerance(double gradient_tolerance) override {
    gradient_tolerance_ = gradient_tolerance;
  }

  /**
   * @brief Check if line search is enabled
   */
  bool isLineSearchEnabled() const { return line_search_enabled_; }

  /**
   * @brief Enable/disable line search
   */
  void setLineSearchEnabled(bool enabled) { line_search_enabled_ = enabled; }

  /**
   * @brief Set iteration logger callback
   */
  void setIterationLogger(
      std::function<void(int, double, double, double, double)> logger)
      override {
    iteration_logger_ = std::move(logger);
  }

 private:
  int max_iterations_;
  double tolerance_;
  double gradient_tolerance_;
  bool line_search_enabled_;
  std::function<void(int, double, double, double, double)> iteration_logger_;
};

}  // namespace metada::base::optimization