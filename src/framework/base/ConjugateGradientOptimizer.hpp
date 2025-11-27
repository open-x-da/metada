#pragma once

#include <algorithm>
#include <cmath>
#include <functional>

#include "OptimizerBase.hpp"

namespace metada::base::optimization {

/**
 * @brief Framework-independent Conjugate Gradient optimizer
 *
 * @details Conjugate Gradient method for optimization. Uses the Polak-Ribière
 * formula for computing the conjugate directions with automatic restart.
 *
 * This implementation is completely independent of the metada framework
 * and can be used in any C++ project.
 */
class ConjugateGradientOptimizer : public OptimizerBase {
 public:
  /**
   * @brief Constructor with parameters
   */
  explicit ConjugateGradientOptimizer(int max_iterations = 100,
                                      double tolerance = 1e-8,
                                      double gradient_tolerance = 1e-6,
                                      bool line_search_enabled = true,
                                      int restart_interval = 0)
      : max_iterations_(max_iterations),
        tolerance_(tolerance),
        gradient_tolerance_(gradient_tolerance),
        line_search_enabled_(line_search_enabled),
        restart_interval_(restart_interval) {}

  /**
   * @brief Minimize the cost function using Conjugate Gradient
   */
  Result minimize(const std::vector<double>& initial_guess,
                  const CostFunction& cost_func,
                  const GradientFunction& gradient_func,
                  std::vector<double>& solution) override {
    Result result;
    solution = initial_guess;

    // Initial evaluation
    double current_cost = cost_func(solution);
    auto gradient = gradient_func(solution);
    auto prev_gradient = gradient;

    // Initial search direction is negative gradient
    auto search_direction = gradient;
    scaleVector(search_direction, -1.0);

    double gradient_norm = computeNorm(gradient);
    result.initial_cost = current_cost;
    result.cost_evaluations = 1;
    result.gradient_evaluations = 1;

    for (int iter = 0; iter < max_iterations_; ++iter) {
      double previous_cost = current_cost;
      // Line search along search direction
      double step_size =
          lineSearch(solution, search_direction, cost_func, gradient_func,
                     current_cost, line_search_enabled_);

      // Update solution
      std::vector<double> new_solution = solution;
      addScaledVector(new_solution, search_direction, step_size);

      // Evaluate new cost and gradient
      double new_cost = cost_func(new_solution);
      prev_gradient = gradient;
      gradient = gradient_func(new_solution);
      result.cost_evaluations++;
      result.gradient_evaluations++;

      // Compute beta for conjugate direction (Polak-Ribière formula)
      std::vector<double> gradient_diff = gradient;
      addScaledVector(gradient_diff, prev_gradient, -1.0);

      double beta = dotProduct(gradient, gradient_diff) /
                    (dotProduct(prev_gradient, prev_gradient) + 1e-12);
      beta = std::max(0.0, beta);  // Polak-Ribière+ modification

      // Check for restart (every restart_interval iterations or if beta is too
      // small)
      bool restart = false;
      if (restart_interval_ > 0 && (iter + 1) % restart_interval_ == 0) {
        restart = true;
      }
      if (beta < 1e-12) {
        restart = true;
      }

      // Update search direction
      if (restart) {
        // Restart with steepest descent
        search_direction = gradient;
        scaleVector(search_direction, -1.0);
      } else {
        // Conjugate direction update
        scaleVector(search_direction, beta);
        addScaledVector(search_direction, gradient, -1.0);
      }

      // Update state and check convergence
      double cost_change = std::abs(new_cost - current_cost);
      solution = std::move(new_solution);
      current_cost = new_cost;
      gradient_norm = computeNorm(gradient);

      result.iterations++;

      if (iteration_logger_) {
        iteration_logger_(result.iterations, previous_cost, current_cost,
                          gradient_norm, step_size);
      }

      // Check convergence
      if (checkConvergence(cost_change, gradient_norm, result.iterations,
                           max_iterations_, tolerance_, gradient_tolerance_,
                           result.convergence_reason)) {
        result.converged = (result.iterations < max_iterations_);
        break;
      }
    }

    result.final_cost = current_cost;
    result.gradient_norm = gradient_norm;
    result.cost_reduction = result.initial_cost - result.final_cost;

    return result;
  }

  std::string getName() const override { return "Conjugate Gradient"; }
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
   * @brief Get restart interval
   */
  int getRestartInterval() const { return restart_interval_; }

  /**
   * @brief Set restart interval (0 = no automatic restart)
   */
  void setRestartInterval(int interval) { restart_interval_ = interval; }

  void setIterationLogger(
      std::function<void(int, double, double, double, double)> logger)
      override {
    iteration_logger_ = std::move(logger);
  }

  /**
   * @brief Check if line search is enabled
   */
  bool isLineSearchEnabled() const { return line_search_enabled_; }

  /**
   * @brief Enable/disable line search
   */
  void setLineSearchEnabled(bool enabled) { line_search_enabled_ = enabled; }

 private:
  int max_iterations_;
  double tolerance_;
  double gradient_tolerance_;
  bool line_search_enabled_;
  int restart_interval_;
  std::function<void(int, double, double, double, double)> iteration_logger_;
};

}  // namespace metada::base::optimization