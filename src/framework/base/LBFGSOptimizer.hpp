#pragma once

#include <algorithm>
#include <cmath>

#include "OptimizerBase.hpp"

namespace metada::base::optimization {

/**
 * @brief Framework-independent L-BFGS optimizer
 *
 * @details Limited-memory Broyden-Fletcher-Goldfarb-Shanno algorithm
 * for quasi-Newton optimization. Uses two-loop recursion to approximate
 * the inverse Hessian without storing the full matrix.
 *
 * This implementation is completely independent of the metada framework
 * and can be used in any C++ project.
 */
class LBFGSOptimizer : public OptimizerBase {
 public:
  /**
   * @brief Constructor with parameters
   */
  explicit LBFGSOptimizer(int max_iterations = 100, double tolerance = 1e-8,
                          double gradient_tolerance = 1e-6,
                          size_t memory_size = 10,
                          bool line_search_enabled = true)
      : max_iterations_(max_iterations),
        tolerance_(tolerance),
        gradient_tolerance_(gradient_tolerance),
        memory_size_(memory_size),
        line_search_enabled_(line_search_enabled) {}

  /**
   * @brief Minimize the cost function using L-BFGS
   */
  Result minimize(const std::vector<double>& initial_guess,
                  const CostFunction& cost_func,
                  const GradientFunction& gradient_func,
                  std::vector<double>& solution) override {
    Result result;
    solution = initial_guess;

    // L-BFGS history storage
    std::vector<std::vector<double>> s_vectors, y_vectors;
    std::vector<double> rho_values;

    // Initial evaluation
    double current_cost = cost_func(solution);
    auto gradient = gradient_func(solution);
    double gradient_norm = computeNorm(gradient);

    result.initial_cost = current_cost;
    result.cost_evaluations = 1;
    result.gradient_evaluations = 1;

    for (int iter = 0; iter < max_iterations_; ++iter) {
      // Compute search direction using L-BFGS two-loop recursion
      auto search_direction =
          computeLBFGSDirection(gradient, s_vectors, y_vectors, rho_values);

      // Line search
      double step_size =
          lineSearch(solution, search_direction, cost_func, gradient_func,
                     current_cost, line_search_enabled_);

      // Update solution
      std::vector<double> new_solution = solution;
      addScaledVector(new_solution, search_direction, step_size);

      // Evaluate new cost and gradient
      double new_cost = cost_func(new_solution);
      auto new_gradient = gradient_func(new_solution);
      result.cost_evaluations++;
      result.gradient_evaluations++;

      // Update L-BFGS vectors
      std::vector<double> s_k(solution.size());
      std::vector<double> y_k(gradient.size());

      for (size_t i = 0; i < solution.size(); ++i) {
        s_k[i] = new_solution[i] - solution[i];
      }

      for (size_t i = 0; i < gradient.size(); ++i) {
        y_k[i] = new_gradient[i] - gradient[i];
      }

      double sy_dot = dotProduct(s_k, y_k);

      if (sy_dot > 1e-12) {  // Curvature condition
        s_vectors.push_back(std::move(s_k));
        y_vectors.push_back(std::move(y_k));
        rho_values.push_back(1.0 / sy_dot);

        // Maintain memory limit
        if (s_vectors.size() > memory_size_) {
          s_vectors.erase(s_vectors.begin());
          y_vectors.erase(y_vectors.begin());
          rho_values.erase(rho_values.begin());
        }
      }

      // Update state
      double cost_change = std::abs(new_cost - current_cost);
      solution = std::move(new_solution);
      current_cost = new_cost;
      gradient = std::move(new_gradient);
      gradient_norm = computeNorm(gradient);

      result.iterations++;

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

  std::string getName() const override { return "L-BFGS"; }
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
   * @brief Get L-BFGS memory size
   */
  size_t getMemorySize() const { return memory_size_; }

  /**
   * @brief Set L-BFGS memory size
   */
  void setMemorySize(size_t memory_size) { memory_size_ = memory_size; }

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
  size_t memory_size_;
  bool line_search_enabled_;

  /**
   * @brief Compute L-BFGS search direction using two-loop recursion
   */
  std::vector<double> computeLBFGSDirection(
      const std::vector<double>& gradient,
      const std::vector<std::vector<double>>& s_vectors,
      const std::vector<std::vector<double>>& y_vectors,
      const std::vector<double>& rho_values) const {
    if (s_vectors.empty()) {
      // No history, use steepest descent
      auto direction = gradient;
      scaleVector(direction, -1.0);
      return direction;
    }

    auto q = gradient;
    std::vector<double> alpha(s_vectors.size());

    // First loop (backward)
    for (int i = static_cast<int>(s_vectors.size()) - 1; i >= 0; --i) {
      alpha[i] = rho_values[i] * dotProduct(s_vectors[i], q);
      addScaledVector(q, y_vectors[i], -alpha[i]);
    }

    // Initial Hessian approximation (identity scaled)
    double gamma = 1.0;
    if (!s_vectors.empty()) {
      const auto& s_k = s_vectors.back();
      const auto& y_k = y_vectors.back();
      gamma = dotProduct(s_k, y_k) / (dotProduct(y_k, y_k) + 1e-12);
    }

    scaleVector(q, gamma);
    auto r = q;

    // Second loop (forward)
    for (size_t i = 0; i < s_vectors.size(); ++i) {
      double beta = rho_values[i] * dotProduct(y_vectors[i], r);
      addScaledVector(r, s_vectors[i], alpha[i] - beta);
    }

    // Return negative for descent direction
    scaleVector(r, -1.0);
    return r;
  }
};

}  // namespace metada::base::optimization