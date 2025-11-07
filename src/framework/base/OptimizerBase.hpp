#pragma once

#include <functional>
#include <string>
#include <vector>

namespace metada::base::optimization {

/**
 * @brief Framework-independent base class for optimization algorithms
 *
 * @details This abstract base class defines a pure mathematical interface
 * for optimization algorithms. It works with standard library types and
 * has no dependencies on the metada framework adapters, making it
 * completely reusable in any project.
 */
class OptimizerBase {
 public:
  /**
   * @brief Cost function type definition
   *
   * Function that takes a parameter vector and returns the cost function value
   */
  using CostFunction = std::function<double(const std::vector<double>&)>;

  /**
   * @brief Gradient function type definition
   *
   * Function that takes a parameter vector and returns the gradient vector
   */
  using GradientFunction =
      std::function<std::vector<double>(const std::vector<double>&)>;

  /**
   * @brief Optimization result structure
   */
  struct Result {
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

  /** @brief Virtual destructor */
  virtual ~OptimizerBase() = default;

  /**
   * @brief Minimize the cost function starting from initial guess
   *
   * @param initial_guess Initial parameter vector
   * @param cost_func Cost function to minimize
   * @param gradient_func Gradient function
   * @param solution Output parameter vector containing the optimized solution
   * @return Optimization result information
   */
  virtual Result minimize(const std::vector<double>& initial_guess,
                          const CostFunction& cost_func,
                          const GradientFunction& gradient_func,
                          std::vector<double>& solution) = 0;

  /**
   * @brief Get the optimizer name
   */
  virtual std::string getName() const = 0;

  /**
   * @brief Get maximum iterations
   */
  virtual int getMaxIterations() const = 0;

  /**
   * @brief Set maximum iterations
   */
  virtual void setMaxIterations(int max_iterations) = 0;

  /**
   * @brief Get cost function tolerance
   */
  virtual double getTolerance() const = 0;

  /**
   * @brief Set cost function tolerance
   */
  virtual void setTolerance(double tolerance) = 0;

  /**
   * @brief Get gradient tolerance
   */
  virtual double getGradientTolerance() const = 0;

  /**
   * @brief Set gradient tolerance
   */
  virtual void setGradientTolerance(double gradient_tolerance) = 0;

  /**
   * @brief Set iteration logger callback
   *
   * @details This callback is invoked at each iteration with iteration
   * information. The signature is: (iteration, previous_cost, current_cost,
   * gradient_norm, step_size). Optimizers that support iteration logging
   * should override this method.
   *
   * @param logger Function to call at each iteration. If null, iteration
   * logging is disabled.
   */
  virtual void setIterationLogger(
      std::function<void(int, double, double, double, double)> logger) {
    // Default implementation: no-op (optimizers can override)
    (void)logger;
  }

 protected:
  /**
   * @brief Check convergence criteria
   */
  bool checkConvergence(double cost_change, double gradient_norm, int iteration,
                        int max_iterations, double tolerance,
                        double gradient_tolerance, std::string& reason) const {
    if (std::abs(cost_change) < tolerance) {
      reason = "Cost function tolerance reached";
      return true;
    }

    if (gradient_norm < gradient_tolerance) {
      reason = "Gradient tolerance reached";
      return true;
    }

    if (iteration >= max_iterations) {
      reason = "Maximum iterations reached";
      return true;
    }

    return false;
  }

  /**
   * @brief Compute vector norm
   */
  double computeNorm(const std::vector<double>& vec) const {
    double norm = 0.0;
    for (double val : vec) {
      norm += val * val;
    }
    return std::sqrt(norm);
  }

  /**
   * @brief Compute dot product of two vectors
   */
  double dotProduct(const std::vector<double>& a,
                    const std::vector<double>& b) const {
    if (a.size() != b.size()) {
      throw std::invalid_argument("Vector sizes must match for dot product");
    }

    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
      result += a[i] * b[i];
    }
    return result;
  }

  /**
   * @brief Scale vector by scalar
   */
  void scaleVector(std::vector<double>& vec, double scale) const {
    for (double& val : vec) {
      val *= scale;
    }
  }

  /**
   * @brief Add scaled vector: a += scale * b
   */
  void addScaledVector(std::vector<double>& a, const std::vector<double>& b,
                       double scale) const {
    if (a.size() != b.size()) {
      throw std::invalid_argument("Vector sizes must match for addition");
    }

    for (size_t i = 0; i < a.size(); ++i) {
      a[i] += scale * b[i];
    }
  }

  /**
   * @brief Simple line search with backtracking
   */
  double lineSearch(const std::vector<double>& x,
                    const std::vector<double>& direction,
                    const CostFunction& cost_func,
                    const GradientFunction& gradient_func, double current_cost,
                    bool enabled = true) const {
    if (!enabled) {
      return 1.0;  // No line search, use unit step
    }

    const double c1 = 1e-4;  // Armijo condition parameter
    const double rho = 0.5;  // Backtracking factor
    const int max_line_search_iter = 20;

    // Compute directional derivative
    auto gradient = gradient_func(x);
    double directional_deriv = dotProduct(gradient, direction);

    double step_size = 1.0;

    for (int i = 0; i < max_line_search_iter; ++i) {
      // Try current step size
      std::vector<double> test_x = x;
      addScaledVector(test_x, direction, step_size);

      double test_cost = cost_func(test_x);

      // Armijo condition
      if (test_cost <= current_cost + c1 * step_size * directional_deriv) {
        return step_size;
      }

      // Reduce step size
      step_size *= rho;
    }

    return step_size;  // Return last tried step size
  }
};

}  // namespace metada::base::optimization