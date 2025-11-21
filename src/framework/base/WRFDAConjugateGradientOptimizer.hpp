#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "OptimizerBase.hpp"

namespace metada::base::optimization {

/**
 * @brief WRFDA-aligned Conjugate Gradient optimizer
 *
 * @details This optimizer implements the exact workflow from WRFDA's
 * `da_minimise_cg.inc`, matching the incremental variational inner-loop
 * minimization algorithm. Key features:
 *
 * - Computes gradient along search direction (not at new solution)
 * - Always uses cost approximation during iterations: J ≈ J0 + 0.5 * <ghat0,
 * xhat>
 * - Supports gradient reorthonormalization (modified Gram-Schmidt)
 * - Uses WRFDA's step size formula: step = rrmold / apdotp
 *
 * Cost function evaluation strategy:
 * - Initial iteration: Compute exact J0 for baseline
 * - During minimization: Use approximation (never compute exact J during
 * iterations)
 * - Final iteration: Compute exact J for final result
 *
 * This matches WRFDA's workflow:
 * 1. Initialize: Compute cost and gradient at xhat=0
 * 2. Iterate: Compute gradient along search direction, update,
 * reorthonormalize, approximate cost
 * 3. Finalize: Compute exact final cost and gradient
 */
class WRFDAConjugateGradientOptimizer : public OptimizerBase {
 public:
  /**
   * @brief Constructor with parameters
   *
   * @param max_iterations Maximum number of iterations
   * @param tolerance Cost function tolerance
   * @param gradient_tolerance Gradient norm tolerance
   * @param orthonorm_gradient If true, enable gradient reorthonormalization
   * (matches WRFDA's orthonorm_gradient=true)
   *
   * @note Cost approximation is always enabled during minimization iterations.
   * Exact cost is only computed at initialization and finalization.
   */
  explicit WRFDAConjugateGradientOptimizer(int max_iterations = 100,
                                           double tolerance = 1e-8,
                                           double gradient_tolerance = 1e-6,
                                           bool orthonorm_gradient = true)
      : max_iterations_(max_iterations),
        tolerance_(tolerance),
        gradient_tolerance_(gradient_tolerance),
        orthonorm_gradient_(orthonorm_gradient) {}

  /**
   * @brief Minimize the cost function using WRFDA-aligned CG
   *
   * @details This implements the exact workflow from WRFDA's da_minimise_cg:
   *
   * [1.0] Initialization:
   *   - Compute initial cost J(xhat=0)
   *   - Compute initial gradient ghat = ∇J(xhat=0)
   *   - Set initial search direction phat = -ghat
   *
   * [2.0] CG Iteration Loop:
   *   For each iteration:
   *   1. Compute gradient along search direction: fhat = ∇J(phat)
   *   2. Compute step size: step = rrmold / apdotp where apdotp = <fhat, phat>
   *   3. Update: ghat = ghat + step * fhat, xhat = xhat + step * phat
   *   4. Reorthonormalize gradient (if enabled)
   *   5. Update search direction: phat = -ghat + β*phat
   *   6. Approximate cost: J ≈ J0 + 0.5 * <ghat0, xhat> (never compute exact J)
   *
   * [3.0] Finalization:
   *   - Compute exact final cost
   *   - Compute exact final gradient
   */
  Result minimize(const std::vector<double>& initial_guess,
                  const CostFunction& cost_func,
                  const GradientFunction& gradient_func,
                  std::vector<double>& solution) override {
    Result result;
    solution = initial_guess;

    //-------------------------------------------------------------------------
    // [1.0] Initialization
    //-------------------------------------------------------------------------

    // Compute initial cost
    double j0_total = cost_func(solution);
    result.cost_evaluations = 1;

    // Compute initial gradient
    std::vector<double> ghat = gradient_func(solution);
    result.gradient_evaluations = 1;

    // Store initial gradient for cost approximation
    std::vector<double> ghat0 = ghat;

    // Compute initial gradient norm squared
    double rrmold = dotProduct(ghat, ghat);
    double j_grad_norm_target = std::sqrt(rrmold);

    if (rrmold == 0.0) {
      // Zero gradient, already at minimum
      result.converged = true;
      result.iterations = 0;
      result.initial_cost = j0_total;
      result.final_cost = j0_total;
      result.gradient_norm = 0.0;
      result.cost_reduction = 0.0;
      result.convergence_reason = "Zero initial gradient";
      return result;
    }

    // Initial search direction: phat = -ghat
    std::vector<double> phat = ghat;
    scaleVector(phat, -1.0);

    // Orthonormalization setup
    std::vector<std::vector<double>> qhat;
    if (orthonorm_gradient_) {
      qhat.reserve(max_iterations_ + 1);
      std::vector<double> q0 = ghat;
      scaleVector(q0, 1.0 / j_grad_norm_target);
      qhat.push_back(std::move(q0));
    }

    result.initial_cost = j0_total;
    double current_cost = j0_total;

    // Log initial iteration
    if (iteration_logger_) {
      iteration_logger_(0, j0_total, j0_total, j_grad_norm_target, 0.0);
    }

    //-------------------------------------------------------------------------
    // [2.0] CG Iteration Loop
    //-------------------------------------------------------------------------

    for (int iter = 1; iter <= max_iterations_; ++iter) {
      if (rrmold == 0.0) {
        result.converged = true;
        result.convergence_reason = "Zero gradient norm";
        break;
      }

      // Step 1: Compute gradient along search direction
      // In WRFDA: fhat = ∇J(phat) where phat is treated as a control variable
      // This computes the gradient at the point phat (the search direction)
      // Note: In WRFDA's incremental variational, this computes the gradient
      // of the cost function evaluated at the search direction vector
      std::vector<double> fhat = gradient_func(phat);
      result.gradient_evaluations++;

      // Step 2: Compute step size using WRFDA formula
      // apdotp = <∇J(phat), phat> (curvature along search direction)
      double apdotp = dotProduct(fhat, phat);
      double step = 0.0;
      if (apdotp > 0.0) {
        step = rrmold / apdotp;
      }

      // Step 3: Update gradient and solution
      // ghat = ghat + step * fhat
      // xhat = xhat + step * phat
      addScaledVector(ghat, fhat, step);
      addScaledVector(solution, phat, step);

      // Step 4: Orthonormalize gradient (modified Gram-Schmidt)
      if (orthonorm_gradient_) {
        for (int i = static_cast<int>(qhat.size()) - 1; i >= 0; --i) {
          double gdot = dotProduct(ghat, qhat[i]);
          addScaledVector(ghat, qhat[i], -gdot);
        }
      }

      // Step 5: Compute new gradient norm squared
      double rrmnew = dotProduct(ghat, ghat);
      double rrmnew_norm = std::sqrt(rrmnew);

      // Step 6: Compute Polak-Ribière beta
      double ratio = 0.0;
      if (rrmold > 0.0) {
        ratio = rrmnew / rrmold;
      }

      // Store normalized gradient for reorthonormalization
      if (orthonorm_gradient_) {
        std::vector<double> q_new = ghat;
        scaleVector(q_new, 1.0 / rrmnew_norm);
        qhat.push_back(std::move(q_new));
      }

      // Step 7: Update search direction: phat = -ghat + β*phat
      scaleVector(phat, ratio);
      std::vector<double> neg_ghat = ghat;
      scaleVector(neg_ghat, -1.0);
      addScaledVector(phat, neg_ghat, 1.0);

      // Step 8: Compute cost using approximation
      // Cost approximation: J ≈ J0 + 0.5 * <ghat0, xhat>
      // We never compute exact J during minimization iterations to avoid
      // expensive forward operator evaluations. Exact J is only computed at
      // initialization (for J0) and finalization (for final result).
      current_cost = j0_total + 0.5 * dotProduct(ghat0, solution);

      // Update for next iteration
      rrmold = rrmnew;

      result.iterations = iter;

      // Log iteration
      if (iteration_logger_) {
        iteration_logger_(iter, current_cost, current_cost, rrmnew_norm, step);
      }

      // Step 9: Check convergence
      // Convergence: ||ghat|| < eps * ||ghat0||
      // For now, use gradient_tolerance_ as eps
      double eps = gradient_tolerance_;
      if (rrmnew_norm < eps * j_grad_norm_target) {
        result.converged = true;
        result.convergence_reason = "Gradient tolerance reached";
        break;
      }

      if (iter >= max_iterations_) {
        result.converged = false;
        result.convergence_reason = "Maximum iterations reached";
        break;
      }
    }

    //-------------------------------------------------------------------------
    // [3.0] Finalization
    //-------------------------------------------------------------------------

    // Compute exact final cost
    double final_cost = cost_func(solution);
    result.cost_evaluations++;

    // Compute exact final gradient
    std::vector<double> final_gradient = gradient_func(solution);
    result.gradient_evaluations++;

    double final_gradient_norm =
        std::sqrt(dotProduct(final_gradient, final_gradient));

    result.final_cost = final_cost;
    result.gradient_norm = final_gradient_norm;
    result.cost_reduction = result.initial_cost - result.final_cost;

    return result;
  }

  std::string getName() const override { return "WRFDA Conjugate Gradient"; }

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
   * @brief Check if gradient reorthonormalization is enabled
   */
  bool isOrthonormGradientEnabled() const { return orthonorm_gradient_; }

  /**
   * @brief Enable/disable gradient reorthonormalization
   */
  void setOrthonormGradientEnabled(bool enabled) {
    orthonorm_gradient_ = enabled;
  }

  void setIterationLogger(
      std::function<void(int, double, double, double, double)> logger)
      override {
    iteration_logger_ = std::move(logger);
  }

 private:
  int max_iterations_;
  double tolerance_;
  double gradient_tolerance_;
  bool orthonorm_gradient_;
  std::function<void(int, double, double, double, double)> iteration_logger_;
};

}  // namespace metada::base::optimization
