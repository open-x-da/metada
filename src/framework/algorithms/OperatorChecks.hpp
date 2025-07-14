#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "CostFunction.hpp"
#include "Increment.hpp"
#include "Logger.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Check the tangent linear and adjoint consistency for an ObsOperator
 * using the inner product test.
 *
 * This function verifies the correctness of the tangent linear (TL) and adjoint
 * (AD) implementations of an observation operator H by checking the following
 * mathematical property:
 *
 * \f[
 *   \langle H dx, dy \rangle = \langle dx, H^T dy \rangle
 * \f]
 *
 * where:
 *   - \f$ dx \f$ is a random state increment
 *   - \f$ dy \f$ is a random observation increment
 *   - \f$ H dx \f$ is the result of applying the tangent linear operator
 *   - \f$ H^T dy \f$ is the result of applying the adjoint operator
 *
 * The function computes the relative error between the two inner products and
 * returns true if the error is less than the specified tolerance.
 *
 * @tparam BackendTag The backend type tag
 * @param obs_op The observation operator
 * @param state The state vector
 * @param obs The observation vector
 * @param tol The tolerance for the relative error (default: 1e-6)
 * @return true if the check passes, false otherwise
 */
template <typename BackendTag>
bool checkObsOperatorTLAD(const ObsOperator<BackendTag>& obs_op,
                          const State<BackendTag>& state,
                          const Observation<BackendTag>& obs,
                          double tol = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  // 1. Create random state increment dx and obs increment dy
  auto dx = Increment<BackendTag>::createFromEntity(state);
  dx.randomize();

  std::vector<double> dy(obs.size());
  for (auto& v : dy) v = (double(rand()) / RAND_MAX - 0.5);

  // 2. Apply tangent linear: H dx
  auto Hdx = obs_op.applyTangentLinear(dx, state, obs);

  // 3. Apply adjoint: H^T dy
  auto HTdy = obs_op.applyAdjoint(dy, state, obs);

  // 4. Compute inner products
  double a = std::inner_product(Hdx.begin(), Hdx.end(), dy.begin(), 0.0);
  double b = dx.dot(HTdy);

  double rel_error =
      std::abs(a - b) / (std::max(std::abs(a), std::abs(b)) + 1e-12);

  logger.Info() << "TL/AD check: <Hdx, dy> = " << std::setprecision(13)
                << std::scientific << a << ", <dx, H^T dy> = " << b
                << ", rel error = " << rel_error;

  return rel_error < tol;
}

/**
 * @brief Check the tangent linear implementation of an ObsOperator using Taylor
 * expansion and finite differences.
 *
 * This function verifies the correctness of the tangent linear (TL)
 * implementation of an observation operator H by comparing the TL result to a
 * finite difference (FD) approximation using Taylor expansion:
 *
 * For a small perturbation \f$ \epsilon \f$ and increment \f$ dx \f$:
 *
 * \f[
 *   H(x + \epsilon dx) \approx H(x) + \epsilon H dx
 * \f]
 *
 * The finite difference approximation is:
 *
 * \f[
 *   \frac{H(x + \epsilon dx) - H(x)}{\epsilon} \approx H dx
 * \f]
 *
 * The function computes the relative error between the TL and FD results for
 * multiple values of \f$ \epsilon \f$, and checks that the error decreases with
 * decreasing \f$ \epsilon \f$ and that the convergence rate is close to 1
 * (first-order accuracy).
 *
 * @tparam BackendTag The backend type tag
 * @param obs_op The observation operator
 * @param state The state vector
 * @param obs The observation vector
 * @param tol The tolerance for the final relative error (default: 1e-6)
 * @param epsilons The set of perturbation sizes to use (default: {1e-3, 1e-4,
 * 1e-5, 1e-6, 1e-7})
 * @return true if the check passes, false otherwise
 */
template <typename BackendTag>
bool checkObsOperatorTangentLinear(
    const ObsOperator<BackendTag>& obs_op, const State<BackendTag>& state,
    const Observation<BackendTag>& obs, double tol = 1e-6,
    const std::vector<double>& epsilons = std::vector<double>{1e-3, 1e-4, 1e-5,
                                                              1e-6, 1e-7}) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  // 1. Create random state increment
  auto dx = Increment<BackendTag>::createFromEntity(state);
  dx.randomize();
  // Normalize dx to unit norm
  double dx_norm = dx.norm();
  if (dx_norm > 0.0) {
    dx /= dx_norm;
  }

  // 2. Apply tangent linear: H dx
  auto Hdx_tl = obs_op.applyTangentLinear(dx, state, obs);

  // 3. Compute Taylor expansion with multiple perturbation sizes
  constexpr double min_epsilon = 1e-10;
  std::vector<double> errors;
  std::vector<double> used_epsilons;
  std::vector<double> convergence_rates;

  auto y0 = obs_op.apply(state, obs);

  for (size_t i = 0; i < epsilons.size(); ++i) {
    double epsilon = epsilons[i];
    if (epsilon < min_epsilon) continue;  // Ignore too-small epsilons

    // Compute perturbed state: x + ε·dx
    auto state_perturbed = state.clone();
    state_perturbed += dx * epsilon;

    // Compute finite difference: (H(x + ε·dx) - H(x)) / ε
    auto y1 = obs_op.apply(state_perturbed, obs);

    std::vector<double> Hdx_fd(y0.size());
    for (size_t j = 0; j < y0.size(); ++j) {
      Hdx_fd[j] = (y1[j] - y0[j]) / epsilon;
    }

    // Compute relative error for this epsilon
    double tl_norm = 0.0, fd_norm = 0.0, diff_norm = 0.0;

    for (size_t j = 0; j < Hdx_tl.size(); ++j) {
      double tl_val = Hdx_tl[j];
      double fd_val = Hdx_fd[j];
      double diff = tl_val - fd_val;

      tl_norm += tl_val * tl_val;
      fd_norm += fd_val * fd_val;
      diff_norm += diff * diff;
    }

    tl_norm = std::sqrt(tl_norm);
    fd_norm = std::sqrt(fd_norm);
    diff_norm = std::sqrt(diff_norm);

    double rel_error = diff_norm / (std::max(tl_norm, fd_norm) + 1e-12);
    errors.push_back(rel_error);
    used_epsilons.push_back(epsilon);

    logger.Info() << "ε = " << std::scientific << std::setprecision(1)
                  << epsilon << ", rel error = " << std::setprecision(13)
                  << std::scientific << rel_error;
  }

  // 4. Compute convergence rate (should be O(ε) for first-order accuracy)
  for (size_t i = 1; i < errors.size(); ++i) {
    double rate = std::log(errors[i - 1] / errors[i]) /
                  std::log(used_epsilons[i - 1] / used_epsilons[i]);
    convergence_rates.push_back(rate);
    logger.Info() << "Convergence rate " << i << ": " << std::setprecision(3)
                  << rate;
  }

  // 5. Check if errors decrease and convergence rate is reasonable
  bool errors_decrease = true;
  bool convergence_good = true;

  for (size_t i = 1; i < errors.size(); ++i) {
    if (errors[i] >= errors[i - 1]) {
      errors_decrease = false;
    }
  }

  for (size_t i = 0; i < convergence_rates.size(); ++i) {
    // Expect convergence rate around 1.0 (first-order accuracy)
    // Allow some tolerance for numerical issues
    if (convergence_rates[i] < 0.5 || convergence_rates[i] > 2.0) {
      convergence_good = false;
    }
  }

  // 6. Final check: use the smallest epsilon for the main tolerance check
  double final_error = errors.empty() ? 0.0 : errors.back();
  bool tol_check = final_error < tol;

  logger.Info() << "Tangent Linear check: final error = "
                << std::setprecision(13) << std::scientific << final_error
                << ", errors decrease = "
                << (errors_decrease ? "true" : "false")
                << ", convergence good = "
                << (convergence_good ? "true" : "false");

  return tol_check && errors_decrease && convergence_good;
}

// Gradient check for CostFunction
// Returns true if the check passes (relative error < tol)
template <typename BackendTag>
bool checkCostFunctionGradient(const CostFunction<BackendTag>& cost_func,
                               const State<BackendTag>& state,
                               double tol = 1e-6, double epsilon = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  auto grad = Increment<BackendTag>::createFromEntity(state);
  cost_func.gradient(state, grad);

  auto dx = Increment<BackendTag>::createFromEntity(state);
  dx.randomize();

  auto state_perturbed = state.clone();
  state_perturbed += dx * epsilon;

  double J0 = cost_func.evaluate(state);
  double J1 = cost_func.evaluate(state_perturbed);

  double fd = (J1 - J0) / epsilon;
  double analytic = grad.dot(dx);

  double rel_error = std::abs(fd - analytic) /
                     (std::max(std::abs(fd), std::abs(analytic)) + 1e-12);

  logger.Info() << "Gradient check: FD = " << std::setprecision(13)
                << std::scientific << fd << ", analytic = " << analytic
                << ", rel error = " << rel_error;

  return rel_error < tol;
}

}  // namespace metada::framework