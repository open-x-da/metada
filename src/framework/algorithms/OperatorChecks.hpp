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
 * @brief Check the tangent linear and adjoint consistency for ObsOperators
 * using the inner product test.
 *
 * This function verifies the correctness of the tangent linear (TL) and adjoint
 * (AD) implementations of observation operators H by checking the following
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
 * @param obs_operators Vector of observation operators
 * @param state The state vector
 * @param observations Vector of observations
 * @param tol The tolerance for the relative error (default: 1e-6)
 * @return true if the check passes, false otherwise
 */
template <typename BackendTag>
bool checkObsOperatorTLAD(
    const std::vector<ObsOperator<BackendTag>>& obs_operators,
    const State<BackendTag>& state,
    const std::vector<Observation<BackendTag>>& observations,
    double tol = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  if (obs_operators.size() != observations.size()) {
    logger.Error() << "Number of observation operators ("
                   << obs_operators.size()
                   << ") must match number of observations ("
                   << observations.size() << ")";
    return false;
  }

  // 1. Create random state increment dx and obs increment dy
  auto dx = Increment<BackendTag>::createFromEntity(state);
  dx.randomize();

  // Create observation increments for all observations
  std::vector<std::vector<double>> dy_vectors;
  for (const auto& obs : observations) {
    std::vector<double> dy(obs.size());
    for (auto& v : dy) v = (double(rand()) / RAND_MAX - 0.5);
    dy_vectors.push_back(std::move(dy));
  }

  // 2. Apply tangent linear: H dx for all operators
  std::vector<std::vector<double>> Hdx_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto Hdx = obs_operators[i].applyTangentLinear(dx, state, observations[i]);
    Hdx_vectors.push_back(std::move(Hdx));
  }

  // 3. Apply adjoint: H^T dy for all operators
  std::vector<Increment<BackendTag>> HTdy_increments;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto HTdy =
        obs_operators[i].applyAdjoint(dy_vectors[i], state, observations[i]);
    HTdy_increments.push_back(std::move(HTdy));
  }

  // 4. Compute inner products for all operators
  double total_a = 0.0, total_b = 0.0;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    double a = std::inner_product(Hdx_vectors[i].begin(), Hdx_vectors[i].end(),
                                  dy_vectors[i].begin(), 0.0);
    double b = dx.dot(HTdy_increments[i]);
    total_a += a;
    total_b += b;
  }

  double rel_error = std::abs(total_a - total_b) /
                     (std::max(std::abs(total_a), std::abs(total_b)) + 1e-12);

  logger.Info() << "TL/AD check: <Hdx, dy> = " << std::setprecision(13)
                << std::scientific << total_a << ", <dx, H^T dy> = " << total_b
                << ", rel error = " << rel_error;

  return rel_error < tol;
}

/**
 * @brief Check the tangent linear implementation of ObsOperators using Taylor
 * expansion and finite differences.
 *
 * This function verifies the correctness of the tangent linear (TL)
 * implementation of observation operators H by comparing the TL result to a
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
 * @param obs_operators Vector of observation operators
 * @param state The state vector
 * @param observations Vector of observations
 * @param tol The tolerance for the final relative error (default: 1e-6)
 * @param epsilons The set of perturbation sizes to use (default: {1e-3, 1e-4,
 * 1e-5, 1e-6, 1e-7})
 * @return true if the check passes, false otherwise
 */
template <typename BackendTag>
bool checkObsOperatorTangentLinear(
    const std::vector<ObsOperator<BackendTag>>& obs_operators,
    const State<BackendTag>& state,
    const std::vector<Observation<BackendTag>>& observations, double tol = 1e-6,
    const std::vector<double>& epsilons = std::vector<double>{1e-3, 1e-4, 1e-5,
                                                              1e-6, 1e-7}) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  if (obs_operators.size() != observations.size()) {
    logger.Error() << "Number of observation operators ("
                   << obs_operators.size()
                   << ") must match number of observations ("
                   << observations.size() << ")";
    return false;
  }

  // 1. Create random state increment
  auto dx = Increment<BackendTag>::createFromEntity(state);
  dx.randomize();
  // Normalize dx to unit norm
  double dx_norm = dx.norm();
  if (dx_norm > 0.0) {
    dx /= dx_norm;
  }

  // 2. Apply tangent linear: H dx for all operators
  std::vector<std::vector<double>> Hdx_tl_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto Hdx_tl =
        obs_operators[i].applyTangentLinear(dx, state, observations[i]);
    Hdx_tl_vectors.push_back(std::move(Hdx_tl));
  }

  // 3. Compute Taylor expansion with multiple perturbation sizes
  constexpr double min_epsilon = 1e-10;
  std::vector<double> errors;
  std::vector<double> used_epsilons;
  std::vector<double> convergence_rates;

  // Get initial observations for all operators
  std::vector<std::vector<double>> y0_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto y0 = obs_operators[i].apply(state, observations[i]);
    y0_vectors.push_back(std::move(y0));
  }

  for (size_t i = 0; i < epsilons.size(); ++i) {
    double epsilon = epsilons[i];
    if (epsilon < min_epsilon) continue;  // Ignore too-small epsilons

    // Compute perturbed state: x + ε·dx
    auto state_perturbed = state.clone();
    state_perturbed += dx * epsilon;

    // Compute finite difference for all operators: (H(x + ε·dx) - H(x)) / ε
    std::vector<std::vector<double>> Hdx_fd_vectors;
    for (size_t j = 0; j < obs_operators.size(); ++j) {
      auto y1 = obs_operators[j].apply(state_perturbed, observations[j]);

      std::vector<double> Hdx_fd(y0_vectors[j].size());
      for (size_t k = 0; k < y0_vectors[j].size(); ++k) {
        Hdx_fd[k] = (y1[k] - y0_vectors[j][k]) / epsilon;
      }
      Hdx_fd_vectors.push_back(std::move(Hdx_fd));
    }

    // Compute relative error for this epsilon (aggregate across all operators)
    double total_tl_norm = 0.0, total_fd_norm = 0.0, total_diff_norm = 0.0;

    for (size_t j = 0; j < obs_operators.size(); ++j) {
      for (size_t k = 0; k < Hdx_tl_vectors[j].size(); ++k) {
        double tl_val = Hdx_tl_vectors[j][k];
        double fd_val = Hdx_fd_vectors[j][k];
        double diff = tl_val - fd_val;

        total_tl_norm += tl_val * tl_val;
        total_fd_norm += fd_val * fd_val;
        total_diff_norm += diff * diff;
      }
    }

    total_tl_norm = std::sqrt(total_tl_norm);
    total_fd_norm = std::sqrt(total_fd_norm);
    total_diff_norm = std::sqrt(total_diff_norm);

    double rel_error =
        total_diff_norm / (std::max(total_tl_norm, total_fd_norm) + 1e-12);
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