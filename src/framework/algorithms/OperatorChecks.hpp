#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "ControlVariable.hpp"
#include "ControlVariableBackend.hpp"
#include "CostFunction.hpp"
#include "IdentityControlVariableBackend.hpp"
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

  // 1. Create random control variable and obs increment dy
  // Use identity backend for these checks (control = increment)
  IdentityControlVariableBackend<BackendTag> identity_backend;
  auto control_var =
      identity_backend.createControlVariable(state.geometry()->backend());
  control_var.randomize();

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
    auto Hdx = obs_operators[i].applyTangentLinear(control_var, state,
                                                   observations[i]);
    Hdx_vectors.push_back(std::move(Hdx));
  }

  // 3. Apply adjoint: H^T dy for all operators
  std::vector<ControlVariable<BackendTag>> HTdy_controls;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto HTdy =
        obs_operators[i].applyAdjoint(dy_vectors[i], state, observations[i]);
    HTdy_controls.push_back(std::move(HTdy));
  }

  // 4. Compute inner products for all operators
  double total_a = 0.0, total_b = 0.0;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    double a = std::inner_product(Hdx_vectors[i].begin(), Hdx_vectors[i].end(),
                                  dy_vectors[i].begin(), 0.0);
    double b = control_var.dot(HTdy_controls[i]);
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
 * @brief Unified check for the tangent linear implementation of ObsOperators
 * using Taylor expansion and finite differences.
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
 * For each epsilon, the function checks both the Taylor expansion and FD forms
 * for all outputs. If only one epsilon is provided, it returns true if all
 * outputs pass the tolerance. If multiple, it also checks that the FD error
 * decreases and the convergence rate is close to 1.
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

  constexpr double min_epsilon = 1e-10;
  std::vector<double> errors;
  std::vector<double> used_epsilons;
  std::vector<double> convergence_rates;

  // 1. Create random control variable
  // Use identity backend for these checks (control = increment)
  IdentityControlVariableBackend<BackendTag> identity_backend;
  auto control_var =
      identity_backend.createControlVariable(state.geometry()->backend());
  control_var.randomize();

  // Also create increment for state perturbation
  auto dx = identity_backend.createIncrement(state.geometry()->backend());
  identity_backend.controlToIncrement(control_var, dx);

  // 2. Apply tangent linear: H dx for all operators
  std::vector<std::vector<double>> Hdx_tl_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto Hdx_tl = obs_operators[i].applyTangentLinear(control_var, state,
                                                      observations[i]);
    Hdx_tl_vectors.push_back(std::move(Hdx_tl));
  }

  // 3. Get initial observations for all operators
  std::vector<std::vector<double>> y0_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto y0 = obs_operators[i].apply(state, observations[i]);
    y0_vectors.push_back(std::move(y0));
  }

  // 4. Loop over epsilons (single or multiple)
  for (size_t i = 0; i < epsilons.size(); ++i) {
    double epsilon = epsilons[i];
    if (epsilon < min_epsilon) continue;  // Ignore too-small epsilons

    // Compute perturbed state: x + ε·dx
    auto state_perturbed = state.clone();
    state_perturbed += dx * epsilon;

    // Compute H(x + ε·dx) for all operators
    std::vector<std::vector<double>> y1_vectors;
    for (size_t j = 0; j < obs_operators.size(); ++j) {
      auto y1 = obs_operators[j].apply(state_perturbed, observations[j]);
      y1_vectors.push_back(std::move(y1));
    }

    // For each operator and output, check both Taylor expansion and FD forms
    double total_tl_norm = 0.0, total_fd_norm = 0.0, total_diff_norm = 0.0;
    for (size_t j = 0; j < obs_operators.size(); ++j) {
      for (size_t k = 0; k < y0_vectors[j].size(); ++k) {
        double taylor_expected =
            y0_vectors[j][k] + Hdx_tl_vectors[j][k] * epsilon;
        double taylor_actual = y1_vectors[j][k];
        double taylor_abs_err = std::abs(taylor_actual - taylor_expected);
        if (taylor_abs_err > tol) {
          logger.Error() << "Taylor expansion check failed at output " << k
                         << ": |" << taylor_actual << " - " << taylor_expected
                         << "| = " << taylor_abs_err << " > tol = " << tol;
        }

        // FD check
        double fd_approx = (y1_vectors[j][k] - y0_vectors[j][k]) / epsilon;
        double tl_val = Hdx_tl_vectors[j][k];
        double diff = tl_val - fd_approx;
        total_tl_norm += tl_val * tl_val;
        total_fd_norm += fd_approx * fd_approx;
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

  // 5. If only one epsilon, just check all outputs passed
  if (errors.size() == 1) {
    if (errors[0] < tol) {
      logger.Info() << "Single-epsilon Taylor/FD check passed for all outputs.";
      return true;
    } else {
      logger.Error() << "Single-epsilon check failed: rel error = "
                     << errors[0];
      return false;
    }
  }

  // Check if all errors are at or below machine precision (linear operator
  // case)
  bool all_errors_machine_prec = std::all_of(
      errors.begin(), errors.end(), [](double e) { return e < 1e-12; });
  if (all_errors_machine_prec) {
    logger.Info() << "All errors are at or below machine precision. Operator "
                     "is likely linear; skipping convergence rate check.";
    return true;
  }

  // 6. Compute convergence rate (should be O(ε) for first-order accuracy)
  for (size_t i = 1; i < errors.size(); ++i) {
    double rate = std::log(errors[i - 1] / errors[i]) /
                  std::log(used_epsilons[i - 1] / used_epsilons[i]);
    convergence_rates.push_back(rate);
    logger.Info() << "Convergence rate " << i << ": " << std::setprecision(3)
                  << rate;
  }

  // 7. Check if errors decrease and convergence rate is reasonable
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

  // 8. Final check: use the smallest epsilon for the main tolerance check
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
                               double tol = 1e-6, double epsilon = 1e-2) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  auto grad =
      Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
  cost_func.gradient(state, grad);

  auto dx =
      Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
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

// Improved gradient check with variable-specific epsilon
// Returns true if the check passes (relative error < tol)
template <typename BackendTag>
bool checkCostFunctionGradientImproved(
    const CostFunction<BackendTag>& cost_func, const State<BackendTag>& state,
    double tol = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  auto grad =
      Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
  cost_func.gradient(state, grad);

  // Create perturbation with variable-specific scaling
  auto dx =
      Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
  dx.randomize();

  // Scale perturbation by variable magnitudes for better numerical stability
  // This ensures epsilon is appropriate for each variable's scale
  double state_norm = state.norm();
  if (state_norm > 0) {
    dx *= (state_norm * 1e-6);  // Scale by 1e-6 of state magnitude
  } else {
    dx *= 1e-6;  // Fallback for zero state
  }

  auto state_perturbed = state.clone();
  state_perturbed += dx;

  double J0 = cost_func.evaluate(state);
  double J1 = cost_func.evaluate(state_perturbed);

  double dx_norm = dx.norm();
  double fd = (J1 - J0) / dx_norm;
  double analytic = grad.dot(dx) / dx_norm;

  double rel_error = std::abs(fd - analytic) /
                     (std::max(std::abs(fd), std::abs(analytic)) + 1e-12);

  logger.Info() << "Improved gradient check: FD = " << std::setprecision(13)
                << std::scientific << fd << ", analytic = " << analytic
                << ", rel error = " << rel_error
                << ", perturbation norm = " << dx_norm;

  return rel_error < tol;
}

// Comprehensive gradient check with multiple random directions
// Returns true if all checks pass (relative error < tol)
template <typename BackendTag>
bool checkCostFunctionGradientMultipleDirections(
    const CostFunction<BackendTag>& cost_func, const State<BackendTag>& state,
    size_t num_directions = 10, double tol = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  logger.Info() << "Testing gradient with " << num_directions
                << " random directions...";

  std::vector<double> rel_errors;
  std::vector<double> fd_values;
  std::vector<double> analytic_values;

  for (size_t i = 0; i < num_directions; ++i) {
    auto grad =
        Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
    cost_func.gradient(state, grad);

    // Create perturbation with variable-specific scaling
    auto dx =
        Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
    dx.randomize();

    // Scale perturbation by variable magnitudes for better numerical stability
    double state_norm = state.norm();
    if (state_norm > 0) {
      dx *= (state_norm * 1e-6);  // Scale by 1e-6 of state magnitude
    } else {
      dx *= 1e-6;  // Fallback for zero state
    }

    auto state_perturbed = state.clone();
    state_perturbed += dx;

    double J0 = cost_func.evaluate(state);
    double J1 = cost_func.evaluate(state_perturbed);

    double dx_norm = dx.norm();
    double fd = (J1 - J0) / dx_norm;
    double analytic = grad.dot(dx) / dx_norm;

    double rel_error = std::abs(fd - analytic) /
                       (std::max(std::abs(fd), std::abs(analytic)) + 1e-12);

    rel_errors.push_back(rel_error);
    fd_values.push_back(fd);
    analytic_values.push_back(analytic);

    logger.Info() << "Direction " << (i + 1) << "/" << num_directions
                  << ": FD = " << std::setprecision(13) << std::scientific << fd
                  << ", analytic = " << analytic
                  << ", rel error = " << rel_error;
  }

  // Statistical analysis
  double max_error = *std::max_element(rel_errors.begin(), rel_errors.end());
  double min_error = *std::min_element(rel_errors.begin(), rel_errors.end());
  double mean_error =
      std::accumulate(rel_errors.begin(), rel_errors.end(), 0.0) /
      rel_errors.size();

  // Calculate standard deviation
  double variance = 0.0;
  for (double error : rel_errors) {
    variance += (error - mean_error) * (error - mean_error);
  }
  double std_dev = std::sqrt(variance / rel_errors.size());

  logger.Info() << "Gradient check statistics:";
  logger.Info() << "  Max relative error: " << std::setprecision(13)
                << std::scientific << max_error;
  logger.Info() << "  Min relative error: " << min_error;
  logger.Info() << "  Mean relative error: " << mean_error;
  logger.Info() << "  Std dev relative error: " << std_dev;
  logger.Info() << "  Pass rate: "
                << std::count_if(rel_errors.begin(), rel_errors.end(),
                                 [tol](double e) { return e < tol; })
                << "/" << num_directions;

  return max_error < tol;
}

// Test specific directions (e.g., unit vectors along each variable)
template <typename BackendTag>
bool checkCostFunctionGradientUnitDirections(
    const CostFunction<BackendTag>& cost_func, const State<BackendTag>& state,
    double tol = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  logger.Info() << "Testing gradient with unit vector directions...";

  auto grad =
      Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
  cost_func.gradient(state, grad);

  std::vector<double> rel_errors;
  size_t state_size = state.size();

  // Create dx and a temporary state for setting unit vectors
  auto dx =
      Increment<BackendTag>::createFromGeometry(state.geometry()->backend());
  auto temp_state = state.clone();
  temp_state.zero();

  // Test unit vectors along each dimension
  for (size_t i = 0; i < std::min(state_size, size_t(10));
       ++i) {  // Limit to first 10 dimensions
    // Create unit vector in direction i using temporary state
    temp_state.zero();
    auto* data = temp_state.template getDataPtr<double>();
    data[i] = 1.0;

    // Transfer to increment
    dx.backend().transferFromState(temp_state.backend());

    // Scale the unit vector
    double state_norm = state.norm();
    if (state_norm > 0) {
      dx *= (state_norm * 1e-6);
    } else {
      dx *= 1e-6;
    }

    auto state_perturbed = state.clone();
    state_perturbed += dx;

    double J0 = cost_func.evaluate(state);
    double J1 = cost_func.evaluate(state_perturbed);

    double dx_norm = dx.norm();
    double fd = (J1 - J0) / dx_norm;
    double analytic = grad.dot(dx) / dx_norm;

    double rel_error = std::abs(fd - analytic) /
                       (std::max(std::abs(fd), std::abs(analytic)) + 1e-12);

    rel_errors.push_back(rel_error);

    logger.Info() << "Unit direction " << i
                  << ": FD = " << std::setprecision(13) << std::scientific << fd
                  << ", analytic = " << analytic
                  << ", rel error = " << rel_error;
  }

  double max_error = *std::max_element(rel_errors.begin(), rel_errors.end());
  logger.Info() << "Unit direction test - Max relative error: "
                << std::setprecision(13) << std::scientific << max_error;

  return max_error < tol;
}

}  // namespace metada::framework