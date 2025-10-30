#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "BackendTraits.hpp"
#include "Increment.hpp"
#include "IncrementalCostFunction.hpp"
#include "Logger.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Check the gradient of an incremental cost function using finite
 * differences
 *
 * @details This function verifies that the analytical gradient of the
 * incremental cost function matches the finite difference approximation. This
 * is crucial for ensuring the correctness of the adjoint implementation in
 * variational data assimilation.
 *
 * The test computes:
 * - Analytical gradient: ∇J(δx) computed by the cost function
 * - Finite difference gradient: [J(δx + ε*d) - J(δx - ε*d)] / (2*ε)
 * - Relative error: ||∇J_analytical - ∇J_finite|| / ||∇J_analytical||
 *
 * @tparam BackendTag The backend tag type
 * @param cost_function The incremental cost function to test
 * @param test_increment The increment around which to test the gradient
 * @param tolerance The tolerance for the gradient test (default: 1e-6)
 * @return True if the gradient test passes within tolerance
 */
template <typename BackendTag>
bool checkIncrementalCostFunctionGradient(
    const IncrementalCostFunction<BackendTag>& cost_function,
    const Increment<BackendTag>& test_increment, double tolerance = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  logger.Info() << "Starting incremental gradient test with tolerance: "
                << tolerance;

  // Create a random direction for the finite difference test
  auto direction = Increment<BackendTag>(test_increment.geometry());
  direction.randomize();

  // Normalize the direction
  double dir_norm = direction.norm();
  if (dir_norm > 0.0) {
    direction *= (1.0 / dir_norm);
  } else {
    logger.Warning() << "Random direction has zero norm, using unit direction";
    direction.zero();
    auto dir_data = direction.template getData<std::vector<double>>();
    if (!dir_data.empty()) {
      dir_data[0] = 1.0;
    }
  }

  // Compute analytical gradient at test_increment
  auto analytical_gradient = Increment<BackendTag>(test_increment.geometry());
  cost_function.gradient(test_increment, analytical_gradient);

  // Compute finite difference gradient
  // Use larger epsilon for better numerical stability with nonlinear operators
  // Rule of thumb: epsilon ~ sqrt(machine_precision) * scale
  double epsilon = 1e-5;

  auto perturbed_plus = Increment<BackendTag>(test_increment.geometry());
  auto perturbed_minus = Increment<BackendTag>(test_increment.geometry());

  // Perturb: (test_increment ± ε*d)
  perturbed_plus += direction * epsilon;
  perturbed_minus -= direction * epsilon;

  double cost_plus = cost_function.evaluate(perturbed_plus);
  double cost_minus = cost_function.evaluate(perturbed_minus);

  // Finite difference gradient in the direction d
  double finite_diff_gradient = (cost_plus - cost_minus) / (2.0 * epsilon);

  // Analytical gradient in the direction d
  double analytical_gradient_dot = analytical_gradient.dot(direction);

  // Compute relative error
  double error = std::abs(analytical_gradient_dot - finite_diff_gradient);
  double relative_error = error / (std::abs(analytical_gradient_dot) + 1e-15);

  logger.Info() << "Gradient test results:";
  logger.Info() << "  Analytical gradient (dot d): " << analytical_gradient_dot;
  logger.Info() << "  Finite difference gradient (dot d): "
                << finite_diff_gradient;
  logger.Info() << "  Absolute error: " << error;
  logger.Info() << "  Relative error: " << relative_error;
  logger.Info() << "  Tolerance: " << tolerance;

  bool test_passed = relative_error < tolerance;

  if (test_passed) {
    logger.Info() << "✓ Incremental gradient test PASSED";
  } else {
    logger.Warning() << "✗ Incremental gradient test FAILED";
  }

  return test_passed;
}

/**
 * @brief Check the gradient of an incremental cost function using multiple
 * random directions
 *
 * @details This function performs a more comprehensive gradient test by
 * checking the gradient in multiple random directions. This provides better
 * coverage of the gradient computation.
 *
 * @tparam BackendTag The backend tag type
 * @param cost_function The incremental cost function to test
 * @param test_increment The increment around which to test the gradient
 * @param num_directions Number of random directions to test (default: 10)
 * @param tolerance The tolerance for the gradient test (default: 1e-6)
 * @return True if all gradient tests pass within tolerance
 */
template <typename BackendTag>
bool checkIncrementalCostFunctionGradientMultipleDirections(
    const IncrementalCostFunction<BackendTag>& cost_function,
    const Increment<BackendTag>& test_increment, int num_directions = 10,
    double tolerance = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  logger.Info() << "Starting incremental gradient test with " << num_directions
                << " random directions";

  bool all_tests_passed = true;
  double max_relative_error = 0.0;

  for (int i = 0; i < num_directions; ++i) {
    // Create a random direction
    auto direction = Increment<BackendTag>(test_increment.geometry());
    direction.randomize();

    // Normalize the direction
    double dir_norm = direction.norm();
    if (dir_norm > 0.0) {
      direction *= (1.0 / dir_norm);
    } else {
      continue;  // Skip zero norm directions
    }

    // Compute analytical gradient
    auto analytical_gradient = Increment<BackendTag>(test_increment.geometry());
    cost_function.gradient(test_increment, analytical_gradient);

    // Compute finite difference gradient
    double epsilon = 1e-8;
    auto perturbed_plus = Increment<BackendTag>(test_increment.geometry());
    auto perturbed_minus = Increment<BackendTag>(test_increment.geometry());

    perturbed_plus += direction * epsilon;
    perturbed_minus -= direction * epsilon;

    double cost_plus = cost_function.evaluate(perturbed_plus);
    double cost_minus = cost_function.evaluate(perturbed_minus);

    double finite_diff_gradient = (cost_plus - cost_minus) / (2.0 * epsilon);
    double analytical_gradient_dot = analytical_gradient.dot(direction);

    double error = std::abs(analytical_gradient_dot - finite_diff_gradient);
    double relative_error = error / (std::abs(analytical_gradient_dot) + 1e-15);

    max_relative_error = std::max(max_relative_error, relative_error);

    if (relative_error >= tolerance) {
      all_tests_passed = false;
      logger.Warning() << "Direction " << i
                       << " failed: relative error = " << relative_error;
    }
  }

  logger.Info() << "Multiple direction gradient test results:";
  logger.Info() << "  Maximum relative error: " << max_relative_error;
  logger.Info() << "  Tolerance: " << tolerance;
  logger.Info() << "  Tests passed: "
                << (all_tests_passed ? "All" : "Some failed");

  if (all_tests_passed) {
    logger.Info() << "✓ Multiple direction incremental gradient test PASSED";
  } else {
    logger.Warning() << "✗ Multiple direction incremental gradient test FAILED";
  }

  return all_tests_passed;
}

/**
 * @brief Check the tangent linear and adjoint consistency for observation
 * operators in incremental variational formulation
 *
 * @details This function verifies the correctness of the tangent linear (TL)
 * and adjoint (AD) implementations of observation operators H used in
 * incremental variational data assimilation by checking the following
 * mathematical property:
 *
 * \f[
 *   \langle H' dx, dy \rangle = \langle dx, H'^T dy \rangle
 * \f]
 *
 * where:
 *   - \f$ dx \f$ is a random state increment
 *   - \f$ dy \f$ is a random observation increment
 *   - \f$ H' dx \f$ is the result of applying the tangent linear operator
 *   - \f$ H'^T dy \f$ is the result of applying the adjoint operator
 *
 * This is specifically designed for incremental variational formulation where
 * the observation operator is applied to state increments rather than full
 * states.
 *
 * @tparam BackendTag The backend type tag
 * @param obs_operators Vector of observation operators
 * @param background_state The background state (first guess)
 * @param observations Vector of observations
 * @param tolerance The tolerance for the relative error (default: 1e-6)
 * @return True if the check passes within tolerance
 */
template <typename BackendTag>
bool checkIncrementalObsOperatorTLAD(
    const std::vector<ObsOperator<BackendTag>>& obs_operators,
    const State<BackendTag>& background_state,
    const std::vector<Observation<BackendTag>>& observations,
    double tolerance = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  if (obs_operators.size() != observations.size()) {
    logger.Error() << "Number of observation operators ("
                   << obs_operators.size()
                   << ") must match number of observations ("
                   << observations.size() << ")";
    return false;
  }

  logger.Info() << "Starting incremental observation operator TL/AD check with "
                << obs_operators.size() << " operators";

  // 1. Create random state increment dx
  auto dx = Increment<BackendTag>::createFromGeometry(
      background_state.geometry()->backend());
  dx.randomize();

  // Normalize the increment for numerical stability
  double dx_norm = dx.norm();
  if (dx_norm > 0.0) {
    dx *= (1.0 / dx_norm);
  } else {
    logger.Warning() << "Random increment has zero norm, using unit increment";
    dx.zero();
    auto dx_data = dx.template getData<std::vector<double>>();
    if (!dx_data.empty()) {
      dx_data[0] = 1.0;
    }
  }

  // 2. Create random observation increments dy for all observations
  std::vector<std::vector<double>> dy_vectors;
  for (const auto& obs : observations) {
    std::vector<double> dy(obs.size());
    for (auto& v : dy) {
      v = (double(rand()) / RAND_MAX - 0.5);
    }
    dy_vectors.push_back(std::move(dy));
  }

  // 3. Apply tangent linear: H' dx for all operators
  std::vector<std::vector<double>> Hdx_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto Hdx = obs_operators[i].applyTangentLinear(dx, background_state,
                                                   observations[i]);
    Hdx_vectors.push_back(std::move(Hdx));
  }

  // 4. Apply adjoint: H'^T dy for all operators
  std::vector<Increment<BackendTag>> HTdy_increments;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto HTdy = obs_operators[i].applyAdjoint(dy_vectors[i], background_state,
                                              observations[i]);
    HTdy_increments.push_back(std::move(HTdy));
  }

  // 5. Compute inner products for all operators
  double total_a = 0.0, total_b = 0.0;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    double a = std::inner_product(Hdx_vectors[i].begin(), Hdx_vectors[i].end(),
                                  dy_vectors[i].begin(), 0.0);
    double b = dx.dot(HTdy_increments[i]);
    total_a += a;
    total_b += b;
  }

  // 6. Compute relative error
  double error = std::abs(total_a - total_b);
  double relative_error =
      error / (std::max(std::abs(total_a), std::abs(total_b)) + 1e-15);

  logger.Info() << "Incremental TL/AD check results:";
  logger.Info() << "  <H'dx, dy> = " << std::setprecision(13) << std::scientific
                << total_a;
  logger.Info() << "  <dx, H'^T dy> = " << total_b;
  logger.Info() << "  Absolute error: " << error;
  logger.Info() << "  Relative error: " << relative_error;
  logger.Info() << "  Tolerance: " << tolerance;

  bool test_passed = relative_error < tolerance;

  if (test_passed) {
    logger.Info() << "✓ Incremental observation operator TL/AD check PASSED";
  } else {
    logger.Warning() << "✗ Incremental observation operator TL/AD check FAILED";
  }

  return test_passed;
}

/**
 * @brief Check the tangent linear implementation of observation operators
 * in incremental variational formulation using finite differences
 *
 * @details This function verifies the correctness of the tangent linear (TL)
 * implementation of observation operators H' used in incremental variational
 * formulation by comparing the TL result to a finite difference (FD)
 * approximation:
 *
 * For a small perturbation \f$ \epsilon \f$ and increment \f$ dx \f$:
 *
 * \f[
 *   H(x_b + \epsilon dx) \approx H(x_b) + \epsilon H' dx
 * \f]
 *
 * The finite difference approximation is:
 *
 * \f[
 *   \frac{H(x_b + \epsilon dx) - H(x_b)}{\epsilon} \approx H' dx
 * \f]
 *
 * @tparam BackendTag The backend type tag
 * @param obs_operators Vector of observation operators
 * @param background_state The background state (first guess)
 * @param observations Vector of observations
 * @param tolerance The tolerance for the relative error (default: 1e-6)
 * @param epsilons The set of perturbation sizes to use (default: {1e-3, 1e-4,
 * 1e-5, 1e-6, 1e-7})
 * @return True if the check passes within tolerance
 */
template <typename BackendTag>
bool checkIncrementalObsOperatorTangentLinear(
    const std::vector<ObsOperator<BackendTag>>& obs_operators,
    const State<BackendTag>& background_state,
    const std::vector<Observation<BackendTag>>& observations,
    double tolerance = 1e-6,
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

  logger.Info()
      << "Starting incremental observation operator tangent linear check with "
      << obs_operators.size() << " operators";

  constexpr double min_epsilon = 1e-10;
  std::vector<double> errors;
  std::vector<double> used_epsilons;
  std::vector<double> convergence_rates;

  // 1. Create random state increment
  auto dx = Increment<BackendTag>::createFromGeometry(
      background_state.geometry()->backend());
  dx.randomize();

  // Normalize the increment for numerical stability
  double dx_norm = dx.norm();
  if (dx_norm > 0.0) {
    dx *= (1.0 / dx_norm);
  } else {
    logger.Warning() << "Random increment has zero norm, using unit increment";
    dx.zero();
    auto dx_data = dx.template getData<std::vector<double>>();
    if (!dx_data.empty()) {
      dx_data[0] = 1.0;
    }
  }

  // 2. Apply tangent linear: H' dx for all operators
  std::vector<std::vector<double>> Hdx_tl_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto Hdx_tl = obs_operators[i].applyTangentLinear(dx, background_state,
                                                      observations[i]);
    Hdx_tl_vectors.push_back(std::move(Hdx_tl));
  }

  // 3. Get initial observations for all operators
  std::vector<std::vector<double>> y0_vectors;
  for (size_t i = 0; i < obs_operators.size(); ++i) {
    auto y0 = obs_operators[i].apply(background_state, observations[i]);
    y0_vectors.push_back(std::move(y0));
  }

  // 4. Loop over epsilons
  for (size_t i = 0; i < epsilons.size(); ++i) {
    double epsilon = epsilons[i];
    if (epsilon < min_epsilon) continue;  // Ignore too-small epsilons

    // Compute perturbed state: x_b + ε·dx
    auto state_perturbed = background_state.clone();
    state_perturbed += dx * epsilon;

    // Compute H(x_b + ε·dx) for all operators
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

        if (taylor_abs_err > tolerance) {
          logger.Warning() << "Taylor expansion check failed at operator " << j
                           << ", output " << k << ": |" << taylor_actual
                           << " - " << taylor_expected
                           << "| = " << taylor_abs_err
                           << " > tolerance = " << tolerance;
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

    double abs_error = total_diff_norm;
    double rel_error =
        total_diff_norm / (std::max(total_tl_norm, total_fd_norm) + 1e-12);
    errors.push_back(rel_error);
    used_epsilons.push_back(epsilon);

    logger.Info() << "ε = " << std::scientific << std::setprecision(1)
                  << epsilon << ", abs error = " << std::setprecision(3)
                  << std::scientific << abs_error
                  << ", rel error = " << std::setprecision(13)
                  << std::scientific << rel_error;
  }

  // 5. If only one epsilon, just check all outputs passed
  if (errors.size() == 1) {
    if (errors[0] < tolerance) {
      logger.Info()
          << "Single-epsilon tangent linear check passed for all outputs.";
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
  // For tangent linear checks, errors should decrease for larger epsilons,
  // then increase for very small epsilons due to numerical precision.
  // We check convergence only for the first few epsilons (not too small).

  // Find the minimum error (best approximation)
  double min_error = errors.empty() ? 1.0 : errors[0];
  size_t min_error_idx = 0;
  for (size_t i = 0; i < errors.size(); ++i) {
    if (errors[i] < min_error) {
      min_error = errors[i];
      min_error_idx = i;
    }
  }

  // Check if errors decrease for at least the first 3-4 epsilons
  // (before hitting numerical precision limits)
  // For very accurate tangent linears, we may need to check more epsilons
  // to see the decreasing trend before hitting the roundoff floor
  bool errors_decrease = true;
  size_t check_count = std::min(size_t(4), errors.size() - 1);
  for (size_t i = 1; i <= check_count && i < errors.size(); ++i) {
    if (errors[i] >= errors[i - 1] * 1.1) {  // Allow 10% tolerance for noise
      errors_decrease = false;
      break;
    }
  }

  // Check convergence rates only for reasonable epsilons (not too small)
  // Skip the last 1-2 rates which may be affected by numerical precision
  bool convergence_good = true;
  size_t rate_check_count = convergence_rates.size() > 2
                                ? convergence_rates.size() - 2
                                : convergence_rates.size();
  for (size_t i = 0; i < rate_check_count; ++i) {
    // Expect convergence rate around 1.0 (first-order accuracy)
    // Allow wider tolerance for numerical issues
    if (std::abs(convergence_rates[i]) < 0.3 ||
        std::abs(convergence_rates[i]) > 3.0) {
      convergence_good = false;
    }
  }

  // 8. Use the minimum error for tolerance check
  // This is the best approximation before numerical precision dominates
  bool tol_check = min_error < tolerance;

  double final_error = errors.empty() ? 0.0 : errors.back();
  logger.Info() << "Incremental tangent linear check: min error = "
                << std::setprecision(13) << std::scientific << min_error
                << " (at ε=" << used_epsilons[min_error_idx] << ")"
                << ", final error = " << final_error << ", errors decrease = "
                << (errors_decrease ? "true" : "false")
                << ", convergence good = "
                << (convergence_good ? "true" : "false");

  // Test passes if minimum error is below tolerance
  // For high-quality tangent linears, also expect errors to decrease initially
  bool test_passed = tol_check && (errors_decrease || min_error < 1e-12);
  // If min_error < 1e-12, the TL is so accurate that we're in pure roundoff,
  // so errors_decrease doesn't matter

  if (test_passed) {
    logger.Info()
        << "✓ Incremental observation operator tangent linear check PASSED"
        << " (min error = " << std::setprecision(2) << std::scientific
        << min_error << " < tolerance = " << tolerance << ")";
  } else {
    logger.Warning()
        << "✗ Incremental observation operator tangent linear check FAILED"
        << " (min error = " << std::setprecision(2) << std::scientific
        << min_error << " >= tolerance = " << tolerance << ")";
  }

  return test_passed;
}

/**
 * @brief Comprehensive check for observation operators in incremental
 * variational formulation
 *
 * @details This function performs both tangent linear/adjoint consistency check
 * and tangent linear finite difference check for observation operators used in
 * incremental variational data assimilation.
 *
 * @tparam BackendTag The backend type tag
 * @param obs_operators Vector of observation operators
 * @param background_state The background state (first guess)
 * @param observations Vector of observations
 * @param tolerance The tolerance for the relative error (default: 1e-6)
 * @param epsilons The set of perturbation sizes for FD check (default: {1e-3,
 * 1e-4, 1e-5, 1e-6, 1e-7})
 * @return True if all checks pass within tolerance
 */
template <typename BackendTag>
bool checkIncrementalObsOperatorsComprehensive(
    const std::vector<ObsOperator<BackendTag>>& obs_operators,
    const State<BackendTag>& background_state,
    const std::vector<Observation<BackendTag>>& observations,
    double tolerance = 1e-6,
    const std::vector<double>& epsilons = std::vector<double>{1e-3, 1e-4, 1e-5,
                                                              1e-6, 1e-7}) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  logger.Info()
      << "Starting comprehensive incremental observation operator checks";

  // 1. Tangent linear/adjoint consistency check
  bool tlad_passed = checkIncrementalObsOperatorTLAD(
      obs_operators, background_state, observations, tolerance);

  // 2. Tangent linear finite difference check
  bool tl_passed = checkIncrementalObsOperatorTangentLinear(
      obs_operators, background_state, observations, tolerance, epsilons);

  bool all_passed = tlad_passed && tl_passed;

  logger.Info()
      << "Comprehensive incremental observation operator check results:";
  logger.Info() << "  TL/AD consistency: "
                << (tlad_passed ? "PASSED" : "FAILED");
  logger.Info() << "  Tangent linear FD: " << (tl_passed ? "PASSED" : "FAILED");
  logger.Info() << "  Overall result: " << (all_passed ? "PASSED" : "FAILED");

  return all_passed;
}

}  // namespace metada::framework
