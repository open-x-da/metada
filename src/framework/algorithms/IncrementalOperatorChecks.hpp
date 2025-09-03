#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include "BackendTraits.hpp"
#include "Increment.hpp"
#include "IncrementalCostFunction.hpp"
#include "Logger.hpp"

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
  auto direction = Increment<BackendTag>(test_increment.state());
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

  // Compute analytical gradient
  auto analytical_gradient = Increment<BackendTag>(test_increment.state());
  cost_function.gradient(test_increment, analytical_gradient);

  // Compute finite difference gradient
  double epsilon = 1e-8;
  auto perturbed_plus = Increment<BackendTag>(test_increment.state());
  auto perturbed_minus = Increment<BackendTag>(test_increment.state());

  // δx + ε*d
  perturbed_plus += direction * epsilon;
  // δx - ε*d
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
    auto direction = Increment<BackendTag>(test_increment.state());
    direction.randomize();

    // Normalize the direction
    double dir_norm = direction.norm();
    if (dir_norm > 0.0) {
      direction *= (1.0 / dir_norm);
    } else {
      continue;  // Skip zero norm directions
    }

    // Compute analytical gradient
    auto analytical_gradient = Increment<BackendTag>(test_increment.state());
    cost_function.gradient(test_increment, analytical_gradient);

    // Compute finite difference gradient
    double epsilon = 1e-8;
    auto perturbed_plus = Increment<BackendTag>(test_increment.state());
    auto perturbed_minus = Increment<BackendTag>(test_increment.state());

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

}  // namespace metada::framework
