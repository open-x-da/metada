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

// Tangent linear check for ObsOperator using Taylor expansion with multiple
// perturbation sizes Returns true if the check passes (relative error < tol and
// convergence order is correct)
template <typename BackendTag>
bool checkObsOperatorTangentLinear(const ObsOperator<BackendTag>& obs_op,
                                   const State<BackendTag>& state,
                                   const Observation<BackendTag>& obs,
                                   double tol = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  // 1. Create random state increment
  auto dx = Increment<BackendTag>::createFromEntity(state);
  dx.randomize();

  // 2. Apply tangent linear: H dx
  auto Hdx_tl = obs_op.applyTangentLinear(dx, state, obs);

  // 3. Compute Taylor expansion with multiple perturbation sizes
  std::vector<double> epsilons = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
  std::vector<double> errors;
  std::vector<double> convergence_rates;

  auto y0 = obs_op.apply(state, obs);

  for (size_t i = 0; i < epsilons.size(); ++i) {
    double epsilon = epsilons[i];

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

    logger.Info() << "ε = " << std::scientific << std::setprecision(1)
                  << epsilon << ", rel error = " << std::setprecision(13)
                  << std::scientific << rel_error;
  }

  // 4. Compute convergence rate (should be O(ε) for first-order accuracy)
  for (size_t i = 1; i < errors.size(); ++i) {
    double rate = std::log(errors[i - 1] / errors[i]) /
                  std::log(epsilons[i - 1] / epsilons[i]);
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
  double final_error = errors.back();
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