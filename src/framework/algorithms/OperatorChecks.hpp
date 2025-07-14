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

// Tangent linear check for ObsOperator using finite differences
// Returns true if the check passes (relative error < tol)
template <typename BackendTag>
bool checkObsOperatorTangentLinear(const ObsOperator<BackendTag>& obs_op,
                                   const State<BackendTag>& state,
                                   const Observation<BackendTag>& obs,
                                   double tol = 1e-6, double epsilon = 1e-6) {
  Logger<BackendTag>& logger = Logger<BackendTag>::Instance();

  // 1. Create random state increment
  auto dx = Increment<BackendTag>::createFromEntity(state);
  dx.randomize();

  // 2. Apply tangent linear: H dx
  auto Hdx_tl = obs_op.applyTangentLinear(dx, state, obs);

  // 3. Compute finite difference approximation
  auto state_perturbed = state.clone();
  state_perturbed += dx * epsilon;

  auto y0 = obs_op.apply(state, obs);
  auto y1 = obs_op.apply(state_perturbed, obs);

  std::vector<double> Hdx_fd(y0.size());
  for (size_t i = 0; i < y0.size(); ++i) {
    Hdx_fd[i] = (y1[i] - y0[i]) / epsilon;
  }

  // 4. Compare tangent linear with finite difference
  double tl_norm = 0.0, fd_norm = 0.0, diff_norm = 0.0;

  for (size_t i = 0; i < Hdx_tl.size(); ++i) {
    double tl_val = Hdx_tl[i];
    double fd_val = Hdx_fd[i];
    double diff = tl_val - fd_val;

    tl_norm += tl_val * tl_val;
    fd_norm += fd_val * fd_val;
    diff_norm += diff * diff;
  }

  tl_norm = std::sqrt(tl_norm);
  fd_norm = std::sqrt(fd_norm);
  diff_norm = std::sqrt(diff_norm);

  double rel_error = diff_norm / (std::max(tl_norm, fd_norm) + 1e-12);

  logger.Info() << "Tangent Linear check: TL norm = " << std::setprecision(13)
                << std::scientific << tl_norm << ", FD norm = " << fd_norm
                << ", diff norm = " << diff_norm
                << ", rel error = " << rel_error;

  return rel_error < tol;
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