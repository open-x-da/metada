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