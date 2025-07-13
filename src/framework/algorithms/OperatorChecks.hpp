#pragma once

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "CostFunction.hpp"
#include "Increment.hpp"
#include "ObsOperator.hpp"
#include "Observation.hpp"
#include "State.hpp"

namespace metada::framework {

// Utility to randomize an Increment (should be implemented for your Increment
// class)
template <typename IncrementT>
void randomize_increment(IncrementT& inc) {
  for (size_t i = 0; i < inc.size(); ++i) {
    inc[i] = (double(rand()) / RAND_MAX - 0.5);
  }
}

// Utility to compute dot product between two Increments or between Increment
// and State (should be implemented for your Increment class)
template <typename IncrementT>
double dot_product(const IncrementT& a, const IncrementT& b) {
  double result = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    result += a[i] * b[i];
  }
  return result;
}

// TL/AD check for ObsOperator
// Returns true if the check passes (relative error < tol)
template <typename BackendTag>
bool checkObsOperatorTLAD(const ObsOperator<BackendTag>& obs_op,
                          const State<BackendTag>& state,
                          const Observation<BackendTag>& obs, double tol = 1e-6,
                          double epsilon = 1e-6, int verbose = 1) {
  // 1. Create random state increment dx and obs increment dy
  auto dx = Increment<State<BackendTag>>::createFromEntity(state);
  randomize_increment(dx);

  std::vector<double> dy(obs.size());
  for (auto& v : dy) v = (double(rand()) / RAND_MAX - 0.5);

  // 2. Apply tangent linear: H dx
  auto Hdx = obs_op.applyTangentLinear(dx, state, obs);

  // 3. Apply adjoint: H^T dy
  auto HTdy = obs_op.applyAdjoint(dy, state);

  // 4. Compute inner products
  double a = std::inner_product(Hdx.begin(), Hdx.end(), dy.begin(), 0.0);
  double b = dot_product(dx, HTdy);

  double rel_error =
      std::abs(a - b) / (std::max(std::abs(a), std::abs(b)) + 1e-12);

  if (verbose) {
    std::cout << "TL/AD check: <Hdx, dy> = " << a << ", <dx, H^T dy> = " << b
              << ", rel error = " << rel_error << std::endl;
  }

  return rel_error < tol;
}

// Gradient check for CostFunction
// Returns true if the check passes (relative error < tol)
template <typename BackendTag>
bool checkCostFunctionGradient(const CostFunction<BackendTag>& cost_func,
                               const State<BackendTag>& state,
                               double tol = 1e-6, double epsilon = 1e-6,
                               int verbose = 1) {
  auto grad = Increment<State<BackendTag>>::createFromEntity(state);
  cost_func.gradient(state, grad);

  auto dx = Increment<State<BackendTag>>::createFromEntity(state);
  randomize_increment(dx);

  auto state_perturbed = state.clone();
  state_perturbed += dx * epsilon;

  double J0 = cost_func.evaluate(state);
  double J1 = cost_func.evaluate(state_perturbed);

  double fd = (J1 - J0) / epsilon;
  double analytic = dot_product(grad, dx);

  double rel_error = std::abs(fd - analytic) /
                     (std::max(std::abs(fd), std::abs(analytic)) + 1e-12);

  if (verbose) {
    std::cout << "Gradient check: FD = " << fd << ", analytic = " << analytic
              << ", rel error = " << rel_error << std::endl;
  }

  return rel_error < tol;
}

}  // namespace metada::framework