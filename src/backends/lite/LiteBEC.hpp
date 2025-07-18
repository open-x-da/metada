#pragma once

#include <cmath>
#include <vector>

#include "LiteState.hpp"

namespace metada::backends::lite {

/**
 * @brief Lite background error covariance backend for concrete testing
 *
 * Implements a diagonal B matrix with unit variances for simplicity
 */
class LiteBEC {
 public:
  LiteBEC() = default;

  // Diagonal B matrix (unit variances for simplicity)
  std::vector<double> variances_;

  // Constructor with config
  template <typename ConfigBackend>
  LiteBEC([[maybe_unused]] const ConfigBackend& config) {
    // Initialize with unit variances for testing
    variances_ = {1.0, 1.0, 1.0};
  }

  // Quadratic form: x^T * B^(-1) * x
  double computeQuadraticFormDiagonal(const LiteState& state) const {
    double result = 0.0;
    const double* data = static_cast<const double*>(state.getData());
    for (size_t i = 0; i < state.size(); ++i) {
      result += data[i] * data[i] / variances_[i];
    }
    return result;
  }

  // Apply inverse: B^(-1) * x
  void applyInverseDiagonal(const LiteState& state, LiteState& result) const {
    const double* input_data = static_cast<const double*>(state.getData());
    double* output_data = static_cast<double*>(result.getData());
    for (size_t i = 0; i < state.size(); ++i) {
      output_data[i] = input_data[i] / variances_[i];
    }
  }

  // Apply forward: B * x
  void applyDiagonal(const LiteState& state, LiteState& result) const {
    const double* input_data = static_cast<const double*>(state.getData());
    double* output_data = static_cast<double*>(result.getData());
    for (size_t i = 0; i < state.size(); ++i) {
      output_data[i] = input_data[i] * variances_[i];
    }
  }

  // Apply square root: B^(1/2) * x
  void applySquareRootDiagonal(const LiteState& state,
                               LiteState& result) const {
    const double* input_data = static_cast<const double*>(state.getData());
    double* output_data = static_cast<double*>(result.getData());
    for (size_t i = 0; i < state.size(); ++i) {
      output_data[i] = input_data[i] * std::sqrt(variances_[i]);
    }
  }

  // Ensemble-based methods (stub implementations)
  double computeQuadraticFormEnsemble(const LiteState& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  void applyInverseEnsemble(const LiteState& state, LiteState& result) const {
    applyInverseDiagonal(state, result);
  }

  void applyEnsemble(const LiteState& state, LiteState& result) const {
    applyDiagonal(state, result);
  }

  void applySquareRootEnsemble(const LiteState& state,
                               LiteState& result) const {
    applySquareRootDiagonal(state, result);
  }

  // Parametric methods (stub implementations)
  double computeQuadraticFormParametric(const LiteState& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  void applyInverseParametric(const LiteState& state, LiteState& result) const {
    applyInverseDiagonal(state, result);
  }

  void applyParametric(const LiteState& state, LiteState& result) const {
    applyDiagonal(state, result);
  }

  void applySquareRootParametric(const LiteState& state,
                                 LiteState& result) const {
    applySquareRootDiagonal(state, result);
  }

  // Hybrid methods (stub implementations)
  double computeQuadraticFormHybrid(const LiteState& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  void applyInverseHybrid(const LiteState& state, LiteState& result) const {
    applyInverseDiagonal(state, result);
  }

  void applyHybrid(const LiteState& state, LiteState& result) const {
    applyDiagonal(state, result);
  }

  void applySquareRootHybrid(const LiteState& state, LiteState& result) const {
    applySquareRootDiagonal(state, result);
  }

  // Full matrix methods (stub implementations)
  double computeQuadraticFormFull(const LiteState& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  void applyInverseFull(const LiteState& state, LiteState& result) const {
    applyInverseDiagonal(state, result);
  }

  void applyFull(const LiteState& state, LiteState& result) const {
    applyDiagonal(state, result);
  }

  void applySquareRootFull(const LiteState& state, LiteState& result) const {
    applySquareRootDiagonal(state, result);
  }

  // Test helper methods
  void setVariances(const std::vector<double>& vars) { variances_ = vars; }
};

}  // namespace metada::backends::lite