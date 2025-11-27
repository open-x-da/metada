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
  template <typename T>
  double computeQuadraticFormDiagonal(const T& state) const {
    double result = 0.0;
    auto data = state.getData();
    for (size_t i = 0; i < data.size(); ++i) {
      result += data[i] * data[i] / variances_[i];
    }
    return result;
  }

  // Apply inverse: B^(-1) * x
  template <typename T>
  void applyInverseDiagonal(const T& state, T& result) const {
    auto input_data = state.getData();
    auto result_data = result.getData();
    for (size_t i = 0; i < input_data.size(); ++i) {
      result_data[i] = input_data[i] / variances_[i];
    }
    // Update result with modified data
    if constexpr (requires { result.setData(result_data); }) {
      result.setData(result_data);
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
  template <typename T>
  double computeQuadraticFormEnsemble(const T& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  template <typename T>
  void applyInverseEnsemble(const T& state, T& result) const {
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
  template <typename T>
  double computeQuadraticFormParametric(const T& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  template <typename T>
  void applyInverseParametric(const T& state, T& result) const {
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
  template <typename T>
  double computeQuadraticFormHybrid(const T& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  template <typename T>
  void applyInverseHybrid(const T& state, T& result) const {
    applyInverseDiagonal(state, result);
  }

  void applyHybrid(const LiteState& state, LiteState& result) const {
    applyDiagonal(state, result);
  }

  void applySquareRootHybrid(const LiteState& state, LiteState& result) const {
    applySquareRootDiagonal(state, result);
  }

  // Full matrix methods (stub implementations)
  template <typename T>
  double computeQuadraticFormFull(const T& state) const {
    return computeQuadraticFormDiagonal(state);
  }

  template <typename T>
  void applyInverseFull(const T& state, T& result) const {
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