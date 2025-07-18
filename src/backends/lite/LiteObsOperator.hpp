#pragma once

#include <string>
#include <vector>

#include "LiteObs.hpp"
#include "LiteState.hpp"

namespace metada::backends::lite {

/**
 * @brief Lite observation operator backend for concrete testing
 *
 * Implements the same linear observation operator as LiteObs
 *
 * This class is designed to comply with the ObsOperatorBackendImpl concept.
 */
class LiteObsOperator {
 public:
  // Deleted default constructor (concept compliance)
  LiteObsOperator() = delete;

  // Deleted copy constructor and copy assignment operator (concept compliance)
  LiteObsOperator(const LiteObsOperator&) = delete;
  LiteObsOperator& operator=(const LiteObsOperator&) = delete;

  LiteObsOperator(LiteObsOperator&&) noexcept = default;
  LiteObsOperator& operator=(LiteObsOperator&&) noexcept = default;

  // Explicit templated config constructor (concept compliance)
  template <typename ConfigBackend>
  explicit LiteObsOperator([[maybe_unused]] const ConfigBackend& config) {
    required_state_vars_ = {"temperature", "pressure", "humidity"};
    required_obs_vars_ = {"temperature", "pressure"};
  }

  // Simple 2x3 observation operator matrix (same as LiteObs)
  static constexpr double H[2][3] = {{1.0, 0.5, 0.0}, {0.0, 1.0, 0.5}};

  // Forward operator: y = H*x
  std::vector<double> apply(const LiteState& state,
                            [[maybe_unused]] const LiteObs& obs) const {
    std::vector<double> result(2);
    const double* state_data = static_cast<const double*>(state.getData());

    for (int i = 0; i < 2; ++i) {
      result[i] = 0.0;
      for (int j = 0; j < 3; ++j) {
        result[i] += H[i][j] * state_data[j];
      }
    }
    return result;
  }

  // Tangent linear: dy = H*dx
  std::vector<double> applyTangentLinear(
      const LiteState& state_increment,
      [[maybe_unused]] const LiteState& reference_state,
      const LiteObs& obs) const {
    // For linear operator, tangent linear is same as forward
    return apply(state_increment, obs);
  }

  // Adjoint: dx = H^T*dy
  void applyAdjoint(const std::vector<double>& obs_increment,
                    [[maybe_unused]] const LiteState& reference_state,
                    LiteState& state_increment,
                    [[maybe_unused]] const LiteObs& obs) const {
    double* state_data = static_cast<double*>(state_increment.getData());

    for (int i = 0; i < 3; ++i) {
      state_data[i] = 0.0;
      for (int j = 0; j < 2; ++j) {
        state_data[i] += H[j][i] * obs_increment[j];
      }
    }
  }

  // Required variables
  const std::vector<std::string>& getRequiredStateVars() const {
    return required_state_vars_;
  }

  const std::vector<std::string>& getRequiredObsVars() const {
    return required_obs_vars_;
  }

  // Capability checks
  bool supportsLinearization() const {
    return true;  // Linear operator supports tangent linear/adjoint
  }

  bool isLinear() const {
    return true;  // This is a linear operator
  }

  bool isInitialized() const {
    return true;  // Always initialized for lite backend
  }

  // Initialize
  template <typename ConfigBackend>
  void initialize([[maybe_unused]] const ConfigBackend& config) {
    // Stub implementation
  }

 private:
  std::vector<std::string> required_state_vars_;
  std::vector<std::string> required_obs_vars_;
};

}  // namespace metada::backends::lite