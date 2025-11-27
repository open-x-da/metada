#pragma once

#include <string>
#include <vector>

#include "LiteIncrement.hpp"
#include "LiteObs.hpp"
#include "LiteState.hpp"

namespace metada::backends::lite {

/**
 * @brief Lite observation operator backend for concrete testing
 *
 * Implements the same linear observation operator as LiteObs
 *
 * This class is designed to comply with the ObsOperatorBackendImpl concept.
 *
 * @tparam ControlVariableBackend Type of control variable backend
 */
template <typename ControlVariableBackend>
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
  explicit LiteObsOperator(
      [[maybe_unused]] const ConfigBackend& config,
      [[maybe_unused]] const ControlVariableBackend& control_backend)
      : control_backend_(nullptr) {
    required_state_vars_ = {"temperature", "pressure", "humidity"};
    required_obs_vars_ = {"temperature", "pressure"};
  }

  // Overload: accept config only (control backend not required for lite)
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

  // Tangent linear: dy = H*dx (operates on increments)
  std::vector<double> applyTangentLinear(
      const LiteIncrement& increment,
      [[maybe_unused]] const LiteState& reference_state,
      const LiteObs& obs) const {
    // For linear operator, tangent linear is same as forward
    auto data = increment.getData();
    std::vector<double> result(obs.size());
    for (int i = 0; i < 2; ++i) {
      result[i] = 0.0;
      for (int j = 0; j < 3; ++j) {
        result[i] += H[i][j] * data[j];
      }
    }
    return result;
  }

  // Adjoint: dx = H^T*dy (outputs an increment)
  void applyAdjoint(const std::vector<double>& obs_increment,
                    [[maybe_unused]] const LiteState& reference_state,
                    LiteIncrement& increment_result,
                    [[maybe_unused]] const LiteObs& obs) const {
    std::vector<double> result(3, 0.0);

    for (int i = 0; i < 3; ++i) {
      result[i] = 0.0;
      for (int j = 0; j < 2; ++j) {
        result[i] += H[j][i] * obs_increment[j];
      }
    }

    // Update increment with computed result
    increment_result.update(result.data(), nullptr, nullptr, nullptr, nullptr);
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

  // No control-space methods needed for lite backend (identity control)

 private:
  std::vector<std::string> required_state_vars_;
  std::vector<std::string> required_obs_vars_;
  const ControlVariableBackend* control_backend_ = nullptr;
};

}  // namespace metada::backends::lite