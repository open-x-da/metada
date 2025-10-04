#pragma once

#include "LiteState.hpp"

namespace metada::backends::lite {

/**
 * @brief Lite model backend for concrete testing
 *
 * Implements a simple linear model: x(t+1) = M*x(t) where M is identity
 *
 * This class is designed to comply with the ModelBackendImpl concept.
 * - Deleted default constructor
 * - Deleted copy constructor and copy assignment operator
 * - Config constructor (templated for flexibility)
 */
class LiteModel {
 public:
  // Deleted default constructor (concept compliance)
  LiteModel() = delete;

  // Deleted copy constructor and copy assignment operator (concept compliance)
  LiteModel(const LiteModel&) = delete;
  LiteModel& operator=(const LiteModel&) = delete;

  // Templated config constructor (for concept compliance and flexibility)
  template <typename ConfigBackend>
  explicit LiteModel([[maybe_unused]] const ConfigBackend& config) {
    // Initialize as identity model
  }

  // Forward model: x(t+1) = M*x(t)
  void run(const LiteState& initial_state, LiteState& final_state) const {
    const double* input_data =
        static_cast<const double*>(initial_state.getData());
    double* output_data = static_cast<double*>(final_state.getData());

    // Identity model for simplicity
    for (size_t i = 0; i < initial_state.size(); ++i) {
      output_data[i] = input_data[i];
    }
  }

  // Tangent linear model: dx(t+1) = M*dx(t)
  void runTangentLinear([[maybe_unused]] const LiteState& reference_initial,
                        [[maybe_unused]] const LiteState& reference_final,
                        const LiteState& initial_perturbation,
                        LiteState& final_perturbation) const {
    // For identity model, tangent linear is same as forward
    run(initial_perturbation, final_perturbation);
  }

  // Adjoint model: dx(t) = M^T*dx(t+1)
  void runAdjoint([[maybe_unused]] const LiteState& initial_state,
                  [[maybe_unused]] const LiteState& final_state,
                  const LiteState& adjoint_forcing,
                  LiteState& adjoint_result) const {
    // For identity model, adjoint is same as forward
    run(adjoint_forcing, adjoint_result);
  }

  // Capability checks
  bool supportsAdjoint() const { return true; }

  bool supportsTangentLinear() const { return true; }

  bool isLinear() const {
    return true;  // Identity model is linear
  }

  bool isInitialized() const {
    return true;  // Always initialized for lite backend
  }

  // Initialize
  template <typename ConfigBackend>
  void initialize([[maybe_unused]] const ConfigBackend& config) {
    // Stub implementation
  }

  // Reset
  void reset() {
    // Stub implementation
  }

  // Finalize
  void finalize() {
    // Stub implementation
  }
};

}  // namespace metada::backends::lite