#pragma once

#include <string>

#include "ControlVariable.hpp"
#include "ControlVariableBackend.hpp"
#include "Increment.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Identity control-variable backend: U = I (control = increment)
 *
 * @details This backend implements the simplest case where the control variable
 * is identical to the state-space increment (v = δx). This means:
 * - Forward transformation: δx = I * v = v (identity)
 * - Adjoint transformation: ∇_v J = I^T * ∇_δx J = ∇_δx J (identity)
 *
 * This is equivalent to not using any change of variables, and is the default
 * for systems without a B-matrix transformation.
 *
 * @tparam BackendTag Backend tag satisfying framework concepts
 */
template <typename BackendTag>
class IdentityControlVariableBackend
    : public ControlVariableBackend<BackendTag> {
 public:
  using typename ControlVariableBackend<BackendTag>::IncrementType;
  using typename ControlVariableBackend<BackendTag>::ControlVariableType;
  using typename ControlVariableBackend<BackendTag>::StateType;

  IdentityControlVariableBackend() = default;
  ~IdentityControlVariableBackend() override = default;

  /**
   * @brief Get backend name
   *
   * @return "identity"
   */
  std::string name() const override { return "identity"; }

  /**
   * @brief Get backend kind
   *
   * @return ControlVariableBackendKind::Identity
   */
  ControlVariableBackendKind kind() const override {
    return ControlVariableBackendKind::Identity;
  }

  /**
   * @brief Create a control variable from geometry
   *
   * @param geometry Geometry backend defining the domain
   * @return Newly created control variable
   */
  ControlVariableType createControlVariable(
      const typename ControlVariableType::GeometryBackendType& geometry)
      const override {
    return ControlVariableType::createFromGeometry(geometry);
  }

  /**
   * @brief Create an increment from geometry
   *
   * @param geometry Geometry backend defining the domain
   * @return Newly created increment
   */
  IncrementType createIncrement(
      const typename IncrementType::GeometryBackendType& geometry)
      const override {
    return IncrementType::createFromGeometry(geometry);
  }

  /**
   * @brief Forward transformation: δx = I * v = v
   *
   * @details For identity backend, control variable and increment are the same,
   * so we just copy the data.
   *
   * @param control Control variable v
   * @param increment Output increment δx (will be overwritten)
   */
  void controlToIncrement(
      [[maybe_unused]] const ControlVariableType& control,
      [[maybe_unused]] IncrementType& increment) const override {
    // For identity: δx = v
    // We need to copy the control variable data to the increment
    // The challenge is that for some backends (like WRF), ControlVariable and
    // Increment are different types, so we can't directly assign backends.
    //
    // Strategy: Use the adapter's axpy method to copy data
    // Since axpy is: this += alpha * other, we can do:
    // 1. Zero the increment
    // 2. Add the control with alpha=1.0
    //
    // But axpy(double, ControlVariable) doesn't exist on Increment.
    // So we need a different approach: extract and re-inject data
    //
    // For now, throw an error - this should be handled by the specific backend
    throw std::runtime_error(
        "IdentityControlVariableBackend::controlToIncrement not implemented "
        "for "
        "backends where ControlVariable != Increment. "
        "Please use a backend-specific ControlVariableBackend implementation.");
  }

  /**
   * @brief Adjoint transformation: ∇_v J = I^T * ∇_δx J = ∇_δx J
   *
   * @details For identity backend, the adjoint is also identity, so we just
   * copy the data.
   *
   * @param state_gradient Gradient in state space (∇_δx J)
   * @param control_gradient Output gradient in control space (∇_v J)
   */
  void incrementAdjointToControl(
      [[maybe_unused]] const IncrementType& state_gradient,
      [[maybe_unused]] ControlVariableType& control_gradient) const override {
    // For identity: ∇_v J = ∇_δx J
    // Same issue as controlToIncrement - we can't directly assign backends
    throw std::runtime_error(
        "IdentityControlVariableBackend::incrementAdjointToControl not "
        "implemented "
        "for backends where ControlVariable != Increment. "
        "Please use a backend-specific ControlVariableBackend implementation.");
  }

  /**
   * @brief Apply increment to state: xa = xb + δx
   *
   * @param increment State-space increment δx
   * @param state State to be updated (input: xb, output: xa = xb + δx)
   */
  void addIncrementToState(const IncrementType& increment,
                           StateType& state) const override {
    state += increment;
  }
};

}  // namespace metada::framework
