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
  void controlToIncrement(const ControlVariableType& control,
                          IncrementType& increment) const override {
    // For identity: δx = v (just copy data)
    auto control_data = control.template getData<std::vector<double>>();
    increment.backend().setFromVector(control_data);
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
      const IncrementType& state_gradient,
      ControlVariableType& control_gradient) const override {
    // For identity: ∇_v J = ∇_δx J (just copy data)
    auto gradient_data = state_gradient.template getData<std::vector<double>>();
    control_gradient.setFromVector(gradient_data);
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
