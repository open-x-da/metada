#pragma once

#include <string>

#include "ControlVariableBackend.hpp"

namespace metada::framework {

/**
 * @brief Control-variable backend that treats grid%xa increments as the control
 *        vector (current behavior).
 */
template <typename BackendTag>
class IdentityControlVariableBackend
    : public ControlVariableBackend<BackendTag> {
 public:
  using typename ControlVariableBackend<BackendTag>::IncrementType;
  using typename ControlVariableBackend<BackendTag>::GeometryBackendType;
  using typename ControlVariableBackend<BackendTag>::StateType;

  IdentityControlVariableBackend() = default;
  ~IdentityControlVariableBackend() override = default;

  IncrementType createIncrement(
      const GeometryBackendType& geometry) const override {
    return IncrementType::createFromGeometry(geometry);
  }

  void addIncrementToState(const IncrementType& increment,
                           StateType& state) const override {
    state += increment;
  }

  std::string name() const override { return "identity"; }

  typename ControlVariableBackend<BackendTag>::ControlVector
  createControlVector(const GeometryBackendType& geometry) const override {
    IncrementType temp = IncrementType::createFromGeometry(geometry);
    auto data = temp.template getData<std::vector<double>>();
    return std::vector<double>(data.size(), 0.0);
  }

  void convertStateIncrementToControl(
      const IncrementType& increment,
      typename ControlVariableBackend<BackendTag>::ControlVector& control)
      const override {
    auto data = increment.template getData<std::vector<double>>();
    control = std::move(data);
  }

  void convertControlToStateIncrement(
      const typename ControlVariableBackend<BackendTag>::ControlVector& control,
      IncrementType& increment) const override {
    increment.backend().setFromVector(control);
  }

  void convertStateGradientToControl(
      const IncrementType& state_gradient,
      typename ControlVariableBackend<BackendTag>::ControlVector&
          control_gradient) const override {
    convertStateIncrementToControl(state_gradient, control_gradient);
  }

  void convertControlGradientToState(
      const typename ControlVariableBackend<BackendTag>::ControlVector&
          control_gradient,
      IncrementType& state_gradient) const override {
    convertControlToStateIncrement(control_gradient, state_gradient);
  }
};

}  // namespace metada::framework
