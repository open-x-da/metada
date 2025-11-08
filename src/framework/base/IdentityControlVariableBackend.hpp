#pragma once

#include <string>

#include "ControlVariableBackend.hpp"

namespace metada::framework {

/**
 * @brief Control-variable backend that treats grid%xa increments as the control
 *        vector (current behavior).
 */
template <typename BackendTag>
class IdentityControlVariableBackend : public ControlVariableBackend<BackendTag> {
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
};

}  // namespace metada::framework


