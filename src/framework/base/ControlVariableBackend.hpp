#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Increment.hpp"
#include "NonCopyable.hpp"
#include "State.hpp"

namespace metada::framework {

/**
 * @brief Abstract interface describing control-variable backend behavior.
 *
 * @details Concrete implementations encapsulate how control vectors are stored,
 *          created, and mapped back to model space. This abstraction allows the
 *          incremental variational workflow to remain agnostic to the chosen
 *          control-variable representation (e.g., direct grid%xa or WRFDA CV5).
 *
 * @tparam BackendTag Backend tag satisfying framework concepts.
 */
template <typename BackendTag>
class ControlVariableBackend : public NonCopyable {
 public:
  using IncrementType = Increment<BackendTag>;
  using StateType = State<BackendTag>;
  using GeometryBackendType = typename IncrementType::GeometryBackendType;
  using ControlVector = std::vector<double>;

  ControlVariableBackend() = default;
  virtual ~ControlVariableBackend() = default;

  /**
   * @brief Create a control increment matching the provided geometry.
   */
  virtual IncrementType createIncrement(
      const GeometryBackendType& geometry) const = 0;

  /**
   * @brief Apply an increment to the provided state (xa = xb + δx).
   */
  virtual void addIncrementToState(const IncrementType& increment,
                                   StateType& state) const = 0;

  /**
   * @brief Human-readable backend name.
   */
  virtual std::string name() const = 0;

  /**
   * @brief Create a control-vector container sized for the provided geometry.
   */
  virtual ControlVector createControlVector(
      const GeometryBackendType& geometry) const = 0;

  /**
   * @brief Convert a state-space increment to control-variable space.
   */
  virtual void convertStateIncrementToControl(const IncrementType& increment,
                                              ControlVector& control) const = 0;

  /**
   * @brief Convert control-variable data to a state-space increment.
   */
  virtual void convertControlToStateIncrement(
      const ControlVector& control, IncrementType& increment) const = 0;

  /**
   * @brief Convert a state-space gradient to control-variable space.
   */
  virtual void convertStateGradientToControl(
      const IncrementType& state_gradient,
      ControlVector& control_gradient) const = 0;

  /**
   * @brief Convert a control-variable gradient to state-space.
   */
  virtual void convertControlGradientToState(
      const ControlVector& control_gradient,
      IncrementType& state_gradient) const = 0;
};

}  // namespace metada::framework
