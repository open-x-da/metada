#pragma once

#include <memory>
#include <string>
#include <vector>

namespace metada::framework {

/**
 * @brief Supported control-variable backend implementations.
 */
enum class ControlVariableBackendKind {
  Identity,  ///< Direct grid%xa control vector (current default)
  WrfdaCv5   ///< WRFDA CV5 control vector
};

// Forward declarations - these will be properly defined in adapter headers
template <typename BackendTag>
class Increment;

template <typename BackendTag>
class ControlVariable;

template <typename BackendTag>
class State;

/**
 * @brief Abstract interface for control-variable backend operations
 *
 * @details In incremental variational data assimilation, the control variable v
 * is the actual optimization variable, while the increment δx is in state
 * space. The relationship between them is:
 *
 *   δx = U v        (forward transformation)
 *   v^T = δx^T U^T  (adjoint transformation)
 *
 * where U is typically related to B^(1/2) (square root of background error
 * covariance matrix). For example:
 * - Identity backend: U = I, so v = δx (direct representation)
 * - WRFDA CV5 backend: U involves balance relationships and vertical/horizontal
 *   transforms
 *
 * This abstraction allows the incremental variational workflow to remain
 * agnostic to the chosen control-variable representation.
 *
 * @tparam BackendTag Backend tag satisfying framework concepts
 */
template <typename BackendTag>
class ControlVariableBackend {
 public:
  using IncrementType = Increment<BackendTag>;
  using ControlVariableType = ControlVariable<BackendTag>;
  using StateType = State<BackendTag>;

  ControlVariableBackend() = default;
  virtual ~ControlVariableBackend() = default;

  // Delete copy and move operations
  ControlVariableBackend(const ControlVariableBackend&) = delete;
  ControlVariableBackend& operator=(const ControlVariableBackend&) = delete;
  ControlVariableBackend(ControlVariableBackend&&) = delete;
  ControlVariableBackend& operator=(ControlVariableBackend&&) = delete;

  /**
   * @brief Human-readable backend name
   *
   * @return Name of the control variable backend (e.g., "identity",
   * "wrfda_cv5")
   */
  virtual std::string name() const = 0;

  /**
   * @brief Get the backend kind identifier
   *
   * @return Backend kind enum value
   */
  virtual ControlVariableBackendKind kind() const = 0;

  /**
   * @brief Create a control variable from geometry
   *
   * @param geometry Geometry backend defining the domain
   * @return Newly created control variable
   */
  virtual ControlVariableType createControlVariable(
      const typename ControlVariableType::GeometryBackendType& geometry)
      const = 0;

  /**
   * @brief Create an increment from geometry
   *
   * @param geometry Geometry backend defining the domain
   * @return Newly created increment
   */
  virtual IncrementType createIncrement(
      const typename IncrementType::GeometryBackendType& geometry) const = 0;

  /**
   * @brief Forward transformation: δx = U v
   *
   * @details Transform control variable to state-space increment. This is used
   * during cost function evaluation to obtain the increment from the control
   * variable.
   *
   * @param control Control variable in control space
   * @param increment Output state-space increment (will be overwritten)
   */
  virtual void controlToIncrement(const ControlVariableType& control,
                                  IncrementType& increment) const = 0;

  /**
   * @brief Adjoint transformation: ∇_v J = U^T ∇_δx J
   *
   * @details Transform state-space gradient to control-space gradient. This is
   * used during gradient computation to transform the gradient from state space
   * to control space.
   *
   * @param state_gradient Gradient in state space (∇_δx J)
   * @param control_gradient Output gradient in control space (∇_v J)
   */
  virtual void incrementAdjointToControl(
      const IncrementType& state_gradient,
      ControlVariableType& control_gradient) const = 0;

  /**
   * @brief Apply increment to state: xa = xb + δx
   *
   * @details Add a state-space increment to a state to produce the analysis
   * state. This is used at the end of the minimization to compute the final
   * analysis state from the background and the optimal increment.
   *
   * @param increment State-space increment δx
   * @param state State to be updated (input: xb, output: xa = xb + δx)
   */
  virtual void addIncrementToState(const IncrementType& increment,
                                   StateType& state) const = 0;
};

}  // namespace metada::framework
