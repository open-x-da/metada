/**
 * @file ITimeStepper.hpp
 * @brief Interface defining the capability for time-stepping models
 * @ingroup repr
 *
 * @details
 * This header provides the capability interface for models that support
 * explicit time stepping, typically used by traditional numerical models.
 * Models that implement this interface can be stepped forward in time
 * with explicit control over the time step.
 */

#pragma once

#include <string>

namespace metada::framework {

// Forward declarations
class IState;

/**
 * @brief Capability interface for time-stepping models
 *
 * @details
 * This interface defines the capabilities required for models that
 * support explicit time stepping. Models that implement this capability
 * can be advanced in time with explicit control over the time step.
 *
 * This capability is typically used by:
 * - Traditional numerical weather prediction models
 * - Fluid dynamics models
 * - Climate models
 * - Any model with explicit time integration
 *
 * Implementations should ensure:
 * - Consistent time step management
 * - Stability of the time stepping scheme
 * - Appropriate handling of time-dependent processes
 */
class ITimeStepper {
 public:
  /**
   * @brief Virtual destructor
   */
  virtual ~ITimeStepper() = default;

  /**
   * @brief Step the model forward by one time step
   *
   * @param currentState The current state of the model
   * @param nextState The next state after one time step (output parameter)
   * @throws std::runtime_error if the model step fails
   */
  virtual void step(const IState& currentState, IState& nextState) = 0;

  /**
   * @brief Get the model's time step size
   *
   * @return The time step size in seconds
   */
  virtual double getTimeStep() const = 0;

  /**
   * @brief Set the model's time step size
   *
   * @param dt The new time step size in seconds
   * @throws std::invalid_argument if the time step is invalid
   */
  virtual void setTimeStep(double dt) = 0;

  /**
   * @brief Get the maximum stable time step for the current model state
   *
   * @param state The current model state
   * @return The maximum stable time step in seconds
   */
  virtual double getMaxStableTimeStep(const IState& state) = 0;

  /**
   * @brief Check if the model uses an adaptive time step
   *
   * @return true if the model uses an adaptive time step, false otherwise
   */
  virtual bool usesAdaptiveTimeStep() const = 0;

  /**
   * @brief Enable or disable adaptive time stepping
   *
   * @param enable true to enable adaptive time stepping, false to disable
   */
  virtual void setAdaptiveTimeStep(bool enable) = 0;

  /**
   * @brief Get the current model time
   *
   * @return The current model time in seconds
   */
  virtual double getCurrentTime() const = 0;

  /**
   * @brief Set the current model time
   *
   * @param time The new model time in seconds
   */
  virtual void setCurrentTime(double time) = 0;

  /**
   * @brief Get the time stepping scheme used by the model
   *
   * @return The name of the time stepping scheme (e.g., "RK4", "Leapfrog",
   * etc.)
   */
  virtual std::string getTimeSteppingScheme() const = 0;
};

}  // namespace metada::framework