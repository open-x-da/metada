/**
 * @file TimeStepperImpl.hpp
 * @brief Implementation of the time stepping capability for models
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains the implementation of the ITimeStepper interface
 * that delegates to a backend model implementation. It provides time stepping
 * functionality for physical numerical models.
 */

#pragma once

#include "../../interfaces/ITimeStepper.hpp"

namespace metada::framework {

/**
 * @brief Implementation of the time stepping capability
 *
 * @details
 * This class implements the ITimeStepper interface by delegating
 * to the backend model implementation. It provides time stepping
 * functionality for physical numerical models.
 *
 * @tparam Backend The backend type implementing the model
 */
template <typename Backend>
class TimeStepperImpl : public ITimeStepper {
 private:
  Backend& backend_;  ///< Reference to the backend model

 public:
  /**
   * @brief Construct a new TimeStepperImpl
   *
   * @param backend Reference to the backend model
   */
  explicit TimeStepperImpl(Backend& backend) : backend_(backend) {}

  /**
   * @brief Step the model forward by one time step
   *
   * @param currentState The current state of the model
   * @param nextState The next state after one time step (output parameter)
   */
  void step(const IState& currentState, IState& nextState) override {
    backend_.step(currentState, nextState);
  }

  /**
   * @brief Get the model's time step size
   *
   * @return The time step size in seconds
   */
  double getTimeStep() const override { return backend_.getTimeStep(); }

  /**
   * @brief Set the model's time step size
   *
   * @param dt The new time step size in seconds
   */
  void setTimeStep(double dt) override { backend_.setTimeStep(dt); }

  /**
   * @brief Get the maximum stable time step for the current model state
   *
   * @param state The current model state
   * @return The maximum stable time step in seconds
   */
  double getMaxStableTimeStep(const IState& state) override {
    return backend_.getMaxStableTimeStep(state);
  }

  /**
   * @brief Check if the model uses an adaptive time step
   *
   * @return true if the model uses an adaptive time step, false otherwise
   */
  bool usesAdaptiveTimeStep() const override {
    return backend_.usesAdaptiveTimeStep();
  }

  /**
   * @brief Enable or disable adaptive time stepping
   *
   * @param enable true to enable adaptive time stepping, false to disable
   */
  void setAdaptiveTimeStep(bool enable) override {
    backend_.setAdaptiveTimeStep(enable);
  }

  /**
   * @brief Get the current model time
   *
   * @return The current model time in seconds
   */
  double getCurrentTime() const override { return backend_.getCurrentTime(); }

  /**
   * @brief Set the current model time
   *
   * @param time The new model time in seconds
   */
  void setCurrentTime(double time) override { backend_.setCurrentTime(time); }

  /**
   * @brief Get the time stepping scheme used by the model
   *
   * @return The name of the time stepping scheme
   */
  std::string getTimeSteppingScheme() const override {
    return backend_.getTimeSteppingScheme();
  }
};

}  // namespace metada::framework