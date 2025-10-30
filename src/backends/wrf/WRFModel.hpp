/**
 * @file WRFModel.hpp
 * @brief WRF model backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <netcdf>
#include <string>

#include "DateTime.hpp"
#include "Duration.hpp"

namespace metada::backends::wrf {

// Forward declaration
template <typename ConfigBackend, typename GeometryBackend>
class WRFState;

/**
 * @brief WRF model backend implementation
 *
 * @details
 * This class implements a model backend for the WRF (Weather Research and
 * Forecasting) model. It provides methods for time stepping and integration of
 * the WRF dynamical core.
 */
template <typename ConfigBackend, typename StateBackend>
class WRFModel {
 public:
  /**
   * @brief Default constructor is deleted
   */
  WRFModel() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  WRFModel(const WRFModel&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  WRFModel& operator=(const WRFModel&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing WRF model options
   */
  explicit WRFModel(const ConfigBackend& config);

  /**
   * @brief Move constructor
   *
   * @param other WRF model backend to move from
   */
  WRFModel(WRFModel&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other WRF model backend to move from
   * @return Reference to this model after assignment
   */
  WRFModel& operator=(WRFModel&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~WRFModel();

  /**
   * @brief Initialize the model with a configuration
   *
   * @param config Configuration containing model options
   */
  void initialize(const ConfigBackend& config);

  /**
   * @brief Reset the model to its initial state
   */
  void reset();

  /**
   * @brief Finalize the model, releasing resources
   */
  void finalize();

  /**
   * @brief Run the model from start time to end time
   *
   * @param initialState Initial state of the model
   * @param finalState Final state after model integration (output)
   * @throws std::runtime_error If model run fails
   */
  void run(const StateBackend& initialState, StateBackend& finalState) const;

  /**
   * @brief Check if the model is initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Check if the model supports adjoint computations
   *
   * @return True if adjoint is supported, false otherwise
   */
  bool supportsAdjoint() const {
    return false;
  }  // TODO: Implement adjoint for WRF

  /**
   * @brief Run adjoint model integration
   *
   * @details This is a placeholder implementation for the adjoint model.
   * In a full implementation, this would run the adjoint of the forward model
   * to compute gradients for variational data assimilation.
   *
   * @param initial_state Initial state for adjoint integration
   * @param final_state Final state for adjoint integration
   * @param forcing_increment Forcing increment for adjoint
   * @param result_increment Resulting increment from adjoint integration
   */
  template <typename StateOrIncrement>
  void runAdjoint([[maybe_unused]] const StateBackend& initial_state,
                  [[maybe_unused]] const StateBackend& final_state,
                  [[maybe_unused]] const StateOrIncrement& forcing_increment,
                  [[maybe_unused]] StateOrIncrement& result_increment) const {
    // Placeholder implementation - in practice this would be a complex
    // adjoint model integration

    // For now, just copy the forcing increment to the result
    // This is obviously not a correct adjoint implementation, but allows
    // the code to compile and run for testing purposes
    if constexpr (std::is_same_v<StateOrIncrement, StateBackend>) {
      result_increment = std::move(*(forcing_increment.clone()));
    } else {
      // For increment, use axpy to copy
      result_increment.zero();
      result_increment.axpy(1.0, forcing_increment);
    }

    // In a real implementation, this would:
    // 1. Load the trajectory from the forward model run
    // 2. Integrate the adjoint model backward in time
    // 3. Apply the adjoint of each model component
    // 4. Accumulate the gradient information

    std::cout << "WRF adjoint model run (placeholder implementation)"
              << std::endl;
  }

 private:
  /**
   * @brief Run a single time step of the model
   *
   * @param inState Input state
   * @param outState Output state after time step
   * @param dt Time step size
   */
  void timeStep(const StateBackend& inState, StateBackend& outState,
                Duration dt) const;

  /**
   * @brief Apply boundary conditions to the model state
   *
   * @param state State to apply boundary conditions to
   */
  void applyBoundaryConditions(StateBackend& state) const;

  // Model configuration
  bool initialized_ = false;

  // Time control
  DateTime startTime_;
  mutable DateTime currentTime_;  // Allow modification in const methods
  DateTime endTime_;
  mutable Duration timeStep_;  // Allow modification in const methods
  bool output_history_;
  std::string history_file_;
  Duration history_frequency_;

  // Physics options
  bool enableMicrophysics_ = true;
  bool enableRadiation_ = true;
  bool enablePBL_ = true;
  bool enableLSM_ = true;

  // Dynamical core options
  std::string advectionScheme_ = "WENO";
  float diffusionCoefficient_ = 0.0;
};

// Constructor implementation with ConfigBackend and StateBackend
template <typename ConfigBackend, typename StateBackend>
WRFModel<ConfigBackend, StateBackend>::WRFModel(const ConfigBackend& config)
    : initialized_(false) {
  // Initialize the model with the provided config
  initialize(config);
}

// Move constructor implementation
template <typename ConfigBackend, typename StateBackend>
WRFModel<ConfigBackend, StateBackend>::WRFModel(
    WRFModel<ConfigBackend, StateBackend>&& other) noexcept
    : initialized_(other.initialized_),
      currentTime_(other.currentTime_),
      timeStep_(other.timeStep_),
      output_history_(other.output_history_),
      history_file_(other.history_file_),
      history_frequency_(other.history_frequency_),
      enableMicrophysics_(other.enableMicrophysics_),
      enableRadiation_(other.enableRadiation_),
      enablePBL_(other.enablePBL_),
      enableLSM_(other.enableLSM_),
      advectionScheme_(std::move(other.advectionScheme_)),
      diffusionCoefficient_(other.diffusionCoefficient_) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.currentTime_ = startTime_;
}

// Move assignment operator implementation
template <typename ConfigBackend, typename StateBackend>
WRFModel<ConfigBackend, StateBackend>&
WRFModel<ConfigBackend, StateBackend>::operator=(
    WRFModel<ConfigBackend, StateBackend>&& other) noexcept {
  if (this != &other) {
    initialized_ = other.initialized_;
    currentTime_ = other.currentTime_;
    timeStep_ = other.timeStep_;
    output_history_ = other.output_history_;
    history_file_ = other.history_file_;
    history_frequency_ = other.history_frequency_;
    enableMicrophysics_ = other.enableMicrophysics_;
    enableRadiation_ = other.enableRadiation_;
    enablePBL_ = other.enablePBL_;
    enableLSM_ = other.enableLSM_;
    advectionScheme_ = std::move(other.advectionScheme_);
    diffusionCoefficient_ = other.diffusionCoefficient_;

    // Reset the moved-from object
    other.initialized_ = false;
    other.currentTime_ = startTime_;
  }
  return *this;
}

// Destructor implementation
template <typename ConfigBackend, typename StateBackend>
WRFModel<ConfigBackend, StateBackend>::~WRFModel() {
  if (initialized_) {
    try {
      finalize();
    } catch (const std::exception& e) {
      std::cerr << "Error during WRF model finalization: " << e.what()
                << std::endl;
    }
  }
}

/**
 * @brief Initialize the model with a configuration
 *
 * @param config Configuration containing model options
 */
template <typename ConfigBackend, typename StateBackend>
void WRFModel<ConfigBackend, StateBackend>::initialize(
    const ConfigBackend& config) {
  if (initialized_) {
    return;  // Already initialized
  }

  try {
    // Time control options
    auto time_control = config.CreateSubsection("time_control");
    startTime_ = DateTime(time_control.Get("start_datetime").asString());

    // Handle forecast length and end time calculations
    Duration forecast_length;
    if (time_control.HasKey("forecast_length")) {
      // Use forecast_length to determine end time
      forecast_length =
          Duration(time_control.Get("forecast_length").asString());
      endTime_ = startTime_ + forecast_length;
    } else if (time_control.HasKey("end_datetime")) {
      // Use explicit end_datetime
      endTime_ = DateTime(time_control.Get("end_datetime").asString());
      forecast_length = endTime_ - startTime_;
    } else {
      throw std::runtime_error(
          "Either forecast_length or end_datetime must be specified");
    }

    // Time step configuration
    timeStep_ = Duration(time_control.Get("time_step").asString());

    // Output history settings
    if (time_control.HasKey("output_history")) {
      output_history_ = time_control.Get("output_history").asBool();
    }

    // If output history is enabled, set the history file and frequency
    if (output_history_) {
      history_file_ = time_control.HasKey("history_file")
                          ? time_control.Get("history_file").asString()
                          : "wrfout";

      history_frequency_ =
          time_control.HasKey("history_frequency")
              ? Duration(time_control.Get("history_frequency").asString())
              : Duration("360s");
    }

    // Physics options
    auto physics = config.CreateSubsection("physics");
    enableMicrophysics_ = physics.Get("microphysics").asBool();
    enableRadiation_ = physics.Get("radiation").asBool();
    enablePBL_ = physics.Get("pbl").asBool();
    enableLSM_ = physics.Get("lsm").asBool();

    // Dynamical core options
    auto dynamics = config.CreateSubsection("dynamics");
    advectionScheme_ = dynamics.Get("advection_scheme").asString();
    diffusionCoefficient_ = dynamics.Get("diffusion").asFloat();

    // Initialization successful
    initialized_ = true;
    currentTime_ = startTime_;

    std::cout << "WRF model initialized with time step: " << timeStep_
              << std::endl;
    std::cout << "Forecast period: " << startTime_ << " to " << endTime_
              << " (duration: " << forecast_length << ")" << std::endl;

  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("WRF model initialization failed: ") +
                             e.what());
  }
}

// Reset implementation
template <typename ConfigBackend, typename StateBackend>
void WRFModel<ConfigBackend, StateBackend>::reset() {
  if (!initialized_) {
    throw std::runtime_error("Cannot reset uninitialized WRF model");
  }

  // Reset model state to initial conditions
  currentTime_ = startTime_;

  std::cout << "WRF model reset to initial state" << std::endl;
}

// Finalize implementation
template <typename ConfigBackend, typename StateBackend>
void WRFModel<ConfigBackend, StateBackend>::finalize() {
  if (!initialized_) {
    return;  // Nothing to finalize
  }

  // Release any resources
  initialized_ = false;

  std::cout << "WRF model finalized" << std::endl;
}

// Run implementation
template <typename ConfigBackend, typename StateBackend>
void WRFModel<ConfigBackend, StateBackend>::run(
    const StateBackend& initialState, StateBackend& finalState) const {
  if (!initialized_) {
    throw std::runtime_error("Cannot run uninitialized WRF model");
  }

  if (endTime_ <= startTime_) {
    throw std::runtime_error("End time must be greater than start time");
  }

  try {
    // Clone the initial state to use as working state
    auto workingState = initialState.clone();

    // Set current time to start time
    currentTime_ = startTime_;

    // Calculate the total simulation time
    Duration totalSimTime = endTime_ - startTime_;

    // Use a temporary state for the time stepping process
    auto tempState = initialState.clone();

    std::cout << "Starting WRF model integration from t=" << startTime_
              << " to t=" << endTime_ << std::endl;

    // Main time stepping loop
    while (currentTime_ < endTime_) {
      // Ensure we don't overshoot the end time
      if (currentTime_ + timeStep_ > endTime_) {
        timeStep_ = endTime_ - currentTime_;
      }

      // Perform a single time step
      timeStep(*workingState, *tempState, timeStep_);

      // Swap states for next iteration
      std::swap(workingState, tempState);

      // Update current time
      currentTime_ += timeStep_;

      // Optional: print progress
      // Check if we're near an hourly mark (report every simulated hour)
      Duration hoursSinceStart = currentTime_ - startTime_;
      if (hoursSinceStart.totalSeconds() % 3600 < timeStep_.totalSeconds()) {
        // Calculate progress as percentage
        double elapsedSeconds = (currentTime_ - startTime_).totalSeconds();
        double totalSeconds = totalSimTime.totalSeconds();
        double progress = (elapsedSeconds / totalSeconds) * 100.0;

        std::cout << "WRF model integration: " << std::fixed
                  << std::setprecision(1) << progress
                  << "% complete (t=" << currentTime_ << ")" << std::endl;
      }
    }

    // Copy the final state
    finalState = std::move(*workingState);

    std::cout << "WRF model integration completed successfully" << std::endl;

  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("WRF model run failed: ") + e.what());
  }
}

// Time step implementation
template <typename ConfigBackend, typename StateBackend>
void WRFModel<ConfigBackend, StateBackend>::timeStep(
    const StateBackend& inState, StateBackend& outState,
    [[maybe_unused]] Duration dt) const {
  // Clone input state to output state first
  outState = std::move(*(inState.clone()));

  try {
    // Get variable names from the state
    const auto& varNames = inState.getVariableNames();

    // Apply model physics to each variable
    // This is a simplified implementation - a real WRF model would have
    // complex physics and dynamics calculations here

    // 1. Apply advection
    if (advectionScheme_ == "WENO") {
      // WENO advection scheme (simplified)
      for ([[maybe_unused]] const auto& varName : varNames) {
        // TODO: Implement setActiveVariable or refactor to work without it
        // outState.setActiveVariable(varName);
        // Perform advection calculations here
      }
    } else if (advectionScheme_ == "RK3") {
      // Runge-Kutta 3 advection scheme (simplified)
      for ([[maybe_unused]] const auto& varName : varNames) {
        // TODO: Implement setActiveVariable or refactor to work without it
        // outState.setActiveVariable(varName);
        // Perform advection calculations here
      }
    } else {
      // Default FD4 advection scheme (simplified)
      for ([[maybe_unused]] const auto& varName : varNames) {
        // TODO: Implement setActiveVariable or refactor to work without it
        // outState.setActiveVariable(varName);
        // Perform advection calculations here
      }
    }

    // 2. Apply physics parameterizations

    // Microphysics (if enabled)
    if (enableMicrophysics_) {
      // Apply microphysics to relevant variables
      if (std::find(varNames.begin(), varNames.end(), "QVAPOR") !=
          varNames.end()) {
        // TODO: Implement setActiveVariable or refactor to work without it
        // outState.setActiveVariable("QVAPOR");
        // Apply microphysics to water vapor
      }
    }

    // Radiation (if enabled)
    if (enableRadiation_) {
      // Apply radiation physics
      if (std::find(varNames.begin(), varNames.end(), "T") != varNames.end()) {
        // TODO: Implement setActiveVariable or refactor to work without it
        // outState.setActiveVariable("T");
        // Apply radiation effects to temperature
      }
    }

    // Planetary Boundary Layer (if enabled)
    if (enablePBL_) {
      // Apply PBL physics
      // Would update wind, temperature, and moisture in the boundary layer
    }

    // Land Surface Model (if enabled)
    if (enableLSM_) {
      // Apply land-surface interactions
      // Would update surface fluxes, soil temperature, and moisture
    }

    // 3. Apply diffusion (if coefficient > 0)
    if (diffusionCoefficient_ > 0.0) {
      for ([[maybe_unused]] const auto& varName : varNames) {
        // TODO: Implement setActiveVariable or refactor to work without it
        // outState.setActiveVariable(varName);
        // Apply diffusion to variable
      }
    }

    // 4. Apply boundary conditions
    applyBoundaryConditions(outState);

  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Error in time step calculation: ") +
                             e.what());
  }
}

// Apply boundary conditions implementation
template <typename ConfigBackend, typename StateBackend>
void WRFModel<ConfigBackend, StateBackend>::applyBoundaryConditions(
    StateBackend& state) const {
  // Apply appropriate boundary conditions to each variable
  // For demonstration purposes, we'll just apply simple conditions

  const auto& varNames = state.getVariableNames();

  for ([[maybe_unused]] const auto& varName : varNames) {
    // TODO: Implement setActiveVariable or refactor to work without it
    // state.setActiveVariable(varName);

    // In a real implementation, this would apply periodic, open, or fixed
    // boundary conditions depending on the variable and domain configuration
  }
}

}  // namespace metada::backends::wrf