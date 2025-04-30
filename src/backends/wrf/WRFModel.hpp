/**
 * @file WRFModel.hpp
 * @brief WRF model backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <memory>
#include <netcdf>
#include <string>
#include <unordered_map>
#include <vector>

#include "DateTime.hpp"
#include "Duration.hpp"

namespace metada::backends::wrf {

// Forward declaration
template <typename ConfigBackend>
class WRFState;

using core::DateTime;
using core::Duration;
/**
 * @brief WRF model backend implementation
 *
 * @details
 * This class implements a model backend for the WRF (Weather Research and
 * Forecasting) model. It provides methods for time stepping and integration of
 * the WRF dynamical core.
 */
template <typename ConfigBackend>
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
   * @param config Configuration containing WRF model parameters
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
   * @param config Configuration containing model parameters
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
   * @brief Get a model parameter
   *
   * @param name Parameter name
   * @return Parameter value as string
   * @throws std::out_of_range If parameter doesn't exist
   */
  std::string getParameter(const std::string& name) const;

  /**
   * @brief Set a model parameter
   *
   * @param name Parameter name
   * @param value Parameter value
   * @throws std::out_of_range If parameter doesn't exist or cannot be set
   */
  void setParameter(const std::string& name, const std::string& value);

  /**
   * @brief Run the model from start time to end time
   *
   * @param initialState Initial state of the model
   * @param finalState Final state after model integration (output)
   * @param startTime Start time
   * @param endTime End time
   * @throws std::runtime_error If model run fails
   */
  void run(const WRFState<ConfigBackend>& initialState,
           WRFState<ConfigBackend>& finalState, DateTime startTime,
           DateTime endTime);

  /**
   * @brief Check if the model is initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

 private:
  /**
   * @brief Run a single time step of the model
   *
   * @param inState Input state
   * @param outState Output state after time step
   * @param dt Time step size
   */
  void timeStep(const WRFState<ConfigBackend>& inState,
                WRFState<ConfigBackend>& outState, Duration dt);

  /**
   * @brief Apply boundary conditions to the model state
   *
   * @param state State to apply boundary conditions to
   */
  void applyBoundaryConditions(WRFState<ConfigBackend>& state) const;

  // Model configuration
  bool initialized_ = false;
  std::unordered_map<std::string, std::string> parameters_;

  // Time control
  DateTime startTime_;
  DateTime currentTime_;
  Duration timeStep_;
  Duration forecast_length_;

  // Physics options
  bool enableMicrophysics_ = true;
  bool enableRadiation_ = true;
  bool enablePBL_ = true;
  bool enableLSM_ = true;

  // Dynamical core options
  std::string advectionScheme_ = "WENO";
  float diffusionCoefficient_ = 0.0;
};

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
WRFModel<ConfigBackend>::WRFModel(const ConfigBackend& config)
    : initialized_(false) {
  // Initialize the model with the provided config
  initialize(config);
}

// Move constructor implementation
template <typename ConfigBackend>
WRFModel<ConfigBackend>::WRFModel(WRFModel<ConfigBackend>&& other) noexcept
    : initialized_(other.initialized_),
      parameters_(std::move(other.parameters_)),
      currentTime_(other.currentTime_),
      timeStep_(other.timeStep_),
      enableMicrophysics_(other.enableMicrophysics_),
      enableRadiation_(other.enableRadiation_),
      enablePBL_(other.enablePBL_),
      enableLSM_(other.enableLSM_),
      advectionScheme_(std::move(other.advectionScheme_)),
      diffusionCoefficient_(other.diffusionCoefficient_) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.parameters_.clear();
  other.currentTime_ = startTime_;
}

// Move assignment operator implementation
template <typename ConfigBackend>
WRFModel<ConfigBackend>& WRFModel<ConfigBackend>::operator=(
    WRFModel<ConfigBackend>&& other) noexcept {
  if (this != &other) {
    initialized_ = other.initialized_;
    parameters_ = std::move(other.parameters_);
    currentTime_ = other.currentTime_;
    timeStep_ = other.timeStep_;
    enableMicrophysics_ = other.enableMicrophysics_;
    enableRadiation_ = other.enableRadiation_;
    enablePBL_ = other.enablePBL_;
    enableLSM_ = other.enableLSM_;
    advectionScheme_ = std::move(other.advectionScheme_);
    diffusionCoefficient_ = other.diffusionCoefficient_;

    // Reset the moved-from object
    other.initialized_ = false;
    other.parameters_.clear();
    other.currentTime_ = startTime_;
  }
  return *this;
}

// Destructor implementation
template <typename ConfigBackend>
WRFModel<ConfigBackend>::~WRFModel() {
  if (initialized_) {
    try {
      finalize();
    } catch (const std::exception& e) {
      std::cerr << "Error during WRF model finalization: " << e.what()
                << std::endl;
    }
  }
}

// Initialize implementation with ConfigBackend
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::initialize(const ConfigBackend& config) {
  if (initialized_) {
    return;  // Already initialized
  }

  try {
    // Read model configuration from the config backend
    startTime_ = DateTime(config.Get("start_datetime").asString());
    forecast_length_ = Duration(config.Get("forecast_length").asString());
    timeStep_ = Duration(config.Get("time_step").asString());

    // Read physics options
    enableMicrophysics_ = config.Get("physics.microphysics").asBool();
    enableRadiation_ = config.Get("physics.radiation").asBool();
    enablePBL_ = config.Get("physics.pbl").asBool();
    enableLSM_ = config.Get("physics.lsm").asBool();

    // Read dynamics options
    advectionScheme_ = config.Get("dynamics.advection_scheme").asString();
    diffusionCoefficient_ = config.Get("dynamics.diffusion").asFloat();

    // Store all parameters for later access
    parameters_["time_step"] = timeStep_.toString();
    parameters_["enable_microphysics"] = enableMicrophysics_ ? "true" : "false";
    parameters_["enable_radiation"] = enableRadiation_ ? "true" : "false";
    parameters_["enable_pbl"] = enablePBL_ ? "true" : "false";
    parameters_["enable_lsm"] = enableLSM_ ? "true" : "false";
    parameters_["advection_scheme"] = advectionScheme_;
    parameters_["diffusion_coefficient"] =
        std::to_string(diffusionCoefficient_);

    // Additional parameters from config
    try {
      auto paramNames = config.Get("parameters").asVectorString();
      for (const auto& name : paramNames) {
        try {
          parameters_[name] = config.Get("parameters." + name).asString();
        } catch (const std::exception& e) {
          std::cerr << "Warning: Failed to read parameter '" << name
                    << "': " << e.what() << std::endl;
        }
      }
    } catch (const std::exception&) {
      // No additional parameters specified, continue
      std::cout << "No additional parameters specified, continuing..."
                << std::endl;
    }

    // Initialization successful
    initialized_ = true;
    currentTime_ = startTime_;

    std::cout << "WRF model initialized with time step: " << timeStep_
              << std::endl;

  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("WRF model initialization failed: ") +
                             e.what());
  }
}

// Reset implementation
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::reset() {
  if (!initialized_) {
    throw std::runtime_error("Cannot reset uninitialized WRF model");
  }

  // Reset model state to initial conditions
  currentTime_ = startTime_;

  std::cout << "WRF model reset to initial state" << std::endl;
}

// Finalize implementation
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::finalize() {
  if (!initialized_) {
    return;  // Nothing to finalize
  }

  // Release any resources
  parameters_.clear();
  initialized_ = false;

  std::cout << "WRF model finalized" << std::endl;
}

// Get parameter implementation
template <typename ConfigBackend>
std::string WRFModel<ConfigBackend>::getParameter(
    const std::string& name) const {
  try {
    return parameters_.at(name);
  } catch (const std::out_of_range&) {
    throw std::out_of_range("Parameter not found: " + name);
  }
}

// Set parameter implementation
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::setParameter(const std::string& name,
                                           const std::string& value) {
  // Special handling for certain parameters that affect model behavior
  if (name == "time_step") {
    try {
      timeStep_ = Duration(value);
    } catch (const std::exception&) {
      throw std::runtime_error("Invalid time step value: " + value);
    }
  } else if (name == "enable_microphysics") {
    enableMicrophysics_ = (value == "true" || value == "1");
  } else if (name == "enable_radiation") {
    enableRadiation_ = (value == "true" || value == "1");
  } else if (name == "enable_pbl") {
    enablePBL_ = (value == "true" || value == "1");
  } else if (name == "enable_lsm") {
    enableLSM_ = (value == "true" || value == "1");
  } else if (name == "advection_scheme") {
    if (value != "WENO" && value != "RK3" && value != "FD4") {
      throw std::runtime_error("Invalid advection scheme: " + value);
    }
    advectionScheme_ = value;
  } else if (name == "diffusion_coefficient") {
    try {
      diffusionCoefficient_ = std::stod(value);
    } catch (const std::exception&) {
      throw std::runtime_error("Invalid diffusion coefficient: " + value);
    }
  }

  // Store the parameter value
  parameters_[name] = value;
}

// Run implementation
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::run(const WRFState<ConfigBackend>& initialState,
                                  WRFState<ConfigBackend>& finalState,
                                  DateTime startTime, DateTime endTime) {
  if (!initialized_) {
    throw std::runtime_error("Cannot run uninitialized WRF model");
  }

  if (endTime <= startTime) {
    throw std::runtime_error("End time must be greater than start time");
  }

  try {
    // Clone the initial state to use as working state
    auto workingState = initialState.clone();

    // Set current time to start time
    currentTime_ = startTime;

    // Calculate the total simulation time
    Duration totalSimTime = endTime - startTime;

    // Use a temporary state for the time stepping process
    auto tempState = initialState.clone();

    std::cout << "Starting WRF model integration from t=" << startTime
              << " to t=" << endTime << " seconds" << std::endl;

    // Main time stepping loop
    while (currentTime_ < endTime) {
      // Ensure we don't overshoot the end time
      if (currentTime_ + timeStep_ > endTime) {
        timeStep_ = endTime - currentTime_;
      }

      // Perform a single time step
      timeStep(*workingState, *tempState, timeStep_);

      // Swap states for next iteration
      std::swap(workingState, tempState);

      // Update current time
      currentTime_ += timeStep_;

      // Optional: print progress
      // Check if we're near an hourly mark (report every simulated hour)
      Duration hoursSinceStart = currentTime_ - startTime;
      if (hoursSinceStart.totalSeconds() % 3600 < timeStep_.totalSeconds()) {
        // Calculate progress as percentage
        double elapsedSeconds = (currentTime_ - startTime).totalSeconds();
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
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::timeStep(const WRFState<ConfigBackend>& inState,
                                       WRFState<ConfigBackend>& outState,
                                       [[maybe_unused]] Duration dt) {
  // Clone input state to output state first
  outState = std::move(*inState.clone());

  try {
    // Get variable names from the state
    const auto& varNames = inState.getVariableNames();

    // Apply model physics to each variable
    // This is a simplified implementation - a real WRF model would have
    // complex physics and dynamics calculations here

    // 1. Apply advection
    if (advectionScheme_ == "WENO") {
      // WENO advection scheme (simplified)
      for (const auto& varName : varNames) {
        outState.setActiveVariable(varName);
        // Perform advection calculations here
      }
    } else if (advectionScheme_ == "RK3") {
      // Runge-Kutta 3 advection scheme (simplified)
      for (const auto& varName : varNames) {
        outState.setActiveVariable(varName);
        // Perform advection calculations here
      }
    } else {
      // Default FD4 advection scheme (simplified)
      for (const auto& varName : varNames) {
        outState.setActiveVariable(varName);
        // Perform advection calculations here
      }
    }

    // 2. Apply physics parameterizations

    // Microphysics (if enabled)
    if (enableMicrophysics_) {
      // Apply microphysics to relevant variables
      if (std::find(varNames.begin(), varNames.end(), "QVAPOR") !=
          varNames.end()) {
        outState.setActiveVariable("QVAPOR");
        // Apply microphysics to water vapor
      }
    }

    // Radiation (if enabled)
    if (enableRadiation_) {
      // Apply radiation physics
      if (std::find(varNames.begin(), varNames.end(), "T") != varNames.end()) {
        outState.setActiveVariable("T");
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
      for (const auto& varName : varNames) {
        outState.setActiveVariable(varName);
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
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::applyBoundaryConditions(
    WRFState<ConfigBackend>& state) const {
  // Apply appropriate boundary conditions to each variable
  // For demonstration purposes, we'll just apply simple conditions

  const auto& varNames = state.getVariableNames();

  for (const auto& varName : varNames) {
    state.setActiveVariable(varName);

    // In a real implementation, this would apply periodic, open, or fixed
    // boundary conditions depending on the variable and domain configuration
  }
}

}  // namespace metada::backends::wrf