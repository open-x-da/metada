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

namespace metada::backends::wrf {

// Forward declaration
template <typename ConfigBackend>
class WRFState;

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
   * @param startTime Start time in seconds
   * @param endTime End time in seconds
   * @throws std::runtime_error If model run fails
   */
  void run(const WRFState<ConfigBackend>& initialState,
           WRFState<ConfigBackend>& finalState, double startTime,
           double endTime);

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
   * @param dt Time step size in seconds
   */
  void timeStep(const WRFState<ConfigBackend>& inState,
                WRFState<ConfigBackend>& outState, double dt);

  /**
   * @brief Calculate adaptive time step based on CFL condition
   *
   * @param state Current model state
   * @param maxDt Maximum allowed time step
   * @return Calculated time step size
   */
  double calculateTimeStep(const WRFState<ConfigBackend>& state,
                           double maxDt) const;

  /**
   * @brief Apply boundary conditions to the model state
   *
   * @param state State to apply boundary conditions to
   */
  void applyBoundaryConditions(WRFState<ConfigBackend>& state) const;

  // Model configuration
  bool initialized_ = false;
  std::unordered_map<std::string, std::string> parameters_;

  // Model state
  double currentTime_ = 0.0;
  double timeStep_ = 0.0;

  // Physics options
  bool enableMicrophysics_ = true;
  bool enableRadiation_ = true;
  bool enablePBL_ = true;
  bool enableLSM_ = true;

  // Dynamical core options
  std::string advectionScheme_ = "WENO";
  double diffusionCoefficient_ = 0.0;
  double cflNumber_ = 0.9;
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
      diffusionCoefficient_(other.diffusionCoefficient_),
      cflNumber_(other.cflNumber_) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.parameters_.clear();
  other.currentTime_ = 0.0;
  other.timeStep_ = 0.0;
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
    cflNumber_ = other.cflNumber_;

    // Reset the moved-from object
    other.initialized_ = false;
    other.parameters_.clear();
    other.currentTime_ = 0.0;
    other.timeStep_ = 0.0;
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
    timeStep_ = config.Get("time_step").asFloat();  // Default: 60 seconds

    // Read physics options
    enableMicrophysics_ = config.Get("physics.microphysics").asBool();
    enableRadiation_ = config.Get("physics.radiation").asBool();
    enablePBL_ = config.Get("physics.pbl").asBool();
    enableLSM_ = config.Get("physics.lsm").asBool();

    // Read dynamics options
    advectionScheme_ = config.Get("dynamics.advection_scheme").asString();
    diffusionCoefficient_ = config.Get("dynamics.diffusion").asFloat();
    cflNumber_ = config.Get("dynamics.cfl").asFloat();

    // Store all parameters for later access
    parameters_["time_step"] = std::to_string(timeStep_);
    parameters_["enable_microphysics"] = enableMicrophysics_ ? "true" : "false";
    parameters_["enable_radiation"] = enableRadiation_ ? "true" : "false";
    parameters_["enable_pbl"] = enablePBL_ ? "true" : "false";
    parameters_["enable_lsm"] = enableLSM_ ? "true" : "false";
    parameters_["advection_scheme"] = advectionScheme_;
    parameters_["diffusion_coefficient"] =
        std::to_string(diffusionCoefficient_);
    parameters_["cfl_number"] = std::to_string(cflNumber_);

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
    }

    // Initialization successful
    initialized_ = true;
    currentTime_ = 0.0;

    std::cout << "WRF model initialized with time step: " << timeStep_
              << " seconds" << std::endl;

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
  currentTime_ = 0.0;

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
      timeStep_ = std::stod(value);
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
  } else if (name == "cfl_number") {
    try {
      double cfl = std::stod(value);
      if (cfl <= 0.0 || cfl > 1.0) {
        throw std::runtime_error("CFL number must be between 0 and 1");
      }
      cflNumber_ = cfl;
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Invalid CFL number: ") + e.what());
    }
  }

  // Store the parameter value
  parameters_[name] = value;
}

// Run implementation
template <typename ConfigBackend>
void WRFModel<ConfigBackend>::run(const WRFState<ConfigBackend>& initialState,
                                  WRFState<ConfigBackend>& finalState,
                                  double startTime, double endTime) {
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
    double totalSimTime = endTime - startTime;

    // Use a temporary state for the time stepping process
    auto tempState = initialState.clone();

    std::cout << "Starting WRF model integration from t=" << startTime
              << " to t=" << endTime << " seconds" << std::endl;

    // Main time stepping loop
    while (currentTime_ < endTime) {
      // Calculate adaptive time step if needed
      double dt = timeStep_;
      if (cflNumber_ > 0.0) {
        dt = calculateTimeStep(*workingState, timeStep_);
      }

      // Ensure we don't overshoot the end time
      if (currentTime_ + dt > endTime) {
        dt = endTime - currentTime_;
      }

      // Perform a single time step
      timeStep(*workingState, *tempState, dt);

      // Swap states for next iteration
      std::swap(workingState, tempState);

      // Update current time
      currentTime_ += dt;

      // Optional: print progress
      if (std::fmod(currentTime_, 3600.0) <
          dt) {  // Report every simulated hour
        double progress = (currentTime_ - startTime) / totalSimTime * 100.0;
        std::cout << "WRF model integration: " << std::fixed
                  << std::setprecision(1) << progress
                  << "% complete (t=" << currentTime_ << "s)" << std::endl;
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
                                       double dt) {
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

// Calculate time step implementation
template <typename ConfigBackend>
double WRFModel<ConfigBackend>::calculateTimeStep(
    const WRFState<ConfigBackend>& state, double maxDt) const {
  // This would calculate an appropriate time step based on the CFL condition
  // For simplicity, we'll just return the maximum time step here

  // In a real implementation, we would:
  // 1. Find the maximum wind speed in the domain
  // 2. Find the minimum grid spacing
  // 3. Calculate dt = cflNumber_ * (minimum grid spacing / maximum wind speed)
  // 4. Return min(dt, maxDt)

  // For demonstration purposes, we'll use a simplified calculation
  try {
    // Get wind components if available
    double maxWindSpeed = 10.0;  // Default assumption in m/s

    const auto& varNames = state.getVariableNames();
    if (std::find(varNames.begin(), varNames.end(), "U") != varNames.end() &&
        std::find(varNames.begin(), varNames.end(), "V") != varNames.end()) {
      // In a real implementation, find the maximum wind speed from U and V
      // For now, we'll use a dummy calculation
      maxWindSpeed = 10.0 + 5.0 * std::sin(currentTime_ /
                                           3600.0);  // Varies between 5-15 m/s
    }

    // Assume a grid spacing of 3 km
    double minGridSpacing = 3000.0;  // meters

    // Calculate time step using CFL condition
    double dt = cflNumber_ * (minGridSpacing / maxWindSpeed);

    // Limit to maximum time step
    return std::min(dt, maxDt);

  } catch (const std::exception&) {
    // Fall back to default time step on error
    return maxDt;
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