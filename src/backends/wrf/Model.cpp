/**
 * @file Model.cpp
 * @brief Implementation of the WRF model backend
 * @ingroup backends
 * @author Metada Framework Team
 */

#include "Model.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "State.hpp"

namespace metada::backends::wrf {

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
Model::Model(const ConfigBackend& config) : initialized_(false) {
  // Initialize the model with the provided config
  initialize(config);
}

// Move constructor implementation
Model::Model(Model&& other) noexcept
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
Model& Model::operator=(Model&& other) noexcept {
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
Model::~Model() {
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
void Model::initialize(const ConfigBackend& config) {
  if (initialized_) {
    return;  // Already initialized
  }

  try {
    // Read model configuration from the config backend
    timeStep_ = config.getDouble("wrf.time_step", 60.0);  // Default: 60 seconds

    // Read physics options
    enableMicrophysics_ = config.getBool("wrf.physics.microphysics", true);
    enableRadiation_ = config.getBool("wrf.physics.radiation", true);
    enablePBL_ = config.getBool("wrf.physics.pbl", true);
    enableLSM_ = config.getBool("wrf.physics.lsm", true);

    // Read dynamics options
    advectionScheme_ =
        config.getString("wrf.dynamics.advection_scheme", "WENO");
    diffusionCoefficient_ = config.getDouble("wrf.dynamics.diffusion", 0.0);
    cflNumber_ = config.getDouble("wrf.dynamics.cfl", 0.9);

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
      auto paramNames = config.getStringVector("wrf.parameters");
      for (const auto& name : paramNames) {
        try {
          parameters_[name] = config.getString("wrf.parameters." + name, "");
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
void Model::reset() {
  if (!initialized_) {
    throw std::runtime_error("Cannot reset uninitialized WRF model");
  }

  // Reset model state to initial conditions
  currentTime_ = 0.0;

  std::cout << "WRF model reset to initial state" << std::endl;
}

// Finalize implementation
void Model::finalize() {
  if (!initialized_) {
    return;  // Nothing to finalize
  }

  // Release any resources
  parameters_.clear();
  initialized_ = false;

  std::cout << "WRF model finalized" << std::endl;
}

// Get parameter implementation
std::string Model::getParameter(const std::string& name) const {
  try {
    return parameters_.at(name);
  } catch (const std::out_of_range&) {
    throw std::out_of_range("Parameter not found: " + name);
  }
}

// Set parameter implementation
void Model::setParameter(const std::string& name, const std::string& value) {
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
void Model::run(const WRFStateBackend& initialState,
                WRFStateBackend& finalState, double startTime, double endTime) {
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
void Model::timeStep(const WRFStateBackend& inState, WRFStateBackend& outState,
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
double Model::calculateTimeStep(const WRFStateBackend& state,
                                double maxDt) const {
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
void Model::applyBoundaryConditions(WRFStateBackend& state) const {
  // Apply appropriate boundary conditions to each variable
  // For demonstration purposes, we'll just apply simple conditions

  const auto& varNames = state.getVariableNames();

  for (const auto& varName : varNames) {
    state.setActiveVariable(varName);

    // In a real implementation, this would apply periodic, open, or fixed
    // boundary conditions depending on the variable and domain configuration
  }
}

// Explicit template instantiations for known config backend types
template Model::Model(const metada::backends::yaml::YAMLConfigBackend& config);
template void Model::initialize(
    const metada::backends::yaml::YAMLConfigBackend& config);

}  // namespace metada::backends::wrf