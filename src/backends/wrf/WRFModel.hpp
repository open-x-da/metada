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

}  // namespace metada::backends::wrf