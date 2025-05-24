/**
 * @file MACOMModel.hpp
 * @brief MACOM model backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <memory>
#include <netcdf>
#include <string>
#include <unordered_map>
#include <vector>

#include "include/MACOMlogging.hpp"

// Include MACOM specific interfaces
// #include "MACOMState.hpp"
#include "include/MACOMFortranInterface.hpp"  // Changed from MACOMFortranInterface.hpp

// Metada framework utilities (if needed)
#include "DateTime.hpp"
#include "Duration.hpp"

namespace metada::backends::macom {

template <typename ConfigBackend, typename GeometryBackend>
class MACOMState;

/**
 * @brief MACOM model backend implementation.
 * @details This class wraps the Fortran-based MACOM ocean model,
 *          making it compatible with the Metada framework via
 * MACOMFortranInterface.
 */
template <typename ConfigBackend, typename StateBackend>
class MACOMModel {
 public:
  /**
   * @brief Default constructor is deleted
   */
  MACOMModel() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  MACOMModel(const MACOMModel&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  MACOMModel& operator=(const MACOMModel&) = delete;

  /**
   * @brief Constructor that takes a configuration backend.
   * @param config Configuration containing MACOM model options.
   * @param mpi_comm MPI communicator to be used by the model.
   */
  explicit MACOMModel(const ConfigBackend& config);

  /**
   * @brief Move constructor.
   * @param other MACOM model backend to move from.
   */
  MACOMModel(MACOMModel&& other) noexcept;

  /**
   * @brief Move assignment operator.
   * @param other MACOM model backend to move from.
   * @return Reference to this model after assignment.
   */
  MACOMModel& operator=(MACOMModel&& other) noexcept;

  /**
   * @brief Destructor.
   */
  ~MACOMModel();

  /**
   * @brief Initialize the model with its configuration.
   *        This will initialize MPI (if not already done externally),
   *        read the namelist, and initialize model components via
   * FortranInterface.
   * @param config Configuration (already passed to constructor, can be used for
   * re-init or verification).
   */
  void initialize(const ConfigBackend& config);

  /**
   * @brief Reset the model to its initial state (placeholder).
   */
  void reset();

  /**
   * @brief Finalize the model, releasing resources.
   *        This will finalize model components and MPI (if initialized by this
   * model) via FortranInterface.
   */
  void finalize();

  /**
   * @brief Run the model.
   *        For now, this will call the Fortran model's main run loop or a
   * series of steps.
   * @param initialState Initial state of the model (currently unused
   * placeholder).
   * @param finalState Final state after model integration (output, currently
   * unused placeholder).
   * @throws std::runtime_error If model run fails.
   */
  void run(const StateBackend& initialState, StateBackend& finalState);

  /**
   * @brief Check if the model is initialized.
   * @return True if initialized, false otherwise.
   */
  bool isInitialized() const { return initialized_; }

 private:
  /**
   * @brief Run a single time step of the model
   *
   * @param inState Input state
   * @param outState Output state after time step
   * @param dt Time step size (in seconds)
   */
  void timeStep(const StateBackend& inState, StateBackend& outState,
                Duration dt);

  void applyBoundaryConditions(StateBackend& state) const;

  // Model configuration
  bool initialized_ = false;

  // Model parameters
  std::string configFile_;  // Path to namelist.nmefc_macom
  std::string outputDir_;   // Output directory

  // Time control
  DateTime startTime_;    // Start time in seconds
  DateTime currentTime_;  // Current model time in seconds
  DateTime endTime_;      // End time in seconds
  Duration timeStep_;     // Time step in seconds

  // Fortran interface
  std::unique_ptr<MACOMFortranInterface> fortranInterface_;
  int mpi_rank_;

  // Dynamical core options
  std::string advectionScheme_ = "WENO";
};

// Constructor implementation with ConfigBackend
template <typename ConfigBackend, typename StateBackend>
MACOMModel<ConfigBackend, StateBackend>::MACOMModel(const ConfigBackend& config)
    : initialized_(false) {
  // Create Fortran interface
  fortranInterface_ = std::make_unique<MACOMFortranInterface>();
  initialize(config);
}

template <typename ConfigBackend, typename StateBackend>
MACOMModel<ConfigBackend, StateBackend>::MACOMModel(MACOMModel&& other) noexcept
    : initialized_(other.initialized_),
      configFile_(std::move(other.configFile_)),
      outputDir_(std::move(other.outputDir_)) {
  other.initialized_ = false;
  other.currentTime_ = startTime_;
}

template <typename ConfigBackend, typename StateBackend>
MACOMModel<ConfigBackend, StateBackend>&
MACOMModel<ConfigBackend, StateBackend>::operator=(
    MACOMModel<ConfigBackend, StateBackend>&& other) noexcept {
  if (this != &other) {
    initialized_ = other.initialized_;
    configFile_ = std::move(other.configFile_);
    outputDir_ = std::move(other.outputDir_);

    // Reset the moved-from object
    other.initialized_ = false;
    other.currentTime_ = startTime_;
  }
  return *this;
}

template <typename ConfigBackend, typename StateBackend>
MACOMModel<ConfigBackend, StateBackend>::~MACOMModel() {
  if (initialized_) {
    try {
      finalize();
    } catch (const std::exception& e) {
      std::cerr << "Error during MACOM model finalization: " << e.what()
                << std::endl;
    }
  }
}

template <typename ConfigBackend, typename StateBackend>
void MACOMModel<ConfigBackend, StateBackend>::initialize(
    const ConfigBackend& config) {
  // Access the nested configuration structure correctly
  try {
    // First, check if we have a time_control subsection in the config
    if (config.HasKey("time_control")) {
      // Direct access to time_control subsection
      auto time_control = config.CreateSubsection("time_control");
      startTime_ = DateTime(time_control.Get("start_datetime").asString());
      endTime_ = DateTime(time_control.Get("end_datetime").asString());
      timeStep_ = Duration(time_control.Get("time_step").asString());

      MACOM_LOG_INFO(
          "MACOMModel",
          "Forecast period: " + time_control.Get("start_datetime").asString() +
              " to " + time_control.Get("end_datetime").asString());
    } else {
      MACOM_LOG_WARNING("MACOMModel",
                        "Configuration does not contain time_control or "
                        "model.time_control sections");
    }

    if (initialized_) {
      MACOM_LOG_WARNING(
          "MACOMModel",
          "Already initialized. Finalize first to re-initialize.");
      return;
    }

    MACOM_LOG_INFO("MACOMModel", "Initializing...");

    // Initialize MPI through Fortran interface
    fortranInterface_->initializeMPI();
    mpi_rank_ = fortranInterface_->getRank();
    MACOM_LOG_INFO("MACOMModel", "MPI Initialized via Fortran. Model Rank: " +
                                     std::to_string(mpi_rank_));

    // Read namelist
    fortranInterface_->readNamelist();
    MACOM_LOG_INFO("MACOMModel", "Fortran Namelist Read.");

    // Initialize model components
    fortranInterface_->initializeModelComponents();
    MACOM_LOG_INFO("MACOMModel", "Fortran model components initialized.");

    initialized_ = true;
    MACOM_LOG_INFO("MACOMModel", "Initialization complete.");
  } catch (const std::exception& e) {
    initialized_ = false;
    throw std::runtime_error(std::string("MACOMModel initialization failed: ") +
                             e.what());
  }
}

template <typename ConfigBackend, typename StateBackend>
void MACOMModel<ConfigBackend, StateBackend>::reset() {}

template <typename ConfigBackend, typename StateBackend>
void MACOMModel<ConfigBackend, StateBackend>::finalize() {
  if (!initialized_ || !fortranInterface_) {
    return;
  }
  MACOM_LOG_INFO("MACOMModel", "Finalizing...");

  try {
    fortranInterface_->finalizeModelComponents();
    MACOM_LOG_INFO("MACOMModel", "Fortran model components finalized.");

    fortranInterface_->finalizeMPI();
    MACOM_LOG_INFO("MACOMModel", "Fortran MPI finalized.");

  } catch (const std::exception& e) {
    MACOM_LOG_ERROR("MACOMModel",
                    "Error during finalization: " + std::string(e.what()));
  }
  initialized_ = false;
  mpi_rank_ = -1;
  MACOM_LOG_INFO("MACOMModel", "Finalization complete.");
}

template <typename ConfigBackend, typename StateBackend>
void MACOMModel<ConfigBackend, StateBackend>::run(
    const StateBackend& initialState, StateBackend& finalState) {
  // TODO: Implement the run method
  if (!initialized_) {
    throw std::runtime_error("Cannot run uninitialized MACOM model");
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

    // // Use a temporary state for the time stepping process
    // auto tempState = initialState.clone();

    // Main time stepping loop
    while (currentTime_ < endTime_) {
      // Ensure we don't overshoot the end time
      if (currentTime_ + timeStep_ > endTime_) {
        timeStep_ = endTime_ - currentTime_;
      }

      // // Perform a single time step
      // timeStep(*workingState, *tempState, timeStep_);

      // // Swap states for next iteration
      // std::swap(workingState, tempState);

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

        MACOM_LOG_INFO("MACOMModel", "MACOM model integration: " +
                                         std::to_string(progress) + ")");
      }
    }

    // Copy the final state
    finalState = std::move(*workingState);

    MACOM_LOG_INFO("MACOMModel",
                   "MACOM model integration completed successfully");

  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("MACOM model run failed: ") +
                             e.what());
  }
}

template <typename ConfigBackend, typename StateBackend>
void MACOMModel<ConfigBackend, StateBackend>::timeStep(
    const StateBackend& inState, StateBackend& outState,
    [[maybe_unused]] Duration dt) {
  // TODO: Implement the time step method
  outState = std::move(*inState.clone());

  try {
    // Get variable names from the state
    const auto& varNames = inState.getVariableNames();
    if (advectionScheme_ == "WENO") {
      // WENO advection scheme (simplified)
      for (const auto& varName : varNames) {
        outState.setActiveVariable(varName);
        // Perform advection calculations here
      }
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("MACOM model time step failed: ") +
                             e.what());
  }
}

}  // namespace metada::backends::macom