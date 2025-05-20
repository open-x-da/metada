/**
 * @file MACOMModel.hpp
 * @brief MACOM model backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <memory>
#include <stdexcept>  // For std::runtime_error
#include <string>
#include <vector>

// Include MACOM specific interfaces
#include "MACOMState.hpp"
#include "include/MACOMFortranInterface.hpp"  // Changed from FortranInterface.hpp

// Metada framework utilities (if needed)
// #include "DateTime.hpp"
// #include "Duration.hpp"

namespace metada::backends::macom {

template <typename ConfigBackend>
class MACOMState;

/**
 * @brief MACOM model backend implementation.
 * @details This class wraps the Fortran-based MACOM ocean model,
 *          making it compatible with the Metada framework via
 * MACOMFortranInterface.
 */
template <typename ConfigBackend>
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
  void initialize([[maybe_unused]] const ConfigBackend& config);

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
  void run([[maybe_unused]] const MACOMState<ConfigBackend>& initialState,
           [[maybe_unused]] MACOMState<ConfigBackend>& finalState);

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
  void timeStep([[maybe_unused]] const MACOMState<ConfigBackend>& inState,
                [[maybe_unused]] MACOMState<ConfigBackend>& outState,
                [[maybe_unused]] double dt);

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
};

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
MACOMModel<ConfigBackend>::MACOMModel(const ConfigBackend& config)
    : initialized_(false) {
  // Create Fortran interface
  fortranInterface_ = std::make_unique<MACOMFortranInterface>();
  initialize(config);
}

template <typename ConfigBackend>
MACOMModel<ConfigBackend>::MACOMModel(MACOMModel&& other) noexcept
    : initialized_(other.initialized_),
      configFile_(std::move(other.configFile_)),
      outputDir_(std::move(other.outputDir_)) {
  other.initialized_ = false;
  other.currentTime_ = startTime_;
}

template <typename ConfigBackend>
MACOMModel<ConfigBackend>& MACOMModel<ConfigBackend>::operator=(
    MACOMModel&& other) noexcept {
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

template <typename ConfigBackend>
MACOMModel<ConfigBackend>::~MACOMModel() {
  if (initialized_) {
    try {
      finalize();
    } catch (const std::exception& e) {
      std::cerr << "Error during MACOM model finalization: " << e.what()
                << std::endl;
    }
  }
}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::initialize(const ConfigBackend& config) {
  std::string start_datetime = config.Get("start_datetime").asString();
  std::string end_datetime = config.Get("end_datetime").asString();
  std::string time_step = config.Get("time_step").asString();

  if (initialized_) {
    std::cout << "[MACOMModel] Start time: " << start_datetime << std::endl;
    std::cout << "[MACOMModel] End time: " << end_datetime << std::endl;
    std::cout << "[MACOMModel] Time step: " << time_step << std::endl;

    std::cout
        << "[MACOMModel] Already initialized. Finalize first to re-initialize."
        << std::endl;
    return;
  }

  std::cout << "[MACOMModel] Initializing..." << std::endl;
  try {
    // Initialize MPI through Fortran interface
    fortranInterface_->initializeMPI();
    mpi_rank_ = fortranInterface_->getRank();
    std::cout << "[MACOMModel] MPI Initialized via Fortran. Model Rank: "
              << mpi_rank_ << std::endl;

    // Read namelist
    fortranInterface_->readNamelist();
    std::cout << "[MACOMModel] Fortran Namelist Read." << std::endl;

    // Initialize model components
    fortranInterface_->initializeModelComponents();
    std::cout << "[MACOMModel] Fortran model components initialized."
              << std::endl;

    initialized_ = true;
    std::cout << "[MACOMModel] Initialization complete." << std::endl;

  } catch (const std::exception& e) {
    initialized_ = false;
    throw std::runtime_error(
        std::string("[MACOMModel] Initialization failed: ") + e.what());
  }
}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::reset() {}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::finalize() {
  if (!initialized_ || !fortranInterface_) {
    return;
  }
  std::cout << "[MACOMModel] Finalizing..." << std::endl;
  try {
    fortranInterface_->finalizeModelComponents();
    std::cout << "[MACOMModel] Fortran model components finalized."
              << std::endl;

    fortranInterface_->finalizeMPI();
    std::cout << "[MACOMModel] Fortran MPI finalized." << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "[MACOMModel] Error during finalization: " << e.what()
              << std::endl;
  }
  initialized_ = false;
  mpi_rank_ = -1;
  std::cout << "[MACOMModel] Finalization complete." << std::endl;
}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::run(
    [[maybe_unused]] const MACOMState<ConfigBackend>& initialState,
    [[maybe_unused]] MACOMState<ConfigBackend>& finalState) {}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::timeStep(
    [[maybe_unused]] const MACOMState<ConfigBackend>& inState,
    [[maybe_unused]] MACOMState<ConfigBackend>& outState,
    [[maybe_unused]] double dt) {}

// // Apply boundary conditions implementation
// template <typename ConfigBackend>
// void MACOMModel<ConfigBackend>::applyBoundaryConditions(
//     MACOMState<ConfigBackend>& state) const {
//   // Apply appropriate boundary conditions to each variable
//   // For demonstration purposes, we'll just apply simple conditions

//   const auto& varNames = state.getVariableNames();

//   for (const auto& varName : varNames) {
//     state.setActiveVariable(varName);

//     // In a real implementation, this would apply periodic, open, or fixed
//     // boundary conditions depending on the variable and domain configuration
//   }
}  // namespace metada::backends::macom