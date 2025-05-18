/**
 * @file MACOMModel.cpp
 * @brief Implementation of the MACOM model backend.
 * @ingroup backends_macom
 * @author Metada Framework Team
 */

#include "../include/MACOMModel.hpp"

#include <iostream>  // For std::cout, std::cerr

namespace metada::backends::macom {

template <typename ConfigBackend>
MACOMModel<ConfigBackend>::MACOMModel(const ConfigBackend& config,
                                      int dummy_mpi_param)
    : config_(config), initialized_(false), mpi_rank_(-1) {
  std::cout << "[MACOMModel] Constructor called with dummy int: "
            << dummy_mpi_param << std::endl;
  < std::endl;
  mm related setup if needed, but FortranInterface handles its
  // conversion.
}

template <typename ConfigBackend>
MACOMModel<ConfigBackend>::MACOMModel(MACOMModel&& other) noexcept
    : config_(other.config_),  // Const reference, so just copy/alias it.
      initialized_(other.initialized_),
      mpi_rank_(other.mpi_rank_),
      fortranInterface_(std::move(other.fortranInterface_)) {
  other.initialized_ = false;
  other.mpi_rank_ = -1;
  std::cout << "[MACOMModel] Move constructor called." << std::endl;
}

template <typename ConfigBackend>
MACOMModel<ConfigBackend>& MACOMModel<ConfigBackend>::operator=(
    MACOMModel&& other) noexcept {
  if (this != &other) {
    // config_ is a const reference, cannot be reassigned.
    // Ensure consistent behavior if your design allows reassigning MACOMModel
    // instances. Typically, models might not be copy/move assignable if they
    // hold significant state or external (e.g., Fortran) resources directly
    // tied to one instance. For this example, we allow moving the interface and
    // status.
    initialized_ = other.initialized_;
    mpi_rank_ = other.mpi_rank_;
    fortranInterface_ = std::move(other.fortranInterface_);

    other.initialized_ = false;
    other.mpi_rank_ = -1;
  }
  std::cout << "[MACOMModel] Move assignment operator called." << std::endl;
  return *this;
}

template <typename ConfigBackend>
MACOMModel<ConfigBackend>::~MACOMModel() {
  std::cout << "[MACOMModel] Destructor called." << std::endl;
  try {
    if (initialized_ && fortranInterface_) {
      finalize();  // Attempt to finalize if it was initialized
    }
  } catch (const std::exception& e) {
    std::cerr << "[MACOMModel] Exception during finalize in destructor: "
              << e.what() << std::endl;
  } catch (...) {
    std::cerr << "[MACOMModel] Unknown exception during finalize in destructor."
              << std::endl;
  }
  // The unique_ptr fortranInterface_ will be automatically destroyed here,
  // and its destructor will handle its own cleanup (like Fortran MPI finalize).
}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::initialize(
    const ConfigBackend& /*config_param*/) {
  // The config_param could be used to re-validate or update, but here we use
  // the member config_.
  if (initialized_) {
    std::cout
        << "[MACOMModel] Already initialized. Finalize first to re-initialize."
        << std::endl;
    return;
  }
  if (!fortranInterface_) {
    throw std::runtime_error(
        "[MACOMModel] Fortran interface is not available for initialization.");
  }

  std::cout << "[MACOMModel] Initializing..." << std::endl;
  try {
    // 1. Initialize MPI through the Fortran interface
    //    This also sets up the rank known to Fortran.
    fortranInterface_->initializeMPI();
    mpi_rank_ = fortranInterface_->getRank();  // Get rank from Fortran side
    std::cout << "[MACOMModel] MPI Initialized via Fortran. Model Rank: "
              << mpi_rank_ << std::endl;

    // 2. Read Namelist (Fortran configuration)
    //    This should be done after MPI init if namelist reading is
    //    rank-dependent (e.g. rank 0 reads)
    fortranInterface_->readNamelist();
    std::cout << "[MACOMModel] Fortran Namelist Read." << std::endl;

    // 3. Initialize core model components in Fortran
    fortranInterface_->initializeModelComponents();
    std::cout << "[MACOMModel] Fortran model components initialized."
              << std::endl;

    // Placeholder: Extract time control parameters if needed from config_ or
    // Fortran Example: if (config_.HasKey("macom_max_iterations")) {
    //     max_fortran_iterations_ =
    //     config_.Get("macom_max_iterations").asInt();
    // }
    // current_fortran_iteration_ = 0; // Or get from Fortran if restart

    initialized_ = true;
    std::cout << "[MACOMModel] Initialization complete." << std::endl;

  } catch (const std::exception& e) {
    initialized_ = false;  // Ensure not marked as initialized on error
    // Attempt to clean up MPI if this class initiated it and an error occurred
    // mid-init.
    if (fortranInterface_ && fortranInterface_->isMPIInitialized()) {
      try {
        fortranInterface_->finalizeMPI();
      } catch (const std::exception& e_mpi) {
        std::cerr
            << "[MACOMModel] Exception during MPI finalize after init error: "
            << e_mpi.what() << std::endl;
      }
    }
    throw std::runtime_error(
        std::string("[MACOMModel] Initialization failed: ") + e.what());
  }
}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::reset() {
  if (!initialized_ || !fortranInterface_) {
    throw std::runtime_error(
        "[MACOMModel] Cannot reset: Model not initialized or interface "
        "missing.");
  }
  std::cout << "[MACOMModel] Resetting model (placeholder)." << std::endl;
  // Actual reset logic would involve:
  // 1. Calling a Fortran routine to reset its internal state via
  // fortranInterface_.
  // 2. Re-initializing MACOMState objects if they hold data that needs reset.
  // For now, this is a conceptual placeholder.
  // Example: fortranInterface_->callFortranReset();
}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::finalize() {
  if (!initialized_ || !fortranInterface_) {
    // std::cout << "[MACOMModel] Finalize called on non-initialized or already
    // finalized model, or interface missing." << std::endl;
    return;
  }
  std::cout << "[MACOMModel] Finalizing..." << std::endl;
  try {
    fortranInterface_->finalizeModelComponents();
    std::cout << "[MACOMModel] Fortran model components finalized."
              << std::endl;

    // MPI finalization is handled by MACOMFortranInterface destructor if it
    // initialized MPI. If MACOMModel explicitly called initializeMPI on the
    // interface, it could also call finalizeMPI here. However, to avoid double
    // finalization, it's safer if the interface manages its own MPI lifecycle.
    // If MACOMModel wants to ensure MPI is finalized NOW, it could call it:
    // if (fortranInterface_->isMPIInitialized()) { // Check if it was
    // initialized by the interface
    //     fortranInterface_->finalizeMPI();
    //     std::cout << "[MACOMModel] Fortran MPI Finalized via explicit call."
    //     << std::endl;
    // }

  } catch (const std::exception& e) {
    std::cerr << "[MACOMModel] Error during finalization: " << e.what()
              << std::endl;
    // Don't rethrow if called from destructor, but this finalize can be called
    // elsewhere.
  }
  initialized_ = false;  // Mark as not initialized
  mpi_rank_ = -1;
  std::cout << "[MACOMModel] Finalization sequence in MACOMModel complete."
            << std::endl;
}

template <typename ConfigBackend>
void MACOMModel<ConfigBackend>::run(
    [[maybe_unused]] const MACOMState<ConfigBackend>& initialState,
    [[maybe_unused]] MACOMState<ConfigBackend>& finalState) {
  if (!initialized_ || !fortranInterface_) {
    throw std::runtime_error(
        "[MACOMModel] Cannot run: Model not initialized or interface missing.");
  }
  std::cout << "[MACOMModel] Running model..." << std::endl;

  // TODO: Implement data transfer from initialState to Fortran model if needed.
  // This would involve MACOMState having methods to get data (e.g., as arrays)
  // and MACOMFortranInterface having methods to set Fortran arrays (not yet
  // defined).

  // --- Example of a simple loop ---
  // This assumes c_macom_run_model_step in Fortran handles one iteration,
  // or if it handles the whole loop, then this C++ loop is not needed.
  // For now, let's assume forecast.cpp or a similar top-level application
  // might control the number of steps or overall duration. Here, we might just
  // call the main Fortran execution routine once if it handles its own loop. If
  // your Fortran `callCSP2()` runs the whole simulation, then just call that
  // once. If `c_macom_run_model_step` is truly one step, you need a loop here
  // or in forecast.cpp.

  // For this example, let's assume `c_macom_run_model_step` is called for one
  // main execution block. If you have a specific number of iterations defined
  // (e.g. from namelist) int start_iter = 0; // Or get from config/Fortran int
  // max_iter = 1;   // Or get from config/Fortran, e.g. max_fortran_iterations_
  // For a single call to a potentially looping Fortran routine:
  try {
    fortranInterface_->runModelStep(0);  // Passing 0 as a placeholder iteration
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::string("[MACOMModel] Error during model run: ") + e.what());
  }

  // If c_macom_run_model_step is literally one step:
  // for (int i = start_iter; i < max_iter; ++i) {
  //     std::cout << "[MACOMModel] Executing model step " << i << std::endl;
  //     fortranInterface_->runModelStep(i);
  // }

  // MPI Barrier might be inside the Fortran c_macom_run_model_step if it's a
  // collective operation. Or, if it represents the full simulation, a barrier
  // might be after it completes. MPI_Barrier(MPI_COMM_WORLD); // Example
  // barrier, ensure this matches Fortran logic

  // TODO: Implement data transfer from Fortran model to finalState.

  std::cout << "[MACOMModel] Model run finished." << std::endl;
}

template <typename ConfigBackend>
bool MACOMModel<ConfigBackend>::isInitialized() const {
  return initialized_;
}

// Explicit template instantiations if needed for your build system, or if
// ConfigBackend is fixed. This is often needed if the .cpp file is compiled
// separately. Replace `YourActualConfigBackend` with the concrete type(s) you
// use, e.g., `backends::config::YamlConfig`. If you support multiple, you might
// need to list them or use other techniques. Example: template class
// MACOMModel<metada::backends::config::YamlConfig>; #ifdef CONFIG_BACKEND_JSON
// // Assuming you have defines for which configs are active template class
// MACOMModel<metada::backends::config::JsonConfig>; #endif

}  // namespace metada::backends::macom
