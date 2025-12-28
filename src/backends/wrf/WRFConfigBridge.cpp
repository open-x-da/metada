/**
 * @file WRFConfigBridge.cpp
 * @brief Implementation of WRF configuration and domain management bridge
 * @ingroup backends
 */

#include "WRFConfigBridge.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <system_error>

#ifdef METADA_USE_MPI
#include <mpi.h>
#endif

namespace metada::backends::wrf {

namespace fs = std::filesystem;

// Initialize static member
bool WRFConfigManager::wrfda_modules_initialized_ = false;

void setupWRFMPICommunicator() {
#ifdef METADA_USE_MPI
  // Check if MPI is initialized
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);

  if (mpi_initialized) {
    // Get MPI_COMM_WORLD and pass it to WRF's domain manager
    // This allows WRFDA to use the same MPI communicator as METADA
    MPI_Comm world_comm = MPI_COMM_WORLD;

    // Convert C communicator to Fortran integer handle using standard MPI
    // function MPI_Comm_c2f is the correct and portable way to convert MPI
    // communicators between C and Fortran. This handles different MPI
    // implementations correctly, including Open MPI where valid handles may be
    // 0.
    MPI_Fint f_handle = MPI_Comm_c2f(world_comm);
    MPI_Fint f_null = MPI_Comm_c2f(MPI_COMM_NULL);

    // Validate conversion result - compare against MPI_COMM_NULL handle
    // Note: In Open MPI, f_handle can be 0 for valid communicators, so we
    //       must compare against f_null rather than checking for 0.
    if (f_handle == f_null) {
      throw std::runtime_error(
          "Failed to convert MPI_COMM_WORLD to Fortran communicator. "
          "Communicator conversion returned null handle. "
          "Check MPI installation and initialization.");
    }

    std::cout << "Setting WRF domain manager communicator from METADA MPI..."
              << std::endl;
    // Note: The Fortran function wrf_set_dm_communicator_from_metada_ will
    //       check STUBMPI (set during WRFDA compilation). If WRFDA was built
    //       without MPI support (STUBMPI defined), the Fortran function will
    //       be a no-op.
    wrf_set_dm_communicator_from_metada_(static_cast<int>(f_handle));
    std::cout << "WRF domain manager communicator set successfully"
              << std::endl;
  } else {
    std::cout << "MPI not initialized - skipping WRF communicator setup"
              << std::endl;
  }
#else
  // Stub: do nothing when MPI is disabled in METADA
  // Note: Even if METADA_USE_MPI is OFF, WRFDA might have been built with MPI.
  //       The Fortran function will handle this via STUBMPI check.
  std::cout << "MPI disabled in METADA - skipping WRF communicator setup"
            << std::endl;
#endif
}

void WRFConfigManager::initializeWRFDAModules() {
  if (!wrfda_modules_initialized_) {
    std::cout << "Initializing WRFDA modules..." << std::endl;

    // Phase 1: Core modules initialization (configuration and constants only)
    wrfda_init_modules_(1);

    // Initialize WRFU time manager (must be between phase 1 and 2)
    std::cout << "Initializing WRFU time utilities..." << std::endl;
    wrfda_wrfu_initialize_();

    // Phase 2: Advanced modules initialization
    wrfda_init_modules_(2);

    wrfda_modules_initialized_ = true;

    // Verify initialization
    if (!wrfda_is_initialized_()) {
      throw std::runtime_error("WRFDA module initialization failed");
    }

    std::cout << "WRFDA modules initialized successfully" << std::endl;
  }
}

WRFConfigManager::WRFConfigManager(int domain_id, bool allocate_domain)
    : domain_id_(domain_id),
      domain_allocated_(false),
      trace_session_started_(false) {
  // Step 0: Setup MPI communicator for WRF/WRFDA (if MPI is enabled)
  // This must be done before WRFDA modules are initialized so that
  // WRFDA can use the communicator for parallel operations
  try {
    setupWRFMPICommunicator();
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::string("Failed to setup WRF MPI communicator: ") + e.what());
  }

  // Step 1: Initialize WRFDA modules (once per process)
  try {
    initializeWRFDAModules();
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::string("Failed to initialize WRFDA modules: ") + e.what());
  }

  // Step 2: Read namelist.input configuration
  // This populates the module-level model_config_rec
  // Note: namelist.input values MUST match the NetCDF file global attributes
  try {
    std::cout << "Reading namelist.input configuration..." << std::endl;
#ifdef METADA_USE_MPI
    // Use parallel initialization sequence (matches WRFDA's
    // da_wrfvar_init1.inc) This ensures config is read on rootproc only and
    // broadcast to all processes
    wrf_initial_config_parallel_();
#else
    // Serial mode: use simple initialization
    wrf_initial_config_();
#endif

    // Step 2a: Copy configuration from model_config_rec to da_control module
    // This replicates the logic from da_wrfvar_init1.inc that includes
    // config_assigns.inc
    std::cout << "Copying configuration to da_control module..." << std::endl;
    copy_config_to_da_control();

    // Step 2b: Validate configuration for common conflicts
    // This replicates the sanity checks from da_solve.inc
    std::cout << "Validating WRFDA configuration..." << std::endl;
    int error_code = 0;
    validate_wrfda_config(&error_code);
    if (error_code != 0) {
      throw std::runtime_error(
          "WRFDA configuration validation failed with error code: " +
          std::to_string(error_code) + ". Check output above for details.");
    }
    std::cout << "WRFDA configuration validated successfully" << std::endl;

    if (wrfda_trace_is_enabled()) {
      std::cout
          << "Trace output requested; ensuring ./trace directory exists..."
          << std::endl;
      std::error_code trace_ec;
      fs::create_directories("trace", trace_ec);
      if (trace_ec) {
        throw std::runtime_error("Failed to create trace directory 'trace': " +
                                 trace_ec.message());
      }
      wrfda_trace_initialize();
      trace_session_started_ = true;
      std::cout << "WRFDA tracing initialized (outputs under ./trace)"
                << std::endl;
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::string("Failed to initialize WRFDA configuration: ") + e.what());
  }

  // Step 3: Allocate and initialize WRFDA domain (following
  // da_wrfvar_init2.inc)
  if (allocate_domain) {
    std::cout << "Allocating WRFDA domain following standard workflow..."
              << std::endl;
    int ierr = wrfda_alloc_and_init_domain_(domain_id_);
    if (ierr != 0) {
      throw std::runtime_error("Failed to allocate WRFDA domain, error code: " +
                               std::to_string(ierr));
    }

    // Verify domain was allocated
    if (!wrfda_head_grid_allocated_()) {
      throw std::runtime_error(
          "WRFDA head_grid was not allocated successfully");
    }

    domain_allocated_ = true;
    std::cout << "WRFDA domain allocated and initialized successfully"
              << std::endl;
  }

  initialized_ = true;
  std::cout << "WRFDA configuration initialized for domain " << domain_id_
            << std::endl;
}

WRFConfigManager::~WRFConfigManager() {
  if (trace_session_started_) {
    wrfda_trace_finalize();
    trace_session_started_ = false;
  }
}

void* WRFConfigManager::getConfigFlagsPtr() const {
  return wrf_get_config_flags_ptr_();
}

size_t WRFConfigManager::getConfigFlagsSize() const {
  return wrf_get_config_flags_size_();
}

void* WRFConfigManager::getGridPtr() const {
  return wrfda_get_head_grid_ptr_();
}

}  // namespace metada::backends::wrf
