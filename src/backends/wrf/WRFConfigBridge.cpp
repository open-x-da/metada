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

namespace metada::backends::wrf {

namespace fs = std::filesystem;

// Initialize static member
bool WRFConfigManager::wrfda_modules_initialized_ = false;

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
    wrf_initial_config_();

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
