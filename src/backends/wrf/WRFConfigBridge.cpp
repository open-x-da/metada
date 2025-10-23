/**
 * @file WRFConfigBridge.cpp
 * @brief Implementation of WRF configuration and domain management bridge
 * @ingroup backends
 */

#include "WRFConfigBridge.hpp"

#include <iostream>
#include <stdexcept>

namespace metada::backends::wrf {

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
    : domain_id_(domain_id), domain_allocated_(false) {
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
