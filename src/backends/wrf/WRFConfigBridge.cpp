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
bool WRFConfigManager::wrf_modules_initialized_ = false;

void WRFConfigManager::initializeWRFModules() {
  if (!wrf_modules_initialized_) {
    std::cout << "Initializing WRF modules..." << std::endl;

    // Phase 1: Core modules initialization
    wrf_init_modules_(1);

    // Phase 2: Advanced modules initialization
    wrf_init_modules_(2);

    wrf_modules_initialized_ = true;

    // Verify initialization
    if (!wrf_is_initialized_()) {
      throw std::runtime_error("WRF module initialization failed");
    }

    std::cout << "WRF modules initialized successfully" << std::endl;
  }
}

WRFConfigManager::WRFConfigManager(int domain_id, bool allocate_domain)
    : domain_id_(domain_id), domain_allocated_(false) {
  // Step 1: Initialize WRF modules (once per process)
  try {
    initializeWRFModules();
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Failed to initialize WRF modules: ") +
                             e.what());
  }

  // Step 2: Read namelist.input configuration
  // This populates the module-level model_config_rec
  // Note: These values are mainly for WRF internal consistency.
  // The actual grid geometry comes from the NetCDF file global attributes.
  try {
    std::cout << "Reading namelist.input configuration..." << std::endl;
    wrf_initial_config_();
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::string("Failed to initialize WRF configuration: ") + e.what());
  }

  // Step 3: Extract domain-specific configuration into module-level
  // config_flags_
  try {
    std::cout << "Extracting configuration for domain " << domain_id_ << "..."
              << std::endl;
    wrf_model_to_grid_config_(&domain_id_);
    initialized_ = true;
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::string("Failed to extract domain configuration: ") + e.what());
  }

  // Step 4: Allocate WRF domain structure (if requested)
  if (allocate_domain) {
    std::cout << "Allocating WRF domain structure for domain " << domain_id_
              << "..." << std::endl;

    int ierr = wrf_alloc_domain_(domain_id_);
    if (ierr != 0) {
      throw std::runtime_error(
          "Failed to allocate WRF domain structure for domain " +
          std::to_string(domain_id_) + ", error code: " + std::to_string(ierr));
    }

    // Verify domain was allocated
    if (!wrf_domain_exists_(domain_id_)) {
      throw std::runtime_error("WRF domain " + std::to_string(domain_id_) +
                               " was not allocated successfully");
    }

    domain_allocated_ = true;
    std::cout << "WRF domain " << domain_id_ << " allocated successfully"
              << std::endl;
  }
}

WRFConfigManager::~WRFConfigManager() {
  // Deallocate WRF domain if we allocated it
  if (domain_allocated_) {
    try {
      std::cout << "Deallocating WRF domain " << domain_id_ << "..."
                << std::endl;
      wrf_dealloc_domain_(domain_id_);
      domain_allocated_ = false;
    } catch (...) {
      // Suppress exceptions in destructor
      std::cerr << "Warning: Exception during WRF domain deallocation"
                << std::endl;
    }
  }
}

bool WRFConfigManager::domainExists() const {
  return wrf_domain_exists_(domain_id_);
}

void* WRFConfigManager::getConfigFlagsPtr() const {
  return wrf_get_config_flags_ptr_();
}

size_t WRFConfigManager::getConfigFlagsSize() const {
  return wrf_get_config_flags_size_();
}

void* WRFConfigManager::getGridPtr() const {
  return wrf_get_grid_ptr_(domain_id_);
}

}  // namespace metada::backends::wrf
