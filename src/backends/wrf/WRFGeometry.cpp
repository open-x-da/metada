/**
 * @file WRFGeometry.cpp
 * @brief WRF geometry backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#include "WRFGeometry.hpp"

#include <netcdf>
#include <stdexcept>
#include <string>

namespace metada::backends::wrf {

// Move constructor implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>::WRFGeometry(WRFGeometry&& other) noexcept
    : wrfFilename_(std::move(other.wrfFilename_)),
      timestamp_(std::move(other.timestamp_)),
      initialized_(other.initialized_),
      longitude_(std::move(other.longitude_)),
      latitude_(std::move(other.latitude_)),
      elevation_(std::move(other.elevation_)),
      nx_(other.nx_),
      ny_(other.ny_),
      nz_(other.nz_),
      periodicX_(other.periodicX_),
      periodicY_(other.periodicY_),
      periodicZ_(other.periodicZ_) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.nx_ = 0;
  other.ny_ = 0;
  other.nz_ = 0;
}

// Move assignment operator implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>& WRFGeometry<ConfigBackend>::operator=(
    WRFGeometry&& other) noexcept {
  if (this != &other) {
    wrfFilename_ = std::move(other.wrfFilename_);
    timestamp_ = std::move(other.timestamp_);
    initialized_ = other.initialized_;
    longitude_ = std::move(other.longitude_);
    latitude_ = std::move(other.latitude_);
    elevation_ = std::move(other.elevation_);
    nx_ = other.nx_;
    ny_ = other.ny_;
    nz_ = other.nz_;
    periodicX_ = other.periodicX_;
    periodicY_ = other.periodicY_;
    periodicZ_ = other.periodicZ_;

    // Reset the moved-from object
    other.initialized_ = false;
    other.nx_ = 0;
    other.ny_ = 0;
    other.nz_ = 0;
  }
  return *this;
}

// Clone implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend> WRFGeometry<ConfigBackend>::clone() const {
  // Use the specialized constructor
  WRFGeometry<ConfigBackend> clone(wrfFilename_, timestamp_, periodicX_,
                                   periodicY_, periodicZ_);

  // Copy data members
  clone.longitude_ = longitude_;
  clone.latitude_ = latitude_;
  clone.elevation_ = elevation_;
  clone.nx_ = nx_;
  clone.ny_ = ny_;
  clone.nz_ = nz_;
  clone.initialized_ = initialized_;

  return clone;
}

// Periodicity checks implementation
template <typename ConfigBackend>
bool WRFGeometry<ConfigBackend>::isPeriodicX() const {
  return periodicX_;
}

template <typename ConfigBackend>
bool WRFGeometry<ConfigBackend>::isPeriodicY() const {
  return periodicY_;
}

template <typename ConfigBackend>
bool WRFGeometry<ConfigBackend>::isPeriodicZ() const {
  return periodicZ_;
}

// isInitialized implementation
template <typename ConfigBackend>
bool WRFGeometry<ConfigBackend>::isInitialized() const {
  return initialized_;
}

// Get accessor methods implementation
template <typename ConfigBackend>
const xt::xarray<double>& WRFGeometry<ConfigBackend>::getLongitude() const {
  return longitude_;
}

template <typename ConfigBackend>
const xt::xarray<double>& WRFGeometry<ConfigBackend>::getLatitude() const {
  return latitude_;
}

template <typename ConfigBackend>
const xt::xarray<double>& WRFGeometry<ConfigBackend>::getElevation() const {
  return elevation_;
}

// Halo exchange implementation
template <typename ConfigBackend>
void WRFGeometry<ConfigBackend>::haloExchangeImpl(void* state_ptr) {
  // Implementation depends on state type and halo exchange requirements
  // This is a placeholder, actual implementation would transfer data between
  // neighboring domains or processes
  if (!initialized_) {
    throw std::runtime_error(
        "Cannot perform halo exchange on uninitialized geometry");
  }

  // Cast state_ptr to appropriate type and perform exchange
  // This is implementation-specific and would depend on the WRF model's needs
}

template <typename ConfigBackend>
WRFGeometry<ConfigBackend>::WRFGeometry(const std::string& fn,
                                        const std::string& ts, bool px, bool py,
                                        bool pz)
    : wrfFilename_(fn),
      timestamp_(ts),
      initialized_(false),
      periodicX_(px),
      periodicY_(py),
      periodicZ_(pz) {
  // The rest will be populated by clone() method
}

}  // namespace metada::backends::wrf