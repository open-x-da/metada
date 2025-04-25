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

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
WRFGeometry::WRFGeometry(const ConfigBackend& config)
    : wrfFilename_(config.getString("wrf.input_file", "")),
      timestamp_(config.getString("wrf.timestamp", "0000-00-00_00:00:00")),
      initialized_(false),
      periodicX_(config.getBool("wrf.periodic_x", false)),
      periodicY_(config.getBool("wrf.periodic_y", false)),
      periodicZ_(config.getBool("wrf.periodic_z", false)) {
  if (wrfFilename_.empty()) {
    throw std::runtime_error(
        "WRF input file path not specified in configuration");
  }

  // Load geometry data from WRF NetCDF file
  loadGeometryData(wrfFilename_, timestamp_);
  initialized_ = true;
}

// Move constructor implementation
WRFGeometry::WRFGeometry(WRFGeometry&& other) noexcept
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
WRFGeometry& WRFGeometry::operator=(WRFGeometry&& other) noexcept {
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
WRFGeometry WRFGeometry::clone() const {
  WRFGeometry clone_geometry(*this);  // Use the copy constructor
  return clone_geometry;
}

// totalGridSize implementation
std::size_t WRFGeometry::totalGridSize() const {
  return nx_ * ny_ * nz_;
}

// Periodicity checks implementation
bool WRFGeometry::isPeriodicX() const {
  return periodicX_;
}

bool WRFGeometry::isPeriodicY() const {
  return periodicY_;
}

bool WRFGeometry::isPeriodicZ() const {
  return periodicZ_;
}

// isInitialized implementation
bool WRFGeometry::isInitialized() const {
  return initialized_;
}

// Get accessor methods implementation
const xt::xarray<double>& WRFGeometry::getLongitude() const {
  return longitude_;
}

const xt::xarray<double>& WRFGeometry::getLatitude() const {
  return latitude_;
}

const xt::xarray<double>& WRFGeometry::getElevation() const {
  return elevation_;
}

// Halo exchange implementation
void WRFGeometry::haloExchangeImpl(void* state_ptr) {
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

// Private helper method to load geometry data from NetCDF file
void WRFGeometry::loadGeometryData(const std::string& filename,
                                   const std::string& timestamp) {
  try {
    // Open NetCDF file
    netCDF::NcFile wrf_file(filename, netCDF::NcFile::read);

    // Check if file is open
    if (!wrf_file.isNull()) {
      // Read grid dimensions
      auto dim_west_east = wrf_file.getDim("west_east");
      auto dim_south_north = wrf_file.getDim("south_north");
      auto dim_bottom_top = wrf_file.getDim("bottom_top");

      if (dim_west_east.isNull() || dim_south_north.isNull() ||
          dim_bottom_top.isNull()) {
        throw std::runtime_error("Missing required dimensions in WRF file");
      }

      nx_ = dim_west_east.getSize();
      ny_ = dim_south_north.getSize();
      nz_ = dim_bottom_top.getSize();

      // Read geographic data variables
      auto var_longitude = wrf_file.getVar("XLONG");
      auto var_latitude = wrf_file.getVar("XLAT");
      auto var_terrain = wrf_file.getVar("HGT");

      if (var_longitude.isNull() || var_latitude.isNull() ||
          var_terrain.isNull()) {
        throw std::runtime_error(
            "Missing required geographic variables in WRF file");
      }

      // Read time index for the requested timestamp if available
      size_t time_idx = 0;  // Default to first time step
      auto times_var = wrf_file.getVar("Times");
      if (!times_var.isNull()) {
        // Find matching timestamp (implementation dependent on WRF file format)
        // This is a placeholder for actual timestamp matching logic
      }

      // Allocate arrays for geographic data
      std::vector<size_t> start = {time_idx, 0, 0};
      std::vector<size_t> count = {1, ny_, nx_};

      // Read longitude data
      std::vector<double> longitude_data(nx_ * ny_);
      var_longitude.getVar(start, count, longitude_data.data());
      longitude_ =
          xt::reshape_view(xt::adapt(longitude_data, {ny_, nx_}), {ny_, nx_});

      // Read latitude data
      std::vector<double> latitude_data(nx_ * ny_);
      var_latitude.getVar(start, count, latitude_data.data());
      latitude_ =
          xt::reshape_view(xt::adapt(latitude_data, {ny_, nx_}), {ny_, nx_});

      // Read terrain/elevation data
      std::vector<double> elevation_data(nx_ * ny_);
      var_terrain.getVar(start, count, elevation_data.data());
      elevation_ =
          xt::reshape_view(xt::adapt(elevation_data, {ny_, nx_}), {ny_, nx_});

    } else {
      throw std::runtime_error("Failed to open WRF file: " + filename);
    }
  } catch (const netCDF::exceptions::NcException& e) {
    throw std::runtime_error("NetCDF error while loading WRF geometry data: " +
                             std::string(e.what()));
  } catch (const std::exception& e) {
    throw std::runtime_error("Error loading WRF geometry data: " +
                             std::string(e.what()));
  }
}

// Explicit template instantiations for known config backend types
// Add additional instantiations as needed for different config backends
template WRFGeometry::WRFGeometry(
    const metada::backends::yaml::YAMLConfigBackend& config);

}  // namespace metada::backends::wrf