/**
 * @file WRFGeometry.hpp
 * @brief WRF geometry backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <memory>
#include <netcdf>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(_WIN32) || defined(__APPLE__)
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#else
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#endif

// Forward declarations
namespace metada::backends::wrf {
template <typename ConfigBackend>
class WRFGeometryIterator;
template <typename ConfigBackend>
class WRFGeometryConstIterator;
}  // namespace metada::backends::wrf

namespace metada::backends::wrf {

/**
 * @brief WRF geometry backend implementation
 *
 * @details
 * This class implements a geometry backend for the WRF model. It reads
 * geographical data from a WRF NetCDF file and provides methods for
 * traversing the grid and managing boundary conditions.
 */
template <typename ConfigBackend>
class WRFGeometry {
 public:
  // Iterator type aliases
  using iterator = WRFGeometryIterator<ConfigBackend>;
  using const_iterator = WRFGeometryConstIterator<ConfigBackend>;

  /**
   * @brief Default constructor is deleted
   */
  WRFGeometry() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  WRFGeometry(const WRFGeometry&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  WRFGeometry& operator=(const WRFGeometry&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing WRF file path and timestamp
   */
  explicit WRFGeometry(const ConfigBackend& config)
      : wrfFilename_(config.Get("input_file").asString()),
        initialized_(false),
        periodicX_(false),
        periodicY_(false),
        periodicZ_(false) {
    if (wrfFilename_.empty()) {
      throw std::runtime_error(
          "WRF input file path not specified in configuration");
    }

    // Load geometry data from WRF NetCDF file
    loadGeometryData(wrfFilename_);
    initialized_ = true;
  };

  /**
   * @brief Move constructor
   *
   * @param other WRF geometry backend to move from
   */
  WRFGeometry(WRFGeometry&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other WRF geometry backend to move from
   * @return Reference to this geometry after assignment
   */
  WRFGeometry& operator=(WRFGeometry&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~WRFGeometry() = default;

  /**
   * @brief Clone this geometry
   *
   * @return A new WRF geometry backend with the same state
   */
  WRFGeometry clone() const;

  /**
   * @brief Get iterator to the beginning of the grid
   *
   * @return Iterator pointing to the first grid point
   */
  iterator begin();

  /**
   * @brief Get iterator to the end of the grid
   *
   * @return Iterator pointing past the last grid point
   */
  iterator end();

  /**
   * @brief Get const iterator to the beginning of the grid
   *
   * @return Const iterator pointing to the first grid point
   */
  const_iterator begin() const;

  /**
   * @brief Get const iterator to the end of the grid
   *
   * @return Const iterator pointing past the last grid point
   */
  const_iterator end() const;

  /**
   * @brief Get the total number of grid points
   *
   * @return Total number of grid points in the geometry
   */
  std::size_t totalGridSize() const { return nx_ * ny_ * nz_; }

  /**
   * @brief Check if the geometry is periodic in X dimension
   *
   * @return True if periodic in X, false otherwise
   */
  bool isPeriodicX() const;

  /**
   * @brief Check if the geometry is periodic in Y dimension
   *
   * @return True if periodic in Y, false otherwise
   */
  bool isPeriodicY() const;

  /**
   * @brief Check if the geometry is periodic in Z dimension
   *
   * @return True if periodic in Z, false otherwise
   */
  bool isPeriodicZ() const;

  /**
   * @brief Perform halo exchange on a state
   *
   * @param state The state on which to perform halo exchange
   */
  template <typename StateBackend>
  void haloExchange(StateBackend& state);

  /**
   * @brief Halo exchange implementation for the backend
   *
   * @param state_ptr Pointer to the state backend
   */
  void haloExchangeImpl(void* state_ptr);

  /**
   * @brief Check if geometry is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const;

  // WRF specific methods

  /**
   * @brief Get the longitude array
   *
   * @return Reference to the longitude array
   */
  const xt::xarray<double>& getLongitude() const;

  /**
   * @brief Get the latitude array
   *
   * @return Reference to the latitude array
   */
  const xt::xarray<double>& getLatitude() const;

  /**
   * @brief Get the elevation array
   *
   * @return Reference to the elevation array
   */
  const xt::xarray<double>& getElevation() const;

 private:
  /**
   * @brief Load geometry data from the WRF NetCDF file
   *
   * @param filename Path to the WRF NetCDF file
   */
  void loadGeometryData(const std::string& filename) {
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
          // Find matching timestamp (implementation dependent on WRF file
          // format) This is a placeholder for actual timestamp matching logic
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
      throw std::runtime_error(
          "NetCDF error while loading WRF geometry data: " +
          std::string(e.what()));
    } catch (const std::exception& e) {
      throw std::runtime_error("Error loading WRF geometry data: " +
                               std::string(e.what()));
    }
  };

  // Special constructor for cloning
  WRFGeometry(const std::string& fn, bool px, bool py, bool pz);

  // WRF NetCDF file
  std::string wrfFilename_;
  bool initialized_ = false;

  // Geographic data
  xt::xarray<double> longitude_;  // Longitude points
  xt::xarray<double> latitude_;   // Latitude points
  xt::xarray<double> elevation_;  // Elevation/terrain height

  // Grid dimensions
  std::size_t nx_ = 0;  // X dimension (west-east)
  std::size_t ny_ = 0;  // Y dimension (south-north)
  std::size_t nz_ = 0;  // Z dimension (vertical levels)

  // Periodicity flags
  bool periodicX_ = false;
  bool periodicY_ = false;
  bool periodicZ_ = false;

  // Friend declaration for iterator
  friend class WRFGeometryIterator<ConfigBackend>;
};

// Template method implementation
template <typename ConfigBackend>
template <typename StateBackend>
void WRFGeometry<ConfigBackend>::haloExchange(StateBackend& state) {
  haloExchangeImpl(static_cast<void*>(&state));
}

}  // namespace metada::backends::wrf

// Include the iterator implementation first
#include "WRFGeometryIterator.hpp"

// Then define the iterator methods
namespace metada::backends::wrf {

template <typename ConfigBackend>
inline typename WRFGeometry<ConfigBackend>::iterator
WRFGeometry<ConfigBackend>::begin() {
  return iterator(this, 0);
}

template <typename ConfigBackend>
inline typename WRFGeometry<ConfigBackend>::iterator
WRFGeometry<ConfigBackend>::end() {
  return iterator(this, totalGridSize());
}

template <typename ConfigBackend>
inline typename WRFGeometry<ConfigBackend>::const_iterator
WRFGeometry<ConfigBackend>::begin() const {
  return const_iterator(this, 0);
}

template <typename ConfigBackend>
inline typename WRFGeometry<ConfigBackend>::const_iterator
WRFGeometry<ConfigBackend>::end() const {
  return const_iterator(this, totalGridSize());
}

}  // namespace metada::backends::wrf