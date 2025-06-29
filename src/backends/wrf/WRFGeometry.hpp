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

#include "PointObservation.hpp"

// --- Begin: GridDimensionInfo struct ---
struct GridDimensionInfo {
  std::string name;  // e.g., "west_east", "west_east_stag"
  bool is_staggered = false;
  size_t size = 0;
  std::vector<double> coordinates;
};
// --- End: GridDimensionInfo struct ---

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
 *
 * It provides access to all possible coordinate arrays and sizes for both
 * staggered and unstaggered axes, but does not know which variable uses which
 * grid.
 */
template <typename ConfigBackend>
class WRFGeometry {
 public:
  // Iterator type aliases
  using iterator = WRFGeometryIterator<ConfigBackend>;
  using const_iterator = WRFGeometryConstIterator<ConfigBackend>;
  using Location = metada::framework::Location;

  // --- Begin: Additions for GeometryBackendImpl concept compliance ---
  using value_type = Location;
  using reference = Location;
  using const_reference = const Location;
  using pointer = Location*;
  using const_pointer = const Location*;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  // Iteration
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }

  // Size information
  size_type size() const { return totalGridSize(); }
  bool empty() const { return size() == 0; }
  size_type max_size() const { return size(); }

  // Element access
  reference operator[](size_type idx) { return getLocation(idx); }
  const_reference operator[](size_type idx) const { return getLocation(idx); }
  reference at(size_type idx) { return getLocation(idx); }
  const_reference at(size_type idx) const { return getLocation(idx); }
  reference front() { return getLocation(0); }
  const_reference front() const { return getLocation(0); }
  reference back() { return getLocation(size() - 1); }
  const_reference back() const { return getLocation(size() - 1); }
  // --- End: Additions for GeometryBackendImpl concept compliance ---

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
  explicit WRFGeometry(const ConfigBackend& config);

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

  // --- Axis/coordinate accessors using GridDimensionInfo ---
  const GridDimensionInfo& x_info() const { return x_info_; }
  const GridDimensionInfo& y_info() const { return y_info_; }
  const GridDimensionInfo& z_info() const { return z_info_; }
  const GridDimensionInfo& x_stag_info() const { return x_stag_info_; }
  const GridDimensionInfo& y_stag_info() const { return y_stag_info_; }
  const GridDimensionInfo& z_stag_info() const { return z_stag_info_; }
  // --- End axis/coordinate accessors ---

  /**
   * @brief Get grid point as Location object
   * @param i X coordinate (west-east)
   * @param j Y coordinate (south-north)
   * @param k Z coordinate (bottom-top)
   * @return Location object with grid coordinates
   */
  Location getLocation(size_t i, size_t j, size_t k = 0) const {
    return Location(static_cast<int>(i), static_cast<int>(j),
                    static_cast<int>(k));
  }

  /**
   * @brief Get grid point as Location object from linear index
   * @param index Linear index into the grid
   * @return Location object with grid coordinates
   */
  Location getLocation(size_t index) const {
    if (index >= totalGridSize()) {
      throw std::out_of_range("Grid index out of range");
    }

    // Calculate 3D indices from linear index
    const size_t nx = x_info_.size;
    const size_t ny = y_info_.size;

    // Using row-major order: index = k*nx*ny + j*nx + i
    size_t k = index / (nx * ny);
    const size_t remainder = index % (nx * ny);
    size_t j = remainder / nx;
    size_t i = remainder % nx;

    return getLocation(i, j, k);
  }

  /**
   * @brief Get geographic location as Location object
   * @param i X coordinate (west-east)
   * @param j Y coordinate (south-north)
   * @return Location object with geographic coordinates
   */
  Location getGeographicLocation(size_t i, size_t j) const {
    if (i >= x_info_.size || j >= y_info_.size) {
      throw std::out_of_range("Grid coordinates out of range");
    }

    double lon = x_info_.coordinates[i];
    double lat = y_info_.coordinates[j];
    double level = 0.0;  // Default level, could be enhanced to get actual level

    return Location(lat, lon, level);
  }

  /**
   * @brief Get total grid size
   * @return Total number of grid points
   */
  size_t totalGridSize() const {
    return x_info_.size * y_info_.size * z_info_.size;
  }

  size_t x_dim() const { return x_info_.size; }
  size_t y_dim() const { return y_info_.size; }

 private:
  /**
   * @brief Load geometry data from the WRF NetCDF file
   *
   * @param filename Path to the WRF NetCDF file
   */
  void loadGeometryData(const std::string& filename);

  // Special constructor for cloning
  WRFGeometry(const std::string& fn, bool px, bool py, bool pz);

  // --- Begin: GridDimensionInfo members ---
  GridDimensionInfo x_info_, y_info_, z_info_;
  GridDimensionInfo x_stag_info_, y_stag_info_, z_stag_info_;
  // --- End: GridDimensionInfo members ---

  // WRF NetCDF file
  std::string wrfFilename_;
  bool initialized_ = false;

  // Periodicity flags
  bool periodicX_ = false;
  bool periodicY_ = false;
  bool periodicZ_ = false;

  // Friend declaration for iterator
  friend class WRFGeometryIterator<ConfigBackend>;
};

// Move constructor implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>::WRFGeometry(WRFGeometry&& other) noexcept
    : wrfFilename_(std::move(other.wrfFilename_)),
      initialized_(other.initialized_),
      x_info_(std::move(other.x_info_)),
      y_info_(std::move(other.y_info_)),
      z_info_(std::move(other.z_info_)),
      x_stag_info_(std::move(other.x_stag_info_)),
      y_stag_info_(std::move(other.y_stag_info_)),
      z_stag_info_(std::move(other.z_stag_info_)),
      periodicX_(other.periodicX_),
      periodicY_(other.periodicY_),
      periodicZ_(other.periodicZ_) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.x_info_.size = 0;
  other.y_info_.size = 0;
  other.z_info_.size = 0;
}

// Move assignment operator implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>& WRFGeometry<ConfigBackend>::operator=(
    WRFGeometry&& other) noexcept {
  if (this != &other) {
    wrfFilename_ = std::move(other.wrfFilename_);
    initialized_ = other.initialized_;
    x_info_ = std::move(other.x_info_);
    y_info_ = std::move(other.y_info_);
    z_info_ = std::move(other.z_info_);
    x_stag_info_ = std::move(other.x_stag_info_);
    y_stag_info_ = std::move(other.y_stag_info_);
    z_stag_info_ = std::move(other.z_stag_info_);
    periodicX_ = other.periodicX_;
    periodicY_ = other.periodicY_;
    periodicZ_ = other.periodicZ_;

    // Reset the moved-from object
    other.initialized_ = false;
    other.x_info_.size = 0;
    other.y_info_.size = 0;
    other.z_info_.size = 0;
  }
  return *this;
}

// Clone implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend> WRFGeometry<ConfigBackend>::clone() const {
  // Use the specialized constructor
  WRFGeometry<ConfigBackend> clone(wrfFilename_, periodicX_, periodicY_,
                                   periodicZ_);

  // Copy data members
  clone.x_info_ = x_info_;
  clone.y_info_ = y_info_;
  clone.z_info_ = z_info_;
  clone.x_stag_info_ = x_stag_info_;
  clone.y_stag_info_ = y_stag_info_;
  clone.z_stag_info_ = z_stag_info_;
  clone.initialized_ = initialized_;

  return clone;
}

// ConfigBackend constructor implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>::WRFGeometry(const ConfigBackend& config)
    : wrfFilename_(config.Get("file").asString()),
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
}

// Implementation of loadGeometryData
template <typename ConfigBackend>
void WRFGeometry<ConfigBackend>::loadGeometryData(const std::string& filename) {
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

      x_info_.size = dim_west_east.getSize();
      y_info_.size = dim_south_north.getSize();
      z_info_.size = dim_bottom_top.getSize();

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
      std::vector<size_t> count = {1, y_info_.size, x_info_.size};

      // Read longitude data
      std::vector<double> longitude_data(x_info_.size * y_info_.size);
      var_longitude.getVar(start, count, longitude_data.data());
      x_info_.coordinates = longitude_data;

      // Read latitude data
      std::vector<double> latitude_data(x_info_.size * y_info_.size);
      var_latitude.getVar(start, count, latitude_data.data());
      y_info_.coordinates = latitude_data;

      // Read terrain/elevation data
      std::vector<double> elevation_data(x_info_.size * y_info_.size);
      var_terrain.getVar(start, count, elevation_data.data());
      z_info_.coordinates = elevation_data;

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
  return iterator(this, x_info_.size * y_info_.size * z_info_.size);
}

template <typename ConfigBackend>
inline typename WRFGeometry<ConfigBackend>::const_iterator
WRFGeometry<ConfigBackend>::begin() const {
  return const_iterator(this, 0);
}

template <typename ConfigBackend>
inline typename WRFGeometry<ConfigBackend>::const_iterator
WRFGeometry<ConfigBackend>::end() const {
  return const_iterator(this, x_info_.size * y_info_.size * z_info_.size);
}

}  // namespace metada::backends::wrf