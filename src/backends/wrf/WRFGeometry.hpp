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
#include "WRFGeometryIterator.hpp"

// --- Begin: GridDimensionInfo struct ---
struct GridDimensionInfo {
  std::string name;  // e.g., "west_east", "west_east_stag"
  bool is_staggered = false;
  size_t size = 0;
  std::vector<double> coordinates;  // 1D coordinate array for this dimension

  // 2D geographic coordinate arrays (for longitude/latitude grids)
  std::vector<double> longitude_2d;  // 2D longitude array (flattened)
  std::vector<double> latitude_2d;   // 2D latitude array (flattened)
  size_t nx = 0;                     // X-dimension size for 2D arrays
  size_t ny = 0;                     // Y-dimension size for 2D arrays

  // Helper methods
  bool has_2d_coords() const {
    return !longitude_2d.empty() && !latitude_2d.empty() && nx > 0 && ny > 0;
  }

  void resize_2d_coords(size_t nx_size, size_t ny_size) {
    nx = nx_size;
    ny = ny_size;
    longitude_2d.resize(nx * ny);
    latitude_2d.resize(nx * ny);
  }
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
  iterator begin() { return iterator(this, 0); }

  /**
   * @brief Get iterator to the end of the grid
   *
   * @return Iterator pointing past the last grid point
   */
  iterator end() { return iterator(this, totalGridSize()); }

  /**
   * @brief Get const iterator to the beginning of the grid
   *
   * @return Const iterator pointing to the first grid point
   */
  const_iterator begin() const { return const_iterator(this, 0); }

  /**
   * @brief Get const iterator to the end of the grid
   *
   * @return Const iterator pointing past the last grid point
   */
  const_iterator end() const { return const_iterator(this, totalGridSize()); }

  // --- Axis/coordinate accessors using GridDimensionInfo ---
  size_t x_dim() const { return unstaggered_grid_.nx; }
  size_t y_dim() const { return unstaggered_grid_.ny; }
  size_t z_dim() const {
    return unstaggered_grid_.size /
           (unstaggered_grid_.nx * unstaggered_grid_.ny);
  }
  size_t x_stag_dim() const { return u_staggered_grid_.nx; }
  size_t y_stag_dim() const { return v_staggered_grid_.ny; }

  // Legacy accessors for compatibility - now return proper dimension info
  const GridDimensionInfo& unstaggered_info() const {
    return unstaggered_grid_;
  }
  const GridDimensionInfo& u_staggered_info() const {
    return u_staggered_grid_;
  }
  const GridDimensionInfo& v_staggered_info() const {
    return v_staggered_grid_;
  }
  // --- End axis/coordinate accessors ---

  // --- Begin: Coordinate grid accessors ---
  const GridDimensionInfo& unstaggered_grid() const {
    return unstaggered_grid_;
  }
  const GridDimensionInfo& u_staggered_grid() const {
    return u_staggered_grid_;
  }
  const GridDimensionInfo& v_staggered_grid() const {
    return v_staggered_grid_;
  }
  // --- End: Coordinate grid accessors ---

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
    const size_t nx = x_dim();
    const size_t ny = y_dim();

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
    if (i >= x_dim() || j >= y_dim()) {
      throw std::out_of_range("Grid coordinates out of range");
    }

    double lon = unstaggered_grid_.longitude_2d[j * x_dim() + i];
    double lat = unstaggered_grid_.latitude_2d[j * x_dim() + i];
    double level = 0.0;  // Default level, could be enhanced to get actual level

    return Location(lat, lon, level);
  }

  /**
   * @brief Get total grid size
   * @return Total number of grid points
   */
  size_t totalGridSize() const { return x_dim() * y_dim() * z_dim(); }

 private:
  /**
   * @brief Load geometry data from the WRF NetCDF file
   *
   * @param filename Path to the WRF NetCDF file
   */
  void loadGeometryData(const std::string& filename);

  // --- Begin: GridDimensionInfo members ---
  GridDimensionInfo unstaggered_grid_;
  GridDimensionInfo u_staggered_grid_;
  GridDimensionInfo v_staggered_grid_;
  // --- End: GridDimensionInfo members ---

  // WRF NetCDF file
  std::string wrfFilename_;
  bool initialized_ = false;

  // Friend declaration for iterator
  friend class WRFGeometryIterator<ConfigBackend>;
};

// Move constructor implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>::WRFGeometry(WRFGeometry&& other) noexcept
    : wrfFilename_(std::move(other.wrfFilename_)),
      initialized_(other.initialized_),
      unstaggered_grid_(std::move(other.unstaggered_grid_)),
      u_staggered_grid_(std::move(other.u_staggered_grid_)),
      v_staggered_grid_(std::move(other.v_staggered_grid_)) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.unstaggered_grid_.size = 0;
  other.u_staggered_grid_.size = 0;
  other.v_staggered_grid_.size = 0;
}

// Move assignment operator implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>& WRFGeometry<ConfigBackend>::operator=(
    WRFGeometry&& other) noexcept {
  if (this != &other) {
    wrfFilename_ = std::move(other.wrfFilename_);
    initialized_ = other.initialized_;
    unstaggered_grid_ = std::move(other.unstaggered_grid_);
    u_staggered_grid_ = std::move(other.u_staggered_grid_);
    v_staggered_grid_ = std::move(other.v_staggered_grid_);

    // Reset the moved-from object
    other.initialized_ = false;
    other.unstaggered_grid_.size = 0;
    other.u_staggered_grid_.size = 0;
    other.v_staggered_grid_.size = 0;
  }
  return *this;
}

// Clone implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend> WRFGeometry<ConfigBackend>::clone() const {
  return *this;
}

// ConfigBackend constructor implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>::WRFGeometry(const ConfigBackend& config)
    : wrfFilename_(config.Get("file").asString()), initialized_(false) {
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
    netCDF::NcFile wrf_file(filename, netCDF::NcFile::read);
    if (!wrf_file.isNull()) {
      // Read dimension information
      auto dim_west_east = wrf_file.getDim("west_east");
      auto dim_south_north = wrf_file.getDim("south_north");
      auto dim_bottom_top = wrf_file.getDim("bottom_top");
      auto dim_west_east_stag = wrf_file.getDim("west_east_stag");
      auto dim_south_north_stag = wrf_file.getDim("south_north_stag");
      auto dim_bottom_top_stag = wrf_file.getDim("bottom_top_stag");

      // Setup unstaggered grid
      unstaggered_grid_.name = "unstaggered";
      unstaggered_grid_.is_staggered = false;
      unstaggered_grid_.size = dim_west_east.getSize() *
                               dim_south_north.getSize() *
                               dim_bottom_top.getSize();
      unstaggered_grid_.resize_2d_coords(dim_west_east.getSize(),
                                         dim_south_north.getSize());

      // Setup U-staggered grid (staggered in X direction)
      u_staggered_grid_.name = "u_staggered";
      u_staggered_grid_.is_staggered = true;
      u_staggered_grid_.size = dim_west_east_stag.getSize() *
                               dim_south_north.getSize() *
                               dim_bottom_top.getSize();
      u_staggered_grid_.resize_2d_coords(dim_west_east_stag.getSize(),
                                         dim_south_north.getSize());

      // Setup V-staggered grid (staggered in Y direction)
      v_staggered_grid_.name = "v_staggered";
      v_staggered_grid_.is_staggered = true;
      v_staggered_grid_.size = dim_west_east.getSize() *
                               dim_south_north_stag.getSize() *
                               dim_bottom_top.getSize();
      v_staggered_grid_.resize_2d_coords(dim_west_east.getSize(),
                                         dim_south_north_stag.getSize());

      // Read unstaggered coordinate arrays (XLONG, XLAT)
      auto var_longitude = wrf_file.getVar("XLONG");
      auto var_latitude = wrf_file.getVar("XLAT");
      if (var_longitude.isNull() || var_latitude.isNull()) {
        throw std::runtime_error("Missing XLONG/XLAT in WRF file");
      }
      std::vector<size_t> start = {0, 0, 0};
      std::vector<size_t> count = {1, unstaggered_grid_.ny,
                                   unstaggered_grid_.nx};
      var_longitude.getVar(start, count, unstaggered_grid_.longitude_2d.data());
      var_latitude.getVar(start, count, unstaggered_grid_.latitude_2d.data());

      // Read U-staggered coordinate arrays (XLONG_U, XLAT_U)
      auto var_longitude_u = wrf_file.getVar("XLONG_U");
      auto var_latitude_u = wrf_file.getVar("XLAT_U");
      if (!var_longitude_u.isNull() && !var_latitude_u.isNull()) {
        std::vector<size_t> count_u = {1, u_staggered_grid_.ny,
                                       u_staggered_grid_.nx};
        var_longitude_u.getVar(start, count_u,
                               u_staggered_grid_.longitude_2d.data());
        var_latitude_u.getVar(start, count_u,
                              u_staggered_grid_.latitude_2d.data());
      }

      // Read V-staggered coordinate arrays (XLONG_V, XLAT_V)
      auto var_longitude_v = wrf_file.getVar("XLONG_V");
      auto var_latitude_v = wrf_file.getVar("XLAT_V");
      if (!var_longitude_v.isNull() && !var_latitude_v.isNull()) {
        std::vector<size_t> count_v = {1, v_staggered_grid_.ny,
                                       v_staggered_grid_.nx};
        var_longitude_v.getVar(start, count_v,
                               v_staggered_grid_.longitude_2d.data());
        var_latitude_v.getVar(start, count_v,
                              v_staggered_grid_.latitude_2d.data());
      }

      // Populate 1D coordinate arrays for dimensional indexing
      // Unstaggered grid coordinates
      unstaggered_grid_.coordinates.resize(unstaggered_grid_.nx);
      for (size_t i = 0; i < unstaggered_grid_.nx; ++i) {
        unstaggered_grid_.coordinates[i] = static_cast<double>(i);
      }

      // U-staggered grid coordinates
      u_staggered_grid_.coordinates.resize(u_staggered_grid_.nx);
      for (size_t i = 0; i < u_staggered_grid_.nx; ++i) {
        u_staggered_grid_.coordinates[i] = static_cast<double>(i) - 0.5;
      }

      // V-staggered grid coordinates
      v_staggered_grid_.coordinates.resize(v_staggered_grid_.nx);
      for (size_t i = 0; i < v_staggered_grid_.nx; ++i) {
        v_staggered_grid_.coordinates[i] = static_cast<double>(i);
      }

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