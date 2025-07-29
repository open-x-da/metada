/**
 * @file WRFGeometry.hpp
 * @brief WRF geometry backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <netcdf>
#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>

#include "Location.hpp"
#include "WRFGeometryIterator.hpp"

/**
 * @brief Grid dimension and coordinate information for WRF grids
 *
 * @details This structure stores comprehensive information about a specific
 * WRF grid type, including dimensions, coordinate arrays, and utility methods
 * for accessing grid properties. It supports both staggered and unstaggered
 * grids used in WRF's Arakawa-C grid system.
 *
 * The structure handles:
 * - 2D geographic coordinate arrays (longitude/latitude)
 * - Vertical coordinate arrays (eta levels)
 * - Grid dimension information and metadata
 * - Utility methods for coordinate validation and manipulation
 */
struct GridDimensionInfo {
  ///@{ @name Grid Metadata
  std::string
      name;  ///< Grid name identifier (e.g., "unstaggered", "u_staggered")
  bool is_staggered = false;  ///< Whether this grid uses staggered coordinates
  size_t size = 0;            ///< Total 3D grid size (nx * ny * nz)
  ///@}

  ///@{ @name Horizontal Grid Information
  size_t nx = 0;  ///< X-dimension size (west-east)
  size_t ny = 0;  ///< Y-dimension size (south-north)
  ///@}

  ///@{ @name Geographic Coordinate Arrays
  std::vector<double>
      longitude_2d;                 ///< 2D longitude array (flattened, degrees)
  std::vector<double> latitude_2d;  ///< 2D latitude array (flattened, degrees)
  ///@}

  ///@{ @name Vertical Grid Information
  size_t nz = 0;  ///< Z-dimension size (bottom-top)
  std::vector<double>
      vertical_coords;  ///< Vertical coordinates (ZNU or ZNW eta levels)
  ///@}

  ///@{ @name Coordinate Validation Methods
  /**
   * @brief Check if 2D geographic coordinates are available
   *
   * @details Verifies that the grid has valid 2D coordinate arrays
   * for geographic transformations. This includes both longitude
   * and latitude arrays with proper dimensions.
   *
   * @return bool True if 2D coordinates are available and valid
   */
  bool has_2d_coords() const {
    return !longitude_2d.empty() && !latitude_2d.empty() && nx > 0 && ny > 0;
  }

  /**
   * @brief Check if vertical coordinates are available
   *
   * @details Verifies that the grid has valid vertical coordinate
   * information for 3D interpolations and level calculations.
   *
   * @return bool True if vertical coordinates are available and valid
   */
  bool has_vertical_coords() const {
    return !vertical_coords.empty() && nz > 0;
  }
  ///@}

  ///@{ @name Grid Resize Methods
  /**
   * @brief Resize 2D coordinate arrays for new grid dimensions
   *
   * @details Allocates memory for 2D coordinate arrays based on the
   * specified horizontal grid dimensions. This method should be called
   * before loading coordinate data from NetCDF files.
   *
   * @param[in] nx_size New X-dimension size (west-east)
   * @param[in] ny_size New Y-dimension size (south-north)
   */
  void resize_2d_coords(size_t nx_size, size_t ny_size) {
    nx = nx_size;
    ny = ny_size;
    longitude_2d.resize(nx * ny);
    latitude_2d.resize(nx * ny);
  }

  /**
   * @brief Resize vertical coordinate array for new grid dimensions
   *
   * @details Allocates memory for vertical coordinate arrays based on
   * the specified vertical grid dimension. This method should be called
   * before loading vertical level data from NetCDF files.
   *
   * @param[in] nz_size New Z-dimension size (bottom-top)
   */
  void resize_vertical_coords(size_t nz_size) {
    nz = nz_size;
    vertical_coords.resize(nz);
  }
  ///@}
};

// Forward declarations
namespace metada::backends::wrf {
template <typename ConfigBackend>
class WRFGeometryIterator;
}  // namespace metada::backends::wrf

namespace metada::backends::wrf {

using framework::CoordinateSystem;

/**
 * @brief WRF geometry backend implementation for meteorological grids
 *
 * @details This class implements a comprehensive geometry backend for the
 * Weather Research and Forecasting (WRF) model. It manages multiple staggered
 * grid systems used in WRF's Arakawa-C grid configuration and provides access
 * to geographic coordinates, grid dimensions, and iteration capabilities.
 *
 * Key features:
 * - Support for multiple grid types (unstaggered, U-staggered, V-staggered,
 * W-staggered)
 * - Geographic coordinate transformations (latitude/longitude)
 * - Vertical coordinate handling (eta levels)
 * - STL-compatible container interface with iterators
 * - NetCDF file reading for WRF grid data
 * - Grid dimension validation and coordinate interpolation
 *
 * Grid Types Supported:
 * - Unstaggered grid: Mass points (T, P, RH, etc.)
 * - U-staggered grid: U-wind component (staggered in X direction)
 * - V-staggered grid: V-wind component (staggered in Y direction)
 * - W-staggered grid: W-wind component (staggered in Z direction)
 *
 * The class provides efficient iteration over grid points and supports
 * both linear indexing and 3D coordinate access patterns commonly used
 * in meteorological applications.
 *
 * @tparam ConfigBackend Configuration backend type providing file paths and
 * settings
 *
 * @see GridDimensionInfo
 * @see WRFGeometryIterator
 * @see framework::Location
 */
template <typename ConfigBackend>
class WRFGeometry {
 public:
  ///@{ @name Type Definitions
  /**
   * @brief Iterator and location type aliases
   *
   * These type definitions provide convenient access to iterator types
   * and location objects used throughout the geometry interface.
   */
  using iterator =
      WRFGeometryIterator<ConfigBackend>;  ///< Mutable iterator type
  using const_iterator =
      WRFGeometryIterator<ConfigBackend>;  ///< Const iterator type
  using Location = framework::Location;    ///< Location object type

  /**
   * @brief Standard container type aliases for STL compatibility
   *
   * These type definitions enable the WRFGeometry to be used with STL
   * algorithms and provide a consistent interface following standard container
   * conventions.
   */
  using value_type = Location;  ///< Type of elements (Location objects)
  using reference = Location;   ///< Reference to element type
  using const_reference = const Location;  ///< Const reference to element type
  using pointer = Location*;               ///< Pointer to element type
  using const_pointer = const Location*;   ///< Const pointer to element type
  using size_type = std::size_t;           ///< Type for sizes and indices
  using difference_type = std::ptrdiff_t;  ///< Type for iterator differences
  ///@}

  ///@{ @name Container Interface Methods
  /**
   * @brief Get const iterator to beginning (alternative interface)
   * @return const_iterator Const iterator to first grid point
   */
  const_iterator cbegin() const { return begin(); }

  /**
   * @brief Get const iterator to end (alternative interface)
   * @return const_iterator Const iterator past last grid point
   */
  const_iterator cend() const { return end(); }

  /**
   * @brief Get total number of grid points
   * @return size_type Total grid size (same as totalGridSize())
   */
  size_type size() const { return totalGridSize(); }

  /**
   * @brief Check if geometry contains no grid points
   * @return bool True if grid is empty, false otherwise
   */
  bool empty() const { return size() == 0; }

  /**
   * @brief Get maximum possible grid size
   * @return size_type Maximum size (same as current size for WRF grids)
   */
  size_type max_size() const { return size(); }

  /**
   * @brief Access grid point by linear index (unchecked)
   * @param[in] idx Linear index into the grid
   * @return reference Location object at the specified index
   */
  reference operator[](size_type idx) { return getLocation(idx); }

  /**
   * @brief Access grid point by linear index (unchecked, const)
   * @param[in] idx Linear index into the grid
   * @return const_reference Const Location object at the specified index
   */
  const_reference operator[](size_type idx) const { return getLocation(idx); }

  /**
   * @brief Access grid point by linear index (bounds-checked)
   * @param[in] idx Linear index into the grid
   * @return reference Location object at the specified index
   * @throws std::out_of_range If index is out of bounds
   */
  reference at(size_type idx) { return getLocation(idx); }

  /**
   * @brief Access grid point by linear index (bounds-checked, const)
   * @param[in] idx Linear index into the grid
   * @return const_reference Const Location object at the specified index
   * @throws std::out_of_range If index is out of bounds
   */
  const_reference at(size_type idx) const { return getLocation(idx); }

  /**
   * @brief Get first grid point
   * @return reference Location object at the first grid point
   */
  reference front() { return getLocation(0); }

  /**
   * @brief Get first grid point (const)
   * @return const_reference Const Location object at the first grid point
   */
  const_reference front() const { return getLocation(0); }

  /**
   * @brief Get last grid point
   * @return reference Location object at the last grid point
   */
  reference back() { return getLocation(size() - 1); }

  /**
   * @brief Get last grid point (const)
   * @return const_reference Const Location object at the last grid point
   */
  const_reference back() const { return getLocation(size() - 1); }
  ///@}

  ///@{ @name Construction and Destruction
  /**
   * @brief Default constructor is deleted
   *
   * @details Default construction is disabled because WRF geometry requires
   * a valid configuration specifying the NetCDF file path and grid parameters.
   */
  WRFGeometry() = delete;

  /**
   * @brief Copy constructor is deleted
   *
   * @details Copy construction is disabled to prevent expensive copying of
   * large coordinate arrays. Use move semantics or create new instances
   * from configuration instead.
   */
  WRFGeometry(const WRFGeometry&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   *
   * @details Copy assignment is disabled to prevent expensive copying of
   * large coordinate arrays. Use move semantics instead.
   */
  WRFGeometry& operator=(const WRFGeometry&) = delete;

  /**
   * @brief Construct WRF geometry from configuration
   *
   * @details Creates a new WRF geometry by reading grid information from
   * the NetCDF file specified in the configuration. The constructor loads
   * all grid types (unstaggered, U-staggered, V-staggered, W-staggered)
   * and their associated coordinate arrays.
   *
   * @param[in] config Configuration containing WRF file path and settings
   * @throws std::runtime_error If file path is missing or file cannot be read
   * @throws netCDF::exceptions::NcException If NetCDF operations fail
   *
   * @note The configuration must contain a "file" key with the WRF NetCDF file
   * path
   */
  explicit WRFGeometry(const ConfigBackend& config);

  /**
   * @brief Move constructor for efficient transfer
   *
   * @details Transfers ownership of all internal data structures without
   * copying large coordinate arrays. The moved-from object is reset to
   * an uninitialized state.
   *
   * @param[in] other WRF geometry backend to move from
   */
  WRFGeometry(WRFGeometry&& other) noexcept;

  /**
   * @brief Move assignment operator for efficient transfer
   *
   * @details Transfers ownership of all internal data structures without
   * copying large coordinate arrays. The moved-from object is reset to
   * an uninitialized state.
   *
   * @param[in] other WRF geometry backend to move from
   * @return WRFGeometry& Reference to this geometry after assignment
   */
  WRFGeometry& operator=(WRFGeometry&& other) noexcept;

  /**
   * @brief Destructor with automatic cleanup
   *
   * @details Automatically releases all allocated memory for coordinate
   * arrays and grid information. No explicit cleanup is required.
   */
  ~WRFGeometry() = default;
  ///@}

  ///@{ @name Geometry Cloning
  /**
   * @brief Clone this geometry (not implemented)
   *
   * @details Cloning is not implemented for WRF geometry due to the complexity
   * of duplicating NetCDF file connections and large coordinate arrays.
   * Create new instances from configuration instead.
   *
   * @return WRFGeometry A new WRF geometry backend with the same state
   * @throws std::runtime_error Always throws as cloning is not supported
   *
   * @note Use the configuration constructor to create equivalent geometries
   */
  WRFGeometry clone() const;
  ///@}

  ///@{ @name Grid Iteration
  /**
   * @brief Get iterator to the beginning of the grid
   *
   * @details Returns a mutable iterator pointing to the first grid point.
   * The iterator traverses the grid in row-major order: X (west-east) first,
   * then Y (south-north), then Z (bottom-top).
   *
   * @return iterator Iterator pointing to the first grid point
   */
  iterator begin() { return iterator(this, 0); }

  /**
   * @brief Get iterator past the end of the grid
   *
   * @details Returns a mutable iterator pointing past the last grid point.
   * This iterator serves as the termination condition for range-based loops.
   *
   * @return iterator Iterator pointing past the last grid point
   */
  iterator end() { return iterator(this, totalGridSize()); }

  /**
   * @brief Get const iterator to the beginning of the grid
   *
   * @details Returns a const iterator pointing to the first grid point.
   * The iterator traverses the grid in row-major order: X (west-east) first,
   * then Y (south-north), then Z (bottom-top).
   *
   * @return const_iterator Const iterator pointing to the first grid point
   */
  const_iterator begin() const { return const_iterator(this, 0); }

  /**
   * @brief Get const iterator past the end of the grid
   *
   * @details Returns a const iterator pointing past the last grid point.
   * This iterator serves as the termination condition for range-based loops.
   *
   * @return const_iterator Const iterator pointing past the last grid point
   */
  const_iterator end() const { return const_iterator(this, totalGridSize()); }
  ///@}

  ///@{ @name Grid Dimension Access
  /**
   * @brief Get X-dimension size for unstaggered grid
   *
   * @details Returns the number of grid points in the west-east direction
   * for the unstaggered (mass point) grid.
   *
   * @return size_t X-dimension size (west-east)
   */
  size_t x_dim() const { return unstaggered_grid_.nx; }

  /**
   * @brief Get Y-dimension size for unstaggered grid
   *
   * @details Returns the number of grid points in the south-north direction
   * for the unstaggered (mass point) grid.
   *
   * @return size_t Y-dimension size (south-north)
   */
  size_t y_dim() const { return unstaggered_grid_.ny; }

  /**
   * @brief Get Z-dimension size for unstaggered grid
   *
   * @details Returns the number of grid points in the bottom-top direction
   * for the unstaggered (mass point) grid.
   *
   * @return size_t Z-dimension size (bottom-top)
   */
  size_t z_dim() const {
    return unstaggered_grid_.size /
           (unstaggered_grid_.nx * unstaggered_grid_.ny);
  }

  /**
   * @brief Get X-dimension size for U-staggered grid
   *
   * @details Returns the number of grid points in the west-east direction
   * for the U-staggered grid (U-wind component).
   *
   * @return size_t X-dimension size for U-staggered grid
   */
  size_t x_stag_dim() const { return u_staggered_grid_.nx; }

  /**
   * @brief Get Y-dimension size for V-staggered grid
   *
   * @details Returns the number of grid points in the south-north direction
   * for the V-staggered grid (V-wind component).
   *
   * @return size_t Y-dimension size for V-staggered grid
   */
  size_t y_stag_dim() const { return v_staggered_grid_.ny; }

  /**
   * @brief Get Z-dimension size for W-staggered grid
   *
   * @details Returns the number of grid points in the bottom-top direction
   * for the W-staggered grid (W-wind component).
   *
   * @return size_t Z-dimension size for W-staggered grid
   */
  size_t z_stag_dim() const {
    return w_staggered_grid_.size /
           (w_staggered_grid_.nx * w_staggered_grid_.ny);
  }
  ///@}

  ///@{ @name Grid Information Access
  /**
   * @brief Get unstaggered grid information
   *
   * @details Returns a reference to the GridDimensionInfo structure containing
   * complete information about the unstaggered (mass point) grid, including
   * dimensions, coordinate arrays, and metadata.
   *
   * @return const GridDimensionInfo& Reference to unstaggered grid info
   */
  const GridDimensionInfo& unstaggered_info() const {
    return unstaggered_grid_;
  }

  /**
   * @brief Get U-staggered grid information
   *
   * @details Returns a reference to the GridDimensionInfo structure containing
   * complete information about the U-staggered grid (U-wind component),
   * including dimensions, coordinate arrays, and metadata.
   *
   * @return const GridDimensionInfo& Reference to U-staggered grid info
   */
  const GridDimensionInfo& u_staggered_info() const {
    return u_staggered_grid_;
  }

  /**
   * @brief Get V-staggered grid information
   *
   * @details Returns a reference to the GridDimensionInfo structure containing
   * complete information about the V-staggered grid (V-wind component),
   * including dimensions, coordinate arrays, and metadata.
   *
   * @return const GridDimensionInfo& Reference to V-staggered grid info
   */
  const GridDimensionInfo& v_staggered_info() const {
    return v_staggered_grid_;
  }

  /**
   * @brief Get W-staggered grid information
   *
   * @details Returns a reference to the GridDimensionInfo structure containing
   * complete information about the W-staggered grid (W-wind component),
   * including dimensions, coordinate arrays, and metadata.
   *
   * @return const GridDimensionInfo& Reference to W-staggered grid info
   */
  const GridDimensionInfo& w_staggered_info() const {
    return w_staggered_grid_;
  }
  ///@}

  ///@{ @name Coordinate Access Methods
  /**
   * @brief Get geographic coordinates for a 3D grid point in the unstaggered
   * grid
   *
   * @details Converts 3D grid indices to geographic coordinates using the
   * unstaggered grid's coordinate arrays. Returns latitude, longitude in
   * degrees and vertical level (eta coordinate).
   *
   * @param[in] i X coordinate index (west-east, 0-based)
   * @param[in] j Y coordinate index (south-north, 0-based)
   * @param[in] k Z coordinate index (bottom-top, 0-based)
   * @return std::tuple<double, double, double> Tuple of (latitude, longitude,
   * level)
   * @throws std::out_of_range If grid coordinates are out of bounds
   *
   * @note Returns (0.0, 0.0, 0.0) if geographic coordinates are unavailable
   * @note This method specifically works with the unstaggered grid
   */
  std::tuple<double, double, double> getGeographicCoords(size_t i, size_t j,
                                                         size_t k = 0) const {
    if (i >= x_dim() || j >= y_dim() || k >= z_dim()) {
      throw std::out_of_range("Grid coordinates out of range");
    }

    // Get geographic coordinates from the unstaggered grid
    if (unstaggered_grid_.has_2d_coords()) {
      size_t idx = j * x_dim() + i;
      if (idx < unstaggered_grid_.longitude_2d.size() &&
          idx < unstaggered_grid_.latitude_2d.size()) {
        double lon = unstaggered_grid_.longitude_2d[idx];
        double lat = unstaggered_grid_.latitude_2d[idx];
        double level = 0.0;

        // Get vertical level if available
        if (unstaggered_grid_.has_vertical_coords() &&
            k < unstaggered_grid_.vertical_coords.size()) {
          level = unstaggered_grid_.vertical_coords[k];
        }

        return std::make_tuple(lat, lon, level);
      }
    }

    // Return default values if geographic coordinates are not available
    return std::make_tuple(0.0, 0.0, 0.0);
  }

  /**
   * @brief Get Location object for a 3D grid point in the unstaggered grid
   *
   * @details Converts 3D grid indices to a Location object containing
   * geographic coordinates from the unstaggered (mass point) grid. Prefers
   * geographic coordinates when available, falls back to grid coordinates
   * otherwise.
   *
   * @param[in] i X coordinate index (west-east, 0-based)
   * @param[in] j Y coordinate index (south-north, 0-based)
   * @param[in] k Z coordinate index (bottom-top, 0-based)
   * @return Location Location object with coordinates
   * @throws std::out_of_range If grid coordinates are out of bounds
   *
   * @note Uses GEOGRAPHIC coordinate system when possible, GRID otherwise
   * @note This method specifically works with the unstaggered grid
   */
  Location getLocation(size_t i, size_t j, size_t k = 0) const {
    if (i >= x_dim() || j >= y_dim() || k >= z_dim()) {
      throw std::out_of_range("Grid coordinates out of range");
    }

    // Get geographic coordinates from the unstaggered grid
    if (unstaggered_grid_.has_2d_coords()) {
      size_t idx = j * x_dim() + i;
      if (idx < unstaggered_grid_.longitude_2d.size() &&
          idx < unstaggered_grid_.latitude_2d.size()) {
        double lon = unstaggered_grid_.longitude_2d[idx];
        double lat = unstaggered_grid_.latitude_2d[idx];
        double level = 0.0;

        // Get vertical level if available
        if (unstaggered_grid_.has_vertical_coords() &&
            k < unstaggered_grid_.vertical_coords.size()) {
          level = unstaggered_grid_.vertical_coords[k];
        }

        return Location(lat, lon, level, CoordinateSystem::GEOGRAPHIC);
      }
    }

    // Fallback to grid coordinates if geographic coordinates are not available
    return Location(static_cast<int>(i), static_cast<int>(j),
                    static_cast<int>(k));
  }

  /**
   * @brief Get Location object from linear grid index across all grids
   *
   * @details Converts a linear grid index across all grids to a Location
   * object. The method first determines which grid type the index belongs to,
   * then calculates the appropriate 3D grid indices and obtains geographic
   * coordinates. Grid order: unstaggered, U-staggered, V-staggered,
   * W-staggered.
   *
   * @param[in] global_index Linear index across all grids (0-based)
   * @return Location Location object with coordinates
   * @throws std::out_of_range If linear index is out of bounds
   */
  Location getLocation(size_t global_index) const {
    auto [grid_type, local_index] = getGridTypeAndIndex(global_index);

    // Select the appropriate grid based on grid type
    const GridDimensionInfo* grid_info;
    switch (grid_type) {
      case 0:
        grid_info = &unstaggered_grid_;
        break;
      case 1:
        grid_info = &u_staggered_grid_;
        break;
      case 2:
        grid_info = &v_staggered_grid_;
        break;
      case 3:
        grid_info = &w_staggered_grid_;
        break;
      default:
        throw std::runtime_error("Invalid grid type");
    }

    // Calculate 3D indices from local index
    const size_t nx = grid_info->nx;
    const size_t ny = grid_info->ny;
    const size_t nz = grid_info->nz;

    if (nx == 0 || ny == 0 || nz == 0) {
      throw std::runtime_error("Grid dimensions not properly initialized");
    }

    // Using row-major order: index = k*nx*ny + j*nx + i
    size_t k = local_index / (nx * ny);
    const size_t remainder = local_index % (nx * ny);
    size_t j = remainder / nx;
    size_t i = remainder % nx;

    // Get geographic coordinates if available
    if (grid_info->has_2d_coords()) {
      size_t coord_idx = j * nx + i;
      if (coord_idx < grid_info->longitude_2d.size() &&
          coord_idx < grid_info->latitude_2d.size()) {
        double lon = grid_info->longitude_2d[coord_idx];
        double lat = grid_info->latitude_2d[coord_idx];
        double level = 0.0;

        // Get vertical level if available
        if (grid_info->has_vertical_coords() &&
            k < grid_info->vertical_coords.size()) {
          level = grid_info->vertical_coords[k];
        }

        return Location(lat, lon, level, CoordinateSystem::GEOGRAPHIC);
      }
    }

    // Fallback to grid coordinates with grid type information
    return Location(static_cast<int>(i), static_cast<int>(j),
                    static_cast<int>(k), CoordinateSystem::GRID);
  }
  ///@}

  ///@{ @name Grid Size Information
  /**
   * @brief Get total number of grid points across all grids
   *
   * @details Calculates the total size of all WRF grids by summing the
   * unstaggered, U-staggered, V-staggered, and W-staggered grid sizes.
   * This represents the total number of grid points across all grid types
   * in the WRF grid system.
   *
   * @return size_t Total number of grid points across all grids
   */
  size_t totalGridSize() const {
    return unstaggered_grid_.size + u_staggered_grid_.size +
           v_staggered_grid_.size + w_staggered_grid_.size;
  }

  /**
   * @brief Get size of unstaggered grid only
   *
   * @details Returns the size of just the unstaggered (mass point) grid.
   * This is useful when you need to know the size of a specific grid type.
   *
   * @return size_t Size of unstaggered grid
   */
  size_t unstaggeredGridSize() const { return unstaggered_grid_.size; }

  /**
   * @brief Get size of U-staggered grid only
   *
   * @details Returns the size of just the U-staggered grid.
   *
   * @return size_t Size of U-staggered grid
   */
  size_t uStaggeredGridSize() const { return u_staggered_grid_.size; }

  /**
   * @brief Get size of V-staggered grid only
   *
   * @details Returns the size of just the V-staggered grid.
   *
   * @return size_t Size of V-staggered grid
   */
  size_t vStaggeredGridSize() const { return v_staggered_grid_.size; }

  /**
   * @brief Get size of W-staggered grid only
   *
   * @details Returns the size of just the W-staggered grid.
   *
   * @return size_t Size of W-staggered grid
   */
  size_t wStaggeredGridSize() const { return w_staggered_grid_.size; }
  ///@}

  ///@{ @name Grid Type Identification
  /**
   * @brief Determine which grid type a linear index belongs to
   *
   * @details Given a linear index across all grids, determines which
   * grid type it belongs to and returns the local index within that grid.
   *
   * @param[in] global_index Linear index across all grids
   * @return std::pair<int, size_t> Pair of (grid_type, local_index)
   *         where grid_type: 0=unstaggered, 1=U-staggered, 2=V-staggered,
   * 3=W-staggered
   * @throws std::out_of_range If global_index is out of bounds
   */
  std::pair<int, size_t> getGridTypeAndIndex(size_t global_index) const {
    if (global_index >= totalGridSize()) {
      throw std::out_of_range("Global grid index out of range");
    }

    // Check unstaggered grid
    if (global_index < unstaggered_grid_.size) {
      return {0, global_index};
    }

    // Check U-staggered grid
    size_t offset = unstaggered_grid_.size;
    if (global_index < offset + u_staggered_grid_.size) {
      return {1, global_index - offset};
    }

    // Check V-staggered grid
    offset += u_staggered_grid_.size;
    if (global_index < offset + v_staggered_grid_.size) {
      return {2, global_index - offset};
    }

    // Must be W-staggered grid
    offset += v_staggered_grid_.size;
    return {3, global_index - offset};
  }
  ///@}

 private:
  ///@{ @name Private Implementation Methods
  /**
   * @brief Load geometry data from WRF NetCDF file
   *
   * @details Reads grid dimensions, coordinate arrays, and metadata from
   * the specified WRF NetCDF file. This method initializes all four grid
   * types (unstaggered, U-staggered, V-staggered, W-staggered) and loads
   * their associated coordinate information.
   *
   * The method performs the following operations:
   * - Reads dimension information from NetCDF file
   * - Allocates memory for coordinate arrays
   * - Loads 2D geographic coordinate arrays (XLONG, XLAT, etc.)
   * - Loads vertical coordinate arrays (ZNU, ZNW)
   * - Validates loaded data for consistency
   *
   * @param[in] filename Path to the WRF NetCDF file
   * @throws std::runtime_error If file cannot be opened or read
   * @throws netCDF::exceptions::NcException If NetCDF operations fail
   *
   * @note This method is called during construction
   */
  void loadGeometryData(const std::string& filename);
  ///@}

  ///@{ @name Grid Data Storage
  /**
   * @brief Grid dimension and coordinate information for all WRF grid types
   *
   * These structures store comprehensive information about each grid type
   * used in WRF's Arakawa-C staggered grid system, including dimensions,
   * coordinate arrays, and metadata.
   */
  GridDimensionInfo unstaggered_grid_;  ///< Mass point grid (T, P, RH, etc.)
  GridDimensionInfo u_staggered_grid_;  ///< U-wind grid (X-staggered)
  GridDimensionInfo v_staggered_grid_;  ///< V-wind grid (Y-staggered)
  GridDimensionInfo w_staggered_grid_;  ///< W-wind grid (Z-staggered)
  ///@}

  ///@{ @name File and State Information
  std::string wrfFilename_;   ///< Path to source WRF NetCDF file
  bool initialized_ = false;  ///< Initialization status flag
  ///@}
};

// Move constructor implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend>::WRFGeometry(WRFGeometry&& other) noexcept
    : wrfFilename_(std::move(other.wrfFilename_)),
      initialized_(other.initialized_),
      unstaggered_grid_(std::move(other.unstaggered_grid_)),
      u_staggered_grid_(std::move(other.u_staggered_grid_)),
      v_staggered_grid_(std::move(other.v_staggered_grid_)),
      w_staggered_grid_(std::move(other.w_staggered_grid_)) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.unstaggered_grid_.size = 0;
  other.u_staggered_grid_.size = 0;
  other.v_staggered_grid_.size = 0;
  other.w_staggered_grid_.size = 0;
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
    w_staggered_grid_ = std::move(other.w_staggered_grid_);

    // Reset the moved-from object
    other.initialized_ = false;
    other.unstaggered_grid_.size = 0;
    other.u_staggered_grid_.size = 0;
    other.v_staggered_grid_.size = 0;
    other.w_staggered_grid_.size = 0;
  }
  return *this;
}

// Clone implementation
template <typename ConfigBackend>
WRFGeometry<ConfigBackend> WRFGeometry<ConfigBackend>::clone() const {
  // Clone is not easily implementable without storing the original config
  // For WRF, geometry instances should be created from configuration
  throw std::runtime_error(
      "WRFGeometry::clone() not implemented - create new instance from "
      "configuration instead");
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
      unstaggered_grid_.nz = dim_bottom_top.getSize();

      // Setup U-staggered grid (staggered in X direction)
      u_staggered_grid_.name = "u_staggered";
      u_staggered_grid_.is_staggered = true;
      u_staggered_grid_.size = dim_west_east_stag.getSize() *
                               dim_south_north.getSize() *
                               dim_bottom_top.getSize();
      u_staggered_grid_.resize_2d_coords(dim_west_east_stag.getSize(),
                                         dim_south_north.getSize());
      u_staggered_grid_.nz = dim_bottom_top.getSize();

      // Setup V-staggered grid (staggered in Y direction)
      v_staggered_grid_.name = "v_staggered";
      v_staggered_grid_.is_staggered = true;
      v_staggered_grid_.size = dim_west_east.getSize() *
                               dim_south_north_stag.getSize() *
                               dim_bottom_top.getSize();
      v_staggered_grid_.resize_2d_coords(dim_west_east.getSize(),
                                         dim_south_north_stag.getSize());
      v_staggered_grid_.nz = dim_bottom_top.getSize();

      // Setup W-staggered grid (staggered in Z direction)
      w_staggered_grid_.name = "w_staggered";
      w_staggered_grid_.is_staggered = true;
      w_staggered_grid_.size = dim_west_east.getSize() *
                               dim_south_north.getSize() *
                               dim_bottom_top_stag.getSize();
      w_staggered_grid_.resize_2d_coords(dim_west_east.getSize(),
                                         dim_south_north.getSize());
      w_staggered_grid_.nz = dim_bottom_top_stag.getSize();

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

      // Read vertical coordinate arrays (ZNU for unstaggered, ZNW for
      // staggered)
      auto var_znu = wrf_file.getVar("ZNU");
      auto var_znw = wrf_file.getVar("ZNW");

      if (!var_znu.isNull()) {
        // ZNU: eta values on half (mass) levels - unstaggered vertical
        // coordinates
        std::vector<size_t> znu_start = {0, 0};
        std::vector<size_t> znu_count = {1, dim_bottom_top.getSize()};
        std::vector<double> znu_data(dim_bottom_top.getSize());
        var_znu.getVar(znu_start, znu_count, znu_data.data());

        // Store ZNU coordinates in grids that use unstaggered vertical levels
        unstaggered_grid_.resize_vertical_coords(dim_bottom_top.getSize());
        unstaggered_grid_.vertical_coords = znu_data;

        u_staggered_grid_.resize_vertical_coords(dim_bottom_top.getSize());
        u_staggered_grid_.vertical_coords = znu_data;

        v_staggered_grid_.resize_vertical_coords(dim_bottom_top.getSize());
        v_staggered_grid_.vertical_coords = znu_data;

        // ZNU coordinates loaded successfully
      } else {
        std::cerr << "Warning: ZNU (unstaggered vertical coordinates) not "
                     "found in WRF file"
                  << std::endl;
      }

      if (!var_znw.isNull()) {
        // ZNW: eta values on full (w) levels - staggered vertical coordinates
        std::vector<size_t> znw_start = {0, 0};
        std::vector<size_t> znw_count = {1, dim_bottom_top_stag.getSize()};
        std::vector<double> znw_data(dim_bottom_top_stag.getSize());
        var_znw.getVar(znw_start, znw_count, znw_data.data());

        // Store ZNW coordinates in W-staggered grid
        w_staggered_grid_.resize_vertical_coords(dim_bottom_top_stag.getSize());
        w_staggered_grid_.vertical_coords = znw_data;

        // ZNW coordinates loaded successfully
      } else {
        std::cerr << "Warning: ZNW (staggered vertical coordinates) not found "
                     "in WRF file"
                  << std::endl;
      }

      // Note: 1D coordinate arrays are not populated as they're not used
      // Geographic coordinates are available in the 2D longitude_2d/latitude_2d
      // arrays

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