/**
 * @file WRFState.hpp
 * @brief WRF state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <chrono>
#include <filesystem>
#include <memory>
#include <netcdf>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
#if defined(_WIN32) || defined(__APPLE__)
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#else
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#endif
#include "WRFGridInfo.hpp"

namespace metada::backends::wrf {

/**
 * @brief WRF state backend implementation for meteorological data
 *
 * @details This class implements a state backend for the Weather Research and
 * Forecasting (WRF) model. It provides comprehensive management of
 * meteorological state variables stored in NetCDF files and implements all
 * operations required by the State adapter interface.
 *
 * Key features:
 * - Multi-variable state management with configurable active variable
 * - Support for both 2D and 3D meteorological fields
 * - Automatic grid association based on variable staggering patterns
 * - NetCDF file I/O operations for state persistence
 * - Vector arithmetic operations for data assimilation algorithms
 * - Container-like interface for element access and iteration
 * - Geographic coordinate to grid index transformations
 *
 * The class supports WRF's staggered grid system:
 * - Mass points (unstaggered grid)
 * - U-wind points (X-staggered grid)
 * - V-wind points (Y-staggered grid)
 * - W-wind points (Z-staggered grid)
 *
 * @tparam ConfigBackend Configuration backend type providing file paths and
 * settings
 * @tparam GeometryBackend Geometry backend type providing grid information
 *
 * @see VariableGridInfo
 * @see framework::State
 */
template <typename ConfigBackend, typename GeometryBackend>
class WRFState {
 public:
  ///@{ @name Type Definitions
  /**
   * @brief Standard container type aliases for STL compatibility
   *
   * These type definitions enable the WRFState to be used with STL algorithms
   * and provide a consistent interface following standard container
   * conventions.
   */
  using value_type = double;               ///< Type of elements stored
  using reference = double&;               ///< Reference to element type
  using const_reference = const double&;   ///< Const reference to element type
  using pointer = double*;                 ///< Pointer to element type
  using const_pointer = const double*;     ///< Const pointer to element type
  using size_type = std::size_t;           ///< Type for sizes and indices
  using difference_type = std::ptrdiff_t;  ///< Type for pointer differences
  ///@}

  ///@{ @name Construction and Destruction
  /**
   * @brief Default constructor is deleted
   *
   * @details WRFState requires configuration and geometry information
   * for proper initialization, so default construction is not allowed.
   */
  WRFState() = delete;

  /**
   * @brief Copy constructor is deleted
   *
   * @details Copy construction is disabled to prevent expensive copying
   * of large meteorological datasets. Use clone() for explicit copying.
   */
  WRFState(const WRFState&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   *
   * @details Copy assignment is disabled to prevent expensive copying
   * of large meteorological datasets. Use clone() for explicit copying.
   */
  WRFState& operator=(const WRFState&) = delete;

  /**
   * @brief Primary constructor for WRF state initialization
   *
   * @details Constructs a WRF state by loading meteorological data from
   * a NetCDF file specified in the configuration. The constructor automatically
   * determines variable grid associations and validates dimensions against
   * the provided geometry.
   *
   * @param[in] config Configuration backend containing:
   *                   - "file": Path to WRF NetCDF file
   *                   - "variables": List of variables to load (optional)
   * @param[in] geometry Geometry backend providing WRF grid information
   *
   * @throws std::runtime_error If file path is missing or file cannot be opened
   * @throws std::runtime_error If variable dimensions don't match geometry
   * @throws netCDF::exceptions::NcException If NetCDF file reading fails
   */
  WRFState(const ConfigBackend& config, const GeometryBackend& geometry);

  /**
   * @brief Move constructor for efficient resource transfer
   *
   * @details Transfers ownership of all internal data structures without
   * copying large arrays. The moved-from object is left in a valid but
   * unspecified state.
   *
   * @param[in] other WRF state to move from
   */
  WRFState(WRFState&& other) noexcept;

  /**
   * @brief Move assignment operator for efficient resource transfer
   *
   * @details Transfers ownership of all internal data structures without
   * copying large arrays. The moved-from object is reset to an uninitialized
   * state.
   *
   * @param[in] other WRF state to move from
   * @return Reference to this state after assignment
   */
  WRFState& operator=(WRFState&& other) noexcept;

  /**
   * @brief Destructor with automatic resource cleanup
   *
   * @details Automatically releases all allocated memory and closes
   * any open file handles. No explicit cleanup is required.
   */
  ~WRFState() = default;
  ///@}

  ///@{ @name State Cloning
  /**
   * @brief Create an independent copy of this state
   *
   * @details Creates a deep copy of all state data including variables,
   * dimensions, and configuration. The cloned state is fully independent
   * and can be modified without affecting the original.
   *
   * @return std::unique_ptr<WRFState> Unique pointer to the cloned state
   *
   * @note This operation involves copying potentially large arrays,
   *       so it may be expensive for states with many variables
   */
  std::unique_ptr<WRFState> clone() const;
  ///@}

  ///@{ @name Data Access
  /**
   * @brief Get mutable access to first variable data (concept compliance)
   *
   * @details Returns a void pointer to the raw data array of the first
   * loaded variable. This method is provided for StateBackendImpl concept
   * compliance. For explicit variable access, use getData(variableName).
   *
   * @return void* Pointer to first variable data, or nullptr if no variables
   * @throws std::runtime_error If no variables are loaded
   *
   * @warning Use with caution - no bounds checking is performed
   * @see getData(const std::string&) for explicit variable access
   */
  void* getData();

  /**
   * @brief Get const access to first variable data (concept compliance)
   *
   * @details Returns a const void pointer to the raw data array of the first
   * loaded variable. This method is provided for StateBackendImpl concept
   * compliance. For explicit variable access, use getData(variableName).
   *
   * @return const void* Const pointer to first variable data, or nullptr if no
   * variables
   * @throws std::runtime_error If no variables are loaded
   *
   * @see getData(const std::string&) for explicit variable access
   */
  const void* getData() const;

  /**
   * @brief Get mutable access to specific variable data
   *
   * @details Returns a void pointer to the raw data array of the specified
   * variable. The data is stored as contiguous double values.
   *
   * @param[in] variableName Name of the variable to access
   * @return void* Pointer to variable data, or nullptr if variable doesn't
   * exist
   *
   * @warning Use with caution - no bounds checking is performed
   */
  void* getData(const std::string& variableName);

  /**
   * @brief Get const access to specific variable data
   *
   * @details Returns a const void pointer to the raw data array of the
   * specified variable. The data is stored as contiguous double values.
   *
   * @param[in] variableName Name of the variable to access
   * @return const void* Const pointer to variable data, or nullptr if variable
   * doesn't exist
   */
  const void* getData(const std::string& variableName) const;

  /**
   * @brief Get mutable reference to variable array
   *
   * @details Returns a direct reference to the xtensor array for the specified
   * variable, allowing full array operations and element access.
   *
   * @param[in] variableName Name of the variable to access
   * @return xt::xarray<double>& Reference to the variable's data array
   * @throws std::out_of_range If variable doesn't exist
   */
  xt::xarray<double> getVariable(const std::string& variableName) {
    const auto& dims = dimensions_.at(variableName);
    size_t offset = variable_offsets_.at(variableName);
    size_t size = 1;
    for (size_t dim : dims) size *= dim;

    // Create a view into the flattened data
    return xt::adapt(flattened_data_.data() + offset, dims);
  }

  /**
   * @brief Get const reference to variable array
   *
   * @details Returns a direct const reference to the xtensor array for the
   * specified variable, allowing read-only array operations and element access.
   *
   * @param[in] variableName Name of the variable to access
   * @return const xt::xarray<double>& Const reference to the variable's data
   * array
   * @throws std::out_of_range If variable doesn't exist
   */
  xt::xarray<double> getVariable(const std::string& variableName) const {
    const auto& dims = dimensions_.at(variableName);
    size_t offset = variable_offsets_.at(variableName);
    size_t size = 1;
    for (size_t dim : dims) size *= dim;

    // Create a view into the flattened data
    return xt::adapt(flattened_data_.data() + offset, dims);
  }
  ///@}

  ///@{ @name Variable Management
  /**
   * @brief Get list of all available variable names
   *
   * @details Returns a vector containing the names of all variables
   * loaded from the WRF NetCDF file.
   *
   * @return const std::vector<std::string>& Reference to vector of variable
   * names
   */
  const std::vector<std::string>& getVariableNames() const;

  /**
   * @brief Check if a variable exists in the state
   *
   * @details Checks whether a variable with the given name has been loaded
   * into this state instance.
   *
   * @param[in] variableName Name of the variable to check
   * @return bool True if variable exists, false otherwise
   */
  bool hasVariable(const std::string& variableName) const {
    return variable_offsets_.find(variableName) != variable_offsets_.end();
  }

  /**
   * @brief Get size of a specific variable
   *
   * @details Returns the total number of elements in the specified variable.
   *
   * @param[in] variableName Name of the variable
   * @return size_t Number of elements in the variable
   * @throws std::out_of_range If variable doesn't exist
   */
  size_t getVariableSize(const std::string& variableName) const {
    const auto& dims = dimensions_.at(variableName);
    size_t size = 1;
    for (size_t dim : dims) size *= dim;
    return size;
  }

  /**
   * @brief Get dimensions of a specific variable
   *
   * @details Returns the dimensions of the specified variable as a vector.
   * For 3D variables: [Z, Y, X] dimensions
   * For 2D variables: [Y, X] dimensions
   *
   * @param[in] variableName Name of the variable
   * @return const std::vector<size_t>& Reference to the variable's dimensions
   * @throws std::out_of_range If variable doesn't exist
   */
  const std::vector<size_t>& getVariableDimensions(
      const std::string& variableName) const {
    return dimensions_.at(variableName);
  }
  ///@}

  ///@{ @name Size
  /**
   * @brief Get the total number of elements across all variables
   *
   * @details Calculates the sum of all elements in all loaded variables.
   * This represents the total size of the complete state vector.
   *
   * @return size_t Total number of elements across all variables
   */
  size_t size() const;
  ///@}

  ///@{ @name Observation Operator Support
  /**
   * @brief Extract data for observation operators
   *
   * @details This method extracts the essential data that observation operators
   * need from the state, optimized for efficient access. It returns a structure
   * containing raw pointers to the data arrays and metadata.
   *
   * @return Structure containing observation operator data
   */
  struct ObsOperatorData {
    // Grid dimensions
    int nx, ny, nz;

    // State variables (as raw pointers for efficiency)
    const double* u = nullptr;
    const double* v = nullptr;
    const double* t = nullptr;
    const double* q = nullptr;
    const double* psfc = nullptr;
    const double* ph = nullptr;   // Geopotential perturbation
    const double* phb = nullptr;  // Base state geopotential
    const double* hf = nullptr;   // Height field (calculated from PH and PHB)
    const double* hgt = nullptr;  // Terrain height (from HGT variable)

    // Grid metadata
    const double* lats_2d = nullptr;
    const double* lons_2d = nullptr;
    const double* levels = nullptr;

    // Variable dimensions (for staggered grids)
    struct VariableDims {
      int nx, ny, nz;
      bool is_staggered = false;
    };

    VariableDims u_dims, v_dims, t_dims, q_dims, psfc_dims, ph_dims, phb_dims,
        hf_dims, hgt_dims;
  };

  ObsOperatorData getObsOperatorData() const;
  ///@}

  ///@{ @name Total State Vector Access
  /**
   * @brief Access element in total state vector (unchecked)
   *
   * @details Provides access to elements across all variables as if they were
   * concatenated into a single vector. Variables are accessed in the order
   * they appear in variableNames_. No bounds checking is performed.
   *
   * @param[in] idx Linear index into the total state vector
   * @return reference Reference to the element at the specified index
   *
   * @warning No bounds checking - undefined behavior if idx >= size()
   * @see at(size_type) for bounds-checked access
   */
  reference operator[](size_type idx) { return flattened_data_[idx]; }

  /**
   * @brief Access element in total state vector (unchecked, const)
   *
   * @details Provides const access to elements across all variables as if they
   * were concatenated into a single vector. Variables are accessed in the order
   * they appear in variableNames_. No bounds checking is performed.
   *
   * @param[in] idx Linear index into the total state vector
   * @return const_reference Const reference to the element at the specified
   * index
   *
   * @warning No bounds checking - undefined behavior if idx >= size()
   * @see at(size_type) for bounds-checked access
   */
  const_reference operator[](size_type idx) const {
    return flattened_data_[idx];
  }

  /**
   * @brief Access element in total state vector (bounds-checked)
   *
   * @details Provides access to elements across all variables as if they were
   * concatenated into a single vector with bounds checking.
   *
   * @param[in] idx Linear index into the total state vector
   * @return reference Reference to the element at the specified index
   * @throws std::out_of_range If idx >= size()
   */
  reference at(size_type idx) {
    if (idx >= size()) {
      throw std::out_of_range("Index out of range");
    }
    return operator[](idx);
  }

  /**
   * @brief Access element in total state vector (bounds-checked, const)
   *
   * @details Provides const access to elements across all variables as if they
   * were concatenated into a single vector with bounds checking.
   *
   * @param[in] idx Linear index into the total state vector
   * @return const_reference Const reference to the element at the specified
   * index
   * @throws std::out_of_range If idx >= size()
   */
  const_reference at(size_type idx) const {
    if (idx >= size()) {
      throw std::out_of_range("Index out of range");
    }
    return operator[](idx);
  }

  /**
   * @brief Get reference to first element of total state vector
   *
   * @return reference Reference to the first element
   * @warning Undefined behavior if state is empty
   */
  reference front() { return flattened_data_.front(); }

  /**
   * @brief Get const reference to first element of total state vector
   *
   * @return const_reference Const reference to the first element
   * @warning Undefined behavior if state is empty
   */
  const_reference front() const { return flattened_data_.front(); }

  /**
   * @brief Get reference to last element of total state vector
   *
   * @return reference Reference to the last element
   * @warning Undefined behavior if state is empty
   */
  reference back() { return flattened_data_.back(); }

  /**
   * @brief Get const reference to last element of total state vector
   *
   * @return const_reference Const reference to the last element
   * @warning Undefined behavior if state is empty
   */
  const_reference back() const { return flattened_data_.back(); }

  /**
   * @brief Access element at geographic location using first variable's grid
   *
   * @details Converts geographic coordinates to grid indices using the first
   * variable's grid type and returns a reference to the corresponding element.
   * Uses nearest-neighbor interpolation to find the closest grid point.
   *
   * @param[in] loc Location containing geographic coordinates
   * @return reference Reference to element at the nearest grid point in first
   * variable
   * @throws std::runtime_error If geographic coordinates are unavailable or no
   * variables loaded
   * @throws std::out_of_range If first variable doesn't exist
   *
   * @note Uses first variable's grid for coordinate conversion
   * @note For explicit variable selection, use at(variableName, location)
   */
  double& at(const framework::Location& loc) {
    if (variableNames_.empty()) {
      throw std::runtime_error("No variables loaded for geographic access");
    }
    return at(variableNames_[0], loc);
  }

  /**
   * @brief Access element at geographic location using first variable's grid
   * (const)
   *
   * @details Converts geographic coordinates to grid indices using the first
   * variable's grid type and returns a const reference to the corresponding
   * element. Uses nearest-neighbor interpolation to find the closest grid
   * point.
   *
   * @param[in] loc Location containing geographic coordinates
   * @return const_reference Const reference to element at the nearest grid
   * point in first variable
   * @throws std::runtime_error If geographic coordinates are unavailable or no
   * variables loaded
   * @throws std::out_of_range If first variable doesn't exist
   *
   * @note Uses first variable's grid for coordinate conversion
   * @note For explicit variable selection, use at(variableName, location)
   */
  const double& at(const framework::Location& loc) const {
    if (variableNames_.empty()) {
      throw std::runtime_error("No variables loaded for geographic access");
    }
    return at(variableNames_[0], loc);
  }
  ///@}

  ///@{ @name Element Access by Variable Name
  /**
   * @brief Access element by variable name and linear index (bounds-checked)
   *
   * @details Provides access to elements in the specified variable with
   * bounds checking. Throws an exception if the index is out of range.
   *
   * @param[in] variableName Name of the variable to access
   * @param[in] idx Linear index into the variable's data
   * @return reference Reference to the element at the specified index
   * @throws std::out_of_range If variable doesn't exist or idx >= variable size
   */
  reference at(const std::string& variableName, size_type idx) {
    size_t offset = variable_offsets_.at(variableName);
    const auto& dims = dimensions_.at(variableName);
    size_t var_size = 1;
    for (size_t dim : dims) var_size *= dim;

    if (idx >= var_size) {
      throw std::out_of_range("Index out of range for variable: " +
                              variableName);
    }

    return flattened_data_[offset + idx];
  }

  /**
   * @brief Access element by variable name and linear index (bounds-checked,
   * const)
   *
   * @details Provides const access to elements in the specified variable with
   * bounds checking. Throws an exception if the index is out of range.
   *
   * @param[in] variableName Name of the variable to access
   * @param[in] idx Linear index into the variable's data
   * @return const_reference Const reference to the element at the specified
   * index
   * @throws std::out_of_range If variable doesn't exist or idx >= variable size
   */
  const_reference at(const std::string& variableName, size_type idx) const {
    size_t offset = variable_offsets_.at(variableName);
    const auto& dims = dimensions_.at(variableName);
    size_t var_size = 1;
    for (size_t dim : dims) var_size *= dim;

    if (idx >= var_size) {
      throw std::out_of_range("Index out of range for variable: " +
                              variableName);
    }

    return flattened_data_[offset + idx];
  }
  ///@}

  ///@{ @name State Initialization
  /**
   * @brief Set all variable values to zero
   *
   * @details Resets all elements in all loaded variables to zero.
   * This is commonly used in data assimilation algorithms for
   * initializing perturbations or increments.
   */
  void zero();
  ///@}

  ///@{ @name Vector Arithmetic Operations
  /**
   * @brief Calculate dot product with another state
   *
   * @details Computes the scalar dot product between this state and another
   * state by summing element-wise products across all variables. Both states
   * must have compatible variable sets and dimensions.
   *
   * @param[in] other State to calculate dot product with
   * @return double Scalar dot product result
   * @throws std::runtime_error If states have incompatible variables or
   * dimensions
   *
   * @note Used in data assimilation algorithms for covariance calculations
   */
  double dot(const WRFState& other) const;

  /**
   * @brief Calculate L2 norm of this state
   *
   * @details Computes the L2 (Euclidean) norm by taking the square root
   * of the sum of squares of all elements across all variables.
   *
   * @return double L2 norm value (always non-negative)
   *
   * @note Used for convergence testing and error quantification
   */
  double norm() const;

  /**
   * @brief Add another state to this one (element-wise)
   *
   * @details Performs element-wise addition: this += other. Both states
   * must have compatible variable sets and dimensions.
   *
   * @param[in] other State to add to this state
   * @throws std::runtime_error If states have incompatible variables or
   * dimensions
   */
  void add(const WRFState& other);

  /**
   * @brief Subtract another state from this one (element-wise)
   *
   * @details Performs element-wise subtraction: this -= other. Both states
   * must have compatible variable sets and dimensions.
   *
   * @param[in] other State to subtract from this state
   * @throws std::runtime_error If states have incompatible variables or
   * dimensions
   */
  void subtract(const WRFState& other);

  /**
   * @brief Multiply all elements by a scalar value
   *
   * @details Performs element-wise scalar multiplication: this *= scalar.
   * All elements in all variables are multiplied by the given value.
   *
   * @param[in] scalar Value to multiply all elements by
   */
  void multiply(double scalar);
  ///@}

  ///@{ @name State Comparison
  /**
   * @brief Check if this state equals another within tolerance
   *
   * @details Compares all elements of all variables between two states
   * using a small numerical tolerance (1e-10). States are considered
   * equal if maximum absolute difference is below tolerance.
   *
   * @param[in] other State to compare with
   * @return bool True if states are numerically equal, false otherwise
   */
  bool equals(const WRFState& other) const;
  ///@}

  ///@{ @name State Status
  /**
   * @brief Check if state has been properly initialized
   *
   * @details Returns true if the state has successfully loaded data from
   * a NetCDF file and is ready for use.
   *
   * @return bool True if initialized and ready for use, false otherwise
   */
  bool isInitialized() const;
  ///@}

  ///@{ @name Geographic Coordinate Access
  /**
   * @brief Access element at geographic location (mutable)
   *
   * @details Converts geographic coordinates (latitude, longitude, level) to
   * grid indices and returns a reference to the corresponding element in the
   * specified variable. Uses nearest-neighbor interpolation to find the closest
   * grid point.
   *
   * @param[in] variableName Name of the variable to access
   * @param[in] loc Location containing geographic coordinates
   * @return reference Reference to element at the nearest grid point
   * @throws std::runtime_error If geographic coordinates are unavailable
   * @throws std::out_of_range If variable doesn't exist
   *
   * @note For production use, consider implementing bilinear interpolation
   */
  double& at(const std::string& variableName, const framework::Location& loc) {
    // Convert geographic coordinates to grid indices
    auto [lat, lon, level] = loc.getGeographicCoords();

    // Get the appropriate geometry grid for the specified variable
    const auto& geometry_grid = getVariableGeometryGrid(variableName);

    // Find the closest grid point to the geographic coordinates
    // This is a simplified approach - in practice, you might want bilinear
    // interpolation
    size_t i, j, k;

    // Find i, j indices from 2D coordinates
    if (geometry_grid.has_2d_coords()) {
      // Find closest grid point by searching through the coordinate arrays
      double min_dist = std::numeric_limits<double>::max();
      size_t min_idx = 0;

      for (size_t idx = 0; idx < geometry_grid.longitude_2d.size(); ++idx) {
        double grid_lon = geometry_grid.longitude_2d[idx];
        double grid_lat = geometry_grid.latitude_2d[idx];

        double dist = std::sqrt((lon - grid_lon) * (lon - grid_lon) +
                                (lat - grid_lat) * (lat - grid_lat));

        if (dist < min_dist) {
          min_dist = dist;
          min_idx = idx;
        }
      }

      // Convert linear index to i, j
      j = min_idx / geometry_grid.nx;
      i = min_idx % geometry_grid.nx;
    } else {
      throw std::runtime_error(
          "Geographic coordinates not available for WRF state access");
    }

    // Find k index from vertical level
    if (geometry_grid.has_vertical_coords()) {
      // Find closest vertical level
      double min_dist = std::numeric_limits<double>::max();
      k = 0;

      for (size_t z = 0; z < geometry_grid.vertical_coords.size(); ++z) {
        double dist = std::abs(level - geometry_grid.vertical_coords[z]);
        if (dist < min_dist) {
          min_dist = dist;
          k = z;
        }
      }
    } else {
      k = 0;  // Default to first level if no vertical coordinates
    }

    // Access the array with proper indexing
    size_t offset = variable_offsets_.at(variableName);
    const auto& dims = dimensions_.at(variableName);

    // Calculate linear index from 3D coordinates
    size_t linear_idx;
    if (dims.size() == 3) {
      linear_idx = k * dims[1] * dims[2] + j * dims[2] + i;  // 3D: [Z, Y, X]
    } else {
      linear_idx = j * dims[1] + i;  // 2D: [Y, X]
    }

    return flattened_data_[offset + linear_idx];
  }

  /**
   * @brief Access element at geographic location (const)
   *
   * @details Converts geographic coordinates (latitude, longitude, level) to
   * grid indices and returns a const reference to the corresponding element in
   * the specified variable. Uses nearest-neighbor interpolation to find the
   * closest grid point.
   *
   * @param[in] variableName Name of the variable to access
   * @param[in] loc Location containing geographic coordinates
   * @return const_reference Const reference to element at the nearest grid
   * point
   * @throws std::runtime_error If geographic coordinates are unavailable
   * @throws std::out_of_range If variable doesn't exist
   *
   * @note For production use, consider implementing bilinear interpolation
   */
  const double& at(const std::string& variableName,
                   const framework::Location& loc) const {
    // Convert geographic coordinates to grid indices
    auto [lat, lon, level] = loc.getGeographicCoords();

    // Get the appropriate geometry grid for the specified variable
    const auto& geometry_grid = getVariableGeometryGrid(variableName);

    // Find the closest grid point to the geographic coordinates
    size_t i, j, k;

    // Find i, j indices from 2D coordinates
    if (geometry_grid.has_2d_coords()) {
      // Find closest grid point by searching through the coordinate arrays
      double min_dist = std::numeric_limits<double>::max();
      size_t min_idx = 0;

      for (size_t idx = 0; idx < geometry_grid.longitude_2d.size(); ++idx) {
        double grid_lon = geometry_grid.longitude_2d[idx];
        double grid_lat = geometry_grid.latitude_2d[idx];

        double dist = std::sqrt((lon - grid_lon) * (lon - grid_lon) +
                                (lat - grid_lat) * (lat - grid_lat));

        if (dist < min_dist) {
          min_dist = dist;
          min_idx = idx;
        }
      }

      // Convert linear index to i, j
      j = min_idx / geometry_grid.nx;
      i = min_idx % geometry_grid.nx;
    } else {
      throw std::runtime_error(
          "Geographic coordinates not available for WRF state access");
    }

    // Find k index from vertical level
    if (geometry_grid.has_vertical_coords()) {
      // Find closest vertical level
      double min_dist = std::numeric_limits<double>::max();
      k = 0;

      for (size_t z = 0; z < geometry_grid.vertical_coords.size(); ++z) {
        double dist = std::abs(level - geometry_grid.vertical_coords[z]);
        if (dist < min_dist) {
          min_dist = dist;
          k = z;
        }
      }
    } else {
      k = 0;  // Default to first level if no vertical coordinates
    }

    // Access the array with proper indexing
    size_t offset = variable_offsets_.at(variableName);
    const auto& dims = dimensions_.at(variableName);

    // Calculate linear index from 3D coordinates
    size_t linear_idx;
    if (dims.size() == 3) {
      linear_idx = k * dims[1] * dims[2] + j * dims[2] + i;  // 3D: [Z, Y, X]
    } else {
      linear_idx = j * dims[1] + i;  // 2D: [Y, X]
    }

    return flattened_data_[offset + linear_idx];
  }
  ///@}

  ///@{ @name File I/O Operations
  /**
   * @brief Save state data to NetCDF file by copying original and updating
   * variables
   *
   * @details Copies the original WRF NetCDF file and updates only the variables
   * that have been loaded and potentially modified. This preserves all original
   * file structure, dimensions, attributes, and variables not explicitly
   * loaded. The method is more efficient and maintains WRF file compatibility.
   *
   * @param[in] filename Path where the NetCDF file should be created
   * @throws std::runtime_error If file creation or copying fails
   * @throws netCDF::exceptions::NcException If NetCDF operations fail
   *
   * @note The output file will contain:
   *       - All original variables and attributes from source file
   *       - Updated data for loaded variables
   *       - Preserved file structure and metadata
   */
  void saveToFile(const std::string& filename) const {
    try {
      // First, copy the original file to preserve all structure
      std::filesystem::copy_file(
          wrfFilename_, filename,
          std::filesystem::copy_options::overwrite_existing);

      // Open the copied file for modification
      netCDF::NcFile nc_file(filename, netCDF::NcFile::write);

      if (nc_file.isNull()) {
        throw std::runtime_error("Failed to open copied NetCDF file: " +
                                 filename);
      }

      // Update only the variables that were loaded (and potentially modified)
      for (const auto& varName : variableNames_) {
        auto nc_var = nc_file.getVar(varName);

        if (!nc_var.isNull()) {
          size_t offset = variable_offsets_.at(varName);

          // Write updated data from flattened storage
          nc_var.putVar(flattened_data_.data() + offset);
        } else {
          // Variable doesn't exist in original file, skip it
          // This could happen if we loaded variables not in the original file
          std::cerr << "Warning: Variable " << varName
                    << " not found in original file, skipping update"
                    << std::endl;
        }
      }

      // Add metadata about the update
      nc_file.putAtt("metada_updated", "true");

      // Convert timestamp to string for NetCDF compatibility
      auto now = std::chrono::system_clock::now();
      auto timestamp = std::chrono::duration_cast<std::chrono::seconds>(
                           now.time_since_epoch())
                           .count();
      nc_file.putAtt("metada_update_time", std::to_string(timestamp));

      // Add list of updated variables as a single string
      std::string updated_vars;
      for (size_t i = 0; i < variableNames_.size(); ++i) {
        if (i > 0) updated_vars += ", ";
        updated_vars += variableNames_[i];
      }
      nc_file.putAtt("metada_updated_variables", updated_vars);

      nc_file.close();

    } catch (const std::filesystem::filesystem_error& e) {
      throw std::runtime_error("File system error while copying WRF file: " +
                               std::string(e.what()));
    } catch (const netCDF::exceptions::NcException& e) {
      throw std::runtime_error("NetCDF error while saving WRF state: " +
                               std::string(e.what()));
    } catch (const std::exception& e) {
      throw std::runtime_error("Error saving WRF state: " +
                               std::string(e.what()));
    }
  }
  ///@}

  ///@{ @name Geometry Access
  /**
   * @brief Get reference to the associated geometry backend
   *
   * @details Returns a const reference to the geometry backend that provides
   * grid information and coordinate transformations for this state.
   *
   * @return const GeometryBackend& Reference to the geometry backend
   */
  const GeometryBackend& geometry() const { return geometry_; }
  ///@}

  ///@{ @name Variable-to-Grid Mapping
  /**
   * @brief Get grid information for a specific variable
   *
   * @details Returns the VariableGridInfo structure containing dimension names,
   * staggering configuration, and grid type association for the specified
   * variable.
   *
   * @param[in] variableName Name of the variable to query
   * @return const VariableGridInfo& Reference to grid information structure
   * @throws std::out_of_range If the specified variable doesn't exist
   */
  const VariableGridInfo& getVariableGridInfo(
      const std::string& variableName) const {
    auto it = variable_grid_info_.find(variableName);
    if (it == variable_grid_info_.end()) {
      throw std::out_of_range("Variable grid info not found: " + variableName);
    }
    return it->second;
  }

  /**
   * @brief Get the appropriate geometry grid for a variable
   *
   * @details Returns a reference to the geometry grid (unstaggered,
   * U-staggered, V-staggered, or W-staggered) that corresponds to the specified
   * variable's grid type. This enables proper coordinate transformations and
   * interpolations.
   *
   * @param[in] variableName Name of the variable
   * @return const auto& Reference to the appropriate GridDimensionInfo from
   * WRFGeometry
   * @throws std::out_of_range If variable doesn't exist
   * @throws std::runtime_error If grid type is unknown
   */
  const auto& getVariableGeometryGrid(const std::string& variableName) const {
    const auto& grid_info = getVariableGridInfo(variableName);

    switch (grid_info.grid_type) {
      case VariableGridInfo::GridType::UNSTAGGERED:
        return geometry_.unstaggered_info();
      case VariableGridInfo::GridType::U_STAGGERED:
        return geometry_.u_staggered_info();
      case VariableGridInfo::GridType::V_STAGGERED:
        return geometry_.v_staggered_info();
      case VariableGridInfo::GridType::W_STAGGERED:
        return geometry_.w_staggered_info();
      default:
        throw std::runtime_error("Unknown grid type for variable: " +
                                 variableName);
    }
  }

  /**
   * @brief Get all variables associated with a specific grid type
   *
   * @details Searches through all loaded variables and returns a list of those
   * that use the specified grid type (unstaggered, U-staggered, V-staggered, or
   * W-staggered). Useful for performing grid-type-specific operations.
   *
   * @param[in] gridType The grid type to search for
   * @return std::vector<std::string> Vector of variable names using the
   * specified grid type
   */
  std::vector<std::string> getVariablesByGridType(
      VariableGridInfo::GridType gridType) const {
    std::vector<std::string> result;
    for (const auto& [varName, gridInfo] : variable_grid_info_) {
      if (gridInfo.grid_type == gridType) {
        result.push_back(varName);
      }
    }
    return result;
  }

  /**
   * @brief Validate that variable dimensions match geometry grid dimensions
   *
   * @details Checks if the dimensions of a specific variable are consistent
   * with the expected dimensions from the associated geometry grid. This
   * validation ensures data integrity and proper grid associations.
   *
   * @param[in] variableName Name of the variable to validate
   * @return bool True if dimensions match expected values, false otherwise
   */
  bool validateVariableGeometry(const std::string& variableName) const {
    try {
      const auto& grid_info = getVariableGridInfo(variableName);
      const auto& var_dims = dimensions_.at(variableName);

      // Check if 2D/3D dimensions match
      if (var_dims.size() >= 2) {
        size_t expected_nx, expected_ny, expected_nz;
        switch (grid_info.grid_type) {
          case VariableGridInfo::GridType::UNSTAGGERED:
            expected_nx = geometry_.x_dim();
            expected_ny = geometry_.y_dim();
            expected_nz = geometry_.z_dim();
            break;
          case VariableGridInfo::GridType::U_STAGGERED:
            expected_nx = geometry_.x_stag_dim();
            expected_ny = geometry_.y_dim();
            expected_nz = geometry_.z_dim();
            break;
          case VariableGridInfo::GridType::V_STAGGERED:
            expected_nx = geometry_.x_dim();
            expected_ny = geometry_.y_stag_dim();
            expected_nz = geometry_.z_dim();
            break;
          case VariableGridInfo::GridType::W_STAGGERED:
            expected_nx = geometry_.x_dim();
            expected_ny = geometry_.y_dim();
            expected_nz = geometry_.z_stag_dim();
            break;
          default:
            return false;
        }

        // For 2D variables: [Y, X]
        // For 3D variables: [Z, Y, X]
        if (var_dims.size() == 2) {
          size_t y_idx = 0;
          size_t x_idx = 1;
          return (var_dims[x_idx] == expected_nx &&
                  var_dims[y_idx] == expected_ny);
        } else if (var_dims.size() == 3) {
          size_t z_idx = 0;
          size_t y_idx = 1;
          size_t x_idx = 2;
          return (var_dims[x_idx] == expected_nx &&
                  var_dims[y_idx] == expected_ny &&
                  var_dims[z_idx] == expected_nz);
        }
      }
      return false;
    } catch (const std::exception&) {
      return false;
    }
  }
  ///@}

  ///@{ @name Stream Output
  /**
   * @brief Stream insertion operator for debugging and logging
   *
   * @details Outputs a human-readable representation of the WRF state including
   * initialization status, source filename, loaded variables, and total data
   * size. Useful for debugging and logging purposes.
   *
   * @param[in,out] os Output stream to write to
   * @param[in] state WRFState instance to output
   * @return std::ostream& Reference to the output stream for chaining
   *
   * @note Output format: WRFState{initialized: <bool>, filename: "<path>",
   *                      variables: [...], total_size: <num>}
   */
  friend std::ostream& operator<<(std::ostream& os, const WRFState& state) {
    os << "WRFState{";
    os << "initialized: " << (state.initialized_ ? "true" : "false");
    os << ", filename: \"" << state.wrfFilename_ << "\"";
    os << ", variables: [";

    for (size_t i = 0; i < state.variableNames_.size(); ++i) {
      if (i > 0) os << ", ";
      os << "\"" << state.variableNames_[i] << "\"";
    }
    os << "]";

    os << ", total_size: " << state.size();

    // Add statistics if data is available
    if (!state.flattened_data_.empty()) {
      auto min_it = std::min_element(state.flattened_data_.begin(),
                                     state.flattened_data_.end());
      auto max_it = std::max_element(state.flattened_data_.begin(),
                                     state.flattened_data_.end());
      double sum = std::accumulate(state.flattened_data_.begin(),
                                   state.flattened_data_.end(), 0.0);
      double mean = sum / state.flattened_data_.size();

      os << ", min: " << *min_it;
      os << ", max: " << *max_it;
      os << ", mean: " << mean;
    }

    os << "}";
    return os;
  }
  ///@}

 private:
  ///@{ @name Private Construction and Cloning
  /**
   * @brief Private constructor for cloning operations
   *
   * @details Special constructor used internally for cloning that skips
   * the expensive NetCDF file loading process. Data is copied directly
   * from the source state instead.
   *
   * @param[in] config Configuration backend reference
   * @param[in] geometry Geometry backend reference
   * @param[in] skip_loading Flag to skip file loading (must be true)
   */
  WRFState(const ConfigBackend& config, const GeometryBackend& geometry,
           bool skip_loading);

  /**
   * @brief Factory method for creating cloneable state instances
   *
   * @details Creates a new WRFState instance ready for cloning operations.
   * The instance is created with file loading disabled to avoid unnecessary
   * I/O operations during the cloning process.
   *
   * @param[in] config Configuration backend reference
   * @param[in] geometry Geometry backend reference
   * @return std::unique_ptr<WRFState> Unique pointer to the new state instance
   */
  static std::unique_ptr<WRFState> createForCloning(
      const ConfigBackend& config, const GeometryBackend& geometry) {
    return std::unique_ptr<WRFState>(new WRFState(config, geometry, true));
  }
  ///@}

  ///@{ @name Private Utility Methods
  /**
   * @brief Load meteorological data from WRF NetCDF file
   *
   * @details Reads variable data from the specified NetCDF file, determines
   * grid associations based on dimension names, and validates dimensions
   * against the geometry backend. Automatically detects variable staggering
   * patterns and associates them with appropriate grid types.
   *
   * @param[in] filename Path to the WRF NetCDF file
   * @param[in] variables List of variable names to load
   * @throws std::runtime_error If file cannot be opened or read
   * @throws netCDF::exceptions::NcException If NetCDF operations fail
   */
  void loadStateData(const std::string& filename,
                     const std::vector<std::string>& variables);

  /**
   * @brief Check compatibility between two WRF states
   *
   * @details Verifies that two states have the same variable names and
   * dimensions, enabling safe arithmetic operations between them.
   * Compatibility is required for add, subtract, dot product, and
   * comparison operations.
   *
   * @param[in] other State to check compatibility with
   * @return bool True if states are compatible, false otherwise
   */
  bool isCompatible(const WRFState& other) const;

  ///@}

  ///@{ @name Member Variables
  const ConfigBackend& config_;      ///< Reference to configuration backend
  const GeometryBackend& geometry_;  ///< Reference to geometry backend

  /// @name File Information
  ///@{
  std::string wrfFilename_;   ///< Path to source WRF NetCDF file
  bool initialized_ = false;  ///< Initialization status flag
  ///@}

  /// @name State Data Storage
  ///@{
  std::vector<double> flattened_data_;  ///< Single flattened state vector
  std::unordered_map<std::string, std::vector<size_t>>
      dimensions_;                          ///< Variable dimensions
  std::vector<std::string> variableNames_;  ///< List of loaded variables
  std::unordered_map<std::string, VariableGridInfo>
      variable_grid_info_;  ///< Variable-to-grid mappings
  std::unordered_map<std::string, size_t>
      variable_offsets_;  ///< Starting offsets for each variable in flattened
                          ///< data
  ///@}
  ///@}
};

// Constructor implementation with ConfigBackend and GeometryBackend
template <typename ConfigBackend, typename GeometryBackend>
WRFState<ConfigBackend, GeometryBackend>::WRFState(
    const ConfigBackend& config, const GeometryBackend& geometry)
    : config_(config),
      geometry_(geometry),
      wrfFilename_(config.Get("file").asString()),
      initialized_(false) {
  if (wrfFilename_.empty()) {
    throw std::runtime_error(
        "WRF input file path not specified in configuration");
  }

  // Get variables to load from config
  std::vector<std::string> variables;

  // Try to get variables list from config, otherwise use defaults
  try {
    variables = config.Get("variables").asVectorString();
  } catch (const std::exception&) {
    std::cerr << "Warning: No variables specified in configuration"
              << std::endl;
  }

  // Load state data from WRF NetCDF file
  loadStateData(wrfFilename_, variables);

  initialized_ = true;
}

// Private constructor for cloning (skips file loading)
template <typename ConfigBackend, typename GeometryBackend>
WRFState<ConfigBackend, GeometryBackend>::WRFState(
    const ConfigBackend& config, const GeometryBackend& geometry,
    [[maybe_unused]] bool skip_loading)
    : config_(config),
      geometry_(geometry),
      wrfFilename_(""),  // Will be set by clone
      initialized_(false) {
  // Skip file loading - data will be copied by clone method
}

// Move constructor implementation
template <typename ConfigBackend, typename GeometryBackend>
WRFState<ConfigBackend, GeometryBackend>::WRFState(
    WRFState<ConfigBackend, GeometryBackend>&& other) noexcept
    : config_(other.config_),
      geometry_(other.geometry_),
      wrfFilename_(std::move(other.wrfFilename_)),
      initialized_(other.initialized_),
      flattened_data_(std::move(other.flattened_data_)),
      dimensions_(std::move(other.dimensions_)),
      variableNames_(std::move(other.variableNames_)),
      variable_grid_info_(std::move(other.variable_grid_info_)),
      variable_offsets_(std::move(other.variable_offsets_)) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.variableNames_.clear();
}

// Move assignment operator implementation
template <typename ConfigBackend, typename GeometryBackend>
WRFState<ConfigBackend, GeometryBackend>&
WRFState<ConfigBackend, GeometryBackend>::operator=(
    WRFState<ConfigBackend, GeometryBackend>&& other) noexcept {
  if (this != &other) {
    // config_ and geometry_ are references, so no assignment needed
    wrfFilename_ = std::move(other.wrfFilename_);
    initialized_ = other.initialized_;
    flattened_data_ = std::move(other.flattened_data_);
    dimensions_ = std::move(other.dimensions_);
    variableNames_ = std::move(other.variableNames_);
    variable_grid_info_ = std::move(other.variable_grid_info_);
    variable_offsets_ = std::move(other.variable_offsets_);

    // Reset the moved-from object
    other.initialized_ = false;
    other.variableNames_.clear();
  }
  return *this;
}

// Clone implementation
template <typename ConfigBackend, typename GeometryBackend>
std::unique_ptr<WRFState<ConfigBackend, GeometryBackend>>
WRFState<ConfigBackend, GeometryBackend>::clone() const {
  // Create a new WRFState with the same config and geometry, skipping file
  // loading
  auto cloned = createForCloning(config_, geometry_);

  // Copy all the state data directly without reloading from file
  cloned->wrfFilename_ = this->wrfFilename_;
  cloned->initialized_ = this->initialized_;

  // Deep copy the variables to avoid sharing memory
  cloned->flattened_data_ = this->flattened_data_;
  cloned->dimensions_ = this->dimensions_;
  cloned->variableNames_ = this->variableNames_;
  cloned->variable_grid_info_ = this->variable_grid_info_;
  cloned->variable_offsets_ = this->variable_offsets_;

  return cloned;
}

// Data access implementation (concept compliance)
template <typename ConfigBackend, typename GeometryBackend>
void* WRFState<ConfigBackend, GeometryBackend>::getData() {
  if (!initialized_ || variableNames_.empty()) {
    return nullptr;
  }
  return flattened_data_.data();
}

// Const data access implementation (concept compliance)
template <typename ConfigBackend, typename GeometryBackend>
const void* WRFState<ConfigBackend, GeometryBackend>::getData() const {
  if (!initialized_ || variableNames_.empty()) {
    return nullptr;
  }
  return flattened_data_.data();
}

// Data access implementation (explicit variable)
template <typename ConfigBackend, typename GeometryBackend>
void* WRFState<ConfigBackend, GeometryBackend>::getData(
    const std::string& variableName) {
  if (!initialized_) {
    return nullptr;
  }

  try {
    size_t offset = variable_offsets_.at(variableName);
    return flattened_data_.data() + offset;
  } catch (const std::out_of_range&) {
    return nullptr;
  }
}

// Const data access implementation (explicit variable)
template <typename ConfigBackend, typename GeometryBackend>
const void* WRFState<ConfigBackend, GeometryBackend>::getData(
    const std::string& variableName) const {
  if (!initialized_) {
    return nullptr;
  }

  try {
    size_t offset = variable_offsets_.at(variableName);
    return flattened_data_.data() + offset;
  } catch (const std::out_of_range&) {
    return nullptr;
  }
}

// Get variable names implementation
template <typename ConfigBackend, typename GeometryBackend>
const std::vector<std::string>&
WRFState<ConfigBackend, GeometryBackend>::getVariableNames() const {
  return variableNames_;
}

// Size implementation
template <typename ConfigBackend, typename GeometryBackend>
size_t WRFState<ConfigBackend, GeometryBackend>::size() const {
  return flattened_data_.size();
}

// Zero implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::zero() {
  std::fill(flattened_data_.begin(), flattened_data_.end(), 0.0);
}

// Dot product implementation
template <typename ConfigBackend, typename GeometryBackend>
double WRFState<ConfigBackend, GeometryBackend>::dot(
    const WRFState<ConfigBackend, GeometryBackend>& other) const {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for dot product");
  }

  double result = 0.0;

  // Compute dot product directly on flattened data
  for (size_t i = 0; i < flattened_data_.size(); ++i) {
    result += flattened_data_[i] * other.flattened_data_[i];
  }

  return result;
}

// Norm implementation
template <typename ConfigBackend, typename GeometryBackend>
double WRFState<ConfigBackend, GeometryBackend>::norm() const {
  double sumSquares = 0.0;

  // Compute norm directly on flattened data
  for (double value : flattened_data_) {
    sumSquares += value * value;
  }

  return std::sqrt(sumSquares);
}

// Equals implementation
template <typename ConfigBackend, typename GeometryBackend>
bool WRFState<ConfigBackend, GeometryBackend>::equals(
    const WRFState<ConfigBackend, GeometryBackend>& other) const {
  if (!isCompatible(other)) {
    return false;
  }

  // Check if flattened data is equal within tolerance
  for (size_t i = 0; i < flattened_data_.size(); ++i) {
    if (std::abs(flattened_data_[i] - other.flattened_data_[i]) > 1e-10) {
      return false;
    }
  }

  return true;
}

// Add implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::add(
    const WRFState<ConfigBackend, GeometryBackend>& other) {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for addition");
  }

  // Add directly on flattened data
  for (size_t i = 0; i < flattened_data_.size(); ++i) {
    flattened_data_[i] += other.flattened_data_[i];
  }
}

// Subtract implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::subtract(
    const WRFState<ConfigBackend, GeometryBackend>& other) {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for subtraction");
  }

  // Subtract directly on flattened data
  for (size_t i = 0; i < flattened_data_.size(); ++i) {
    flattened_data_[i] -= other.flattened_data_[i];
  }
}

// Multiply implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::multiply(double scalar) {
  // Multiply flattened data by the scalar
  for (double& value : flattened_data_) {
    value *= scalar;
  }
}

// Is initialized implementation
template <typename ConfigBackend, typename GeometryBackend>
bool WRFState<ConfigBackend, GeometryBackend>::isInitialized() const {
  return initialized_;
}

// Private helper to load state data
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::loadStateData(
    const std::string& filename, const std::vector<std::string>& variables) {
  try {
    // Open NetCDF file
    netCDF::NcFile wrf_file(filename, netCDF::NcFile::read);

    if (!wrf_file.isNull()) {
      // Clear any existing data
      flattened_data_.clear();
      dimensions_.clear();
      variableNames_.clear();
      variable_grid_info_.clear();
      variable_offsets_.clear();

      // Load each requested variable
      for (const auto& varName : variables) {
        auto var = wrf_file.getVar(varName);

        if (!var.isNull()) {
          // Get variable dimensions
          const auto& varDims = var.getDims();
          std::vector<size_t> dims;

          // --- Begin: Fill VariableGridInfo ---
          VariableGridInfo grid_info;
          size_t dim_offset = 0;
          if (!varDims.empty() && varDims[0].getName() == "Time") {
            dim_offset = 1;
          }

          // Assume WRF variable dimensions are ordered as [Z, Y, X] after Time
          if (varDims.size() >= dim_offset + 2) {  // At least 2D (Y, X)
            if (varDims.size() >= dim_offset + 3) {
              // 3D variable: [Z, Y, X]
              grid_info.z_dim_name = varDims[dim_offset + 0].getName();
              grid_info.y_dim_name = varDims[dim_offset + 1].getName();
              grid_info.x_dim_name = varDims[dim_offset + 2].getName();
            } else {
              // 2D variable: [Y, X]
              grid_info.y_dim_name = varDims[dim_offset + 0].getName();
              grid_info.x_dim_name = varDims[dim_offset + 1].getName();
              grid_info.z_dim_name = "";  // No Z dimension
            }

            // Determine staggering based on dimension names
            grid_info.z_staggered =
                !grid_info.z_dim_name.empty() &&
                (grid_info.z_dim_name.find("_stag") != std::string::npos);
            grid_info.y_staggered =
                (grid_info.y_dim_name.find("_stag") != std::string::npos);
            grid_info.x_staggered =
                (grid_info.x_dim_name.find("_stag") != std::string::npos);

            // Determine grid type based on staggering pattern
            grid_info.grid_type = grid_info.determineGridType();

            // Validate grid type against WRF geometry
            try {
              switch (grid_info.grid_type) {
                case VariableGridInfo::GridType::UNSTAGGERED:
                  // Check if dimensions match unstaggered grid
                  if (grid_info.x_dim_name != "west_east" ||
                      grid_info.y_dim_name != "south_north") {
                    std::cerr << "Warning: Variable " << varName
                              << " classified as unstaggered but has "
                                 "unexpected dimensions: "
                              << grid_info.x_dim_name << ", "
                              << grid_info.y_dim_name << std::endl;
                  }
                  break;
                case VariableGridInfo::GridType::U_STAGGERED:
                  // Check if dimensions match U-staggered grid
                  if (grid_info.x_dim_name != "west_east_stag" ||
                      grid_info.y_dim_name != "south_north") {
                    std::cerr << "Warning: Variable " << varName
                              << " classified as U-staggered but has "
                                 "unexpected dimensions: "
                              << grid_info.x_dim_name << ", "
                              << grid_info.y_dim_name << std::endl;
                  }
                  break;
                case VariableGridInfo::GridType::V_STAGGERED:
                  // Check if dimensions match V-staggered grid
                  if (grid_info.x_dim_name != "west_east" ||
                      grid_info.y_dim_name != "south_north_stag") {
                    std::cerr << "Warning: Variable " << varName
                              << " classified as V-staggered but has "
                                 "unexpected dimensions: "
                              << grid_info.x_dim_name << ", "
                              << grid_info.y_dim_name << std::endl;
                  }
                  break;
                case VariableGridInfo::GridType::W_STAGGERED:
                  // Check if dimensions match W-staggered grid
                  if (grid_info.x_dim_name != "west_east" ||
                      grid_info.y_dim_name != "south_north" ||
                      grid_info.z_dim_name != "bottom_top_stag") {
                    std::cerr << "Warning: Variable " << varName
                              << " classified as W-staggered but has "
                                 "unexpected dimensions: "
                              << grid_info.x_dim_name << ", "
                              << grid_info.y_dim_name << ", "
                              << grid_info.z_dim_name << std::endl;
                  }
                  break;
              }
            } catch (const std::exception& e) {
              std::cerr << "Warning: Could not validate grid type for variable "
                        << varName << ": " << e.what() << std::endl;
            }

            // Debug output
            std::cout << "Variable " << varName
                      << " -> Grid type: " << grid_info.getGridTypeString()
                      << " (dims: " << grid_info.x_dim_name << ", "
                      << grid_info.y_dim_name;
            if (!grid_info.z_dim_name.empty()) {
              std::cout << ", " << grid_info.z_dim_name;
            }
            std::cout << ")" << std::endl;
          } else {
            // Variable has fewer than 2 dimensions, treat as unstaggered
            grid_info.grid_type = VariableGridInfo::GridType::UNSTAGGERED;
            std::cerr << "Warning: Variable " << varName
                      << " has fewer than 2 dimensions, treating as unstaggered"
                      << std::endl;
          }

          variable_grid_info_[varName] = grid_info;
          // --- End: Fill VariableGridInfo ---

          // Skip time dimension if present
          bool hasTimeDim =
              (varDims.size() > 0 && varDims[0].getName() == "Time");
          size_t startDim = hasTimeDim ? 1 : 0;

          for (size_t i = startDim; i < varDims.size(); ++i) {
            dims.push_back(varDims[i].getSize());
          }

          // Store dimensions
          dimensions_[varName] = dims;

          // Calculate total size
          size_t totalSize = 1;
          for (size_t dim : dims) {
            totalSize *= dim;
          }

          // Prepare start and count vectors for reading
          std::vector<size_t> start(varDims.size(), 0);
          std::vector<size_t> count(varDims.size(), 1);

          if (hasTimeDim) {
            start[0] = 0;  // or time_idx if you want to support multiple times
          }

          for (size_t i = startDim; i < varDims.size(); ++i) {
            count[i] = varDims[i].getSize();
          }

          // Allocate memory for the data
          std::vector<double> data(totalSize);

          // Read data
          var.getVar(start, count, data.data());

          // Calculate offset for this variable in flattened data
          size_t offset = flattened_data_.size();
          variable_offsets_[varName] = offset;

          // Append to flattened data
          flattened_data_.insert(flattened_data_.end(), data.begin(),
                                 data.end());

          // Store dimensions
          dimensions_[varName] = dims;
          variableNames_.push_back(varName);
        } else {
          std::cerr << "Warning: Variable not found in WRF file: " << varName
                    << std::endl;
        }
      }

      // Print summary of variable-to-grid associations
      if (!variableNames_.empty()) {
        std::cout << *this;
      }
    } else {
      throw std::runtime_error("Failed to open WRF file: " + filename);
    }
  } catch (const netCDF::exceptions::NcException& e) {
    throw std::runtime_error("NetCDF error while loading WRF state data: " +
                             std::string(e.what()));
  } catch (const std::exception& e) {
    throw std::runtime_error("Error loading WRF state data: " +
                             std::string(e.what()));
  }
}

// Helper to check compatibility
template <typename ConfigBackend, typename GeometryBackend>
bool WRFState<ConfigBackend, GeometryBackend>::isCompatible(
    const WRFState<ConfigBackend, GeometryBackend>& other) const {
  // Check if both states have the same variable names
  if (variableNames_.size() != other.variableNames_.size()) {
    return false;
  }

  // Check if each variable has the same dimensions
  for (const auto& varName : variableNames_) {
    // Check if variable exists in other state
    if (std::find(other.variableNames_.begin(), other.variableNames_.end(),
                  varName) == other.variableNames_.end()) {
      return false;
    }

    // Check dimensions
    const auto& thisDims = dimensions_.at(varName);
    const auto& otherDims = other.dimensions_.at(varName);

    if (thisDims != otherDims) {
      return false;
    }
  }

  return true;
}

// Implementation of getObsOperatorData
template <typename ConfigBackend, typename GeometryBackend>
typename WRFState<ConfigBackend, GeometryBackend>::ObsOperatorData
WRFState<ConfigBackend, GeometryBackend>::getObsOperatorData() const {
  ObsOperatorData data;

  // Set grid dimensions
  data.nx = static_cast<int>(geometry_.x_dim());
  data.ny = static_cast<int>(geometry_.y_dim());
  data.nz = static_cast<int>(geometry_.z_dim());

  // Get state variables
  data.u = static_cast<const double*>(getData("U"));
  data.v = static_cast<const double*>(getData("V"));
  data.t = static_cast<const double*>(getData("T"));
  data.q = static_cast<const double*>(getData("QVAPOR"));
  data.psfc = static_cast<const double*>(getData("PSFC"));
  data.ph = static_cast<const double*>(getData("PH"));
  data.phb = static_cast<const double*>(getData("PHB"));
  data.hgt = static_cast<const double*>(getData("HGT"));

  // Get grid metadata
  const auto& gi = geometry_.unstaggered_info();
  data.lats_2d = gi.latitude_2d.empty() ? nullptr : gi.latitude_2d.data();
  data.lons_2d = gi.longitude_2d.empty() ? nullptr : gi.longitude_2d.data();
  data.levels =
      gi.vertical_coords.empty() ? nullptr : gi.vertical_coords.data();

  // Set variable dimensions and staggered flags
  auto set_dims = [&](const std::string& var_name,
                      ObsOperatorData::VariableDims& dims) {
    try {
      const auto& var_dims = getVariableDimensions(var_name);
      if (var_dims.size() >= 3) {
        dims.nz = static_cast<int>(var_dims[0]);
        dims.ny = static_cast<int>(var_dims[1]);
        dims.nx = static_cast<int>(var_dims[2]);
      } else if (var_dims.size() == 2) {
        dims.ny = static_cast<int>(var_dims[0]);
        dims.nx = static_cast<int>(var_dims[1]);
        dims.nz = 1;
      } else {
        dims.nx = dims.ny = dims.nz = 1;
      }

      // Check if staggered
      const auto& grid_info = getVariableGridInfo(var_name);
      dims.is_staggered =
          (grid_info.grid_type != VariableGridInfo::GridType::UNSTAGGERED);
    } catch (...) {
      dims.nx = dims.ny = dims.nz = 1;
      dims.is_staggered = false;
    }
  };

  set_dims("U", data.u_dims);
  set_dims("V", data.v_dims);
  set_dims("T", data.t_dims);
  set_dims("QVAPOR", data.q_dims);
  set_dims("PSFC", data.psfc_dims);
  set_dims("PH", data.ph_dims);
  set_dims("PHB", data.phb_dims);
  set_dims("HGT", data.hgt_dims);

  // Validate required variables are present - fail fast if missing
  if (!data.u) {
    throw std::runtime_error(
        "Required variable 'U' (U-wind) not found in WRF state");
  }
  if (!data.v) {
    throw std::runtime_error(
        "Required variable 'V' (V-wind) not found in WRF state");
  }
  if (!data.t) {
    throw std::runtime_error(
        "Required variable 'T' (Temperature) not found in WRF state");
  }
  if (!data.q) {
    throw std::runtime_error(
        "Required variable 'QVAPOR' (Water vapor) not found in WRF state");
  }
  if (!data.psfc) {
    throw std::runtime_error(
        "Required variable 'PSFC' (Surface pressure) not found in WRF state");
  }
  if (!data.ph) {
    throw std::runtime_error(
        "Required variable 'PH' (Geopotential perturbation) not found in WRF "
        "state");
  }
  if (!data.phb) {
    throw std::runtime_error(
        "Required variable 'PHB' (Base state geopotential) not found in WRF "
        "state");
  }
  if (!data.hgt) {
    throw std::runtime_error(
        "Required variable 'HGT' (Terrain height) not found in WRF state");
  }

  // Calculate height field: hf = (ph + phb) / gravity
  // WRF gravity constant is typically 9.81 m/s²
  constexpr double gravity = 9.81;  // m/s²

  // PH, PB, and HF are vertically staggered variables with nz+1 levels
  const int staggered_nz = data.nz + 1;

  // Allocate memory for height field calculation
  static thread_local std::vector<double> height_field_storage;
  const size_t total_size =
      static_cast<size_t>(data.nx) * data.ny * staggered_nz;
  height_field_storage.resize(total_size);

  // Calculate height field for each grid point (including staggered levels)
  for (size_t k = 0; k < static_cast<size_t>(staggered_nz); ++k) {
    for (size_t j = 0; j < static_cast<size_t>(data.ny); ++j) {
      for (size_t i = 0; i < static_cast<size_t>(data.nx); ++i) {
        const size_t idx = k * data.ny * data.nx + j * data.nx + i;
        height_field_storage[idx] = (data.ph[idx] + data.phb[idx]) / gravity;
      }
    }
  }

  // Set height field pointer and dimensions
  data.hf = height_field_storage.data();
  data.hf_dims = data.ph_dims;  // Height field has same dimensions as PH/PHB
  data.hf_dims.nz =
      staggered_nz;  // Override with correct staggered vertical dimension

  return data;
}

}  // namespace metada::backends::wrf