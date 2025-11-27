/**
 * @file WRFState.hpp
 * @brief WRF state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <netcdf>
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
#include "wrfda/WRFDAStateTransforms.hpp"

// WRFDA C bindings for domain initialization
extern "C" {
// NEW: WRFDA First Guess Initialization Functions
// These functions implement the refactored WRFState initialization strategy
// that leverages WRFDA's proven initialization pipeline
int wrfda_load_first_guess(void* grid_ptr, const char* filename,
                           int filename_len);

int wrfda_extract_background_state(void* grid_ptr, double* u, double* v,
                                   double* t, double* q, double* psfc,
                                   double* p, double* ph, double* phb,
                                   double* hgt, double* lats, double* lons,
                                   int* nx, int* ny, int* nz);

int wrfda_extract_additional_fields(void* grid_ptr, double* w, double* mu,
                                    double* mub, double* pb, double* t_init);

// WARNING: DO NOT USE IN DESTRUCTORS!
// da_control module variables (c1f, c2f, etc.) are module-level and shared
// across all WRFState instances. Only call during final WRFDA shutdown.
void wrfda_cleanup_vertical_coords();

// LEGACY: Old domain construction methods (kept for backward compatibility)
int wrfda_construct_domain_from_arrays(
    const int* nx, const int* ny, const int* nz, const double* u,
    const double* v, const double* t, const double* q, const double* psfc,
    const double* ph, const double* phb, const double* hf, const double* hgt,
    const double* p, const double* pb, const double* lats2d,
    const double* lons2d);

int wrfda_init_domain_from_wrf_fields(
    const int* nx, const int* ny, const int* nz, const double* u,
    const double* v, const double* w, const double* t, const double* mu,
    const double* mub, const double* p, const double* pb, const double* ph,
    const double* phb, const double* xlat, const double* xlong,
    const double* ht, const double* znu, const double* znw, const double* dn,
    const double* dnw, const double* rdnw, const double* rdn,
    const double* p_top, const double* t_init, const double* moist,
    const int* num_moist, const double* psfc, const int* start_year,
    const int* start_month, const int* start_day, const int* start_hour);
}

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
  using reference = double&;              ///< Reference to element type
  using const_reference = const double&;  ///< Const reference to element type
  using size_type = std::size_t;          ///< Type for sizes and indices
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
   * @brief Destructor
   *
   * @details Uses default destructor. WRFDA module-level variables
   * (da_control::c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h) are NOT cleaned
   * up here as they are shared across all WRFState instances. WRFDA
   * manages these internally during finalization.
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
   * @brief Get mutable access to core state variables data (concept compliance)
   *
   * @details Returns a void pointer to the raw data array of core state
   * variables (U, V, T, QVAPOR, PSFC) stored contiguously at the beginning.
   * This method is provided for StateBackendImpl concept compliance.
   * For explicit variable access, use getData(variableName).
   *
   * @return void* Pointer to core state variables data, or nullptr if no core
   * variables
   * @throws std::runtime_error If no core variables are loaded
   *
   * @warning Use with caution - no bounds checking is performed
   * @see getData(const std::string&) for explicit variable access
   */
  void* getData();

  /**
   * @brief Get const access to core state variables data (concept compliance)
   *
   * @details Returns a const void pointer to the raw data array of core state
   * variables (U, V, T, QVAPOR, PSFC) stored contiguously at the beginning.
   * This method is provided for StateBackendImpl concept compliance.
   * For explicit variable access, use getData(variableName).
   *
   * @return const void* Const pointer to core state variables data, or nullptr
   * if no core variables
   * @throws std::runtime_error If no core variables are loaded
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

  ///@{ @name State Data Access
  /**
   * @brief Extract core state data for direct access
   *
   * @details This method extracts only the core state variables from the WRF
   * state, providing raw pointers for efficient access. This is a lightweight
   * alternative to ObsOperatorData when only state variables are needed.
   *
   * @return Structure containing core state data
   */
  struct StateData {
    // State variables (as raw pointers for efficiency)
    const double* u = nullptr;     // U-wind component
    const double* v = nullptr;     // V-wind component
    const double* t = nullptr;     // Temperature
    const double* q = nullptr;     // Specific humidity
    const double* psfc = nullptr;  // Surface pressure

    // Core state data locations in flattened_data_
    size_t u_begin = 0, u_end = 0;        // U-wind location range
    size_t v_begin = 0, v_end = 0;        // V-wind location range
    size_t t_begin = 0, t_end = 0;        // Temperature location range
    size_t q_begin = 0, q_end = 0;        // Specific humidity location range
    size_t psfc_begin = 0, psfc_end = 0;  // Surface pressure location range

    // Total core state size (u + v + t + q + psfc)
    size_t core_state_size = 0;
  };

  StateData getStateData() const;

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
    const double* w = nullptr;  // Vertical velocity (W-staggered)
    const double* t = nullptr;
    const double* q = nullptr;
    const double* psfc = nullptr;
    const double* ph = nullptr;   // Geopotential perturbation
    const double* phb = nullptr;  // Base state geopotential
    const double* hf = nullptr;   // Height field (calculated from PH and PHB)
    const double* hgt = nullptr;  // Terrain height (from HGT variable)
    const double* p = nullptr;    // Pressure perturbation
    const double* pb = nullptr;   // Base state pressure
    const double* mu = nullptr;   // Column dry air mass perturbation
    const double* mub = nullptr;  // Base state column dry air mass
    const double* t_init = nullptr;  // Initial temperature field

    // Grid metadata
    const double* lats_2d = nullptr;
    const double* lons_2d = nullptr;
    const double* levels = nullptr;

    // Variable dimensions (for staggered grids)
    struct VariableDims {
      int nx, ny, nz;
      bool is_staggered = false;
    };

    VariableDims u_dims, v_dims, w_dims, t_dims, q_dims, psfc_dims, ph_dims,
        phb_dims, hf_dims, hgt_dims, p_dims, pb_dims, mu_dims, mub_dims,
        t_init_dims;
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

  /**
   * @brief Add increment to this state using WRFDA's proven routines
   *
   * @details Adds analysis increments to background state by delegating to
   * WRFDA's da_transfer_xatowrf, which handles all WRF-specific
   * transformations:
   * - Conversion of specific humidity to mixing ratio
   * - Computation of dry air mass increments
   * - Conversion of temperature to potential temperature
   * - Computation of geopotential height
   * - Arakawa-C grid staggering (A-grid to C-grid conversion)
   * - Positivity constraints (e.g., moisture >= 0)
   * - Recomputation of 2m/10m diagnostic fields
   *
   * This is the fundamental operation for combining analysis increments
   * with background states: x_analysis = x_background + increment
   *
   * @tparam IncrementBackend The increment backend type
   * @param[in] increment Increment to add to this state
   *
   * @note This leverages WRFDA's battle-tested code instead of reimplementing
   * complex WRF-specific physics and grid transformations
   */
  template <typename IncrementBackend>
  void addIncrement(const IncrementBackend& increment);
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
   * @brief Save state data using WRFDA's da_update_firstguess routine
   *
   * @details Delegates to WRFDA's battle-tested NetCDF writer that copies the
   * background (fg) file and updates only the analysis variables touched by
   * assimilation. This ensures complete compatibility with WRF's expected file
   * structure and metadata handling.
   *
   * @param[in] filename Path where the updated analysis file should be written
   * @throws std::runtime_error If filesystem operations or WRFDA call fails
   */
  void saveToFile(const std::string& filename) const {
    try {
      if (filename.empty()) {
        throw std::runtime_error(
            "Output filename for WRF state cannot be empty");
      }
      void* grid_ptr = geometry_.getGridPtr();
      if (!grid_ptr) {
        throw std::runtime_error(
            "WRFDA grid pointer is null; cannot write analysis file");
      }
      std::filesystem::path working_dir = std::filesystem::current_path();
      std::filesystem::path fg_target = working_dir / "fg";
      std::filesystem::path fg_backup = working_dir / "fg.metada.bak";
      bool fg_preexisted = std::filesystem::exists(fg_target);

      auto restoreOriginalFg = [&]() noexcept {
        try {
          if (fg_preexisted) {
            if (std::filesystem::exists(fg_target)) {
              std::filesystem::remove(fg_target);
            }
            if (std::filesystem::exists(fg_backup)) {
              std::filesystem::rename(fg_backup, fg_target);
            }
          } else if (std::filesystem::exists(fg_target)) {
            std::filesystem::remove(fg_target);
          }
        } catch (...) {
        }
      };

      try {
        if (fg_preexisted) {
          std::filesystem::rename(fg_target, fg_backup);
        }
      } catch (const std::filesystem::filesystem_error& e) {
        throw std::runtime_error("Failed to backup existing fg file: " +
                                 std::string(e.what()));
      }

      try {
        std::filesystem::copy_file(
            wrfFilename_, fg_target,
            std::filesystem::copy_options::overwrite_existing);
      } catch (const std::filesystem::filesystem_error& e) {
        restoreOriginalFg();
        throw std::runtime_error("Failed to stage fg file for WRFDA: " +
                                 std::string(e.what()));
      }

      try {
        wrfda::WRFDAStateTransforms::updateFirstGuess(grid_ptr, filename);
      } catch (...) {
        restoreOriginalFg();
        throw;
      }

      restoreOriginalFg();
    } catch (const std::exception& e) {
      throw std::runtime_error(
          "Error saving WRF state via WRFDA da_update_firstguess: " +
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
    os << "WRFState{\n";
    os << "  initialized: " << (state.initialized_ ? "true" : "false") << ",\n";
    os << "  filename: \"" << state.wrfFilename_ << "\",\n";
    os << "  variables: [";

    for (size_t i = 0; i < state.variableNames_.size(); ++i) {
      if (i > 0) os << ", ";
      os << "\"" << state.variableNames_[i] << "\"";
    }
    os << "],\n";

    os << "  total_size: " << state.size();

    // Add statistics if data is available
    if (!state.flattened_data_.empty()) {
      os << ",\n  variable_stats: {\n";

      for (size_t i = 0; i < state.variableNames_.size(); ++i) {
        const std::string& varName = state.variableNames_[i];
        os << "    \"" << varName << "\": {";

        // Get variable offset and dimensions
        auto offset_it = state.variable_offsets_.find(varName);
        auto dims_it = state.dimensions_.find(varName);

        if (offset_it != state.variable_offsets_.end() &&
            dims_it != state.dimensions_.end()) {
          size_t offset = offset_it->second;
          size_t var_size = 1;
          for (size_t dim : dims_it->second) var_size *= dim;

          // Calculate statistics for this variable
          auto var_begin = state.flattened_data_.begin() + offset;
          auto var_end = var_begin + var_size;

          auto var_min_it = std::min_element(var_begin, var_end);
          auto var_max_it = std::max_element(var_begin, var_end);
          double var_sum = std::accumulate(var_begin, var_end, 0.0);
          double var_mean = var_sum / var_size;

          os << "min=" << std::fixed << std::setprecision(6) << *var_min_it;
          os << ", max=" << *var_max_it;
          os << ", mean=" << var_mean;
          os << ", size=" << var_size;
        } else {
          os << "invalid";
        }

        os << "}";
        if (i < state.variableNames_.size() - 1) os << ",";
        os << "\n";
      }

      os << "  }";
    }

    os << "\n}\n";
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
   * @brief Load 1D vertical coordinate variables from WRF NetCDF file
   *
   * @details Loads ZNU, ZNW, DN, DNW, RDNW, RDN arrays that define the
   * vertical coordinate system. These are typically 1D arrays stored
   * as coordinate variables in the NetCDF file.
   *
   * @param[in] filename Path to the WRF NetCDF file
   * @throws std::runtime_error If required coordinate variables are missing
   */
  void loadVerticalCoordinates(const std::string& filename);

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

  /**
   * @brief Initialize WRFDA domain structure from state data
   *
   * @details Constructs the WRFDA domain structure (persistent_grid) from the
   * loaded state data. This includes setting up grid dimensions, allocating
   * state arrays, copying data from C++ arrays to Fortran structures, and
   * initializing WRFDA module variables. This method should be called after
   * the state data is loaded and before any observation operator operations.
   *
   * @throws std::runtime_error If WRFDA domain construction fails
   * @throws std::runtime_error If required variables are missing
   */
  void initializeWRFDADomain();

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

  /// @name Core State Variables (stored contiguously at beginning)
  ///@{
  size_t u_begin_ = 0, u_end_ = 0;        ///< U-wind location range
  size_t v_begin_ = 0, v_end_ = 0;        ///< V-wind location range
  size_t t_begin_ = 0, t_end_ = 0;        ///< Temperature location range
  size_t q_begin_ = 0, q_end_ = 0;        ///< Specific humidity location range
  size_t psfc_begin_ = 0, psfc_end_ = 0;  ///< Surface pressure location range
  size_t core_state_size_ = 0;            ///< Total core state size
  ///@}

  /// @name Vertical Coordinate Arrays (1D)
  ///@{
  std::vector<double> znu_;  ///< Eta values on mass levels (size: nz)
  std::vector<double> znw_;  ///< Eta values on staggered levels (size: nz+1)
  std::vector<double> dn_;   ///< Delta eta on mass levels (size: nz)
  std::vector<double> dnw_;  ///< Delta eta on staggered levels (size: nz+1)
  std::vector<double>
      rdnw_;  ///< Inverse delta eta on staggered levels (size: nz+1)
  std::vector<double> rdn_;  ///< Inverse delta eta on mass levels (size: nz)
  double p_top_ = 0.0;       ///< Pressure at model top (Pa)
  ///@}
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

  // Check if we should use WRFDA-based initialization (recommended)
  bool use_wrfda_init = true;
  try {
    if (config.HasKey("use_wrfda_initialization")) {
      use_wrfda_init = config.Get("use_wrfda_initialization").asBool();
    }
  } catch (const std::exception&) {
    // Default to WRFDA initialization
  }

  if (use_wrfda_init) {
    // NEW APPROACH: Use WRFDA's proven initialization pipeline
    std::cout << "Using WRFDA first guess initialization for: " << wrfFilename_
              << std::endl;

    // Get grid pointer from geometry (allocated via WRFConfigManager)
    void* grid_ptr = geometry.getGridPtr();
    if (!grid_ptr) {
      throw std::runtime_error(
          "Failed to get WRFDA grid pointer. Ensure WRFConfigManager allocated "
          "domain.");
    }

    // Load first guess using WRFDA with pre-allocated head_grid
    int rc = wrfda_load_first_guess(grid_ptr, wrfFilename_.c_str(),
                                    static_cast<int>(wrfFilename_.length()));
    if (rc != 0) {
      throw std::runtime_error("Failed to load first guess with WRFDA, code " +
                               std::to_string(rc));
    }

    // Step 2: Get grid dimensions from geometry (already known from NetCDF)
    const int nx = static_cast<int>(geometry.x_dim());
    const int ny = static_cast<int>(geometry.y_dim());
    const int nz = static_cast<int>(geometry.z_dim());

    std::cout << "Grid dimensions: nx=" << nx << ", ny=" << ny << ", nz=" << nz
              << std::endl;

    // Step 3: Allocate temporary storage for background state extraction
    std::vector<double> u_temp(nx * ny * nz);
    std::vector<double> v_temp(nx * ny * nz);
    std::vector<double> t_temp(nx * ny * nz);
    std::vector<double> q_temp(nx * ny * nz);
    std::vector<double> psfc_temp(nx * ny);
    std::vector<double> p_temp(nx * ny * nz);
    std::vector<double> ph_temp(nx * ny * (nz + 1));
    std::vector<double> phb_temp(nx * ny * (nz + 1));
    std::vector<double> hgt_temp(nx * ny);
    std::vector<double> lats_temp(nx * ny);
    std::vector<double> lons_temp(nx * ny);

    // Extract full background state in a SINGLE call (dimensions already known)
    int nx_out = nx, ny_out = ny, nz_out = nz;  // For verification
    rc = wrfda_extract_background_state(
        grid_ptr, u_temp.data(), v_temp.data(), t_temp.data(), q_temp.data(),
        psfc_temp.data(), p_temp.data(), ph_temp.data(), phb_temp.data(),
        hgt_temp.data(), lats_temp.data(), lons_temp.data(), &nx_out, &ny_out,
        &nz_out);
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to extract background state from WRFDA, code " +
          std::to_string(rc));
    }

    // Verify dimensions match (sanity check)
    if (nx_out != nx || ny_out != ny || nz_out != nz) {
      throw std::runtime_error(
          "Dimension mismatch: Expected (" + std::to_string(nx) + "," +
          std::to_string(ny) + "," + std::to_string(nz) + "), got (" +
          std::to_string(nx_out) + "," + std::to_string(ny_out) + "," +
          std::to_string(nz_out) + ")");
    }

    // Step 4: Extract additional fields
    std::vector<double> w_temp(nx * ny * (nz + 1));
    std::vector<double> mu_temp(nx * ny);
    std::vector<double> mub_temp(nx * ny);
    std::vector<double> pb_temp(nx * ny * nz);
    std::vector<double> t_init_temp(nx * ny * nz);

    rc = wrfda_extract_additional_fields(grid_ptr, w_temp.data(),
                                         mu_temp.data(), mub_temp.data(),
                                         pb_temp.data(), t_init_temp.data());
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to extract additional fields from WRFDA, code " +
          std::to_string(rc));
    }

    // Step 5: Populate internal data structures
    // Core variables (stored contiguously at beginning)
    std::vector<std::string> core_variables = {"U", "V", "T", "QVAPOR", "PSFC"};

    for (const auto& varName : core_variables) {
      if (varName == "U") {
        variable_offsets_["U"] = flattened_data_.size();
        u_begin_ = flattened_data_.size();
        flattened_data_.insert(flattened_data_.end(), u_temp.begin(),
                               u_temp.end());
        u_end_ = flattened_data_.size();
        dimensions_["U"] = {static_cast<size_t>(nz), static_cast<size_t>(ny),
                            static_cast<size_t>(nx)};
        variableNames_.push_back("U");
      } else if (varName == "V") {
        variable_offsets_["V"] = flattened_data_.size();
        v_begin_ = flattened_data_.size();
        flattened_data_.insert(flattened_data_.end(), v_temp.begin(),
                               v_temp.end());
        v_end_ = flattened_data_.size();
        dimensions_["V"] = {static_cast<size_t>(nz), static_cast<size_t>(ny),
                            static_cast<size_t>(nx)};
        variableNames_.push_back("V");
      } else if (varName == "T") {
        variable_offsets_["T"] = flattened_data_.size();
        t_begin_ = flattened_data_.size();
        flattened_data_.insert(flattened_data_.end(), t_temp.begin(),
                               t_temp.end());
        t_end_ = flattened_data_.size();
        dimensions_["T"] = {static_cast<size_t>(nz), static_cast<size_t>(ny),
                            static_cast<size_t>(nx)};
        variableNames_.push_back("T");
      } else if (varName == "QVAPOR") {
        variable_offsets_["QVAPOR"] = flattened_data_.size();
        q_begin_ = flattened_data_.size();
        flattened_data_.insert(flattened_data_.end(), q_temp.begin(),
                               q_temp.end());
        q_end_ = flattened_data_.size();
        dimensions_["QVAPOR"] = {static_cast<size_t>(nz),
                                 static_cast<size_t>(ny),
                                 static_cast<size_t>(nx)};
        variableNames_.push_back("QVAPOR");
      } else if (varName == "PSFC") {
        variable_offsets_["PSFC"] = flattened_data_.size();
        psfc_begin_ = flattened_data_.size();
        flattened_data_.insert(flattened_data_.end(), psfc_temp.begin(),
                               psfc_temp.end());
        psfc_end_ = flattened_data_.size();
        dimensions_["PSFC"] = {static_cast<size_t>(ny),
                               static_cast<size_t>(nx)};
        variableNames_.push_back("PSFC");
      }
    }

    core_state_size_ = flattened_data_.size();

    // Add additional fields (after core state variables)
    variable_offsets_["P"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), p_temp.begin(), p_temp.end());
    dimensions_["P"] = {static_cast<size_t>(nz), static_cast<size_t>(ny),
                        static_cast<size_t>(nx)};
    variableNames_.push_back("P");

    variable_offsets_["PH"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), ph_temp.begin(),
                           ph_temp.end());
    dimensions_["PH"] = {static_cast<size_t>(nz + 1), static_cast<size_t>(ny),
                         static_cast<size_t>(nx)};
    variableNames_.push_back("PH");

    variable_offsets_["PHB"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), phb_temp.begin(),
                           phb_temp.end());
    dimensions_["PHB"] = {static_cast<size_t>(nz + 1), static_cast<size_t>(ny),
                          static_cast<size_t>(nx)};
    variableNames_.push_back("PHB");

    variable_offsets_["HGT"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), hgt_temp.begin(),
                           hgt_temp.end());
    dimensions_["HGT"] = {static_cast<size_t>(ny), static_cast<size_t>(nx)};
    variableNames_.push_back("HGT");

    variable_offsets_["XLAT"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), lats_temp.begin(),
                           lats_temp.end());
    dimensions_["XLAT"] = {static_cast<size_t>(ny), static_cast<size_t>(nx)};
    variableNames_.push_back("XLAT");

    variable_offsets_["XLONG"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), lons_temp.begin(),
                           lons_temp.end());
    dimensions_["XLONG"] = {static_cast<size_t>(ny), static_cast<size_t>(nx)};
    variableNames_.push_back("XLONG");

    variable_offsets_["W"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), w_temp.begin(), w_temp.end());
    dimensions_["W"] = {static_cast<size_t>(nz + 1), static_cast<size_t>(ny),
                        static_cast<size_t>(nx)};
    variableNames_.push_back("W");

    variable_offsets_["MU"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), mu_temp.begin(),
                           mu_temp.end());
    dimensions_["MU"] = {static_cast<size_t>(ny), static_cast<size_t>(nx)};
    variableNames_.push_back("MU");

    variable_offsets_["MUB"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), mub_temp.begin(),
                           mub_temp.end());
    dimensions_["MUB"] = {static_cast<size_t>(ny), static_cast<size_t>(nx)};
    variableNames_.push_back("MUB");

    variable_offsets_["PB"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), pb_temp.begin(),
                           pb_temp.end());
    dimensions_["PB"] = {static_cast<size_t>(nz), static_cast<size_t>(ny),
                         static_cast<size_t>(nx)};
    variableNames_.push_back("PB");

    variable_offsets_["T_INIT"] = flattened_data_.size();
    flattened_data_.insert(flattened_data_.end(), t_init_temp.begin(),
                           t_init_temp.end());
    dimensions_["T_INIT"] = {static_cast<size_t>(nz), static_cast<size_t>(ny),
                             static_cast<size_t>(nx)};
    variableNames_.push_back("T_INIT");

    // Load vertical coordinate arrays (1D) from NetCDF
    // WRFDA doesn't provide these through extraction, so load directly
    loadVerticalCoordinates(wrfFilename_);

    initialized_ = true;

    std::cout << "WRFDA-based initialization completed successfully"
              << std::endl;
    std::cout << "Loaded " << variableNames_.size()
              << " variables, core state size: " << core_state_size_
              << std::endl;
  } else {
    // LEGACY APPROACH: Manual NetCDF loading + WRFDA initialization
    std::cout << "Using legacy initialization for: " << wrfFilename_
              << std::endl;

    // Load state data from WRF NetCDF file
    // Core variables (u, v, t, q, psfc) will be stored contiguously at
    // beginning Only insert if not already present

    // Required fields for da_transfer_wrftoxb
    std::vector<std::string> required_fields = {
        // Geopotential (staggered vertical)
        "PH", "PHB",
        // Terrain height
        "HGT",
        // Pressure perturbation and base state
        "P", "PB",
        // Vertical velocity (staggered vertical)
        "W",
        // Column dry air mass (perturbation and base state)
        "MU", "MUB",
        // Initial temperature field
        "T_INIT",
        // Additional moisture species (beyond QVAPOR)
        "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUPEL"};

    for (const auto& var : required_fields) {
      if (std::find(variables.begin(), variables.end(), var) ==
          variables.end()) {
        variables.push_back(var);
      }
    }

    loadStateData(wrfFilename_, variables);

    // Load vertical coordinate arrays (1D)
    loadVerticalCoordinates(wrfFilename_);

    initialized_ = true;

    // Initialize WRFDA domain structure after state is fully initialized
    try {
      initializeWRFDADomain();
    } catch (const std::exception& e) {
      std::cerr << "Warning: WRFDA domain initialization failed: " << e.what()
                << std::endl;
      std::cerr << "WRFDA observation operators may not work correctly."
                << std::endl;
      // Continue with initialization - WRFDA is optional for basic state
      // operations
    }
  }
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
      variable_offsets_(std::move(other.variable_offsets_)),
      u_begin_(other.u_begin_),
      u_end_(other.u_end_),
      v_begin_(other.v_begin_),
      v_end_(other.v_end_),
      t_begin_(other.t_begin_),
      t_end_(other.t_end_),
      q_begin_(other.q_begin_),
      q_end_(other.q_end_),
      psfc_begin_(other.psfc_begin_),
      psfc_end_(other.psfc_end_),
      core_state_size_(other.core_state_size_) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.variableNames_.clear();
  other.core_state_size_ = 0;
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

    // Copy core state tracking variables
    u_begin_ = other.u_begin_;
    u_end_ = other.u_end_;
    v_begin_ = other.v_begin_;
    v_end_ = other.v_end_;
    t_begin_ = other.t_begin_;
    t_end_ = other.t_end_;
    q_begin_ = other.q_begin_;
    q_end_ = other.q_end_;
    psfc_begin_ = other.psfc_begin_;
    psfc_end_ = other.psfc_end_;
    core_state_size_ = other.core_state_size_;

    // Reset the moved-from object
    other.initialized_ = false;
    other.variableNames_.clear();
    other.core_state_size_ = 0;
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

  // Copy core state tracking variables
  cloned->u_begin_ = this->u_begin_;
  cloned->u_end_ = this->u_end_;
  cloned->v_begin_ = this->v_begin_;
  cloned->v_end_ = this->v_end_;
  cloned->t_begin_ = this->t_begin_;
  cloned->t_end_ = this->t_end_;
  cloned->q_begin_ = this->q_begin_;
  cloned->q_end_ = this->q_end_;
  cloned->psfc_begin_ = this->psfc_begin_;
  cloned->psfc_end_ = this->psfc_end_;
  cloned->core_state_size_ = this->core_state_size_;

  return cloned;
}

// Data access implementation (concept compliance) - returns only core state
// variables
template <typename ConfigBackend, typename GeometryBackend>
void* WRFState<ConfigBackend, GeometryBackend>::getData() {
  if (!initialized_ || core_state_size_ == 0) {
    return nullptr;
  }
  // Return pointer to core state variables (stored at beginning of
  // flattened_data_)
  return flattened_data_.data();
}

// Const data access implementation (concept compliance) - returns only core
// state variables
template <typename ConfigBackend, typename GeometryBackend>
const void* WRFState<ConfigBackend, GeometryBackend>::getData() const {
  if (!initialized_ || core_state_size_ == 0) {
    return nullptr;
  }
  // Return pointer to core state variables (stored at beginning of
  // flattened_data_)
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

// Size implementation - returns size of core state variables only
template <typename ConfigBackend, typename GeometryBackend>
size_t WRFState<ConfigBackend, GeometryBackend>::size() const {
  return core_state_size_;
}

// Zero implementation - only zeros core state variables (U, V, T, QVAPOR, PSFC)
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::zero() {
  // Only zero the core state variables which are stored contiguously at the
  // beginning
  if (core_state_size_ > 0) {
    std::fill(flattened_data_.begin(),
              flattened_data_.begin() + core_state_size_, 0.0);
  }
}

// Dot product implementation - only operates on core state variables
template <typename ConfigBackend, typename GeometryBackend>
double WRFState<ConfigBackend, GeometryBackend>::dot(
    const WRFState<ConfigBackend, GeometryBackend>& other) const {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for dot product");
  }

  double result = 0.0;

  // Compute dot product only on core state variables
  for (size_t i = 0; i < core_state_size_; ++i) {
    result += flattened_data_[i] * other.flattened_data_[i];
  }

  return result;
}

// Norm implementation - only operates on core state variables
template <typename ConfigBackend, typename GeometryBackend>
double WRFState<ConfigBackend, GeometryBackend>::norm() const {
  double sumSquares = 0.0;

  // Compute norm only on core state variables
  for (size_t i = 0; i < core_state_size_; ++i) {
    sumSquares += flattened_data_[i] * flattened_data_[i];
  }

  return std::sqrt(sumSquares);
}

// Equals implementation - only compares core state variables
template <typename ConfigBackend, typename GeometryBackend>
bool WRFState<ConfigBackend, GeometryBackend>::equals(
    const WRFState<ConfigBackend, GeometryBackend>& other) const {
  if (!isCompatible(other)) {
    return false;
  }

  // Check if core state variables are equal within tolerance
  for (size_t i = 0; i < core_state_size_; ++i) {
    if (std::abs(flattened_data_[i] - other.flattened_data_[i]) > 1e-10) {
      return false;
    }
  }

  return true;
}

// Add implementation - only operates on core state variables
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::add(
    const WRFState<ConfigBackend, GeometryBackend>& other) {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for addition");
  }

  // Add only core state variables
  for (size_t i = 0; i < core_state_size_; ++i) {
    flattened_data_[i] += other.flattened_data_[i];
  }
}

// Subtract implementation - only operates on core state variables
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::subtract(
    const WRFState<ConfigBackend, GeometryBackend>& other) {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for subtraction");
  }

  // Subtract only core state variables
  for (size_t i = 0; i < core_state_size_; ++i) {
    flattened_data_[i] -= other.flattened_data_[i];
  }
}

// Multiply implementation - only operates on core state variables
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::multiply(double scalar) {
  // Multiply only core state variables by the scalar
  for (size_t i = 0; i < core_state_size_; ++i) {
    flattened_data_[i] *= scalar;
  }
}

// Add increment implementation
template <typename ConfigBackend, typename GeometryBackend>
template <typename IncrementBackend>
void WRFState<ConfigBackend, GeometryBackend>::addIncrement(
    const IncrementBackend& increment) {
  // Use WRFDA's proven da_transfer_xatowrf routine to add increments
  // This handles all WRF-specific physics and grid transformations:
  // - Conversion of specific humidity to mixing ratio
  // - Computation of dry air mass increments
  // - Conversion of temperature to potential temperature (theta)
  // - Computation of geopotential height for WRF's vertical coordinate
  // - Arakawa-C grid staggering (A-grid to C-grid conversion)
  // - Positivity constraints (ensures moisture >= 0, etc.)
  // - Recomputation of 2m/10m diagnostic fields (T2, Q2, U10, V10, TH2)

  // Step 1: Sync increment to grid%xa (WRFDA's analysis increment workspace)
  // This writes the increment data into the WRFDA grid structure
  increment.syncToGrid();

  // Step 2: Call WRFDA's da_transfer_xatowrf to add increment to background
  // This performs the complete transformation: x_analysis = x_background +
  // increment After this call, the WRF grid structure (grid%u_2, grid%v_2,
  // grid%t_2, etc.) contains the full analysis state ready for WRF model
  // integration
  try {
    wrfda::WRFDAStateTransforms::transferXaToWRF(
        geometry_.getGridPtr(),  // WRF grid structure (contains xb and xa)
        geometry_.getConfigFlagsPtr()  // WRF configuration flags
    );
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::string(
            "Failed to add increment using WRFDA da_transfer_xatowrf: ") +
        e.what());
  }

  // Note: After da_transfer_xatowrf, the WRF grid structure already contains
  // the full analysis state in grid%u_2, grid%v_2, grid%t_2, grid%moist, etc.
  // These are the authoritative state variables that WRFState's getData()
  // methods reference, so no additional extraction is needed.
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

      // Reset core state tracking
      u_begin_ = u_end_ = v_begin_ = v_end_ = t_begin_ = t_end_ = q_begin_ =
          q_end_ = psfc_begin_ = psfc_end_ = 0;
      core_state_size_ = 0;

      // Define core state variables that will be stored contiguously at
      // beginning
      const std::vector<std::string> core_variables = {"U", "V", "T", "QVAPOR",
                                                       "PSFC"};

      // First pass: Load core variables and store them contiguously
      for (const auto& varName : core_variables) {
        if (std::find(variables.begin(), variables.end(), varName) !=
            variables.end()) {
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

            // Assume WRF variable dimensions are ordered as [Z, Y, X] after
            // Time
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
                std::cerr
                    << "Warning: Could not validate grid type for variable "
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
              std::cerr
                  << "Warning: Variable " << varName
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
              start[0] =
                  0;  // or time_idx if you want to support multiple times
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

            // Track core state variable locations
            size_t var_end = flattened_data_.size();
            if (varName == "U") {
              u_begin_ = offset;
              u_end_ = var_end;
            } else if (varName == "V") {
              v_begin_ = offset;
              v_end_ = var_end;
            } else if (varName == "T") {
              t_begin_ = offset;
              t_end_ = var_end;
            } else if (varName == "QVAPOR") {
              q_begin_ = offset;
              q_end_ = var_end;
            } else if (varName == "PSFC") {
              psfc_begin_ = offset;
              psfc_end_ = var_end;
            }

            // Store dimensions
            dimensions_[varName] = dims;
            variableNames_.push_back(varName);
          } else {
            std::cerr << "Warning: Variable not found in WRF file: " << varName
                      << std::endl;
          }
        }
      }

      // Calculate total core state size after loading all core variables
      // Core variables are stored contiguously at the beginning
      core_state_size_ = flattened_data_.size();

      // Second pass: Load non-core variables (diagnostic/background fields)
      for (const auto& varName : variables) {
        // Skip if already loaded as core variable
        if (std::find(core_variables.begin(), core_variables.end(), varName) !=
            core_variables.end()) {
          continue;
        }

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

          // Append to flattened data (after core variables)
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
  data.w = static_cast<const double*>(getData("W"));
  data.t = static_cast<const double*>(getData("T"));
  data.q = static_cast<const double*>(getData("QVAPOR"));
  data.psfc = static_cast<const double*>(getData("PSFC"));
  data.ph = static_cast<const double*>(getData("PH"));
  data.phb = static_cast<const double*>(getData("PHB"));
  data.hgt = static_cast<const double*>(getData("HGT"));
  data.p = static_cast<const double*>(getData("P"));
  data.pb = static_cast<const double*>(getData("PB"));
  data.mu = static_cast<const double*>(getData("MU"));
  data.mub = static_cast<const double*>(getData("MUB"));
  data.t_init = static_cast<const double*>(getData("T_INIT"));

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
  set_dims("W", data.w_dims);
  set_dims("T", data.t_dims);
  set_dims("QVAPOR", data.q_dims);
  set_dims("PSFC", data.psfc_dims);
  set_dims("PH", data.ph_dims);
  set_dims("PHB", data.phb_dims);
  set_dims("HGT", data.hgt_dims);
  set_dims("P", data.p_dims);
  set_dims("PB", data.pb_dims);
  set_dims("MU", data.mu_dims);
  set_dims("MUB", data.mub_dims);
  set_dims("T_INIT", data.t_init_dims);

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
  if (!data.p) {
    throw std::runtime_error(
        "Required variable 'P' (Pressure perturbation) not found in WRF "
        "state");
  }
  if (!data.pb) {
    throw std::runtime_error(
        "Required variable 'PB' (Base state pressure) not found in WRF "
        "state");
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

// Implementation of getStateData
template <typename ConfigBackend, typename GeometryBackend>
typename WRFState<ConfigBackend, GeometryBackend>::StateData
WRFState<ConfigBackend, GeometryBackend>::getStateData() const {
  StateData data;

  // Get state variables from contiguous storage at beginning of
  // flattened_data_
  data.u = flattened_data_.data() + u_begin_;
  data.v = flattened_data_.data() + v_begin_;
  data.t = flattened_data_.data() + t_begin_;
  data.q = flattened_data_.data() + q_begin_;
  data.psfc = flattened_data_.data() + psfc_begin_;

  // Set location ranges
  data.u_begin = u_begin_;
  data.u_end = u_end_;
  data.v_begin = v_begin_;
  data.v_end = v_end_;
  data.t_begin = t_begin_;
  data.t_end = t_end_;
  data.q_begin = q_begin_;
  data.q_end = q_end_;
  data.psfc_begin = psfc_begin_;
  data.psfc_end = psfc_end_;

  // Set total core state size
  data.core_state_size = core_state_size_;

  return data;
}

template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::initializeWRFDADomain() {
  // Use wrfda_init_domain_from_wrf_fields which properly leverages
  // WRFDA's da_transfer_wrftoxb subroutine for complete initialization

  // Get observation operator data which contains all required state variables
  const auto state_data = getObsOperatorData();

  // Get grid dimensions
  const int nx = static_cast<int>(geometry_.x_dim());
  const int ny = static_cast<int>(geometry_.y_dim());
  const int nz = static_cast<int>(geometry_.z_dim());

  // Get pointers to all required WRF fields from state_data
  // Note: getObsOperatorData provides u, v, t, q, psfc, ph, phb, hf, hgt, p,
  // pb, lats_2d, lons_2d For additional fields not in state_data, use data()
  // method
  const double* w_data = state_data.w;
  const double* mu_data = state_data.mu;
  const double* mub_data = state_data.mub;
  const double* t_init_data = state_data.t_init;

  // Get moisture species - for now just use QVAPOR from state_data
  // TODO: Combine all moisture species (QCLOUD, QRAIN, etc.) into moist array
  int num_moist = 1;  // Only QVAPOR for now

  // Get time information from config or WRF file
  // For now, use placeholder values - should come from WRF file metadata
  int start_year = 2000;
  int start_month = 1;
  int start_day = 24;
  int start_hour = 12;

  // Call wrfda_init_domain_from_wrf_fields with all required parameters
  int rc = wrfda_init_domain_from_wrf_fields(
      &nx, &ny, &nz, state_data.u, state_data.v, w_data, state_data.t, mu_data,
      mub_data, state_data.p, state_data.pb, state_data.ph, state_data.phb,
      state_data.lats_2d, state_data.lons_2d, state_data.hgt, znu_.data(),
      znw_.data(), dn_.data(), dnw_.data(), rdnw_.data(), rdn_.data(), &p_top_,
      t_init_data, state_data.q, &num_moist, state_data.psfc, &start_year,
      &start_month, &start_day, &start_hour);

  if (rc != 0) {
    throw std::runtime_error(
        "Failed to initialize WRFDA domain with da_transfer_wrftoxb, code " +
        std::to_string(rc));
  }
}

template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::loadVerticalCoordinates(
    const std::string& filename) {
  try {
    // Open NetCDF file
    netCDF::NcFile wrf_file(filename, netCDF::NcFile::read);

    if (wrf_file.isNull()) {
      throw std::runtime_error("Failed to open WRF file: " + filename);
    }

    // Get vertical dimension size
    auto nz_dim = wrf_file.getDim("bottom_top");
    if (nz_dim.isNull()) {
      throw std::runtime_error(
          "Cannot find 'bottom_top' dimension in WRF file");
    }
    size_t nz = nz_dim.getSize();
    size_t nz_stag = nz + 1;

    // Load ZNU (eta values on mass levels)
    auto znu_var = wrf_file.getVar("ZNU");
    if (!znu_var.isNull()) {
      znu_.resize(nz);
      znu_var.getVar(znu_.data());
    } else {
      throw std::runtime_error("Required variable ZNU not found in WRF file");
    }

    // Load ZNW (eta values on staggered levels)
    auto znw_var = wrf_file.getVar("ZNW");
    if (!znw_var.isNull()) {
      znw_.resize(nz_stag);
      znw_var.getVar(znw_.data());
    } else {
      throw std::runtime_error("Required variable ZNW not found in WRF file");
    }

    // Load DN (delta eta on mass levels)
    auto dn_var = wrf_file.getVar("DN");
    if (!dn_var.isNull()) {
      dn_.resize(nz);
      dn_var.getVar(dn_.data());
    } else {
      throw std::runtime_error("Required variable DN not found in WRF file");
    }

    // Load DNW (delta eta on staggered levels)
    auto dnw_var = wrf_file.getVar("DNW");
    if (!dnw_var.isNull()) {
      dnw_.resize(nz_stag);
      dnw_var.getVar(dnw_.data());
    } else {
      throw std::runtime_error("Required variable DNW not found in WRF file");
    }

    // Load RDNW (inverse delta eta on staggered levels)
    auto rdnw_var = wrf_file.getVar("RDNW");
    if (!rdnw_var.isNull()) {
      rdnw_.resize(nz_stag);
      rdnw_var.getVar(rdnw_.data());
    } else {
      throw std::runtime_error("Required variable RDNW not found in WRF file");
    }

    // Load RDN (inverse delta eta on mass levels)
    auto rdn_var = wrf_file.getVar("RDN");
    if (!rdn_var.isNull()) {
      rdn_.resize(nz);
      rdn_var.getVar(rdn_.data());
    } else {
      throw std::runtime_error("Required variable RDN not found in WRF file");
    }

    // Load P_TOP (pressure at model top)
    auto p_top_var = wrf_file.getVar("P_TOP");
    if (!p_top_var.isNull()) {
      p_top_var.getVar(&p_top_);
    } else {
      throw std::runtime_error("Required variable P_TOP not found in WRF file");
    }

    std::cout << "Loaded vertical coordinates: nz=" << nz
              << ", p_top=" << p_top_ << " Pa" << std::endl;

  } catch (const netCDF::exceptions::NcException& e) {
    throw std::runtime_error("NetCDF error loading vertical coordinates: " +
                             std::string(e.what()));
  }
}

}  // namespace metada::backends::wrf