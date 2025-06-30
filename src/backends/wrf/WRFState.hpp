/**
 * @file WRFState.hpp
 * @brief WRF state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

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
#include "Location.hpp"

namespace metada::backends::wrf {

// --- Begin: VariableGridInfo struct ---
struct VariableGridInfo {
  std::string x_dim_name;
  std::string y_dim_name;
  std::string z_dim_name;
  bool x_staggered = false;
  bool y_staggered = false;
  bool z_staggered = false;

  // Grid type association with WRFGeometry
  enum class GridType {
    UNSTAGGERED,  // Uses unstaggered_grid() from WRFGeometry
    U_STAGGERED,  // Uses u_staggered_grid() from WRFGeometry
    V_STAGGERED,  // Uses v_staggered_grid() from WRFGeometry
    W_STAGGERED   // Uses w_staggered_grid() from WRFGeometry
  };
  GridType grid_type = GridType::UNSTAGGERED;

  // Helper method to determine grid type
  GridType determineGridType() const {
    if (x_staggered && !y_staggered && !z_staggered) {
      return GridType::U_STAGGERED;  // Staggered in X (U wind component)
    } else if (!x_staggered && y_staggered && !z_staggered) {
      return GridType::V_STAGGERED;  // Staggered in Y (V wind component)
    } else if (!x_staggered && !y_staggered && z_staggered) {
      return GridType::W_STAGGERED;  // Staggered in Z (W wind component)
    } else {
      return GridType::UNSTAGGERED;  // Not staggered or multiple staggering
                                     // (rare)
    }
  }

  // Get grid type as string for debugging
  std::string getGridTypeString() const {
    switch (grid_type) {
      case GridType::UNSTAGGERED:
        return "unstaggered";
      case GridType::U_STAGGERED:
        return "u_staggered";
      case GridType::V_STAGGERED:
        return "v_staggered";
      case GridType::W_STAGGERED:
        return "w_staggered";
      default:
        return "unknown";
    }
  }
};
// --- End: VariableGridInfo struct ---

/**
 * @brief WRF state backend implementation
 *
 * @details
 * This class implements a state backend for the WRF model. It manages
 * meteorological state variables stored in NetCDF files and provides
 * operations required by the State adapter.
 */
template <typename ConfigBackend, typename GeometryBackend>
class WRFState {
 public:
  /**
   * @brief Default constructor is deleted
   */
  WRFState() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  WRFState(const WRFState&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  WRFState& operator=(const WRFState&) = delete;

  /**
   * @brief Constructor that takes a configuration backend and geometry backend
   *
   * @param config Configuration containing WRF file path and variables
   * @param geometry Geometry backend providing grid information
   */
  WRFState(const ConfigBackend& config, const GeometryBackend& geometry);

  /**
   * @brief Move constructor
   *
   * @param other WRF state backend to move from
   */
  WRFState(WRFState&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other WRF state backend to move from
   * @return Reference to this state after assignment
   */
  WRFState& operator=(WRFState&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~WRFState() = default;

  /**
   * @brief Clone this state
   *
   * @return Unique pointer to a new identical WRF state backend
   */
  std::unique_ptr<WRFState> clone() const;

  /**
   * @brief Get mutable access to the underlying data
   *
   * @return Void pointer to the active data array
   */
  void* getData();

  /**
   * @brief Get const access to the underlying data
   *
   * @return Const void pointer to the active data array
   */
  const void* getData() const;

  /**
   * @brief Get the active variable name
   *
   * @return Name of the currently active variable
   */
  const std::string& getActiveVariable() const;

  /**
   * @brief Set the active variable
   *
   * @param name Name of the variable to set as active
   * @throws std::out_of_range If variable doesn't exist
   */
  void setActiveVariable(const std::string& name);

  /**
   * @brief Get all variable names in this state
   *
   * @return Vector of variable names
   */
  const std::vector<std::string>& getVariableNames() const;

  /**
   * @brief Get the total size of the state vector
   *
   * @return Total number of elements in the state vector
   */
  size_t size() const;

  /**
   * @brief Set all values to zero
   */
  void zero();

  /**
   * @brief Calculate dot product with another state
   *
   * @param other State to calculate dot product with
   * @return Scalar dot product result
   * @throws std::runtime_error If states are incompatible
   */
  double dot(const WRFState& other) const;

  /**
   * @brief Calculate L2 norm of this state
   *
   * @return L2 norm value
   */
  double norm() const;

  /**
   * @brief Check if this state equals another
   *
   * @param other State to compare with
   * @return True if states are equal, false otherwise
   */
  bool equals(const WRFState& other) const;

  /**
   * @brief Add another state to this one
   *
   * @param other State to add
   * @throws std::runtime_error If states are incompatible
   */
  void add(const WRFState& other);

  /**
   * @brief Subtract another state from this one
   *
   * @param other State to subtract
   * @throws std::runtime_error If states are incompatible
   */
  void subtract(const WRFState& other);

  /**
   * @brief Multiply this state by a scalar
   *
   * @param scalar Value to multiply by
   */
  void multiply(double scalar);

  /**
   * @brief Check if state is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const;

  double& at(const framework::Location& loc) {
    // Convert geographic coordinates to grid indices
    auto [lat, lon, level] = loc.getGeographicCoords();

    // Get the appropriate geometry grid for the active variable
    const auto& geometry_grid = getVariableGeometryGrid(activeVariable_);

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
    auto& arr = variables_.at(activeVariable_);
    if (arr.dimension() == 3) {
      return arr(k, j, i);  // 3D: [Z, Y, X]
    } else {
      return arr(j, i);  // 2D: [Y, X]
    }
  }

  const double& at(const framework::Location& loc) const {
    // Convert geographic coordinates to grid indices
    auto [lat, lon, level] = loc.getGeographicCoords();

    // Get the appropriate geometry grid for the active variable
    const auto& geometry_grid = getVariableGeometryGrid(activeVariable_);

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
    const auto& arr = variables_.at(activeVariable_);
    if (arr.dimension() == 3) {
      return arr(k, j, i);  // 3D: [Z, Y, X]
    } else {
      return arr(j, i);  // 2D: [Y, X]
    }
  }

  // --- Begin: Additions for StateBackendType/Impl concept compliance ---
  using value_type = double;
  using reference = double&;
  using const_reference = const double&;
  using pointer = double*;
  using const_pointer = const double*;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  // --- End: Additions for StateBackendType/Impl concept compliance ---

  // --- Begin: More additions for StateBackendType/Impl concept compliance ---
  bool empty() const { return size() == 0; }
  size_type max_size() const { return size(); }

  reference operator[](size_type idx) {
    return variables_.at(activeVariable_).flat(idx);
  }
  const_reference operator[](size_type idx) const {
    return variables_.at(activeVariable_).flat(idx);
  }
  reference at(size_type idx) {
    return variables_.at(activeVariable_).flat(idx);
  }
  const_reference at(size_type idx) const {
    return variables_.at(activeVariable_).flat(idx);
  }
  reference front() { return variables_.at(activeVariable_).flat(0); }
  const_reference front() const {
    return variables_.at(activeVariable_).flat(0);
  }
  reference back() { return variables_.at(activeVariable_).flat(size() - 1); }
  const_reference back() const {
    return variables_.at(activeVariable_).flat(size() - 1);
  }
  // --- End: More additions for StateBackendType/Impl concept compliance ---

  // --- Begin: Add saveToFile for StateBackendImpl concept compliance ---
  void saveToFile(const std::string& filename) const {
    try {
      // Create NetCDF file
      netCDF::NcFile nc_file(filename, netCDF::NcFile::replace);

      if (nc_file.isNull()) {
        throw std::runtime_error("Failed to create NetCDF file: " + filename);
      }

      // Save each variable to the NetCDF file
      for (const auto& varName : variableNames_) {
        const auto& var_data = variables_.at(varName);
        const auto& var_dims = dimensions_.at(varName);

        // Create dimensions
        std::vector<netCDF::NcDim> nc_dims;
        for (size_t i = 0; i < var_dims.size(); ++i) {
          std::string dim_name = "dim_" + std::to_string(i) + "_" + varName;
          nc_dims.push_back(nc_file.addDim(dim_name, var_dims[i]));
        }

        // Create variable
        netCDF::NcVar nc_var =
            nc_file.addVar(varName, netCDF::ncDouble, nc_dims);

        // Write data
        nc_var.putVar(var_data.data());
      }

      // Add global attributes
      nc_file.putAtt("title", "WRF State Data");
      nc_file.putAtt("source", "Metada Framework");
      nc_file.putAtt("active_variable", activeVariable_);
      nc_file.putAtt("total_size", netCDF::ncInt, size());

      nc_file.close();

    } catch (const netCDF::exceptions::NcException& e) {
      throw std::runtime_error("NetCDF error while saving WRF state: " +
                               std::string(e.what()));
    } catch (const std::exception& e) {
      throw std::runtime_error("Error saving WRF state: " +
                               std::string(e.what()));
    }
  }
  // --- End: Add saveToFile for StateBackendImpl concept compliance ---

  const GeometryBackend& geometry() const { return geometry_; }

  // --- Variable-to-grid mapping ---
  std::unordered_map<std::string, VariableGridInfo> variable_grid_info_;

  /**
   * @brief Get grid information for a specific variable
   * @param variableName Name of the variable
   * @return VariableGridInfo structure containing grid association
   * @throws std::out_of_range If variable doesn't exist
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
   * @param variableName Name of the variable
   * @return Reference to the GridDimensionInfo from WRFGeometry
   * @throws std::out_of_range If variable doesn't exist
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
   * @param gridType The grid type to search for
   * @return Vector of variable names using the specified grid type
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
   * @param variableName Name of the variable to validate
   * @return True if dimensions match, false otherwise
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

  /**
   * @brief Stream operator for WRFState
   * @param os Output stream
   * @param state WRFState to output
   * @return Reference to the output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const WRFState& state) {
    os << "WRFState{";
    os << "initialized: " << (state.initialized_ ? "true" : "false");
    os << ", filename: \"" << state.wrfFilename_ << "\"";
    os << ", active_variable: \"" << state.activeVariable_ << "\"";
    os << ", variables: [";

    for (size_t i = 0; i < state.variableNames_.size(); ++i) {
      if (i > 0) os << ", ";
      os << "\"" << state.variableNames_[i] << "\"";
    }
    os << "]";

    os << ", total_size: " << state.size();
    os << "}";
    return os;
  }

  // --- End variable-to-grid mapping ---

 private:
  /**
   * @brief Private constructor for cloning (skips file loading)
   *
   * @param config Configuration backend
   * @param geometry Geometry backend
   * @param skip_loading Flag to skip file loading (for cloning)
   */
  WRFState(const ConfigBackend& config, const GeometryBackend& geometry,
           bool skip_loading);

  /**
   * @brief Static factory method for creating a clone
   * @param config Configuration backend
   * @param geometry Geometry backend
   * @return Unique pointer to a new WRFState ready for cloning
   */
  static std::unique_ptr<WRFState> createForCloning(
      const ConfigBackend& config, const GeometryBackend& geometry) {
    return std::unique_ptr<WRFState>(new WRFState(config, geometry, true));
  }

  /**
   * @brief Load state data from the WRF NetCDF file
   *
   * @param filename Path to the WRF NetCDF file
   * @param variables List of variables to load
   */
  void loadStateData(const std::string& filename,
                     const std::vector<std::string>& variables);

  /**
   * @brief Verify that two states have compatible variables and dimensions
   *
   * @param other State to check compatibility with
   * @return True if compatible, false otherwise
   */
  bool isCompatible(const WRFState& other) const;

  const ConfigBackend& config_;
  const GeometryBackend& geometry_;

  // WRF NetCDF file information
  std::string wrfFilename_;
  bool initialized_ = false;

  // State data
  std::unordered_map<std::string, xt::xarray<double>> variables_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
  std::vector<std::string> variableNames_;
  std::string activeVariable_;  // Currently active variable for data access
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

  // Set the active variable to the first one if available
  if (!variables.empty()) {
    activeVariable_ = variables[0];
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
      variables_(std::move(other.variables_)),
      dimensions_(std::move(other.dimensions_)),
      variableNames_(std::move(other.variableNames_)),
      activeVariable_(std::move(other.activeVariable_)) {
  // Reset the moved-from object
  other.initialized_ = false;
  other.variableNames_.clear();
  other.activeVariable_.clear();
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
    variables_ = std::move(other.variables_);
    dimensions_ = std::move(other.dimensions_);
    variableNames_ = std::move(other.variableNames_);
    activeVariable_ = std::move(other.activeVariable_);

    // Reset the moved-from object
    other.initialized_ = false;
    other.variableNames_.clear();
    other.activeVariable_.clear();
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
  cloned->variables_ = this->variables_;
  cloned->dimensions_ = this->dimensions_;
  cloned->variableNames_ = this->variableNames_;
  cloned->activeVariable_ = this->activeVariable_;
  cloned->variable_grid_info_ = this->variable_grid_info_;

  return cloned;
}

// Data access implementation
template <typename ConfigBackend, typename GeometryBackend>
void* WRFState<ConfigBackend, GeometryBackend>::getData() {
  if (!initialized_ || activeVariable_.empty()) {
    return nullptr;
  }

  try {
    return variables_.at(activeVariable_).data();
  } catch (const std::out_of_range&) {
    return nullptr;
  }
}

// Const data access implementation
template <typename ConfigBackend, typename GeometryBackend>
const void* WRFState<ConfigBackend, GeometryBackend>::getData() const {
  if (!initialized_ || activeVariable_.empty()) {
    return nullptr;
  }

  try {
    return variables_.at(activeVariable_).data();
  } catch (const std::out_of_range&) {
    return nullptr;
  }
}

// Get active variable implementation
template <typename ConfigBackend, typename GeometryBackend>
const std::string& WRFState<ConfigBackend, GeometryBackend>::getActiveVariable()
    const {
  return activeVariable_;
}

// Set active variable implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::setActiveVariable(
    const std::string& name) {
  if (variables_.find(name) == variables_.end()) {
    throw std::out_of_range("Variable not found: " + name);
  }
  activeVariable_ = name;
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
  size_t totalSize = 0;

  // Sum the sizes of all variables
  for (const auto& [name, data] : variables_) {
    totalSize += data.size();
  }

  return totalSize;
}

// Zero implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::zero() {
  for (auto& [name, data] : variables_) {
    data.fill(0.0);
  }
}

// Dot product implementation
template <typename ConfigBackend, typename GeometryBackend>
double WRFState<ConfigBackend, GeometryBackend>::dot(
    const WRFState<ConfigBackend, GeometryBackend>& other) const {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for dot product");
  }

  double result = 0.0;

  // Sum dot products of all variables
  for (const auto& varName : variableNames_) {
    const auto& thisVar = variables_.at(varName);
    const auto& otherVar = other.variables_.at(varName);

    // Use xtensor to compute element-wise multiplication and sum
    result += xt::sum(thisVar * otherVar)();
  }

  return result;
}

// Norm implementation
template <typename ConfigBackend, typename GeometryBackend>
double WRFState<ConfigBackend, GeometryBackend>::norm() const {
  double sumSquares = 0.0;

  // Sum squares of all variables
  for (const auto& [name, data] : variables_) {
    sumSquares += xt::sum(xt::square(data))();
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

  // Check if all variables have the same values
  for (const auto& varName : variableNames_) {
    const auto& thisVar = variables_.at(varName);
    const auto& otherVar = other.variables_.at(varName);

    // Check if arrays are equal within a small tolerance
    auto diff = xt::abs(thisVar - otherVar);
    double maxDiff = xt::amax(diff)();

    if (maxDiff > 1e-10) {
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

  // Add each variable
  for (const auto& varName : variableNames_) {
    auto& thisVar = variables_.at(varName);
    const auto& otherVar = other.variables_.at(varName);

    thisVar += otherVar;
  }
}

// Subtract implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::subtract(
    const WRFState<ConfigBackend, GeometryBackend>& other) {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for subtraction");
  }

  // Subtract each variable
  for (const auto& varName : variableNames_) {
    auto& thisVar = variables_.at(varName);
    const auto& otherVar = other.variables_.at(varName);

    thisVar -= otherVar;
  }
}

// Multiply implementation
template <typename ConfigBackend, typename GeometryBackend>
void WRFState<ConfigBackend, GeometryBackend>::multiply(double scalar) {
  // Multiply each variable by the scalar
  for (auto& [name, data] : variables_) {
    data *= scalar;
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
      variables_.clear();
      dimensions_.clear();
      variableNames_.clear();
      variable_grid_info_.clear();

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

          // Create xtensor array with proper shape
          xt::xarray<double> xdata = xt::adapt(data, dims);

          // Store the variable
          variables_[varName] = std::move(xdata);
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

}  // namespace metada::backends::wrf