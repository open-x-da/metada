/**
 * @file MACOMState.hpp
 * @brief MACOM state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <fstream>
#include <iostream>
#include <memory>
#include <netcdf>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "MACOMParallel.hpp"
#include "include/MACOMlogging.hpp"
// #include "MACOMGeometry.hpp"

// #include <xtensor/containers/xadapt.hpp>
// #include <xtensor/containers/xarray.hpp>

// namespace metada::backends::macom {
// template <typename ConfigBackend>
// class MACOMGeometry;
// }

namespace metada::backends::macom {

/**
 * @brief MACOM state backend implementation
 *
 * @details
 * This class implements a state backend for the MACOM model. It manages
 * oceanographic state variables and provides operations required by the State
 * adapter.
 */
template <typename ConfigBackend, typename GeometryBackend>
class MACOMState {
 public:
  // =============================================================================
  // FRAMEWORK CONCEPTS REQUIRED INTERFACES
  // Required by StateBackendImpl concept
  // =============================================================================

  // --- Resource management (required by framework) ---

  /**
   * @brief Default constructor is deleted (required by framework)
   */
  MACOMState() = delete;

  /**
   * @brief Copy constructor is deleted (required by framework)
   */
  MACOMState(const MACOMState&) = delete;

  /**
   * @brief Copy assignment operator is deleted (required by framework)
   */
  MACOMState& operator=(const MACOMState&) = delete;

  /**
   * @brief Constructor that takes a configuration backend (required by
   * framework)
   * @param config Configuration containing MACOM file path and variables
   * @param geometry Geometry object containing grid information
   */
  MACOMState(const ConfigBackend& config, const GeometryBackend& geometry);

  /**
   * @brief CopySkipInitialization constructor
   */
  MACOMState(const ConfigBackend& config, const GeometryBackend& geometry,
             bool CopySkipInitialization);

  /**
   * @brief Move constructor (required by framework)
   * @param other MACOM state backend to move from
   */
  MACOMState(MACOMState&& other) noexcept;

  /**
   * @brief Move assignment operator (required by framework)
   * @param other MACOM state backend to move from
   * @return Reference to this state after assignment
   */
  MACOMState& operator=(MACOMState&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~MACOMState() = default;

  // --- Data access interface (required by framework) ---

  /**
   * @brief Get mutable access to the underlying data (required by framework)
   * @return Void pointer to the active data array
   */
  void* getData();

  /**
   * @brief Get const access to the underlying data (required by framework)
   * @return Const void pointer to the active data array
   */
  const void* getData() const;

  // --- Variable management interface (required by framework) ---

  /**
   * @brief Get variable names (required by framework)
   * @return Const reference to vector of variable names
   */
  const std::vector<std::string>& getVariableNames() const;

  /**
   * @brief Get total size of state vector (required by framework)
   * @return Number of elements in the active variable
   */
  size_t size() const {
    return variables_.empty() ? 0 : variables_.at(activeVariable_).size();
  }

  // --- Index-based access (required by BackgroundErrorCovariance) ---

  /**
   * @brief Access element by index (required by BackgroundErrorCovariance)
   * @param index Index of the element to access
   * @return Reference to the element
   */
  double& operator[](size_t index) {
    if (variables_.empty() || activeVariable_.empty()) {
      throw std::runtime_error("No active variable set");
    }
    if (index >= variables_.at(activeVariable_).size()) {
      throw std::out_of_range("Index out of range");
    }
    return variables_.at(activeVariable_)[index];
  }

  /**
   * @brief Access element by index (const version)
   * @param index Index of the element to access
   * @return Const reference to the element
   */
  const double& operator[](size_t index) const {
    if (variables_.empty() || activeVariable_.empty()) {
      throw std::runtime_error("No active variable set");
    }
    if (index >= variables_.at(activeVariable_).size()) {
      throw std::out_of_range("Index out of range");
    }
    return variables_.at(activeVariable_)[index];
  }

  // --- Output operator (required by framework) ---

  /**
   * @brief Output operator for MACOMState
   * @param os Output stream
   * @param state MACOMState to output
   * @return Output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const MACOMState& state) {
    os << "MACOMState{";
    os << "initialized: " << (state.initialized_ ? "true" : "false");
    os << ", active_variable: \"" << state.activeVariable_ << "\"";
    os << ", variables: [";
    for (size_t i = 0; i < state.variableNames_.size(); ++i) {
      if (i > 0) os << ", ";
      os << "\"" << state.variableNames_[i] << "\"";
    }
    os << "]";
    os << ", size: " << state.size();
    os << ", nlpb: " << state.nlpb_;
    os << ", nk: " << state.nk_;
    os << "}";
    return os;
  }

  // --- Cloning interface (required by framework) ---

  /**
   * @brief Clone this state (required by framework)
   * @return Unique pointer to a new identical MACOM state backend
   */
  std::unique_ptr<MACOMState> clone() const;

  // --- Vector arithmetic interface (required by framework) ---

  /**
   * @brief Zero all elements (required by framework)
   */
  void zero();

  /**
   * @brief Add another state to this one (required by framework)
   * @param other State to add
   */
  void add(const MACOMState& other);

  /**
   * @brief Subtract another state from this one (required by framework)
   * @param other State to subtract
   */
  void subtract(const MACOMState& other);

  /**
   * @brief Multiply this state by a scalar (required by framework)
   * @param scalar Scalar value to multiply by
   */
  void multiply(double scalar);

  /**
   * @brief Add random perturbations to the state variables
   * @param perturbation_magnitude Standard deviation of Gaussian perturbations
   * @param seed Random seed (optional, uses current time if not provided)
   */
  void addPerturbation(double perturbation_magnitude, unsigned int seed = 0);

  /**
   * @brief Directly modify state variables at a specific location by a fixed
   * value
   * @param lat Latitude of the modification location (degrees)
   * @param lon Longitude of the modification location (degrees)
   * @param delta_value Value to add to the state variables (negative to
   * subtract)
   */
  void modifyValueAtLocation(double lat, double lon, double delta_value);

  /**
   * @brief Compute dot product with another state (required by framework)
   * @param other State to compute dot product with
   * @return Dot product value
   */
  double dot(const MACOMState& other) const;

  /**
   * @brief Compute norm of this state (required by framework)
   * @return Norm value
   */
  double norm() const;

  // --- Comparison interface (required by framework) ---

  /**
   * @brief Check equality with another state (required by framework)
   * @param other State to compare with
   * @return True if states are equal, false otherwise
   */
  bool equals(const MACOMState& other) const;

  // --- File I/O interface (required by framework) ---

  /**
   * @brief Save state to file (required by framework)
   * @param filename Output file path
   */
  void saveToFile(const std::string& filename) const {
    if (!initialized_) {
      throw std::runtime_error("MACOMState not initialized");
    }

    // Determine output format based on filename extension
    std::string ext = filename.substr(filename.find_last_of(".") + 1);

    if (ext == "nc") {
      saveToNetCDF(filename);
    } else {
      saveToText(filename);
    }
  }

  /**
   * @brief Set template file for NetCDF output
   * @param template_file Path to template NetCDF file
   */
  void setTemplateFile(const std::string& template_file) {
    template_file_ = template_file;
    MACOM_LOG_INFO("MACOMState", "Template file set to: " + template_file);
  }

  // =============================================================================
  // IDENTITY OBS OPERATOR COMPATIBILITY INTERFACES
  // Required for compatibility with IdentityObsOperator
  // =============================================================================

  /**
   * @brief Get access to the associated geometry (required by
   * IdentityObsOperator)
   * @return Const reference to the geometry backend
   */
  const GeometryBackend& geometry() const { return geometry_; }

  /**
   * @brief Access data at grid coordinates (required by IdentityObsOperator)
   *
   * @param coord Grid coordinates as a pair (i, j)
   * @return Reference to data at the specified coordinates
   */
  double& at(const std::pair<int, int>& coord) {
    if (!initialized_ || activeVariable_.empty()) {
      throw std::runtime_error(
          "MACOMState not initialized or no active variable set");
    }

    try {
      auto& var_data = variables_.at(activeVariable_);
      size_t x_dim = geometry_.getNlpb();
      size_t index = coord.second * x_dim + coord.first;

      if (index >= var_data.size()) {
        throw std::out_of_range("Coordinate access index out of bounds");
      }

      return var_data[index];
    } catch (const std::out_of_range& e) {
      throw std::runtime_error("Variable '" + activeVariable_ +
                               "' not found in state");
    }
  }

  /**
   * @brief Access data at grid coordinates (const version)
   *
   * @param coord Grid coordinates as a pair (i, j)
   * @return Const reference to data at the specified coordinates
   */
  const double& at(const std::pair<int, int>& coord) const {
    if (!initialized_ || activeVariable_.empty()) {
      throw std::runtime_error(
          "MACOMState not initialized or no active variable set");
    }

    try {
      const auto& var_data = variables_.at(activeVariable_);
      size_t x_dim = geometry_.getNlpb();
      size_t index = coord.second * x_dim + coord.first;

      if (index >= var_data.size()) {
        throw std::out_of_range("Coordinate access index out of bounds");
      }

      return var_data[index];
    } catch (const std::out_of_range& e) {
      throw std::runtime_error("Variable '" + activeVariable_ +
                               "' not found in state");
    }
  }

  /**
   * @brief Access data at geographic location (required by framework)
   *
   * @param location Geographic location
   * @return Reference to data at the nearest grid point
   */
  double& at(const framework::Location& location) {
    if (!initialized_ || activeVariable_.empty()) {
      throw std::runtime_error(
          "MACOMState not initialized or no active variable set");
    }

    // Convert geographic location to grid coordinates
    if (location.getCoordinateSystem() ==
        framework::CoordinateSystem::GEOGRAPHIC) {
      auto [lat, lon, level] = location.getGeographicCoords();

      // Find nearest grid point
      auto nearest_point = geometry_.findNearestGridPoint(lon, lat);

      // Access data at the nearest grid point
      try {
        auto& var_data = variables_.at(activeVariable_);
        if (nearest_point.index >= var_data.size()) {
          throw std::out_of_range("Grid point index out of bounds");
        }

        return var_data[nearest_point.index];
      } catch (const std::out_of_range& e) {
        throw std::runtime_error("Variable '" + activeVariable_ +
                                 "' not found in state");
      }
    } else {
      throw std::runtime_error(
          "MACOMState only supports geographic coordinate locations");
    }
  }

  /**
   * @brief Access data at geographic location (const version)
   *
   * @param location Geographic location
   * @return Const reference to data at the nearest grid point
   */
  const double& at(const framework::Location& location) const {
    if (!initialized_ || activeVariable_.empty()) {
      throw std::runtime_error(
          "MACOMState not initialized or no active variable set");
    }

    // Convert geographic location to grid coordinates
    if (location.getCoordinateSystem() ==
        framework::CoordinateSystem::GEOGRAPHIC) {
      auto [lat, lon, level] = location.getGeographicCoords();

      // Find nearest grid point
      auto nearest_point = geometry_.findNearestGridPoint(lon, lat);

      // Access data at the nearest grid point
      try {
        const auto& var_data = variables_.at(activeVariable_);
        if (nearest_point.index >= var_data.size()) {
          throw std::out_of_range("Grid point index out of bounds");
        }

        return var_data[nearest_point.index];
      } catch (const std::out_of_range& e) {
        throw std::runtime_error("Variable '" + activeVariable_ +
                                 "' not found in state");
      }
    } else {
      throw std::runtime_error(
          "MACOMState only supports geographic coordinate locations");
    }
  }

  // =============================================================================
  // MACOM SPECIFIC FUNCTIONALITY
  // These are MACOM-specific methods beyond framework requirements
  // =============================================================================

  // --- MACOM variable management ---

  /**
   * @brief Get active variable name
   * @return Const reference to active variable name
   */
  const std::string& getActiveVariable() const;

  /**
   * @brief Set active variable
   * @param name Variable name to set as active
   */
  void setActiveVariable(const std::string& name);

  /**
   * @brief Get dimensions for a specific variable
   * @param name Variable name
   * @return Const reference to vector of dimensions
   */
  const std::vector<size_t>& getDimensions(const std::string& name) const;

  // --- MACOM grid interface ---

  /**
   * @brief Set the associated geometry object
   * @param geometry Pointer to a MACOMGeometry object
   */
  void getDimensionsFromGeometry(const GeometryBackend* geometry);

  // --- MACOM data retrieval ---

  /**
   * @brief Get values at nearest grid points for multiple query locations
   * @param query_lons Vector of query longitudes
   * @param query_lats Vector of query latitudes
   * @param query_depths Vector of query depths
   * @param var_name Variable name to query
   * @param use_horizontal_interp Whether to use horizontal interpolation
   * @return Vector of interpolated values
   */
  std::vector<double> getValuesAtNearestPoints(
      const std::vector<double>& query_lons,
      const std::vector<double>& query_lats,
      const std::vector<double>& query_depths, const std::string& var_name,
      bool use_horizontal_interp = false) const;

  // --- MACOM mode checking ---

  /**
   * @brief Check if running in Fortran mode
   * @return True if in Fortran mode, false if in C++ mode
   */
  bool isFortranMode() const {
    return backends::macom::MACOMParallel::getInstance().isFortranMode();
  }

  // --- MACOM initialization ---

  /**
   * @brief Initialize the state with configuration and geometry
   * @param config Configuration containing MACOM file path and variables
   * @param geometry Geometry object containing grid information
   * @return bool True if initialization successful, false otherwise
   */
  bool CPPInitialization(const ConfigBackend& config,
                         const GeometryBackend& geometry);

  /**
   * @brief Initialize the state in Fortran mode
   * @return bool True if initialization successful, false otherwise
   */
  bool FortranInitialization();

 private:
  /**
   * @brief Load variable dimensions from NetCDF file
   *
   * @param ncFile NetCDF file handle
   */
  void loadVariableDimensions(netCDF::NcFile& ncFile);

  /**
   * @brief Load variable arrays from NetCDF file
   *
   * @param ncFile NetCDF file handle
   */
  void loadVariableArrays(netCDF::NcFile& ncFile,
                          const std::vector<std::string>& variables);

  /**
   * @brief Check compatibility with another state for operations
   * (same variables and dimensions)
   * @param other The other MACOMState to compare with
   * @return True if compatible, false otherwise
   */
  bool isCompatible(const MACOMState& other) const;

  /**
   * @brief Save state to NetCDF format based on template file
   * @param filename Output NetCDF file path
   */
  void saveToNetCDF(const std::string& filename) const;

  /**
   * @brief Save state to text format (backup option)
   * @param filename Output text file path
   */
  void saveToText(const std::string& filename) const;

  // Reference to configuration
  const ConfigBackend& config_;
  const GeometryBackend& geometry_;

  // Variable dimensions
  std::size_t nlpb_ = 0;  // Number of grid points
  std::size_t nk_ = 0;    // Number of vertical levels

  // Variable data
  std::vector<double> u;  // u-velocity
  std::vector<double> v;  // v-velocity
  std::vector<double> t;  // temperature
  std::vector<double> s;  // salinity
  std::vector<double> w;  // w-velocity

  // State information
  std::string inputFile_;      // Input data file path
  std::string timestamp_;      // Timestamp of the data
  std::string template_file_;  // Template file for NetCDF output
  std::unique_ptr<MACOMFortranInterface> fortranInterface_;
  bool initialized_ = false;
  bool CopySkipInitialization = true;

  // Data storage (this would be filled by the actual implementation)
  std::unordered_map<std::string, std::vector<double>> variables_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
  std::vector<std::string> variableNames_;  // Available variables
  std::string activeVariable_;              // Currently active variable

};  // namespace metada::backends::macom

// ConfigBackend constructor implementation and GeometryBackend constructor
template <typename ConfigBackend, typename GeometryBackend>
MACOMState<ConfigBackend, GeometryBackend>::MACOMState(
    const ConfigBackend& config, const GeometryBackend& geometry)
    : config_(config), geometry_(geometry) {
  if (isFortranMode()) {
    // MACOM_LOG_INFO("MACOMState", "Running in Fortran mode");
    if (!FortranInitialization()) {
      throw std::runtime_error(
          "Failed to initialize MACOM state in Fortran mode");
    }
  } else {
    MACOM_LOG_INFO("MACOMState", "Running in C++ mode");
    if (!CPPInitialization(config, geometry)) {
      throw std::runtime_error("Failed to initialize MACOM state");
    }
  }
}

// ConfigBackend constructor implementation and GeometryBackend constructor
template <typename ConfigBackend, typename GeometryBackend>
MACOMState<ConfigBackend, GeometryBackend>::MACOMState(
    const ConfigBackend& config, const GeometryBackend& geometry,
    bool CopySkipInitialization)
    : config_(config), geometry_(geometry) {
  if (!CopySkipInitialization) {
    CPPInitialization(config, geometry);
  }
}

// Constructor implementation with ConfigBackend
template <typename ConfigBackend, typename GeometryBackend>
MACOMState<ConfigBackend, GeometryBackend>::MACOMState(
    MACOMState<ConfigBackend, GeometryBackend>&& other) noexcept
    : config_(other.config_),
      geometry_(other.geometry_),
      inputFile_(std::move(other.inputFile_)),
      initialized_(other.initialized_),
      variables_(std::move(other.variables_)),
      dimensions_(std::move(other.dimensions_)),
      variableNames_(std::move(other.variableNames_)),
      activeVariable_(std::move(other.activeVariable_)),
      nlpb_(other.nlpb_),
      nk_(other.nk_) {
  other.initialized_ = false;
  other.variableNames_.clear();
  other.activeVariable_.clear();
}

// Move assignment operator implementation
template <typename ConfigBackend, typename GeometryBackend>
MACOMState<ConfigBackend, GeometryBackend>&
MACOMState<ConfigBackend, GeometryBackend>::operator=(
    MACOMState<ConfigBackend, GeometryBackend>&& other) noexcept {
  if (this != &other) {
    inputFile_ = std::move(other.inputFile_);
    initialized_ = other.initialized_;
    variables_ = std::move(other.variables_);
    dimensions_ = std::move(other.dimensions_);
    variableNames_ = std::move(other.variableNames_);
    activeVariable_ = std::move(other.activeVariable_);

    other.initialized_ = false;
    other.variableNames_.clear();
    other.activeVariable_.clear();
  }
  return *this;
}

template <typename ConfigBackend, typename GeometryBackend>
std::unique_ptr<MACOMState<ConfigBackend, GeometryBackend>>
MACOMState<ConfigBackend, GeometryBackend>::clone() const {
  auto cloned = std::make_unique<MACOMState<ConfigBackend, GeometryBackend>>(
      config_, geometry_, CopySkipInitialization);

  cloned->initialized_ = this->initialized_;
  cloned->inputFile_ = this->inputFile_;
  cloned->variables_ = this->variables_;
  cloned->activeVariable_ = this->activeVariable_;
  cloned->variableNames_ = this->variableNames_;
  cloned->dimensions_ = this->dimensions_;

  cloned->nlpb_ = this->nlpb_;
  cloned->nk_ = this->nk_;

  return cloned;
}

template <typename ConfigBackend, typename GeometryBackend>
void* MACOMState<ConfigBackend, GeometryBackend>::getData() {
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
const void* MACOMState<ConfigBackend, GeometryBackend>::getData() const {
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
const std::string&
MACOMState<ConfigBackend, GeometryBackend>::getActiveVariable() const {
  return activeVariable_;
}

// Set active variable implementation
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::setActiveVariable(
    const std::string& name) {
  if (variables_.find(name) == variables_.end()) {
    throw std::out_of_range("Variable not found: " + name);
  }
  activeVariable_ = name;
}

// Get variable names implementation
template <typename ConfigBackend, typename GeometryBackend>
const std::vector<std::string>&
MACOMState<ConfigBackend, GeometryBackend>::getVariableNames() const {
  return variableNames_;
}

// Get dimensions implementation
template <typename ConfigBackend, typename GeometryBackend>
const std::vector<size_t>&
MACOMState<ConfigBackend, GeometryBackend>::getDimensions(
    const std::string& name) const {
  try {
    // auto it = dimensions_.find(name);
    // if (it != dimensions_.end()) {
    //   return it->second;
    // }

    // if (geometry_ptr_ && geometry_ptr_->isInitialized()) {
    //   std::vector<size_t> dims;

    //   if (name == "u" || name == "v" || name == "t" || name == "s") {
    //     dims = {geometry_ptr_->nlpb_, geometry_ptr_->nk_};
    //   } else if (name == "w") {
    //     dims = {geometry_ptr_->nlpb_, geometry_ptr_->nkp1_};
    //   } else {
    //     throw std::out_of_range("Unknown variable: " + name);
    //   }

    //   dimensions_[name] = dims;
    //   return dimensions_[name];
    // }

    throw std::out_of_range("No dimension information available for: " + name);
    return dimensions_.at(name);
  } catch (const std::out_of_range&) {
    throw std::out_of_range("Variable not found: " + name);
  }
}
// Zero implementation
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::zero() {
  for (auto& [name, data] : variables_) {
    std::fill(data.begin(), data.end(), 0.0);
  }
}

// Dot product implementation
template <typename ConfigBackend, typename GeometryBackend>
double MACOMState<ConfigBackend, GeometryBackend>::dot(
    const MACOMState<ConfigBackend, GeometryBackend>& other) const {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for dot product");
  }

  double result = 0.0;

  // // Sum dot products of all variables
  // for (const auto& varName : variableNames_) {
  //   const auto& thisVar = variables_.at(varName);
  //   const auto& otherVar = other.variables_.at(varName);

  //   // Use xtensor to compute element-wise multiplication and sum
  //   result += xt::sum(thisVar * otherVar)();
  // }

  return result;
}

// Norm implementation
template <typename ConfigBackend, typename GeometryBackend>
double MACOMState<ConfigBackend, GeometryBackend>::norm() const {
  double sumSquares = 0.0;

  return std::sqrt(sumSquares);
}

// Equals implementation
template <typename ConfigBackend, typename GeometryBackend>
bool MACOMState<ConfigBackend, GeometryBackend>::equals(
    const MACOMState<ConfigBackend, GeometryBackend>& other) const {
  if (!isCompatible(other)) {
    return false;
  }

  // // Check if all variables have the same values
  // for (const auto& varName : variableNames_) {
  //   const auto& thisVar = variables_.at(varName);
  //   const auto& otherVar = other.variables_.at(varName);

  //   // Check if arrays are equal within a small tolerance
  //   auto diff = xt::abs(thisVar - otherVar);
  //   double maxDiff = xt::amax(diff)();

  //   if (maxDiff > 1e-10) {
  //     return false;
  //   }
  // }

  return true;
}

template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::add(const MACOMState& other) {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for addition");
  }

  for (const auto& [name, data] : variables_) {
    if (other.variables_.find(name) != other.variables_.end()) {
      const auto& other_data = other.variables_.at(name);
      if (data.size() == other_data.size()) {
        for (size_t i = 0; i < data.size(); ++i) {
          variables_[name][i] += other_data[i];
        }
      }
    }
  }
}

template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::subtract(
    const MACOMState<ConfigBackend, GeometryBackend>& other) {
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for subtraction");
  }

  for (const auto& [name, data] : variables_) {
    if (other.variables_.find(name) != other.variables_.end()) {
      const auto& other_data = other.variables_.at(name);
      if (data.size() == other_data.size()) {
        for (size_t i = 0; i < data.size(); ++i) {
          variables_[name][i] -= other_data[i];
        }
      }
    }
  }
}

template <typename ConfigBackend, typename GeometryBackend>
bool MACOMState<ConfigBackend, GeometryBackend>::isCompatible(
    const MACOMState& other) const {
  // Check if both states are initialized
  if (!initialized_ || !other.initialized_) {
    return false;
  }

  // Check if they have the same variables
  if (variables_.size() != other.variables_.size()) {
    return false;
  }

  // Check if all variables exist in both states and have the same dimensions
  for (const auto& [name, data] : variables_) {
    auto other_it = other.variables_.find(name);
    if (other_it == other.variables_.end()) {
      return false;  // Variable not found in other state
    }

    if (data.size() != other_it->second.size()) {
      return false;  // Different sizes for the same variable
    }
  }

  // Check if dimensions match
  if (nlpb_ != other.nlpb_ || nk_ != other.nk_) {
    return false;
  }

  return true;
}

template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::multiply(double scalar) {
  // Multiply each variable by the scalar
  for (auto& pair : variables_) {
    // Assuming variables_ stores xt::xarray<double> or similar that supports *=
    // If it's std::vector<double>, loop through elements:
    for (double& val : pair.second) {
      val *= scalar;
    }
  }
}

// Add perturbation implementation
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::addPerturbation(
    double perturbation_magnitude, unsigned int seed) {
  if (!initialized_) {
    throw std::runtime_error("MACOMState not initialized");
  }

  // Initialize random number generator
  std::random_device rd;
  std::mt19937 gen(seed == 0 ? rd() : seed);
  std::normal_distribution<double> dist(0.0, perturbation_magnitude);

  // Add perturbations to all variables
  for (auto& [name, data] : variables_) {
    for (auto& value : data) {
      value += dist(gen);
    }
  }

  std::cout << "Added Gaussian perturbations (std=" << perturbation_magnitude
            << ") to " << variables_.size() << " variables with seed=" << seed
            << std::endl;
}

// Implementation of loadVariableDimensions
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::loadVariableDimensions(
    netCDF::NcFile& ncFile) {
  auto readDimension = [&ncFile](const std::string& name, std::size_t& value) {
    // MACOM_LOG_INFO("MACOMState", "Attempting to read dimension: " + name);
    auto dim = ncFile.getDim(name);
    if (dim.isNull()) {
      // MACOM_LOG_ERROR("MACOMState",
      //                 "Dimension '" + name + "' not found in grid file");
      throw std::runtime_error("Dimension '" + name +
                               "' not found in grid file");
    }
    // MACOM_LOG_INFO("MACOMState", "Found dimension: " + name);
    value = dim.getSize();
    // MACOM_LOG_INFO("MACOMState",
    //                "Read size for " + name + ": " + std::to_string(value));
  };

  // Read all variable dimensions
  readDimension("nlpb", nlpb_);
  readDimension("nk", nk_);

  // MACOM_LOG_INFO("MACOMState", "Loaded variable dimensions:");
  // MACOM_LOG_INFO("MACOMState", "  nlpb=" + std::to_string(nlpb_) +
  //                                  ", nk=" + std::to_string(nk_));
}

// Implementation of loadVariableArrays
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::loadVariableArrays(
    netCDF::NcFile& ncFile, const std::vector<std::string>& variables_to_load) {
  // This part resizes individual u,v,t,s,w members - keeping as is for now,
  // but note that primary data storage is in variables_ map.
  // for (const auto& variable_name : variables_to_load) {
  //   if (variable_name == "u") {
  //     u.resize(nlpb_ * nk_);
  //   } else if (variable_name == "v") {
  //     v.resize(nlpb_ * nk_);
  //   } else if (variable_name == "t") {
  //     t.resize(nlpb_ * nk_);
  //   } else if (variable_name == "s") {
  //     s.resize(nlpb_ * nk_);
  //   } else if (variable_name == "w") {
  //     // 'w' might have different vertical dimension (e.g., nkp1_ from
  //     geometry)
  //     // For now, using nk_ consistent with other variables from file.
  //     w.resize(nlpb_ * nk_);
  //   }
  // }

  auto readVarLambda = [&ncFile](const std::string& name,
                                 std::vector<double>& data_vec) {
    auto var = ncFile.getVar(name);
    if (var.isNull()) {
      throw std::runtime_error("Variable '" + name +
                               "' not found in initial data file");
    }
    // Data vector should be resized before var.getVar() is called.
    // The size comes from this->nlpb_ and this->nk_ read from the file.
    // data_vec.resize(this->nlpb_ * this->nk_); // This is now done before
    // calling lambda
    var.getVar(data_vec.data());
  };

  // // Clear existing data in variables_ map and dimensions_ map to avoid
  // issues
  // if this is a reload
  variables_.clear();
  dimensions_.clear();

  // Read grid data
  for (const auto& variable_name : variables_to_load) {
    // Ensure vector exists in map and is correctly sized
    // Using file's global nlpb_ and nk_ for all variables.
    // A more robust solution would query each NetCDF variable for its specific
    // dimensions.
    variables_[variable_name].resize(this->nlpb_ * this->nk_);
    readVarLambda(variable_name, variables_[variable_name]);
    dimensions_[variable_name] = {
        this->nlpb_, this->nk_};  // Store dimensions for this variable

    // // Output sample data: take 10 values from the middle of variable data
    // std::stringstream sample_ss;
    // sample_ss << "Sample values for " << variable_name << ": ";

    // size_t total_size = variables_[variable_name].size();
    // size_t mid_point = total_size / 2;
    // size_t start_idx = (total_size <= 10) ? 0 : (mid_point - 5);
    // size_t end_idx = std::min(start_idx + 10, total_size);

    // for (size_t i = start_idx; i < end_idx; ++i) {
    //   sample_ss << variables_[variable_name][i];
    //   if (i < end_idx - 1) sample_ss << ", ";
    // }

    // MACOM_LOG_INFO("MACOMState", sample_ss.str());

    // if (!variables_[variable_name].empty()) {
    //   double min_val = variables_[variable_name][0];
    //   double max_val = variables_[variable_name][0];
    //   double sum = 0.0;

    //   for (const auto& val : variables_[variable_name]) {
    //     min_val = std::min(min_val, val);
    //     max_val = std::max(max_val, val);
    //     sum += val;
    //   }

    //   double avg = sum / total_size;
    //   std::stringstream stats_ss;
    //   stats_ss << "Statistics for " << variable_name << ": min=" << min_val
    //            << ", max=" << max_val << ", avg=" << avg
    //            << ", count=" << total_size;

    //   MACOM_LOG_INFO("MACOMState", stats_ss.str());
    // }
  }

  // // Replace original output statements
  // if (!variables_to_load.empty()) {
  //   MACOM_LOG_INFO("MACOMState", "Successfully loaded " +
  //                                    std::to_string(variables_to_load.size())
  //                                    + " variable arrays");
  // } else {
  //   MACOM_LOG_WARNING("MACOMState", "No variables specified to load");
  // }
  // Note: this->variableNames_ and this->activeVariable_ are set in the
  // constructor after this function successfully completes.
}

// Implementation of loadVariableData - THIS METHOD IS NOW REMOVED
// template <typename ConfigBackend, typename GeometryBackend>
// void MACOMState<ConfigBackend, GeometryBackend>::loadVariableData(
// ...

template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::getDimensionsFromGeometry(
    const GeometryBackend* geometry) {
  if (!geometry) {
    // MACOM_LOG_ERROR("MACOMState",
    //                 "getDimensionsFromGeometry: null pointer passed");
    throw std::runtime_error("getDimensionsFromGeometry: null pointer passed");
  }
  if (!geometry->isInitialized()) {
    // MACOM_LOG_ERROR("MACOMState",
    //                 "getDimensionsFromGeometry: geometry not initialized");
    throw std::runtime_error(
        "getDimensionsFromGeometry: geometry not initialized");
  }
  // now pull the numbers out of the geometry
  nlpb_ = geometry->getNlpb();
  nk_ = geometry->getNk();
  // (if you need nkp1 you can grab that too)
}

template <typename ConfigBackend, typename GeometryBackend>
std::vector<double>
MACOMState<ConfigBackend, GeometryBackend>::getValuesAtNearestPoints(
    const std::vector<double>& query_lons,
    const std::vector<double>& query_lats,
    const std::vector<double>& query_depths, const std::string& var_name,
    bool use_horizontal_interp) const {
  if (!initialized_) {
    throw std::runtime_error(
        "MACOMState is not initialized. Cannot get values.");
  }
  if (!geometry_.isInitialized()) {
    throw std::runtime_error(
        "MACOMGeometry is not initialized. Cannot get values via geometry.");
  }
  if (query_lons.size() != query_lats.size() ||
      query_lons.size() != query_depths.size()) {
    throw std::invalid_argument(
        "Longitude, latitude and depth arrays must have the same size for "
        "batch get.");
  }

  auto var_it = variables_.find(var_name);
  if (var_it == variables_.end()) {
    throw std::runtime_error("Variable '" + var_name +
                             "' not found in MACOMState.");
  }
  const std::vector<double>& var_data = var_it->second;

  std::vector<metada::backends::macom::GeoPoint> nearest_points;
  if (var_name == "u") {
    nearest_points =
        geometry_.findNearestGridPointsBatchW(query_lons, query_lats);
  } else if (var_name == "v") {
    nearest_points =
        geometry_.findNearestGridPointsBatchS(query_lons, query_lats);
  } else {
    nearest_points =
        geometry_.findNearestGridPointsBatch(query_lons, query_lats);
  }

  // // Print grid point information
  // MACOM_LOG_INFO("MACOMState",
  //                "Found nearest grid points for " + var_name + ":");
  // for (size_t i = 0; i < nearest_points.size(); ++i) {
  //   MACOM_LOG_INFO(
  //       "MACOMState",
  //       "Point " + std::to_string(i) +
  //           ": index=" + std::to_string(nearest_points[i].index) +
  //           ", lon=" + std::to_string(nearest_points[i].lon) +
  //           ", lat=" + std::to_string(nearest_points[i].lat) +
  //           ", distance=" + std::to_string(nearest_points[i].distance) + "
  //           km");
  // }

  // Get nearest vertical points
  std::vector<metada::backends::macom::VerticalPoint> vertical_points =
      geometry_.findNearestVerticalPointsBatch(query_depths);

  std::vector<double> result_values(query_lons.size());

  for (size_t i = 0; i < nearest_points.size(); ++i) {
    size_t center_idx = nearest_points[i].index;
    const auto& vp = vertical_points[i];

    if (vp.is_outside) {
      result_values[i] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    // Horizontal interpolation
    if (use_horizontal_interp) {
      // Get neighbor indices (assume tw_, te_, tn_, ts_ are available and same
      // size as grid)
      size_t idx_west = geometry_.getTw()[center_idx];
      size_t idx_east = geometry_.getTe()[center_idx];
      size_t idx_north = geometry_.getTn()[center_idx];
      size_t idx_south = geometry_.getTs()[center_idx];

      // Collect indices
      std::vector<size_t> idxs = {center_idx, idx_west, idx_east, idx_north,
                                  idx_south};

      // Compute horizontal distances (haversine or Euclidean)
      std::vector<double> dists(5);
      for (size_t j = 0; j < 5; ++j) {
        // Get lon/lat for each neighbor
        double lon_j, lat_j;
        if (var_name == "u") {
          lon_j = geometry_.getLonW()[idxs[j]];
          lat_j = geometry_.getLatW()[idxs[j]];
        } else if (var_name == "v") {
          lon_j = geometry_.getLonS()[idxs[j]];
          lat_j = geometry_.getLatS()[idxs[j]];
        } else {
          lon_j = geometry_.getLonC()[idxs[j]];
          lat_j = geometry_.getLatC()[idxs[j]];
        }

        dists[j] = metada::backends::macom::MACOMGeometryIterator<
            ConfigBackend>::haversineDistance(query_lons[i], query_lats[i],
                                              lon_j, lat_j);
        // Avoid zero distance (for center point)
        if (dists[j] < 1e-8) dists[j] = 1e-8;
      }

      // Compute inverse distance weights
      double wsum = 0.0;
      std::vector<double> weights(5);
      for (size_t j = 0; j < 5; ++j) {
        weights[j] = 1.0 / dists[j] * dists[j];
        wsum += weights[j];
      }
      for (size_t j = 0; j < 5; ++j) {
        weights[j] /= wsum;
      }

      // Vertical interpolation for each neighbor
      double interp_val = 0.0;
      bool all_land = true;  // Flag to check if all points are land
      for (size_t j = 0; j < 5; ++j) {
        size_t lower_index = idxs[j] * this->nk_ + vp.lower_index;
        size_t upper_index = idxs[j] * this->nk_ + vp.upper_index;
        if (lower_index >= var_data.size() || upper_index >= var_data.size()) {
          throw std::out_of_range("Batch: Calculated indices out of bounds.");
        }
        double lower_value = var_data[lower_index];
        double upper_value = var_data[upper_index];

        // Check if both levels are zero (land)
        if (std::abs(lower_value) > 1e-10 || std::abs(upper_value) > 1e-10) {
          all_land = false;
        }

        double v_interp =
            lower_value + vp.interp_coef * (upper_value - lower_value);
        interp_val += weights[j] * v_interp;
      }

      // If all points are land, set to -99999
      result_values[i] = all_land ? -99999.0 : interp_val;
    }
    // Nearest point method (default)
    else {
      size_t lower_index = center_idx * this->nk_ + vp.lower_index;
      size_t upper_index = center_idx * this->nk_ + vp.upper_index;
      if (lower_index >= var_data.size() || upper_index >= var_data.size()) {
        throw std::out_of_range("Batch: Calculated indices out of bounds.");
      }
      double lower_value = var_data[lower_index];
      double upper_value = var_data[upper_index];

      // Check if both levels are zero (land)
      if (std::abs(lower_value) < 1e-10 && std::abs(upper_value) < 1e-10) {
        result_values[i] = -99999.0;  // Land point
      } else {
        result_values[i] =
            lower_value + vp.interp_coef * (upper_value - lower_value);
      }
    }
  }
  return result_values;
}

template <typename ConfigBackend, typename GeometryBackend>
bool MACOMState<ConfigBackend, GeometryBackend>::CPPInitialization(
    const ConfigBackend& config, const GeometryBackend& geometry) {
  try {
    MACOM_LOG_INFO("MACOMState",
                   "Initializing state with CPPInitialization...");

    // Get input file path from config (support both "file" and "input_file" for
    // backward compatibility)
    if (config.HasKey("file")) {
      inputFile_ = config.Get("file").asString();
    } else if (config.HasKey("input_file")) {
      inputFile_ = config.Get("input_file").asString();
    } else {
      MACOM_LOG_ERROR("MACOMState",
                      "Input file path not specified in configuration (neither "
                      "'file' nor 'input_file' found)");
      return false;
    }

    if (inputFile_.empty()) {
      MACOM_LOG_ERROR("MACOMState",
                      "Input file path is empty in configuration");
      return false;
    }

    // // Get timestamp from file if available
    // timestamp_ = config.Get("timestamp").asString();

    // Get dimensions from geometry
    getDimensionsFromGeometry(&geometry);

    // Open and validate NetCDF file
    netCDF::NcFile ncFile(inputFile_, netCDF::NcFile::read);
    if (ncFile.isNull()) {
      MACOM_LOG_ERROR("MACOMState",
                      "Failed to open NetCDF file: " + inputFile_);
      return false;
    }

    // // Step 1: Load grid dimensions from the state file
    // loadVariableDimensions(ncFile);  // Sets this->nlpb_ and this->nk_

    // Step 2: Compare dimensions with geometry (if needed)
    // MACOM_LOG_INFO("MACOMState", "Dimensions from State File: nlpb_ = " +
    //                                  std::to_string(this->nlpb_) +
    //                                  ", nk_ = " + std::to_string(this->nk_));

    // Step 3: Load variable arrays
    // Get variables to load from config
    std::vector<std::string> variables;
    variables = config.Get("variables").asVectorString();
    loadVariableArrays(ncFile, variables);

    // Populate variableNames_ and activeVariable_
    this->variableNames_ = variables;
    if (!this->variableNames_.empty()) {
      this->activeVariable_ = this->variableNames_[0];
      // MACOM_LOG_INFO("MACOMState",
      //                "Active variable set to: " + this->activeVariable_);
    } else {
      MACOM_LOG_WARNING("MACOMState", "No variables loaded or specified.");
    }

    initialized_ = true;

    // --- Test code for getValuesAtNearestPoints ---
    // Check if testing is enabled via configuration
    bool enable_testing = false;
    try {
    } catch (const std::exception& e) {
      // If config key doesn't exist, default to false (no testing)
      enable_testing = false;
    }

    if (enable_testing) {
      try {
        MACOM_LOG_INFO(
            "MACOMState",
            "Testing getValuesAtNearestPoints (enabled via config)...");

        // Test points
        std::vector<double> test_lons = {160.0, 161.0, 162.0};
        std::vector<double> test_lats = {36.0, 37.0, 38.0};
        std::vector<double> test_depths = {100.0, 200.0, 300.0};

        // Test each variable
        for (const auto& var_name : variables) {
          if (variables_.find(var_name) != variables_.end()) {
            auto values = getValuesAtNearestPoints(test_lons, test_lats,
                                                   test_depths, var_name, true);

            MACOM_LOG_INFO("MACOMState",
                           "Test results for variable " + var_name + ":");
            for (size_t i = 0; i < values.size(); ++i) {
              MACOM_LOG_INFO("MACOMState",
                             "  Point " + std::to_string(i) + " (" +
                                 std::to_string(test_lons[i]) + ", " +
                                 std::to_string(test_lats[i]) + ", " +
                                 std::to_string(test_depths[i]) +
                                 "): " + std::to_string(values[i]));
            }
          }
        }
      } catch (const std::exception& e) {
        MACOM_LOG_WARNING("MACOMState",
                          "Test code for getValuesAtNearestPoints failed: " +
                              std::string(e.what()));
        // Note: Don't set initialized_ = false here as test failure
        // shouldn't prevent normal operation
      }
    } else {
      MACOM_LOG_INFO("MACOMState",
                     "getValuesAtNearestPoints testing is disabled "
                     "(enable_getvalues_test=false)");
    }

    MACOM_LOG_INFO("MACOMState", "State initialization completed successfully");
    return true;

  } catch (const netCDF::exceptions::NcException& e) {
    MACOM_LOG_ERROR("MACOMState", "NetCDF error during initialization: " +
                                      std::string(e.what()));
    return false;
  } catch (const std::exception& e) {
    MACOM_LOG_ERROR("MACOMState",
                    "Error during initialization: " + std::string(e.what()));
    return false;
  }
}

template <typename ConfigBackend, typename GeometryBackend>
bool MACOMState<ConfigBackend, GeometryBackend>::FortranInitialization() {
  // MACOM_LOG_INFO("MACOMState", "Initializing state in Fortran mode...");

  auto& parallel = MACOMParallel::getInstance();
  if (!parallel.isFortranMode()) {
    MACOM_LOG_ERROR("MACOMState", "Not in Fortran mode");
    return false;
  }

  try {
    // Get the initialized Fortran interface from MACOMParallel
    // fortranInterface_ = std::make_unique<MACOMFortranInterface>();
    auto& parallelInterface = parallel.getFortranInterface();
    if (!parallelInterface) {
      MACOM_LOG_ERROR("MACOMState", "Fortran interface not available");
      return false;
    }

    // Share the Fortran interface
    fortranInterface_ = std::move(parallelInterface);

    // Only compute processes should initialize mitice
    if (fortranInterface_->getRank() < fortranInterface_->getCompProcs()) {
      // // Get config flags from Fortran
      // bool mitice_on, restart_in, assim_in;
      // int init_iter, max_iter;
      // fortranInterface_->getConfigFlags(mitice_on, restart_in, assim_in,
      // init_iter, max_iter);

      // if (mitice_on) {
      //   MACOM_LOG_INFO("MACOMState", "Initializing mitice components...");
      //   fortranInterface_->initializeMitice();
      // }

      fortranInterface_->initializeMitice();
      fortranInterface_->sendInfoToIO();
      fortranInterface_->openMiscRunInfo();
      fortranInterface_->initCSP();
      fortranInterface_->miticeInitAll();
      fortranInterface_->restartAndAssim();
      fortranInterface_->runCspStep();
    } else {
      fortranInterface_->CspIoMain();
      // MACOM_LOG_INFO(
      //     "MACOMState",
      //     std::string(
      //         "Skipping mitice initialization for I/O process, rank = ") +
      //         std::to_string(fortranInterface_->getRank()));
    }

    fortranInterface_->finalizeMitice();

    initialized_ = true;
    // MACOM_LOG_INFO("MACOMState", "Fortran initialization completed
    // successfully");
    return true;
  } catch (const std::exception& e) {
    MACOM_LOG_ERROR("MACOMState", "Error during Fortran initialization: " +
                                      std::string(e.what()));
    return false;
  }
}

// Implementation of saveToText method
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::saveToText(
    const std::string& filename) const {
  MACOM_LOG_INFO("MACOMState", "Saving state to text file: " + filename);

  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file for writing: " + filename);
  }

  // Write header information
  file << "# MACOM State File (Text Format)\n";
  file << "# Input file: " << inputFile_ << "\n";
  file << "# Active variable: " << activeVariable_ << "\n";
  file << "# Grid dimensions: nlpb=" << nlpb_ << ", nk=" << nk_ << "\n";
  file << "# Variables: ";
  for (const auto& var : variableNames_) {
    file << var << " ";
  }
  file << "\n\n";

  // Write variable data
  for (const auto& [name, data] : variables_) {
    file << "# Variable: " << name << " (size=" << data.size() << ")\n";
    for (size_t i = 0; i < data.size(); ++i) {
      file << data[i];
      if (i < data.size() - 1) file << " ";
      // Add line breaks every 10 values for readability
      if ((i + 1) % 10 == 0) file << "\n";
    }
    file << "\n\n";
  }

  file.close();
  MACOM_LOG_INFO("MACOMState", "Successfully saved state to: " + filename);
}

// Implementation of saveToNetCDF method
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::saveToNetCDF(
    const std::string& filename) const {
  MACOM_LOG_INFO("MACOMState", "Saving state to NetCDF file: " + filename);

  // Use template file if available, otherwise use input file
  std::string source_file =
      template_file_.empty() ? inputFile_ : template_file_;

  if (source_file.empty()) {
    throw std::runtime_error(
        "No template or input file available for NetCDF output");
  }

  MACOM_LOG_INFO("MACOMState", "Using template file: " + source_file);

  // Copy template file structure using C++ file operations
  std::ifstream source(source_file, std::ios::binary);
  if (!source.is_open()) {
    throw std::runtime_error("Failed to open source file: " + source_file);
  }

  std::ofstream dest(filename, std::ios::binary);
  if (!dest.is_open()) {
    throw std::runtime_error("Failed to create destination file: " + filename);
  }

  dest << source.rdbuf();
  source.close();
  dest.close();

  // Open the copied file for writing
  netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
  if (ncFile.isNull()) {
    throw std::runtime_error("Failed to open NetCDF file for writing: " +
                             filename);
  }

  // Update variable data
  for (const auto& [varName, varData] : variables_) {
    auto ncVar = ncFile.getVar(varName);
    if (!ncVar.isNull() && varData.size() > 0) {
      ncVar.putVar(varData.data());
      MACOM_LOG_INFO("MACOMState", "Updated variable: " + varName);
    } else {
      MACOM_LOG_WARNING("MACOMState",
                        "Variable not found in NetCDF file: " + varName);
    }
  }

  ncFile.close();
  MACOM_LOG_INFO("MACOMState", "Successfully saved NetCDF state");
}

// Simple value modification - now modifies ALL grid points for testing
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::modifyValueAtLocation(
    double lat, double lon, double delta_value) {
  if (!initialized_) {
    throw std::runtime_error("MACOMState not initialized");
  }

  // Modify values for the active variable only
  if (activeVariable_.empty()) {
    throw std::runtime_error("No active variable set");
  }

  auto it = variables_.find(activeVariable_);
  if (it == variables_.end()) {
    throw std::runtime_error("Active variable '" + activeVariable_ +
                             "' not found");
  }

  auto& data = it->second;
  size_t modified_count = 0;

  // FOR TESTING: Modify ALL grid points by the delta value
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] += delta_value;
    modified_count++;
  }

  std::string message =
      "TESTING MODE: Modified ALL " + std::to_string(modified_count) +
      " grid points of variable '" + activeVariable_ + "' by " +
      std::to_string(delta_value) + " (originally requested for location " +
      std::to_string(lat) + ", " + std::to_string(lon) + ")";
  std::cout << message << std::endl;
}

}  // namespace metada::backends::macom