/**
 * @file MACOMState.hpp
 * @brief MACOM state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <iostream>
#include <memory>
#include <netcdf>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "MACOMParallel.hpp"
#include "include/MACOMlogging.hpp"
// #include "MACOMGeometry.hpp"

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>

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
  /**
   * @brief Set the associated geometry object
   *
   * @param geometry Pointer to a MACOMGeometry object
   */
  void getDimensionsFromGeometry(const GeometryBackend* geometry);

  /**
   * @brief Default constructor is deleted
   */
  MACOMState() = delete;

  /**
   * @brief Copy constructor is deleted
   */
  MACOMState(const MACOMState&) = delete;

  /**
   * @brief Copy assignment operator is deleted
   */
  MACOMState& operator=(const MACOMState&) = delete;

  /**
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing MACOM file path and variables
   */
  MACOMState(const ConfigBackend& config, const GeometryBackend& geometry);

  /**
   * @brief CopySkipInitialization
   */
  MACOMState(const ConfigBackend& config, const GeometryBackend& geometry,
             bool CopySkipInitialization);

  /**
   * @brief Move constructor
   *
   * @param other MACOM state backend to move from
   */
  MACOMState(MACOMState&& other) noexcept;

  /**
   * @brief Move assignment operator
   *
   * @param other MACOM state backend to move from
   * @return Reference to this state after assignment
   */
  MACOMState& operator=(MACOMState&& other) noexcept;

  /**
   * @brief Destructor
   */
  ~MACOMState() = default;

  /**
   * @brief Clone this state
   *
   * @return Unique pointer to a new identical MACOM state backend
   */
  std::unique_ptr<MACOMState> clone() const;

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
   * @brief Get the dimensions of a variable
   *
   * @param name Name of the variable
   * @return Vector of dimension sizes
   * @throws std::out_of_range If variable doesn't exist
   */
  const std::vector<size_t>& getDimensions(const std::string& name) const;

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
  double dot(const MACOMState& other) const;

  /**
   * @brief Calculate norm of the state
   *
   * @return Norm of the state
   */
  double norm() const;

  /**
   * @brief Check if this state equals another
   *
   * @param other State to compare with
   * @return True if states are equal, false otherwise
   */
  bool equals(const MACOMState& other) const;
  /**
   * @brief Add another state to this one
   *
   * @param other State to add
   * @throws std::runtime_error If states are incompatible
   */
  void add(const MACOMState& other);

  /**
   * @brief Subtract another state from this one
   *
   * @param other State to subtract
   * @throws std::runtime_error If states are incompatible
   */
  void subtract(const MACOMState& other);

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
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Get the values of a variable at the nearest grid points for multiple
   * locations (first vertical layer).
   *
   * @param query_lons Vector of longitudes (degrees)
   * @param query_lats Vector of latitudes (degrees)
   * @param query_depths Vector of depths (meters)
   * @param var_name Name of the variable
   * @param use_horizontal_interp Boolean to use horizontal interpolation
   * @return Vector of variable values at the nearest grid points (first layer)
   * @throws std::runtime_error if not initialized or variable/points not found
   */
  std::vector<double> getValuesAtNearestPoints(
      const std::vector<double>& query_lons,
      const std::vector<double>& query_lats,
      const std::vector<double>& query_depths, const std::string& var_name,
      bool use_horizontal_interp = false) const;

  bool isFortranMode() const {
    return backends::macom::MACOMParallel::getInstance().isFortranMode();
  }

  /**
   * @brief Initialize the state with configuration and geometry
   *
   * @param config Configuration containing MACOM file path and variables
   * @param geometry Geometry object containing grid information
   * @return bool True if initialization successful, false otherwise
   */
  bool CPPInitialization(const ConfigBackend& config,
                         const GeometryBackend& geometry);

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
  // void loadVariableData(const std::string& filename,  // Removed this method
  //                       const std::vector<std::string>& variables);

  bool isCompatible(const MACOMState& other) const;

  // Reference to configuration
  const ConfigBackend& config_;
  const GeometryBackend& geometry_;

  // Variable dimensions
  std::size_t nlpb_ = 0;  // Number of grid points
  std::size_t nk_ = 0;    // Number of vertical levels

  std::size_t nlpb_grid = 0;  // Number of grid points
  std::size_t nk_grid = 0;    // Number of vertical levels

  // Variable data
  std::vector<double> u;  // u-velocity
  std::vector<double> v;  // v-velocity
  std::vector<double> t;  // temperature
  std::vector<double> s;  // salinity
  std::vector<double> w;  // w-velocity

  // State information
  std::string inputFile_;  // Input data file path
  std::string timestamp_;  // Timestamp of the data
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
    MACOM_LOG_INFO("MACOMState", "Running in Fortran mode");
    initialized_ = true;
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
      activeVariable_(std::move(other.activeVariable_)) {
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
  // cloned->variables_ = this->variables_;
  cloned->activeVariable_ = this->activeVariable_;
  cloned->variableNames_ = this->variableNames_;
  cloned->dimensions_ = this->dimensions_;

  // cloned->nlpb_ = this->nlpb_;
  // cloned->nk_ = this->nk_;
  // cloned->nlpb_grid = this->nlpb_grid;
  // cloned->nk_grid = this->nk_grid;
  // cloned->u = this->u;
  // cloned->v = this->v;
  // cloned->t = this->t;
  // cloned->s = this->s;
  // cloned->w = this->w;

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
  // for (auto& [name, data] : variables_) {
  //   data.fill(0.0);
  // }
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
void MACOMState<ConfigBackend, GeometryBackend>::add(
    [[maybe_unused]] const MACOMState& other) {}

template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::subtract(
    const MACOMState<ConfigBackend, GeometryBackend>& other) {
  // Subtract each variable
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for subtraction");
  }
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
  for (const auto& variable_name : variables_to_load) {
    if (variable_name == "u") {
      u.resize(nlpb_ * nk_);
    } else if (variable_name == "v") {
      v.resize(nlpb_ * nk_);
    } else if (variable_name == "t") {
      t.resize(nlpb_ * nk_);
    } else if (variable_name == "s") {
      s.resize(nlpb_ * nk_);
    } else if (variable_name == "w") {
      // 'w' might have different vertical dimension (e.g., nkp1_ from geometry)
      // For now, using nk_ consistent with other variables from file.
      w.resize(nlpb_ * nk_);
    }
  }

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

  // Clear existing data in variables_ map and dimensions_ map to avoid issues
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

    // // 输出样本数据：从变量数据中间取10个值
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

  // // 替换原来的输出语句
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
  nlpb_grid = geometry->getNlpb();
  nk_grid = geometry->getNk();
  // (if you need nkp1 you can grab that too)
  // MACOM_LOG_INFO("MACOMState",
  //                "State bound to geometry: nlpb=" + std::to_string(nlpb_grid)
  //                +
  //                    " nk=" + std::to_string(nk_grid));
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

  // Print grid point information
  MACOM_LOG_INFO("MACOMState",
                 "Found nearest grid points for " + var_name + ":");
  for (size_t i = 0; i < nearest_points.size(); ++i) {
    MACOM_LOG_INFO(
        "MACOMState",
        "Point " + std::to_string(i) +
            ": index=" + std::to_string(nearest_points[i].index) +
            ", lon=" + std::to_string(nearest_points[i].lon) +
            ", lat=" + std::to_string(nearest_points[i].lat) +
            ", distance=" + std::to_string(nearest_points[i].distance) + " km");
  }

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

    // Get input file path from config
    inputFile_ = config.Get("input_file").asString();
    if (inputFile_.empty()) {
      MACOM_LOG_ERROR("MACOMState",
                      "Input file path not specified in configuration");
      return false;
    }

    // Get timestamp from file if available
    timestamp_ = config.Get("timestamp").asString();

    // Get dimensions from geometry
    getDimensionsFromGeometry(&geometry);

    // MACOM_LOG_INFO("MACOMState", "Dimensions from Geometry: nlpb_grid = " +
    //                          std::to_string(nlpb_grid) +
    //                          ", nk_grid = " + std::to_string(nk_grid));

    // Open and validate NetCDF file
    netCDF::NcFile ncFile(inputFile_, netCDF::NcFile::read);
    if (ncFile.isNull()) {
      MACOM_LOG_ERROR("MACOMState",
                      "Failed to open NetCDF file: " + inputFile_);
      return false;
    }

    // // Step 1: Load grid dimensions from the state file
    // loadVariableDimensions(ncFile);  // Sets this->nlpb_ and this->nk_

    // Step 2: Compare dimensions with geometry
    // MACOM_LOG_INFO("MACOMState", "Dimensions from State File: nlpb_ = " +
    //                                  std::to_string(this->nlpb_) +
    //                                  ", nk_ = " + std::to_string(this->nk_));
    // if (this->nlpb_ != this->nlpb_grid || this->nk_ != this->nk_grid) {
    //   std::ostringstream error_msg;
    //   error_msg
    //       << "Dimension mismatch between MACOMState file and
    //       MACOMGeometry:\n"
    //       << "  State File dimensions: nlpb=" << this->nlpb_
    //       << ", nk=" << this->nk_ << "\n"
    //       << "  Geometry dimensions:   nlpb=" << this->nlpb_grid
    //       << ", nk=" << this->nk_grid;
    //   MACOM_LOG_ERROR("MACOMState", error_msg.str());
    //   throw std::runtime_error(error_msg.str());
    // }
    // MACOM_LOG_INFO("MACOMState",
    //                "State file dimensions match geometry dimensions.");
    nlpb_ = nlpb_grid;
    nk_ = nk_grid;

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

    // --- Test code for getValuesAtNearestPoints ---
    try {
      MACOM_LOG_INFO("MACOMState", "Testing getValuesAtNearestPoints...");

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
    }

    initialized_ = true;
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

}  // namespace metada::backends::macom