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
   * @brief Initialize variable data from NetCDF file
   *
   * @param filename Path to the NetCDF variable data file
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
    : config_(config), geometry_(geometry), initialized_(false) {
  inputFile_ = config.Get("input_file").asString();
  if (inputFile_.empty()) {
    // MACOM_LOG_ERROR(
    //     "MACOMState",
    //     "MACOM state input file path not specified in configuration");
    throw std::runtime_error(
        "MACOM state input file path not specified in configuration");
  }

  timestamp_ = config.Get("timestamp").asString();

  // Get variables to load from config
  std::vector<std::string> vars_to_load_from_config;
  try {
    vars_to_load_from_config = config.Get("variables").asVectorString();
  } catch (const std::exception&) {
    // MACOM_LOG_WARNING("MACOMState",
    //                   "No variables specified in configuration, attempting to
    //                   " "load all from file.");
  }

  getDimensionsFromGeometry(&geometry_);  // Sets nlpb_grid, nk_grid

  // MACOM_LOG_INFO("MACOMState", "Dimensions from Geometry: nlpb_grid = " +
  //                                  std::to_string(nlpb_grid) +
  //                                  ", nk_grid = " + std::to_string(nk_grid));

  try {
    // MACOM_LOG_INFO(
    //     "MACOMState",
    //     "Attempting to open NetCDF file for state data: " + inputFile_);
    netCDF::NcFile ncFile(inputFile_, netCDF::NcFile::read);
    if (ncFile.isNull()) {
      throw std::runtime_error("Failed to open NetCDF state file: " +
                               inputFile_);
    }
    // MACOM_LOG_INFO("MACOMState", "Successfully opened NetCDF state file.");

    // Step 1: Load grid dimensions from the state file
    loadVariableDimensions(ncFile);  // Sets this->nlpb_ and this->nk_

    // Step 2: Compare dimensions with geometry
    // MACOM_LOG_INFO("MACOMState", "Dimensions from State File: nlpb_ = " +
    //                                  std::to_string(this->nlpb_) +
    //                                  ", nk_ = " + std::to_string(this->nk_));
    if (this->nlpb_ != this->nlpb_grid || this->nk_ != this->nk_grid) {
      std::ostringstream error_msg;
      error_msg
          << "Dimension mismatch between MACOMState file and MACOMGeometry:\n"
          << "  State File dimensions: nlpb=" << this->nlpb_
          << ", nk=" << this->nk_ << "\n"
          << "  Geometry dimensions:   nlpb=" << this->nlpb_grid
          << ", nk=" << this->nk_grid;
      MACOM_LOG_ERROR("MACOMState", error_msg.str());
      throw std::runtime_error(error_msg.str());
    }
    // MACOM_LOG_INFO("MACOMState",
    //                "State file dimensions match geometry dimensions.");

    // Step 3: Load variable arrays
    loadVariableArrays(ncFile, vars_to_load_from_config);

    // Populate variableNames_ and activeVariable_
    this->variableNames_ = vars_to_load_from_config;
    if (!this->variableNames_.empty()) {
      this->activeVariable_ = this->variableNames_[0];
      // MACOM_LOG_INFO("MACOMState",
      //                "Active variable set to: " + this->activeVariable_);
    } else {
      MACOM_LOG_WARNING("MACOMState", "No variables loaded or specified.");
    }

    initialized_ = true;
    // MACOM_LOG_INFO("MACOMState",
    //                "MACOMState initialized successfully from " + inputFile_);

  } catch (const netCDF::exceptions::NcException& e) {
    MACOM_LOG_ERROR(
        "MACOMState",
        "NetCDF error while initializing MACOMState: " + std::string(e.what()));
    throw std::runtime_error("NetCDF error while initializing MACOMState: " +
                             std::string(e.what()));
  } catch (const std::exception& e) {
    MACOM_LOG_ERROR("MACOMState",
                    "Error initializing MACOMState: " + std::string(e.what()));
    throw std::runtime_error("Error initializing MACOMState: " +
                             std::string(e.what()));
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
      config_, geometry_);
  cloned->initialized_ = this->initialized_;
  cloned->inputFile_ = this->inputFile_;
  cloned->variables_ = this->variables_;
  cloned->activeVariable_ = this->activeVariable_;
  cloned->variableNames_ = this->variableNames_;
  cloned->dimensions_ = this->dimensions_;
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

    // 输出样本数据：从变量数据中间取10个值
    std::stringstream sample_ss;
    sample_ss << "Sample values for " << variable_name << ": ";

    size_t total_size = variables_[variable_name].size();
    size_t mid_point = total_size / 2;
    size_t start_idx = (total_size <= 10) ? 0 : (mid_point - 5);
    size_t end_idx = std::min(start_idx + 10, total_size);

    for (size_t i = start_idx; i < end_idx; ++i) {
      sample_ss << variables_[variable_name][i];
      if (i < end_idx - 1) sample_ss << ", ";
    }

    // 使用新的日志系统输出样本数据
    MACOM_LOG_INFO("MACOMState", sample_ss.str());

    // 也可以添加简单的统计信息
    if (!variables_[variable_name].empty()) {
      double min_val = variables_[variable_name][0];
      double max_val = variables_[variable_name][0];
      double sum = 0.0;

      for (const auto& val : variables_[variable_name]) {
        min_val = std::min(min_val, val);
        max_val = std::max(max_val, val);
        sum += val;
      }

      double avg = sum / total_size;
      std::stringstream stats_ss;
      stats_ss << "Statistics for " << variable_name << ": min=" << min_val
               << ", max=" << max_val << ", avg=" << avg
               << ", count=" << total_size;

      MACOM_LOG_INFO("MACOMState", stats_ss.str());
    }
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

}  // namespace metada::backends::macom