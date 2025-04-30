/**
 * @file WRFState.hpp
 * @brief WRF state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

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

namespace metada::backends::wrf {

/**
 * @brief WRF state backend implementation
 *
 * @details
 * This class implements a state backend for the WRF model. It manages
 * meteorological state variables stored in NetCDF files and provides
 * operations required by the State adapter.
 */
template <typename ConfigBackend>
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
   * @brief Constructor that takes a configuration backend
   *
   * @param config Configuration containing WRF file path and variables
   */
  explicit WRFState(const ConfigBackend& config);

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

  /**
   * @brief Set data for a specific variable
   *
   * @param values Vector of values to set
   * @param variableName Optional variable name (uses active variable if not
   * specified)
   * @throws std::runtime_error If variable doesn't exist or dimensions don't
   * match
   */
  void setData(const std::vector<double>& values,
               const std::string& variableName = "");

  /**
   * @brief Get data for a specific variable
   *
   * @param variableName Optional variable name (uses active variable if not
   * specified)
   * @return const xt::xarray<double>& Reference to the data array
   * @throws std::runtime_error If variable doesn't exist
   */
  const xt::xarray<double>& getData(const std::string& variableName = "") const;

 private:
  /**
   * @brief Load state data from the WRF NetCDF file
   *
   * @param filename Path to the WRF NetCDF file
   * @param variables List of variables to load
   * @param timestamp Timestamp to read from the file
   */
  void loadStateData(const std::string& filename,
                     const std::vector<std::string>& variables,
                     const std::string& timestamp);

  /**
   * @brief Verify that two states have compatible variables and dimensions
   *
   * @param other State to check compatibility with
   * @return True if compatible, false otherwise
   */
  bool isCompatible(const WRFState& other) const;

  const ConfigBackend& config_;

  // WRF NetCDF file information
  std::string wrfFilename_;
  std::string timestamp_;
  bool initialized_ = false;

  // State data
  std::unordered_map<std::string, xt::xarray<double>> variables_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
  std::vector<std::string> variableNames_;
  std::string activeVariable_;  // Currently active variable for data access
};

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
WRFState<ConfigBackend>::WRFState(const ConfigBackend& config)
    : config_(config),
      wrfFilename_(config.Get("input_file").asString()),
      timestamp_(config.Get("timestamp").asString()),
      initialized_(false) {
  if (wrfFilename_.empty()) {
    throw std::runtime_error(
        "WRF input file path not specified in configuration");
  }

  // Get variables to load from config
  std::vector<std::string> defaultVars = {"T", "U", "V", "QVAPOR", "P"};
  std::vector<std::string> variables;

  // Try to get variables list from config, otherwise use defaults
  try {
    variables = config.Get("wrf.variables").asVectorString();
  } catch (const std::exception&) {
    variables = defaultVars;
  }

  // Set the active variable to the first one if available
  if (!variables.empty()) {
    activeVariable_ = variables[0];
  }

  // Load state data from WRF NetCDF file
  loadStateData(wrfFilename_, variables, timestamp_);
  initialized_ = true;
}

// Move constructor implementation
template <typename ConfigBackend>
WRFState<ConfigBackend>::WRFState(WRFState<ConfigBackend>&& other) noexcept
    : config_(std::move(other.config_)),
      wrfFilename_(std::move(other.wrfFilename_)),
      timestamp_(std::move(other.timestamp_)),
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
template <typename ConfigBackend>
WRFState<ConfigBackend>& WRFState<ConfigBackend>::operator=(
    WRFState<ConfigBackend>&& other) noexcept {
  if (this != &other) {
    wrfFilename_ = std::move(other.wrfFilename_);
    timestamp_ = std::move(other.timestamp_);
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
template <typename ConfigBackend>
std::unique_ptr<WRFState<ConfigBackend>> WRFState<ConfigBackend>::clone()
    const {
  auto cloned = std::make_unique<WRFState<ConfigBackend>>(config_);
  // Copy all the state data to the cloned object
  cloned->wrfFilename_ = this->wrfFilename_;
  cloned->timestamp_ = this->timestamp_;
  cloned->initialized_ = this->initialized_;
  cloned->variables_ = this->variables_;
  cloned->dimensions_ = this->dimensions_;
  cloned->variableNames_ = this->variableNames_;
  cloned->activeVariable_ = this->activeVariable_;
  return cloned;
}

// Data access implementation
template <typename ConfigBackend>
void* WRFState<ConfigBackend>::getData() {
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
template <typename ConfigBackend>
const void* WRFState<ConfigBackend>::getData() const {
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
template <typename ConfigBackend>
const std::string& WRFState<ConfigBackend>::getActiveVariable() const {
  return activeVariable_;
}

// Set active variable implementation
template <typename ConfigBackend>
void WRFState<ConfigBackend>::setActiveVariable(const std::string& name) {
  if (variables_.find(name) == variables_.end()) {
    throw std::out_of_range("Variable not found: " + name);
  }
  activeVariable_ = name;
}

// Get variable names implementation
template <typename ConfigBackend>
const std::vector<std::string>& WRFState<ConfigBackend>::getVariableNames()
    const {
  return variableNames_;
}

// Get dimensions implementation
template <typename ConfigBackend>
const std::vector<size_t>& WRFState<ConfigBackend>::getDimensions(
    const std::string& name) const {
  try {
    return dimensions_.at(name);
  } catch (const std::out_of_range&) {
    throw std::out_of_range("Variable not found: " + name);
  }
}

// Zero implementation
template <typename ConfigBackend>
void WRFState<ConfigBackend>::zero() {
  for (auto& [name, data] : variables_) {
    data.fill(0.0);
  }
}

// Dot product implementation
template <typename ConfigBackend>
double WRFState<ConfigBackend>::dot(
    const WRFState<ConfigBackend>& other) const {
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
template <typename ConfigBackend>
double WRFState<ConfigBackend>::norm() const {
  double sumSquares = 0.0;

  // Sum squares of all variables
  for (const auto& [name, data] : variables_) {
    sumSquares += xt::sum(xt::square(data))();
  }

  return std::sqrt(sumSquares);
}

// Equals implementation
template <typename ConfigBackend>
bool WRFState<ConfigBackend>::equals(
    const WRFState<ConfigBackend>& other) const {
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
template <typename ConfigBackend>
void WRFState<ConfigBackend>::add(const WRFState<ConfigBackend>& other) {
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
template <typename ConfigBackend>
void WRFState<ConfigBackend>::subtract(const WRFState<ConfigBackend>& other) {
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
template <typename ConfigBackend>
void WRFState<ConfigBackend>::multiply(double scalar) {
  // Multiply each variable by the scalar
  for (auto& [name, data] : variables_) {
    data *= scalar;
  }
}

// Is initialized implementation
template <typename ConfigBackend>
bool WRFState<ConfigBackend>::isInitialized() const {
  return initialized_;
}

// Set data implementation
template <typename ConfigBackend>
void WRFState<ConfigBackend>::setData(const std::vector<double>& values,
                                      const std::string& variableName) {
  // Use active variable if none specified
  std::string varName = variableName.empty() ? activeVariable_ : variableName;

  if (varName.empty() || variables_.find(varName) == variables_.end()) {
    throw std::runtime_error("Variable not found: " + varName);
  }

  auto& data = variables_.at(varName);
  const auto& dims = dimensions_.at(varName);

  // Calculate total size from dimensions
  size_t totalSize = 1;
  for (size_t dim : dims) {
    totalSize *= dim;
  }

  // Check if input data size matches expected size
  if (values.size() != totalSize) {
    throw std::runtime_error("Data size mismatch for variable: " + varName);
  }

  // Copy data to the array
  std::copy(values.begin(), values.end(), data.begin());
}

// Get data implementation
template <typename ConfigBackend>
const xt::xarray<double>& WRFState<ConfigBackend>::getData(
    const std::string& variableName) const {
  // Use active variable if none specified
  std::string varName = variableName.empty() ? activeVariable_ : variableName;

  if (varName.empty() || variables_.find(varName) == variables_.end()) {
    throw std::runtime_error("Variable not found: " + varName);
  }

  return variables_.at(varName);
}

// Private helper to load state data
template <typename ConfigBackend>
void WRFState<ConfigBackend>::loadStateData(
    const std::string& filename, const std::vector<std::string>& variables,
    const std::string& timestamp) {
  try {
    // Open NetCDF file
    netCDF::NcFile wrf_file(filename, netCDF::NcFile::read);

    if (!wrf_file.isNull()) {
      // Read dimensions
      auto dim_west_east = wrf_file.getDim("west_east");
      auto dim_south_north = wrf_file.getDim("south_north");
      auto dim_bottom_top = wrf_file.getDim("bottom_top");

      if (dim_west_east.isNull() || dim_south_north.isNull() ||
          dim_bottom_top.isNull()) {
        throw std::runtime_error("Missing required dimensions in WRF file");
      }

      // Find time index for the requested timestamp
      size_t time_idx = 0;  // Default to first time step
      auto times_var = wrf_file.getVar("Times");
      if (!times_var.isNull()) {
        // Code to find matching timestamp would go here
        // For simplicity, we'll use index 0
      }

      // Clear any existing data
      variables_.clear();
      dimensions_.clear();
      variableNames_.clear();

      // Load each requested variable
      for (const auto& varName : variables) {
        auto var = wrf_file.getVar(varName);

        if (!var.isNull()) {
          // Get variable dimensions
          const auto& varDims = var.getDims();
          std::vector<size_t> dims;

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
            start[0] = time_idx;
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
template <typename ConfigBackend>
bool WRFState<ConfigBackend>::isCompatible(
    const WRFState<ConfigBackend>& other) const {
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