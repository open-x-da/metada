/**
 * @file WRFState.cpp
 * @brief Implementation of the WRF state backend
 * @ingroup backends
 * @author Metada Framework Team
 */

#include "WRFState.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#if defined(_WIN32) || defined(WIN32)
#include <xtensor/core/xmath.hpp>
#include <xtensor/views/xview.hpp>
#else
#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#endif

namespace metada::backends::wrf {

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
WRFState<ConfigBackend>::WRFState(const ConfigBackend& config)
    : config_(config),
      wrfFilename_(config.Get("wrf.input_file").asString()),
      timestamp_(config.Get("wrf.timestamp").asString()),
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
    variables = config.Get("wrf.variables").as;
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
    : wrfFilename_(std::move(other.wrfFilename_)),
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

      size_t nx = dim_west_east.getSize();
      size_t ny = dim_south_north.getSize();
      size_t nz = dim_bottom_top.getSize();

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

/*
 * TODO: Template instantiations for config backends need to be fixed
 * There is a namespace/include issue with YamlConfig that needs to be addressed
 */

}  // namespace metada::backends::wrf