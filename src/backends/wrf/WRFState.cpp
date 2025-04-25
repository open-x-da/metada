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
#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>

namespace metada::backends::wrf {

// Constructor implementation with ConfigBackend
template <typename ConfigBackend>
WRFState::WRFState(const ConfigBackend& config)
    : wrfFilename_(config.getString("wrf.input_file", "")),
      timestamp_(config.getString("wrf.timestamp", "0000-00-00_00:00:00")),
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
    variables = config.getStringVector("wrf.variables");
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
WRFState::WRFState(WRFState&& other) noexcept
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
WRFState& WRFState::operator=(WRFState&& other) noexcept {
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
std::unique_ptr<WRFState> WRFState::clone() const {
  // Use copy constructor for cloning (private constructor would be better in a
  // real implementation)
  return std::make_unique<WRFState>(*this);
}

// Data access implementation
void* WRFState::getData() {
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
const void* WRFState::getData() const {
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
const std::string& WRFState::getActiveVariable() const {
  return activeVariable_;
}

// Set active variable implementation
void WRFState::setActiveVariable(const std::string& name) {
  if (variables_.find(name) == variables_.end()) {
    throw std::out_of_range("Variable not found: " + name);
  }
  activeVariable_ = name;
}

// Get variable names implementation
const std::vector<std::string>& WRFState::getVariableNames() const {
  return variableNames_;
}

// Get dimensions implementation
const std::vector<size_t>& WRFState::getDimensions(
    const std::string& name) const {
  try {
    return dimensions_.at(name);
  } catch (const std::out_of_range&) {
    throw std::out_of_range("Variable not found: " + name);
  }
}

// Zero implementation
void WRFState::zero() {
  for (auto& [name, data] : variables_) {
    data.fill(0.0);
  }
}

// Dot product implementation
double WRFState::dot(const WRFState& other) const {
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
double WRFState::norm() const {
  double sumSquares = 0.0;

  // Sum squares of all variables
  for (const auto& [name, data] : variables_) {
    sumSquares += xt::sum(xt::square(data))();
  }

  return std::sqrt(sumSquares);
}

// Equals implementation
bool WRFState::equals(const WRFState& other) const {
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
void WRFState::add(const WRFState& other) {
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
void WRFState::subtract(const WRFState& other) {
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
void WRFState::multiply(double scalar) {
  // Multiply each variable by the scalar
  for (auto& [name, data] : variables_) {
    data *= scalar;
  }
}

// Is initialized implementation
bool WRFState::isInitialized() const {
  return initialized_;
}

// Set data implementation
void WRFState::setData(const std::vector<double>& values,
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
const xt::xarray<double>& WRFState::getData(
    const std::string& variableName) const {
  // Use active variable if none specified
  std::string varName = variableName.empty() ? activeVariable_ : variableName;

  if (varName.empty() || variables_.find(varName) == variables_.end()) {
    throw std::runtime_error("Variable not found: " + varName);
  }

  return variables_.at(varName);
}

// Private helper to load state data
void WRFState::loadStateData(const std::string& filename,
                             const std::vector<std::string>& variables,
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
bool WRFState::isCompatible(const WRFState& other) const {
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

// Explicit template instantiations for known config backend types
// Add additional instantiations as needed for different config backends
template WRFState::WRFState(
    const metada::backends::yaml::YAMLConfigBackend& config);

}  // namespace metada::backends::wrf