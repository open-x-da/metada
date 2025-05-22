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
#include <string>
#include <unordered_map>
#include <vector>

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
  // /**
  //  * @brief Set the associated geometry object
  //  *
  //  * @param geometry Pointer to a MACOMGeometry object
  //  */
  // void setGeometry(const Geometry_MACOM* geometry) {
  //   if (!geometry) {
  //     throw std::runtime_error("setGeometry: null pointer passed");
  //   }
  //   if (!geometry->isInitialized()) {
  //     throw std::runtime_error("setGeometry: geometry not initialized");
  //   }
  //   geometry_ptr_ = geometry;
  //   // now pull the numbers out of the geometry
  //   nlpb_ = geometry_ptr_->getNlpb();
  //   nk_ = geometry_ptr_->getNk();
  //   // (if you need nkp1 you can grab that too)
  //   std::cout << "State bound to geometry: nlpb=" << nlpb_ << " nk=" << nk_
  //             << "\n";
  // }

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
  void multiply();

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
  void loadVariableData(const std::string& filename,
                        const std::vector<std::string>& variables);

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
  std::string input_filename = config.Get("input_file").asString();
  if (input_filename.empty()) {
    throw std::runtime_error(
        "MACOM state input file path not specified in configuration");
  }
  inputFile_ = input_filename;

  std::string timestamp_ = config.Get("timestamp").asString();

  // Get variables to load from config
  std::vector<std::string> variables;

  // Try to get variables list from config, otherwise use defaults
  try {
    variables = config.Get("variables").asVectorString();
  } catch (const std::exception&) {
    std::cerr << "Warning: No variables specified in configuration"
              << std::endl;
  }

  std::cout << "nlpb_grid = " << nlpb_grid << ", nk_grid = " << nk_grid
            << std::endl;

  loadVariableData(input_filename, variables);

  std::cout << "MACOMState initialized = " << initialized_ << std::endl;

  // variableNames_ = variables;

  // // Set the active variable to the first one if available
  // if (!variables.empty()) {
  //   activeVariable_ = variables[0];
  // }
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
    [[maybe_unused]] const MACOMState& other) {
  // 实现框架
}

template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::subtract(
    const MACOMState<ConfigBackend, GeometryBackend>& other) {
  // Subtract each variable
  if (!isCompatible(other)) {
    throw std::runtime_error("States are incompatible for subtraction");
  }
}

template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::multiply() {
  // Multiply each variable by the scalar
  // for (auto& [name, data] : variables_) {
  //   data *= scalar;
  // }
}

// Implementation of loadVariableDimensions
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::loadVariableDimensions(
    netCDF::NcFile& ncFile) {
  auto readDimension = [&ncFile](const std::string& name, std::size_t& value) {
    std::cout << "Attempting to read dimension: " << name << std::endl;
    auto dim = ncFile.getDim(name);
    if (dim.isNull()) {
      throw std::runtime_error("Dimension '" + name +
                               "' not found in grid file");
    }
    std::cout << "Found dimension: " << name << std::endl;
    value = dim.getSize();
    std::cout << "Read size for " << name << ": " << value << std::endl;
  };

  // Read all variable dimensions
  readDimension("nlpb", nlpb_);
  readDimension("nk", nk_);

  std::cout << "Loaded variable dimensions:" << std::endl;
  std::cout << "  nlpb=" << nlpb_ << ", nk=" << nk_ << std::endl;
}

// Implementation of loadVariableArrays
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::loadVariableArrays(
    netCDF::NcFile& ncFile, const std::vector<std::string>& variables) {
  for (const auto& variable : variables) {
    if (variable == "u") {
      u.resize(nlpb_ * nk_);
    } else if (variable == "v") {
      v.resize(nlpb_ * nk_);
    } else if (variable == "t") {
      t.resize(nlpb_ * nk_);
    } else if (variable == "s") {
      s.resize(nlpb_ * nk_);
    } else if (variable == "w") {
      w.resize(nlpb_ * nk_);
    }
  }

  // Read variables
  auto readVar = [&ncFile](const std::string& name, std::vector<double>& data) {
    auto var = ncFile.getVar(name);
    if (var.isNull()) {
      throw std::runtime_error("Variable '" + name +
                               "' not found in initial data file");
    }
    var.getVar(data.data());
  };

  // Read grid data
  for (const auto& variable : variables) {
    readVar(variable, variables_[variable]);
  }

  std::cout << "Loaded variable arrays" << std::endl;
}

// Implementation of loadVariableData
template <typename ConfigBackend, typename GeometryBackend>
void MACOMState<ConfigBackend, GeometryBackend>::loadVariableData(
    const std::string& filename, const std::vector<std::string>& variables) {
  try {
    std::cout << "Attempting to open NetCDF file: " << filename << std::endl;

    // Open the NetCDF file
    netCDF::NcFile ncFile(filename, netCDF::NcFile::read);

    if (ncFile.isNull()) {
      throw std::runtime_error("Failed to open NetCDF file: " + filename);
    }

    std::cout << "Successfully opened NetCDF file" << std::endl;

    // Step 1: Load grid dimensions
    loadVariableDimensions(ncFile);

    // Step 2: Load grid arrays
    loadVariableArrays(ncFile, variables);

    // Close the file
    ncFile.close();

    initialized_ = true;

    std::cout << "Successfully initialized variable data from " << filename
              << std::endl;

  } catch (const netCDF::exceptions::NcException& e) {
    throw std::runtime_error("NetCDF error while reading variable data file: " +
                             std::string(e.what()));
  } catch (const std::exception& e) {
    throw std::runtime_error("Error initializing variable data: " +
                             std::string(e.what()));
  }
}

}  // namespace metada::backends::macom