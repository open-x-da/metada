/**
 * @file MACOMState.hpp
 * @brief MACOM state backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <algorithm>  // For std::equal, std::find
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>    // For std::inner_product (potentially for dot product)
#include <stdexcept>  // For std::runtime_error
#include <string>
#include <unordered_map>
#include <vector>

#include "MACOMGeometry.hpp"

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
template <typename ConfigBackend>
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
  explicit MACOMState(const ConfigBackend& config);

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
  void setActiveVariable([[maybe_unused]] const std::string& name);

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
   * @brief Add another state to this one
   *
   * @param other State to add
   * @throws std::runtime_error If states are incompatible
   */
  void add([[maybe_unused]] const MACOMState& other);

  /**
   * @brief Subtract another state from this one
   *
   * @param other State to subtract
   * @throws std::runtime_error If states are incompatible
   */
  void subtract([[maybe_unused]] const MACOMState& other);

  /**
   * @brief Multiply this state by a scalar
   *
   * @param scalar Value to multiply by
   */
  void multiply([[maybe_unused]] double scalar);

  /**
   * @brief Check if state is properly initialized
   *
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Calculate the norm of this state
   *
   * @return Norm of the state
   */
  double norm() const {
    double sum_sq = 0.0;
    for (const auto& pair : variables_) {
      for (double val : pair.second) {
        sum_sq += val * val;
      }
    }
    return std::sqrt(sum_sq);
  }

  /**
   * @brief 计算与另一个状态的点积
   *
   * @param other 另一个状态
   * @return double 点积结果
   */
  double dot(const MACOMState& other) const {
    if (variableNames_.size() != other.variableNames_.size()) {
      // Or handle more gracefully if some variables might be missing
      throw std::runtime_error(
          "MACOMState::dot: States have different number of variables.");
    }
    double total_dot_product = 0.0;
    for (const auto& name : variableNames_) {
      const auto& vec_this = this->variables_.at(name);
      const auto& vec_other = other.variables_.at(name);
      if (vec_this.size() != vec_other.size()) {
        throw std::runtime_error("MACOMState::dot: Variable '" + name +
                                 "' has different sizes.");
      }
      total_dot_product += std::inner_product(vec_this.begin(), vec_this.end(),
                                              vec_other.begin(), 0.0);
    }
    return total_dot_product;
  }

  /**
   * @brief 检查两个状态是否相等
   *
   * @param other 要比较的另一个状态
   * @return bool 如果相等返回true，否则返回false
   */
  bool equals(const MACOMState& other) const {
    if (initialized_ != other.initialized_ ||
        activeVariable_ != other.activeVariable_ ||
        variableNames_.size() != other.variableNames_.size() ||
        variables_.size() != other.variables_.size() ||
        dimensions_.size() != other.dimensions_.size()) {
      return false;
    }

    // Check variable names
    std::vector<std::string> vn1 = variableNames_;
    std::vector<std::string> vn2 = other.variableNames_;
    std::sort(vn1.begin(), vn1.end());
    std::sort(vn2.begin(), vn2.end());
    if (vn1 != vn2) return false;

    for (const auto& name : variableNames_) {
      if (variables_.at(name).size() != other.variables_.at(name).size() ||
          dimensions_.at(name) != other.dimensions_.at(name)) {
        return false;
      }
      // Approx equals for floating point data
      const auto& data_this = variables_.at(name);
      const auto& data_other = other.variables_.at(name);
      for (size_t i = 0; i < data_this.size(); ++i) {
        if (std::abs(data_this[i] - data_other[i]) >
            1e-9) {  // Tolerance for float comparison
          return false;
        }
      }
    }
    return true;
  }

 private:
  using Geometry_MACOM = MACOMGeometry<ConfigBackend>;
  const Geometry_MACOM* geometry_ptr_ = nullptr;
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
   * @brief Initialize grid from NetCDF file
   *
   * @param filename Path to the NetCDF grid file
   */
  void loadVariableData(const std::string& filename,
                        const std::vector<std::string>& variables);

  // Variable dimensions
  std::size_t nlpb_ = 0;  // Number of grid points
  std::size_t nk_ = 0;    // Number of vertical levels

  std::size_t nlpb_grid = 0;  // Number of grid points
  std::size_t nk_grid = 0;    // Number of vertical levels
  const ConfigBackend* config_ptr_ = nullptr;

  // Variable data
  std::vector<double> u;  // u-velocity
  std::vector<double> v;  // v-velocity
  std::vector<double> t;  // temperature
  std::vector<double> s;  // salinity
  std::vector<double> w;  // w-velocity

  // Reference to configuration
  const ConfigBackend& config_;

  // State information
  bool initialized_ = false;
  std::string inputFile_;                   // Input data file path
  std::string timestamp_;                   // Timestamp of the data
  std::vector<std::string> variableNames_;  // Available variables
  std::string activeVariable_;              // Currently active variable

  // Data storage (this would be filled by the actual implementation)
  std::unordered_map<std::string, std::vector<double>> variables_;
  std::unordered_map<std::string, std::vector<size_t>> dimensions_;
};  // namespace metada::backends::macom

// ConfigBackend constructor implementation
template <typename ConfigBackend>
MACOMState<ConfigBackend>::MACOMState(const ConfigBackend& config)
    : config_(config), initialized_(false) {
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
template <typename ConfigBackend>
MACOMState<ConfigBackend>::MACOMState(MACOMState&& other) noexcept
    : config_(other.config_),
      initialized_(other.initialized_),
      inputFile_(std::move(other.inputFile_)),
      variableNames_(std::move(other.variableNames_)),
      activeVariable_(std::move(other.activeVariable_)),
      variables_(std::move(other.variables_)),
      dimensions_(std::move(other.dimensions_)) {
  other.initialized_ = false;
}

template <typename ConfigBackend>
MACOMState<ConfigBackend>& MACOMState<ConfigBackend>::operator=(
    MACOMState&& other) noexcept {
  if (this != &other) {
    initialized_ = other.initialized_;
    inputFile_ = std::move(other.inputFile_);
    variableNames_ = std::move(other.variableNames_);
    activeVariable_ = std::move(other.activeVariable_);
    variables_ = std::move(other.variables_);
    dimensions_ = std::move(other.dimensions_);

    other.initialized_ = false;
  }
  return *this;
}

template <typename ConfigBackend>
std::unique_ptr<MACOMState<ConfigBackend>> MACOMState<ConfigBackend>::clone()
    const {
  auto cloned = std::make_unique<MACOMState>(config_);
  cloned->initialized_ = this->initialized_;
  cloned->inputFile_ = this->inputFile_;
  cloned->variableNames_ = this->variableNames_;
  cloned->activeVariable_ = this->activeVariable_;
  cloned->variables_ = this->variables_;
  cloned->dimensions_ = this->dimensions_;
  return cloned;
}

template <typename ConfigBackend>
void* MACOMState<ConfigBackend>::getData() {
  // 基本实现框架
  return nullptr;
}

template <typename ConfigBackend>
const void* MACOMState<ConfigBackend>::getData() const {
  // 基本实现框架
  return nullptr;
}

template <typename ConfigBackend>
const std::string& MACOMState<ConfigBackend>::getActiveVariable() const {
  return activeVariable_;
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::setActiveVariable(
    [[maybe_unused]] const std::string& name) {
  // 实现框架
}

template <typename ConfigBackend>
const std::vector<std::string>& MACOMState<ConfigBackend>::getVariableNames()
    const {
  return variableNames_;
}

template <typename ConfigBackend>
const std::vector<size_t>& MACOMState<ConfigBackend>::getDimensions(
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

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::zero() {
  // 基本实现框架
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::add([[maybe_unused]] const MACOMState& other) {
  // 实现框架
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::subtract(
    [[maybe_unused]] const MACOMState& other) {
  // 实现框架
}

template <typename ConfigBackend>
void MACOMState<ConfigBackend>::multiply([[maybe_unused]] double scalar) {
  // 实现框架
}

// Implementation of loadVariableDimensions
template <typename ConfigBackend>
void MACOMState<ConfigBackend>::loadVariableDimensions(netCDF::NcFile& ncFile) {
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
template <typename ConfigBackend>
void MACOMState<ConfigBackend>::loadVariableArrays(
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
template <typename ConfigBackend>
void MACOMState<ConfigBackend>::loadVariableData(
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